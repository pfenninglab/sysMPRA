#!/usr/bin/python

#Use multiple alignment files to get the correct alignment between the cbp mouse enhancers and human

import fileinput
import subprocess
import sys
import os
from os import path
import glob
import re
import math;
import io;
import gzip;

from datetime import datetime



#import urllib.request
#from urllib.request import urlopen
#import xml.etree.ElementTree as ET
#from xml.etree.ElementTree import parse

#from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
#from Bio.Alphabet.IUPAC import unambiguous_dna, ambiguous_dna
#from Bio import SeqIO
#from Bio.SeqRecord import SeqRecord


#perl -p -e 's/\r//g' startingAnnot/ensembl/targetProteins.ensp.2.1.1.txt > startingAnnot/ensembl/targetProteins.ensp.2.1.2.txt

#Get information from deSeq2 about significance of enhancer difference
#Intersect with mouse enhancers (must be greater than 50bp to count)
#Get the 0.9 quantiles over each region
#Annotate the enhancer with the closest gene
#Annotate the file with chromatin state
#concat.allInfo.2.1.py brings everything together into one file

#awk '{a[$4]++} END {print length(a)}' mergeMark.3.3.bed, print the number of unique original ids

#pip3 install biopython --user

################################################
#   Initialization                  #
################################################

#initBedStr = "use BEDTools";
#print initBedStr
#subprocess.Popen(initBedStr, shell=True).wait();

#srun --partition=pfen2 --mem=8G --pty bash
#module load python36
#cd /home/apfennin/projects/covidTest/orthGene/


################################################
#   Functions                         #
################################################





################################################
#   Version control  #
################################################
#Version 1.1.1 - Does basic quantification of MPRA barcodes
#Version 1.1.2 - Updated to loop through all fastq files


################################################
#   PrepCode   #
################################################

#Get a subset of the GTF file that only include CDS entries
#OLd
#awk -vOFS='\t' -vFS='\t' '{if ($3 == "CDS") {print $0}}' startingAnnot/GRCh38_latest_genomic.gtf > startingAnnot/GRCh38_latest_genomic_onlyCds.gtf

#Copy test files to run

#head /projects/MPRA/MPRA/MPRAi/MPRAi_v2_NovaSeq_Genewiz/1-3-R-11/1-3-R-11_R1_001.fastq -n 1000 > input/1-3-R-11_R1_001.test.fastq

#head /projects/MPRA/MPRA/MPRAi/MPRAi_v2_NovaSeq_Genewiz/1-3-D-19/1-3-D-19_R1_001.fastq -n 1000 > input/1-3-D-19_R1_001.test.fastq

#grep TCAGAACATGAGACTC /projects/MPRA/MPRA/MPRAi/MPRAi_batch1mice_1-208906698/FASTQ_Generation_2020-11-06_19_51_52Z-338176598/MPRA_Array1_sequences

################################################
#   Command Line Run   #
################################################

#Location
#/projects/MPRA/andreas/21_06_arrayProc

#module load python36
#cd /projects/MPRA/andreas/21_06_arrayProc


#python3 code/orthoGene_200M.2.1.1.py startingAnnot/GRCh38_latest_genomic_onlyCds.gtf startingAnnot/refSeq2Chr.04-06-20.1.txt NP_068576.1 zoonomia/zoonomiaGuide_4-6-20_2.csv /data/pfenninggroup/align/mam200/alignment_files/200m-v1.hal xSpecProtBed

#python3 arrayProc.1.1.1.py input/Sample_index.txt /projects/MPRA/MPRA/MPRAi/MPRAi_v2_NovaSeq_Genewiz/ /projects/MPRA/MPRA/MPRAi/MPRAi_batch1mice_1-208906698/FASTQ_Generation_2020-11-06_19_51_52Z-338176598/MPRA_Array1_sequences

#python3 arrayProc.1.1.2.py input/Sample_index.txt /projects/MPRA/MPRA/MPRAi/MPRAi_v2_NovaSeq_Genewiz/ /projects/MPRA/MPRA/MPRAi/MPRAi_batch1mice_1-208906698/FASTQ_Generation_2020-11-06_19_51_52Z-338176598/MPRA_Array1_sequences

################################################
#   Parameters   #
################################################


#Read in
indexMapFn = str(sys.argv[1]); #The file that contains index information
seqDir = str(sys.argv[2]); #The directory that contains the fastqs
barcodeMapFn = str(sys.argv[3]); #The map between barcodes and enhancers



################################################
#   Program Control   #
################################################


#Get genome stats
loopDetailedB = False;

loopFastB = True;

writeCountsB = True;


################################################
#   Read in basic information from the annotation files  #
################################################

######### Create a list of barcodes and their associated index  ###############

rcCharD = {};
rcCharD['A'] = 'T';
rcCharD['T'] = 'A';
rcCharD['C'] = 'G';
rcCharD['G'] = 'C';
rcCharD['N'] = 'N';

def rcSeqFc(seq) :
    seq2 = "";
    for curChar in seq :
        seq2 = rcCharD[curChar] + seq2;
    #seq3 = seq2[::-1]
    return(seq2);

barcode2indexD = {}; #Maps barcodes to the index number
barcodeRc2indexD = {}; #Maps barcode revComp to the index number

curRow = 0;
curIndex = 0;

for curLine in fileinput.input([barcodeMapFn]):
    curLineP = curLine.rstrip().split("\t");
    if(curRow > 0 and curLineP[0] != '') :
        barcode2indexD[curLineP[3]] = curIndex;
        barcodeRc2indexD[rcSeqFc(curLineP[3])] = curIndex;
        curIndex = curIndex + 1;

        #print(rcSeqFc(curLineP[3]));


    curRow = curRow + 1;

barcodeArrayLen = curIndex;


#print(barcode2indexD["GCCTATTAACTCACTA"]);
#print("Results - " + str(barcodeRc2indexD["GAGTCTCATGTTCTGA"]));

######### Read in the sample information  ###############

samp2numD = {}; #maps samples to how many there currently are
base2InfoD = {}; #maps the base numbers to sample and replicate information
                 #

for curLine in fileinput.input([indexMapFn]): #Read in list of things for gene distance
    curLineP = curLine.rstrip().split("\t");
    if curLineP[1] in samp2numD :
        samp2numD[curLineP[1]] = samp2numD[curLineP[1]] + 1;
    else :
        samp2numD[curLineP[1]] = 1;

    curId = curLineP[1] + "_" + str(samp2numD[curLineP[1]]);

    base2InfoD[curLineP[0]] = [curId,curLineP[0],curLineP[1],samp2numD[curLineP[1]]];


#print(samp2numD["1-2"]);


######### Get the set of all sequenced experiments  ###############


fastqDirsV = os.listdir(seqDir);
#print(fastqDirsV);

dir2InfoD = {}; #id, directory, (DNA/RNA), baseId, base, tissue, replicate number, fastqR1, fastqR2
#indexed by id
idSetV = []; #A list of the IDs included

for curDir in  fastqDirsV :

    curDirSplitV = curDir.split("-");

    if len(curDirSplitV) > 1 :

        curBase = curDirSplitV[0] + "-" + curDirSplitV[1];
        curBaseV = base2InfoD[curBase];

        curMolec = curDirSplitV[2];

        curFn1 = seqDir + curDir + "/" + curDir + "_R1_001.fastq";
        curFn2 = seqDir + curDir + "/" + "_R1_001.fastq"


        curId = curBaseV[0] + "_" + curMolec;
        curBaseV = [curId, curDir, curMolec]  + curBaseV + [curFn1,curFn2];

        dir2InfoD[curId] = curBaseV;

        idSetV.append(curId);


print(curBaseV);


################################################
#   Loop through the fast q file performing calcs   #
################################################

#For R1, reverse complement




############### Set parameters ###########################

reSiteRc1 = "TCTAGAGGTACC"; #Perfect match
reSiteRc2 = "TCTAGACGTACC"; #Second match
reSiteRc3 = "TCTAGAAGTACC"; #Third match
reSiteRc4 = "TCTAGATGTACC"; #Fourth match

reSiteV = [reSiteRc1,reSiteRc2,reSiteRc3,reSiteRc4];


#primer1Rc = "CGACGCTCTTCCGATCT";


############### Iterate through rows of a fastq file - Detailed #############

statSummaryM = [];

for curExp in idSetV :

    curFastQFn = dir2InfoD[curExp][7];
    #curFastQFn = "input/1-3-R-11_R1_001.test.fastq";
    #curFastQFn = "input/1-3-D-19_R1_001.test.fastq";
    #curCountFn = "counts/counts.STR_1_D_test.1.1.csv";
    curCountFn = "counts/counts." + dir2InfoD[curExp][0] + ".1.1.csv";

    print(curFastQFn);
    print(curCountFn);

    barcodeCountV = [0] * barcodeArrayLen;

    if loopDetailedB :

        reSiteLocV = []; #Location of restriction enzyme site matches
        reSiteLocLaxV = []; #Location of restriction enzyme site matches, relaxed

        eightyFourD = {}; #Dictionary counting sequences at position 84;

        bcStatusV = []; #Vector with whether a good barcode was found
                        #-1 no RE site; 0 = no bc hit; 1=good bc hit


        curRow = 0;


        for curLine in fileinput.input([curFastQFn]):

            if curRow % 4 == 1 : #Only take certain rows of the fastq
                #print(curLine);

                #Restriction Enzyme site match
                curResLoc = curLine.find(reSiteRc1);
                #print(curResLoc);
                reSiteLocV.append(curResLoc);

                #Sequences location at position 84 (looking for RE)
                curPos84 = curLine[84:96];
                #print(curPos84);
                if curPos84 in eightyFourD :
                    eightyFourD[curPos84] = eightyFourD[curPos84] + 1;
                else :
                    eightyFourD[curPos84] = 1;

                #Find RE site by looping through possibilities, prioritize best match
                curResLoc = -1;
                for curRe in reSiteV[::-1] :
                    newPos = curLine.find(curRe);
                    if newPos > 0 :
                        curResLoc = newPos;
                reSiteLocLaxV.append(curResLoc);

                #Use the updated RE location to extract the barcode (revComp in sequence)
                if curResLoc > 0 :
                    curBcRcSeq = curLine[(curResLoc-16):curResLoc]
                    #print(curBcRcSeq);

                    if curBcRcSeq in barcodeRc2indexD :
                        #print(barcodeRc2indexD[curBcRcSeq])
                        barcodeCountV[barcodeRc2indexD[curBcRcSeq]] = barcodeCountV[barcodeRc2indexD[curBcRcSeq]] + 1
                        bcStatusV.append(1);
                    else :
                        bcStatusV.append(0);
                else :
                    bcStatusV.append(-1);

            curRow = curRow + 1;



        ### Print out various stats ####
        reSiteLocV.sort();
        print(reSiteLocV);

        reSiteLocLaxV.sort();
        print(reSiteLocLaxV);

        numNoRe = sum(1 for item in bcStatusV if item==(-1));
        numNoHit = sum(1 for item in bcStatusV if item==(0));
        numHit = sum(1 for item in bcStatusV if item==(1));

        print("RE Found Perc = " + str(numHit/len(bcStatusV)));

        #for curKey in eightyFourD :
        #    print(curKey + " - " + str(eightyFourD[curKey]))

        #print(",".join(map(str,barcodeCountV)));


    ############### Iterate through rows of a fastq file - Fast/efficient #############

    if loopFastB :

        now = datetime.now()

        current_time = now.strftime("%H:%M:%S")
        print("Starting Current Time =", current_time)

        reSiteLocLaxV = []; #Location of restriction enzyme site matches, relaxed
        eightyFourD = {}; #Dictionary counting sequences at position 84;

        bcStatusV = []; #Vector with whether a good barcode was found
                        #-1 no RE site; 0 = no bc hit; 1=good bc hit

        curRow = 0;

        for curLine in fileinput.input([curFastQFn]):

            if curRow % 4 == 1 : #Only take certain rows of the fastq
                #print(curLine);

                curResLoc = -1;

                curPos84 = curLine[84:96];

                if curPos84 in reSiteV :
                    curResLoc = 84;
                else :
                    curPos83 = curLine[83:95];
                    if curPos83 in reSiteV :
                        curResLoc = 83;
                    else :
                        curPos82 = curLine[82:94];
                        if curPos82 in reSiteV :
                            curResLoc = 82;
                        else :
                            curPos85 = curLine[85:97];
                            if curPos85 in reSiteV :
                                curResLoc = 85;
                            else : #Find RE site by looping through possibilities, prioritize best match
                                for curRe in reSiteV[::-1] :
                                    newPos = curLine.find(curRe);
                                    if newPos > 0 :
                                        curResLoc = newPos;

                reSiteLocLaxV.append(curResLoc);

                #Use the updated RE location to extract the barcode (revComp in sequence)
                if curResLoc > 0 :
                    curBcRcSeq = curLine[(curResLoc-16):curResLoc]
                    #print(curBcRcSeq);

                    if curBcRcSeq in barcodeRc2indexD :
                        #print(barcodeRc2indexD[curBcRcSeq])
                        barcodeCountV[barcodeRc2indexD[curBcRcSeq]] = barcodeCountV[barcodeRc2indexD[curBcRcSeq]] + 1
                        bcStatusV.append(1);
                    else :
                        bcStatusV.append(0);
                else :
                    bcStatusV.append(-1);

            curRow = curRow + 1;

        ### Print out various stats ####

        #reSiteLocLaxV.sort();
        #print(reSiteLocLaxV);

        numNoRe = sum(1 for item in bcStatusV if item==(-1));
        numNoHit = sum(1 for item in bcStatusV if item==(0));
        numHit = sum(1 for item in bcStatusV if item==(1));

        print("BC Found Perc = " + str(numHit/len(bcStatusV)));
        print("No Hit Perc = " + str(numNoHit/len(bcStatusV)));
        print("No RE Perc = " + str(numNoRe/len(bcStatusV)));

        statSummaryM.append([curExp,numHit,numNoHit,numNoRe])


        #for curKey in eightyFourD :
        #    print(curKey + " - " + str(eightyFourD[curKey]))

        #print(",".join(map(str,barcodeCountV)));

        now = datetime.now()

        current_time = now.strftime("%H:%M:%S")
        print("Starting Current Time =", current_time)


    if writeCountsB :

        curCountF = open(curCountFn,"w");
        curCountF.write(",".join(map(str,barcodeCountV)) + "\n");
        curCountF.close();



print(statSummaryM);

endProgram();



################################################
#   Conclusions Notes   #
################################################

#Some RE sites aren't found because they're not perfect
#Almost all RE sites are in positions 82,83,84,85; most in 84
    #Could be used to find imperfect sites that still have barcode
#Primer Site 2 is a complete mess,
    #often starts with TCA, not TCT, very short
#Real at position 84 TCTAGAGGTACC (98)
    #Common alternative at 84 TCTAGACGTACC (83), TCTAGAAGTACC (11), TCTAGATGTACC (4)
    #Also some shifted one base
#Generally 75-80% of sequences have a reverse barcode hit
#Numbers are better and there is less RE wobble for hit
