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
#Version 2.1.1 - Update to work in parallel per file
#              - Read in GZ fastq


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

#python3 arrayProc.2.1.1.py input/ input/MPRA_Array1_sequences 105DLiver 1-3-R-11_R1_001.test.fastq countsV2/ qcV2/

#python3 arrayProc.2.1.1.py /projects/MPRA/MPRA/MPRAi/MPRAi_v2_NovaSeq_Genewiz2/ input/MPRA_Array1_sequences 105DLiver 105DLiver_R1_001.fastq.gz countsV2/ qcV2/

################################################
#   Parameters   #
################################################


#Read in general information
#indexMapFn = str(sys.argv[1]); #The file that contains index information
seqDir = str(sys.argv[1]); #The directory that contains the fastqs
barcodeMapFn = str(sys.argv[2]); #The map between barcodes and enhancers

#Read in general information
sampName = str(sys.argv[3]); #105DLiver
fastqFn = str(sys.argv[4]); #105DLiver_R1_001.fastq.gz
countsDir = str(sys.argv[5]); #countsV2
qcDir = str(sys.argv[6]); #qcV2


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


################################################
#   Perform calcs for a specific fastq file   #
################################################
#For R1, reverse complement

outCountFn = countsDir + "counts." + sampName + ".2.1.1.csv";
outQcFn = qcDir + "qc." + sampName + ".2.1.1.csv";

fastqFnFull = seqDir + fastqFn;


############### Set parameters ###########################

reSiteRc1 = "TCTAGAGGTACC"; #Perfect match
reSiteRc2 = "TCTAGACGTACC"; #Second match
reSiteRc3 = "TCTAGAAGTACC"; #Third match
reSiteRc4 = "TCTAGATGTACC"; #Fourth match

reSiteV = [reSiteRc1,reSiteRc2,reSiteRc3,reSiteRc4];


#primer1Rc = "CGACGCTCTTCCGATCT";

#curExp sampName
#curCountFn outCountFn
#curFastQFn fastqFn

############### Iterate through rows of a fastq file - Detailed #############


barcodeCountV = [0] * barcodeArrayLen;

if loopDetailedB :

    reSiteLocV = []; #Location of restriction enzyme site matches
    reSiteLocLaxV = []; #Location of restriction enzyme site matches, relaxed

    eightyFourD = {}; #Dictionary counting sequences at position 84;

    bcStatusV = []; #Vector with whether a good barcode was found
                    #-1 no RE site; 0 = no bc hit; 1=good bc hit


    curRow = 0;


    for curLine in fileinput.input([fastqFnFull]):

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

    #for curLine in fileinput.input([fastqFnFull]):
    with gzip.open(fastqFnFull,'rt') as f:

        for curLine in f:
        #print('got line', line)

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


        #for curKey in eightyFourD :
        #    print(curKey + " - " + str(eightyFourD[curKey]))

        #print(",".join(map(str,barcodeCountV)));

        now = datetime.now()

        current_time = now.strftime("%H:%M:%S")
        print("Starting Current Time =", current_time)


if writeCountsB :

    curCountF = open(outCountFn,"w");
    curCountF.write(",".join(map(str,barcodeCountV)) + "\n");
    curCountF.close();

    curQcF = open(outQcFn,"w");
    curQcF.write(",".join(map(str,[sampName,numHit,numNoHit,numNoRe])) + "\n");
    curQcF.close();


endProgram();

if False :
    curCountF = open("sumStats.mpra.csv","w");
    for curLineV in summaryStatsM :
        curCountF.write(",".join(map(str,curLineV)) + "\n");
    curCountF.close();

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
