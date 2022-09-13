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


curFastQFn = dir2InfoD["STR_1_R"][7];
#curFastQFn = "input/1-3-R-11_R1_001.test.fastq";
#curFastQFn = "input/1-3-D-19_R1_001.test.fastq";
#curCountFn = "counts/counts.STR_1_D_test.1.1.csv";
curCountFn = "counts/counts." + dir2InfoD["STR_1_R"][0] + ".1.1.csv";
print(curFastQFn);

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


################################################
#   Loop through the fast q file performing calcs   #
################################################

#Load the matrix of zoonomia species annotations
halInfoM = [];

genomeSourceD = {};
genomeFastaD = {};
referenceGenomeFn = "";
spec2indD = {};

namingConventionMasterD = {"GenBank_to_RefSeq":"GenbankToRefseqName","GenBank_to_UCSC":"GenbankToUCSCName","GenBank_to_UCSC":"GenbankToUCSCName","Sequence-Name_to_GenBank":"SequenceNameToGenbankName"};

curRow = 0;

for curLine in fileinput.input([xspecGuideFn]): #Read in list of things for gene distance
    curLineP = curLine.rstrip().split(",");
    if(curRow > 0 and curLineP[0] != '') :
        #halInfoM.append(curLineP);

        curSpec = curLineP[1];

        halInfoM.append(curLineP)
        genomeSourceD[curSpec] = curLineP[9];

        if(curSpec == referenceSpecies) :
            referenceGenomeFn = curLineP[11];

        genomeFastaD[curSpec] = curLineP[11];
        spec2indD[curSpec] = curRow - 1;


    curRow = curRow + 1;

#print(halInfoM);
print(referenceGenomeFn);


#Load the matrix of pairwise distances
specDistD = {}; #Dictionary object the stores the pairwise distances
curRow = 0;
specDistOrderP = []; #Order of species in matrix, not stars at index 1

for curLine in fileinput.input([specDistFn]): #Read in list of things for gene distance
    curLineP = curLine.rstrip().split(",");
    if curRow == 0 :
        specDistOrderP = curLineP;

    if(curRow > 0 and curLineP[1] != '') :
        curSpec1 = curLineP[0];

        for curInd in range(len(curLineP)-1) :
            curSpec2 = specDistOrderP[curInd + 1];
            curVal = curLineP[curInd+1];
            if(curVal != '') :
                curKey = curSpec1 + "-" + curSpec2;
                specDistD[curKey] = float(curVal);
                #print(curKey + " : " + curVal);
                curKey = curSpec2 + "-" + curSpec1;
                specDistD[curKey] = float(curVal);

    curRow = curRow + 1;


################################################
#   Extract reference coding sequence location from GTF   #
################################################

############# Read in refSeq chromosomes coding ###################

print("Starting to read in GTF File\n");

refSeq2ChrD = {}; #If needed, a way to change the reference species chromosome names


############# Read in refSeq chromosomes coding ###################

#proteinSrch = re.compile("protein_id")

protNoVersion = useRefProtId.split(".")[0];
#print(protNoVersion);
curProtBedFn = "protBed/protBed.%s.%s.bed" % (useRefProtId,fileSuffixStr);
refProtBed = []; #Bed of the reference protein

refProtGeneEns = "";

if gtf2BedB :

    curProtBedF = open(curProtBedFn,"w");
    with gzip.open(proteinGtfFn,'rt') as fInput:

        for curLine in fInput : #Read in list of things for gene distance
            if not curLine.startswith("#") :
                curLineP = curLine.rstrip().split("\t");

                if curLineP[2] == "CDS" :
                    #print(curLineP[8]);
                    #curRe = re.search('protein_id \"(w+)\"',curLineP[8]);
                    #print(curLineP);
                    #stophere();
                    #curRe = re.search('protein_id ["](\w+)[.][(\w+)]["]',curLineP[8]);
                    curRe = re.search('protein_id ["](\w+)["]',curLineP[8]);
                    if curRe != None :
                        curProt = curRe.group(1);
                        #curProt2 = re.findall(r'"([^"]*)"', curProt)[0];
                        curProt2 = curProt;
                        #print(curProt);
                        if curProt2 == useRefProtId : #Only if it matches the right protein
                            #print(curLineP[8]);

                            #Get the exon number
                            curRe2 = re.search('exon_number ["](\w+)["]',curLineP[8]);
                            curExon = curRe2.group(0);
                            curExon2 = re.findall(r'"([^"]*)"', curExon)[0];
                            curExon3 = "e_" + curExon2;
                            curId = curProt2 + "." + curExon3;

                            curScore = str(int(curLineP[4]) - int(curLineP[3]) + 1);

                            curChr = curLineP[0];

                            curRe3 = re.search('gene_id ["](\w+)["]',curLineP[8]);
                            refProtGeneEns = curRe3.group(1);

                            #Adjust the nucleotide position by 1 because of indexing
                            if curChr in refSeq2ChrD :
                                bedOutP = [refSeq2ChrD[curChr],str(int(curLineP[3])-1),curLineP[4],curId,curScore,curLineP[6]];

                            else :
                                bedOutP = ["chr" + curChr ,str(int(curLineP[3])-1),curLineP[4],curId,curScore,curLineP[6]];


                            print(bedOutP);
                            curProtBedF.write("\t".join(bedOutP) + "\n");
                            refProtBed.append(bedOutP);


    curProtBedF.close();


print (refProtGeneEns);

################################################
#  Get human protein sequence #
################################################

curRow = 0;
referenceAaLen = 0;
referenceStrand = "";
exon2indD = {}; #Map of the exon name to the index in the table
exonStatsColP = ["chrVal","hgStart","hgStop","exonId","len_hg","strand"]; #Columns for the bed files\

exonStatsM = []; #Statistics for the reference exon
exonNamesP = []; #Names for the reference exon

anchorProteinD = {};
anchorNucsD = {};


### Build the BED file information for the reference species ####
for curLine in fileinput.input([curProtBedFn]): #Read in stats for human protein bed file
    curLineP = curLine.rstrip().split("\t");

    curExonName = curLineP[3];
    referenceStrand = curLineP[5];

    exonStatsM.append(curLineP);
    exonNamesP.append(curLineP[3]);
    exon2indD[curLineP[3]] = curRow;

    referenceAaLen = referenceAaLen + int(curLineP[2]) - int(curLineP[1]);
    curRow = curRow + 1;

referenceAaLen = referenceAaLen/3;
print(referenceAaLen);

### Get the sequence for the reference species ####
getSeqStr = "bedtools getfasta -fi %s -bed %s" % (referenceGenomeFn,curProtBedFn);
curSpecCdsSeq = Seq("");
process = subprocess.Popen(getSeqStr, shell=True, stdout=subprocess.PIPE);

### Pull out the sequences from the BEDs ###
for curIoLine in io.TextIOWrapper(process.stdout, encoding="utf-8"):
    #print(line.rstrip());
    if(not curIoLine.startswith(">")) :
        newSpecCdsSeq = Seq(curIoLine.rstrip(),ambiguous_dna);
        if(referenceStrand == "-") :
            newSpecCdsSeq=newSpecCdsSeq.reverse_complement();
        curSpecCdsSeq = curSpecCdsSeq + newSpecCdsSeq;
    else :
        curIoLine2 = curIoLine.rstrip();


referenceProtSeq = curSpecCdsSeq.translate();

anchorProteinD[referenceSpecies] = referenceProtSeq;
anchorNucsD[referenceSpecies] = curSpecCdsSeq;

print(referenceProtSeq);


################################################
#   Mapping to anchor species   #
################################################

curRow = 0;

annotFd = "startingAnnot/ensembl/annotationGTFs/";
orthoFd = "startingAnnot/ensembl/orthology/";
anchorFd = "anchorBeds/";

anchorSpecM = []; #Species Name, Annotation File, Orthology File

print("Mapping to anchor species\n");

#New ENSEMBLID, Orthology Type, Number of protein matches, best prot id, length overlap from hal, length missing from hal
anchorStatsM = []; #Basic information about the mapping of each anchor species

#Which anchor proteins to use based on which mapped across genomes
useAnchorP = [referenceSpecies];

#Which bed file to use an anchor species
useAnchorBedD = {};
useAnchorBedD[referenceSpecies] = curProtBedFn;


#Loop through the anchor species and get the best proteins
for curLine in fileinput.input([anchorSpecFn]): #Read in list of things for gene distance
    curLineP = curLine.rstrip().split(",");

    if(curRow > 0) :
        curSpecName = curLineP[0];
        curAnnotFn = annotFd + curLineP[1];
        curOrthoFn = orthoFd + curLineP[2];
        anchorSpecM.append(curLineP);
    curRow = curRow + 1;

#print(anchorSpecM);

allAnchorP = [];

if anchorOverB :

    for curLineP in anchorSpecM :

        curSpecName = curLineP[0];
        curAnnotFn = annotFd + curLineP[1];
        curOrthoFn = orthoFd + curLineP[2];

        allAnchorP.append(curSpecName);

        curSpecLineP = halInfoM[spec2indD[curSpecName]]; #Load the information for this species

        #New ENSEMBLID, Orthology Type, Number of protein matches, best prot id, length overlap from hal, length missing from hal
        curAnchorStatP = ["NA"] * 7;



        curEnsChrNameB = False; #Whether the ENSEMBL names need to be converted
        curEnsChrNameFn = "startingAnnot/ensembl/ens2ncbi-%s.csv" % (curSpecName);
        curEnsChrNameD = {};

        if path.exists(curEnsChrNameFn) :
            curEnsChrNameB = True;

            for curLine2 in fileinput.input([curEnsChrNameFn]): #Read in conversion of chromosome names
                curLine2P = curLine2.rstrip().split(",");

                curEnsChrName = curLine2P[0];
                curEnsChrName = curEnsChrName.replace("Chromosome ","");

                curEnsNewName = "";

                if(curSpecLineP[8] == "UCSC") :
                    curEnsNewName = "chr" + curEnsChrName;
                if(curSpecLineP[8] == "GenBank") :
                    curEnsNewName = curLine2P[1];
                if(curSpecLineP[8] == "RefSeq") :
                    curEnsNewName = curLine2P[2];


                curEnsChrNameD[curEnsChrName] = curEnsNewName

        #print(curSpecLineP);

        print(curLineP);

        #Run HAL liftover to this species to be able to benchmark other proteins
        #hloOutbedFn = "%s/protBed.%s.%s.2.1.1.bed" % (outBedDir,curSpecName,useRefProtId);
        hloOutbedFn = "%s/protBed.hal.%s.%s.%s.bed" % (outBedDir,curSpecName,useRefProtId,fileSuffixStr);
        callHLOStr = "halLiftover --bedType 6 %s Homo_sapiens %s %s %s" % (halFn,curProtBedFn,curSpecName,hloOutbedFn);
        #print(callHLOStr);
        print(curSpecName);
        subprocess.Popen(callHLOStr, shell=True).wait();
        #stophere();

        #Read in HAL liftover results in curXSpecHalBedP
        curXSpecHalBedP = [];
        for curLine2 in fileinput.input([hloOutbedFn]): #Read in results of HAL Liftover
            curLine2P = curLine2.rstrip().split("\t");
            curXSpecHalBedP.append(curLine2P);
        #print(curXSpecHalBedP);

        #Get the orthologous ENSEMBL ID
        curRow2 = 0;
        curXSpecGeneEns = "NA";
        for curLine2 in fileinput.input([curOrthoFn]): #Read in table of orthologous sequences
            curLine2P = curLine2.rstrip().split("\t");
            if(curLine2P[0] == refProtGeneEns) :
                if(curLine2P[3] == "ortholog_one2one") : #This is where criteria are set based on orthology for using a gene
                    curXSpecGeneEns = curLine2P[2];
                else :
                    curAnchorStatP[1] = curLine2P[3]
            curRow2 = curRow2 + 1;

        curAnchorStatP[0] = curXSpecGeneEns;

        #print(curXSpecGeneEns);

        if(curXSpecGeneEns != "NA") :
            curAnchorStatP[1] = "ortholog_one2one";
            #Scan through the GTP file for every protein corresponding to that gene
            curXSpecBedD = {}; #A dictionary indexed by the proteinID where every element is a different bed file
            curXSpecProtP = []; #The different protein IDs for that gene
            with gzip.open(curAnnotFn,'rt') as fInput:

                for curLine2 in fInput : #Read in list of things for gene distance
                    if not curLine2.startswith("#") :
                        curLine2P = curLine2.rstrip().split("\t");

                        if curLine2P[2] == "CDS" :
                            curRe = re.search('gene_id ["](\w+)["]',curLine2P[8]);
                            if curRe != None :
                                curGeneId = curRe.group(1);
                                if curGeneId == curXSpecGeneEns : #Only if it matches the right geneID
                                    #print(curLineP[8]);

                                    curRe3 = re.search('protein_id ["](\w+)["]',curLine2P[8]);
                                    curProt2 = curRe3.group(1);

                                    #Get the exon number
                                    curRe2 = re.search('exon_number ["](\w+)["]',curLine2P[8]);
                                    curExon = curRe2.group(1);
                                    curExon3 = "e_" + curExon;
                                    curId = curProt2 + "." + curExon3;

                                    curScore = str(int(curLine2P[4]) - int(curLine2P[3]) + 1); #Score is sequence length

                                    curChr = curLine2P[0];

                                    if(curEnsChrNameB and curChr in curEnsChrNameD) : #Fix the chromosomes
                                        curChr2 = curEnsChrNameD[curChr];
                                    else :
                                        if (genomeSourceD[curSpecName] == "UCSC") :
                                            curChr2 = "chr" + curChr;
                                        else :
                                            curChr2 = curChr

                                    bedOutP = [curChr2 ,str(int(curLine2P[3])-1),str(int(curLine2P[4])),curId,curScore,curLine2P[6]];

                                    if curProt2 in curXSpecBedD :
                                        curXSpecBedD[curProt2].append(bedOutP);
                                    else :
                                        curXSpecBedD[curProt2] = [bedOutP];
                                        curXSpecProtP.append(curProt2)


            #print(curXSpecBedD[curProt2]);
            #print(curXSpecProtP);

            curAnchorStatP[2] = len(curXSpecProtP);

            if( len(curXSpecProtP) > 0) : #Only scan through different proteins if there are any found in the anchor species

                #Loop through each unique protein
                xSpecProtLenLapP = [0] * len(curXSpecProtP);
                xSpecProtLenMissingP = [0] * len(curXSpecProtP);

                curProtInd = 0;
                for curXSpecProt in curXSpecProtP :

                    curXSpecEnsBedP = curXSpecBedD[curXSpecProt];
                    #curLapResP =  calcBedLap(curXSpecHalBedP,curXSpecEnsBedP); #Get the overlap of the two beds [overlap,missing]
                    curLapResP =  calcBedLap(curXSpecEnsBedP,curXSpecHalBedP); #Get the overlap of the two beds [overlap,missing]
                    xSpecProtLenLapP[curProtInd] = curLapResP[0];
                    xSpecProtLenMissingP[curProtInd] = curLapResP[1];

                    curProtInd = curProtInd + 1;


                #### Choose the best of the ENSEMBL proteins as the referece #####
                curBestProtScore = -1000000;
                curBestProInd = -1;

                for curProtInd in range(len(curXSpecProtP)) :
                    curProtScore = xSpecProtLenLapP[curProtInd] - xSpecProtLenMissingP[curProtInd]; # Score each base that is different
                    if curProtScore > curBestProtScore :
                        curBestProInd = curProtInd;
                        curBestProtScore = curProtScore;

                curBestXspecProt = curXSpecProtP[curBestProInd]; #The best ENSEMBL protein
                curBestXspecBedP = curXSpecBedD[curBestXspecProt]; #The BED of the best ENSEMBL protein

                #print(curBestProInd);
                print(curBestXspecProt);

                curAnchorStatP[3] = curBestXspecProt;
                curAnchorStatP[4] = xSpecProtLenLapP[curBestProInd];
                curAnchorStatP[5] = xSpecProtLenMissingP[curBestProInd];

                ensOutbedFn = "%s/protBed.%s.%s.halLiftChr.%s.bed" % (anchorFd,curSpecName,useRefProtId,fileSuffixStr);

                curProtBedF = open(ensOutbedFn,"w");
                for curProtBedP in curBestXspecBedP :
                    curProtBedF.write("\t".join(curProtBedP) + "\n");

                    curStrand = curProtBedP[5];

                curProtBedF.close();

                #### Create the naming convension dictionary ####
                skipNameB = False;

                if(curSpecLineP[10] == "N") : #If name convension doesn't match

                    namingConvKey = curSpecLineP[8] + "_to_" + curSpecLineP[9];
                    namingConvValue = namingConventionMasterD[namingConvKey];

                    #Read in a dictionary that will transfer the naming conventions
                    curChrNameFn = "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/200MammalsFastas/%s_%s.txt" % (namingConvValue,curSpecName);

                    if not path.exists(curChrNameFn) : #adjust naming convension name if you can't find it
                        curChrNameFn = "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/200MammalsFastas/%s_%s" % (namingConvValue,curSpecName);

                    if path.exists(curChrNameFn) :

                        curChrNameD = {};
                        for curLine in fileinput.input([curChrNameFn]): #Read in list of things for gene distance
                            curLineP = curLine.rstrip().split("\t");
                            curChrNameD[curLineP[0]] = curLineP[1];
                    else :
                        skipNameB = True;

                ### Convert the BED file to read the from the FASTA ###
                ensOutbedFastaFn = "%s/protBed.%s.%s.fastaChr.%s.bed" % (anchorFd,curSpecName,useRefProtId,fileSuffixStr);
                curProtBedF = open(ensOutbedFastaFn,"w");
                for curProtBedP in curBestXspecBedP :

                    #Update the naming convention
                    if(curSpecLineP[10] == "N" and not skipNameB) :
                        curProtBedP[0] = curChrNameD[curProtBedP[0]]; # Change the chromosome naming convension

                    curProtBedF.write("\t".join(curProtBedP) + "\n");

                curProtBedF.close();

                ### Get the sequence for the anchor species ####
                curSpecFasta = curSpecLineP[11];
                if( not path.exists(curSpecFasta) ) :
                    print("FASTA NOT FOUND - " + curSpecFasta );
                else :

                    getSeqStr = "bedtools getfasta -fi %s -bed %s" % (curSpecFasta,ensOutbedFastaFn);
                    curSpecCdsSeq = Seq("");
                    process = subprocess.Popen(getSeqStr, shell=True, stdout=subprocess.PIPE);

                    ### Pull out the sequences from the BEDs ###
                    for curIoLine in io.TextIOWrapper(process.stdout, encoding="utf-8"):
                        #print(line.rstrip());
                        if(not curIoLine.startswith(">")) :
                            newSpecCdsSeq = Seq(curIoLine.rstrip(),ambiguous_dna);
                            if(curStrand == "-") :
                                newSpecCdsSeq=newSpecCdsSeq.reverse_complement();
                            curSpecCdsSeq = curSpecCdsSeq + newSpecCdsSeq;
                        else :
                            curIoLine2 = curIoLine.rstrip();


                    #print(curSpecCdsSeq);
                    #print();
                    curSpecCdsProt = curSpecCdsSeq.translate();
                    print(curSpecCdsProt);

                    anchorProteinD[curSpecName] = curSpecCdsProt;
                    anchorNucsD[curSpecName] = curSpecCdsSeq;
                    curAnchorStatP[6] = len(curSpecCdsProt);
                    useAnchorP.append(curSpecName);
                    useAnchorBedD[curSpecName] = ensOutbedFn;
        anchorStatsM.append(curAnchorStatP);



#print(anchorStatsM);
##New ENSEMBLID, Orthology Type, Number of protein matches, best prot id, length overlap from hal, length missing from hal
anchorNumMFn = "output/exonStats/anchorSpecSummary.%s.%s.txt" % (useRefProtId,fileSuffixStr);
anchorInfoP = ["ensGene","orthoType","numProtMatches","ensProt","halOverlapBp","halMissingBp","aaLen"];

print(anchorStatsM);
print(useAnchorP[1:]);

writeTable(allAnchorP,anchorInfoP,anchorStatsM,anchorNumMFn);



################################################
#   Based on which anchor proteins map, calculate the best match   #
#   Also gather exon stats from anchor species
################################################

exonAnchorStatsD = {}; #Exon stats, each species is an entry, each element is a  matrix
exonAnchorNamesD = {}; #Exon names, each species is an entry, each element is a vector
exon2indDD = {}; #A dictionary of dictionaries, mapping exon names to indices in the matrix

for curSpecLineP in halInfoM :

    curUseAnchor = "anchor";
    curSpecName = curSpecLineP[1];

    curAnchorDist = 1.1;

    if(not curSpecName in  useAnchorP) : #Don't search for anchor is species is an anchor
        for curAnchorSpec in useAnchorP : #Loop through anchor species and get the lowest distance
            subAnchorDist = specDistD[curSpecName + "-" + curAnchorSpec];
            if(subAnchorDist < curAnchorDist) :
                curAnchorDist = subAnchorDist;
                curUseAnchor = curAnchorSpec;

    else : #Gather exon stats if the species is an anchor

        curAnchorDist = 0;

        if(curSpecName == referenceSpecies) :
            exonAnchorStatsD[curSpecName] = exonStatsM;
            exonAnchorNamesD[curSpecName] = exonNamesP;
            exon2indDD[curSpecName] = exon2indD;

        else :

            curExonStatsM = []; #Statistics for the reference exon
            curExonNamesP = []; #Names for the reference exon
            curExon2indD = {}; #Maps exons to indices

            #curAchorBedFn = "%s/protBed.%s.%s.halLiftChr.2.1.1.bed" % (anchorFd,curSpecName,useRefProtId);
            #curAchorBedFn = "%s/protBed.%s.%s.halLiftChr.2.1.1.bed" % (anchorFd,curSpecName,useRefProtId);
            #ensOutbedFastaFn = "%s/protBed.%s.%s.fastaChr.2.1.1.bed" % (anchorFd,curSpecName,useRefProtId);
            curAchorBedFn = useAnchorBedD[curSpecName];

            curRow = 0;
            ### Build the BED file information for the reference species ####
            for curLine in fileinput.input([curAchorBedFn]): #Read in stats for human protein bed file
                curLineP = curLine.rstrip().split("\t");

                curExonName = curLineP[3];
                curReferenceStrand = curLineP[5];

                curExonStatsM.append(curLineP);
                curExonNamesP.append(curLineP[3]);
                curExon2indD[curLineP[3]] = curRow;

                #referenceAaLen = referenceAaLen + int(curLineP[2]) - int(curLineP[1]);
                curRow = curRow + 1;

            exonAnchorStatsD[curSpecName] = curExonStatsM;
            exonAnchorNamesD[curSpecName] = curExonNamesP;
            exon2indDD[curSpecName] = curExon2indD;

    curSpecLineP.append(curUseAnchor);
    curSpecLineP.append(curAnchorDist);

    #print(curSpecLineP);

print("Print Index Info");

for curKey in exon2indDD :
    print(curKey);
    for curKey2 in exon2indDD[curKey] :
        print("\t" + curKey2);


#stophere();


################################################
#   Code to run halLiftover on bed file for a protein   #
################################################

#print(halInfoM);

if halOverB :

    print("Running HAL Liftover\n");

    for curSpecLineP in halInfoM :

        curSpecName = curSpecLineP[1];

        if( not curSpecName in useAnchorP) :

            curAnchorSpec = curSpecLineP[15];
            curAchorBedFn = useAnchorBedD[curAnchorSpec];

            hloOutbedFn = "%s/protBed.%s.%s.%s.bed" % (outBedDir,curSpecName,useRefProtId,fileSuffixStr);

            callHLOStr = "halLiftover --bedType 6 %s %s %s %s %s" % (halFn,curAnchorSpec,curAchorBedFn,curSpecName,hloOutbedFn);
            #print(callHLOStr);

            print(curSpecName + " from " + curAnchorSpec);

            subprocess.Popen(callHLOStr, shell=True).wait();

            #stophere();


################################################
#   Post-processing on halLiftover results   #
################################################

print("Processing HAL Liftover stats\n"); #Don't bother doing for anchor species

halSpeciesNamesP = []; #The names of the species that have summary information and have been mapped over
spec2HalIndD = {}; #The indices of the species in the hal summary
halSpecBedM = []; #An array of matrices where each element is BED file for that species

#Items to handle decision to merge exons
mergedExonsP = []; #Indexed by species_exonName, list of exons for printing
halSpecBedMergeD = {}; #Indexed by species_exonName, dictionary where every element is merged line of bed file
halSpecBedSepD = {}; #Indexed by species_exonName, dictionary where every element is a matrix of separated elements


#Different matrices for each anchor species
exonNucNumDD = {}; #Dictionary of the number of nucleotides per exon
exonFragNumDD = {}; #Dictionary of number of fragments per exon
exonScaffoldsNumDD = {}; #Dictionary of number of scaffolds per exon

#Initialize the dictionary
for curAnchor in useAnchorP :
    exonNucNumDD[curAnchor] = {};
    exonFragNumDD[curAnchor] = {};
    exonScaffoldsNumDD[curAnchor] = {};

#fragNumP = []; #Vector of number of fragments per gene
#scafNumP = []; #Vector of scaffolds per gene

#curSpecName,curFragNum,curScafNum,curAANum,curMissingExon,curMultFragExon,curMultScafExon;
specSummaryM = []; #Summary information for each species


curRow = 0;
for curSpecLineP in halInfoM :
#for curSpecLineP in halInfoM[1:10] :
    curSpecName = curSpecLineP[1];
    print(curSpecName);

    hloOutbedFn = "%s/protBed.%s.%s.%s.bed" % (outBedDir,curSpecName,useRefProtId,fileSuffixStr); #Bed file of mapping for current species


    if(path.exists(hloOutbedFn) and not curSpecName in useAnchorP and not curSpecLineP[8] == '') :

        curAnchorSpec = curSpecLineP[15];


        halSpeciesNamesP.append(curSpecName);
        exonStatsColP.append(curSpecName);


        #### Gather statistics on the exons #######
        curAnchorExonStatsM = exonAnchorStatsD[curAnchorSpec];

        curExonNucNumP = [0] * len(curAnchorExonStatsM); #Vector, number of nucleotides per exon
        curExonFragNumP = [0] * len(curAnchorExonStatsM); #Vector, number of fragments per exon
        curExonScafListP = [[] for i in range(len(curAnchorExonStatsM))]; #Matrix of scaffolds, one per exon
        curExonScafNumP = [0] * len(curAnchorExonStatsM); #Vector, number of scaffolds per exon

        curAllScafList = []; #Matrix of scaffolds, one per gene

        curAANum = 0;

        for curLine in fileinput.input([hloOutbedFn]) : #Read in bed file for current protein
            curLineP = curLine.rstrip().split("\t");
            curExon = curLineP[3];
            curExonIndex = exon2indDD[curAnchorSpec][curExon];

            curExonNucNumP[curExonIndex] = curExonNucNumP[curExonIndex] + int(curLineP[2]) - int(curLineP[1]);
            curExonFragNumP[curExonIndex] = curExonFragNumP[curExonIndex] + 1;
            curExonScafListP[curExonIndex].append(curLineP[0]);
            curAllScafList.append(curLineP[0]);
            #curAANum = curAANum + int(curLineP[2]) - int(curLineP[1]) + 1;
            curAANum = curAANum + (int(curLineP[2]) - int(curLineP[1]));

        for curExonNum in range(len(curAnchorExonStatsM)) :
            curExonScafNumP[curExonNum] = len(set(curExonScafListP[curExonNum]));

        #exonNucNumM.append(curExonNucNumP);
        #exonFragNumM.append(curExonFragNumP);
        #exonScaffoldsNumM.append(curExonScafNumP);
        exonNucNumDD[curAnchorSpec][curSpecName] = curExonNucNumP;
        exonFragNumDD[curAnchorSpec][curSpecName] = curExonFragNumP;
        exonScaffoldsNumDD[curAnchorSpec][curSpecName] = curExonScafNumP;

        curFragNum = sum(curExonFragNumP);
        curScafNum = len(set(curAllScafList))
        curAANum = curAANum/3;

        curMissingExon = sum(frag == 0 for frag in curExonFragNumP);
        curMultFragExon = sum(frag > 1 for frag in curExonFragNumP);
        curMultScafExon = sum(scaf > 1 for scaf in curExonScafNumP);

        ##### Produce a new bed file #######
        curSpecFasta = curSpecLineP[11];

        #Create the naming convension dictionary
        skipNameB = False;
        if(curSpecLineP[10] == "N") :

            namingConvKey = curSpecLineP[8] + "_to_" + curSpecLineP[9];
            namingConvValue = namingConventionMasterD[namingConvKey];
            #print(namingConvValue);

            #Read in a dictionary that will transfer the naming conventions
            curChrNameFn = "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/200MammalsFastas/%s_%s.txt" % (namingConvValue,curSpecName);

            if not path.exists(curChrNameFn) :
                curChrNameFn = "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/200MammalsFastas/%s_%s" % (namingConvValue,curSpecName);

            if path.exists(curChrNameFn) :

                curChrNameD = {};
                for curLine in fileinput.input([curChrNameFn]): #Read in list of things for gene distance
                    curLineP = curLine.rstrip().split("\t");
                    curChrNameD[curLineP[0]] = curLineP[1];
            else :
                skipNameB = True;

        #Transfer naming conventions and bed coordinate shifts to create a new bed file
        hloOutbed2Fn = "%s/protBedHalName.%s.%s.%s.bed" % (outBedDir,curSpecName,useRefProtId,fileSuffixStr);
        hloOutbed2F = open(hloOutbed2Fn,"w");

        curStrand = "NA";
        #prevExonName = ""; #Only subtract from coordinate if there is a new exon
        newExonCoordM = [["NA","NA","NA","NA","NA","NA"] for i in range(len(exonAnchorNamesD[curAnchorSpec]))]; #Matrix (one row per exon)
        newExonCoordD = {}; #Dictorionary that links the FASTA header back to the line

        curExonMapD = {}; #Dictionry indexed by exonName. Every element is matrix of all BED mathces.

        for curLine in fileinput.input([hloOutbedFn]): # Read in the Output of HAL liftover
            curLineP = curLine.rstrip().split("\t");
            curStrand = curLineP[5];
            curExonName = curLineP[3];
            curExonIndex = exon2indDD[curAnchorSpec][curExonName];

            #print(curLineP);

            #Update the naming convention
            if(curSpecLineP[10] == "N" and not skipNameB) :
                curLineP[0] = curChrNameD[curLineP[0]]; # Change the chromosome naming convension

            if(not curExonName in curExonMapD) :
                curExonMapD[curExonName] = [];

            curExonMapD[curExonName].append(curLineP);
            newExonInfoP = curLineP.copy();

            #Create new record if exon hasn't been observed
            if(newExonCoordM[curExonIndex][0] == "NA") :
                newExonCoordM[curExonIndex] = newExonInfoP;
            else :
                if(curLineP[1] < newExonCoordM[curExonIndex][1]) :
                    newExonCoordM[curExonIndex][1] = curLineP[1];
                if(curLineP[2] > newExonCoordM[curExonIndex][2]) :
                    newExonCoordM[curExonIndex][2] = curLineP[2];

        #Print out the final exon results
        myCurExonInd = 0;

        #print(curExonMapD["ENSCHIP00000016062.e_3"]);

        #print(newExonCoordM);

        curNumExonSep = 0;
        curNumExonMerge = 0;
        curNumExonSimp = 0;
        curAANumNew = 0;

        curFinalBedM = [];

        for curExonP in newExonCoordM :
            if(not curExonP[0] == "NA") :
                #print("exon name - " + curExonP[3]);
                curExonCoordId = ">%s:%s-%s" % (curExonP[0],curExonP[1],curExonP[2]);
                newExonCoordD[curExonCoordId] = curExonP;
                anchorNucNum = int(curExonP[4]); #Nucleotide length of the anchor
                origNucNum = curExonNucNumP[myCurExonInd]; #Nucletoide length of the separated
                newNucNum = int(curExonP[2]) - int(curExonP[1]); #Nucleotide length of the merged

                if(len(curExonMapD[curExonP[3]]) > 1) :

                    #Decide whether to collapse the parts of the exon into a larger piece
                    curCollapseExonB = False; #Default is not to collapse


                    #print("\tAnchNuc " + str(anchorNucNum) + "\tOrigNuc " + str(origNucNum) + "\tmergeNuc " + str(newNucNum));

                    if(not curExonScafNumP[myCurExonInd] > 1) : #Don't collapse exons for multiple scaffolds
                        #print("\tonly one scaffold");
                        if(newNucNum < (3*anchorNucNum)) : #Don't collapse if you more than triple in size from anchor by collapsing
                            #print("\tnot too long");
                            if( ((anchorNucNum % 3) == (origNucNum % 3)) and not ((anchorNucNum % 3) == (newNucNum % 3)) ) : #Collpase unless you lose inframe exons
                                #print("\tsameModulo");
                                curCollapseExonB = False;
                            else :
                                curCollapseExonB = True;

                    if(curCollapseExonB) :
                        curAANumNew = curAANumNew + newNucNum;
                        hloOutbed2F.write("\t".join(curExonP) + "\n");
                        curNumExonMerge = curNumExonMerge + 1;
                        curFinalBedM.append(curExonP); #Print merged coordinates

                        #Save merged exon information to go back if merged doesn't work well
                        curSpecExonName = curSpecName + "_" + curExonP[3];
                        mergedExonsP.append(curSpecExonName);
                        halSpecBedMergeD[curSpecExonName] = curExonP; #Store merged coordinates
                        halSpecBedSepD[curSpecExonName] = curExonMapD[curExonP[3]]; #Store separated coordates

                        #print("Exon");
                        #print(curExonP);
                    else :

                        curAANumNew = curAANumNew + origNucNum;

                        for curSubExonP in curExonMapD[curExonP[3]] :
                            #print("SubExon");
                            #print(curSubExonP);
                            hloOutbed2F.write("\t".join(curSubExonP) + "\n");
                            curNumExonSep = curNumExonSep + 1;
                            curFinalBedM.append(curSubExonP);

                else :
                    curAANumNew = curAANumNew + origNucNum;
                    hloOutbed2F.write("\t".join(curExonP) + "\n");
                    curNumExonSimp = curNumExonSimp + 1;
                    curFinalBedM.append(curExonP);



            #prevExonName = curExonName;
            myCurExonInd = myCurExonInd + 1;

        hloOutbed2F.close();

        curAANumNew = curAANumNew/3

        specSummaryP = [curSpecName,curFragNum,curScafNum,curAANum,curMissingExon,curMultFragExon,curMultScafExon,curNumExonSep,curNumExonMerge,curNumExonSimp,curAANumNew];
        specSummaryM.append(specSummaryP);

        spec2HalIndD[curSpecName] = curRow;
        halSpecBedM.append(curFinalBedM)
        curRow = curRow + 1;

        #stophere();

#rowNamesP, colnamesP, matrixM, fileNameStr
#print(exonNamesP);
#print(speciesNamesP);
#print(len(exonNucCountM));
#print(len(exonNucCountM[1]));

#exonNucNumDD[curAnchrSpec][curSpecName] = curExonNucNumP;
#exonFragNumDD[curAnchrSpec][curSpecName] = curExonFragNumP;
#exonScaffoldsNumDD[curAnchrSpec][curSpecName] = curExonScafNumP;

for curAnchorSpec in useAnchorP :

    curPrintSpecP = [];
    exonNucNumM = [];
    exonFragNumM = [];
    exonScaffoldsNumM = [];

    for specKey in exonNucNumDD[curAnchorSpec] :
        curPrintSpecP.append(specKey);
        exonNucNumM.append(exonNucNumDD[curAnchorSpec][specKey]);
        exonFragNumM.append(exonFragNumDD[curAnchorSpec][specKey]);
        exonScaffoldsNumM.append(exonScaffoldsNumDD[curAnchorSpec][specKey]);


    exonNucNumMFn = "output/exonStats/exonNucNumM.anchor_%s.%s.%s.txt" % (curAnchorSpec,useRefProtId,fileSuffixStr);
    writeTable(curPrintSpecP,exonNamesP,exonNucNumM,exonNucNumMFn);

    exonFragNumMFn = "output/exonStats/exonFragNumM.anchor_%s.%s.%s.txt" % (curAnchorSpec,useRefProtId,fileSuffixStr);
    writeTable(curPrintSpecP,exonNamesP,exonFragNumM,exonFragNumMFn);

    exonScaffoldsNumMFn = "output/exonStats/exonScaffoldsNumM.anchor_%s.%s.%s.txt" % (curAnchorSpec,useRefProtId,fileSuffixStr);
    writeTable(curPrintSpecP,exonNamesP,exonScaffoldsNumM,exonScaffoldsNumMFn);

specSummaryMFn = "output/exonStats/specSummaryM.%s.%s.txt" % (useRefProtId,fileSuffixStr);
specSummaryColP = ["curSpecName","curFragNum","curScafNum","curAANum","curMissingExon","curMultFragExon","curMultScafExon","finalNumExonSep","finalNumExonMerge","finalNumExonSimp","finalNumAA"];
writeTable(halSpeciesNamesP,specSummaryColP,specSummaryM,specSummaryMFn);

#stophere();

################################################
#   Get the genome/protein sequence   #
################################################
#Adjust the BED file to be be correct for pulling out sequence

#orthProteinFn = "output/orthProteins.%s.1.1.4.fasta" % (useRefProtId);

specSummary2M = []; #Updated matrix of summary across species
speciesNamesWriteP = []; #Which species to print out the information for

if getSeqB:

    print("Loading the sequence from FASTA files\n");
    numHighQuality = 0;

    orthProteinP = [];
    orthNucsP = [];

    #Add the protien information from the anchor species
    for curAnchorSpec in useAnchorP :
        curProtId = "%s-%s" % (useRefProtId,curSpecName);
        curAnchorProt = anchorProteinD[curAnchorSpec];
        curAnchorProtSeqR = SeqRecord(curAnchorProt,id=curProtId,description="pfenning_anchor_zoonomiaCACTUS");
        orthProteinP.append(curAnchorProtSeqR);

        curAnchorNucs = anchorNucsD[curAnchorSpec];
        curAnchorNucsSeqR = SeqRecord(curAnchorNucs,id=curProtId,description="pfenning_anchor_zoonomiaCACTUS");
        orthNucsP.append(curAnchorNucsSeqR);


    #Loop through the HAL files
    for curSpecLineP in halInfoM :
    #for curSpecLineP in halInfoM[1:10] :

        #print(curSpecLineP);

        curSpecName = curSpecLineP[1];
        curAnchorSpec = curSpecLineP[15];
        curProtId = "%s-%s" % (useRefProtId,curSpecName);

        hloOutbed2Fn = "%s/protBedHalName.%s.%s.%s.bed" % (outBedDir,curSpecName,useRefProtId,fileSuffixStr);
        curSpecFasta = curSpecLineP[11];

        curFastaB = True;

        if( not path.exists(curSpecFasta) ) :
            print("FASTA NOT FOUND - " + curSpecFasta + " - " + curSpecName);

        if(curSpecName in halSpeciesNamesP) :

            #0=curSpecName,1=curFragNum,2=curScafNum,3=curAANum,4=curMissingExon,5=curMultFragExon,6=curMultScafExon;
            curSpecSummaryP = specSummaryM[spec2HalIndD[curSpecName]]; #Load the stats for the gene

            #Test to see if naming convention can make the chromosomes work
            #Test to see if the species is in the stats list
            if(not curSpecLineP[8].startswith("N/A") and not curSpecLineP[8] == "" and path.exists(curSpecFasta)) :

                ### Pull out the sequences from the BEDs ###
                #getSeqStr = "bedtools getfasta -fi %s -bed %s" % (curSpecFasta,hloOutbedFn);
                getSeqStr = "bedtools getfasta -fi %s -bed %s" % (curSpecFasta,hloOutbed2Fn);
                print(getSeqStr);

                #### Read in sequences from the command line processing ####
                curSpecCdsSeq = Seq("");
                curExonSeq = Seq(""); # If in use, building up the sequence separated exons
                #if(curSpecName == "Mus_musculus") :

                curHalBedM = halSpecBedM[spec2HalIndD[curSpecName]];
                curBedPos = 0; #The current position in the bed file

                curFixExon = 0;
                curChooseMerge = 0;
                curChooseSep = 0;
                curSimp = 0;

                #print(len(curHalBedM));
                print("Starting " + curSpecName)

                process = subprocess.Popen(getSeqStr, shell=True, stdout=subprocess.PIPE);
                for curIoLine in io.TextIOWrapper(process.stdout, encoding="utf-8"):
                    #print(line.rstrip());
                    if(not curIoLine.startswith(">")) :

                        newSpecCdsSeq = Seq(curIoLine.rstrip(),ambiguous_dna);
                        curStrand = curHalBedM[curBedPos][5];
                        #print("Seq 1 - " + newSpecCdsSeq);


                        ### If the exon is split apart mulitple into fragments, order could be upended ####
                        #Test to see if the next exon is that same. If it is, adjust where you add the sequence
                        #If the at the of the file, or exon name changes, add to sequence
                        curEndExonB = False;
                        #print(curBedPos);

                        if( (curBedPos+2) > len(curHalBedM)) :
                            curEndExonB = True;
                        else :
                            if (curHalBedM[curBedPos][3] != curHalBedM[curBedPos+1][3]) : #New exon or at the end
                                curEndExonB = True;

                        if(curEndExonB) :
                            if( len(curExonSeq) > 0) : #If there is sequence stored in the exon, add it
                                #print("\tExon Seq Add = ");
                                #print(curExonSeq);
                                #print(newSpecCdsSeq);
                                curExonSeq = curExonSeq + newSpecCdsSeq;

                            testSpecCdsSeq = newSpecCdsSeq;
                            #print("Seq 2 - " + testSpecCdsSeq);

                            if(curStrand == "-") :
                                testSpecCdsSeq = testSpecCdsSeq.reverse_complement();

                            ### Adjust the nucleotide position of the merged sequence based on frame shifts ###
                            curAnchorNucNum = int(curHalBedM[curBedPos][4]); #Nucleotide length of the anchor
                            curMergeNucNum = int(len(testSpecCdsSeq)); #Nucleotide lenght of the merge sequence
                            curMergeNucDiff = curAnchorNucNum - curMergeNucNum;
                            curNumMergeBad = 0;
                            if( not (curMergeNucDiff % 3) == 0 ) : #Try to make adjustments if exon frame changes
                                curReviseTestSeq1 = "";
                                curReviseTestSeq2 = "";
                                curFixExon = curFixExon + 1;
                                if ( (curMergeNucDiff % 3) == 1) : # Remove two nucleotides to make things work (from begin or end)
                                    curReviseTestSeq1 = testSpecCdsSeq[2:];
                                    curReviseTestSeq2 = testSpecCdsSeq[:-2];
                                if ( (curMergeNucDiff % 3) == 2) : # Remove two nucleotides to make things work (from begin or end)
                                    curReviseTestSeq1 = testSpecCdsSeq[1:];
                                    curReviseTestSeq2 = testSpecCdsSeq[:-1];
                                curTmpCds1 = curSpecCdsSeq + curReviseTestSeq1;
                                curTmpCds2 = curSpecCdsSeq + curReviseTestSeq2;
                                curTmpProt1 = curTmpCds1.translate();
                                curTmpProt2 = curTmpCds2.translate();
                                curNumBad1 = str(curTmpProt1).count("X") + str(curTmpProt1).count("*");
                                curNumBad2 = str(curTmpProt2).count("X") + str(curTmpProt2).count("*");

                                #If either scores are equal, it likely means that removing from end is best
                                if(curNumBad1 < curNumBad2) : #Remove from beginning if you reduce stop codons
                                    testSpecCdsSeq = curReviseTestSeq1;
                                    curNumMergeBad = curNumBad1;
                                else :
                                    testSpecCdsSeq = curReviseTestSeq2;
                                    curNumMergeBad = curNumBad2;

                            curExon = curHalBedM[curBedPos][3];
                            curExonIden = curSpecName + "_" + curExon;

                            if curExonIden in mergedExonsP : #If exon had been merged, pull up original sequence

                                curExonSepBedM = halSpecBedSepD[curExonIden];
                                curSepSeq = "";
                                for curSepExonRowP in curExonSepBedM :
                                    curAdjStart = int(curSepExonRowP[1]) - int(curHalBedM[curBedPos][1]);
                                    curAdjStop = int(curSepExonRowP[2]) - int(curHalBedM[curBedPos][1]);
                                    curSepSeq  = curSepSeq + curExonSeq[curAdjStart:curAdjStop]; #Subset non rev-comp

                                curSepSeqTest = curSepSeq;
                                if(curStrand == "-") :
                                    curSepSeqTest=curSepSeq.reverse_complement();

                                ### Adjust the nucleotide position of the merged sequence based on frame shifts ###
                                #curAnchorNucNum = int(curHalBedM[curBedPos][4]); #Nucleotide length of the anchor
                                curSepNucNum = int(len(curSepSeqTest)); #Nucleotide lenght of the merge sequence
                                curSepNucSepDiff = curAnchorNucNum - curSepNucNum;
                                curNumSepBad = 0;
                                if( not (curSepNucSepDiff % 3) == 0 ) : #Try to make adjustments if exon frame changes
                                    curReviseTestSeq1 = "";
                                    curReviseTestSeq2 = "";
                                    if ( (curSepNucSepDiff % 3) == 1) : # Remove two nucleotides to make things work (from begin or end)
                                        curReviseTestSeq1 = curSepSeqTest[2:];
                                        curReviseTestSeq2 = curSepSeqTest[:-2];
                                    if ( (curSepNucSepDiff % 3) == 2) : # Remove two nucleotides to make things work (from begin or end)
                                        curReviseTestSeq1 = curSepSeqTest[1:];
                                        curReviseTestSeq2 = curSepSeqTest[:-1];
                                    curTmpCds1 = curSpecCdsSeq + curReviseTestSeq1;
                                    curTmpCds2 = curSpecCdsSeq + curReviseTestSeq2;
                                    curTmpProt1 = curTmpCds1.translate();
                                    curTmpProt2 = curTmpCds2.translate();
                                    curNumSepBad1 = str(curTmpProt1).count("X") + str(curTmpProt1).count("*");
                                    curNumSepBad2 = str(curTmpProt2).count("X") + str(curTmpProt2).count("*");
                                    #If either scores are equal, it likely means that removing from end is best
                                    if(curNumSepBad1 < curNumSepBad2) : #Remove from beginning if you reduce stop codons
                                        curSepSeqTest = curReviseTestSeq1;
                                        curNumSepBad = curNumSepBad1;
                                    else :
                                        curSepSeqTest = curReviseTestSeq2;
                                        curNumSepBad = curNumSepBad2;

                                ### Compare the separated and new sequence ####
                                #Default is merged if at all possible, decide whether to separate out again
                                useSepB = False; #By default use the merged sequence;

                                curMergeCds = curSpecCdsSeq + testSpecCdsSeq; #Add the the reverse complemented sequence to test
                                curMergeProt = curMergeCds.translate();
                                curMergeNumBad =  str(curMergeProt).count("X") + str(curMergeProt).count("*");

                                curSepCds = curSpecCdsSeq + curSepSeqTest; #Add the the reverse complemented sequence to test
                                curSepProt = curSepCds.translate();
                                curSepNumBad =  str(curSepProt).count("X") + str(curSepProt).count("*");

                                #curAnchorNucNum
                                #print("\tMergeAnchorDiff: " + curMergeNucDiff)
                                #print("\tSepAnchorDiff: " + curSepNucSepDiff)

                                #print("\tTargetMod3: " + (curAnchorNucNum % 3));
                                #print("\tMergeMod3: " + str(len(testSpecCdsSeq) % 3));
                                #print("\tSepMod3: " + (curSepSeqTest % 3));


                                if(curSepNumBad <  curMergeNumBad) : #Assume that you are moving onto the next exon and add the sequence
                                    useSepB = True;
                                    curSpecCdsSeq = curSpecCdsSeq + curSepSeqTest;
                                    #print("Adding exon to full sequence - use separated");
                                    #print("Seq 3.1 - " + curSepSeqTest);
                                    curChooseSep = curChooseSep + 1;
                                else :
                                    curSpecCdsSeq = curSpecCdsSeq + testSpecCdsSeq;
                                    #print("Adding exon to full sequence - use merged");
                                    #print("Seq 3.2 - " + testSpecCdsSeq);
                                    curChooseMerge = curChooseMerge + 1;

                                curExonSeq = Seq(""); #
                            else :
                                #curExonSeq = curExonSeq + newSpecCdsSeq; #If not new exon or end, just store seq
                                #curExonSeq = curExonSeq + curSepSeqTest;
                                #print("\tMergeMod3_2: " + str(len(testSpecCdsSeq) % 3));
                                curSpecCdsSeq = curSpecCdsSeq + testSpecCdsSeq; #If exon hasn't been merged, just add the seq test
                                #print("Adding exon to full sequence - no merge check");
                                #print("Seq 3.3 - " + testSpecCdsSeq);
                                curSimp = curSimp + 1;
                        else :
                            curExonSeq = curExonSeq + newSpecCdsSeq; #If not new exon or end, just store seq
                            #print("Adding Sequence to Exon, not support");
                            #print("Adding Exon Seq - " + newSpecCdsSeq);



                        #print(newSpecCdsSeq);
                        #curSpecCdsSeq = curSpecCdsSeq + newSpecCdsSeq;
                        curBedPos = curBedPos + 1;
                    else :

                        #print(curIoLine.rstrip())
                        curIoLine2 = curIoLine.rstrip();
                        #print(newExonCoordD[curIoLine2])


                curSpecProtSeq = curSpecCdsSeq.translate();
                print(curSpecProtSeq);
                curNumAaBroken = str(curSpecProtSeq).count("*");
                curNumAaAmbig = str(curSpecProtSeq).count("X");
                #curAaLen = len(curSpecProtSeq);
                curAaLen = len(curSpecCdsSeq)/3;

                #0=curSpecName,1=curFragNum,2=curScafNum,3=curAANum,4=curMissingExon,5=curMultFragExon,6=curMultScafExon;
                #curSpecSummaryP = specSummaryM[spec2indD[curSpecName]]; #Load the stats for the gene
                #curFixExon = 0,curChooseMerge = 0,curChooseSep = 0,curSimp = 0;

                curSpecSummaryP.append(curAaLen); #Add amino acid length based on final sequence
                curSpecSummaryP.append(curNumAaBroken); #Add the number that is broken
                curSpecSummaryP.append(curNumAaAmbig); #Add the number that is ambiguous
                curSpecSummaryP.append(curFixExon); #Add the where the merged exon is fixed
                curSpecSummaryP.append(curChooseMerge); #Add the number where the merged exon is chosed
                curSpecSummaryP.append(curChooseSep); #Add the number where the sepearated exons are chosen
                curSpecSummaryP.append(curSimp); #Add the number where the full exon is added without conflict
                curSpecSummaryP.append(curAnchorSpec); #Add the anchor species for this gene
                curSpecSummaryP.append(specDistD[curSpecName + "-" + curAnchorSpec]); #Add the anchor distance species for this gene
                curSpecSummaryP.append(specDistD[curSpecName + "-" + referenceSpecies]); #Add the reference distance species for this gene


                ##### Set criteria for including in protein alignment ########
                printCurProtB = True;

                curAnchorProtLen = len(anchorProteinD[curAnchorSpec]);

                if curNumAaBroken > curAaLen*0.01 : #Do not print protein if > 1/100 amino acids are stop codons
                    printCurProtB = False;

                minAALen = curAnchorProtLen - (curAnchorProtLen*.3);
                maxAAlen = curAnchorProtLen + (curAnchorProtLen*.3);

                if curAaLen < minAALen or curAaLen > maxAAlen  : #Do not print protein if the length is too different
                    printCurProtB = False;

                #print([curNumAaBroken,curAaLen,curSpecSummaryP[6],curSpecSummaryP[4]]);

                if(printCurProtB):
                    print("%s - High Quality" % (curSpecName));
                    numHighQuality = numHighQuality + 1;

                    curSpecProtSeqR = SeqRecord(curSpecProtSeq,id=curProtId,description="pfenning_HAL_zoonomiaCACTUS");
                    orthProteinP.append(curSpecProtSeqR);

                    curSpecNucsSeqR = SeqRecord(curSpecCdsSeq,id=curProtId,description="pfenning_HAL_zoonomiaCACTUS");
                    orthNucsP.append(curSpecNucsSeqR);

                else :
                    print("%s - Low Quality lengthAA-%f broken-%i" % (curSpecName,curAaLen,curNumAaBroken));



                curSpecSummaryP.append(printCurProtB)
                specSummary2M.append(curSpecSummaryP);
                speciesNamesWriteP.append(curSpecName);



            else :

                curSpecSummaryP.append(-1); #Add amino acid length based on final sequence
                curSpecSummaryP.append(-1); #Add the number that is broken
                curSpecSummaryP.append(-1); #Add the number that is ambiguous
                curSpecSummaryP.append(-1); #Add the where the merged exon is fixed
                curSpecSummaryP.append(-1); #Add the number where the merged exon is chosed
                curSpecSummaryP.append(-1); #Add the number where the sepearated exons are chosen
                curSpecSummaryP.append(-1); #Add the number where the full exon is added without conflict
                curSpecSummaryP.append(curAnchorSpec); #Add the anchor species for this gene
                curSpecSummaryP.append(specDistD[curSpecName + "-" + curAnchorSpec]); #Add the anchor distance species for this gene
                curSpecSummaryP.append(specDistD[curSpecName + "-" + referenceSpecies]); #Add the reference distance species for this gene

                curSpecSummaryP.append(False); #Add the number that is ambiguous
                specSummary2M.append(curSpecSummaryP);
                speciesNamesWriteP.append(curSpecName);



    print("Number Correct " + str(numHighQuality));

    if(printProtB) :

        orthProteinFn = "output/fastaRes/orthProteins.%s.%s.fasta" % (useRefProtId,fileSuffixStr);
        SeqIO.write(orthProteinP, orthProteinFn, "fasta");

        orthNucsFn = "output/fastaRes/orthNucleotides.%s.%s.fasta" % (useRefProtId,fileSuffixStr);
        SeqIO.write(orthNucsP, orthNucsFn, "fasta");

        specSummaryM2Fn = "output/fastaRes/specSummary2M.%s.%s.txt" % (useRefProtId,fileSuffixStr);
        specSummaryCol2P = specSummaryColP + ["curAaLen2","curAaStopCodon","curNumAaAmbig","curFixExon","curChooseMerge","curChooseSep","curSimp","achorSpecies","anchorDistance","referenceDistance","curHighQual"];

        writeTable(speciesNamesWriteP,specSummaryCol2P,specSummary2M,specSummaryM2Fn);

################################################
#   Post-run test code   #
################################################

#Score == length, so count the total length of the mapped sequence
#awk -vOFS='\t' -vFS='\t' '{sum += $5} END {print sum/3}' protBed/protBed.NP_068576.1.1.1.1.bed
#awk -vOFS='\t' -vFS='\t' '{sum += $3-$2 + 1} END {print sum/3}' protBed/protBed.NP_001358344.1.1.1.1.bed

#hg38.fa
#bedtools getfasta -fi /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/200MammalsFastas/hg38.fa -bed xSpecProtBed/protBedHalName.Homo_sapiens.NP_068576.1.1.1.1.bed


################################################
#   Post-run test use HALPER   #
################################################

#Run halLiftover
#halLiftover --bedType 6 /data/pfenninggroup/align/mam200/alignment_files/200m-v1.hal Homo_sapiens /home/apfennin/projects/covidTest/orthGene/protBed/protBed.NP_068576.1.1.1.1.bed Mus_musculus /home/apfennin/projects/covidTest/orthGene/xSpecProtBed/protBed.NP_068576.Mus_musculus.1.1.1.1.bed

#awk -vOFS='' -vFS='\t' '{print $1,":",$2,"-",$3}' /home/apfennin/projects/covidTest/orthGene/xSpecProtBed/protBed.NP_068576.Mus_musculus.1.1.1.1.bed

#awk -vOFS='\t' -vFS='\t' '{sum += $3-$2 + 1} END {print sum/3}' /home/apfennin/projects/covidTest/orthGene/xSpecProtBed/protBed.NP_068576.Mus_musculus.1.1.1.1.bed

#awk -vOFS='\t' -vFS='\t' '{sum += $3-$2 + 1} END {print sum/3}' xSpecProtBed/protBed.NP_068576.Mus_musculus.1.1.1.1.bed
#awk -vOFS='\t' -vFS='\t' '{sum += $3-$2 + 1} END {print sum/3}' xSpecProtBed/protBed.Ziphius_cavirostris.NP_068576.1.1.1.1.bed
#awk -vOFS='\t' -vFS='\t' '{sum += $3-$2 + 1} END {print sum/3}' xSpecProtBed/protBed.Acinonyx_jubatus.NP_068576.1.1.1.1.bed

########## Post-test run - fasta ############

#bedtools getfasta -fi /projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/200MammalsFastas/mm10.fa -bed xSpecProtBed/protBed.NP_068576.Mus_musculus.1.1.1.1.bed






############# Output Columns ##############

#ID
#TssEns
#TssSym
#TssDist
#TssProtEns
#TssProtSym
#TssProtDist
#Location
#FullEns
#FullSym
#FullDist
#FullProtEns
#FullProtSym
#FullProtDist
