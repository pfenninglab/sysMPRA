#!/bin/bash

#SBATCH -a 2-35
##SBATCH -a 1-1
#SBATCH -p pfen1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -o /projects/MPRA/andreas/21_06_arrayProc/tmp/logs/bcCount.slurm.out.%A.%a.txt
#SBATCH -e /projects/MPRA/andreas/21_06_arrayProc/tmp/logs/bcCount.slurm.err.%A.%a.txt
##SBATCH -o /projects/MPRA/andreas/21_06_arrayProc/tmp/logs/bcCount.slurm.out.%j.txt
##SBATCH -e /projects/MPRA/andreas/21_06_arrayProc/tmp/logs/bcCount.slurm.err.%j.txt
##SBATCH -t 0-24:00:00

cd /projects/MPRA/andreas/21_06_arrayProc/

#/home/apfennin/meme/bin/ame -oc bindingSite/test/wilcInc -control out15/fasta/regSeq.verm14_hg38.5.2.1.fa -bgformat 0 -scoring max -method ranksum out15/fasta/regSeq.verm14_hg38_prgInc.5.2.1.fa /home/apfennin/tools/tfbs/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme
#/home/apfennin/meme/bin/ame -oc bindingSite/test/fishInc -control out15/fasta/regSeq.verm14_hg38.5.2.1.fa -bgformat 0 -scoring max -method fisher out15/fasta/regSeq.verm14_hg38_prgInc.5.2.1.fa /home/apfennin/tools/tfbs/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme

MYFILE=/projects/MPRA/andreas/21_06_arrayProc/input/sampV2.names.2.1.txt
MYID=$SLURM_ARRAY_TASK_ID

name=$(awk "NR==${MYID}" $MYFILE)

echo $name

python3 arrayProc.2.1.1.py /projects/MPRA/MPRA/MPRAi/MPRAi_v2_NovaSeq_Genewiz2/ input/MPRA_Array1_sequences $name "$name"_R1_001.fastq.gz countsV2/ qcV2/

#bamToBed -i cantin2016/bam/$name.raw.bam | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c -f > cantin2016/tag/$name.raw.tag.gz

#zcat cantin2016/tag/$name.raw.tag.gz | sort -k1,1 -k2,2n -k3,3n > cantin2016/tag/$name.raw.sort.tag

#intersectBed -sorted -g cantin2016/taeGut2.sort.chrom.sizes  -c -a peaks/tagGut2.4sort.bed -b cantin2016/tag/$name.raw.sort.tag > tmp/counts/cant.$name.counts.1.bed

#$name

#sbatch scripts/mpraBcCount.2.1.1.sh
