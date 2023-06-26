#! /bin/bash

#Script Romain Coppee
#Creation data: 09/02/2020
#Last modification: 09/07/2020

####--------General Goal: 

###############################################################################################################

#Adding MINIMAP2 to the PATH environment
export PATH=/home/virologie/Téléchargements/minimap2/:$PATH


#location of PORECHOP software
#PORECHOP=/home/virologie/Documents/Porechop

#Location of the Reference genome
REF_GEN=/home/virologie/Documents/sars_cov_reference/sars_cov_2_ref_genome.fasta

# #Location of the fastq files
# FASTQ_REP=/media/virologie/My_Passport/SIID_pour_Romain/FASTQ_BAM/FASTQ_GZ

# #Location of the work directory
# MY_REP=/media/virologie/My_Passport/Romane/Data_bichat

#Location of the fastq files
FASTQ_REP=/media/virologie/My_Passport/Romane/Bichat_accurate/fastq_7

#Location of the work directory
MY_REP=/media/virologie/My_Passport/Romane/Bichat_accurate/Seq_7

#FILES contains the fastq file for each sample
#cd $FASTQ_REP
FILES_FASTQ=*[0-9].fastq.gz

#FILES contains the fastq file filtered and trimmed by nanofilt for each sample
#cd $MY_REP
FILES_FASTQ_FILTER=*filtered.fastq.gz

#FILES contains the fixed, indexed BAM file for each sample
FILES_FIX=*.fix

#FILES contains the SAM file for each sample
FILES_SAM=*.sam

#FILES contains the BAM file for each sample
FILES_BAM=*.bam

#FILES contains the indexed BAM file for each sample
FILES_BAI=*.bai

FILES_PILEUP=*.pileup


FILES_TSV=*.tsv

###############################################################################################################
#####-------Goal 1: 

##### Regroup all the fastq.qz into one per sample 
while read SAMPLE; do
  echo $SAMPLE
  BARCODE=$(echo $SAMPLE | cut -d ";" -f 2)
  NAME=$(echo $SAMPLE | cut -d ";" -f 1)
  echo $BARCODE
  echo $NAME
    
  cat $FASTQ_REP/${BARCODE}/*.fastq.gz > $MY_REP/$NAME.fastq.gz
done < $MY_REP/samplelist_seq7.csv


#cd $FASTQ_REP
cd $MY_REP

for f in `ls $FILES_FASTQ`
do
    current_name=$(basename $f .fastq.gz)
    echo "$current_name"
    gunzip -c $f | python3 NanoFilt -q 10 --headcrop 10 | gzip > $MY_REP/${current_name}_filtered.fastq.gz
    echo "NANOFILT PROCESSED $f"
done

#cd $MY_REP

for f in `ls $FILES_FASTQ_FILTER`
do
    current_name=$(basename $f _filtered.fastq.gz)
    echo "$current_name"
    minimap2 -ax map-ont -t 8 $REF_GEN $f > $current_name.sam
    echo "MINIMAP2 PROCESSED $f"
done

for f in `ls $FILES_SAM`
do
    current_name=$(basename $f .sam)
    samtools view -b -S -@ 8 $f > $current_name.bam
    echo "SAM to BAM PROCESSED $f"
    rm $f
done

for f in `ls $FILES_BAM`
do
    current_name=$(basename $f .bam)
    samtools sort -@ 8 $f -o $current_name.sorted.bam
    echo "SORTING PROCESSED $f"
    rm $f
done

for f in `ls $FILES_BAM`
do
    samtools index $f
    echo "INDEXATION PROCESSED $f"
done

for f in `ls $FILES_BAM`
do
    current_name=$(basename $f .sorted.bam)
    samtools stats $f > ${current_name}_stats.txt
    echo "STATS PROCESSED $f"
done

for f in `ls $FILES_BAM`
do
    current_name=$(basename $f .sorted.bam)
    echo $f
    #samtools mpileup --max-depth 8000 -a -f $REF_GEN -r NC_045512.2:21563-25384 $f -o $current_name.pileup
    samtools mpileup --max-depth 8000 -a -B -q 60 -Q 15 -f $REF_GEN $f -o $current_name.pileup
done


for f in `ls $FILES_PILEUP`
do
    current_name=$(basename $f .pileup)
    echo "$current_name"
    perl /home/virologie/Documents/scripts/extract_info_from_pileup.pl -p $f -o ./
    mv info_from_pileup.tsv $current_name.tsv
    awk -v new="$current_name" '$1=new' < $current_name.tsv > $current_name.tsv.fix
done

rm $FILES_PILEUP
for f in `ls $FILES_FIX`
do 
    mv -- "$f" "${f%.fix}"
    echo "fixing PROCESSED $f"
done

for f in `ls $FILES_TSV`
do
    current_name=$(basename $f .tsv)
    echo "$f"
    sed -i '1 i\Id_sample Position n_A n_C n_G n_T n_indel major_nucl complement qual' $f
    awk -v OFS="\t" '$1=$1' $f > $f.fix
done

rm $FILES_TSV
for f in `ls $FILES_FIX`
do 
    mv -- "$f" "${f%.fix}"
    echo "fixing PROCESSED $f"
done