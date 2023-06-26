#! /bin/bash

#Script Romane Guilbaud
#Creation data: 24/02/2023
#Last modification: 24/02/2023

####--------General Goal: Counts the number of reads before and after trimming and filtereing

###############################################################################################################

# #Repertories
# MY_REP=/media/virologie/My_Passport/Romane/SIID-Romane

# FASTQ_BICHAT=/media/virologie/My_Passport/SIID_pour_Romain/FASTQ_BAM/FASTQ_GZ
# FASTQ_BICHAT2=/media/virologie/My_Passport/SIID_pour_Romain/FASTQ_BAM/FASTQ
# BICHAT_REP=/media/virologie/My_Passport/Romane/Data_bichat

FASTQ_BICHAT_ACC=/media/virologie/My_Passport/Romane/Bichat_accurate/Seq_7
BICHAT_ACC_REP=/media/virologie/My_Passport/Romane/Bichat_accurate/Seq_7

# FASTQ_PITIE=/media/virologie/MyPassport/SIID_pitie/J0/fastq/
# PITIE_REP=/media/virologie/My_Passport/Romane/Data_pitie

# Files
# cd $FASTQ_BICHAT
# FILES_FASTQ_B=*[0-9].fastq.gz

# cd $FASTQ_BICHAT_ACC
# FILES_FASTQ_Ba=*[0-9].fastq.gz

# cd $FASTQ_BICHAT2
# FILES_FASTQ_B2=*[0-9].fastq

# cd "$FASTQ_PITIE"
# FILES_FASTQ_P=$(find $FASTQ_PITIE -type f -name "*.fastq.gz")

################################################################################################################
cd $BICHAT_ACC_REP
echo -e "Sample\tNb_raw_reads\tnb_reads_filtered" > nb_reads.txt

### Samples Bichat 

# accurate

while read SAMPLE; do
    echo $SAMPLE
    NAME=$(echo $SAMPLE | cut -d ";" -f 1)
    echo $NAME
    NBREADS="$(zcat $NAME.fastq.gz | wc -l)"
    NBFILTERED="$(zcat ${NAME}_filtered.fastq.gz | wc -l)"
    
    echo -e "$NAME\t"$NBREADS"\t"$NBFILTERED >> nb_reads.txt
    
done < $BICHAT_ACC_REP/samplelist_seq7.csv

# # fastq.gz
# for f in `ls $FILES_FASTQ_B`
# do
# 	cd $FASTQ_BICHAT
#     current_name=$(basename $f .fastq.gz)
#     echo "$current_name"
#     NBREADS="$(zcat $f | wc -l)"
#     cd $BICHAT_REP
#     NBFILTERED="$(zcat ${current_name}_filtered.fastq.gz | wc -l)"

#     cd $MY_REP
# 	echo -e "$current_name\t"$NBREADS"\t"$NBFILTERED >> nb_reads.txt

# done


# # fastq
# for f in `ls $FILES_FASTQ_B2`
# do
#     cd $FASTQ_BICHAT2
#     current_name=$(basename $f .fastq)
#     echo "$current_name"
#     NBREADS="$(wc -l $f )"
#     cd $BICHAT_REP
#     NBFILTERED="$(zcat ${current_name}_filtered.fastq.gz | wc -l)"

#     cd $MY_REP
#     echo -e "$current_name\t"$NBREADS"\t"$NBFILTERED >> nb_reads.txt

# done


# #### Samples Pitie
# # fastq.gz
# for f in `ls $FILES_FASTQ_P`
# do
#     cd $FASTQ_PITIE
#     current_name=$(basename $f .fastq.gz)
#     echo "$current_name"
#     NBREADS="$(zcat $f | wc -l)"
#     cd $PITIE_REP
#     NBFILTERED="$(zcat ${current_name}_filtered.fastq.gz | wc -l)"

#     cd $MY_REP
#     echo -e "$current_name\t"$NBREADS"\t"$NBFILTERED >> nb_reads.txt

# done


