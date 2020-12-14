#!/bin/bash
#
# Script Name: download_ribosomal_models.sh
#
# Author: Helena Cooper
# Last edited: 14/12/2020
#
# Description: This script automates the process of downloading Rfam covariance models and Pfam hidden Markov models.
#

#--------------------------------------------------------------------------------------------------------------------------

### Constants & Set-up

#--------------------------------------------------------------------------------------------------------------------------

IFS=$'\n'                                     # Set delimiter as new line during for loops.
DATE=$(date +%y%m%d )                         # Date parameter.

if command -v cmcalibrate &> /dev/null        # If INFERNAL executables not available, then exit script.
then
    :
else
    echo "Please add the INFERNAL executables to your PATH."
    exit
fi

mkdir $DATE-models/ && cd $DATE-models/

#--------------------------------------------------------------------------------------------------------------------------

# rRNA Covariance Models

#--------------------------------------------------------------------------------------------------------------------------

# 16S rRNA (RF00177)
wget -O RF00177.stk "http://rfam.xfam.org/family/RF00177/alignment?acc=RF00177&format=stockholm&download=1"
cmbuild -F RF00177-16S.cm RF00177.stk
cmcalibrate RF00177-16S.cm
rm RF00177.stk

# 23S rRNA (RF02541)
wget -O RF02541.stk "http://rfam.xfam.org/family/RF02541/alignment?acc=RF02541&format=stockholm&download=1"
cmbuild -F RF02541-23S.cm RF02541.stk
cmcalibrate RF02541-23S.cm
rm RF02541.stk

#--------------------------------------------------------------------------------------------------------------------------

# SSU Hidden Markov Models

#--------------------------------------------------------------------------------------------------------------------------

# SSU uS4 (PF00163)
wget -O PF00163-uS4.hmm "https://pfam.xfam.org/family/PF00163/hmm"

# SSU uS12 (PF00164)
wget -O PF00164-uS12.hmm "https://pfam.xfam.org/family/PF00164/hmm"

# SSU uS7 (PF00177)
wget -O PF00177-uS7.hmm "https://pfam.xfam.org/family/PF00177/hmm"

# SSU uS3 C-terminal (PF00189)
wget -O PF00189-uS3C.hmm "https://pfam.xfam.org/family/PF00189/hmm"

# SSU uS19 (PF00203)
wget -O PF00203-uS19.hmm "https://pfam.xfam.org/family/PF00203/hmm"

# SSU uS14 (PF00253)
wget -O PF00253-uS14.hmm "https://pfam.xfam.org/family/PF00253/hmm"

# SSU uS15 (PF00312)
wget -O PF00312-uS15.hmm "https://pfam.xfam.org/family/PF00312/hmm"

# SSU uS2 (PF00318)
wget -O PF00318-uS2.hmm "https://pfam.xfam.org/family/PF00318/hmm"

# SSU uS5 N-terminal (PF00333)
wget -O PF00333-uS5N.hmm "https://pfam.xfam.org/family/PF00333/hmm"

# SSU uS10 (PF00338)
wget -O PF00338-uS10.hmm "https://pfam.xfam.org/family/PF00388/hmm"

# SSU uS17 (PF00366)
wget -O PF00366-uS17.hmm "https://pfam.xfam.org/family/PF00366/hmm"

# SSU uS9 (PF00380)
wget -O PF00380-uS9.hmm "https://pfam.xfam.org/family/PF00380/hmm"

# SSU uS8 (PF00410)
wget -O PF00410-uS8.hmm "https://pfam.xfam.org/family/PF00410/hmm"

# SSU uS11 (PF00411)
wget -O PF00411-uS11.hmm "https://pfam.xfam.org/family/PF00411/hmm"

# SSU uS13 (PF00416)
wget -O PF00416-uS13.hmm "https://pfam.xfam.org/family/PF00416/hmm"

# SSU uS5 C-terminal (PF03719)
wget -O PF03719-uS5C.hmm "https://pfam.xfam.org/family/PF03719/hmm"

#--------------------------------------------------------------------------------------------------------------------------

# LSU Hidden Markov Models

#--------------------------------------------------------------------------------------------------------------------------

# LSU uL2 (PF00181)
wget -O PF00181-uL2.hmm "https://pfam.xfam.org/family/PF00181/hmm"

# LSU uL22 (PF00237)
wget -O PF00237-uL22.hmm "https://pfam.xfam.org/family/PF00237/hmm"

# LSU uL14 (PF00238)
wget -O PF00238-uL14.hmm "https://pfam.xfam.org/family/PF00238/hmm"

# LSU uL16 (PF00252)
wget -O PF00252-uL16.hmm "https://pfam.xfam.org/family/PF00252/hmm"

# LSU uL23 (PF00276)
wget -O PF00276-uL23.hmm "https://pfam.xfam.org/family/PF00276/hmm"

# LSU uL5 (PF00281)
wget -O PF00281-uL5.hmm "https://pfam.xfam.org/family/PF00281/hmm"

# LSU uL3 (PF00297)
wget -O PF00297-uL3.hmm "https://pfam.xfam.org/family/PF00297/hmm"

# LSU uL11 (PF00298)
wget -O PF00298-uL11.hmm "https://pfam.xfam.org/family/PF00298/hmm"

# LSU uL30 (PF00327)
wget -O PF00327-uL30.hmm "https://pfam.xfam.org/family/PF00327/hmm"

# LSU uL6 (PF00347)
wget -O PF00347-uL6.hmm "https://pfam.xfam.org/family/PF00347/hmm"

# LSU uL10 (PF00466)
wget -O PF00466-uL10.hmm "https://pfam.xfam.org/family/PF00466/hmm"

# LSU uL13 (PF00572)
wget -O PF00572-uL13.hmm "https://pfam.xfam.org/family/PF00572/hmm"

# LSU uL4 (PF00573)
wget -O PF00573-uL4.hmm "https://pfam.xfam.org/family/PF00573/hmm"

# LSU uL1 (PF00687)
wget -O PF00687-uL1.hmm "https://pfam.xfam.org/family/PF00687/hmm"

# LSU uL15 (PF00828)
wget -O PF00282-uL15.hmm "https://pfam.xfam.org/family/PF00828/hmm"

# LSU uL29 (PF00831)
wget -O PF00831-uL29.hmm "https://pfam.xfam.org/family/PF00831/hmm"

# LSU uL18 (PF00861)
wget -O PF00861-uL18.hmm "https://pfam.xfam.org/family/PF00861/hmm"

# LSU uL11 N-terminal (PF03946)
wget -O PF03946-uL11N.hmm "https://pfam.xfam.org/family/PF03946/hmm"

# LSU uL24 (PF17136)
wget -O PF17136-uL24.hmm "https://pfam.xfam.org/family/PF17136/hmm"

#--------------------------------------------------------------------------------------------------------------------------

cd ..
