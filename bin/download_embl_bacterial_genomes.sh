#!/bin/bash
#
# Script Name: download_embl_bacterial_genomes.sh
#
# Author: Helena Cooper
# Last edited: 09/11/2020
#
# Description: Download all bacterial and Haloarcula marismortui genomes from ENA.
#

#--------------------------------------------------------------------------------------------------------------------------

### Constants & Set-up

#--------------------------------------------------------------------------------------------------------------------------

d=$(date +%y%m%d)   # Date parameter.

### Check if required dependency is available before running script.
if command -v sreformat &> /dev/null
then
    :
else
    echo "Please add the sreformat to your PATH."
    exit
fi

mkdir -p $d-bacteria    # Create folder for downloaded genomes.
cd $d-bacteria

#--------------------------------------------------------------------------------------------------------------------------

### Main pipeline

#--------------------------------------------------------------------------------------------------------------------------

### Create summary of downloaded gemomes.
curl -G https://www.ebi.ac.uk/genomes/bacteria.details.txt > $d-bacteria.details.txt

### Fetch emble files.
grep -v ^# $d-bacteria.details.txt | perl -lane 'print "curl -G \42http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/embl/$F[0]\42 > $F[0]\.embl" if (not -s "$F[0]\.embl");' | sh

### Fetch Haloarcula marismortui genome.
curl -G https://www.ebi.ac.uk/ena/browser/api/embl/AY596297.1?download=true > AY596297.1.embl

### Move the contigs & obviously incomplete sequences aside.
mkdir -p contigs
grep ^CO *embl | tr "." "\t" | cut -f 1,2 | uniq | perl -lane 'print "curl -G \42http://www.ebi.ac.uk/ena/data/view/$F[0]&display=fasta&download&filename=$F[0].fasta\42 > $F[0].$F[1].fasta && cp $F[0].$F[1].embl contigs/" if (not -e "$F[0]\.$F[1]\.fasta")' | sh
grep contig $d-bacteria.details.txt | perl -lane 'print "mv $F[0]* contigs/"' | sh

### Generate bacteria.fa.
ls *embl | perl -lane 'if (/(\S+)\.embl/){print "sreformat fasta $F[0]  | tr -d \47;\47 > $1\.fasta" }' | sh
rm -f $d-bacteria.fa

rm -rf AY596297.1.embl    # Doesn't format correctly with bacteria.taxonomy.txt pipe. 

### Summarise taxonomy data.
echo -ne "EMBL_ACC\tSPECIES\tNCBI_CLASSIFICATION" > $d-bacteria.taxonomy.txt
grep ^O *embl | tr ":" "\t" | perl -lane 'if (($F[0] ne $e) or (/U00096.*Escherichia\./)){print "$e\t$s\t$c"; ($e,$s,$c)=("","","");} if(/(\S+)\tOS\s+(.*)/){ $e=$1; $s=$2; }elsif(/(\S+)\tOC\s+(.*)/){$c.=$2;}' >> $d-bacteria.taxonomy.txt

### Add Haloarcula marismortui taxonomy data.
echo -e "AY596297.1.embl\tHaloarcula marismortui AT 43049\tArchaea; Euryarchaeota; Stenosarchaea group; Halobacteria; Halobacteriales; Haloarculaceae; Haloarcula" >> $d-bacteria.taxonomy.txt

#####################################################################

rm -rf *.embl
