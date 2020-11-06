#!/bin/bash
#
# Script Name: ribosome_sequence_analysis.sh
#
# Author: Helena Cooper
# Last edited: 06/11/2020
#
# Description: This script predicts ribsomal RNA or protein sequence from genomes, and then generates a maximum likelihood
#              phylogeny tree.
#

#--------------------------------------------------------------------------------------------------------------------------

### Constants

#--------------------------------------------------------------------------------------------------------------------------

IFS=$'\n'                                     # Set delimiter as new line during for loops.
DATE=$(date +%y%m%d )                         # Date parameter.
SCRIPT_NAME=`basename $0 | sed -e "s/.sh$//"` # Name of the script.

# List of all pathogen species of interest:
PATHOGENS="Acinetobacter_baumannii|Actinomyces_israelii|Bacillus_anthraces|Bacteroides_fragilis|Bartenella_henselae|Bordetella_pertussis|Borrelia_recurrentis|Brucella_melitensis|Burkholderia_cepacia|Campylobacter_jejuni|Capnocytophaga_canimorsus|Cardiobacterium_hominis|Chlamydia_trachomatis|Chlamydophila_pneumonia|Clostridium_perfrigens|Clostridium_difficile|Corynebacterium_diptheriae|Coxiella_burnetii|Cutibacterium_acnes|Ehrlichia_canis|Enterobacter_cloacae|Enterococcus_faecium|Erysipelothrix_rhusiopathiae|Escherichia_coli|Francisella_tularensis|Fusobacterium_necrophorum|Haemophilus_influenzae|Helicobacter_pylori|Kingella_kingae|Klebsiella_pneumoniae|Legionella_pneumophila|Leptospira_interrogans|Listeria_monocytogenes|Moraxella_catarrhalis|Morganella_morganii|Mycobacterium_leprae|Mycoplasma_pneumoniae|Neisseria_gonorrhoeae|Neisseria_menigitis|Nocardia_asteroides|Pastuerella_multocida|Proteus_mirabilis|Providencia_stuartii|Rickettsia_rickettsii|Salmonella_enterica|Serrate_marcescens|Shigella_dysenteriae|Staphylococcus_epidermidis|Streptococcus_pneumoniae|Streptococcus_pyogenes|Streptoccus_viridans|Treponema_pallidum|Ureaplasma_urealyticum|Vibrio_cholerae|Yersinia_pestis|Yersinia_pseudotuberculosis"

# List of all species of interest with solved structures:
STRUCTURES="Deinococcus_radiodurans|Escherichia_coli|Haloarcula_marismortui|Thermus_thermophilus|Enterococcus_faecalis|Mycobacterium_tuberculosis|Staphylococcus_aureus|Pseudomonas_aeruginosa"

#--------------------------------------------------------------------------------------------------------------------------

### Functions

#--------------------------------------------------------------------------------------------------------------------------

### Calculation of float numbers.

calc() { awk "BEGIN{print $*}"; }

#--------------------------------------------------------------------------------------------------------------------------

### Check if executables are on $PATH prior to running main pipeline.

check_executable () {

    ### $1 is an executable.
    ### $2 is the package.

    if command -v $1 &> /dev/null
    then
        :
    else
        echo "Please add the "$2" executables to your PATH."
        exit
    fi
}

#--------------------------------------------------------------------------------------------------------------------------

### Checks if Haloarcula marismortui genome has been downloaded prior to running main pipeline.

check_hm () {

    ### No user input.

    if grep -q "Haloarcula marismortui" $GENOME_FOLDER/*taxonomy.txt
    then
        :
    else
        echo "Please download Haloarcula marismortui genome and add to bacterial genomes folder."
        exit
    fi

}

#--------------------------------------------------------------------------------------------------------------------------

### Ensure barrnap has been altered prior to running main pipeline.

check_barrnap() {

    if command -v barrnap.hacked &> /dev/null
    then
        BARRNAP=barrnap.hacked 
    elif command -v barrnap &> /dev/null
    then
        read -p "Has line 108 of the barrnap script been changed to `my $score = defined $x[13] ? $x[13] : '.';` (y/n)?" ANSWER
        if [[ "$ANSWER" == "y" ]]
        then
            BARRNAP=barrnap
        else
            exit
        fi
    else
        echo "Please add the altered barrnap script (barrnap.hacked) to your PATH."
        exit
    fi

}

#--------------------------------------------------------------------------------------------------------------------------

### Display script help syntax.

help() {

  echo ""
  echo "Usage: ${SCRIPT_NAME}.sh -g genome_folder -m sequence_model -i seqID [ -v | –h | -s | -t ]"
  echo ""
  echo "  Required:"
  echo "    -g         Folder containing downloaded bacterial genomes."
  echo "    -m         CM or HMM for sequence of interest."
  echo "    -i         Custom sequence identifier (16S or uS12)."
  echo ""
  echo "  Optional:"
  echo "    -v         Verbose logging."
  echo "    -h         Script help."
  echo "    -s         Predict ribosomal sequence only."
  echo "    -t         Build ML phylogeny tree only."
  echo ""
  echo "  Dependencies:"
  echo "    INFERNAL v1.1.2"
  echo "    HMMER v3.1"
  echo "    barrnap v0.9"
  echo "    phylip v3.6.7"
  echo "    bedtools v2.29.0"
  echo ""

  return 0

}

#--------------------------------------------------------------------------------------------------------------------------

### Processes Command Line Arguments (and call help if the option has been specified).

process_args() {

  # Need : after anything that requires an option

  while getopts g:m:i:vhst OPTION
  do
    case "$OPTION" in
      g)  GENOME_FOLDER=${OPTARG};;
      m)  SEQ_MODEL=${OPTARG};;
      i)  SEQID=${OPTARG};;
      v)  USE_DEBUG="Y";;
      h)  USE_HELP="Y";;
      s)  SEQUENCE_ONLY="Y";;
      t)  TREE_ONLY="Y";;   
    esac
  done

  shift $(($OPTIND-1))

  # If help option specified, then call help and exit from the script.
  if [[ "$USE_HELP" == "Y" ]] ; then help ; exit ; fi
  
  if [ -z "${GENOME_FOLDER}" ] ; then echo "You must specify a genome folder for -g. See –h for full instructions." ; exit ; fi
  if [ -z "${SEQ_MODEL}" ] ; then echo "You must specify a CM or HMM for -m. See –h for full instructions." ; exit ; fi
  if [ -z "${SEQID}" ] ; then echo "You must specify a sequence identifier for -i. See –h for full instructions." ; exit ; fi

  return 0

}

#--------------------------------------------------------------------------------------------------------------------------

### Predict rRNA sequences and create a multiple sequence alignment.

generate_rRNA () {

    ### No user input.
    
    for GENOME in $FILES   # Parses through each genome file.
    do
        NAME=$( echo $GENOME | sed 's/.fasta//' | rev | cut -d '/' -f 1 | rev)   # Genome Accession Number.
        GEN=$FILE_ID-$NAME    # File identifier during for loop (format: date-seqID-genomeID)

        if [ -z "$( grep -v ">" $GENOME )" ]    # If no sequence in genome file, do not process.
        then
            echo $NAME >> $FILE_ID-removed.txt
        else
            SPECIES_NAME=$( grep "$NAME" $GENOME_FOLDER/*taxonomy.txt | cut -f 2 | tr -d "[]'" | cut -d ' ' -f 1,2 )   

            ### Generate rRNA sequence predictions using altered barrnap script ($BARRNAP defined in check_barrnap).
            
            if [[ "$SPECIES_NAME" == "Haloarcula marismortui" ]]  # Use Archaea model for H. marismortui, otherwise use bacteria model.
            then
                $BARRNAP --quiet --kingdom arc --reject 0.80 $GENOME > $GEN-output.txt
            else
                $BARRNAP --quiet --kingdom bac --reject 0.80 $GENOME > $GEN-output.txt 
            fi

            ### Sort file, and then select the sequence with the largest bitscore.
            
            sort -k6 -n -r $GEN-output.txt > $GEN-barrnap.txt
            CHECK=$(cat $GEN-barrnap.txt | wc -l)

            if (( $(echo "$CHECK == 1" | bc -l) ))    # If barrnap couldn't identify any rRNA above threshold (output only contains file header).
            then
                echo $NAME >> $FILE_ID-removed.txt
            else
                bedtools getfasta -fo $GEN.fa -fi $GENOME -bed $GEN-barrnap.txt   # Extract rRNA from genome file.
                VAR=$( echo $SEQID rRNA | tr " " "_" )    # Sequence type in format for parsing through barrnap output.
                VARIABLE=$( cat $GEN-barrnap.txt | grep -m 1 "$VAR" | cut -f 5 )    # Unique entry from barrnap output.
                if [ -z "$VARIABLE" ]   # If no unique variable, then specific rRNA type being genereated wasn't predicted.
                then
                    echo $NAME >> $FILE_ID-removed.txt
                else
                    grep $VARIABLE -A 1 $GEN.fa > $FOLDER/$SEQID-$NAME.fa   # Extract sequence into its own fasta file.
                    BITSCORE=$( grep $VARIABLE $GEN-barrnap.txt | cut -f 6 )    # Bitscore for extracted sequence.
                    echo $NAME,$SPECIES_NAME,$BITSCORE >> $FOLDER/$FILE_ID-bitscores.csv    # File used for later filtering step.
                fi
            fi
        fi

        rm -rf $GEN*.*    # Remove excess files.

    done
    
    ### Identify species with the most conserved rRNA sequence (highest bitscore - third column).
    # The first pipe sorts the file by ',' delimiter, first by the second column (names in alphabetical order) and then the third (high to low bitscore).
    # The second pipe keeps the first unique hit in the second column (species names), which is the individual for each species with the highest bitscore.

    UNIQUE_IDs=$( sort -t, -k2,2 -k3nr $FOLDER/$FILE_ID-bitscores.csv | awk -F "," '!a[$2]++' | cut -d ',' -f 1 )
    
    ### Create a multiple sequence alignment using filtered species.
    
    for UNIQUE in $UNIQUE_IDs ; do cat $FOLDER/$SEQID-$UNIQUE.fa >> $FOLDER/seqdb-$SEQID ; done
    cmsearch -E 1E-6 -A $FOLDER/$FILE_ID-aln.stk $SEQ_MODEL $FOLDER/seqdb-$SEQID &> /dev/null

} 

#--------------------------------------------------------------------------------------------------------------------------

### Predict ribosomal protein sequences and create a multiple sequence alignment.

generate_protein () {

    ### No user input.

    LENG=$( grep "LENG" $SEQ_MODEL | tr -d "LENG" | tr -d ' ' )   # Sequence length from the HMM. 
    NUMBER=$( calc 0.80*$LENG )   # Expected sequence length (80%).
    THRESHOLD=$( printf "%.0f\n" $NUMBER )   

    for GENOME in $FILES     # Parses through each genome file.
    do
        NAME=$( echo $GENOME | sed 's/.fasta//' | rev | cut -d '/' -f 1 | rev)    # Genome Accession Number.
        GEN=$FILE_ID-$NAME     # File identifier during for loop (format: date-seqID-genomeID)

        if [ -z "$( grep -v ">" $GENOME )" ]     # If no sequence in genome file, don't process.
        then
            echo $NAME >> $FILE_ID-removed.txt
        else
            esl-translate $GENOME > $GEN.pep    # 6-frame translation of ORFs.

            ### Reformat output so that there are no '\n' in the sequences.

            n=$( grep -c ">" $GEN.pep )
            num=$(( n*2 ))
            cat $GEN.pep | awk '/^>/ {printf("\n%s\n",$1);next; } { printf("%s",$1);}  END {printf("\n");}' | tail -"$num" > $GEN-translated.pep

            ### Predict ribosomal protein sequences.
            
            hmmsearch -E 1E-6 --tblout $GEN-results.tbl -o $GEN-output.tbl $SEQ_MODEL $GEN-translated.pep &> /dev/null

            if [ -z "$( grep -v "#" $GEN-results.tbl )" ]   # If no hmmsearch hit, don't process.
            then
                echo $NAME >> $FILE_ID-removed.txt
            else
                INFO=$( grep -v "#" $GEN-results.tbl | tr -s ' ' | cut -d ' ' -f 1,6 | sort -k2 -n | tail -1 )   # Obtain largest bit score if there are multiple hits.
                ID=$( echo $INFO | cut -d ' ' -f 1 )    # ID of the sequence with the best hit to HMM.
                SCORE=$( echo $INFO | cut -d ' ' -f 2 )   # Bitscore of the sequence with the best hit to HMM.

                START=$( sed -n "/>> $ID/,/^$/p" $GEN-output.tbl | grep -m 1 '!' | tr -s ' ' | cut -d ' ' -f 11 )  # Start position.
                END=$( sed -n "/>> $ID/,/^$/p" $GEN-output.tbl | grep -m 1 '!' | tr -s ' ' | cut -d ' ' -f 12 )    # End position.

                SEQUENCE=$( grep -A 1 -w ">$ID" $GEN-translated.pep | tail -1 | cut -c$START-$END )   # Extract ribosomal sequence from esl-translate output.

                LEN=$( expr length $SEQUENCE )    # Length of predicted sequence.
                if [ $LEN -lt $THRESHOLD ]    # If predicted sequence is shorter than threshold, do not process.
                then
                    echo $NAME >> $FILE_ID-removed.txt
                else
                    echo ">$NAME" >> $FOLDER/$SEQID-$NAME.fa    # Create individual fasta file for ribosomal protein sequence.
                    echo $SEQUENCE >> $FOLDER/$SEQID-$NAME.fa
                    SPECIES_NAME=$( grep "$NAME" $GENOME_FOLDER/*taxonomy.txt | cut -f 2 | tr -d "[]'" | cut -d ' ' -f 1,2 )
                    echo $NAME,$SPECIES_NAME,$SCORE >> $FOLDER/$FILE_ID-bitscores.csv    # File used for later filtering step.
                fi
            fi
        fi

        rm -rf $GEN*.*    # Remove excess files.
  
    done

    ### Identify species with the most conserved rRNA sequence (highest bitscore - third column).
    # The first pipe sorts the file by ',' delimiter, first by the second column (names in alphabetical order) and then the third (high to low bitscore).
    # The second pipe keeps the first unique hit in the second column (species names), which is the individual for each species with the highest bitscore.
    
    UNIQUE_IDs=$( sort -t, -k2,2 -k3nr $FOLDER/$FILE_ID-bitscores.csv | awk -F "," '!a[$2]++' | cut -d ',' -f 1 )

    ### Create a multiple sequence alignment using filtered species.
    
    for UNIQUE in $UNIQUE_IDs ; do cat $FOLDER/$SEQID-$UNIQUE.fa >> $FOLDER/seqdb-$SEQID ; done
    hmmalign -o $FOLDER/$FILE_ID-aln.stk $SEQ_MODEL $FOLDER/seqdb-$SEQID &> /dev/null

    DUPLICATES=$( grep -v "#" $FOLDER/$FILE_ID-aln.stk | grep . | cut -d '/' -f 1 | sort | uniq -d )    # List of species that had multiple domain alignments.

    if [ -z "$DUPLICATES" ]   # If no duplicates, skip formatting step.
    then
        :
    else 
        for DUPLICATE in $DUPLICATES 
        do
            ### Format first domain hit.
            
            line_one_first=$( grep "$DUPLICATE" $FOLDER/$FILE_ID-aln.stk | grep -v "#" | head -1 | cut -d '-' -f 1,2 )    # Seq for first domain.
            first_count=$( echo $line_one_first | wc -c )   # Length of sequence in first domain.
            first_count=$( calc $first_count-1 )    # Adjust so cut region is correct.
            line_two_first=$( grep "$DUPLICATE" $FOLDER/$FILE_ID-aln.stk | grep "#=GR" | head -1 | cut -c 1-$first_count )    # Sec Struc for first domain.

            ### Format second domain hit.
            
            line_one_third=$( grep "$DUPLICATE" $FOLDER/$FILE_ID-aln.stk | grep -v "#" | tr -s ' ' | cut -d ' ' -f 2 | tail -1 | rev | cut -d '-' -f 1 | rev ) # Seq for second domain.
            third_count=$( echo $line_one_third | wc -c )   # Length of sequence in second domain.
            third_count=$( calc $third_count-1 )    # Adjust so cut region is correct.
            line_two_third=$( grep "$DUPLICATE" $FOLDER/$FILE_ID-aln.stk | grep "#=GR" | tail -1 | rev | tr -s ' ' | cut -d ' ' -f 1 | cut -c 1-$third_count | rev )     # Sec Struc for second domain.

            ### Format region separting the two domains.
            
            second_count=$( grep "$DUPLICATE" $FOLDER/$FILE_ID-aln.stk | grep -v "#" | tail -1 | cut -c $first_count- | rev | cut -d '-' -f 2- | rev | wc -c )    # Total length of the connecting region.
            second_count=$( calc $second_count+$first_count )   # Adjust so cut region is correct.
            line_one_second=$( grep "$DUPLICATE" $FOLDER/$FILE_ID-aln.stk | grep -v "#" | head -1 | cut -c $first_count-$second_count )   # Seq for connecting region.
            line_two_second=$( grep "$DUPLICATE" $FOLDER/$FILE_ID-aln.stk | grep "#=GR" | head -1 | cut -c $first_count-$second_count )   # Sec Struc for connection region.

            ### Merge alignments into one line.
            
            LINE_ONE=$( echo $line_one_first$line_one_second$line_one_third )   # Reformatted first line (sequence).
            LINE_TWO=$( echo $line_two_first$line_two_second$line_two_third )   # Reformatted second line (secondary structure).

            ### Replace the two domain hits with a singular combined hit.

            grep -v "$DUPLICATE" $FOLDER/$FILE_ID-aln.stk | sed "/^#=GC PP_cons.*/i $LINE_ONE\n$LINE_TWO" > $FILE_ID-interim
            cp $FILE_ID-interim $FOLDER/$FILE_ID-aln.stk

        done
    fi

    rm -rf $FILE_ID-interim

}

#--------------------------------------------------------------------------------------------------------------------------

### Build a maximum likelihood phylogeny tree and incorporate metadata into outtree.

build_ML_tree () {

    ### No user input.

    ### Convert to Stockholm aligment to phylip format and build tree.
    
    esl-reformat phylip $FILE_ID-aln.stk | tr '.' '|' | tr ":" "|" > infile

    if type phylip >/dev/null 2>&1    # How the program is called depends on how it's installed.
    then
        printf 'Y\n' | phylip $PROGRAM
    else
        printf 'Y\n' | $PROGRAM
    fi

    ### Incorporate metadata into philip outtree.

    [ ! -f $FILE_ID-outtree ] && cp outtree $FILE_ID-outtree    # Rename outtree in case there is a formatting error with sed.

    for LINE in $( grep -v "NCBI_CLASSIFICATION" ../$GENOME_FOLDER/*taxonomy.txt )
    do
        ID=$( echo $LINE | cut -f 1 | cut -d '.' -f 1 )   # Genome accession number.
        NAME=$( echo $LINE | cut -f 2 | cut -d ' ' -f 1,2 | tr -d "[]'" | tr " " "_" )    # Species name.
        PHYLA=$( echo $LINE | cut -d ';' -f 2,3 | tr -d ' ' | tr ';' ' ' | tr " " "_" | tr -d '.' | tr '/' '-' )    # Phyla and Class.

        if [[ $PATHOGENS =~ $NAME ]]    # If species is a pathogen of interest, label as "Pathogenic".
        then
            TARGET="Pathogenic"
        elif [[ $STRUCTURES =~ $NAME ]]   # If species is has a solved structure, label as "Structural".
        then
            TARGET="Structural"
        else
            TARGET="Other"    # Otherwise, label as "Other".
        fi
   
        BACTERIA_ID=$( echo $ID $NAME $PHYLA $TARGET | tr " " "_" )   # Formatted metadata so it can be parsed in R.

        sed -i "s/$ID|./$BACTERIA_ID/g" $FILE_ID-outtree    # Replace genome ID with metadata.
   
    done

    sed -i "s/-/./g" $FILE_ID-outtree   # Reformat dashes for easier parsing in R.

}

#--------------------------------------------------------------------------------------------------------------------------

# Define main function/pipeline.

#--------------------------------------------------------------------------------------------------------------------------

main () {

  ### Process main arguments.

  process_args $*

  [[ "$USE_DEBUG" == "Y" ]] && set -xv
  [[ "$SEQ_MODEL" == *"hmm"* ]] && MODEL=protein || MODEL=rrna    # Whether protein or rRNA is being analysed.

  ### Check all dependencies are available.

  check_hm 
  [[ "$MODEL" == "rrna" ]] && check_barrnap
  check_executable hmmsearch HMMER
  check_executable cmsearch INFERNAL
  check_executable esl-translate HMMER/easal
  check_executable bedtools bedtools
  command -v phylip &> /dev/null || check_executable dnaml phylip   # Depends on how phylip was installed.

  ### Format starting files and variables.

  FILE_ID=$DATE-$SEQID
  FOLDER=$FILE_ID-sequences   # Folder that sequences and phylogeny tree will go into.
  FILES=$( ls $GENOME_FOLDER/*.fasta | sort -V ) # List of all genomic sequence files.
  
  [ -f $FILE_ID-removed.txt ] && rm -rf $FILE_ID-removed.txt    # List of genomes with no predicted sequence.
  [ -f $FILE_ID-bitscores.csv ] && rm -rf $FILE_ID-bitscores.csv    # List of all genomes with predicted sequence.
  [ -d $FOLDER/ ] && rm -rf $FOLDER/* || mkdir $FOLDER

  ### Run main pipeline.

  if [[ "$TREE_ONLY" == "Y" ]]    # If only the phylogeny tree has been specified, skip sequence prediction steps.
  then
      :
  elif [[ "$MODEL" == "rrna" ]]   
  then
      generate_rRNA
      PROGRAM=dnaml
      [[ "$SEQUENCE_ONLY" == "Y" ]] && exit  # If only sequence generation has been specified, exit script before tree-building.

  elif [[ "$MODEL" == "protein" ]]
  then
      generate_protein
      PROGRAM=proml
      [[ "$SEQUENCE_ONLY" == "Y" ]] && exit   # If only sequence generation has been specified, exit script before tree-building.
  else
      echo "Check input CM/HMM is valid."   # Since this loop is dependent on the $MODEL variable defined by the second user input.
      exit
  fi

  cd $FOLDER/ && build_ML_tree && cd ../    # Move into sequence folder prior to building phylogeny tree so that script can be run in parallel. 

}

### Run main function.

main $*   

# END
