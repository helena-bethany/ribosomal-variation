# Investigating bacterial ribosomal sequence variation in regards to future structural and antibiotic research.

* [Required Dependencies](#required-dependencies)

* [Downloading the Initial Data](#downloading-the-initial-data)

* [Annotation of Ribosomal Genes and Phylogeny Analysis](#annotation-of-ribosomal-genes-and-phylogeny-analysis)

* [Multidimensional Scaling and Selection of Proposed Structures](#multidimensional-scaling-and-selection-of-proposed-structures)

* [References](#references)

-------------------------------------------------------------------------------------------------

### Required Dependencies

| Software | Program(s) Used | Reference |  
|:------------:|:------------:|:------------:|
|[barrnap 0.9](https://github.com/tseemann/barrnap)| barrnap | - |
|[Bedtools 2.29.0](https://bedtools.readthedocs.io/en/latest/content/installation.html)| bedtools | Quinlan, 2014 |
|[EMBOSS 6.6.0.0](http://emboss.sourceforge.net/what/)| transeq | Rice et al., 2000 |
|[HMMER 3.1](hmmer.org)| hmmsearch; hmmalign; esl-translate ; esl-reformat| - |
|[INFERNAL 1.1.2](http://eddylab.org/infernal/)| cmcalibrate; cmsearch; cmbuild | Nawrocki & Eddy, 2013 |
|[phylip 3.6.7](https://evolution.genetics.washington.edu/phylip.html)| dnaml ; proml |Felenstein, 2009|

-------------------------------------------------------------------------------------------------

### Downloading the Initial Data

The script [download_embl_bacterial_genomes.sh](bin/download_embl_bacterial_genomes.sh) was used to download all 3,758 bacterial genomes, in addition to the Archaean *Haloarcula marismorui* genome (Accession: [AY596297.1](https://www.ebi.ac.uk/ena/browser/view/AY596297)), from the European Nucleotide Archive (Amid et al., 2019; Baliga et al., 2004). The output includes a dated folder (eg: ```201112-bacteria/```) containing the genome sequences (eg: ```AY596297.1.fasta```), a summary of the downloaded genomes (eg: ```201112-bacteria.details.txt```) and a summary of taxonomy data (eg: ```201112-bacteria.taxonomy.txt```).

```
# Example of download_embl_bacterial_genomes.sh usage
bash download_embl_bacterial_genomes.sh
```

The script [download_ribosomal_models.sh](bin/download_ribosomal_models.sh) was used to download covariance models from Rfam v14.3 and hidden Markov models from Pfam v33.1 in order to annotate ribosomal genes (Kalvari et al., 2018; El-Gebali et al., 2019). This script also list the IDs for each downloaded model. The output includes a dated folder (eg: ```201112-models/```) containing the calibrated covariance models (eg: ```RF00177-16S.cm```) and hidden Markov models (eg: ```PF00416-uS13.hmm```).

```
# Example of download_embl_bacterial_genomes.sh usage
bash download_ribosomal_models.sh
```

--------------------------------------------------------------------------------------------------

### Annotation of Ribosomal Genes and Phylogeny Analysis

The script [ribosome_sequence_analysis.sh](bin/ribosome_sequence_analysis) was used to annotate ribosomal genes and to generate a maximum likelihood phylogeny tree for each type of ribosomal gene. There are three required input: the folder containing the downloaded genome files, the appropriate covariance or hidden Markov model and the name of the sequence being annotated. Further options are available and can be listed using help (```-h```). The main outputs include a dated folder (eg: ```201112-16S-sequences/```) containing the annotated ribosomal sequences (eg: ```16S-AY596297.1.fa```)  and a newick outtree file which includes the metadata for each species (eg: ```201112-16S-outtree```).

```
# Example of ribosome_sequence_analysis.sh usage
bash ribosome_sequence_analysis.sh -g 201112-bacteria -m RF00177-16S.cm -i 16S
```

##### Ribosomal RNA

Ribosomal RNA (rRNA) genes were annotated using ```barrnap``` with an 80\% sequence threshold (```--reject 0.80```), which meant that sequences that were not 80\% of the expected length for the rRNA were excluded from analysis. For *H. marismortui*, the archaea model (```--kingdom arc```) was to annotate the rRNAs instead of the bacteria model (```--kingdom bac```). Line 108 of barrnap was edited so that bitscores were reported, rather than the E-value, to allow the rRNA and protein methods to be comparable. The annotated sequences were then filtered to include the paralogue with the hightest bit score, and one sequence was retained per species. To create the multiple sequence alignment input for the phylogeny tree, ```cmsearch``` was used to align the annoated sequences to a rRNA covariance model using an E-value of 1E-6 (```-E 1E-6```) to ensure each sequence only aligned to the model once (Nawrocki & Eddy, 2013). The generated Stockholm file was then converted to phylip using ```esl-reformat``` was used to convert the Stockholm file to phylip, and ```dnaml``` was used to create the maximum likelihood tree (Felsenstein, 2009). If an example has not been provided below, then default parameters were used for the executable.

```
# Altered Line 108 in barrnap
my $score = defined $x[13] ? $x[13] : '.';

# Example of barrnap usage
barrnap --quiet --kingdom bac --reject 0.80 CP009789.1.fasta 

# Example of cmsearch usage
cmsearch -E 1E-6 -A 201112-16S-aln.stk ../RF00177-16S.cm seqdb-16S 
```

##### Ribosomal Protein

Prior to annotaing ribosomal protein genes, the genome files were underwent a six-frame translation using ```transeq``` and open-reading frame predictions using ```esl-translate``` to account for potentially inconsistent or absent genome annotations (Rice et al., 2000). Ribosomal proteins were annotated using ```hmmsearch``` with an 80\% sequence threshold, which meant that sequences that were not 80\% of the expected length for the protein, based on the ```LENG``` variable in the hidden Markov model, were excluded from analysis. The same pipeline was used for *H. marismortui* as universally conserved ribosomal proteins should be conserved across all domains of life. The annotated sequences were then filtered to include all unique domains with the hightest bit score and combined into one sequence, with one protein being retained per species. Large subunit protein uL6 was filtered separately to ensure only one domain was retained, as it contains two identical functional domains which introduce errors in the multiple alignment file (Golden et al., 1993). To create the multiple sequence alignment input for the phylogeny tree, ```hmmalign``` was used to align the annoated sequences to a protein covariance model. The generated Stockholm file was then converted to phylip using ```esl-reformat``` was used to convert the Stockholm file to phylip, and ```proml``` was used to create the maximum likelihood tree (Felsenstein, 2009). If an example has not been provided below, then default parameters were used for the executable.

```
# Example of transeq usage
transeq -sequence CP009789.1.fasta -outseq CP009789.1.pep -frame=6 

# Example of hmmsearch usage
hmmsearch -E 1E-6 -A 201112-16S-aln.stk ../RF00177-16S.cm seqdb-16S 
```

---------------------------------------------------------------------------------------------------

### Multidimensional Scaling and Selection of Proposed Structures

-----------------------------------------------------------------------------------------------------

### References

Amid C, Alako BTF, Balavenkataraman Kadhirvelu V, Burdett T, Burgin J, Fan J, Harrison PW, Holt S, Hussein A, Ivanov E, et al. 2020. The European Nucleotide Archive in 2019. Nucleic Acids Res 48: D70–D76.

Baliga NS, Bonneau R, Facciotti MT, Pan M, Glusman G, Deutsch EW, Shannon P, Chiu Y, Weng RS, Gan RR, et al. 2004. Genome sequence of Haloarcula marismortui: a halophilic archaeon from the Dead Sea. Genome Res 14: 2221–2234.

El-Gebali S, Mistry J, Bateman A, Eddy SR, Luciani A, Potter SC, Qureshi M, Richardson LJ, Salazar GA, Smart A, et al. 2019. The Pfam protein families database in 2019. Nucleic Acids Res 47: D427–D432.

Felsenstein J. 2009. PHYLIP (Phylogeny Inference Package) version 3.7a. Distributed by the author. http://www.evolution.gs.washington.edu/phylip.html. http://www.evolution.gs.washington.edu/phylip.html (Accessed September 7, 2019).

Golden BL, Ramakrishnan V, White SW. 1993. Ribosomal protein L6: structural evidence of gene duplication from a primitive RNA binding protein. EMBO J 12: 4901–4908.

Kalvari I, Argasinska J, Quinones-Olvera N, Nawrocki EP, Rivas E, Eddy SR, Bateman A, Finn RD, Petrov AI. 2018. Rfam 13.0: shifting to a genome-centric resource for non-coding RNA families. Nucleic Acids Res 46: D335–D342.

Nawrocki EP, Eddy SR. 2013. Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics 29: 2933–2935.

Quinlan AR. 2014. BEDTools: The Swiss-Army Tool for Genome Feature Analysis. Curr Protoc Bioinformatics 47: 11.12.1–34.

Rice P, Longden I, Bleasby A. 2000. EMBOSS: the European Molecular Biology Open Software Suite. Trends Genet 16: 276–277.
