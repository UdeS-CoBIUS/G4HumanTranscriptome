Welcome to G4Conserve's documentation!
======================================

This project aims to predict G-quadruplexes (G4) in the whole human's
transcriptome.

Data availability
=================

If the upper link is not owrking here are the main step to retrieve all input files:

* Go to the ensembl `archive Jul 2018 v93`_
* Go to the **biotmart** download tool
* Select the database "Ensembl Genes 93" and the dataset "Human genes (GRCh38.p12)"
* Then, in filter select one chromosomes from **region** (1 to 22 and X Y), you
  need to download one chromosome by one chromosome
* For the **first file type** you need to select the **attribute** Structures
  and select in order : Gene stable ID, Transcript stable ID,
  Chromosome/scaffold name, Gene type, 5' UTR start, 5' UTR end, 3' UTR start,
  3' UTR end, Exon region start (bp), Exon region end (bp), Exon rank in
  transcript, Strand. This file is named **HS_transcript_unspliced_chrY.txt**
  for the Y chromosome.
* The **second file type** contains the **attribute** sequences of Unspliced
  Gene and as header the features : Gene stable ID, Gene start (bp), Gene end
  (bp). This type of file are named **HS_gene_unspliced_chrY.txt**
* The **third file type** is an **attribute** feature to get transcript
  class : Gene stable ID, Transcript stable ID, Chromosome/scaffold name,
  Transcript type. This file type is contained in a special directory
  "transcriptType" and file are named like **transcriptType_chrY**
* The last files are whole chromosomes fasta file, downloaded from the
  FTP : ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/. Only
  **Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz** file are downloaded for each chromosomes.

.. _archive Jul 2018 v93: http://jul2018.archive.ensembl.org/index.html

Prerequisites
-------------

To use all scripts in this project you will need to have python 2.7.

You will also need to install G4RNA Screener from
`Michelle Scott gitlab`_.

.. _Michelle Scott gitlab: http://gitlabscottgroup.med.usherbrooke.ca/J-Michel/g4rna_screener.git

Following is procedure to install G4RNA screener ::

    git clone http://gitlabscottgroup.med.usherbrooke.ca/J-Michel/g4rna_screener.git
    sudo apt install python-pip
    pip install biopython
    pip install numpy
    pip install pandas
    pip install PyBrain
    pip install regex
    pip install scipy
    cd PATH/TO/PYTHON/dist-packages/    or    cd PATH/TO/PYTHON/site-packages/
    sudo -i
    wget https://dev.mysql.com/get/Downloads/Connector-Python/mysql-connector-python-2.1.4.tar.gz
    tar -xzf mysql-connector-python-2.1.4.tar.gz
    cd mysql-connector-python-2.1.4
    python setup.py install

Installing
----------
There is no particular installation because their is only individual home made
scripts. Here is a soft descriptino of them:

* ReplaceInformation.py -> This script aims to create
  junction between exons and also compute the coordinates of introns.
  To do that, this script takes as a file with informations about transcripts.
  This file is avaible in Ensembl and downable with Biomart.
* G4RNA Screener -> used on fasta files with genes and junction.
* all G4Annotation -> filters pG4 and add there location and biotype. Inputs are
  output from G4RNA Screener, the parsed file containing introns coordinates and
  also a file containing transcripts's biotype (available on Ensembl).
* G4Calcul.py -> computes all statistics and enrichment about our pG4.
* Coverage.py -> computes coverage between different data set. Takes as input
  3 data set : pG4, G4RNA G4 and G4 from rG4seq study.
* recurentFunction -> little library of functions that are often used.

Glossary
========

.. toctree::
   :maxdepth: 2

   computeDensities
   countLocation
   countTranscript
   FilterG4DependingOnTr
   G4Annotation
   G4AnnotationCirculaire
   G4AnnotationPIRNA
   G4AnnotationTRNA
   GenerateExonJunctionSequences
   generateShuffle
   getDataFig
   getMainDensities
   readRandom
   recurentFunction
   RunBlast
   readBlast
   readBlastpG4Q

Here is a global pipeline of how to use those scripts : |H|

.. |H| image:: /images/Pipeline.png

Command lines
=============

First step, on a local machine : configuration and data retrival
----------------------------------------------------------------

Go to the folder were you want to make this project or create one. You need
to launch scripts from this folder or you will have to change paths.

.. code-block:: bash
   :linenos:
   :emphasize-lines: 0

   mkdir G4prediction
   cd G4prediction

Then, we configure this folder, retrieve scripts and then data

.. code-block:: bash
  :linenos:
  :emphasize-lines: 0

   mkdir Data
   cd Data/
   for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do mkdir chr$chr; done
   mkdir transcriptType
   mkdir Fasta
   cd Fasta
   wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.*
   rm Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz
   gunzip Homo_sapiens.GRCh38.dna.chromosome.*
   cd ../../
   git clone https://github.com/UdeS-CoBIUS/G4HumanTranscriptome.git

At this point, you can download data, as explained in Data availability. Then
You can launch the first script to generate Exon-Exon sequences, and upload your
data to your cluster of calcul (Here it's graham from compute canada).

.. code-block:: bash
 :linenos:
 :emphasize-lines: 0

   for chr in $(ls Data/ | grep chr); do time python G4Conserve/script/GenerateExonJunctionSequences.py -chr $chr; done
   scp -r -p Data/ userID@graham.computecanada.ca:~/scratch/G4prediction

Second step on cluster : generate shuffled dataset, and predict G4
------------------------------------------------------------------

Here follow all step to create tmp folder, results folder and to predict G4
and annotate pG4. all bash scripts are given in the folder sBash, but path
should be changed depending on your cluster or user ID.

.. code-block:: bash
  :linenos:
  :emphasize-lines: 0

   # create a few folder for tmp files
   for chr in $(ls ~/scratch/G4prediction/Data/ | grep chr); do mkdir ~/scratch/TestRepro/Data/$chr/SplitFile;done
   for chr in $(ls ~/scratch/G4prediction/Data/ | grep chr); do mkdir ~/scratch/TestRepro/Data/$chr/CSVFile;done
   # generate random and then split all fasta file
   for chr in $(ls ~/scratch/G4prediction/Data/ | grep chr); do sbatch ../../Scripts/generateRandom.sh $chr ;done
   for chr in $(ls ~/scratch/G4prediction/Data/ | grep chr); do sbatch ../../Scripts/fasta_enum_sbatch.sh $chr $(ls ~/scratch/TestRepro/Data/$chr | grep 'gene_unspliced'); done
   for chr in $(ls ~/scratch/G4prediction/Data/ | grep chr); do sbatch ../../Scripts/fasta_enum_sbatch.sh $chr HS_transcript_unspliced_$chr_Sequence.txt; done
   for chr in $(ls ~/scratch/G4prediction/Data/ | grep chr); do for file in $(ls ~/scratch/TestRepro/Data/$chr | grep Shuffle); do sbatch ../../Scripts/fasta_enum_sbatch.sh $chr $file; done; done
   # predict G4 with G4RNA screener, it create one task for each splited file
   # with a max of 999 runing at the same time
   sbatch ../../Scripts/submitG4RNAScreener.sh
   # create folder for results
   mkdir ~/scratch/G4prediction/Results
   mkdir ~/scratch/G4prediction/Results/perChromosome
   # filter the output of G4RNA screener and annotate pG4r for both dataset
   for chr in $(ls ~/scratch/G4prediction/Data/ | grep chr); do sbatch ../../Scripts/G4Annotation.sh $chr; done
   for chr in $(ls ~/scratch/G4prediction/Data/ | grep chr); do sbatch ../../Scripts/readRandom.sh $chr; done
   # if needed, SplitFile folders and CSVFile folders can be deleted

Third step on local machine : compute densities and generate figures
--------------------------------------------------------------------

.. code-block:: bash
 :linenos:
 :emphasize-lines: 0

   scp -r -p userID/graham.computecanada.ca:~/scratch/G4prediction/ ./
   # don't forget to go to your main folder
   mkdir Documents/G4prediction/Results/All
   # join all results chromosomes files or annotation file into one for the
   # entire human genome
   for file in $(ls Documents/G4prediction/Results/perChromosome | grep 'G4InTran'); do cat Documents/TestRepro/Results/perChromosome/$file; done >> Documents/TestRepro/Results/All/HS_All_G4InTranscript.txt
   for file in $(ls Results/perChromosome | grep shuffle); do cat Results/perChromosome/$file; done >> Results/All/pG4r_shuffle.csv
   for file in $(ls Data/transcriptType/ | grep -v All); do cat Data/transcriptType/$file; done >> Data/transcriptType/transcriptType_All.txt
   for chr in $(ls Data | grep chr); do for file in $(ls Data/$chr | grep Shuffle); do cat Data/$chr/$file | grep '>'; done; done >> Results/All/Locations.txt
   for chr in $(ls Data/ | grep chr); do cat 'Data/'$chr'/HS_transcript_unspliced_'$chr'_Index.txt'; done >> Data/transcriptType/HS_transcript_unspliced_All.txt
   mkdir Results/Figures
   # create first a table with densities at the location level for all subclasses
   # then compute a table for each figure
   python ./G4Conserve/scripts/getMainDensities.py
   python ./G4Conserve/scripts/getDataFig.py

At this point, you should have all results done. The last script to use is the
rmd one, were a rmarkdown on html was made.

rG4seq dataset comparison
-------------------------

For this part, those scripts are used : RunBlast , readBlast and readBlastpG4Q.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
