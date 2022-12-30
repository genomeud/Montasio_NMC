# lattoinnesto
Code for the paper **"Integrated approach to explore the microbial biodiversity of natural milk culture for Montasio cheese production"**

## Required software
- R (we used R 3.6.1, any newer version should do)
- kraken2 (we used 2.1.2)

## Used data
- Reads can be retrieved from SRA, under accession **PRJNA916949**
- Physicochemical parameter data can be found in the **06_chempar** folder of this repository

## Instructions
- Create a folder for the project e.g. myhome/lattoinnesto
- Download the reads from SRA accession **PRJNA916949** and store them in the project folder, in a subfolder called "01_raw_reads"
- Download the repository in the project folder
- Edit the lines in the **00_lattoinnesto_workflow.sh** file where local variables are stored. Remember to uncomment the variables. 
- Run the **00_lattoinnesto_workflow.sh** file. I suggest you execute the commands step by step by copying and pasting them in the terminal, so that you can spot mistakes immediately. There may be errors due to wrong folder structure (my fault!), or due to missing software (not my fault!)   
