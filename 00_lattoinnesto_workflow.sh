#!/bin/bash


##############################################
#
# Variables you will need to set based on your file system
#
# BASEDIR=/base/project/folder/
# KRAKBIN=/path/to/kraken2/executable
# BRACK_BIN=/path/bracken/executable
# FUNCDIR=/path/to/functions/folder/
#
###############################################


############################################
#
#  _              _              ___  
# | |            | |            |__ \ 
# | | ___ __ __ _| | _____ _ __    ) |
# | |/ / '__/ _` | |/ / _ \ '_ \  / / 
# |   <| | | (_| |   <  __/ | | |/ /_ 
# |_|\_\_|  \__,_|_|\_\___|_| |_|____|
#                                     
############################################
                                     



READDIR=${BASEDIR}/01_raw_reads
KRES=${BASEDIR}/03_kraken2_minikraken
KDB=${BASEDIR}/databases/minikraken2_v2_8GB_201904_UPDATE
mkdir -p $KRES


READPATTERN=_R1_001.fastq.gz
#End of parameters that you likely need to change

#Run kraken on the selected subreads
# PRE Execution
#s1 - classify
for aaa in $READDIR/*${READPATTERN}
do
read1=$(basename $aaa)
read2=${read1/_R1_001.fastq.gz/_R2_001.fastq.gz}
pref=${read1/_R1_001.fastq.gz/}
echo $read1
echo $read2
mkdir -p ${KRES}/logs
if [ -s ${KRES}/${pref}.kraken.report.txt ]
	then
	echo "Kraken already run, i will skip it"	
else
cd ${KRES} 
export TMPDIR=${KRES}
${KRAKBIN}/kraken2 --threads 8 --fastq-input --preload --paired --gzip-compressed --db $KDB $READDIR/${read1} $READDIR/${read2} --output ${KRES}/${pref}.kraken --use-names --report ${KRES}/${pref}.kraken.report.txt
fi
done


#We try using SILVA KRAKEN DB

KRES=${BASEDIR}/03_kraken2_silva
KDB=${BASEDIR}/databases/16S_SILVA138_k2db
mkdir -p $KRES


READPATTERN=_R1_001.fastq.gz
#End of parameters that you likely need to change

#Run kraken on the selected subreads
# PRE Execution
#s1 - classify
for aaa in $READDIR/*${READPATTERN}
do
read1=$(basename $aaa)
read2=${read1/_R1_001.fastq.gz/_R2_001.fastq.gz}
pref=${read1/_R1_001.fastq.gz/}
echo $read1
echo $read2
mkdir -p ${KRES}/logs
if [ -s ${KRES}/${pref}.kraken.report.txt ]
	then
	echo "Kraken already run, i will skip it"	
else
cd ${KRES} 
export TMPDIR=${KRES}
${KRAKBIN}/kraken2 --threads 8 --fastq-input --preload --paired --gzip-compressed --db $KDB $READDIR/${read1} $READDIR/${read2} --output ${KRES}/${pref}.kraken --use-names --report ${KRES}/${pref}.kraken.report.txt
fi
done


###########################################                                      
#  _                    _              
# | |                  | |             
# | |__  _ __ __ _  ___| | _____ _ __  
# | '_ \| '__/ _` |/ __| |/ / _ \ '_ \ 
# | |_) | | | (_| | (__|   <  __/ | | |
# |_.__/|_|  \__,_|\___|_|\_\___|_| |_|
#                                      
###########################################                                      


RESDIR=${BASEDIR}/03_kraken2_silva
THRESHOLD=10
READ_LEN=200
for CLASSIFICATION_LEVEL in S G
do
for READ1 in ${READDIR}/*${READPATTERN}
do
read1=$(basename $READ1)
read1=${read1/_R1_001.fastq.gz/}
echo $read1
KRES=$RESDIR
#Run bracken on kraken results
echo $THRESHOLD
#rm ${KRES}/*_${CLASSIFICATION_LEVEL}.bracken.txt
${BRACK_BIN}/bracken -d ${KDB} -i ${KRES}/${read1}.kraken.report.txt -o ${KRES}/${read1}_${CLASSIFICATION_LEVEL}.bracken.txt -r ${READ_LEN} -l ${CLASSIFICATION_LEVEL} -t ${THRESHOLD}
done
done


################################################
#  _____  _       _           _                     _                      
# |  __ \| |     | |         | |                   | |                     
# | |__) | | ___ | |_    __ _| |__  _   _ _ __   __| | __ _ _ __   ___ ___ 
# |  ___/| |/ _ \| __|  / _` | '_ \| | | | '_ \ / _` |/ _` | '_ \ / __/ _ \
# | |    | | (_) | |_  | (_| | |_) | |_| | | | | (_| | (_| | | | | (_|  __/
# |_|    |_|\___/ \__|  \__,_|_.__/ \__,_|_| |_|\__,_|\__,_|_| |_|\___\___|
#                                                                          
################################################

mkdir -p ${BASEDIR}/04_plots/
mkdir -p ${BASEDIR}/05_tables/
                                                                          
#Using minikraken db
Rscript ${FUNCDIR}/01_summarize_bracken.r \
        -I ${BASEDIR}/03_kraken2_minikraken \
        -O ${BASEDIR}/05_tables/bracken.txt \
        -G ${BASEDIR}/04_plots/species_kraken_barplot.pdf \
        -N 15 --remove_thermo TRUE

Rscript ${FUNCDIR}/01_summarize_bracken.r \
        -I ${BASEDIR}/03_kraken2_minikraken \
        -O ${BASEDIR}/05_tables/bracken.txt \
        -G ${BASEDIR}/04_plots/species_kraken_barplot.pdf \
        -N 15 --remove_thermo FALSE

mkdir -p ${BASEDIR}/04_plots_silva/
mkdir -p ${BASEDIR}/05_tables_silva/

#Using silva db
Rscript ${FUNCDIR}/01a_summarize_genus_bracken.r \
        -I ${BASEDIR}/03_kraken2_silva \
        -O ${BASEDIR}/05_tables_silva/genus_bracken.txt \
        -G ${BASEDIR}/04_plots_silva/genus_kraken_barplot.pdf \
        -N 15 --remove_Strepto TRUE

Rscript ${FUNCDIR}/01a_summarize_genus_bracken.r \
        -I ${BASEDIR}/03_kraken2_silva \
        -O ${BASEDIR}/05_tables_silva/genus_bracken.txt \
        -G ${BASEDIR}/04_plots_silva/genus_kraken_barplot.pdf \
        -N 15 --remove_Strepto FALSE


################################################
#   _____                      _       _   
#  / ____|                    | |     | |  
# | |     ___  _ __ _ __ _ __ | | ___ | |_ 
# | |    / _ \| '__| '__| '_ \| |/ _ \| __|
# | |___| (_) | |  | |  | |_) | | (_) | |_ 
#  \_____\___/|_|  |_|  | .__/|_|\___/ \__|
#                       | |                
#                       |_|                
#
################################################

#Run correlation for all chemical parameters. Do not worry if you get errors when you only have one parameter: in that case there are problems
#with the plot. I will eventually fix it
for INFILE in Acidi.txt aromi.txt pH.txt zuccheri.txt Treatment.txt
do
Rscript ${FUNCDIR}/02_correl.r \
        -A ${BASEDIR}/05_tables/bracken.txt \
        -O ${BASEDIR}/06_chempar/ \
        -C ${BASEDIR}/06_chempar/${INFILE}
done

#Use silva results
mkdir -p ${BASEDIR}/06_chempar_silva
for INFILE in Acidi.txt aromi.txt pH.txt zuccheri.txt Treatment.txt
do
Rscript ${FUNCDIR}/02_correl.r \
       -A ${BASEDIR}/05_tables_silva/genus_bracken.txt \
	   -O ${BASEDIR}/06_chempar_silva/ \
	   -C ${BASEDIR}/06_chempar/${INFILE}
done





################################################
#  _____  _                    _ _         
# |  __ \(_)                  (_) |        
# | |  | |___   _____ _ __ ___ _| |_ _   _ 
# | |  | | \ \ / / _ \ '__/ __| | __| | | |
# | |__| | |\ V /  __/ |  \__ \ | |_| |_| |
# |_____/|_| \_/ \___|_|  |___/_|\__|\__, |
#                                     __/ |
#                                    |___/ 
#
################################################


#Kraken species on minikraken (with and without Streptococcus thermophilus)
Rscript ${FUNCDIR}/05_diversity.r \
       -A ${BASEDIR}/05_tables/bracken_raw.txt -O ${BASEDIR}/05_tables/diversity.txt

Rscript ${FUNCDIR}/05_diversity.r \
	   -R 'Streptococcus thermophilus' \
       -A ${BASEDIR}/05_tables/bracken_raw.txt -O ${BASEDIR}/05_tables/diversity_no_Sth.txt

#Kraken genus on silva (with and withoout Streptococcus)
Rscript ${FUNCDIR}/05_diversity.r \
       -A ${BASEDIR}/05_tables_silva/genus_bracken_raw.txt -O ${BASEDIR}/05_tables_silva/diversity.txt

Rscript ${FUNCDIR}/05_diversity.r \
	   -R 'Streptococcus' \
       -A ${BASEDIR}/05_tables_silva/genus_bracken_raw.txt -O ${BASEDIR}/05_tables_silva/diversity_no_Strepto.txt


################################################
#   _____                      _       _   
#  / ____|                    | |     | |  
# | |     ___  _ __ _ __ _ __ | | ___ | |_ 
# | |    / _ \| '__| '__| '_ \| |/ _ \| __|
# | |___| (_) | |  | |  | |_) | | (_) | |_ 
#  \_____\___/|_|  |_|  | .__/|_|\___/ \__|
#                       | |                
#                       |_|                
#
################################################


#Run correlation for all chemical parameters with diversity parameters, and not with single species. 
#This is the version used for the paper, i.e.. based on species
for INFILE in Acidi.txt aromi.txt pH.txt zuccheri.txt Treatment.txt
do
Rscript ${FUNCDIR}/08_correl_div.r \
       -D ${BASEDIR}/05_tables/diversity.txt \
	   -O ${BASEDIR}/06_chempar/ \
       -C ${BASEDIR}/06_chempar/${INFILE}
done



#Use silva results at the genus level
# for INFILE in Acidi.txt aromi.txt pH.txt zuccheri.txt Treatment.txt
# do
# Rscript ${FUNCDIR}/08_correl_div.r \
       # -D ${BASEDIR}/05_tables_silva/diversity.txt \
	   # -O ${BASEDIR}/06_chempar_silva/ \
	   # -C ${BASEDIR}/06_chempar/${INFILE}
# done


#Plot correlation like shown in Zheng et al 
#Using minikraken
Rscript ${FUNCDIR}/09_correl_network.r \
-A ${BASEDIR}/05_tables/bracken.txt \
-O ${BASEDIR}/06_chempar/ \
-r '' -g ${BASEDIR}/06_chempar/sugar_aroma_netcorr.pdf

Rscript ${FUNCDIR}/09_correl_network.r \
-A ${BASEDIR}/05_tables/bracken.txt \
-O ${BASEDIR}/06_chempar/ \
-r 'Streptococcus thermophilus' -g ${BASEDIR}/06_chempar/sugar_aroma_nothermo_netcorr.pdf

#In this version we remove all streptococci that are not thermophilus 
#This is because they may be missclassified thermo 
#Even if htye aren't, they behave similar to thermo for most properties and do not add information
Rscript ${FUNCDIR}/09_correl_network.r \
-A ${BASEDIR}/05_tables/bracken.txt \
-O ${BASEDIR}/06_chempar/ \
-r 'Streptococcus agalactiae,Streptococcus anginosus,Streptococcus dysgalactiae,Streptococcus infantarius,Streptococcus iniae,Streptococcus marmotae,Streptococcus mutans,Streptococcus parauberis,Streptococcus pneumoniae,Streptococcus pyogenes,Streptococcus sobrinus,Streptococcus suis' \
-g ${BASEDIR}/06_chempar/sugar_aroma_onlythermo_netcorr.pdf

#Using silva at the genus level
Rscript ${FUNCDIR}/09_correl_network.r \
-A ${BASEDIR}/05_tables_silva/genus_bracken.txt \
-O ${BASEDIR}/06_chempar_silva/ \
-r '' -g ${BASEDIR}/06_chempar_silva/sugar_aroma_netcorr.pdf

Rscript ${FUNCDIR}/09_correl_network.r \
-A ${BASEDIR}/05_tables_silva/genus_bracken.txt \
-O ${BASEDIR}/06_chempar_silva/ \
-r 'Streptococcus' -g ${BASEDIR}/06_chempar_silva/sugar_aroma_nostrepto_netcorr.pdf

################################################################################################
#
#  _____                                     _       _   _             
# |  __ \                                   | |     | | (_)            
# | |__) |_ _ _ __    ___ ___  _ __ _ __ ___| | __ _| |_ _  ___  _ __  
# |  ___/ _` | '__|  / __/ _ \| '__| '__/ _ \ |/ _` | __| |/ _ \| '_ \ 
# | |  | (_| | |    | (_| (_) | |  | | |  __/ | (_| | |_| | (_) | | | |
# |_|   \__,_|_|     \___\___/|_|  |_|  \___|_|\__,_|\__|_|\___/|_| |_|
#                                                                      
################################################################################################
                                                                      
Rscript ${FUNCDIR}/10_par_cor.r
