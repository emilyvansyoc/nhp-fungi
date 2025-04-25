#!/bin/bash

### ----- CHECK CONDA ENVIRONMENT -----
# Check if the required conda environment is active
if [[ "$CONDA_DEFAULT_ENV" != "its_pipeline_env" ]]; then
    echo "Error: Please activate the 'its_pipeline_env' conda environment first."
    echo "Run: conda activate its_pipeline_env"
    echo "If you haven't set up the environment yet, run:"
    echo "bash $(dirname "$0")/setup_its_env.sh"
    exit 1
fi

## master pipeline for pulling ITS data from SRA and creating OTUs/ASVs
# EVS 2/2025

## if needed, activate the conda environment within the bash script:
#eval "$(conda shell.bash hook)"
#conda activate its_pipeline_env
#Trap Ctrl+C and exit the script completely -> makes cmd+c quit the entire script 'gracefully'
trap "echo 'Script interrupted! Exiting...'; exit 1" SIGINT
set -e # quit immediately if any command fails

### ----- VARIABLES ----

# accession number
ACC=$1

# print progress message
echo "Running pipeline for accession: $ACC"

# input directory
DIR=$2
mkdir -p $DIR/$ACC

# location to SINTAX UNITE DB
DB=$3

### ----- get runinfo ----
# get runinfo
esearch -db sra -query $ACC | efetch -format runinfo > $DIR/runinfo_$ACC.csv

# get runids - search first for SRR, then try ERR
PREFIXES=("SRR" "ERR")
success=false # flag to error catch

# start loop for prefixes
for PREF in "${PREFIXES[@]}"; do 
    echo "Looking for runs that start with: $PREF"
    if cat $DIR/runinfo_$ACC.csv | cut -f 1 -d ',' | grep "$PREF" > $DIR/runids_$ACC.txt; then
        success=true
        break # exit loop if successfull
    fi

    # if that fails, print this
    echo "Prefix $PREF failed. Trying next..."
done

# if all attempts failed, exit with error message
if ! $success; then
    echo "Error: Cannot find runs that start with SRR or ERR." >&2  # Print to stderr
    exit 1  # Exit with a failure status
fi

# print message
echo "Number of samples detected: "
cat $DIR/runids_$ACC.txt | wc -l

## ---- fastq-dump ----

# make directory
RAW=$DIR/$ACC/rawreads
mkdir -p $RAW

# test if samples are paired or single
if grep -q "PAIRED" $DIR/runinfo_$ACC.csv; then
    # check to make sure that there is not both paired and single end samples, fail if so       
    if grep -q "SINGLE" $DIR/runinfo_$ACC.csv; then
        echo "Error: Both paired and single end samples detected"
        exit 1
    fi  
    echo "Detecting paired end samples"
    # Set boolean for paired-end reads
    is_paired=true
else
    if grep -q "SINGLE" $DIR/runinfo_$ACC.csv; then
        echo "Detecting single end samples"
        is_paired=false
    else
        echo "Error: No paired or single end samples detected"
        exit 1
    fi
fi

# run fastq-dump on all samples
# test to see if all files are already downloaded, skip if so
# run for paired or single end
if [[ "$is_paired" == "true" ]]; then
    cat $DIR/runids_$ACC.txt | \
        awk -v raw="$RAW" '{if (system("test -f " raw"/"$1"_1.fastq") || system("test -f " raw"/"$1"_2.fastq")) print $1}' | \
        parallel -j 10 --will-cite --halt now,fail=1 "fastq-dump --split-files --outdir $RAW {} > $DIR/$ACC/fqdump.log 2>&1"
else
    cat $DIR/runids_$ACC.txt | \
        awk -v raw="$RAW" '{if (system("test -f " raw"/"$1".fastq")) print $1}' | \
        parallel -j 10 --will-cite --halt now,fail=1 "fastq-dump --outdir $RAW {} > $DIR/$ACC/fqdump.log 2>&1"
fi

# # error catch: if no reverse reads found, print warning
# if [ -e "$RAW/*_2.fastq" ]; then
#     echo "Reverse reads not found: are you sure this is paired end?"
# else 
#     echo "Number of paired end reads downloaded:" 
#     ls $RAW | grep fastq | wc -l
# fi

# print stats for the reads
seqkit stat $RAW/*.fastq --quiet -T -j 10 -o $DIR/$ACC/stats_raw_$ACC.txt

### ---- fastqc and multiqc ----

# print progress message
echo "Running FastQC..."

# make directory
FQC=$DIR/$ACC/fastqc_files/
mkdir -p $FQC

# run for paired or single end
if [[ "$is_paired" == "true" ]]; then
    mkdir -p $FQC/forward/
    mkdir -p $FQC/reverse/
    
    # Check if FastQC files exist for all samples
    all_files_exist=true
    while read -r runid; do
        if [ ! -f "$FQC/forward/${runid}_1_fastqc.html" ] || \
           [ ! -f "$FQC/forward/${runid}_1_fastqc.zip" ] || \
           [ ! -f "$FQC/reverse/${runid}_2_fastqc.html" ] || \
           [ ! -f "$FQC/reverse/${runid}_2_fastqc.zip" ]; then
            all_files_exist=false
            break
        fi
    done < "$DIR/runids_$ACC.txt"

    if [[ "$all_files_exist" == "true" ]]; then
        echo "FastQC files already exist for all paired-end samples, skipping..."
    else
        echo "Running FastQC on paired-end samples..."
        cat $DIR/runids_$ACC.txt | \
            parallel -j 10 --will-cite --halt now,fail=1 "
                fastqc --quiet --threads 10 --outdir $FQC/forward/ $RAW/{}_1.fastq
                fastqc --quiet --threads 10 --outdir $FQC/reverse/ $RAW/{}_2.fastq"
    fi

    # Run MultiQC if reports don't exist
    if [ ! -f "$DIR/$ACC/multiqc_fwd_$ACC.html" ] || [ ! -f "$DIR/$ACC/multiqc_rev_$ACC.html" ]; then
        echo "Running MultiQC..."
        multiqc -q $FQC/forward/* -n multiqc_fwd_$ACC -o $DIR/$ACC
        multiqc -q $FQC/reverse/* -n multiqc_rev_$ACC -o $DIR/$ACC
    else
        echo "MultiQC reports already exist, skipping..."
    fi
else
    # Check if FastQC files exist for all single-end samples
    all_files_exist=true
    while read -r runid; do
        if [ ! -f "$FQC/${runid}_fastqc.html" ] || [ ! -f "$FQC/${runid}_fastqc.zip" ]; then
            all_files_exist=false
            break
        fi
    done < "$DIR/runids_$ACC.txt"

    if [[ "$all_files_exist" == "true" ]]; then
        echo "FastQC files already exist for all single-end samples, skipping..."
    else
        echo "Running FastQC on single-end samples..."
        cat $DIR/runids_$ACC.txt | \
            parallel -j 10 --will-cite --halt now,fail=1 "
                fastqc --quiet --threads 10 --outdir $FQC/ $RAW/{}.fastq"
    fi

    # Run MultiQC if report doesn't exist
    if [ ! -f "$DIR/$ACC/multiqc_$ACC.html" ]; then
        echo "Running MultiQC..."
        multiqc -q $FQC/* -n multiqc_$ACC -o $DIR/$ACC
    else
        echo "MultiQC report already exists, skipping..."
    fi
fi

### ---- primer removal ----

# if a primer file is found in the directory, run primer removal
# otherwise skip to next step
# Set cutadapt directory path
CUTAD="$RAW/cutadapt"

# Check if cutadapt directory exists and has files
if [ -d "$CUTAD" ] && [ "$(ls -A $CUTAD)" ]; then
    echo "Cutadapt files already exist, skipping primer removal..."
else
    # run for paired or single end
    if [[ "$is_paired" == "true" ]]; then 
        echo "Running cutadapt for paired-end samples"
        Rscript "amplicon_processing/dada2-cutadapt.R" "$RAW" "$DIR/primers_$ACC.fasta"
    else
        echo "Running cutadapt for single-end samples"
        Rscript "amplicon_processing/dada2-cutadapt_single.R" "$RAW" "$DIR/primers_$ACC.fasta"
    fi
fi

# run stats on filtN reads
seqkit stat $RAW/filtN/*.fastq --quiet -T -j 10 -o $DIR/$ACC/stats_filtN_$ACC.txt

## ---- VSEARCH for dereplication and OTU generation ----

### this creates several subdirectories

# write output to log file
LOG="$DIR/$ACC/vsearch.log"
> $LOG

# print progress message
echo "Running VSEARCH for dereplication and OTU generation, see $LOG for details..."

# location for quality filtered reads
QUAL="$DIR/$ACC/qualfilt/"
mkdir -p $QUAL

# merged reads
MERGE=$DIR/$ACC/merge/
mkdir -p $MERGE

# dereplicated reads
DEREP=$DIR/$ACC/derep/
mkdir -p $DEREP

## run the code VSEARCH code ######

## error catch: if primers were not removed with cutadapt, use filtN files instead
# skip this step if reads are single end
if [[ "$is_paired" == "true" ]]; then
    testfile=($(ls $RAW | head -1 | cut -f 1 -d _ ))
    # if detect paired end, check for primer removed reads
    if [ -e "$CUTAD/${testfile}_1.fastq" ]; then
        echo "VSEARCH is merging the primer-removed reads"
    # merge paired-end reads
    {
        echo "-------------------"
        echo "$(date): Starting merge of primer-removed reads"
        cat $DIR/runids_$ACC.txt | \
            parallel -j 10 --will-cite --halt now,fail=1 "
                {
                    echo \"Processing {}\";
                    vsearch --quiet --fastq_mergepairs $CUTAD/{}_1.fastq --threads 10 --reverse $CUTAD/{}_2.fastq --fastq_minovlen 30 --fastq_allowmergestagger --fastqout $MERGE/{}.fastq --log /dev/stdout;
                } 2>&1"
        } >> $LOG 2>&1

    # print stats for the reads
    seqkit stat $CUTAD/*.fastq --quiet  -T -j 10 -o $DIR/$ACC/stats_cutadapt_$ACC.txt

    else
    echo "VSEARCH merging the reads with Ns removed"
    # merge paired-end reads
    {
        echo "-------------------"
        echo "$(date): Starting merge of N-filtered reads"
        cat $DIR/runids_$ACC.txt | \
            parallel -j 10 --will-cite --halt now,fail=1 "
                {
                    echo \"Processing {}\";
                    vsearch --quiet --fastq_mergepairs $RAW/filtN/{}_1.fastq --threads 10 --reverse $RAW/filtN/{}_2.fastq --fastq_minovlen 30 --fastq_allowmergestagger --fastqout $MERGE/{}.fastq --log /dev/stdout;
                } 2>&1"
        } >> $LOG 2>&1
    fi  
    # run stats
    seqkit stat $MERGE/*.fastq --quiet -T -j 10 -o $DIR/$ACC/stats_merged_$ACC.txt
else    
    echo "Skipping merge step as reads are single end"
fi



# quality filtered merged reads
# if reads are paired end, run quality filter on merged reads
# otherwise run quality filter on filtN reads
if [[ "$is_paired" == "true" ]]; then
    {
        echo "-------------------"
        echo "$(date): Starting quality filtering"
    cat $DIR/runids_$ACC.txt | \
        parallel -j 10 --will-cite --halt now,fail=1 "
            {
                echo \"Processing quality filtering for {}\";
                vsearch --quiet --fastx_filter $MERGE/{}.fastq --fastq_minlen 20 --fastq_truncqual 10 --fastaout $QUAL/{}.fasta --log /dev/stdout;
            } 2>&1"
    } >> $LOG 2>&1
   
else
    # run quality filter on filtN reads
    {
        echo "-------------------"
        echo "$(date): Starting quality filtering"
    cat $DIR/runids_$ACC.txt | \
        parallel -j 10 --will-cite --halt now,fail=1 "
            {
                echo \"Processing quality filtering for {}\";
                vsearch --quiet --fastx_filter $RAW/filtN/{}.fastq --fastq_minlen 20 --fastq_truncqual 10 --fastaout $QUAL/{}.fasta --log /dev/stdout;
            } 2>&1"
    } >> $LOG 2>&1
fi

# run stats
seqkit stat $QUAL/*.fasta --quiet -T -j 10 -o $DIR/$ACC/stats_qualfilt_$ACC.txt

# sample-wise dereplication
{
    echo "-------------------"
    echo "$(date): Starting per-sample dereplication"
    cat $DIR/runids_$ACC.txt | \
        parallel -j 10 --will-cite --halt now,fail=1 "
            {
                echo \"Processing dereplication for {}\";
                vsearch --quiet --derep_fulllength $QUAL/{}.fasta --output $DEREP/{}_derep.fasta --sizeout --relabel {}.derep --log /dev/stdout;
            } 2>&1"
} >> $LOG 2>&1

# merge all samples
cat $DEREP/*_derep.fasta > $DIR/$ACC/all.fasta 

# dereplicate whole datasets
{
    echo "-------------------"
    echo "$(date): Starting final dereplication of all sequences"
    echo "Minimum unique size: 2"
    vsearch --quiet --derep_fulllength $DIR/$ACC/all.fasta --minuniquesize 2 --sizein --sizeout --output $DIR/$ACC/all.derep.fasta --log /dev/stdout
} >> $LOG 2>&1

### NOTE: there is NO masking done in clustering or de novo steps; this can be changed

# cluster
{
    echo "-------------------"
    echo "$(date): Starting clustering at 97% identity"
    vsearch --quiet --cluster_size $DIR/$ACC/all.derep.fasta --threads 10 --id 0.97 --qmask none --sizein --sizeout --centroids $DIR/$ACC/centroids.fasta --log /dev/stdout
} >> $LOG 2>&1

# sort and remove singletons
{
    echo "-------------------"
    echo "$(date): Sorting by size and removing singletons"
    vsearch --quiet --sortbysize $DIR/$ACC/centroids.fasta --sizein --sizeout --minsize 2 --output $DIR/$ACC/sorted.fasta --log /dev/stdout
} >> $LOG 2>&1

# de novo chimeras
{
    echo "-------------------"
    echo "$(date): Detecting de novo chimeras"
    vsearch --quiet --uchime_denovo $DIR/$ACC/sorted.fasta --sizein --sizeout --qmask none --nonchimeras $DIR/$ACC/denovo.nonchimeras.fasta --log /dev/stdout
} >> $LOG 2>&1

# relabel OTUs
{
    echo "-------------------"
    echo "$(date): Relabeling OTUs"
    vsearch --quiet --fastx_filter $DIR/$ACC/denovo.nonchimeras.fasta --sizein --sizeout --relabel OTU_ --fastaout $DIR/$ACC/otus.fasta --log /dev/stdout
} >> $LOG 2>&1

# assign OTUs to sequences
{
    echo "-------------------"
    echo "$(date): Assigning sequences to OTUs"
    vsearch --quiet --usearch_global $DIR/$ACC/all.fasta --threads 10 --db $DIR/$ACC/otus.fasta --id 0.97 --strand plus --sizein --sizeout --qmask none --dbmask none --otutabout $DIR/$ACC/otutab.txt --log /dev/stdout
} >> $LOG 2>&1

# run SINTAX
# print which database is being used
echo "Running SINTAX using database: $DB"
{
    echo "-------------------"
    echo "$(date): Running taxonomic classification with SINTAX (50% cutoff)"
    vsearch --quiet --sintax $DIR/$ACC/otus.fasta --db $DB --tabbedout $DIR/$ACC/sintax50.txt --sintax_cutoff .50 --log /dev/stdout
} >> $LOG 2>&1

# exit without errors
exit 0