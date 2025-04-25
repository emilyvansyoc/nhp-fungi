#!/bin/bash

## compile all samples and run VSEARCH clustering on the global dataset
# EVS 3/10/2025

# activate environment
eval "$(conda shell.bash hook)"
conda activate its_pipeline_env
# Trap Ctrl+C and exit the script completely -> makes cmd+c quit the entire script 'gracefully'
trap "echo 'Script interrupted! Exiting...'; exit 1" SIGINT
set -e # quit immediately if any command fails

# get into the working directory
DIR=$1
cd $DIR

## OTU cutoff to use
OTU_CUTOFF=0.99

# add flags to one just one at a time optionally
skip18s="true"
skipITS1="true"
skipITS2="false"
## print progress message to track which samples are being processed
if [ $skip18s = "false" ]; then
    echo "##### Processing 18S samples #####"
fi
if [ $skipITS1 = "false" ]; then
    echo "##### Processing ITS1 samples #####"
fi
if [ $skipITS2 = "false" ]; then
    echo "##### Processing ITS2 samples #####"
fi

# make a working directory and subdirectory for all files
WORKDIR=$DIR/vsearch_all_primate
mkdir -p $WORKDIR
DIR18S=$WORKDIR/all18S_99
DIRITS1=$WORKDIR/allITS1_99
DIRITS2=$WORKDIR/allITS2_ITS3ITS4_99
mkdir -p $DIR18S
mkdir -p $DIRITS1
mkdir -p $DIRITS2
mkdir -p $DIR18S/allreads/
mkdir -p $DIRITS1/allreads/
mkdir -p $DIRITS2/allreads/

# set location to UNITE database
DB=$2
# master list of IDS
IDS=$DIR/data/all_primate_runids.txt

# make arrays of accessions for ITS and for 18S
acc18=(
    "PRJEB32407"
    "PRJNA694983_18S"
    "PRJNA913825"
    )
accITS1=(
    "PRJNA701719"
    "PRJNA447371"
    "PRJNA593927"
)
#accITS2=(
#    "PRJNA634617"
#    "PRJEB39443"
#    "PRJNA694983_ITS2"
#    "PRJNA686661"
#    "PRJEB37770_Zenodo"
#)
## ITS2 accessions with primers ITS3 and ITS4
accITS2=(
	"PRJEB37770_Zenodo"
	"PRJNA686661"
	"PRJEB39443"
)

# search each directory for either cutadapt or filtN files that match $IDS
# Find all PRJ* directories and search for matching files

        ####18S
if [ $skip18s = "false" ]; then
        for acc in "${acc18[@]}"; do     
            # Check for cutadapt directory first
            if [ -d "${DIR}/${acc}/rawreads/cutadapt" ]; then
                SEARCHDIR="${DIR}/${acc}/rawreads/cutadapt"
                # print progress message
                echo "Processing ${DIR}/${acc}/rawreads/cutadapt"
            # Fall back to filtN if cutadapt doesn't exist
            elif [ -d "${DIR}/${acc}/rawreads/filtN" ]; then
                SEARCHDIR="${DIR}/${acc}/rawreads/filtN"
                # print progress message
                echo "Processing ${DIR}/${acc}/rawreads/filtN"
            else
                echo "Warning: Neither cutadapt nor filtN directory found in ${DIR}/${acc}"
                continue
            fi

            # print the number of files in $SEARCHDIR that have a partial match to $IDS
            echo "Number of files that have a partial match to IDs file: $(ls $SEARCHDIR | grep -f $IDS | wc -l)"

            # get runids file
            runids_file=$(find "$DIR" -type f -name "runids*_${acc}.txt" | head -n 1)

            # Read IDs from file and search for matching fastq files
            while IFS= read -r id; do
                # Find and copy matching files to INDIR
                find "$SEARCHDIR" -type f \( -name "${id}_1.fastq" -o -name "${id}.fastq" \) -exec cp {} "$DIR18S/allreads/" \; #2>/dev/null
            done < "$runids_file"
        done
fi
        ####ITS1
if [ $skipITS1 = "false" ]; then
        for acc in "${accITS1[@]}"; do     
            # Check for cutadapt directory first
            if [ -d "${DIR}/${acc}/rawreads/cutadapt" ]; then
                SEARCHDIR="${DIR}/${acc}/rawreads/cutadapt"
                # print progress message
                echo "Processing ${DIR}/${acc}/rawreads/cutadapt"
            # Fall back to filtN if cutadapt doesn't exist
            elif [ -d "${DIR}/${acc}/rawreads/filtN" ]; then
                SEARCHDIR="${DIR}/${acc}/rawreads/filtN"
                # print progress message
                echo "Processing ${DIR}/${acc}/rawreads/filtN"
            else
                echo "Warning: Neither cutadapt nor filtN directory found in ${DIR}/${acc}"
                continue
            fi

            # print the number of files in $SEARCHDIR that have a partial match to $IDS
            echo "Number of files that have a partial match to IDs file: $(ls $SEARCHDIR | grep -f $IDS | wc -l)"

            # get runids file
            runids_file=$(find "$DIR" -type f -name "runids*_${acc}.txt" | head -n 1)

            # Read IDs from file and search for matching fastq files
            while IFS= read -r id; do
                # Find and copy matching files to INDIR
                find "$SEARCHDIR" -type f \( -name "${id}_1.fastq" -o -name "${id}.fastq" \) -exec cp {} "$DIRITS1/allreads/" \; #2>/dev/null
            done < "$runids_file"
        done
fi
        ####ITS2
if [ $skipITS2 = "false" ]; then
        for acc in "${accITS2[@]}"; do     
            # Check for cutadapt directory first
            if [ -d "${DIR}/${acc}/rawreads/cutadapt" ]; then
                SEARCHDIR="${DIR}/${acc}/rawreads/cutadapt"
                # print progress message
                echo "Processing ${DIR}/${acc}/rawreads/cutadapt"
            # Fall back to filtN if cutadapt doesn't exist
            elif [ -d "${DIR}/${acc}/rawreads/filtN" ]; then
                SEARCHDIR="${DIR}/${acc}/rawreads/filtN"
                # print progress message
                echo "Processing ${DIR}/${acc}/rawreads/filtN"
            else
                echo "Warning: Neither cutadapt nor filtN directory found in ${DIR}/${acc}"
                continue
            fi

            # print the number of files in $SEARCHDIR that have a partial match to $IDS
            echo "Number of files that have a partial match to IDs file: $(ls $SEARCHDIR | grep -f $IDS | wc -l)"

            # get runids file
            runids_file=$(find "$DIR" -type f -name "runids*_${acc}.txt" | head -n 1)

            # Read IDs from file and search for matching fastq files
            while IFS= read -r id; do
                # Find and copy matching files to INDIR
                find "$SEARCHDIR" -type f \( -name "${id}_1.fastq" -o -name "${id}.fastq" \) -exec cp {} "$DIRITS2/allreads/" \; #2>/dev/null
            done < "$runids_file"
        done
fi

##### ITS1 ######
if [ $skipITS1 = "false" ]; then 
	# Process ITS samples
	WORKDIR="$DIRITS1"
	# write output to log file
	LOG="$WORKDIR/vsearch.log"
	> $LOG

	# print progress message
	echo "Running VSEARCH clustering on ITS1 samples, see $LOG for details..."

	# location for quality filtered reads
	QUAL="$WORKDIR/qualfilt_vsearch/"
	mkdir -p $QUAL

	# dereplicated reads
	DEREP=$WORKDIR/derep_vsearch/
	mkdir -p $DEREP

	# loop through each file in $INDIR and do VSEARCH quality filtering


	echo "Quality filtering ITS1 samples..."
	for file in $WORKDIR/allreads/*.fastq; do
		# get basename without extension
		base=$(basename "$file" .fastq)
		vsearch --quiet --fastx_filter "$file" --fastq_minlen 20 --fastq_truncqual 10 --fastaout "$QUAL/${base}.fasta" --log /dev/stdout >> $LOG 2>&1
	done 

	# run stats
	seqkit stat $QUAL/*.fasta -T -j 10 -o $WORKDIR/stats_qualfilt_vsearch.txt

	# run sample-wise dereplication
	echo "Dereplicating ITS1 samples sample-wise..."
	for file in $QUAL/*.fasta; do
		base=$(basename "$file" .fasta)
		vsearch --quiet --derep_fulllength "$file" --output "$DEREP/${base}.derep.fasta" --sizeout --relabel "${base}.derep" --log /dev/stdout >> $LOG 2>&1
	done 


	# merge all samples
	cat $DEREP/*.derep.fasta > $WORKDIR/derep.fasta
	#seqkit scat -I fasta -i fasta -j 4 -f $DEREP > $WORKDIR/derep_all.fasta

	# dereplicate whole dataset
	echo "Dereplicating whole ITS1 dataset..."
	vsearch --derep_fulllength $WORKDIR/derep.fasta --minuniquesize 2 --sizein --sizeout --output $WORKDIR/all.derep.fasta >> $LOG 2>&1

	# cluster
	echo "Clustering ITS1 OTUs at $OTU_CUTOFF identity..."
	vsearch --cluster_size $WORKDIR/all.derep.fasta --threads 10 --id $OTU_CUTOFF --qmask none --sizein --sizeout --centroids $WORKDIR/all.centroids.fasta >> $LOG 2>&1

	# sort and remove singletons
	echo "Sorting and removing ITS1 singletons..."
	vsearch --sortbysize $WORKDIR/all.centroids.fasta --sizein --sizeout --minsize 2 --output $WORKDIR/all.sorted.fasta >> $LOG 2>&1

	# denovo chimeras
	echo "Removing ITS1 denovo chimeras..."
	vsearch --quiet --uchime_denovo $WORKDIR/all.sorted.fasta --sizein --sizeout --qmask none --nonchimeras $WORKDIR/all.denovo.nonchimeras.fasta  >> $LOG 2>&1

	# relabel OTUs
	echo "Relabeling ITS1 OTUs..."
	vsearch --fastx_filter $WORKDIR/all.denovo.nonchimeras.fasta --sizein --sizeout --relabel OTU_ --fastaout $WORKDIR/all.otus.fasta  >> $LOG 2>&1

	# assign OTUs to sequences
	echo "Assigning ITS1 OTUs to sequences..."
	vsearch --usearch_global $WORKDIR/all.derep.fasta --threads 10 --db $WORKDIR/all.otus.fasta --id $OTU_CUTOFF --strand plus --sizein --sizeout --qmask none --dbmask none --otutabout $WORKDIR/all.otutab.txt  >> $LOG 2>&1

	# run SINTAX
	echo "Running SINTAX on ITS1 $DB..."
	vsearch --sintax $WORKDIR/all.otus.fasta --db $DB --tabbedout $WORKDIR/all.sintax50.txt --sintax_cutoff .50 --threads 10  >> $LOG 2>&1

fi


###### 18S ######
if [ $skip18s = "false" ]; then
	# Process 18S samples
	WORKDIR="$DIR18S"
	# write output to log file
	LOG="$WORKDIR/vsearch.log"
	> $LOG

	# print progress message
	echo "Running VSEARCH clustering on 18S samples, see $LOG for details..."

	# location for quality filtered reads
	QUAL="$WORKDIR/qualfilt_vsearch/"
	mkdir -p $QUAL

	# dereplicated reads
	DEREP=$WORKDIR/derep_vsearch/
	mkdir -p $DEREP

	# loop through each file in $INDIR and do VSEARCH quality filtering

 
		echo "Quality filtering 18S samples..."
		for file in $WORKDIR/allreads/*.fastq; do
		# get basename without extension
		base=$(basename "$file" .fastq)
		vsearch --quiet --fastx_filter "$file" --fastq_minlen 20 --fastq_truncqual 10 --fastaout "$QUAL/${base}.fasta" --log /dev/stdout >> $LOG 2>&1
		done 

		# run stats
		seqkit stat $QUAL/*.fasta -T -j 10 -o $WORKDIR/stats_qualfilt_vsearch.txt

		# run sample-wise dereplication
		echo "Dereplicating 18S samples sample-wise..."
		for file in $QUAL/*.fasta; do
			base=$(basename "$file" .fasta)
			vsearch --quiet --derep_fulllength "$file" --output "$DEREP/${base}.derep.fasta" --sizeout --relabel "${base}.derep" --log /dev/stdout >> $LOG 2>&1
		done 
		




	# merge all samples
	cat $DEREP/*.derep.fasta > $WORKDIR/derep.fasta

	# dereplicate whole dataset
	echo "Dereplicating whole 18S dataset..."
	vsearch --derep_fulllength $WORKDIR/derep.fasta --minuniquesize 2 --sizein --sizeout --output $WORKDIR/all.derep.fasta >> $LOG 2>&1

	# cluster
	echo "Clustering 18S OTUs at $OTU_CUTOFF identity..."
	vsearch --cluster_size $WORKDIR/all.derep.fasta --threads 10 --id $OTU_CUTOFF --qmask none --sizein --sizeout --centroids $WORKDIR/all.centroids.fasta >> $LOG 2>&1

	# sort and remove singletons
	echo "Sorting and removing 18S singletons..."
	vsearch --sortbysize $WORKDIR/all.centroids.fasta --sizein --sizeout --minsize 2 --output $WORKDIR/all.sorted.fasta >> $LOG 2>&1

	# denovo chimeras
	echo "Removing 18S denovo chimeras..."
	vsearch --uchime_denovo $WORKDIR/all.sorted.fasta --sizein --sizeout --qmask none --nonchimeras $WORKDIR/all.denovo.nonchimeras.fasta  >> $LOG 2>&1

	# relabel OTUs
	echo "Relabeling 18S OTUs..."
	vsearch --fastx_filter $WORKDIR/all.denovo.nonchimeras.fasta --sizein --sizeout --relabel OTU_ --fastaout $WORKDIR/all.otus.fasta  >> $LOG 2>&1

	# assign OTUs to sequences
	echo "Assigning 18S OTUs to sequences..."
	vsearch  --usearch_global $WORKDIR/all.derep.fasta --threads 10 --db $WORKDIR/all.otus.fasta --id $OTU_CUTOFF --strand plus --sizein --sizeout --qmask none --dbmask none --otutabout $WORKDIR/all.otutab.txt  >> $LOG 2>&1

	# run SINTAX
	echo "Running SINTAX on 18S $DB..."
	vsearch --sintax $WORKDIR/all.otus.fasta --db $DB --tabbedout $WORKDIR/all.sintax50.txt --sintax_cutoff .50  >> $LOG 2>&1

fi

###### ITS2 ######
if [ $skipITS2 = "false" ]; then    
	# process ITS2 samples
	WORKDIR="$DIRITS2"
	# write output to log file
	LOG="$WORKDIR/vsearch.log"
	> $LOG

	# print progress message
	echo "Running VSEARCH clustering on ITS2 samples, see $LOG for details..."

	# location for quality filtered reads
	QUAL="$WORKDIR/qualfilt_vsearch/"
	mkdir -p $QUAL

	# dereplicated reads
	DEREP=$WORKDIR/derep_vsearch/
	mkdir -p $DEREP

	# loop through each file in $INDIR and do VSEARCH quality filtering

		echo "Quality filtering ITS2 samples..."
		for file in $WORKDIR/allreads/*.fastq; do
		# get basename without extension
		base=$(basename "$file" .fastq)
		vsearch --fastx_filter "$file" --fastq_minlen 20 --fastq_truncqual 10 --fastaout "$QUAL/${base}.fasta" --log /dev/stdout >> $LOG 2>&1
		done 

		# run stats
		seqkit stat $QUAL/*.fasta -T -j 10 -o $WORKDIR/stats_qualfilt_vsearch.txt

		# run sample-wise dereplication
		echo "Dereplicating ITS2 samples sample-wise..."
		for file in $QUAL/*.fasta; do
			base=$(basename "$file" .fasta)
			vsearch --derep_fulllength "$file" --output "$DEREP/${base}.derep.fasta" --sizeout --relabel "${base}.derep" --log /dev/stdout >> $LOG 2>&1
		done 

		


	# merge all samples
	cat $DEREP/*.derep.fasta > $WORKDIR/derep.fasta

	# dereplicate whole dataset
	echo "Dereplicating whole ITS2 dataset..."
	vsearch --derep_fulllength $WORKDIR/derep.fasta --minuniquesize 2 --sizein --sizeout --output $WORKDIR/all.derep.fasta >> $LOG 2>&1

	# cluster
	echo "Clustering ITS2 OTUs at $OTU_CUTOFF identity..."
	vsearch --cluster_size $WORKDIR/all.derep.fasta --threads 10 --id $OTU_CUTOFF --qmask none --sizein --sizeout --centroids $WORKDIR/all.centroids.fasta >> $LOG 2>&1

	# sort and remove singletons
	echo "Sorting and removing ITS2 singletons..."
	vsearch --sortbysize $WORKDIR/all.centroids.fasta --sizein --sizeout --minsize 2 --output $WORKDIR/all.sorted.fasta >> $LOG 2>&1

	# denovo chimeras
	echo "Removing ITS2 denovo chimeras..."
	vsearch --uchime_denovo $WORKDIR/all.sorted.fasta --sizein --sizeout --qmask none --nonchimeras $WORKDIR/all.denovo.nonchimeras.fasta  >> $LOG 2>&1

	# relabel OTUs
	echo "Relabeling ITS2 OTUs..."
	vsearch --fastx_filter $WORKDIR/all.denovo.nonchimeras.fasta --sizein --sizeout --relabel OTU_ --fastaout $WORKDIR/all.otus.fasta  >> $LOG 2>&1

	# assign OTUs to sequences
	echo "Assigning ITS2 OTUs to sequences..."
	vsearch --usearch_global $WORKDIR/all.derep.fasta --threads 10 --db $WORKDIR/all.otus.fasta --id $OTU_CUTOFF --strand plus --sizein --sizeout --qmask none --dbmask none --uc $WORKDIR/usearch.uc--otutabout $WORKDIR/all.otutab.txt  >> $LOG 2>&1

	# run SINTAX
	echo "Running SINTAX on ITS2 $DB..."
	vsearch --sintax $WORKDIR/all.otus.fasta --db $DB --tabbedout $WORKDIR/all.sintax50.txt --sintax_cutoff .50  --threads 10 >> $LOG 2>&1

fi

# exit without errors
exit 0


