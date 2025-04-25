# nhp-fungi
Repository for scripts related to Van Syoc et al, under review

## ITS amplicon data processing (`amplicon_processing`)  

The user can first install a conda environment with the necessary softwares to run this module using `setup_its_env.sh`. This module is designed to be fully reproducible and at least partially executable for a different project based on ITS bioinformatics for SRA accession numbers.  

Two main bash scripts for bioinformatics:  
`its-pipeline-pairedandsingle.sh` takes a Bioproject accession number as input and performs a series of steps to retrieve the data, perform quality control and primer removal (if a primer sequence is provided in a fasta file), and runs OTU clustering in VSEARCH with taxonomy in SINTAX in a user-supplied database fasta (must be formatted for sintax). This is built to run on paired or single end data. Right now very little customization is allowed.  

`run-vsearch-allsamples.sh` takes an array of Bioproject accessions and combines the quality filtered files together into a one folder, then runs the VSEARCH OTU clustering on all sequencing files together. This is run after `its-pipeline-pairedandsingle.sh` since it will look for folders with the same naming convention.  

`get-sample-metadata.R` wrangles metadata for each study and combines to one object. 

`setup_ITS2.R` wrangles VSEARCH output and performs rarefaction and quality filtering to produce final OTU table.