# nhp-fungi
Repository for scripts related to Van Syoc et al, under review

#"An agglomeration of Frankenscripts in pursuit of scientific novelty..."

## ITS amplicon data processing (`amplicon_processing`)  

The user can first install a conda environment with the necessary softwares to run this module using `setup_its_env.sh`. This module is designed to be fully reproducible and at least partially executable for a different project based on ITS bioinformatics for SRA accession numbers.  

Two main bash scripts for bioinformatics:  
`its-pipeline-pairedandsingle.sh` takes a Bioproject accession number as input and performs a series of steps to retrieve the data, perform quality control and primer removal (if a primer sequence is provided in a fasta file), and runs OTU clustering in VSEARCH with taxonomy in SINTAX in a user-supplied database fasta (must be formatted for sintax). This is built to run on paired or single end data. Right now very little customization is allowed.  

`run-vsearch-allsamples.sh` takes an array of Bioproject accessions and combines the quality filtered files together into a one folder, then runs the VSEARCH OTU clustering on all sequencing files together. This is run after `its-pipeline-pairedandsingle.sh` since it will look for folders with the same naming convention.  

`get-sample-metadata.R` wrangles metadata for each study and combines to one object. 

`setup_ITS2.R` wrangles VSEARCH output and performs rarefaction and quality filtering to produce final OTU table.

## Topological congruency testing for phylosymbiosis (folder `phylosymbiosis`)  

This is a series of steps to wrangle a phyloseq object into QIIME2 for the beta rarefaction pipeline, then manipulates the resulting newick tree to create a mycobiome dendrogram, and finally conducts statistical testing with the TreeCMP software as originally implemented in [Brooks 2016](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2000225).  

`phylo_to_qiime.R` wrangles a phyloseq object into a biom table for input to QIIME2.  

`qiime2_betadiv.sh` conducts beta diversity rarefaction with QIIME2. NOTE that the output of this (the 'qzv' file) has to be manually input into QIIME2 viewer to download the Newick tree (or you can unzip the file and go into the QIIME file structure, but this is faster IMO). The trees are then manually rooted with an outgroup in Figtree (no code associated - this is a weird step but Figtree roots in a much more intuitive way compared to the ape R Package, and this has a measurable difference in the final topology). These dendrograms are placed in the folder `data/rooted_dendrograms/`.  

`prep_hosttree.R` Reads in the primate phylogeny accessed from Time Tree and does some basic wrangling like renaming cogener branchs and subsetting to the host species of our dataset.

`prep_qiime2_topocong.R` The fungal dendrograms are rotated on the internal branches. TreeCMP implements topological congruency testing that is fairly 'brittle', meaning that it doesn't automatically rotate branches to try and find the best fit, so this is done by hand. 

`topocong_on_qiimedendrogram.sh` Sets up topological congruency testing for each of the dendrograms.  

`run_topocong_script.R` This is the R script called in `topocong_on_qiimedendrogram.sh` that does the heavy lifting. This is a bit of a Russian Dolls situation since it calls several older scripts that call other scripts, etc... In a world with more time, these can all be cleaned and concatenated into one Python or Nextflow file! It's also a bit of a Frankenscript that moves back and forth between R language and bash (perhaps better suited for a Jupyter notebook?). This is a fully executable script to conduct topological congruency testing. It uses the exact same setup as Brooks 2016 but has been converted to run in R. Some scripts are original to the Brooks 2016 paper. This calls `run_randomtrees.sh`  to generate random trees for a null distribution (note that `run_randomtrees.sh` then calls `tree_random.py`, Andy's original Python2 script from Brooks 2016 for generating random trees - to run this, install python2 as outlined in `create-py2-env.sh`). 

## Bray-Curtis distances (folder `pairwise_braycurtis`)  

`pairwise_distances.R` Performs statistical testing for intraspecific and interspecific Bray-Curtis pairwise distances. Other ancillary statistical tests are included and several supplementary tables are generated. These objects are then sourced into `figures/Figure1.R` to create the figures.  

## Human-specific mycobiomes (folder `human_mycobiome`)  

`differential_relabundance.R` calculates some ancillary statistics and runs linear models for differential relative abundance testing.  

## Cophylogeny (folder `cophylogeny`)  

This is another collection of scripts that go back and forth between R and bash. The main script `cophylogeny.R` is designed to call each other script in the folder `cophylogeny_scripts` as it runs. These scripts are not fully executable because the fasta file of sequences is big for GitHub and some of the hominid metadata was shared directly by study authors. The intermediate files are also clunky; multiple trees and text files are written during tree generation for each fungal genus and subsequent wrangling to get it ready for cophylogeny testing. The basic framework can be applied on a phyloseq object with corresponding metadata (like OTU cluster centroids) from VSEARCH.  

`results_parafit_paco.R` wrangles the results and calculates multiple comparisons corrections.  

## Sequence divergence times (folder `seqdivergence`)  

`sequence_divergence.R` uses the framework of Moeller 2016 to calibrate a fungal molecular clock using nucleotide sequence divergence, then estimates hominid speciation based on those clocks. The intermediate fasta files from cophylogenetic testing are not uploaded to GitHub (large), but an image file can be uploaded midway through the script to make the clock estimates executable.



## Figures (folder `figures`)  

These source the analysis scripts and a helper script that sets the colors for each host species to create each figure and supplementary figure. Note that Figure 3 is the only one with paneling and text labels done in Illustrator; all others are top-to-bottom hardcoded in R.

