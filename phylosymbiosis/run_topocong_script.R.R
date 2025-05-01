### run topological congruency from the command line
# EVS 4/2025

# print start time
cat("\n########## Running Topological Congruency ########################  \n Started at: ", format(Sys.time(), "%a %b %d %X %Y"), "\n")

suppressPackageStartupMessages({
  suppressWarnings({
    if (!require(optparse)) stop("Package 'optparse' not found")
    if (!require(tidyverse)) stop("Package 'tidyverse' not found")
    if (!require(ape)) stop("Package 'ape' not found")
  })
})

# parse command line arguments
option_list <- list(
  # make_option(c("--help"),
  #   action = "store_true", default = FALSE,
  #  help = "Show this help message and exit"
  # ),
  make_option(c("-i", "--inputhost"),
              type = "character", default = "",
              help = "Path to host tree, expected in newick format", metavar = "character"
  ),
  make_option(c("-d", "--dendrogram"),
              type = "character", default = "",
              help = "Path to microbiome dendrogram, expected in newick format", metavar = "character"
  ),
  make_option(c("-c", "--treecmp"),
              type = "character", default = "/bin/TreeCmp_v2.0-b76/bin/treeCmp.jar",
              help = "Path to TreeCMP executable", metavar = "character"
  ),
  make_option(c("-m", "--makerandom"),
              type = "logical", default = "true",
              help = "Logical: make random trees? Necessary if path to random tree is not provided", metavar = "logical"
  ),
  make_option(c("-r", "--random"),
              type = "character", default = "",
              help = "Path to random trees for P value calculations, expected in newick format", metavar = "character"
  ),
  make_option(c("-t", "--tests"),
              type = "character", default = "mc rc ns tt mp mt co",
              help = "Tests to run, see TreeCmp documentation for options, expecting one space separated string", metavar = "character"
  ),
  make_option(c("-o", "--output"),
              type = "character", default = "",
              help = "Path to output DIRECTORY for writing multiple output files: FILES WILL OVERWRITE IF THEY ALREADY EXIST", metavar = "character"
  )
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))


# check if required arguments are provided
if (opt$inputhost == "" || opt$dendrogram == "" || opt$output == "") {
  cat("Error: All required arguments are not provided\n")
  cat("Usage: Rscript run_topocong_script.R -i <inputhost> -d <dendrogram> -c <treecmp> -o <output>\n")
  q(status = 1)
}

# check if random tree is provided if makerandom is true
if (opt$makerandom == "false" && opt$random == "") {
  cat("Error: Random tree is required if makerandom is false\n")
  cat("Usage: Rscript run_topocong_script.R -i <inputhost> -d <dendrogram> -c <treecmp> -o <output> -m false -r <randomtree>\n")
  q(status = 1)
}

# create output directory if it doesn't exist
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}
output_path <- gsub("/$", "", opt$output)

# check that path to TreeCmp executable exists
# Check that TreeCmp.jar exists and is executable
if (!file.exists(opt$treecmp)) {
  cat("Error: TreeCmp.jar not found in specified directory\n")
  cat("Usage: Rscript run_topocong_script.R -i <inputhost> -d <dendrogram> -c <treecmp> -o <output> -m false -r <randomtree>\n")
  q(status = 1)
}

# Check if file is executable (has read permissions)
if (!file.access(opt$treecmp, mode = 4) == 0) {
  cat("Error: TreeCmp.jar exists but is not readable/executable\n")
  cat("Usage: Rscript run_topocong_script.R -i <inputhost> -d <dendrogram> -c <treecmp> -o <output> -m false -r <randomtree>\n")
  q(status = 1)
}

##### if makerandom is true, make random trees ####
if (opt$makerandom == TRUE) {
  # print progress message
  cat("Making random trees (this runs for several minutes)...\n")
  # read in host tree
  host_tree <- read.tree(opt$inputhost)
  # get tip labels and save as a vector in comma separated list with no spaces
  tip_labels <- host_tree$tip.label
  tip_labels <- paste(tip_labels, collapse = ",")
  # write tip labels to file
  writeLines(tip_labels, file.path(output_path, "tip_labels_for_random_trees.txt"))
  
  # run the random tree script
  system(paste0(
    "bash phylosymbiosis/run_randomtrees.sh ",
    tip_labels, " ",
    file.path(output_path, "random_trees.tre")
  ))
  
  # save variable for random trees
  rand <- file.path(output_path, "random_trees.tre")
} else {
  # if makerandom is false, use the random tree provided
  rand <- opt$random
}


#### run TreeCMP on random trees versus host tree ####

cat("Running TreeCMP on random trees versus host tree...\n")

system(paste0(
  "java -jar ", opt$treecmp, " -r ", opt$inputhost, " -i ", rand, " -d ", opt$tests, " -N -o ", file.path(output_path, "host_v_random.txt")
))

#### run TreeCMP on host tree versus microbiome dendrogram ####

## fix tip labels if needed
d <- read.tree(opt$dendrogram)
h <- read.tree(opt$inputhost)
if (!all(d$tip.label %in% h$tip.label) || !all(h$tip.label %in% d$tip.label)) {
  # print message
  cat("Tip labels are not identical, fixing...\n")
  # fix tip labels
  d$tip.label <- if_else(str_detect(d$tip.label, "troglodyte"), "pan_troglodytes", d$tip.label)
  d$tip.label <- if_else(str_detect(d$tip.label, "procolobus"), "procolobus", d$tip.label)
  # check
  if (!all(d$tip.label %in% h$tip.label)) {
    cat("Error: Tip labels are not identical and cannot be fixed, please check the tip labels of the host tree and the dendrogram\n")
    q(status = 1)
  } else {
    cat("Tip labels fixed successfully\n")
    write.tree(d, file.path(output_path, "dendrogram_fixedtiplabels.newick"))
    dend <- file.path(output_path, "dendrogram_fixedtiplabels.newick")
  }
} else {
  dend <- opt$dendrogram
}

cat("Running TreeCMP on host tree versus microbiome dendrogram...\n")

system(paste0(
  "java -jar ", opt$treecmp, " -r ", opt$inputhost, " -i ", dend, " -d ", opt$tests, " -N -o ", file.path(output_path, "host_v_dendrogram.txt")
))


#### define function to calculate p values ####

myPval <- function(obs, stoch, normalized = TRUE) {
  ############ get files #################
  # read files
  obstab <- read.table(obs, sep = "\t", header = TRUE)
  sttab <- read.table(stoch, sep = "\t", header = TRUE)
  
  ########### get scores ####################
  cat("\n showing NORMALIZED scores and p values \n ")
  
  # subset to get just normalized data
  obssub <- obstab %>% dplyr::select(ends_with("_toUnifAvg"))
  stsub <- sttab %>% dplyr::select(ends_with("_toUnifAvg"))
  
  outdf <- data.frame()
  ### for each metric, get scores
  for (i in 1:length(colnames(stsub))) {
    # get observed
    obs.score <- obssub[, i]
    # get stochastic scores
    st.sc <- length(which(stsub[, i] <= obs.score))
    
    ## get p values
    p <- st.sc / nrow(stsub)
    
    # print metrics
    cat("\n -------------------------------------------")
    cat("\n total number of stochastic comparisons: ", nrow(stsub))
    cat("\n observed ", colnames(obssub)[i], ": ", obs.score)
    cat("\n p value ", colnames(obssub)[i], ": ", p)
    cat("\n -------------------------------------------")
    
    ## build results dataframe
    res <- data.frame(
      no.comps = nrow(stsub),
      metric = colnames(obssub)[i],
      obsscore = obs.score,
      pval = p,
      meanstoch = summary(stsub)[, i][4],
      minstoch = summary(stsub)[, i][1],
      maxstoch = summary(stsub)[, i][6]
    )
    
    outdf <- rbind(outdf, res)
  }
  
  
  ### return observed scores
  # remove weird text
  outdf <- outdf %>%
    mutate(across(ends_with("stoch"), ~ str_extract(., "(\\d)\\.(\\d){1,10}")))
  
  return(outdf)
}

##### run p value calculations ####
cat("Calculating p values...\n")

pval <- myPval(obs = file.path(output_path, "host_v_dendrogram.txt"), stoch = file.path(output_path, "host_v_random.txt"), normalized = TRUE)

## write output p values to file
write.table(pval, file.path(output_path, "p_values.txt"), sep = "\t", row.names = FALSE)

# print finished message with the time it took to run and system date
cat("\nFinished at: ", format(Sys.time(), "%a %b %d %X %Y"), "\n")
