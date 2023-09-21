# Originally Written by Cameron MacPherson
# Modified by Preston Leung (Take in multiple entrez ID rather than just a single one)

suppressWarnings(
  suppressPackageStartupMessages(
    {
      library(optparse)
      library(dplyr)
      library(purrr)
      library(magrittr)
    }
  )
)

#- Load functions
source("/dockerWorkDir/codes/GeneAnalysis_SKATO_Helper.R")


# set work
# setwd("L:/LovbeskyttetMapper/persimune/Preston/Code_from_computerome/script/")

# set default output directory to save results to
default_dir = getwd()
default_fileName = "GeneAnalysis_SKATO_Results.tsv"
default_filePath = file.path(default_dir,default_fileName)




# - Define the help message so that if we don't know what
#   to put in, we can have this to help out
help_msg = paste(
  "GeneAnalysis_SKATO.R is a tool to help calculate associations of SNP",
  "and groups of SNPs using SKAT-O. Users can edit the param file to",
  "configure the tool."
)


#- Defining options when running tool outside StudioR
option_list = list(
  make_option(c("-p", "--param_file"), action="store", default=NULL, type='character',
              help="parameter file to configure SKAT-O run. Must be in comma separated format."),
  make_option(c("-o", "--out"), action="store", default= default_filePath,
              help="Output tab-separated result. Default name will be current directory with file name: %default.")
)

opt = parse_args(
  OptionParser(
    option_list = option_list,
    description = help_msg
  )
)

#- Define map function | returns a scalar
update.cohort = function(cohort) last.cohort <<- cohort
reset.cohort = function() last.cohort <<- "NONE"
do = function(entrez, window, pidfile, confounders, cohort, outcome, affy_file, pheno_file, geno_file, ignore = 0) {
  #- Inform user
  cat(
    sprintf(
      "Cohort=%s\nOutcome=%s\nEntrez=%s\nWindow=%s\nConfounders=%s\n",
      cohort, 
      outcome, 
      entrez, 
      window, 
      paste0("[",paste0(confounders, collapse="+"),"]")
    )
  )
  message("######################")
  message("# Step 1 - Load data #")
  message("######################")
  #- Load data
  if (cohort != last.cohort) {
    cat(sprintf("loading %s...", cohort))
    cat("\n")    
    loadData(
      affy_file = affy_file,
      pheno_file = pheno_file,
      geno_file = geno_file
    )
    cat("LOADED!\n")
  }
  update.cohort(cohort)
  setDefault.TestMutator(log) # log is the log function
  
  #- Exit early if outcome not in cohort (the column title is specific to each cohort data)
  # e.g. START cohort outcome of rnabase line is named RNA_00. 
  # while in FIRST it is named rnabl
  if (!(outcome %in% colnames(DATA__$PHENOTYPES))) {
    message("Outcome not found in cohort, returning 2.\n")
    return (2)
  }
  
  #- Specify query. pids are patient ids from the cohorts 
  pids = read.table(pidfile)[,1]
  
  message("#######################")
  message("# Step 2 - Setup test #")
  message("#######################")
  
  #- Setup test - here might be able to modify the entrez param for 1+ genes
  setupTest.withEntrez(
    # entrez = entrez, 
    entrez = unlist(stringr::str_split(entrez,"_")),
    controls = confounders, 
    outcome = outcome, 
    pad = window, 
    pids = pids
  )
  
  TEST__$y[is.infinite(TEST__$y)] <<- NA  
  
  message("###########################")
  message("# Step 3 - Run SKATO test #")
  message("###########################")
  #- Run test
  runSKAT.SKATO()
  #- Return p.value
  result = getResultsSKAT.SKATO()$p.value
  #- Inform user
  # cat(sprintf("p-value = %s", result))
  message("# SKATO test Completed!!! #")  
  return(result)
}


helpMsg = function(){
  print_help(OptionParser(
    option_list=option_list,
    description = help_msg
  ))
  message("-------------------------------------------")
  message("Missing needed files. Exit.")
}


####################################################
#- Read in parameters (Change parameter file here) #
####################################################

# & !is.null(opt$data_path_file)
if(!is.null(opt$param_file)){
  if(opt$out != default_filePath){
    default_filePath = opt$out
    message("# Set default file path: ", default_filePath," #")
  }
  message("# Reading SKAT-O parameter file #")
  parameters = read.delim(opt$param_file,sep=",",as.is=TRUE)  
  
  #- Tidy magic to collapse confounders to a vector
  parameters$confounders = parameters %>%
    dplyr::select(dplyr::starts_with("confounder.")) %>%
    purrr::pmap(function (...) {list(...)}) %>%
    purrr::map(function (listOfConfounders) as.character(listOfConfounders))
  
  #- Remove all confounder.# columns
  reset.cohort()
  parameters %<>% dplyr::select(-dplyr::starts_with("confounder."))

  #- Now we map the parameters to run SKAT-O
  results = purrr::pmap_dbl(parameters, do)
  result_formatted = tibble(
    entrez = parameters$entrez,
    cohort = parameters$cohort,
    outcomes = parameters$outcome,
    pvalue = results
  )
  data.table::fwrite(
    result_formatted,
    file = opt$out,
    sep = "\t"
  )
}else{
  helpMsg()
}
