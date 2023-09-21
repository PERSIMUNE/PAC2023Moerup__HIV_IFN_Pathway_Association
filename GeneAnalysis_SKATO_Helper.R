# source("/home/projects/cu_10103/people/cammac/TET2CGA/script/checkDependencies.R", echo = FALSE, verbose = FALSE)
source("/dockerWorkDir/codes/checkDependencies.R", echo = FALSE, verbose = FALSE)
library(magrittr)

#### GLOBALS ----

DATA__ = list(AFFYPROBE = NA, PHENOTYPES=NA, GENOTYPES=NA)
TEST__ = list(type = "TEST__", X=NA, y=NA, Z=NA, results = list())
passthrough = function(a) a
default.mutator = passthrough
`%/%` = function(a,b) paste(a,b,sep="/")
`%.%` = function(a,b) paste(a,b,sep=".")

#### EXTERNALS

# Original Code by Cameron
# INPUT: Single entrez ID
# OUTPUT: a named list with names: chrom, start, end
getGenomicCoordinate.byEntrez = function(entrez) {
  gene = mygene::getGene(entrez, fields = c("genomic_pos_hg19"))[[1]]$genomic_pos_hg19
  
  # if gene co-ordinates not found, put default location for ease of checking
  if (is.null(gene)) return(list(chrom = 1, start = -999, end = -998)) 
  if (is.null(names(gene))) gene = gene[[1]]
  list(chrom = gene$chr, start = min(gene$start, gene$end), end = max(gene$start, gene$end))
}

# New code modded by pres 
# this function is to replace getGenomicCoordinate.byEntrez to cater for 
# 1+ entrez ids given. Now returns a tibble rather than a list with named variables
getGenomicCoordinate.byListOfEntrez = function(entrez) {
  gene = NULL
  for (entrez_id in entrez){
    gene_temp = mygene::getGene(entrez_id, fields = c("genomic_pos_hg19"))[[1]]$genomic_pos_hg19
    
    # if gene co-ordinates not found, put default location for ease of checking
    if (is.null(gene_temp)){
      # gene_temp = tibble(chrom = 1, start = -999, end = -998, strand = 0)
      gene_temp = tibble(chr = 1, start = -999, end = -998)
    }else{
      gene_temp = tibble(
        chr = gene_temp$chr,
        start = min(gene_temp$start, gene_temp$end),
        end = max(gene_temp$start, gene_temp$end)
      )
    }
    # not exactly sure where does this case come into effect, but was in getGenomicCoordinate.byEntrez
    # Here the condition is tested, if it goes through, it'll print a message notifying the case
    # and then does nothing.
    if (is.null(names(gene_temp)) & !is.null(gene_temp)){
      message("Weird case where gene_temp is not null but names(gene_temp) is null.")
      # gene_temp = gene_temp[[1]] <== reactivate once we know what this is doing
    } 
    gene = rbind(gene,gene_temp)
  }
  colnames(gene)[1] = c("chrom")
  return(gene)
}


#### DATA__ functions ----
resetData = function() DATA__ <<- list(AFFYPROBE = NA, PHENOTYPES=NA, GENOTYPES=NA)

loadData = function(affy_file, pheno_file, geno_file){
  cat(paste("Now loading:", pheno_file))
  cat("\n")
  DATA__$PHENOTYPES <<- readRDS(pheno_file)
  cat(paste("Now loading:", affy_file))
  cat("\n")
  DATA__$AFFYPROBE  <<- readRDS(affy_file)
  cat(paste("Now loading:", geno_file))
  cat("\n")
  DATA__$GENOTYPES <<- as.matrix(readRDS(geno_file))
}

overridePHENO = function(path) {
  DATA__$PHENOTYPES <<- readRDS(path)
}

getRS.byGenomicCoordinate = function(query = list(chrom = 1, start = 1, end = 100)) {
  result = DATA__$AFFYPROBE %>% 
    dplyr::filter(
      Chrom == as.numeric(query$chrom) & 
      dplyr::between(Position.HG19, query$start, query$end)
    ) %>% 
    dplyr::select(rsID)
  result[[1]]
}

# Original code by Cameron
# INPUT: entrez id of ONE gene, padding to add +/- regions
# OUTPUT: corresponding rsIDs within the region of the gene given
# by entrez id (in a list form)
getRS.byEntrez = function(entrez, pad = 0) {
  genomicCoordinate = getGenomicCoordinate.byEntrez(entrez) # genomicCoordinate here is a named list
  genomicCoordinate$start = genomicCoordinate$start - pad
  genomicCoordinate$end   = genomicCoordinate$end   + pad
  getRS.byGenomicCoordinate(genomicCoordinate) # a list of rsids
}


# New code modded by pres
# INPUT: a list of entrez ids (multiple genes), padding to add +/- regions
# OUTPUT: a tibble with 2 columns: Entrez id and list of rsids
getRS.byListOfEntrez = function(entrez, pad = 0){
  # genomicCordinate is a tibble with N rows for multiple entrez id
  genomicCoordinate = getGenomicCoordinate.byListOfEntrez(entrez)
  genomicCoordinate$start = genomicCoordinate$start - pad
  genomicCoordinate$end   = genomicCoordinate$end   + pad
  
  x = NULL
  entrez_rs = bind_rows(
    lapply(
      seq(1,nrow(genomicCoordinate)),
      function(i){
        x = tibble(entrez = entrez[i])
        x$rsID = list(getRS.byGenomicCoordinate(genomicCoordinate[i,]))
        return(x)
      }
    )
  )
  return(entrez_rs)
}


# Original code by Cameron
getAffyID.fromRS = function(query) {
  query %<>% rsFormatFilter()
  result = DATA__$AFFYPROBE %>% 
    dplyr::filter(rsID %in% query) %>% 
    dplyr::select(AffyID)
  result[[1]] %>% 
    DATA__probeParityFilter()
}

# new code modded by Pres
getAffyID.fromTibbleRS = function(query) {
  query_rsid = unlist(query$rsID)
  query_rsid%<>% rsFormatFilter()
  result = DATA__$AFFYPROBE %>% 
    dplyr::filter(rsID %in% query_rsid) %>% 
    dplyr::select(AffyID)
  result[[1]] %>%
    DATA__probeParityFilter()
}

# Original code by Cameron
getAffyID.fromEntrez = function(entrez, pad = 0) {
  rsIDs = getRS.byEntrez(entrez, pad)
  getAffyID.fromRS(rsIDs)
}

# New code modded by Pres
getAffyID.fromListOfEntrez = function(entrez, pad = 0) {
  entrez_rsIDs = getRS.byListOfEntrez(entrez, pad)
  getAffyID.fromTibbleRS(entrez_rsIDs)
}


filterGENOTYPES.byENTREZ = function(entrez, pad = 0) {
  # affyIDs = getAffyID.fromEntrez(entrez, pad) # original code by Cameron  
  affyIDs = getAffyID.fromListOfEntrez(entrez, pad) # new code modded by Pres  
  DATA__$GENOTYPES[,affyIDs,drop=FALSE]
}

filterGENOTYPES.byRS = function(rsIDs) {
  affyIDs = getAffyID.fromRS(rsIDs)
  DATA__$GENOTYPES[,affyIDs]
}

filterGENOTYPES.byAFFY = function(affyIDs) {
  affyIDs %<>% DATA__probeParityFilter()
  DATA__$GENOTYPES[,affyIDs]
}

DATA__probeParityFilter = function(affyIDs) {
  affyIDs[affyIDs %in% colnames(DATA__$GENOTYPES)]
}

rsFormatFilter = function(rsIDs) {
  rsIDs[stringr::str_detect(rsIDs, "^rs")]
}

listPHENOTYPES = function() {
  cat(paste0(paste(colnames(DATA__$PHENOTYPES), collapse="\n"),"\n"))
}

getPHENOTYPES = function() colnames(DATA__$PHENOTYPES)

#### TEST__ functions ----
#- These functions handle selection and mutation of genotype/phenotype from DATA__ to TEST__
#- These functions manipulate X, y, Z, W, type, and results within TEST__

resetTest   = function(   ) TEST__ <<- list(type = "TEST__", X=NA, y=NA, Z=NA, W=NA, 
                                            results = list(SKAT = list(),
                                                           PDA  = list(),
                                                           PLA  = list()
                                            ))

getTest    = function(   ) TEST__

setTest = function(obj) {
  if (isTest(obj)) TEST__ <<- obj
  else cat("Restore failed, not a TEST__ object.\n")
}

isTest  = function(obj) {
  if (!is.list(obj)) return (FALSE)
  if (!any("type" %in% names(obj))) return (FALSE)
  if (obj$type == "TEST__") return (TRUE)
  FALSE
}

syncTestPIDs = function(pids = NULL) {
  if (is.null(pids)) pids = rownames(TEST__$Z)
  if (!is.character(pids)) pids %<>% as.character()
  pids %<>% filterPIDs.byTestParity()
  TEST__$Z <<- TEST__$Z[pids,,drop=FALSE]
  TEST__$PHENOTYPES <<- (DATA__$PHENOTYPES %>% tibble::column_to_rownames("PID"))[pids,] %>% tibble::rownames_to_column("PID")
  cat(paste0(length(pids), " patient IDs included.\n"))
}

removeAllPhenotypeNAs = function() {
  remove.these = c(which(rowSums(is.na(TEST__$X)) > 0), which(is.na(TEST__$y)), which(is.infinite(TEST__$y)))
  TEST__$X <<- TEST__$X[-remove.these,]
  TEST__$y <<- TEST__$y[-remove.these]
  TEST__$Z <<- TEST__$Z[-remove.these,]
}

setupTest.withEntrez = function(entrez, controls, outcome, pad = 0, pids = NULL, mutate.outcome = NA, ...) {
  resetTest()
  TEST__$Z <<- filterGENOTYPES.byENTREZ(entrez, pad)
  setup__(controls, outcome, pids, mutate.outcome, ...)
}

setupTest.withRS = function(rsIDs, controls, outcome, pids = NULL, mutate.outcome = NA, ...) {
  resetTest()
  TEST__$Z <<- filterGENOTYPES.byRS(rsIDs)
  setup__(controls, outcome, pids, mutate.outcome, ...)
}

setupTest.withAffy = function(affyIDs, controls, outcome, pids = NULL, mutate.outcome = NA, ...) {
  resetTest()
  TEST__$Z <<- filterGENOTYPES.byAFFY(affyIDs)
  setup__(controls, outcome, pids, mutate.outcome, ...)
}

setup__ = function(controls, outcome, pids, mutate.outcome, na.remove = FALSE, silent = FALSE, ...) {
  syncTestPIDs(pids) #must be run here; needs error handling
  setTestControls(controls)
  setTestOutcome(outcome, mutate.outcome)
  if (na.remove) removeAllPhenotypeNAs()
  if (!silent) cat(paste0(ncol(TEST__$Z), " probes pooled.\n"))
}

setTestControls = function(controls) {  
  TEST__$X <<- as.matrix(TEST__$PHENOTYPES %>% dplyr::select(all_of(dplyr::one_of(controls))))
  cat(paste0(ncol(TEST__$X), " of ", length(controls), " control variables set.\n"))
}

setTestOutcome = function(outcome, override = NA) {
  if (is.na(override)) {
    message("SetTestOutcome -> Path 1.")
    TEST__$y <<- (TEST__$PHENOTYPES %>% dplyr::select(all_of(outcome)))[[1]] %>% default.mutator()    
  } else {
    message("SetTestOutcome -> Path 2.")
    TEST__$y <<- (TEST__$PHENOTYPES %>% dplyr::select(all_of(outcome)))[[1]] %>% override()
  }
  cat(paste0("Outcome is '", outcome,"'.\n"))
}

mutateTestOutcome = function(f) {
  TEST__$y <<- f(TEST__$y)
}

logTestOutcome = function(base=exp(1)) {
  TEST__$y <<- log(TEST__$y, base)
}

ordinateTestOutcome = function(classes = 2, method = "content") {
  if (!is.numeric(TEST__$y)) {
    cat("Outcome not numeric, can not convert to a nominal variable.")
  }
  TEST__$y.ord <<- OneR::bin(TEST__$y, nbin = classes, labels = paste0("L",1:classes), method = method) %>% as.character()
}

#### SKAT functions ----
#- These functions make use of X, y, and Z from TEST__

runSKAT.CommonRare = function() {
  obj = SKAT::SKAT_Null_Model(TEST__$y ~ TEST__$X, out_type="C")
  TEST__$results$SKAT$CommonRare <<- SKAT::SKAT_CommonRare(TEST__$Z, obj, r.corr.common = 0, r.corr.rare = 1, method="C")
}

runSKAT.SKATO = function(kernel = "linear.weighted", beta.shape = c(1,25), impute.geno.method = "fixed", missingness.allowed = 0.15, commonness.allowed = 1, maf.estimate.method = "all", resample.n = 0) {
  obj = SKAT::SKAT_Null_Model(TEST__$y ~ TEST__$X, out_type="C", n.Resampling = resample.n)
  TEST__$results$SKAT$SKATO <<- SKAT::SKAT(TEST__$Z, obj, method = "SKATO", kernel = kernel, 
                                           weights.beta = beta.shape, impute.method = impute.geno.method, 
                                           missing_cutoff = missingness.allowed, max_maf = commonness.allowed, 
                                           estimate_MAF = ifelse(maf.estimate.method == "all",1,2))
}

runSKAT.SKATOBinary = function(kernel = "linear.weighted", method.bin = "Hybrid", resampling = 2*2^-6, seednum = NULL, beta.shape = c(1,25), impute.geno.method = "bestguess", missingness.allowed = 0.15, commonness.allowed = 1, maf.estimate.method = "all", resample.n = 0) {
  obj = SKAT::SKAT_Null_Model(TEST__$y ~ TEST__$X, out_type="D", n.Resampling = resample.n)
  TEST__$results$SKAT$SKATOBinary <<- SKAT::SKATBinary(TEST__$Z, obj, method = "SKATO", kernel = kernel, 
                                           method.bin = method.bin, N.Resampling = resampling, seednum = seednum,
                                           weights.beta = beta.shape, impute.method = impute.geno.method, 
                                           missing_cutoff = missingness.allowed, max_maf = commonness.allowed, 
                                           estimate_MAF = ifelse(maf.estimate.method == "all",1,2))
}

getResultsSKAT.CommonRare  = function() TEST__$results$SKAT$CommonRare
getResultsSKAT.SKATO       = function() TEST__$results$SKAT$SKATO
getResultsSKAT.SKATOBinary = function() TEST__$results$SKAT$SKATOBinary

setDefault.TestMutator = function(f) {
  default.mutator <<- f
}

filterPIDs.byTestParity = function(pids) {
  pids[pids %in% rownames(TEST__$Z)]
}

filterAffy.byTestParity = function(affyIDs) {
  affyIDs[affyIDs %in% colnames(TEST__$Z)]
}

#### PLA functions ----
#- These functions make use of X, y, and Z from TEST__

findControlResiduals.byLinearRegression = function() {
  df = as.data.frame(cbind(TEST__$X, TEST__$y))
  colnames(df)[ncol(df)] = "O"
  TEST__$y.res <<- lm(O ~ ., data = df)$residuals
}

findOrthogonalGenotypes.byPCA = function(discard=c()) {
  if (any(is.na(TEST__$Z))) {
    geno = TEST__$Z
    geno[is.na(geno)] = 0
  } else geno = TEST__$Z
  pca = prcomp(geno)
  if (length(discard) > 0) pca$x = pca$x[,-discard]
  TEST__$Z.pca <<- pca
  TEST__$W <<- pca$x
}

ordinateTestResiduals = function(classes = 2, method = "content") {
  TEST__$y.res.ord <<- OneR::bin(TEST__$y.res, nbin = classes, labels = paste0("L",1:classes), method = method, na.omit = FALSE)
}

runPDA = function(removePCs=c()) {
  findControlResiduals.byLinearRegression()
  findOrthogonalGenotypes.byPCA(removePCs)
  ordinateTestResiduals(classes = 3, method = "length")
  if (nrow(TEST__$W) != length(TEST__$y.res.ord)) {
    cat("Error: Dimensions of ordinals and principal genotypes do not agree. Please try run `removeAllPhenotypeNAs()` before runPDA().\n")
    return ()
  }
  df = cbind(TEST__$W, TEST__$y.res.ord) %>% as.data.frame()
  colnames(df)[ncol(df)] = "O"
  TEST__$results$PDA$lda    <<- MASS::lda(O ~ ., data = df)
  TEST__$results$PDA$lda.cv <<- MASS::lda(O ~ ., data = df, CV = TRUE)
  TEST__$results$PDA$manova <<- manova(as.matrix(df %>% dplyr::select(-O)) ~ O, data = df)
  TEST__$results$PDA$pca    <<- TEST__$Z.pca
}

runPLA = function(removePCs=c()) {
  findControlResiduals.byLinearRegression()
  findOrthogonalGenotypes.byPCA(removePCs)
  if (nrow(TEST__$W) != length(TEST__$y.res)) {
    cat("Error: Dimensions of residuals and principal genotypes do not agree. Please try run `removeAllPhenotypeNAs()` before runPLA().\n")
    return ()
  }
  df = cbind(TEST__$W, TEST__$y.res) %>% as.data.frame()
  colnames(df)[ncol(df)] = "O"
  TEST__$results$PLA$lm  <<- df %>% 
                             dplyr::select(-O) %>%
                             purrr::map( ~lm(df$O ~ .x, data = df)) %>% 
                             purrr::map(summary) %>% 
                             purrr::map_df(broom::glance) %>%
                             tibble::rowid_to_column("PC") %>% 
                             dplyr::arrange(p.value) %>% 
                             round(3)
  TEST__$results$PLA$pca <<- TEST__$Z.pca
}

getResultsPDA = function() TEST__$results$PDA
getResultsPLA = function() TEST__$results$PLA
