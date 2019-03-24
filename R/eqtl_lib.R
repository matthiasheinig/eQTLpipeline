#' A package to facilitate the full pipeline of cis eQTL analysis including
#' permutations on a HPC cluster supported by BatchJobs
#'
#' \tabular{ll}{
#' Package: \tab eQTLpipeline\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1\cr
#' Date: \tab 2015-11-03\cr
#' License: \tab LGPL \cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @name eQTLpipeline
#' @aliases eQTLpipeline
#' @docType package
#' @title full pipeline of eQTL analysis including permutations on HPC systems
#' @docType package
#' @references 
#' @keywords package
#' @import MatrixEQTL
#' @import data.table
#' @import BatchJobs
#' @examples
#' #insertexample()
NULL



#' Basic interface to MatrixEQTL passing all data in files
#'
#' @family eqtl functions
#' @title basic eQTL interface
#' @param expression_file_name filename of a tab separated file with expression
#'        values in a matrix ngene x nsample
#' @param genotype_file_name filename of a tab separated file with genotype
#'        dosage values (usually between 0 and 2) in a matrix nsnp x nsample
#' @param covariates_file_name filename of a tab separated file with covariate
#'        values in a matrix ncovar x nsample
#' @param gene.position data.frame with gene positions in the first 4 columns
#'        named "chrom", "start", "end", "gene_id"
#' @param snp.pos data.frame with SNP positions in the first three columns
#'        named "snp_id", "chrom", "snp_pos"
#' @param prefix name prefix for the outputfilenames
#' @param redo logical, if TRUE the results will always be recomputed, if FALSE
#'        the results will be loaded if they exist (default FALSE)
#' @param threshold significance threshold for eQTLs (default 1e-5)
#' @param compute.all MatrixEQTL reports only significant results, if results
#'        for all tests are required set this to TRUE (default FALSE)
#' @param what flag to determine what should be returned (default "table"). Can
#'        be set to "object" to get the MatrixEQTL result object.
#' @return if what is "table" the function returns the eQTL results as a table
#'         if what is "object" the function returns the MatrixEQTL object.
#' @references
#' @export
eqtl <- function(expression_file_name, genotype_file_name, covariates_file_name, gene.position, snp.pos, prefix, redo=FALSE, threshold=1e-5, compute.all=FALSE, what="table") {
  ## require(MatrixEQTL)
  ## require(data.table)
  
  # we use this check sum to store intermediate results
  if (compute.all) {
    eqtl.file = paste(prefix, "_eQTL_results_R.txt", sep="")
  } else {
    eqtl.file = paste(prefix, "_eQTL_results_R_significant.txt", sep="")
  }
  if (file.exists(eqtl.file) && !redo) {
    eqtl = fread(eqtl.file, sep="\t", stringsAsFactors=F)
    return(eqtl)
  }
  
  dir.create(dirname(prefix), recursive=T)
 
    
  # compute eqtl
  useModel = modelLINEAR ; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

  # Only associations significant at this level will be saved
  if (compute.all) {
    pvOutputThreshold_cis = 1
  } else {
    pvOutputThreshold_cis = threshold
  }
  
  pvOutputThreshold_tra = 0 # only look at cis eqtl
  
  errorCovariance = numeric()
  # errorCovariance = read.table("Sample_Data/errorCovariance.txt")
  cisDist = 1e6

  # load the expression data for matrix eqtl
  gene = SlicedData$new();
  gene$fileDelimiter = "\t"; # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1; # one row of column labels
  gene$fileSkipColumns = 1; # one column of row labels
  gene$fileSliceSize = 2000; # read file in pieces of 2,000 rows
  gene$LoadFile(expression_file_name );

 

  # load the covariate data for matrix eqtl
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t"; # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1; # one row of column labels
  cvrt$fileSkipColumns = 1; # one column of row labels
  cvrt$fileSliceSize = 2000; # read file in one piece
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }  
  

  ## Load genotype data
      
  snps = SlicedData$new();
  snps$fileDelimiter = "\t"; # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1   ; # one row of column labels
  snps$fileSkipColumns = 1; # one col of snp ids
  snps$fileSliceSize = 2000; # read file in pieces of 2,000 rows
  snps$LoadFile(genotype_file_name);
      
  ## Run the analysis
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = "",
    pvOutputThreshold = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = eqtl.file,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snp.pos[,1:3],
    genepos = gene.position[,1:4],
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp=(what != "table"));

  pdf(file=paste(prefix, "_qqplot.pdf", sep=""))
  plot(me)
  dev.off()

  if (what == "table") {
    eqtl = fread(eqtl.file, sep="\t", stringsAsFactors=F)

    ## add information about the snp position to the eqtl data
    eqtl = cbind(eqtl, snp.pos[match(eqtl$SNP, snp.pos[,"snp_id"]),-1], gene.position[match(eqtl$gene, gene.position[,1]),-1])
    eqtl = eqtl[order(eqtl$chrom, eqtl$snp_pos, eqtl$gene),]
    write.table(eqtl, eqtl.file, sep="\t", row.names=F, quote=F)
    
    ## also create a filtered version
    if (compute.all) {
      significant = eqtl[eqtl[["p-value"]] < threshold,]
      sig.file = gsub("_eQTL_results_R.txt$", "_eQTL_results_R_significant.txt", eqtl.file)
      write.table(significant, sig.file, sep="\t", row.names=F, quote=F)
    }
    return(invisible(eqtl))
  } else {
    return(me)
  }
}

## they use min(p) per gene as the test statistic
## then they permute phenotypes to derive the null distribution
## they run 1000 permutations for all genes
## the exit criterion is > 15 permuted min(p) > actual min(p)


#' Basic interface to MatrixEQTL passing expression and covariate data in
#' matrices, and genotypes as a file uses the \code{\link{eqtl}} function
#'
#' @family eqtl functions
#' @title basic eQTL interface
#' @param expr expression values in a matrix ngene x nsample
#' @param covar covariate values in a matrix ncovar x nsample
#' @param genotype_file_name filename of a tab separated file with genotype
#'        dosage values (usually between 0 and 2) in a matrix nsnp x nsample
#' @param gene.position data.frame with gene positions in the first 4 columns
#'        named "chrom", "start", "end", "gene_id"
#' @param snp.pos data.frame with SNP positions in the first three columns
#'        named "snp_id", "chrom", "snp_pos"
#' @param id character indicating an id used for tempfile creation
#' @return the same output as the \code{\link{eqtl}} function
#' @references
#' @export
eqtl.run <- function(expr, covar, genotype_file_name, gene.position, snp.pos, id, ...) {

  ## prepare the input files
  outdir = tempfile(pattern=id)
  cat("creating temp dir for eqtl run", outdir, "\n")
  dir.create(outdir)

  prefix = file.path(outdir, id)
  expr.file = paste(prefix, "_expression.txt", sep="")
  write.table(expr, file=expr.file, sep="\t", quote=F)
  
  if (!is.null(covar)) {
    covar.file = paste(prefix, "_covariates.txt", sep="")
    write.table(covar, file=covar.file, sep="\t", quote=F)
  } else {
    covar.file = NULL
  }

  ## run the eqtl analysis
  eqtls = eqtl(expr.file, genotype_file_name, covar.file, gene.position, snp.pos, prefix, ...)

  return(eqtls)
}

#' Interface to the MatrixEQTL Matrix_eQTL_engine()
#' This is a simplified version for more general (trans)QTL analyses,
#' e.g. GWAS or metaboliteQTL and does not need gene annotations
#'
#' @family eqtl functions
#' @title transQTL interface
#' @param prefix path and prefix for the outputfiles (default NULL)
#' @param genotype_file_name filename of a tab separated file with genotype
#'        dosage values (usually between 0 and 2) in a matrix nsnp x nsample
#' @param vcf is the genotype file a vcf file that needs to be converted?
#'        (default FALSE)
#' @param snp.pos data.frame with SNP positions to be added to matrixEQTL object
#'        in the end. Columns named "snps", "chr", "pos", containing
#'        information in that order (default NULL)
#'        if vcf=TRUE and snp.pos=TRUE | snp.pos=NULL, information is extracted
#'        from vcf file unless save.memory=TRUE.
#'        if snp.pos=FALSE, no annotation is added
#' @param expression_file_name filename of a tab separated file with expression
#'        values in a matrix ngene x nsample.
#'        if expression_file_name is a numeric vector (lenght nsample),
#'        matrix (ngene x nsample) or data.frame (ngene x nsample), the
#'        expression file is created with the name prefix_gene.txt.
#' @param covariates_file_name filename of a tab separated file with covariate
#'        values in a matrix ncovar x nsample
#'        to not use any covariates, set covariates_file_name=character()
#'        if covariates_file_name is a numeric vector (lenght nsample),
#'        matrix (ngene x nsample) or data.frame (ngene x nsample), the
#'        covariate file is created with the name prefix_covariates.txt.
#' @param threshold significance threshold for eQTLs (default 1e-5), set to 1 to
#'        get all QTLs. threshold set to 1 for compute.all=TRUE
#' @param compute.all MatrixEQTL reports only significant results, if results
#'        for all tests are required set this to TRUE (default FALSE)
#' @param save.memory prevent FDR correction (default FALSE). only functional in
#'        the github-version of MatrixEQTL, see also
#'        https://github.com/andreyshabalin/MatrixEQTL/issues/8
#' @param what flag to determine if a table or the MatrixEQTL object should be
#'        returned (default NULL). If unspecified, a table is returned except
#'        either save.memory=TRUE or min.pv.by.genesnp=TRUE, then the MatrixEQTL
#'        result object is returned.
#' @param min.pv.by.genesnp option to report the minimal pvalue per SNP and per
#'        gene in the MatrixEQTL object (default FALSE). If you want to use this
#'        option, make sure to set what="object" to get the results.
#' @param verbose print comments of progress (default TRUE)
#' @param redo logical, if TRUE the results will always be recomputed, if FALSE
#'        the results will be loaded if they exist (default TRUE)
#' @param load.qtls (default FALSE) only relevant for save.memory=TRUE which is
#'        mostly used for big datasets. In that case, QTL results are not
#'        reported in the MatrixEQTL object but can be read in with
#'        load.qtls=TRUE.
#' @param keep.outfile unless compute.all=T or save.memory=T, results are written
#'        to a temp file that is deleted. (default FALSE)
#' 
#' Parameters /options I would like to add:
#' @param id.check check if ids and order of ids match in genotype, expression
#'        and covariate data (default FALSE)
#' @param threshold.cis
#' @param threshold.trans
#'
#' Parameters/options in the old wrapper that do not work currently:
#' @param gene.position data.frame with gene positions in the first 4 columns
#'        named "gene_id", "chrom", "start", "end"
#'        
#'
#' @return if what="table", the data.frame contained in the MatrixEQTL object
#'         me$all$eqtls is returned. If what="object", the whole MatrixEQTL
#'         object is returned.
#' @author Ines Assum (2019-03-21)
#' @references
#' @export
trans.qtl <- function(prefix,
                      genotype_file_name,
                      vcf=FALSE,
                      snp.pos=NULL,
                      expression_file_name,
                      covariates_file_name=character(),
                      threshold=1e-5,
                      compute.all=FALSE,
                      what=NULL,
                      min.pv.by.genesnp=FALSE,
                      save.memory=FALSE,
                      verbose=TRUE,
                      redo=TRUE,
                      load.qtls=FALSE,
                      keep.outfile=FALSE,
                      id.check=FALSE,
                      gene.position=NULL) {

  # TODO:
  # implement id.check:
  # check if ids and order of ids match in genotype, expression and covariate data

  # libraries:
  require(MatrixEQTL)
  require(data.table)
  if(vcf){
    require(dplyr)
    require(VariantAnnotation)
    require(snpStats)
  }
  
  if(is.null(prefix)){
    if(!dir.exists(paste0(getwd(), "/matrixEQTL"))){
      dir.create(paste0(getwd(), "/matrixEQTL"), showWarnings = F, recursive = T)
    }
    prefix <- paste0(getwd(), "/matrixEQTL/test")
  }
  
  if (is.null(what)){
    what <- "table"
    if (save.memory | min.pv.by.genesnp){
      what <- "object"
    }
  }
  
  # Output file name
  if (save.memory | compute.all | keep.outfile){
    output_file_name = paste0(prefix, "_QTL_results.txt");
  } else {
    output_file_name = tempfile();
  }
  
  if (file.exists(output_file_name) && !redo) {
    cat("QTL_results already exist in file\n", output_file_name, "\n",
        "and will not be recomputed because redo=F\n")
    eqtl = fread(output_file_name, sep="\t", stringsAsFactors=F)
    colnames(me$all$eqtls) <- c("snps", "gene", "beta", "statistic", "pvalue")
    return(eqtl)
  }
  
  
  ## Location of the package with the data files.
  base.dir = prefix;
  
  # data prep

  if(vcf){
    vcf.file <- readVcf(genotype_file_name)
    geno.file <- paste0(prefix, "_SNP.txt")
    gt.mat <- t(as(genotypeToSnpMatrix(vcf.file)$genotype, "numeric")) %>% as.data.frame
    write.table(cbind(snpid=rownames(gt.mat), gt.mat),
                row.names = F, col.names = T, quote = F, sep = "\t",
                file = geno.file)
    if(verbose){
      cat("Created genotype file from .vcf file at:\n", geno.file, "\n")
    }
    genotype_file_name <- geno.file
    if(is.null(snp.pos) & vcf==T & isFALSE(save.memory)){
      snp.pos <- T
    }
    if(isTRUE(snp.pos)){
      snp.pos <- rowRanges(vcf.file)
      snp.pos <- data.frame(snps=names(snp.pos),
                            chr=as.numeric(as.character(seqnames(snp.pos))),
                            pos=start(snp.pos),
                            stringsAsFactors = F)
      if(verbose){
        cat("SNP annotations extracted from the .vcf file.\n")
      }
    }
  }
  if(is.data.frame(snp.pos)){
    if(verbose){
      cat("SNP position annotation will be added to the results in the end.\n")
    }
  }
  
  if (is.character(expression_file_name) && length(expression_file_name)>0) {
    if(verbose){
      cat("Expression file will be read from:\n", expression_file_name, "\n")
    }
  } else if (is.data.frame(expression_file_name)){
    expr.file <- paste0(prefix, "_gene.txt")
    write.table(cbind(geneid=rownames(expression_file_name), expression_file_name),
                row.names = F, col.names = T, quote = F, sep = "\t",
                file = expr.file)
    if(verbose){
      dim <- dim(expression_file_name)
      cat("Created expression file from data.frame for ",
          dim[1], " genes by ",
          dim[2], " individuals at:\n", expr.file, "\n")
    }
    expression_file_name <- expr.file
  } else if (is.matrix(expression_file_name)){
    expr.file <- paste0(prefix, "_gene.txt")
    write.table(cbind(geneid=rownames(expression_file_name), expression_file_name),
                row.names = F, col.names = T, quote = F, sep = "\t",
                file = expr.file)
    if(verbose){
      dim <- dim(expression_file_name)
      cat("Created expression file from matrix for ",
          dim[1], " genes by ",
          dim[2], " individuals at:\n", expr.file, "\n")
    }
    expression_file_name <- expr.file
  } else if (is.vector(expression_file_name) && is.numeric(expression_file_name)) {
    expr.file <- paste0(prefix, "_gene.txt")
    write.table(rbind(c("geneid", paste0("sample", 1:length(expression_file_name))),
                      c("gene1", expression_file_name)),
                row.names = F, col.names = F, quote = F, sep = "\t",
                file = expr.file)
    if(verbose){
      dim <- length(expression_file_name)
      cat("Created expression file for a single gene for", dim,
          "individuals at:\n", expr.file, "\n")
    }
    expression_file_name <- expr.file
  } else {
    error("Something went wrong with importing your expression data.")
  }
  
  if (is.character(covariates_file_name) && length(covariates_file_name)>0) {
    if(verbose){
      cat("Covariate file will be read from:\n", covariates_file_name, "\n")
    }
  } else if (is.data.frame(covariates_file_name)){
    cov.file <- paste0(prefix, "_covariate.txt")
    write.table(cbind(covid=rownames(covariates_file_name), covariates_file_name),
                row.names = F, col.names = T, quote = F, sep = "\t",
                file = cov.file)
    if(verbose){
      dim <- dim(covariates_file_name)
      cat("Created covariate file from data.frame for ",
          dim[1], " genes by ",
          dim[2], " individuals at:\n", cov.file, "\n")
    }
    covariates_file_name <- cov.file
  } else if (is.matrix(covariates_file_name)){
    cov.file <- paste0(prefix, "_covariate.txt")
    write.table(cbind(covid=rownames(covariates_file_name), covariates_file_name),
                row.names = F, col.names = T, quote = F, sep = "\t",
                file = cov.file)
    if(verbose){
      dim <- dim(covariates_file_name)
      cat("Created covariate file from matrix for ",
          dim[1], " genes by ",
          dim[2], " individuals at:\n", cov.file, "\n")
    }
    covariates_file_name <- cov.file
  } else if (is.vector(covariates_file_name) && is.numeric(covariates_file_name)) {
    cov.file <- paste0(prefix, "_covariate.txt")
    write.table(rbind(c("covid", paste0("sample", 1:length(covariates_file_name))),
                      c("cov1", covariates_file_name)),
                row.names = F, col.names = F, quote = F, sep = "\t",
                file = cov.file)
    if(verbose){
      dim <- length(covariates_file_name)
      cat("Created covariate file for a single gene for", dim,
          "individuals at:\n", cov.file, "\n")
    }
    covariates_file_name <- cov.file
  } else if (length(covariates_file_name)==0) {
    if(verbose){
      cat("No covariates used in this analysis\n")
    }
  } else {
    error("Something went wrong with importing your covariate data.")
  }

  ## Settings
  
  # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  
  
  # Only associations significant at this level will be saved
  if (compute.all) {
    pvOutputThreshold = 1
  } else {
    pvOutputThreshold = threshold
  }
  
  # Error covariance matrix
  # Set to numeric() for identity.
  errorCovariance = numeric();
  # errorCovariance = read.table("Sample_Data/errorCovariance.txt");
  
  
  ## Load genotype data
  
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(genotype_file_name);
  
  ## Load gene expression data
  
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  ## Load covariates
  
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
  
  ## Run the analysis
  
  me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = verbose,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = min.pv.by.genesnp,
    noFDRsaveMemory = save.memory);
  
  if(save.memory & load.qtls){
    # noFDRsaveMemory logical.
    # Set noFDRsaveMemory = TRUE to save significant gene-SNP pairs directly to
    # the output files, reduce memory footprint and skip FDR calculation.
    # The eQTLs are not recorded in the returned object if noFDRsaveMemory = TRUE.
    
    # let's add them manually
    cat("Read QTL_results from file\n", output_file_name, "\n")
    me$all$eqtls <- data.frame(fread(output_file_name));
    colnames(me$all$eqtls) <- c("snps", "gene", "beta", "statistic", "pvalue")
  }
  
  if(!(compute.all |save.memory)){
    unlink(output_file_name);
  }

  
  ## Results:
  
  cat('Analysis of ', me$all$ntests, 'QTLs done in: ', me$time.in.sec, ' seconds', '\n');
  cat('Detected eQTLs: ', me$all$neqtls, '\n');

  if(is.data.frame(snp.pos)){
    me$all$eqtls <- merge(snp.pos, me$all$eqtls,
                all.y=T)
    me$all$eqtls$chr <- as.character(me$all$eqtls$chr)
    me$all$eqtls$gene <- as.character(me$all$eqtls$gene)
    if(verbose){
      cat("SNP position annotation added to the QTL results.\n")
    }
  }
  
  if(what=="table" & min.pv.by.genesnp==FALSE){
    me <- data.frame(me$all$eqtls, stringsAsFactors = F)
  }
  
  return(me)
}



#' convenience function to get cumulative chromosome lengths for manhattan plots
#'
#' @family annotation functions
#' @title cumulative chromosome lengths
#' @param build genome build to use, options: 'b37' for GRCh37.p13,
#'        'b38' for GRCh38.p10, 'custom' (default 'b37')
#' @param df data.frame containing custom annotations (default NULL)
#' @param build.names names for colums with annotations for custom build,
#'        order: chromosome (chr), length (pos) (default NULL)
#' 
#' @author Ines Assum (2019-03-21)
#' @references
#' @export
add.chr.len.anno <- function(build="b37", df=NULL, build.names=NULL){
  chromosomeLength38 <- c(248956422, 242193529, 198295559, 190214555, 181538259,
                          170805979, 159345973, 145138636, 138394717, 133797422,
                          135086622, 133275309, 114364328, 107043718, 101991189,
                          90338345, 83257441, 80373285, 58617616, 64444167,
                          46709983, 50818468) # bp length from Ensembl GRCh38.p10
  chromosomeLength37 <- c(249250621, 243199373, 198022430, 191154276, 180915260,
                          171115067, 159138663, 146364022, 141213431, 135534747,
                          135006516, 133851895, 115169878, 107349540, 102531392,
                          90354753, 81195210, 78077248, 59128983, 63025520,
                          48129895, 51304566) # bp length from Ensembl GRCh37.p13
  if (build == 'b37') {
    pos <- data.frame(chr=1:22,
                      chrLen=chromosomeLength37)
  } else if (build == 'b38') {
    pos <- data.frame(chr=1:22,
                      chrLen=chromosomeLength38)
  } else if (build == 'custom') {
    library(data.table)
    df[, c("chr", "pos")] <- df[, build.names]
    pos <- as.data.table(df[, c("chr", "pos")])
    pos <- pos[, chrLen := max(pos), by = chr]
    pos <- pos[!duplicated(pos$chr), ]
    pos <- pos[order(pos$chr), ]
  } else {
    stop("Build not supported")
  }
  pos$cumSum2 <- cumsum(pos$chrLen)
  pos$cumSum <- pos$cumSum2-pos$chrLen
  pos$tick <- pos$cumSum + (pos$chrLen)/2
  return(pos)
}


#' convenience function to create ggplot2 manhattan plots of MatrixEQTL objects
#'
#' @family eqtl functions
#' @title make manhattan plots of MatrixEQTL results
#' @param me MatrixEQTL results as list or data.frame
#' @param build genome build to use, options: "b37", "b38", "custom"
#'        (default "b37")
#' @param build.names names for colums with annotations for custom build,
#'        order: chr, length (default NULL)
#' @param snp.pos data.frame with SNP annotations if not in me, SNP by cols:
#'        colnames(snp.pos) <- c("snps", "chr", "pos")
#' @param vcf vcf file to extract SNP annotations from (default NULL)
#' @param trait character vector of gene/trait names to plot
#' 
#' @author Ines Assum (2019-03-21)
#' @references
#' @export
manhattan.qtl <- function(me, build="b37", snp.pos=NULL, trait=NULL, vcf=NULL, build.names=NULL) {
  require(ggplot2)
  require(data.table)

  if (!is.data.frame(me)){
    me <- data.frame(me$all$eqtls,
                     stringsAsFactors = F)
  }
  
  if(!is.null(trait)){
    me <- me[me$gene %in% trait, ]
  }
  if(!is.null(snp.pos)){
    me <- merge(snp.pos, me,
                all.y=T)
  }
  pos <- add.chr.len.anno(build, df=me, build.names)
  me <- merge(me, pos,
              all.x = T)
  me$cumPos <- me$pos + me$cumSum
  qtl.plot <- ggplot(me, aes(x=cumPos, y=-log10(pvalue), col = chr)) +
    geom_point() +
    theme_bw() +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22),
                       guide=FALSE) +
    scale_x_continuous(label = pos$chr, breaks = pos$tick)
  return(qtl.plot)
}


#' Basic interface to obtain minimal P-values per gene from MatrixEQTL passing
#' expression and covariate data in matrices, and genotypes as a file. The
#' function uses \code{\link{eqtl}}.
#'
#' @family eqtl functions
#' @title basic eQTL interface
#' @param expr expression values in a matrix ngene x nsample
#' @param covar covariate values in a matrix ncovar x nsample
#' @param genotype_file_name filename of a tab separated file with genotype
#'        dosage values (usually between 0 and 2) in a matrix nsnp x nsample
#' @param gene.position data.frame with gene positions in the first 4 columns
#'        named "chrom", "start", "end", "gene_id"
#' @param snp.pos data.frame with SNP positions in the first three columns
#'        named "snp_id", "chrom", "snp_pos"
#' @param id character indicating an id used for tempfile creation
#' @return the same output as the \code{\link{eqtl}} function
#' @references
#' @export
eqtl.min.p <- function(expr, ...) {

  ## report the min(p) for each gene
  eqtls = eqtl.run(expr, ..., what="object")
  min.p = eqtls$cis$min.pv.gene

  ## make sure the order is the same as the input expression matrix
  min.p = min.p[rownames(expr)]
  
  return(min.p)
}



#' Extension to the waitForJobs function of the BatchJobs package which shows
#' some strange behaviour when waiting for jobs (database locked)
#' so we need to make it extra failsafe.
#'
#' @family eqtl functions
#' @title basic eQTL interface
#' @param reg BatchJobs registry
#' @param waittime time to wait before updating job status
#' @param nretry number of time to retry getting the job status before throwing #'        an error
myWaitForJobs <- function(reg, waittime=3, nretry=100) {
  success = FALSE
  while (nretry > 0 && !success) {
    status = tryCatch({
      while (TRUE) {
        status = showStatus(reg)
        if (status$done + status$expired == status$n) {
          cat("done\n")
          return(list(success=TRUE, nretry=nretry))
        }
        Sys.sleep(waittime)
      }
      return(list(success=FALSE, nretry=nretry))
    }, error=function(e) {
      cat("Error while waiting for jobs:\n")
      print(e)
      cat("\nnumber of retries left: ", nretry - 1, "\n")
      Sys.sleep(waittime + runif(1, 0, 3))
      return(list(success=FALSE, nretry=nretry - 1))
    })
    success = status$success
    nretry = status$nretry
    cat("success after the tryCatch block:", success, "\n")
    cat("nretry after the tryCatch block:", nretry, "\n")
  }
  

  if (!success) {
    err.msg = paste("Error during batch processing in registry")
    save(envir=sys.frame(), list=ls(envir=sys.frame()), file=file.path(dir, "error_image.RData"))
    stop(err.msg)
  }
}

#' Run eQTL permutation analysis using BatchJobs to distribute jobs on a HPC.
#' This function uses \code{\link{eqtl.run}} function.
#'
#' @family eqtl HPC functions
#' @title identification of eGenes
#' @param expr expression values in a matrix ngene x nsample
#' @param covar covariate values in a matrix ncovar x nsample
#' @param gene.position data.frame with gene positions in the first 4 columns
#'        named "chrom", "start", "end", "gene_id"
#' @param snp.pos data.frame with SNP positions in the first three columns
#'        named "snp_id", "chrom", "snp_pos"
#' @param genotype_file_name filename of a tab separated file with genotype
#'        dosage values (usually between 0 and 2) in a matrix nsnp x nsample
#' @param dir work directory where the BatchJobs registry is stored to. This
#'        path must be accessible by / visible to all compute nodes
#' @param min.perm minimal number of permutations
#' @param max.perm maximal number of permutations
#' @param seed random seed
#' @param exit.criterion number of times a more extreme test statistic has to
#'        be observed before a gene is not permuted further
#' @param actual.min.p if the actual minimal P-value per gene has already been
#'        computed it can be passed
#' @return a table with minimal P-values per gene and the empirical P-values
#' @references GTEx Consortium, Ardlie, K. G., Wright, F. A., & Dermitzakis, E. T. (2015). The Genotype-Tissue Expression (GTEx) pilot analysis: multitissue gene regulation in humans. Science, 348(6235), 648–660. \url{http://doi.org/10.1126/science.1262110}
#' @export
find.eGenes <- function(expr, covar, gene.position, snp.pos, genotype_file_name, dir, min.perm=1000, max.perm=10000, seed=0, exit.criterion=15, actual.min.p=NULL) {
  require(BatchJobs)
  ## initilaize the random generator
  set.seed(seed)
  
  ## check if we need to determine the actual min.p
  if (is.null(actual.min.p)) {
    actual.min.p = eqtl.min.p(expr, covar, genotype_file_name, gene.position, snp.pos, "actual")
  }
  
  
  ## submit batches of min.perm permutation runs to the cluster
  ## evaluate which genes go to the next round
  
  rfun <- function(pcount, new.order, expr, covar, genotype_file_name, gene.position, snp.pos) {
    source("R/eqtl_lib.R")
    id = paste("permutation", pcount, "___", sep="")

    o = new.order[pcount,]
    pexpr = expr[,o]
    pcovar = covar[,o]
    colnames(pexpr) = colnames(expr)
    colnames(pcovar) = colnames(covar)

    min.p = eqtl.min.p(pexpr, pcovar, genotype_file_name, gene.position, snp.pos, id)
    return(min.p)
  }

  
  genes.to.permute = 1:nrow(expr)
  pcount = 0
  count = rep(0, nrow(expr))
  nperm.per.gene = rep(0, nrow(expr))
  
  while (pcount < max.perm && length(genes.to.permute) > 0) {

    cat("\n\n\n\npermutation run", pcount, "analysing", length(genes.to.permute), "genes\n")
    
    ## permute sample labels
    new.order = t(sapply(1:min.perm, function(i) sample(1:ncol(expr), ncol(expr))))

    ## submit jobs
    rdir = file.path(dir, "batchjobs")
    reg = makeRegistry(id="permutations", seed=seed, file.dir=rdir, packages=c("MatrixEQTL", "data.table"))

    resources = list(memory="8G", queue="standard", time="4:00:00", longrun="False")
    more.args = list(new.order=new.order, expr=expr[genes.to.permute,,drop=F], covar=covar, genotype_file_name=genotype_file_name, gene.position=gene.position, snp.pos=snp.pos)
    
    batchMap(reg, rfun, 1:min.perm, more.args=more.args)
    
    submitJobs(reg, resources=resources, job.delay=function(n, i) runif(1, 0, 0.5))

    ## there is a problem with concurrent acess to the sqlite database
    ## so we need to wrap this in a loop with a number of retries
    myWaitForJobs(reg, waittime=30, nretry=100)
    
    
    ## reduce matrix is concatenating results row wise
    min.p = reduceResultsMatrix(reg)
    
    ## remove the registry
    system(paste("rm -rf", rdir))

    ## count for each gene how many times the permuted minp was smaller
    smaller = min.p <= rep(actual.min.p[genes.to.permute], each=nrow(min.p))
    count[genes.to.permute] = count[genes.to.permute] + colSums(smaller)

    ## increment the counter
    pcount = pcount + min.perm

    ## remember for each gene how many permutations were done
    nperm.per.gene[genes.to.permute] = pcount

    ## check which genes need to go to the next round
    genes.to.permute = which(count < exit.criterion)
    save(envir=sys.frame(), list=ls(envir=sys.frame()), file=file.path(dir, "current_image.RData"))
  }

  ## return the results
  tab = data.frame(actual.min.p, count, nperm.per.gene, empirical.p=count/nperm.per.gene)
  return(tab)
}


#' Run eQTL permutation analysis using BatchJobs to distribute jobs on a HPC.
#' This function uses \code{\link{eqtl.run}} function.
#'
#' @family eqtl HPC functions
#' @title identification of eSNPs
#' @param expr expression values in a matrix ngene x nsample
#' @param covar covariate values in a matrix ncovar x nsample
#' @param gene.position data.frame with gene positions in the first 4 columns
#'        named "chrom", "start", "end", "gene_id"
#' @param snp.pos data.frame with SNP positions in the first three columns
#'        named "snp_id", "chrom", "snp_pos"
#' @param genotype_file_name filename of a tab separated file with genotype
#'        dosage values (usually between 0 and 2) in a matrix nsnp x nsample
#' @param dir work directory where the BatchJobs registry is stored to. This
#'        path must be accessible by / visible to all compute nodes
#' @param min.perm minimal number of permutations
#' @param max.perm maximal number of permutations
#' @param seed random seed
#' @param exit.criterion number of times a more extreme test statistic has to
#'        be observed before a gene is not permuted further
#' @param actual.eqtls if the actual eQTLs have already been computed they can
#'        be passed
#' @param blocksize the number of permutations collected at the same time. This
#'        parameter tunes the balance between speed and memory usage. The more
#'        permutation runs are collected in the same block, the faster, but
#'        also more memory consuming.
#' @return a table with eQTL P-values and the empirical P-values for eSNPs
#' @references GTEx Consortium, Ardlie, K. G., Wright, F. A., & Dermitzakis, E. T. (2015). The Genotype-Tissue Expression (GTEx) pilot analysis: multitissue gene regulation in humans. Science, 348(6235), 648–660. \url{http://doi.org/10.1126/science.1262110}
#' @export
find.eSNPs <- function(expr, covar, gene.position, snp.pos, genotype_file_name, dir, threshold, min.perm=1000, max.perm=10000, seed=0, exit.criterion=15, actual.eqtls=NULL, blocksize=10) {

  ## as input we need the permutation threshold which is the empirical min(p)
  ## that corresponds to the FDR threshold across genes

  ## then we compute gene wise nominal p-values that correspond to the
  ## empirical threshold
  
  require(BatchJobs)
  ## initilaize the random generator
  set.seed(seed)
  
  ## check if we need to determine the actual min.p
  if (is.null(actual.eqtls)) {
    actual.eqtls= eqtl.run(expr, covar, genotype_file_name, gene.position, snp.pos, "actual")
  }
  
  
  ## submit batches of min.perm permutation runs to the cluster
  ## evaluate which genes go to the next round
  
  rfun <- function(pcount, new.order, expr, covar, genotype_file_name, gene.position, snp.pos) {
    source("R/eqtl_lib.R")
    id = paste("permutation", pcount, "___", sep="")

    o = new.order[pcount,]
    pexpr = expr[,o]
    pcovar = covar[,o]
    colnames(pexpr) = colnames(expr)
    colnames(pcovar) = colnames(covar)

    min.p = eqtl.min.p(pexpr, pcovar, genotype_file_name, gene.position, snp.pos, id)
    return(min.p)
  }


  ## we do not load all min.perm vectors at once
  ## instead we will aggregate by counting if the permuted min(p) value
  ## for a gene was smaller than the actual p value of a pair
  
  afun <- function(aggr, job, res, actual.pvalue, gene) {
    
    ## from the BatchJobs man package (reduceResults)
    ## Here, 'job' is the
    ## current job descriptor (see 'Job'), 'result' is the current
    ## result object and 'aggr' are the so far aggregated results.
    ## When using 'reduceResults', your function should add the
    ## stuff you want to have from 'job' and 'result' to 'aggr' and
    ## return that

    ## res is the min(p) per gene
    minp = res[gene]
    
    ## aggr is just the sum of permutations that was smaller
    aggr = aggr + as.numeric(minp <= actual.pvalue)

    return(aggr)
  }
  
  pairs.to.permute = 1:nrow(actual.eqtls)
  genes.to.permute = 1:nrow(expr)
  pcount = 0
  count = rep(0, nrow(actual.eqtls))
  nperm.per.gene = rep(0, nrow(actual.eqtls))
  
  while (pcount < max.perm && length(genes.to.permute) > 0) {

    cat("\n\n\n\npermutation run", pcount, "analysing", length(genes.to.permute), "genes\n")
    
    ## permute sample labels
    new.order = t(sapply(1:min.perm, function(i) sample(1:ncol(expr), ncol(expr))))

    ## submit jobs
    rdir = file.path(dir, "batchjobs")
    reg = makeRegistry(id="permutations", seed=seed, file.dir=rdir, packages=c("MatrixEQTL", "data.table"))

    resources = list(memory="8G", queue="standard", time="4:00:00", longrun="False")
    more.args = list(new.order=new.order, expr=expr[genes.to.permute,,drop=F], covar=covar, genotype_file_name=genotype_file_name, gene.position=gene.position, snp.pos=snp.pos)
    
    batchMap(reg, rfun, 1:min.perm, more.args=more.args)
    
    submitJobs(reg, resources=resources, job.delay=function(n, i) runif(1, 0, 0.5))

    ## there is a problem with concurrent acess to the sqlite database
    ## so we need to wrap this in a loop with a number of retries
    myWaitForJobs(reg, waittime=30, nretry=100)
    

    ## somehow reducing results one by one is very slow, so we do it block
    ## wise not to waste too much memory

    gene = actual.eqtls[pairs.to.permute, "gene"]
    actual.pvalue = actual.eqtls[pairs.to.permute, "p-value"]

    ## reshape to the size of the block
    actual.pvalue = matrix(rep(actual.pvalue, each=blocksize), nrow=blocksize)

    ids = findDone(reg)
    nsuccess = length(ids)
    while (length(ids) > 0) {
      ## get the block of results
      min.p = reduceResultsMatrix(reg, ids=head(ids, blocksize), progressbar=FALSE)

      ## reshape to the size of the eqtl results
      min.p = min.p[,gene,drop=F]

      ## check if the block is of full size
      if (nrow(min.p) < blocksize) {
        ## if not reduce the size of the actual pvalues matrix
        actual.pvalue = actual.pvalue[1:nrow(min.p),]
      }
      
      ## count
      smaller = min.p <= actual.pvalue
      count[pairs.to.permute] = count[pairs.to.permute] + colSums(smaller)

      ## remove the job ids that were just processed
      ids = tail(ids, -blocksize)
    }
    
    ## remove the registry
    system(paste("rm -rf", rdir))
    
    ## increment the counter
    pcount = pcount + nsuccess

    ## remember for each gene how many permutations were done
    nperm.per.gene[genes.to.permute] = pcount

    ## check which genes need to go to the next round
    pairs.to.permute = which(count < exit.criterion)
    genes.to.permute = which(rownames(expr) %in% actual.eqtls[pairs.to.permute,"gene"])
    
    save(envir=sys.frame(1), list=ls(envir=sys.frame(1)), file=file.path(dir, "current_image.RData"))
  }

  ## return the results
  nperm.per.gene = nperm.per.gene[match(actual.eqtls[,"gene"], rownames(expr))]
  tab = data.frame(as.data.frame(actual.eqtls), count, nperm.per.gene, empirical.p=count/nperm.per.gene)
  return(tab)
}



#' Get egenes from a esnp table
#'
#' @family eqtl HPC functions
#' @param esnps an esnps table obtained from \code{\link{find.esnps}}
#' @param fdr false discovery rate threshold
#' @return a table of egenes with columns for minimal P-values and emprical
#' P-values
#' @export
get.egenes.from.esnps <- function(esnps, fdr=0.05) {
  require(qvalue)
  
  ## get the gene wise min(p)
  min.p = tapply(1:nrow(esnps), esnps[,"gene"], function(idx) idx[which.min(esnps[idx,"p.value"])])
  egenes = esnps[min.p,c("gene", "p.value", "chrom", "start", "end", "gene_name", "gene_type", "count", "nperm.per.gene", "empirical.p")]
  
  ## compute the fdr from the empirical p values
  qval = qvalue(egenes[,"empirical.p"])
  egenes = cbind(egenes, qvalue=qval$qvalues)

  max.p = max(egenes[egenes[,"qvalue"] < fdr, "empirical.p"])
  
  ## determine the gene specific p value threshold that corresponds to an fdr 5%
  threshold = tapply(1:nrow(esnps), esnps[,"gene"], function(idx) {
    sig = idx[esnps[idx,"empirical.p"] < max.p]
    return(max(c(0, egenes[sig, "p.value"])))
  })
  
  egenes = cbind(egenes, gene.specific.threshold=threshold[egenes[,"gene"]])
  return(egenes)
}

#' Quantile - quantile plot against the uniform distribution
#'
#' @param q P-values
#' @export
qquniform <- function(p) {
  p = sort(p)
  expected = (1:length(p)) / length(p)
  plot(-log10(expected), -log10(p))
}

#' Estimate latent factors using the PEER method
#' @param k the number of latent factors
#' @param expr expression values in a matrix ngene x nsample
#' @param covar measured covariate values in a matrix  nsample x ncovar.
#'        If covar is NULL, no measured covariates will be used.
#' @return the latent factors in a  nsample x k matrix
get.peer.factors <- function(k, expr, covar) {
  require(peer)
  model = PEER()
  PEER_setPhenoMean(model, t(expr))
  PEER_setNk(model, k)
  if (!is.null(covar)) {
    PEER_setCovariates(model, as.matrix(covar))
  }
  PEER_update(model)
  factors = PEER_getX(model)
  weights = PEER_getW(model)
  return(factors)
}


#' Quantile normalization
#'
#' @param x ngenes x nsamples matrix to be normalized
#' @return quantile normalized matrix
#' @export
normalize.quantile <- function(x) {
  x = as.matrix(x)
  o = apply(x, 2, order)
  xsort = x
  for (col in 1:ncol(x)) {
    xsort[,col] = x[o[,col], col]
  }
  means = apply(xsort, 1, mean)
  normalized = matrix(NA, nrow=nrow(x), ncol=ncol(x), dimnames=dimnames(x))
  for (col in 1:ncol(x)) {
    normalized[o[,col], col] = means
  }
  return(normalized)
}
