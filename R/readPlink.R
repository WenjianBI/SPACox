#' SaddlePoint Approximation implementation of a surival analysis (Plink input)
#'
#' A fast and accurate method for a genome-wide survival analysis on a large-scale dataset.
#' @param obj.null an R object returned from function SPACox_Null_Model()
#' @param plink.file character, represents the prefix of PLINK input file.
#' @param output.file character, represents the prefix of output file.
#' @param memory.chunk a numeric value (default: 4, unit=Gb) to specify how much memory is used to store genotype matrix from plink files.
#' @param min.maf a numeric value (default: 0.0001) to specify the cutoff of the minimal MAF. Any SNP with MAF < cutoff will be excluded from the analysis.
#' @param Cutoff a numeric value (Default: 2) to specify the standard deviation cutoff to be used.
#'               If the test statistic lies within the standard deviation cutoff, its p value is calculated based on a normal distribution approximation,
#'               otherwise, its p value is calculated based on a saddlepoint approximation.
#' @param impute.method a character string (default: "fixed") to specify the method to impute missing genotypes.
#'                      "fixed" imputes missing genotypes (NA) by assigning the mean genotype value (i.e. 2p where p is MAF).
#' @param missing.cutoff a numeric value (default: 0.15) to specify the cutoff of the missing rates.
#'                       Any variant with missing rate higher than this cutoff will be excluded from the analysis.
#' @param CovAdj.cutoff a numeric value (default: 5e-5). If the p-value is less than this cutoff, then we would use an additional technic to adjust for covariates.
#' @details To run SPACox, the following two steps are required:
#' \itemize{
#'   \item Step 1. Use function SPACox_Null_Model() to fit a null Cox model.
#'   \item Step 2: Use function SPACox() or SPACox.plink() to calculate p value for each genetic variant.
#' }
#'
#' SPACox uses a hybrid strategy with both saddlepoint approximation and normal distribution approximation.
#' Generally speaking, saddlepoint approximation is more accurate than, but a little slower than, the traditional normal distribution approximation.
#' Hence, when the score statistic is close to 0 (i.e. p-values are not small), we use the normal distribution approximation.
#' And when the score statistic is far away from 0 (i.e. p-values are small), we use the saddlepoint approximation.
#' Argument 'Cutoff' is to specify the standard deviation cutoff.
#'
#' To calibrate the score statistics, SPACox uses martingale residuals which are calculated via R package survival.
#' All extentions (such as strata, ties, left-censoring) supported by package survival could also be used in SPACox.
#' Time-varying covariates are also supported by splitting each subject into several observations.
#' Simulation studies and real data analyses indicate that SPACox works well if one subject corresponds to 2~3 observations.
#' While, if there are more than 4 observations for each subject, SPACox has not been fully evaluated and the results should be carefully intepreted.
#'
#' Sometimes, the order of subjects between phenotype data and genotype data are different, which could lead to some errors.
#' To avoid that, we ask users to specify the IDs of both phenotype data (pIDs) and genotype data (gIDs) when fitting the null model.
#' Users are responsible to check the consistency between pIDs and formula, and the consistency between gIDs and Geno.mtx(plink.file).
#'
#' @return the function outputs to a file (output.file) with the following columns
#' \item{markerID}{marker IDs}
#' \item{MAF}{Minor allele frequencies}
#' \item{missing.rate}{Missing rates}
#' \item{p.value.spa}{p value (recommanded) from a saddlepoint approximation.}
#' \item{p.value.norm}{p value from a normal distribution approximation.}
#' \item{Stat}{score statistics}
#' \item{Var}{estimated variances of the score statistics}
#' \item{z}{z values corresponding to the score statistics}
#' @examples
#' # Simulation phenotype and genotype
#' N = 1000
#' fam.file = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.fam", package = "SPACox")
#' fam.data = read.table(fam.file)
#' IDs = fam.data$V2
#' Phen.mtx = data.frame(ID = IDs,
#'                       event=rbinom(N,1,0.5),
#'                       time=runif(N),
#'                       Cov1=rnorm(N),
#'                       Cov2=rbinom(N,1,0.5))
#' plink.file = gsub("-ext.fam","-ext",fam.file)
#' output.file = gsub("-ext.fam","-ext.output",fam.file)
#'
#' # Attach the survival package so that we can use its function Surv()
#' library(survival)
#' obj.null = SPACox_Null_Model(Surv(time,event)~Cov1+Cov2, data=Phen.mtx,
#'                              pIDs=Phen.mtx$ID, gIDs=fam.data$V2)
#'
#' ## output is written in output.file
#' SPACox.plink(obj.null, plink.file, output.file)
#' @export
#' @import seqminer
#' @import data.table
SPACox.plink = function(obj.null,
                        plink.file,
                        output.file,
                        memory.chunk = 4,
                        Cutoff = 2,
                        impute.method = "fixed",
                        missing.cutoff = 0.15,
                        min.maf = 0.0001,
                        CovAdj.cutoff = 5e-5,
                        G.model = "Add")
{
  ## check plink files input
  bim.file = paste0(plink.file, ".bim")
  bed.file = paste0(plink.file, ".bed")
  fam.file = paste0(plink.file, ".fam")

  if(!file.exists(bim.file)) stop("Could not find paste0(plink.file,'.bim')")
  if(!file.exists(bed.file)) stop("Could not find paste0(plink.file,'.bed')")
  if(!file.exists(fam.file)) stop("Could not find paste0(plink.file,'.fam')")
  if(file.exists(output.file)) stop("'output.file' existed. Please give a different 'output.file' or remove the existing 'output.file'.")

  fam.data = read.table(fam.file)
  bim.data = read.table(bim.file)

  N = nrow(fam.data)
  M = nrow(bim.data)

  print(paste0("Totally ", M, " markers in plink files."))
  if(any(obj.null$gIDs != fam.data$V2))
    stop("any(obj.null$gIDs != fam.data$V2): when fitting null model, please give gIDs as the same order as in the plink files.")

  M.chunk = floor(memory.chunk * 1e9 / 4 / N)  # number of markers in each chunk
  n.chunk = ceiling(M / M.chunk)

  print(paste0("Split all markers into ", n.chunk, " chunks."))
  print(paste0("Each chunk includes less than ", M.chunk, " markers."))

  for(i in 1:n.chunk){
    if(i == n.chunk){
      markerIndex = ((n.chunk-1)*M.chunk+1):M;
    }else{
      markerIndex = 1:M.chunk + (i-1) * M.chunk;
    }
    print(paste0("Analyzing chunk ",i,"/",n.chunk,"."))
    Geno.mtx = seqminer::readPlinkToMatrixByIndex(plink.file, 1:N, markerIndex)
    colnames(Geno.mtx) = bim.data$V2[markerIndex]

    output = SPACox(obj.null, Geno.mtx, Cutoff, impute.method, missing.cutoff,
                    min.maf, CovAdj.cutoff, G.model)

    output = cbind(rownames(output), output)
    colnames(output)[1] = "markerID"

    if(i == 1){
      data.table::fwrite(output, output.file, sep = "\t", append = F, row.names = F, col.names = T)
    }else{
      data.table::fwrite(output, output.file, sep = "\t", append = T, row.names = F, col.names = F)
    }
  }

  print("Analysis Complete.")
  print(Sys.time())

  # return(output)  # output is stored in outFile
}
