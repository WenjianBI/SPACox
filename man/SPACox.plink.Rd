% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readPlink.R
\name{SPACox.plink}
\alias{SPACox.plink}
\title{SaddlePoint Approximation implementation of a surival analysis (Plink input)}
\usage{
SPACox.plink(
  obj.null,
  plink.file,
  output.file,
  memory.chunk = 4,
  Cutoff = 2,
  impute.method = "fixed",
  missing.cutoff = 0.15,
  min.maf = 1e-04,
  CovAdj.cutoff = 5e-05,
  G.model = "Add"
)
}
\arguments{
\item{obj.null}{an R object returned from function SPACox_Null_Model()}

\item{plink.file}{character, represents the prefix of PLINK input file.}

\item{output.file}{character, represents the prefix of output file.}

\item{memory.chunk}{a numeric value (default: 4, unit=Gb) to specify how much memory is used to store genotype matrix from plink files.}

\item{Cutoff}{a numeric value (Default: 2) to specify the standard deviation cutoff to be used.
If the test statistic lies within the standard deviation cutoff, its p value is calculated based on a normal distribution approximation,
otherwise, its p value is calculated based on a saddlepoint approximation.}

\item{impute.method}{a character string (default: "fixed") to specify the method to impute missing genotypes.
"fixed" imputes missing genotypes (NA) by assigning the mean genotype value (i.e. 2p where p is MAF).}

\item{missing.cutoff}{a numeric value (default: 0.15) to specify the cutoff of the missing rates.
Any variant with missing rate higher than this cutoff will be excluded from the analysis.}

\item{min.maf}{a numeric value (default: 0.0001) to specify the cutoff of the minimal MAF. Any SNP with MAF < cutoff will be excluded from the analysis.}

\item{CovAdj.cutoff}{a numeric value (default: 5e-5). If the p-value is less than this cutoff, then we would use an additional technic to adjust for covariates.}
}
\value{
the function outputs to a file (output.file) with the following columns
\item{markerID}{marker IDs}
\item{MAF}{Minor allele frequencies}
\item{missing.rate}{Missing rates}
\item{p.value.spa}{p value (recommanded) from a saddlepoint approximation.}
\item{p.value.norm}{p value from a normal distribution approximation.}
\item{Stat}{score statistics}
\item{Var}{estimated variances of the score statistics}
\item{z}{z values corresponding to the score statistics}
}
\description{
A fast and accurate method for a genome-wide survival analysis on a large-scale dataset.
}
\details{
To run SPACox, the following two steps are required:
\itemize{
  \item Step 1. Use function SPACox_Null_Model() to fit a null Cox model.
  \item Step 2: Use function SPACox() or SPACox.plink() to calculate p value for each genetic variant.
}

SPACox uses a hybrid strategy with both saddlepoint approximation and normal distribution approximation.
Generally speaking, saddlepoint approximation is more accurate than, but a little slower than, the traditional normal distribution approximation.
Hence, when the score statistic is close to 0 (i.e. p-values are not small), we use the normal distribution approximation.
And when the score statistic is far away from 0 (i.e. p-values are small), we use the saddlepoint approximation.
Argument 'Cutoff' is to specify the standard deviation cutoff.

To calibrate the score statistics, SPACox uses martingale residuals which are calculated via R package survival.
All extentions (such as strata, ties, left-censoring) supported by package survival could also be used in SPACox.
Time-varying covariates are also supported by splitting each subject into several observations.
Simulation studies and real data analyses indicate that SPACox works well if one subject corresponds to 2~3 observations.
While, if there are more than 4 observations for each subject, SPACox has not been fully evaluated and the results should be carefully intepreted.

Sometimes, the order of subjects between phenotype data and genotype data are different, which could lead to some errors.
To avoid that, we ask users to specify the IDs of both phenotype data (pIDs) and genotype data (gIDs) when fitting the null model.
Users are responsible to check the consistency between pIDs and formula, and the consistency between gIDs and Geno.mtx(plink.file).
}
\examples{
# Simulation phenotype and genotype
N = 1000
fam.file = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.fam", package = "SPACox")
fam.data = read.table(fam.file)
IDs = fam.data$V2
Phen.mtx = data.frame(ID = IDs,
                      event=rbinom(N,1,0.5),
                      time=runif(N),
                      Cov1=rnorm(N),
                      Cov2=rbinom(N,1,0.5))
plink.file = gsub("-ext.fam","-ext",fam.file)
output.file = gsub("-ext.fam","-ext.output",fam.file)

# Attach the survival package so that we can use its function Surv()
library(survival)
obj.null = SPACox_Null_Model(Surv(time,event)~Cov1+Cov2, data=Phen.mtx,
                             pIDs=Phen.mtx$ID, gIDs=fam.data$V2)

## output is written in output.file
SPACox.plink(obj.null, plink.file, output.file)
}
