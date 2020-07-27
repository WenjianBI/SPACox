#' Fits a NULL model for SPACox
#'
#' Fits a null Cox proportional hazards model and then calculates the empirical cumulant generation function (CGF) of the martingale residuals
#' @param formula a formula to be passed to function coxph(). For more details, please refer to package survival.
#' @param data a data.frame in which to interpret the variables named in the formula
#' @param pIDs a character vector of subject IDs. NOTE: its order should be the same as the subjects order in the formula.
#' @param gIDs a character vector of subject IDs. NOTE: its order should be the same as the subjects order of the Geno.mtx (i.e. the input of the function SPACox()).
#' @param range a two-element numeric vector (default: c(-100,100)) to specify the domain of the empirical CGF.
#' @param length.out a positive integer (default: 9999) for empirical CGF. Larger length.out corresponds to longer calculation time and more accurate estimated empirical CGF.
#' @param ... Other arguments passed to function coxph(). For more details, please refer to package survival.
#' @return an object with a class of "SPACox_NULL_Model".
#' @examples
#' Please check help(SPACox) for a simulated example.
#' @export
#' @import survival
SPACox_Null_Model = function(formula,
                             data=NULL,
                             pIDs=NULL,
                             gIDs=NULL,
                             range=c(-100,100),
                             length.out = 10000,
                             ...)
{
  Call = match.call()

  ### Fit a Cox model
  obj.coxph = coxph(formula, data=data, x=T, ...)

  ### Check input arguments
  obj.check = check_input(pIDs, gIDs, obj.coxph, range)
  p2g = obj.check$p2g
  pIDs = obj.check$pIDs

  ### Get the covariate matrix to adjust for genotype
  mresid = obj.coxph$residuals
  Cova = obj.coxph$x

  X = cbind(1, Cova)
  X.invXX = X %*% solve(t(X)%*%X)
  tX = t(X)

  ### calculate empirical CGF for martingale residuals
  idx0 = qcauchy(1:length.out/(length.out+1))
  idx1 = idx0 * max(range) / max(idx0)

  cumul = NULL
  print("Start calculating empirical CGF for martingale residuals...")
  c = 0
  for(i in idx1){
    c = c+1
    t = i
    e_resid = exp(mresid*t)
    M0 = mean(e_resid)
    M1 = mean(mresid*e_resid)
    M2 = mean(mresid^2*e_resid)
    K0 = log(M0)
    K1 = M1/M0
    K2 = (M0*M2-M1^2)/M0^2
    cumul = rbind(cumul, c(t, K0, K1, K2))
    if(c %% 1000 == 0) print(paste0("Complete ",c,"/",length.out,"."))
  }

  K_org_emp = approxfun(cumul[,1], cumul[,2], rule=2)
  K_1_emp = approxfun(cumul[,1], cumul[,3], rule=2)
  K_2_emp = approxfun(cumul[,1], cumul[,4], rule=2)

  var.resid = var(mresid)

  re=list(resid = mresid,
          var.resid = var.resid,
          K_org_emp = K_org_emp,
          K_1_emp = K_1_emp,
          K_2_emp = K_2_emp,
          Call = Call,
          obj.coxph = obj.coxph,
          tX = tX,
          X.invXX = X.invXX,
          p2g = p2g,
          gIDs = gIDs,
          pIDs = pIDs)

  class(re)<-"SPACox_NULL_Model"
  return(re)
}

#' SaddlePoint Approximation implementation of a surival analysis
#'
#' A fast and accurate method for a genome-wide survival analysis on a large-scale dataset.
#' @param obj.null an R object returned from function SPACox_Null_Model()
#' @param Geno.mtx a numeric genotype matrix with each row as an individual and each column as a genetic variant.
#'                 Column names of genetic variations and row names of subject IDs are required.
#'                 Missng genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
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
#'   \item Step 2: Use function SPACox() to calculate p value for each genetic variant.
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
#' Users are responsible to check the consistency between pIDs and formula, and the consistency between gIDs and Geno.mtx.
#'
#' @return an R matrix with the following columns
#' \item{MAF}{Minor allele frequencies}
#' \item{missing.rate}{Missing rates}
#' \item{p.value.spa}{p value (recommanded) from a saddlepoint approximation.}
#' \item{p.value.norm}{p value from a normal distribution approximation.}
#' \item{Stat}{score statistics}
#' \item{Var}{estimated variances of the score statistics}
#' \item{z}{z values corresponding to the score statistics}
#' @examples
#' # Simulation phenotype and genotype
#' N = 10000
#' nSNP = 1000
#' MAF = 0.1
#' Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                       event=rbinom(N,1,0.5),
#'                       time=runif(N),
#'                       Cov1=rnorm(N),
#'                       Cov2=rbinom(N,1,0.5))
#' Geno.mtx = matrix(rbinom(N*nSNP,2,MAF),N,nSNP)
#'
#' # NOTE: The row and column names of genotype matrix are required.
#' rownames(Geno.mtx) = paste0("IID-",1:N)
#' colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
#' Geno.mtx[1:10,1]=NA   # please use NA for missing genotype
#'
#' # Attach the survival package so that we can use its function Surv()
#' library(survival)
#' obj.null = SPACox_Null_Model(Surv(time,event)~Cov1+Cov2, data=Phen.mtx,
#'                              pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
#' SPACox.res = SPACox(obj.null, Geno.mtx)
#'
#' # we recommand using column of 'p.value.spa' to associate genotype with time-to-event phenotypes
#' head(SPACox.res)
#'
#' ## missing data in response/indicator variables is also supported. Please do not remove pIDs of subjects with missing data, the program will do it.
#' Phen.mtx$event[2] = NA
#' Phen.mtx$Cov1[5] = NA
#' obj.null = SPACox_Null_Model(Surv(time,event)~Cov1+Cov2, data=Phen.mtx,
#'                              pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
#' SPACox.res = SPACox(obj.null, Geno.mtx)
#'
#' # The below is an example code to use survival package
#' coxph(Surv(time,event)~Cov1+Cov2+Geno.mtx[,1], data=Phen.mtx)
#' @export
SPACox = function(obj.null,
                  Geno.mtx,
                  Cutoff = 2,
                  impute.method = "fixed",
                  missing.cutoff = 0.15,
                  min.maf = 0.0001,
                  CovAdj.cutoff = 5e-5,
                  G.model = "Add")
{
  ## check input
  par.list = list(pwd=getwd(),
                  sessionInfo=sessionInfo(),
                  Cutoff=Cutoff,
                  impute.method=impute.method,
                  missing.cutoff=missing.cutoff,
                  min.maf=min.maf,
                  CovAdj.cutoff=CovAdj.cutoff,
                  G.model=G.model)

  check_input1(obj.null, Geno.mtx, par.list)
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))

  ### Prepare the main output data frame
  n.Geno = ncol(Geno.mtx)
  output = matrix(NA, n.Geno, 7)
  colnames(output) = c("MAF","missing.rate","p.value.spa","p.value.norm","Stat","Var","z")
  rownames(output) = colnames(Geno.mtx)

  ### Start analysis
  print("Start Analyzing...")
  print(Sys.time())

  # Cycle for genotype matrix
  for(i in 1:n.Geno){

    g = Geno.mtx[,i]
    output.one.SNP = SPACox.one.SNP(g,
                                    obj.null,
                                    Cutoff,
                                    impute.method,
                                    missing.cutoff,
                                    min.maf,
                                    CovAdj.cutoff,
                                    G.model)
    output[i,] = output.one.SNP
  }

  print("Analysis Complete.")
  print(Sys.time())
  return(output)
}

#' SaddlePoint Approximation implementation of Cox regression surival analysis (One-SNP-version)
#'
#' One-SNP-version SPACox function. This function is to facilitate users that prefer reading and analyzing genotype line-by-line.
#' @param g a numeric genotype vector. Missing genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param others the same as function SPACox. NOTE that we do not check subject order in this one-snp-version !!!
#' @return the same as function SPACox.
#' @export
SPACox.one.SNP = function(g,
                          obj.null,
                          Cutoff = 2,
                          impute.method = "fixed",
                          missing.cutoff = 0.15,
                          min.maf = 0.0001,
                          CovAdj.cutoff = 5e-5,
                          G.model = "Add")
{
  g[g==-9]=NA  # since we add plink input
  ## calculate MAF and update genotype vector
  MAF = mean(g, na.rm=T)/2
  N = length(g)
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/N

  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }

  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
  }

  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)

  if(MAF < min.maf)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA))

  if(!is.null(obj.null$p2g))
    g = g[obj.null$p2g]

  ## Score statistic
  S = sum(g * obj.null$resid)

  ## estimated variance without adjusting for covariates
  G1 = g - 2*MAF   # centered genotype (such that mean=0)
  S.var1 = obj.null$var.resid * sum(G1^2)
  z1 = S/sqrt(S.var1)

  if(abs(z1) < Cutoff){
    pval.norm = pnorm(abs(z1), lower.tail = FALSE)*2
    return(c(MAF, missing.rate, pval.norm, pval.norm, S, S.var1, z1))
  }

  N1set = which(g!=0)  # position of non-zero genotypes
  N0 = N-length(N1set)

  G1norm = G1/sqrt(S.var1)  # normalized genotype (such that sd=1)

  G1N1 = G1norm[N1set]
  G1N0 = -2*MAF/sqrt(S.var1)   # all subjects with g=0 share the same normlized genotype, this is to reduce computation time

  pval1 = GetProb_SPA(obj.null, G1N1, G1N0, N1set, N0, abs(z1), lower.tail = FALSE)
  pval2 = GetProb_SPA(obj.null, G1N1, G1N0, N1set, N0, -abs(z1), lower.tail = TRUE)
  pval = pval1 + pval2

  if(pval[1] > CovAdj.cutoff)
    return(c(MAF, missing.rate, pval, S, S.var1, z1))

  ## estimated variance after adjusting for covariates

  G2 = g - obj.null$X.invXX %*% (obj.null$tX[,N1set,drop=F] %*% g[N1set])
  S.var2 = obj.null$var.resid * sum(G2^2)
  z2 = S/sqrt(S.var2)

  G2norm = G2/sqrt(S.var2)

  N1set = 1:N
  N0 = 0
  G2N1 = G2norm
  G2N0 = 0   # since N0=0, this value actually does not matter

  pval1 = GetProb_SPA(obj.null, G2N1, G2N0, N1set, N0, abs(z2), lower.tail = FALSE)
  pval2 = GetProb_SPA(obj.null, G2N1, G2N0, N1set, N0, -abs(z2), lower.tail = TRUE)
  pval = pval1 + pval2

  return(c(MAF, missing.rate, pval, S, S.var2, z2))
}


GetProb_SPA = function(obj.null, G2NB, G2NA, NBset, N0, q2, lower.tail){

  out = uniroot(K1_adj, c(-20,20), extendInt = "upX",
                G2NB=G2NB, G2NA=G2NA, NBset=NBset,
                N0=N0, q2=q2, obj.null=obj.null)
  zeta = out$root

  k1 = K_org(zeta,  G2NB=G2NB, G2NA=G2NA, NBset=NBset, N0=N0, obj.null=obj.null)
  k2 = K2(zeta,  G2NB=G2NB, G2NA=G2NA, NBset=NBset, N0=N0, obj.null=obj.null)

  temp1 = zeta * q2 - k1

  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}

  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)
  pval.norm = pnorm(q2, lower.tail = lower.tail)

  re = c(pval, pval.norm)
  return(re)
}


K_org = function(t, G2NB, G2NA, NBset, N0, obj.null){

  n.t = length(t)
  out = rep(0,n.t)
  for(i in 1:n.t){
    t1 = t[i]
    t2NA = t1*G2NA
    t2NB = t1*G2NB
    out[i] = N0*obj.null$K_org_emp(t2NA) + sum(obj.null$K_org_emp(t2NB))
  }
  return(out)
}

K1_adj = function(t, G2NB, G2NA, NBset, N0, q2, obj.null)
{
  n.t = length(t)
  out = rep(0,n.t)

  for(i in 1:n.t){
    t1 = t[i]
    t2NA = t1*G2NA
    t2NB = t1*G2NB
    out[i] = N0*G2NA*obj.null$K_1_emp(t2NA) + sum(G2NB*obj.null$K_1_emp(t2NB)) - q2
  }
  return(out)
}

K2 = function(t, G2NB, G2NA, NBset, N0, obj.null)
{
  n.t = length(t)
  out = rep(0,n.t)

  for(i in 1:n.t){
    t1 = t[i]
    t2NA = t1*G2NA
    t2NB = t1*G2NB
    out[i] = N0*G2NA^2*obj.null$K_2_emp(t2NA) + sum(G2NB^2*obj.null$K_2_emp(t2NB))
  }
  return(out)
}


check_input = function(pIDs, gIDs, obj.coxph, range)
{
  if(is.null(pIDs) & is.null(gIDs))
    stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")

  pIDs = as.character(pIDs)
  gIDs = as.character(gIDs)
  if(!is.null(obj.coxph$na.action)){
    posNA = c(obj.coxph$na.action)
    if(any(posNA > length(pIDs)))
      stop("Number of input data is larger than length(pIDs).")
    pIDsNA = pIDs[posNA]

    print(paste0("Due to missing data in response/indicators, ",length(posNA)," entries are removed from analysis."))
    print("If concerned about the power loss, we suggest users impute data first and then use SPACox package.")
    print(head(cbind(posNA=posNA, pIDsNA=pIDsNA)))

    pIDs = pIDs[-1*posNA]  # remove IDs with missing data
  }

  if(any(!is.element(pIDs, gIDs)))
    stop("All elements in pIDs should be also in gIDs.")

  if(anyDuplicated(gIDs)!=0)
    stop("Argument 'gIDs' should not have a duplicated element.")

  if(range[2]!=-1*range[1])
    stop("range[2] should be -1*range[1]")

  mresid = obj.coxph$residuals

  if(length(mresid)!=length(pIDs))
    stop("length(mresid)!=length(pIDs) where mresid is the martingale residuals from coxph() in survival package.")

  p2g = NULL
  if(length(pIDs)!=length(gIDs)){
    p2g = match(pIDs, gIDs)
  }else{
    if(any(pIDs != gIDs))
      p2g = match(pIDs, gIDs)
  }

  return(list(p2g=p2g,pIDs=pIDs))
}

check_input1 = function(obj.null, Geno.mtx, par.list)
{
  if(class(obj.null)!="SPACox_NULL_Model")
    stop("obj.null should be a returned outcome from SPACox_Null_Model()")

  if(any(obj.null$gIDs != rownames(Geno.mtx))) stop("gIDs should be the same as rownames(Geno.mtx).")
  if(is.null(rownames(Geno.mtx))) stop("Row names of 'Geno.mtx' should be given.")
  if(is.null(colnames(Geno.mtx))) stop("Column names of 'Geno.mtx' should be given.")
  if(!is.numeric(Geno.mtx)|!is.matrix(Geno.mtx)) stop("Input 'Geno.mtx' should be a numeric matrix.")

  if(!is.numeric(par.list$min.maf)|par.list$min.maf<0|par.list$min.maf>0.5) stop("Argument 'min.maf' should be a numeric value >= 0 and <= 0.5.")
  if(!is.numeric(par.list$Cutoff)|par.list$Cutoff<0) stop("Argument 'Cutoff' should be a numeric value >= 0.")
  # if(!is.element(par.list$impute.method,c("none","bestguess","random","fixed"))) stop("Argument 'impute.method' should be 'none', 'bestguess', 'random' or 'fixed'.")
  if(!is.element(par.list$impute.method,c("fixed"))) stop("Argument 'impute.method' should be 'fixed'.")
  if(!is.numeric(par.list$missing.cutoff)|par.list$missing.cutoff<0|par.list$missing.cutoff>1) stop("Argument 'missing.cutoff' should be a numeric value between 0 and 1.")
  if(!is.element(par.list$G.model,c("Add","Dom","Rec"))) stop("Argument 'G.model' should be 'Add', 'Dom' or 'Rec'.")
}

