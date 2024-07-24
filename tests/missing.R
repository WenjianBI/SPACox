library(SPACox)
# Simulation phenotype and genotype
N = 10000
nSNP = 1000
MAF = 0.1
Phen.mtx = data.frame(ID = paste0("IID-",1:N),
                      event=rbinom(N,1,0.5),
                      time=runif(N),
                      Cov1=rnorm(N),
                      Cov2=rbinom(N,1,0.5))
Geno.mtx = matrix(rbinom(N*nSNP,2,MAF),N,nSNP)

# NOTE: The row and column names of genotype matrix are required.
rownames(Geno.mtx) = paste0("IID-",1:N)
colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
Geno.mtx[1:10,1]=NA   # please use NA for missing genotype
Geno.mtx[,3]=NA   # set an entire column NA

# Attach the survival package so that we can use its function Surv()
library(survival)
obj.null = SPACox_Null_Model(Surv(time,event)~Cov1+Cov2, data=Phen.mtx,
                             pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
SPACox.res = SPACox(obj.null, Geno.mtx)

# we recommand using column of 'p.value.spa' to associate genotype with time-to-event phenotypes
print(head(SPACox.res))
