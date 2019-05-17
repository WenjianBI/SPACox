## How to use the new codes
1. download the three files to any directory
2. open R and then enter
```{r}
setwd(DIR)    # set working directory to the directory that stores the 3 files
source("coxphf.R")
source("decomposeSurv.r")
system("R CMD SHLIB coxphf_v3.f90")
dyn.load("coxphf_v3.so")
```
3. Now, you can follow the same usage as in R package coxphf to use its function (much faster).

## Why we revise the coxphf package
The original R package coxphf (version 1.13) is very slow, especially when analyzing interval time data whose max(start_time)>min(end_end). The following is an example  
```{r}
library(coxphf)
library(rbenchmark)
N = 3000
t1=runif(N)
t2=runif(N)

Data=data.frame(start=t1,end=t1+t2,
                event=rbinom(N,1,0.1),
                X1=rnorm(N),
                X2=rbinom(N,1,0.5))

a=benchmark(coxphf(Surv(start,end,event)~X1+X2,data=Data),replications = 10)
# a
#                                                     test replications elapsed relative user.self sys.self user.child sys.child
# 1 coxphf(Surv(start, end, event) ~ X1 + X2, data = Data)           10   24.77        1     24.69        0         NA        NA
```

## How we revise the code
1. We update the loop with a smarter order to subtract exp(risk) et al. of subjects whose start_time < current_time_in_loop. The update can decrease computation time from N^2 to 2N.
2. We add some if statement to avoid large amount of unnecessary calculation. This can save computation from 2N to 2K, where K is event number.
3. We use the code in section 'ifastmode=T' to replace the code in section 'ifastmode=F'.

