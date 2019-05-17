## Why we revise the coxphf package
The original R package coxphf (version 1.13) is very slow, especially when analyzing interval time data whose max(start_time)>min(end_end). The following is an example  
 ```{r }
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

## How to use the new codes
