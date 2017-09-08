### filter the KO according to p-value
t2d <- read.table("~/git-code/R/lmm/t2d_fit2.txt", sep = "\t")
#t2d[is.na(t2d)] <- -1
c <- c()
for(i in 1:ncol(t2d)){
    indx <- which(!is.na(t2d[,i]) & t2d[,i] < 0.00001 & t2d[,i] >0)
    if(length(indx) > 0){
        c <- c(c,i)
    }
}
t2d <- t2d[,c]

### filter the bacteria according to p-value
v <- c()
for(i in 1:nrow(t2d)){
    indx <- which(!is.na(t2d[i,]) & t2d[i,] < 0.00001 & t2d[i,] >0)
    if(length(indx) > 0){
        v <- c(v,i)
    }
}
t2d <- t2d[v,]

c <- c()
control <- read.table("~/git-code/R/lmm/control_fit2.txt", sep = "\t")
for(i in 1:ncol(control)){
    indx <- which(!is.na(control[,i]) & control[,i] < 0.00001 & control[,i] >0)
    if(length(indx) > 0){
        c <- c(c,i)
    }
}
control <- control[,c]
### filter the bacteria according to p-value
v <- c()
for(i in 1:nrow(control)){
    indx <- which(!is.na(control[i,]) & control[i,] < 0.00001 & control[i,] >0)
    if(length(indx) > 0){
        v <- c(v,i)
    }
}
control <- control[v,]

#intersect(colnames(t2d),colnames(control))
# KO map to pathway
library(KEGGREST)
result <- c()
k <- c()
for(ko in colnames(t2d)){
    k <- c(k,ko)
    if(length(k)>= 10){
        result <- c(result,names(keggGet(k)[[1]]$PATHWAY))
        k <- c()
    }
}
if(length(k) >0){
    result <- c(result,names(keggGet(k)[[1]]$PATHWAY)) 
}
result <- as.factor(result)


control_result <- c()
k <- c()
for(ko in colnames(control)){
    k <- c(k,ko)
    if(length(k)>= 10){
        control_result <- c(control_result,names(keggGet(k)[[1]]$PATHWAY))
        k <- c()
    }
}
if(length(k) >0){
    control_result <- c(control_result,names(keggGet(k)[[1]]$PATHWAY)) 
}
control_result <- as.factor(control_result)


