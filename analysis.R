### filter the KO according to p-value
t2d <- read.table("~/code/R/lmm/data/t2d_fit.txt", sep = "\t")
t2d[is.na(t2d)] <- -1
c <- c()
for(i in 1:ncol(t2d)){
    indx <- which(t2d[,i] < 0.001 & t2d[,i] >0)
    if(length(indx) > 0){
        c <- c(c,i)
    }
}
t2d <- t2d[,c]

### filter the bacteria according to p-value
v <- c()
for(i in 1:nrow(t2d)){
    indx <- which(t2d[i,] < 0.001 & t2d[i,] >0)
    if(length(indx) > 0){
        v <- c(v,i)
    }
}
t2d <- t2d[v,]

c <- c()
control <- read.table("~/code/R/lmm/data/control_fit.txt", sep = "\t")
control[is.na(control)] <- -1
for(i in 1:ncol(control)){
    indx <- which(control[,i] < 0.001 & control[,i] >0)
    if(length(indx) > 0){
        c <- c(c,i)
    }
}
control <- control[,c]
### filter the bacteria according to p-value
v <- c()
for(i in 1:nrow(control)){
    indx <- which(control[i,] < 0.001 & control[i,] >0)
    if(length(indx) > 0){
        v <- c(v,i)
    }
}
control <- control[v,]


