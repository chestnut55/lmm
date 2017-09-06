### PCA analysis
source("code/R/lmm/read_tables.R")
t2d <- read_ptr_data()$T2D
t2d <- t2d[,which(!apply(t2d,2,FUN = function(x){all(x == 0)}))]
t2d <- t2d[which(!apply(t2d,1,FUN = function(x){all(x == 0)})),]

###psych package
library(psych)
fa.parallel(t(t2d), fa = "pc", n.iter = 1000, show.legend = FALSE, main = "screen plot with parallel analysis")
pc <- principal(t(t2d), nfactors = 4)

###another PCA package
t2d.pr <- princomp(t(t2d),cor = TRUE)
screeplot(t2d.pr, type = "lines")
summary(t2d.pr, loadings = TRUE)

kegg_t2d_data <- read_kegg_data()$T2D

### remove the row and column all is zeros
kegg_t2d_data <- kegg_t2d_data[,which(!apply(kegg_t2d_data,2,FUN = function(x){all(x == 0)}))]
kegg_t2d_data <- kegg_t2d_data[which(!apply(kegg_t2d_data,1,FUN = function(x){all(x == 0)})),]

### the psych package analysis PCA failed because of exception
#library(psych)
#fa.parallel(t(kegg_t2d_data), fa = "pc", n.iter = 100, show.legend = FALSE, main = "screen plot with parallel analysis")

###exception:"'princomp'只能在单位比变量多的情况下使用"
kegg_t2d.pr <- princomp(t(kegg_t2d_data),cor = TRUE)
screeplot(kegg_t2d.pr, type = "lines")
summary(kegg_t2d.pr, loadings = TRUE)