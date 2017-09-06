###############数据预处理 2017.07.06#########################
bgi_SRP008047 <-
  read.csv("~/code/R/linear_mixed_model/data/EBI/mg_project_summary_20170706_SRP008047.csv",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    sep = ","
  )

#bgi_SRP008047 <- bgi_SRP008047[,1:ncol(bgi_SRP008047)-1]


bgi_SRP011011 <-
  read.csv(file = "~/code/R/linear_mixed_model/data/EBI/mg_project_summary_20170706_SRP011011.csv",
           header = TRUE,
           check.names = FALSE,
           stringsAsFactors = FALSE,
           sep = ",")

#bgi_SRP011011 <- bgi_SRP011011[,1:ncol(bgi_SRP011011)-1]

bgi_data <- rbind(bgi_SRP008047, bgi_SRP011011)

#### type 2 diabetes samples, control group and disease group
stageII_samples <-   read.csv(file = "~/code/R/linear_mixed_model/data/T2D/stage II.csv",
                              header = TRUE,
                              check.names = FALSE,
                              sep = ",")
stageII_T2D_samples <- stageII_samples[stageII_samples$`Diabetic (Y or N)` == 'Y',] # T2D group samples
stageII_control_samples <- stageII_samples[stageII_samples$`Diabetic (Y or N)` == 'N',] # contro group samples

bgi_data_stageII_T2D <- bgi_data[bgi_data$`Sample Name` %in% as.vector(stageII_T2D_samples$`Sample ID`),]
bgi_data_stageII_control <- bgi_data[bgi_data$`Sample Name` %in% as.vector(stageII_control_samples$`Sample ID`),]

### PTR sample with run id instead of sample name, T2D and control group samples from PTR
ptr_data <-
  read.table(file = "~/code/R/linear_mixed_model/data/PTR/ptr_qin.txt",
             header = TRUE,
             stringsAsFactors = FALSE,
             sep = "\t")
rownames(ptr_data) <- as.vector(ptr_data[,1])
ptr_data <- ptr_data[,c(-1,-2)]


### change the PTR run id with sample name
ptr_data_stageII_T2D <- ptr_data[,as.vector(bgi_data_stageII_T2D$`Run ID`)]
ptr_data_stageII_control <- ptr_data[,as.vector(bgi_data_stageII_control$`Run ID`)]

### KEGG --- T2D and control group for KEGG
kegg_data <-
  read.table(file = "~/code/R/linear_mixed_model/data/KEGG/KEGG.StageII.relative_abun.txt",
             header = TRUE,
             check.names = FALSE,
             stringsAsFactors = FALSE,
             sep = "\t")
rownames(kegg_data) <- kegg_data[,1]
kegg_data <- kegg_data[,-1]
kegg_data_stageII_T2D <- kegg_data[,as.vector(bgi_data_stageII_T2D$`Sample Name`)]
kegg_data_stageII_control <- kegg_data[,as.vector(bgi_data_stageII_control$`Sample Name`)]

### change the PTR run id with sample name
library(data.table)
setnames(ptr_data_stageII_T2D, as.vector(bgi_data_stageII_T2D$`Run ID`), as.vector(bgi_data_stageII_T2D$`Sample Name`))
setnames(ptr_data_stageII_control, as.vector(bgi_data_stageII_control$`Run ID`), as.vector(bgi_data_stageII_control$`Sample Name`))

# order by sample name
ptr_data_stageII_T2D <- ptr_data_stageII_T2D[,c(sort(colnames(ptr_data_stageII_T2D)))] # T2D
kegg_data_stageII_T2D <- kegg_data_stageII_T2D[,sort(colnames(kegg_data_stageII_T2D))] # T2D
kegg_data_stageII_T2D <- kegg_data_stageII_T2D[rowSums(kegg_data_stageII_T2D) !=0,]

# set the NA value to column mean
ind <- which(is.na(ptr_data_stageII_T2D), arr.ind=TRUE)
ptr_data_stageII_T2D[ind] <- rowMeans(ptr_data_stageII_T2D,  na.rm = TRUE)[ind[,1]]

ptr_data_stageII_control <- ptr_data_stageII_control[,c(sort(colnames(ptr_data_stageII_control)))]  # control
kegg_data_stageII_control <- kegg_data_stageII_control[,sort(colnames(kegg_data_stageII_control))] # control
kegg_data_stageII_control <- kegg_data_stageII_control[rowSums(kegg_data_stageII_control) !=0,]
# set the NA value to column mean
ind <- which(is.na(ptr_data_stageII_control), arr.ind=TRUE)
ptr_data_stageII_control[ind] <- rowMeans(ptr_data_stageII_control,  na.rm = TRUE)[ind[,1]]
#ptr_str <- as.vector(ptr_data$X)
rownames(ptr_data_stageII_T2D) <- rownames(ptr_data)
rownames(ptr_data_stageII_control) <- rownames(ptr_data)
########KEGG KO##############
library(XML)
library(sqldf)
library(KEGGREST)
library(data.table)

#bacts <- as.vector(ptr_data$X)
#orgs <- data.frame(keggList("organism"))
#orgs <- sqldf("select o1.organism,o1.species from orgs o1 INNER JOIN ptr_data p1 ON o1.species == p1.X")

#ko_ids <- as.vector(kegg_data$K)
#row <- length(bacts)
#col <- length(ko_ids)

organism_species <-
  read.table(file = "~/code/R/linear_mixed_model/data/organism_species.csv",
             header = TRUE,
             check.names = FALSE,
             stringsAsFactors = FALSE,
             sep = ",")
### random generate the 0,1 matrix
random_values <- sample(c(0,1), replace=TRUE, size=nrow(organism_species) * nrow(kegg_data))
random_maxtrix <- matrix(random_values, nrow = nrow(organism_species), ncol = nrow(kegg_data), byrow = TRUE)
result_maxtrix <- data.frame(random_maxtrix)
colnames(result_maxtrix) <- colnames(kegg_data)
rownames(result_maxtrix) <- organism_species$organism

##################线性混合模型############
#======install dependencies================
picaplot_dependencies <- c("ggplot2", "knitr", "rmarkdown", "mclust")
install.packages(picaplot_dependencies)

source("http://bioconductor.org/biocLite.R")
biocLite("SummarizedExperiments")

#installing picaplot
install.packages("devtools")
library("devtools")
devtools::install_github("jinhyunju/picaplot")
devtools::install_github("jinhyunju/confeti")

#======install dependencies finished========

#setwd("~/code/R/linear_mixed_model")
#source("confeti_functions.R")
library("picaplot")
library("confeti")
#ica_result <- runICA(kegg_data)
#ica_object <- ica_genotype_test(ica_result, t(ptr_data))
#confeti_results <- get_similarity_mx(ica_object)
#ica_object$confeti <- confeti_results

K <- cor(ptr_data_stageII_T2D)
decomp <- eigen(K, symmetric=TRUE)
tmp  <- t(ptr_data_stageII_T2D)

### the final result for lmm
lmm_result_maxtrix <- data.frame(matrix(nrow = nrow(organism_species), ncol = nrow(kegg_data), byrow = TRUE)) # intercepet
rownames(lmm_result_maxtrix) <- organism_species$organism
colnames(lmm_result_maxtrix) <- rownames(kegg_data)



result_maxtrix_row_names <- rownames(result_maxtrix)
for(i in 1:nrow(kegg_data_stageII_T2D)){
  y <- t(kegg_data_stageII_T2D[i,])
  xx <- c()
  for(j in 1:40){
    if(result_maxtrix[j,i] == 1){
      xx <- c(xx,j)
    }
  }

  if(length(xx) >0){
    formula <- as.formula(paste("y~", paste(paste0(paste0("tmp[,",xx),"]"), collapse="+")))
    #print(formula)
    fit <- lrgpr(formula, decomp)
    p_values <- as.vector(fit$p.values)
    #print(p_values)
    if(length(p_values) > 0){
      p_values <- p_values[2:length(p_values)]
      lmm_result_maxtrix[xx,i] <- p_values
    }
  }
}


#kegg_data <- t(kegg_data)
#ptr_data <- t(ptr_data)
#single_kegg <- as.numeric(kegg_data[1,])
#single_ptr <- as.numeric(ptr_data[1,])
#lmm_fit <- lrgpr(single_kegg~single_ptr, decomp <- svd(confeti_results$Kmx))