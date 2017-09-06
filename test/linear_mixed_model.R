###############数据预处理 2017.07.06#########################
bgi_SRP008047 <-
  read.csv("~/code/R/linear_mixed_model/data/EBI/mg_project_summary_20170706_SRP008047.csv",
    header = TRUE,
    check.names = FALSE,
    sep = ","
  )

#bgi_SRP008047 <- bgi_SRP008047[,1:ncol(bgi_SRP008047)-1]


bgi_SRP011011 <-
  read.csv(file = "~/code/R/linear_mixed_model/data/EBI/mg_project_summary_20170706_SRP011011.csv",
           header = TRUE,
           check.names = FALSE,
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
             sep = "\t")

### change the PTR run id with sample name
ptr_data_stageII_T2D <- ptr_data[,as.vector(bgi_data_stageII_T2D$`Run ID`)]
ptr_data_stageII_control <- ptr_data[,as.vector(bgi_data_stageII_control$`Run ID`)]

### KEGG --- T2D and control group for KEGG
kegg_data <-
  read.table(file = "~/code/R/linear_mixed_model/data/KEGG/KEGG.StageII.relative_abun.txt",
             header = TRUE,
             check.names = FALSE,
             sep = "\t")

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


########KEGG KO##############
library(XML)
library(sqldf)
library(KEGGREST)

bacts <- as.vector(ptr_data$X)
orgs <- data.frame(keggList("organism"))
orgs <- sqldf("select o1.organism,o1.species from orgs o1 INNER JOIN ptr_data p1 ON o1.species == p1.X")

ko_ids <- as.vector(kegg_data$K)

result_maxtrix <- data.frame(matrix(ncol = length(ko_ids), nrow = length(bacts)))
colnames(result_maxtrix) <- ko_ids
rownames(result_maxtrix) <- bacts
#ko_ids <-c("K00001","K00002","K00003","K00004","K00005","K00006","K00007","K00008","K00009","K00010")
#ko_ids <-c("K00005")
ko <- c()
for(ko_id in ko_ids){
  ko <- c(ko,ko_id)
  tmp <- as.character(ko)
  if(length(ko) == 5){
    ####### ko length is 10
    url <- paste("http://www.kegg.jp/kegg-bin/view_ortholog_table?against=bacteria&orthology=",paste(ko,collapse="+"),"&mode=complete",collapse = "",sep = "")
    print(url)
    orthology <- readHTMLTable(url)
    print("readHTMLTable finished")
    if(length(orthology) >0 && !is.null(orthology[[1]])){
      len <- length(orthology[[1]])
      for(i in 4:len){
        orthology_matrix <- data.frame(orthology[[1]][3],orthology[[1]][i],stringsAsFactors = FALSE)
        colnames(orthology_matrix) <- c("organism","ko")
        
        result <- sqldf("select o1.organism as organism,o1.species as species from orgs o1 INNER JOIN orthology_matrix p1 ON o1.organism == p1.organism")
        if(length(result$organism) !=0 || length(result$species) !=0){
          rindex <- which(rownames(result_maxtrix) %in% result$species)
          cindex <- which(colnames(result_maxtrix)== tmp[i -3])
          result_maxtrix[rindex, cindex] <- 1
        }
      }
    }
    ko <- c() ## clean 
  }
}

if(length(ko) >0){
  ####### ko length is 10
  url <- paste("http://www.kegg.jp/kegg-bin/view_ortholog_table?against=bacteria&orthology=",paste(ko,collapse="+"),"&mode=complete",collapse = "",sep = "")
  print(url)
  orthology <- readHTMLTable(url)
  #print("readHTMLTable finished")
  if(length(orthology) >0 && !is.null(orthology[[1]])){
    len <- length(orthology[[1]])
    for(i in 4:len){
      orthology_matrix <- data.frame(orthology[[1]][3],orthology[[1]][i],stringsAsFactors = FALSE)
      colnames(orthology_matrix) <- c("organism","ko")
      
      result <- sqldf("select o1.organism as organism,o1.species as species from orgs o1 INNER JOIN orthology_matrix p1 ON o1.organism == p1.organism")
      if(length(result$organism) !=0 || length(result$species) !=0){
        rindex <- which(rownames(result_maxtrix) %in% result$species)
        cindex <- which(colnames(result_maxtrix) == tmp[i -3])
        result_maxtrix[rindex, cindex] <- 1
      }
    }
  }
  ko <- c() ## clean 
}

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
ica_result <- runICA(kegg_data)
ica_object <- ica_genotype_test(ica_result, t(ptr_data))
confeti_results <- get_similarity_mx(ica_object)
ica_object$confeti <- confeti_results

kegg_data <- t(kegg_data)
ptr_data <- t(ptr_data)
single_kegg <- as.numeric(kegg_data[1,])
single_ptr <- as.numeric(ptr_data[1,])
lmm_fit <- lrgpr(single_kegg~single_ptr, decomp <- svd(confeti_results$Kmx))