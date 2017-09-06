read_organism_species <- function(){
  organism_species <-
    read.table(file = "~/code/R/linear_mixed_model/data/organism_species.csv",
               header = TRUE,
               check.names = FALSE,
               stringsAsFactors = FALSE,
               sep = ",")
  
  return (organism_species)
}


read_stageII_sample <- function(){
  #### type 2 diabetes samples, control group and disease group
  stageII_samples <-   read.csv(file = "~/code/R/linear_mixed_model/data/T2D/stage II.csv",
                                header = TRUE,
                                check.names = FALSE,
                                sep = ",")
  data_result <- list()
  data_result$T2D <- stageII_samples[stageII_samples$`Diabetic (Y or N)` == 'Y',] # T2D group samples
  data_result$control <- stageII_samples[stageII_samples$`Diabetic (Y or N)` == 'N',] # contro group samples
  
  return (data_result)
}


read_bgi_data <- function(){
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
  
  data_result <- list()
  
  samples_data <- read_stageII_sample()
  stageII_T2D_samples <- samples_data$T2D
  stageII_control_samples <- samples_data$control
  
  
  data_result$T2D <- bgi_data[bgi_data$`Sample Name` %in% as.vector(stageII_T2D_samples$`Sample ID`),]
  data_result$control <- bgi_data[bgi_data$`Sample Name` %in% as.vector(stageII_control_samples$`Sample ID`),]
  
  return (data_result)
}



read_ptr_data <- function(){
  ### PTR sample with run id instead of sample name, T2D and control group samples from PTR
  data <-
    read.table(file = "~/code/R/linear_mixed_model/data/PTR/ptr_qin.txt",
               header = TRUE,
               stringsAsFactors = FALSE,
               sep = "\t")
  rownames(data) <- as.vector(data[,1])
  data <- data[,c(-1,-2)]
  
  bgi_data <- read_bgi_data()
  bgi_data_stageII_T2D <- bgi_data$T2D
  bgi_data_stageII_control <- bgi_data$control
  
  ### T2D and control group PTR data
  ptr_data_stageII_T2D <- data[,as.vector(bgi_data_stageII_T2D$`Run ID`)]
  ptr_data_stageII_control <- data[,as.vector(bgi_data_stageII_control$`Run ID`)]
  
  ### change the PTR run id with sample name
  library(data.table)
  setnames(ptr_data_stageII_T2D, as.vector(bgi_data_stageII_T2D$`Run ID`), as.vector(bgi_data_stageII_T2D$`Sample Name`))
  setnames(ptr_data_stageII_control, as.vector(bgi_data_stageII_control$`Run ID`), as.vector(bgi_data_stageII_control$`Sample Name`))
  
  # order by sample name
  ptr_data_stageII_T2D <- ptr_data_stageII_T2D[,c(sort(colnames(ptr_data_stageII_T2D)))] # T2D
  ptr_data_stageII_control <- ptr_data_stageII_control[,c(sort(colnames(ptr_data_stageII_control)))]  # control
  
  # set the NA value to column mean
  ind <- which(is.na(ptr_data_stageII_T2D), arr.ind=TRUE)
  ptr_data_stageII_T2D[ind] <- rowMeans(ptr_data_stageII_T2D,  na.rm = TRUE)[ind[,1]]
  
  ind <- which(is.na(ptr_data_stageII_control), arr.ind=TRUE)
  ptr_data_stageII_control[ind] <- rowMeans(ptr_data_stageII_control,  na.rm = TRUE)[ind[,1]]
  
  #rownames(ptr_data_stageII_T2D) <- rownames(data)
  #rownames(ptr_data_stageII_control) <- rownames(data)
  

  
  data_result <- list()
  data_result$T2D <- ptr_data_stageII_T2D
  data_result$control <- ptr_data_stageII_control
  
  return (data_result)
}

read_kegg_data <- function(){
  ### KEGG --- T2D and control group for KEGG
  kegg_data <-
    read.table(file = "~/code/R/linear_mixed_model/data/KEGG/KEGG.StageII.relative_abun.txt",
               header = TRUE,
               check.names = FALSE,
               stringsAsFactors = FALSE,
               sep = "\t")
  rownames(kegg_data) <- kegg_data[,1]
  kegg_data <- kegg_data[,-1]
  kegg_data <- kegg_data[rowSums(kegg_data) !=0,]
  
  bgi_data <- read_bgi_data()
  bgi_data_stageII_T2D <- bgi_data$T2D
  bgi_data_stageII_control <- bgi_data$control
  
  kegg_data_stageII_T2D <- kegg_data[,as.vector(bgi_data_stageII_T2D$`Sample Name`)]
  kegg_data_stageII_control <- kegg_data[,as.vector(bgi_data_stageII_control$`Sample Name`)]
  
  kegg_data_stageII_T2D <- kegg_data_stageII_T2D[,sort(colnames(kegg_data_stageII_T2D))] # T2D
  kegg_data_stageII_control <- kegg_data_stageII_control[,sort(colnames(kegg_data_stageII_control))] # control
  
  data_result <- list()
  data_result$T2D <- kegg_data_stageII_T2D
  data_result$control <- kegg_data_stageII_control
  
  return (data_result)
}