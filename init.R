##############################################################################
##############################################################################
#########init the intern file to save time####################################
library("picaplot")
library("confeti")
library("KEGGREST")
### 1. read tables
source("~/git-code/R/lmm/read_tables.R")
### 2. construct matrix
source("~/git-code/R/lmm/construct_matrix.R")
### 3. lmm fit
source("~/git-code/R/lmm/lmm_fit.R")
### 4. result analysis
source("~/git-code/R/lmm/analysis.R")
### 5. draw graphs
source("~/git-code/R/lmm/graph.R")

########################PART 1. read tables########################
### organism and species
ptr_data <- read_ptr_data()
kegg_data <- read_kegg_data()


#### modify @2017.09.11 for construct 0-1 matrix
row_names <- rownames(ptr_data$control)
col_names <- rownames(kegg_data$control)
related_matrix <- data.frame(matrix(rep(0,length(row_names) * length(col_names)), nrow = length(row_names),ncol = length(col_names)))
colnames(related_matrix) <- col_names
rownames(related_matrix) <- row_names
related_matrix <- construct_01_matrix(related_matrix)
write.table(x = related_matrix, file = "~/git-code/R/lmm/data/generated/result_01_matrix.txt", sep = "\t")
#### modify @2017.09.11 end

########################PART 2. construct 0-1 matrix########################
#sample_01_matrix <- ramdom_sample_matrix(organism_species_data,cbind(kegg_data$T2D, kegg_data$control))
related_matrix <- retrieve_related_matrix()
related_matrix <- related_matrix[rownames(ptr_data$control),rownames(kegg_data$control)]

#t2d <- fit(ptr_data = ptr_data$T2D,kegg_data = kegg_data$T2D,organism_species = organism_species_data,result_maxtrix = related_matrix)
#write.table(x = t2d, file = "t2d_fit.txt", sep = "\t")
#control <- fit(ptr_data = ptr_data$control,kegg_data = kegg_data$control,organism_species = organism_species_data,result_maxtrix = related_matrix)
#write.table(x = control, file = "control_fit.txt", sep = "\t")
########################PART 3. fit the lmm model########################
t2d <- lmm_fit(ptr_data = ptr_data$T2D,kegg_data = kegg_data$T2D,result_maxtrix = related_matrix)
write.table(x = t2d, file = "~/git-code/R/lmm/data/generated/fit/t2d_fit.txt", sep = "\t")
control <- lmm_fit(ptr_data = ptr_data$control,kegg_data = kegg_data$control,result_maxtrix = related_matrix)
write.table(x = control, file = "~/git-code/R/lmm/data/generated/fit/control_fit.txt", sep = "\t")

########################PART 4. result analysis########################
### filter the KO according to p-value
P_VALUE <- 0.0005
t2d <- remove_no_signicant(t2d, P_VALUE)
control <- remove_no_signicant(control, P_VALUE)

########################ko and species map###################################################
# KO map to pathway , a little slowly
t2d_pathway_result <- map2pathway("~/git-code/R/lmm/data/generated/t2d_pathway",t2d, P_VALUE)
control_pathway_result <- map2pathway("~/git-code/R/lmm/data/generated/control_pathway",control, P_VALUE)

# get KEGG protein access number according to KO number
t2d_gene_set <- convert2genes("~/git-code/R/lmm/data/generated/gsea/t2d/",t2d, P_VALUE)
control_gene_set <- convert2genes("~/git-code/R/lmm/data/generated/gsea/control/",control, P_VALUE)

### convert to ncbi protein id
convert2ncbi_protein_id("~/git-code/R/lmm/data/generated/gsea/t2d/", 
                        "~/git-code/R/lmm/data/generated/protein/t2d/") # t2d
convert2ncbi_protein_id("~/git-code/R/lmm/data/generated/gsea/control/",
                        "~/git-code/R/lmm/data/generated/protein/control/") # control

t2d_ko_list <- get_significant_ko("~/git-code/R/lmm/data/generated/ko/t2d/",t2d, P_VALUE)
control_ko_list <- get_significant_ko("~/git-code/R/lmm/data/generated/ko/control/",control, P_VALUE)
