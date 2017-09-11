source("~/git-code/R/lmm/read_tables.R")
source("~/git-code/R/lmm/construct_matrix.R")
source("~/git-code/R/lmm/lmm_fit.R")
source("~/git-code/R/lmm/analysis.R")

### organism and species
organism_species_data <- read_organism_species()
ptr_data <- read_ptr_data()
kegg_data <- read_kegg_data()

#sample_01_matrix <- ramdom_sample_matrix(organism_species_data,cbind(kegg_data$T2D, kegg_data$control))
related_matrix <- retrieve_related_matrix()
related_matrix <- related_matrix[rownames(ptr_data$control),rownames(kegg_data$control)]

library("picaplot")
library("confeti")

#t2d <- fit(ptr_data = ptr_data$T2D,kegg_data = kegg_data$T2D,organism_species = organism_species_data,result_maxtrix = related_matrix)
#write.table(x = t2d, file = "t2d_fit.txt", sep = "\t")
#control <- fit(ptr_data = ptr_data$control,kegg_data = kegg_data$control,organism_species = organism_species_data,result_maxtrix = related_matrix)
#write.table(x = control, file = "control_fit.txt", sep = "\t")

t2d <- lmm_fit(ptr_data = ptr_data$T2D,kegg_data = kegg_data$T2D,result_maxtrix = related_matrix)
write.table(x = t2d, file = "~/git-code/R/lmm/data/generated/fit/t2d_fit.txt", sep = "\t")
control <- lmm_fit(ptr_data = ptr_data$control,kegg_data = kegg_data$control,result_maxtrix = related_matrix)
write.table(x = control, file = "~/git-code/R/lmm/data/generated/fit/control_fit.txt", sep = "\t")


############result analysis###########################
library(KEGGREST)
### filter the KO according to p-value
P_VALUE <- 0.00001

t2d <- read.table("~/git-code/R/lmm/data/generated/fit/t2d_fit.txt", sep = "\t")
control <- read.table("~/git-code/R/lmm/data/generated/fit/control_fit.txt", sep = "\t")

t2d <- remove_no_signicant(t2d, P_VALUE)
control <- remove_no_signicant(control, P_VALUE)

#intersect(colnames(t2d),colnames(control))
# KO map to pathway
result <- map2pathway(t2d)
control_result <- map2pathway(control)

t2d_gene_set <- convet2genes("~/git-code/R/lmm/data/generated/gsea/t2d/",t2d, P_VALUE)
control_gene_set <- convet2genes("~/git-code/R/lmmdata/generated/gsea/control/",control, P_VALUE)

#########################################################



