source("~/git-code/R/lmm/read_tables.R")
source("~/git-code/R/lmm/construct_matrix.R")
source("~/git-code/R/lmm/lmm_fit.R")

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
write.table(x = t2d, file = "t2d_fit2.txt", sep = "\t")
control <- lmm_fit(ptr_data = ptr_data$control,kegg_data = kegg_data$control,result_maxtrix = related_matrix)
write.table(x = control, file = "control_fit2.txt", sep = "\t")



