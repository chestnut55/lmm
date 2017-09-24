library("picaplot")
library("confeti")
library("KEGGREST")
library("stringr") ## str_count method
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

### filter the KO according to p-value
P_VALUE <- 0.0005

t2d <- read.table("~/git-code/R/lmm/data/generated/fit/t2d_fit.txt", sep = "\t")
control <- read.table("~/git-code/R/lmm/data/generated/fit/control_fit.txt", sep = "\t")

t2d <- remove_no_signicant(t2d, P_VALUE)
control <- remove_no_signicant(control, P_VALUE)

########################pathway and species map###################################################
t2d_pathway_result <- read.table("~/git-code/R/lmm/data/generated/t2d_pathway", sep = "\t")
control_pathway_result <- read.table("~/git-code/R/lmm/data/generated/control_pathway", sep = "\t")

t2d_count_species_in_pathway <- count_species_in_pathway(t2d_pathway_result)
control_count_species_in_pathway <- count_species_in_pathway(control_pathway_result)
remove(t2d_pathway_result)
remove(control_pathway_result)
#######################################end###################################################
####################################sampling#################################################
t2d_sampling_count <- sampling("~/git-code/R/lmm/data/generated/t2d_sampling", t2d_count_species_in_pathway)
control_sampling_count <- sampling("~/git-code/R/lmm/data/generated/control_sampling", control_count_species_in_pathway)
######################draw graphs for KO and species##############################
draw_bip_network(t2d, p_value = P_VALUE)
draw_bip_network(control, p_value = P_VALUE)

######################draw weighted graphs for pathway and species##############################
t2d <- read.table("~/git-code/R/lmm/data/generated/t2d_pathway", sep = "\t")
t2d[is.na(t2d)] <- 0
control <- read.table("~/git-code/R/lmm/data/generated/control_pathway", sep = "\t")
control[is.na(control)] <- 0
draw_weight_bip_network(t2d)
draw_weight_bip_network(control)