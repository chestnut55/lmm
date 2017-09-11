#### construct 0,1 matrix for bacteria PTR and kegg abundance
#### column: KEGG abundance
#### row: bacteria PTR
####
#### this method run only one time, call retrieve_related_matrix
####
#construct_related_matrix <- function(organism, kegg_ko){
#  result <- data.frame(matrix(nrow = length(organism),ncol = length(kegg_ko), byrow = TRUE), stringsAsFactors = FALSE)
#  for(ko_id in kegg_ko){
#   url <- paste0("http://www.kegg.jp/kegg-bin/view_ortholog_table?against=bacteria&orthology=",ko_id,"&mode=complete")
#    print(url)
#   orthology <- readHTMLTable(url)
#   if(length(orthology) >0 && !is.null(orthology[[1]])){
#     orthology_matrix <- data.frame(orthology[[1]][3],orthology[[1]][4],stringsAsFactors = FALSE)
#     colnames(orthology_matrix) <- c("organism","ko")
#     
#     orthology_matrix <- subset(orthology_matrix, orthology_matrix$organism %in% organism)
#     na.omit(orthology_matrix)
#     result[c(orthology_matrix$organism),ko_id] <- 1
#    }
#  }
#  return (result)
#}

construct_related_matrix <- function() {
  ### PTR sample with run id instead of sample name, T2D and control group samples from PTR
  ptr_data <-
    read.table(
      file = "~/git-code/R/lmm/data/PTR/ptr_qin.txt",
      header = TRUE,
      stringsAsFactors = FALSE,
      sep = "\t"
    )
  rownames(ptr_data) <- as.vector(ptr_data[, 1])
  ptr_data <- ptr_data[, c(-1, -2)]
  
  ### KEGG --- T2D and control group for KEGG
  kegg_data <-
    read.table(
      file = "~/git-code/R/lmm/data/KEGG/KEGG.StageII.relative_abun.txt",
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      sep = "\t"
    )
  rownames(kegg_data) <- kegg_data[, 1]
  kegg_data <- kegg_data[, -1]
  
  library(XML)
  library(sqldf)
  library(KEGGREST)
  
  
  result_maxtrix <-
  data.frame(matrix(ncol = length(rownames(kegg_data)), nrow = length(rownames(ptr_data))))
  colnames(result_maxtrix) <- rownames(kegg_data)
  rownames(result_maxtrix) <- rownames(ptr_data)
  orgs <-
    read.table(
      file = "~/git-code/R/lmm/data/organism_species.csv",
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      sep = ","
    )
  #ko_ids <-c("K00001","K00002","K00003","K00004","K00005","K00006","K00007","K00008","K00009","K00010")
  #ko_ids <-c("K00005")
    ko_ids <- c(colnames(result_maxtrix))
    ko <- c()
    for (ko_id in ko_ids) {
      ko <- c(ko, ko_id)
      tmp <- as.character(ko)
      if (length(ko) == 10) {
        ####### ko length is 10
        url <-
          paste(
            "http://www.kegg.jp/kegg-bin/view_ortholog_table?against=bacteria&orthology=",
            paste(ko, collapse = "+"),
            "&mode=complete",
            collapse = "",
            sep = ""
          )
        #print(url)
        orthology <- readHTMLTable(url)
        if (length(orthology) > 0 && !is.null(orthology[[1]])) {
          for (i in 4:length(orthology[[1]])) {
            orthology_matrix <-
              data.frame(orthology[[1]][3], orthology[[1]][i], stringsAsFactors = FALSE)
            colnames(orthology_matrix) <- c("organism", "ko")
            
            result <-
              sqldf(
                "select o1.organism as organism,o1.species as species from orgs o1 INNER JOIN orthology_matrix p1 ON o1.organism == p1.organism"
              )
            if (length(result$organism) != 0 ||
                length(result$species) != 0) {
              rindex <- which(rownames(result_maxtrix) %in% result$organism)
              cindex <- which(colnames(result_maxtrix) == tmp[i - 3])
              result_maxtrix[rindex, cindex] <- 1
            }
          }
        }
        ko <- c() ## clean
      }
    }
    
    if (length(ko) > 0) {
      ####### ko length is 10
      url <-
        paste(
          "http://www.kegg.jp/kegg-bin/view_ortholog_table?against=bacteria&orthology=",
          paste(ko, collapse = "+"),
          "&mode=complete",
          collapse = "",
          sep = ""
        )
      print(url)
      orthology <- readHTMLTable(url)
      #print("readHTMLTable finished")
      if (length(orthology) > 0 && !is.null(orthology[[1]])) {
        for (i in 4:length(orthology[[1]])) {
          orthology_matrix <-
            data.frame(orthology[[1]][3], orthology[[1]][i], stringsAsFactors = FALSE)
          colnames(orthology_matrix) <- c("organism", "ko")
          
          result <-
            sqldf(
              "select o1.organism as organism,o1.species as species from orgs o1 INNER JOIN orthology_matrix p1 ON o1.organism == p1.organism"
            )
          if (length(result$organism) != 0 ||
              length(result$species) != 0) {
            rindex <- which(rownames(result_maxtrix) %in% result$organism)
            cindex <- which(colnames(result_maxtrix) == tmp[i - 3])
            result_maxtrix[rindex, cindex] <- 1
          }
        }
      }
      ko <- c() ## clean
    }
    write.table(result_maxtrix, "result_matrix01.txt", sep = "\t")
}


### construct matrix through keggGet API
### m_matrix is species * ko matrix
construct_01_matrix <- function(m_matrix){
    all_correct_ko <- intersect(get_all_ko(), colnames(m_matrix))
    for(ko in all_correct_ko){
        genes_set <- keggGet(ko)[[1]]$GENES
        for(i in 1:length(genes_set)){
            ul <- unlist(strsplit(genes_set[i],":"))
            row_indx <- match(tolower(ul[1]), rownames(m_matrix))
            if(!is.na(row_indx)){
                col_indx <- match(ko, colnames(m_matrix))
                m_matrix[row_indx, col_indx] <- 1
                break
            }
        }
    }
    return (m_matrix)
}


retrieve_related_matrix <- function(){
  related_matrix <- read.table(file = "~/git-code/R/lmm/data/generated/result_01_matrix.txt", sep = "\t")
  return (related_matrix)
}

### generate random sample matrix for test purpose
ramdom_sample_matrix <- function(organism_species, kegg_data){
  ### random generate the 0,1 matrix
  random_values <- sample(c(0,1), replace=TRUE, size=nrow(organism_species) * nrow(kegg_data))
  random_maxtrix <- matrix(random_values, nrow = nrow(organism_species), ncol = nrow(kegg_data), byrow = TRUE)
  result_maxtrix <- data.frame(random_maxtrix)
  colnames(result_maxtrix) <- colnames(kegg_data)
  rownames(result_maxtrix) <- organism_species$organism
  
  return (result_maxtrix)
}

### get all ko , otherwise if the ko no exist, the methods keggGet will throw exception
###
get_all_ko <- function(){
    kk <- keggList("ko")
    result <- c()
    for(k in names(kk)){
        ko <- unlist(strsplit(k,":"))[2]
        result <- c(result, ko)
    }
    return(result)
}