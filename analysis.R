### remove the row or column without any significant level
### data : t2d or control
remove_no_signicant <- function(data, P_VALUE){
    c <- c()
    for(i in 1:ncol(data)){
        indx <- which(!is.na(data[,i]) & data[,i] < P_VALUE & data[,i] >0)
        if(length(indx) > 0){
            c <- c(c,i)
        }
    }
    data <- data[,c]
    
    ### filter the bacteria according to p-value
    v <- c()
    for(i in 1:nrow(data)){
        indx <- which(!is.na(data[i,]) & data[i,] < P_VALUE & data[i,] >0)
        if(length(indx) > 0){
            v <- c(v,i)
        }
    }
    data <- data[v,]
    
    return (data)
}
### keggGet or keggFind all OK, keggFind is simple but a little slowly
convert2genes <- function(path, species_ko, P_VALUE){
    genes_set <- c()
    rownames <- rownames(species_ko)
    colnames <- colnames(species_ko)
    for(i in 1:nrow(species_ko)){
        indx <- which(!is.na(species_ko[i,]) & species_ko[i,] < P_VALUE & species_ko[i,] >0)
        row_name <- rownames[i]
        for(k in colnames[indx]){
            genes <- names(keggFind("genes",c(row_name,k)))
            if(length(genes) > 0){
                for(j in 1:length(genes)){
                    genes_set <- c(genes_set, genes[j])
                }   
            }
        }
        write.table(x = genes_set, file = paste0(path,row_name), sep = "\n", row.names = FALSE, col.names = FALSE)
        genes_set <- c()
    }
}

### map KO of species to pathway
map2pathway <- function(path, data, P_VALUE){
    result <- list()
   
    rownames <- rownames(data)
    colnames <- colnames(data)
    for(i in 1:nrow(data)){
        pathway_set <- c()
        indx <- which(!is.na(data[i,]) & data[i,] < P_VALUE & data[i,] >0)
        row_name <- rownames[i]
        for(k in colnames[indx]){
            pathway <- names(keggGet(k)[[1]]$PATHWAY)
            if(length(pathway) > 0){
                pathway_set <- c(pathway_set, pathway)
            }
        }
        if(length(pathway_set) >0){
            write.table(x = pathway_set, file = paste0(path,row_name), sep = "\n", row.names = FALSE, col.names = FALSE)
            result[[row_name]] <- table(pathway_set)
        }
    }
    return(result)
}

to_pathway_matrix <- function(filename, list){
    column_names <- c()
    row_names <- names(list)
    for(i in 1:length(row_names)){
         t <- list[[row_names[i]]]
         column_names <- c(column_names, names(t))
    }
    column_names <- unique(column_names)
    result <- data.frame(matrix(nrow = length(row_names), ncol = length(column_names)))
    colnames(result) <- column_names
    rownames(result) <- row_names
    
    for(i in 1:length(row_names)){
        df <- as.data.frame(list[[row_names[i]]])
        for(j in 1:nrow(df)){
            col_name <- as.character(df[j,1])
            result[c(row_names[i]), c(col_name)] <- df[j,2]
        }
    }
    
    result <- result[ , order(names(result))]
    write.table(x = result, file = filename, sep = "\t")
    return (result)
}

get_significant_ko <- function(write_path, data, p_value){
    result <- list()
    for(i in 1:ncol(data)){
       indx <- which(data[i,]<p_value & data[i,]>0)
       result[[rownames(data)[i]]] <- colnames(data)[indx]
       write.table(x = colnames(data)[indx], file = paste0(write_path,rownames(data)[i]), sep = "\n", row.names = FALSE, col.names = FALSE)
    }
    return (result)
}

### add #2017.09.12
### convert KEGG access number to ncbi protein id
convert2ncbi_protein_id <- function(file_in_path, file_out_path){
    t2d_files <- list.files(path = file_in_path, all.files = FALSE, full.names = FALSE, recursive = FALSE, include.dirs = FALSE)
    lapply(t2d_files, function(x) {
        file_in_name <- paste0(file_in_path, x)
        print(file_in_name)
        kegg_genes_id <- read.table(file = file_in_name, sep = "\n", header = FALSE)
        set <- c()
        for(gene in kegg_genes_id){
            ncbi_protein_id <- unname(keggConv("ncbi-proteinid",gene))
            set <- c(set, ncbi_protein_id)
        }
        file_out_name <- paste0(file_out_path, x)
        write.table(x = set, file = file_out_name, sep = "\n", row.names = FALSE, col.names = FALSE)
    })
}
