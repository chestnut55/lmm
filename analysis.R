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

### map KO of species to pathway:species -> KO -> pathway
map2pathway <- function(filename, data, P_VALUE){
    result <- c()
    rownames <- rownames(data)
    colnames <- colnames(data)
    for(i in 1:nrow(data)){
        value <- c()
        indx <- which(!is.na(data[i,]) & data[i,] < P_VALUE & data[i,] >0)
        for(k in colnames[indx]){
            pathway <- names(keggGet(k)[[1]]$PATHWAY)
            if(length(pathway) >0){
                for(j in 1:length(pathway)){
                    value <- rbind(value, c(rownames[i], k, pathway[j]))
                }
            }
        }
        if(length(value) >0){
            result <- rbind(result, value)
        }
        #write.table(x = value, file = paste0(path,rownames[i]), sep = "\t", row.names = FALSE, col.names = FALSE)
    }
    
    write.table(x = result, file = filename, sep = "\t", row.names = FALSE, col.names = FALSE)
    return (result)
}

### count the species in the pathway
### list : species as key, pathway list as value
count_species_in_pathway <- function(data){
    species <- as.character(unique(data[,1]))
    pathways <- as.character(unique(data[,3]))
    matrix <- data.frame(matrix(nrow = length(species), ncol = length(pathways)))
    colnames(matrix) <- pathways
    rownames(matrix) <- species
    
    for(pathway in pathways){
        for(sp in species){
            matrix[c(sp),c(pathway)] <- length(which(data[,3] == pathway & data[,1] == sp))
        }
    }
    
    return (matrix)
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

### matrix : species * pathway
### the value in each matrix cell is the number of species in the pathway
sampling <-function(filename, matrix){
    r_names <- rownames(matrix)
    c_names <- colnames(matrix)
    r_sums <- rowSums(matrix)
    for(i in 1:length(r_names)){
        ### random matrix
        m <- c()
        for( j in 1:10000){
            s <- sample(c_names, r_sums[[i]], replace = TRUE)
            m <- rbind(m,s)
        }
        ### count
        indx <- which(matrix[i,] >0)
        for(l in 1:length(indx)){
            v <- matrix[i, indx[l]] # the number of species in the pathway
            colname <- c_names[l]
            count <- 0
            for(k in 1:nrow(m)){
                len <- length(which(str_count(m[k,], colname) == 1))
                if(len >= v){
                    count <- count + 1
                }
            }
            matrix[i, indx[l]] <- count
        }
    }
    write.table(x = matrix, file = filename, sep = "\t")
    return (matrix)
}
