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
convet2genes <- function(path, species_ko, P_VALUE){
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
                    ll <- unlist(strsplit(genes[j],":"))
                    if(length(ll) >1){
                        genes_set <- c(genes_set, ll[2])
                    }
                }   
            }
        }
        write.table(x = genes_set, file = paste0(path,row_name), sep = "\n", row.names = FALSE, col.names = FALSE)
        genes_set <- c()
    }
}

map2pathway <- function(data){
    # KO map to pathway
    result <- c()
    k <- c()
    for(ko in colnames(data)){
        k <- c(k,ko)
        if(length(k)>= 10){
            result <- c(result,names(keggGet(k)[[1]]$PATHWAY))
            k <- c()
        }
    }
    if(length(k) >0){
        result <- c(result,names(keggGet(k)[[1]]$PATHWAY)) 
    }
    result <- as.factor(result)
    
    return (result)
}

get_significant_ko <- function(write_path, data, p_value){
    for(i in 1:ncol(data)){
       indx <- which(data[i,]<p_value & data[i,]>0)
       write.table(x = colnames(data)[indx], file = paste0(write_path,rownames(data)[i]), sep = "\n", row.names = FALSE, col.names = FALSE)
    }
}
