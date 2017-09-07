### fit the linear mixed model
### ptr_data : t2d or control group
### kegg_data : t2d or control group
### organism_species : organism:species pairs
### result_maxtrix : 0-1 matrix for ptr and kegg 
fit <- function(ptr_data, kegg_data, organism_species, result_maxtrix){
    K <- cor(ptr_data) #####right???????
    decomp <- eigen(K, symmetric=TRUE)
    tmp  <- t(ptr_data)
    
    ### the final result for lmm
    lmm_result_maxtrix <- data.frame(matrix(nrow = nrow(organism_species), ncol = nrow(kegg_data), byrow = TRUE)) 
    rownames(lmm_result_maxtrix) <- organism_species$organism
    colnames(lmm_result_maxtrix) <- rownames(kegg_data)
    
    result_maxtrix_row_names <- rownames(result_maxtrix)
    for(i in 1:nrow(kegg_data)){
        y <- t(kegg_data[i,])
        xx <- c()
        for(j in 1:40){
            if(!is.na(result_maxtrix[j,i]) & result_maxtrix[j,i] == 1){
                xx <- c(xx,j)
            }
        }
        
        if(length(xx) >0){
            formula <- as.formula(paste("y~", paste(paste0(paste0("tmp[,",xx),"]"), collapse="+")))
            #print(formula)
            fit <- lrgpr(formula, decomp)
            #fit <- lm(formula,data = cbind(y,tmp)) # fit the linear model
            p_values <- as.vector(fit$p.values)
            #print(p_values)
            if(length(p_values) > 0){
                p_values <- p_values[2:length(p_values)] ### remove the intercepet column
                lmm_result_maxtrix[xx,i] <- p_values
            }
        }
    }
    return (lmm_result_maxtrix)
}

fit2 <- function(ptr_data, kegg_data, organism_species, result_maxtrix){
    decomp <- eigen(cor(ptr_data), symmetric=TRUE)
    species  <- t(ptr_data)
    
    ### the final result for lmm
    lmm_result_maxtrix <- data.frame(matrix(nrow = nrow(organism_species), ncol = nrow(kegg_data), byrow = TRUE)) 
    rownames(lmm_result_maxtrix) <- organism_species$organism
    colnames(lmm_result_maxtrix) <- rownames(kegg_data)
    
    result_maxtrix_row_names <- rownames(result_maxtrix)
    for(i in 1:nrow(kegg_data)){
        y <- t(kegg_data[i,])
        xx <- c()
        for(j in 1:40){
            if(!is.na(result_maxtrix[j,i]) & result_maxtrix[j,i] == 1){
                xx <- c(xx,j)
            }
        }
        
        if(length(xx) >0){
            formula <- as.formula(paste("y~", paste(paste0(paste0("species[,",xx),"]"), collapse="+")))
            #print(formula)
            fit <- lrgpr(formula, decomp)
            #fit <- lm(formula,data = cbind(y,tmp)) # fit the linear model
            p_values <- as.vector(fit$p.values)
            #print(p_values)
            if(length(p_values) > 0){
                p_values <- p_values[2:length(p_values)] ### remove the intercepet column
                lmm_result_maxtrix[xx,i] <- p_values
            }
        }
    }
    return (lmm_result_maxtrix)
}

### fit the linear mixed model
### ptr_data : t2d or control group: species * samples
### kegg_data : t2d or control group: ko * samples
### result_maxtrix : 0-1 matrix for ptr and kegg: species * ko
lmm_fit <- function(ptr_data, kegg_data, result_maxtrix){
    result_maxtrix[is.na(result_maxtrix)] <- 0 # replace the NA with 0
    decomp <- svd(scale(t(ptr_data)))
    #K <- cor(ptr_data) # covariance matrix
    #decomp <- eigen(K, symmetric=TRUE)
    # y is n * 1 vector, n is number of samples
    # x is n * 1 vector, n is number of samples
    # α is 1 * m vector, m is number of species
    # β is m * 1 vector, m is number of species
    # formula: y = x*α*β
    # first, get the design matrix x*α
    for(i in 1:nrow(kegg_data)){
        y <- t(kegg_data[i,])
        ko <- rownames(kegg_data)[i]
        x <- t(ptr_data)
        #α <- result_maxtrix[,c(ko)]
        α <- which(result_maxtrix[,c(ko)] == 1)
        if(length(α) > 0){
            formula <- as.formula(paste("y~", paste(paste0(paste0("x[,",α),"]"), collapse="+")))
            fit <- lrgpr(formula, decomp)
            p_values <- as.vector(fit$p.values)
            if(length(p_values) > 0){
                p_values <- p_values[2:length(p_values)] ### remove the intercepet column
                result_maxtrix[α,c(ko)] <- p_values
            }
        }
    }
    return (result_maxtrix)
}