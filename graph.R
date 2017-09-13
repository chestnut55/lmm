library(multiplex)
draw_bip_network <- function(data, p_value){
    for(i in 1: ncol(data)){
        index <- which(data[,i]>0&data[,i]<p_value)
        data[index, i] <- 1
        data[-index,i] <- 0
    }
    bmgraph(data, layout = "force", seed = 123, cex = 1.5,tcex = .5, pch = c(19, 19), lwd = 2,  
            vcol = 2:3, ecol = 8, rot = 65)
}
