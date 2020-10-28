#### Function to compute the Goulard & Voltz (1992) algorithm
#### Goulard, M., & Voltz, M. (1992).
#### Linear coregionalization model: tools for estimation and choice of cross-variogram matrix.
#### Mathematical Geology, 24(3), 269-286.

#### Author of the code: Jonas Alcaina-Mateos (jonas.alcaina@upf.edu)

# build a square matrix from a vector 
buildMat <- function(x) {
    n <- (sqrt(1 + (8*length(x)))-1) / 2
    m <- matrix(NA, n, n)
    m[lower.tri(m, diag=TRUE)] <- x
    m <- t(m)
    m[lower.tri(m, diag=TRUE)] <- x
    return(m) 
}

# positive semi-definite function
posDef <- function(x) {
    q <- eigen(buildMat(x))
    d <- q$values
    d[d < 0] <- 1e-10   # small positive number 
    res <- q$vectors %*% diag(d) %*% t(q$vectors)
    return(res[lower.tri(res, diag=T)])
}

# Goulard & Voltz (1992) algorithm
fitGV <- function(semvar, model) {
    # prepare data 
    semvar <- split(semvar, f=semvar$id)
    model <- model[1:length(semvar)] # eliminate non-models
    semvar <- semvar[names(model)]
    # matrices and vectors 
    sills <- outer(1:length(model), 1:nrow(model[[1]]), FUN=Vectorize(function(i,j) model[[i]]$psill[j]))
    yvect <- sapply(1:length(semvar), function(i) matrix(semvar[[i]]$gamma, ncol=1)) 
    values <- sapply(1:nrow(model[[1]]), function(x) variogramLine(vgm(1, as.character(model[[1]][x,]$model), model[[1]][x,]$range), dist_vector = semvar[[1]]$dist)$gamma)
    w <- sapply(1:length(semvar), function(x) semvar[[x]]$np)
    # iterate (max. iteration: 500)
    WSS <- c(sum(w * (yvect - (values %*% t(sills)))^2))
    for (it in 1:500) {
        for (s in 1:ncol(sills)) {
            aux1 <- (yvect - (values[,-s] %*% t(sills[,-s])))
            aux2 <- w * apply(aux1, 2, function(x) values[,s]*x)
            aux3 <- colSums(aux2) / colSums(w * values[,s]^2)
            sills[,s] <- posDef(aux3)
        }
        WSS[it+1] <- sum(w * (yvect - (values %*% t(sills)))^2)
        # stop criterion (10e-4)
        if(it > 1 & abs(WSS[it] - WSS[it+1]) < 10^-4) break 
    }
    # arrange the model
    for (i in 1:length(model)) model[[i]]$psill <- sills[i,] 
    return(model)
}
