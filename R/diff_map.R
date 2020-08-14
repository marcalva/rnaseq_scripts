
#' Guassian kernel formula
#'
#' @param x Numeric vector for which to calculate guassian kernel.
#' @param sig Sigma squared value.
#'
#' @return Numeric vector giving guassian kernels
gk <- function(x, sig = 1){
    exp(-0.5 * (x^2/sig))
}   

#' Calculate gaussian kernel for a distance matrix
#'
#' Calculates the guassian kernel for each row in \code{d} using 
#' the width (variance) given by \code{sig}, or if \code{kn} is 
#' specified, using an adaptive kernel. The width (variance) of 
#' the adaptive guassian kernel is given by the squared distance 
#' to the \code{kn}th closest neighbor.
#'
#' @param d Numeric vector of distance values
#' @param sig Sigma squared value
#' @param kn Integer giving the kth neighbor from which to specify the 
#'  width of the guassian kernel.
#' 
#' @return Numeric value giving gaussian kernel.
gkd <- function(d, sig = 1, kn = NULL){
    if (kn > ncol(d)){
        stop("kn must be less than or equal to the number of columns in d")
    }

    if (is.null(kn)){
        r <- apply(d, c(1,2), gk, sig)
    } else {
        r <- apply(d, 1, function(x){
                   sig <- x[kn]^2
                   gk(x, sig) })
        r <- t(r)
    }
    return(r)
}

# spase matrix
#   i is row index of corresponding x value
#   p is starting index in x
#   x are values

#' Convert nearest neighbor object to sparse matrix
#'
#' @param nn Nearest neighbor object returned by RANN's nn2 function.
#'
#' @return A sparse matrix giving the distance values between pairs.
nn2sparse <- function(nn){
    idx <- nn$nn.idx[,1:ncol(nn$nn.idx)]
    dists <- nn$nn.dists[,1:ncol(nn$nn.dists)]
    i <- sapply(1:nrow(idx), rep, ncol(idx))
    i <- as.vector(i) # row index
    j <- as.numeric(t(idx)) # column index
    x <- as.numeric(t(dists))
    d <- c(nrow(idx), nrow(idx))
    mat <- sparseMatrix(i = i, j = j, x = x, dims = d)
    return(mat)
}

#' Calculates an adjacency matrix.
#'
#' First calculates the \code{k} nearest neighbors and the euclidean distances 
#' for each data point using the 
#' \code{\link[RANN]{nn2}} function from the RANN package. Then calculates 
#' the guassian kernel between these pairs to get teh adjacency matrix. 
#' By default, it calculates an 
#' adaptive gaussian kernel by setting the variance to the squared distance 
#' to the \code{kn} nearest neighbor for each data point. The parameter 
#' \code{kn} can be specified, or set to 1/3 of \code{k} (default). 
#' Alternatively, the width (variance) can be set with \code{sig}.
#'
#' @param x An n by m matrix, corresponding to n data points and 
#'  m dimensions (e.g. genes or PCs).
#' @param k Number of nearests neighbors to calculate distances up 
#'  to.
#' @param Asym Whether to symmetrize the nearest neighbos graph 
#'  by adding the transpose
#' @param sig Sigma squared value of gaussian kernel.
#' @param kn Set the standard deviation of the gaussian kernel 
#'  for each data point to the distance to the kn nearest neighbor.
#'  Sets kn to 1/3 of k by default.
#' @param weighted Boolean, indicating whether to weigh each row 
#'  by the sum of adjacencies.
#' 
#' @return Sparse adjacency matrix.
get_adj <- function(x, 
                    k = 200, 
                    Asym = TRUE, 
                    sig = NULL, 
                    kn = NULL, 
                    weighted = TRUE){
    require(RANN)
    require(RSpectra)
    require(Matrix)

    message("Getting nearest neigbors...")
    x <- as.matrix(x)
    nns <- RANN::nn2(data = x, k = k)

    # Convert nn distances to adjacency with guassian kernel
    if (!is.null(sig)){
        nns$nn.dists <- gkd(nns$nn.dists, sig = sig, kn = NULL)
    } else {
        if (is.null(kn)) kn <- floor(k/3)
        nns$nn.dists <- gkd(nns$nn.dists, sig = sig, kn = kn)
    }

    # Convert nn object to sparse matrix
    A <- nn2sparse(nn = nns)

    # Symmetrize
    if (Asym)
        A <- A + t(A)

    # Weighted
    if (weighted){
        w <- 1 / Matrix::rowSums(A)
        W <- Matrix::Diagonal(x = w, n = ncol(A))
        A <- W %*% A
    }

    return(A)
}

#' Calculate diffusion map
#'
#' @param A An adjacency matrix.
#' @param a Multiply the eigenvalues by this constant. This value 
#'  can be increased to model a diffusion process with a steps.
#' @param p Number of eigenvectors to return.
#'
#' @return The object returned by the eigs function from RSpectra, 
#'  removing the first trivial constant eigenvector/eigenvalue.
diff_map <- function(A, a = 1, p = 10){
    require(RANN)
    require(RSpectra)
    require(Matrix)

    message("Running eigendecomposition...")
    ed <- RSpectra::eigs(A = A, k = p+1, type = "LR")
    ed$values <- Re(ed$values)[2:length(ed$values)]
    ed$vectors <- Re(ed$vectors)[,2:ncol(ed$vectors)]
    ed$values <- ed$values ^ a
    message("Done")
    return(ed)
}

