
#' Read counts from STARSolo to a column sparse matrix.
#'
#' @param p is the path to the STAR solo output that contains the matrix.mtx,
#'   features.tsv and barcodes.tsv files
#' @param gene_col is the column index to use for the gene names in the 
#'   barcodes.tsv file
#' @param sep is the separator if gene names are not unique, 
#'   the argument to make.unique function
read_solo_counts <- function(p, gene_col = 1, sep = "."){
    require(Matrix)
    mtx_file <- file.path(p, "matrix.mtx")
    genes_file <- file.path(p, "features.tsv")
    barcode_file <- file.path(p, "barcodes.tsv")

    counts <- Matrix::readMM(mtx_file)
    barcode_names <- readLines(barcode_file)
    genes <- read.delim(genes_file, 
                        header = FALSE, 
                        sep = "\t", 
                        stringsAsFactors = FALSE)


    colnames(counts) <- make.unique(barcode_names, sep=sep)
    rownames(counts) <- make.unique(genes[,gene_col], sep=sep)

    counts <- as(counts, "CsparseMatrix")

    return(counts)
}

set_breaks_10 <- function(x){
	xmax <- x[2]
	bk <- 10
	brks <- c(bk)
	while (bk < xmax){
		bk <- bk * 10
		brks <- c(brks, bk)
	}
	return(brks)
}

#' Create a barcode rank plot
#'
#' @param counts A sparse count matrix with genes in rows and droplets in 
#'  columns. Must contain all droplets.
#' @param title Title of the plot
#' @param ret Logical indicating whether to return a ggplot object
#'
#' @return Nothing. If return=TRUE, then return a ggplot object
#' @import ggplot2
#' @importFrom scales comma
#' @export
barcode_rank_plot <- function(sce, title = "Barcode-rank", ret = TRUE){
    require(Matrix)
    require(ggplot2)
    require(scales)
    tc <- Matrix::colSums(counts)
    tc <- sort(tc, decreasing=TRUE)
    ranks <- seq(length(tc))
    dups <- duplicated(tc)
    tc[dups] <- NA
    ranks[dups] <- NA
    datf <- data.frame("Rank"=ranks, "Count"=tc)
    datf <- datf[!is.na(datf[,2]),,drop=FALSE]
    datf <- datf[datf[,"Count"] != 0,]
    p <- ggplot(datf, aes_string(x = "Rank", y = "Count")) +
        geom_point() +
        scale_x_continuous(trans='log10', breaks=set_breaks_10, labels=comma) +
        scale_y_continuous(name="Droplet size", trans='log10', labels=comma)  +
        theme_minimal() + 
        theme(plot.title=element_text(hjust=0.5),
              axis.text.x=element_text(angle=45, hjust=1)) +
        ggtitle(title)
    if (ret) return(p)
    else print(p)
}

#' Get proportions from the columns of a data frame
#'
#' @param x A data frame
#' @param c1 Column 1, which will be the rows of the proportions
#' @param c2 Column 2, which will be the columns of the proportions
#' @param prop Whether to return proportions (TRUE) or counts (FALSE). Default is TRUE.
#'
#' @return data frame
get_props <- function(x, c1="orig.ident", c2="seurat_clusters", prop=TRUE){
    datf <- table(x[,c2], x[,c1])
    if (prop){
        datf <- sweep(datf, 2, colSums(datf), "/")
    }
    datf <- as.data.frame.array(datf)
    return(datf)
}

#' Collapse expression along cell types
#'
#' @param x Sparse matrix of counts, the raw gene-barcode matrix
#' @param identities Character or factor vector of cell type identities. 
#'  Each element corresponds to the column in x.
#' @param counts Whether to return counts or normalized counts.
#' @param scale_size If normalizing, scale so column sums equal this size.
#' @param scale_cells Whether to scale so column sums are all equal.
#' @param logt Whether to log1p transform the counts.
#'
#' @return matrix, with rows corresponding to genes and columns 
#'  corresponding to identities.
collapse = function(x, 
                    identities, 
                    counts=FALSE, 
                    scale_size=1e4, 
                    scale_cells=TRUE,
                    logt=TRUE){
    require(Matrix)
    if (is.null(attr(class(x), "package")) || 
        attr(class(x), "package") != "Matrix"){
        x <- Matrix::Matrix(x)
    }
    ngenes <- nrow(x)
    if (class(identities) != "factor"){
        identities <- factor(identities)
    }
    celltypes <- levels(identities)
    ncells <- table(identities)[celltypes]
    panel <- matrix(0, nrow=ngenes, ncol=length(celltypes))
    rownames(panel) <- rownames(x)
    colnames(panel) <- celltypes
    for (celltype in celltypes){
        panel[,celltype] <- Matrix::rowSums(x[,identities == celltype,drop=FALSE])
    }

    if (counts) return(panel)

    cs <- colSums(panel)
    cs[cs == 0] <- 1
    panel <- sweep(panel, 2, cs, "/")
    
    if (scale_cells){
        panel <- scale_size*panel
    }
    
    if (logt){
        panel <- log1p(panel)
    }

    return(panel)
}

#' Collapse counts across 2 factors
#' 
#' 
collapse2d <- function(x, 
                       ix1, 
                       ix2){
    require(Matrix)

    if (is.null(attr(class(x), "package")) || 
        attr(class(x), "package") != "Matrix"){
        x <- Matrix::Matrix(x)
    }

    if (length(ix1) != ncol(x) || length(ix2) != ncol(x))
        stop("i1 and i2 must have length equal to number of columns in x")

    i1 <- rownames(x)
    i1n <- rownames(x)
    i2 <- factor(ix1)
    i2n <- levels(i2)
    i3 <- factor(ix2)
    i3n <- levels(i3)

    panel <- array(data = NA, 
                   dim = c(length(i1n), length(i2n), length(i3n)), 
                   dimnames = list(i1n, i2n, i3n))

    for (k in i3n){
        for (j in i2n){
            colkeep <- i2 == j & i3 == k
            panel[,j,k] <- Matrix::rowSums(x[,colkeep,drop=FALSE])
        }
    }

    return(panel)
}






