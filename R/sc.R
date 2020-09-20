
#' Read counts from STARSolo to a column sparse matrix.
#'
#' @param p is the path to the STAR solo output that contains the matrix.mtx,
#'   features.tsv and barcodes.tsv files
#' @param gene_col is the column index to use for the gene names in the 
#'   barcodes.tsv file
#' @param sep is the separator if gene names are not unique, 
#'   the argument to make.unique function
read_solo_counts <- function(p, gene_col = 1, sep = "."){
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

