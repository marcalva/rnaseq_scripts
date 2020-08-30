
#' Save a data frame to file
save2file <- function(datf, 
                      dirout = "./", 
                      fn = "datf.txt", 
                      gz = FALSE, 
                      row.names = TRUE, 
                      col.names = NA, 
                      quote = FALSE, 
                      sep = "\t"){
    dir.create(dirout, showWarnings = FALSE, recursive = TRUE)
    fn <- file.path(dirout, fn)
    if (gz){
        fn <- paste0(fn, ".gz")
        gz1 <- gzfile(fn, "w")
        write.table(datf, 
                    gz1, 
                    row.names = row.names, 
                    col.names = col.names, 
                    quote = quote, 
                    sep = sep)
        close(gz1)
    } else {
        write.table(datf, 
                    fn, 
                    row.names = row.names, 
                    col.names = col.names, 
                    quote = quote, 
                    sep = sep)
    }
}

#=================================================
#=================================================

runedger <- function(yf, ds, prior.count = 0.125, p.adjust = "fdr"){
    require(edgeR)
    fit <- glmFit(yf, ds, prior.count = prior.count, coef = ncol(yf$design))
    lrt <- glmLRT(fit, coef = ncol(fit$design))
    diffTable <- lrt$table

    diffTable$p_adj <- p.adjust(diffTable$PValue, method = p.adjust)
    diffTable$FoldChange <- 2^(diffTable$logFC)
    diffTable$AverageCPM <- 2^(diffTable$logCPM)
    diffTable <- diffTable[,c("logFC", "FoldChange", "logCPM", "AverageCPM", "LR", "PValue", "p_adj")]
    diffTable <- diffTable[order(diffTable$PValue),]
    return(list("diffTable" = diffTable, "fit" = lrt))
}

runlimma <- function(yf, ds, prior.count = 0.125, lfc = 0, 
                     p.adjust = "fdr"){
    require(edgeR)
    logCPM <- cpm(yf, log = TRUE, prior.count = prior.count)

    fit <- lmFit(logCPM, ds)
    fit <- eBayes(fit, trend=TRUE)
    diffTable <- topTable(fit, lfc = lfc, coef=ncol(ds), number = nrow(fit), 
                          adjust.method = p.adjust)
    
    diffTable$FoldChange = 2^(diffTable$logFC)
    diffTable = diffTable[,c("logFC", "FoldChange", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
    colnames(diffTable)[ colnames(diffTable) == "P.Value" ] <- "PValue"
    colnames(diffTable)[ colnames(diffTable) == "p_adj" ] <- "p_adj"
    return(list("diffTable" = diffTable, "fit" = fit))
}

#=================================================
#=================================================

run_de_all <- function(counts, pheno, covars = NULL, gene_info = NULL,
                       de_method = "edgeR", lfc = 0,  p.adjust = "fdr", 
                       plots = FALSE, 
                       nm_method = "TMM", prior.count = 0.125, min_cpm = 0, 
                       sfrac = 0.9, dirout = "", verbose = TRUE, gz = TRUE){
    require(edgeR)

    s <- intersect(colnames(counts), rownames(pheno))
    if (!is.null(covars)){
        s <- intersect(s, rownames(covars))
        covars <- covars[s,,drop=FALSE]
    }
    pheno <- pheno[s,,drop=FALSE]
    counts <- counts[,s,drop=FALSE]

    if (!is.null(covars)){
        dat <- as.data.frame(covars)
        for (i in colnames(pheno)){
            dat[,i] <- pheno[,i]
        }
    } else {
        dat <- pheno
    }


    y <- DGEList(counts = counts, genes = gene_info)
    keep <- rowSums(cpm(y) > min_cpm) >= (sfrac * ncol(y))
    y <- y[keep, , keep.lib.sizes=FALSE]

    y <- calcNormFactors(y, method = nm_method)

    if (plots){
        fn <- paste0(dirout, "mds.edgeR.pdf")
        pdf(fn)
        plotMDS(y)
        dev.off()
    }

    for (p in colnames(pheno)){
        if (verbose){
            message("Testing ", p)
        }
        if (is.null(covars)){
            ds <- pheno[,p,drop=FALSE]
        } else {
            ds <- covars
            ds[,p] <- pheno[,p]
        }
        ds <- model.frame(ds)
        design <- model.matrix(terms(ds), ds)

        y <- DGEList(counts = counts[,rownames(design)], genes = gene_info)
        keep <- rowSums(cpm(y) > min_cpm) >= (sfrac * ncol(y))
        y <- y[keep, , keep.lib.sizes=FALSE]

        if (verbose){
            message("testing ", nrow(y), " genes across ", ncol(y), " samples")
        }

        y <- calcNormFactors(y, method = nm_method)
        y <- estimateDisp(y, as.data.frame(ds), robust=TRUE)

        if (de_method == "edgeR"){
            ret <- runedger(y, design, prior.count = prior.count, p.adjust = p.adjust)
            de <- ret$diffTable
            if (plots){
                pdf(paste0(dirout, p, ".bcv.pdf"))
                edgeR::plotBCV(y)
                dev.off()
            }
            deg <- cbind(gene_info[rownames(de),], de)

            fn <- paste0(dirout, p, ".edgeR.txt")
            if (gz){
                fn <- paste0(fn, ".gz")
                gz1 <- gzfile(fn, "w")
                write.table(deg, gz1, row.names = T, col.names = NA, quote = F, sep = "\t")
                close(gz1)
            } else {
                write.table(deg, fn, row.names = T, col.names = NA, quote = F, sep = "\t")
            }
        } else if (de_method == "limma"){
            ret <- runlimma(y, design, prior.count = prior.count, lfc = lfc, 
                            p.adjust = p.adjust)
            fit <- ret$fit
            de <- ret$diffTable
            if (plots){
                pdf(paste0(dirout, p, ".SA.pdf"))
                limma::plotSA(fit)
                dev.off()
            }
            deg <- cbind(gene_info[rownames(de),], de)

            fn <- paste0(dirout, p, ".limma.txt")
            if (gz){
                fn <- paste0(fn, ".gz")
                gz1 <- gzfile(fn, "w")
                write.table(deg, gz1, row.names = T, col.names = NA, quote = F, sep = "\t")
                close(gz1)
            } else {
                write.table(deg, fn, row.names = T, col.names = NA, quote = F, sep = "\t")
            }
        }
        else{
            stop(sQuote(de_method), " must be one of ", sQuote(edgeR)," or ", sQuote(limma))
        }
    }
}


#=================================================
#=================================================

#' Run DE using a linear model
#'
#' Given a gene by sample matrix \code{counts} and a design matrix
#' \code{design}, run a linear model for each gene. Returns the 
#' coefficients and p-value for the terms specified by \code{cf} 
#' in the design matrix. The 
#' model is run with \code{lm(x ~ 0 + design)}.
#'
#' @param counts Normalized expression data. Genes are in rows and 
#'  samples are in columns.
#' @param design Design matrix. The last column is the variable of 
#'  interest. Note this design matrix must include the intercept term, 
#'  as the lm function is called without an intercept.
#' @param cf Coefficients to return. By default, returns the last term 
#'  in the design matrix.
#' @param p_adj Adjust p values using this method from \link{p.adjust}.
#' @param p_filter Remove genes with an adjusted p-value greater than this 
#'  value.
#' @param gene_info Data frame cotnaining gene annotations to append to 
#'  results.
#' @param save2file Boolean indicating whether to save results to file.
#' @param dirout Output directory.
#' @param fn File name.
#' @param gz Boolean indicating whether to gzip the output file.
#' 
#' @return A data frame with one gene per row.
run_de_lm <- function(counts, 
                      design, 
                      cf = NULL, 
                      p_adj = "fdr", 
                      gene_info = NULL,
                      save2file = TRUE, 
                      dirout = "./", 
                      fn = "de.lm.txt", 
                      gz = TRUE){

    if (nrow(design) != ncol(counts)){
        stop("Number of rows in design must be the same as number of columns in counts")
    }

    if (is.null(cf)) cf <- ncol(design)
    c_names <- sapply(cf, function(i){
                      c(paste0(colnames(design)[i], "_effect"), 
                        paste0(colnames(design)[i], "_p")) })
    detable <- apply(counts, 1, function(x){
                     lmr <- lm(x ~ 0 + design)
                     lmrs <- summary(lmr)
                     coefs <- lmrs$coefficients
                     ret <- as.vector(t(coefs[cf, c(1,4)]))
                     names(ret) <- c_names
                     return(ret) })
    detable <- t(detable)
    detable <- as.data.frame(detable)

    # Adjust p-values
    for (i in cf){
        cname <- paste0(colnames(design)[i], "_p")
        cname2 <- paste0(colnames(design)[i], "_p.adj")
        detable[,cname2] <- p.adjust(detable[,cname], method = p_adj)
    }

    # Re-order columns
    reo <- c()
    for (i in cf){
        addto <- c(paste0(colnames(design)[i], "_effect"), 
                   paste0(colnames(design)[i], "_p"),
                   paste0(colnames(design)[i], "_p.adj"))
        reo <- c(reo, addto)
    }
    detable <- detable[,reo]

    # Add anotations
    if (!is.null(gene_info))
        detable <- cbind(gene_info[rownames(detable),], detable)

    return(detable)
}



