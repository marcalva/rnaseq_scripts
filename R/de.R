

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


