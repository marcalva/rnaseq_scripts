
#' Analyze powers that maximize scale free topology
#' 
#' @param x Normalized expression data. Rows correspond to samples and 
#'  columns correspond to genes.
#' @param powers Powers to test for scale-free topology.
#' @param plot_top Boolean, whether to output a plot
#' @param plot_dir Plot directory
#' @param plot_file Plot file name
#' @param nThreads Number of threads.
#' @param ... Additional parameters to pass to pickSoftThreshold.
#'  For example, signed networks can be constructed by adding 
#'  networkType = "signed".
#'
#' @return The list returned by pickSoftThreshold 
choose_power <- function(x, 
                         powers = c(c(1:10), seq(from = 12, to = 20, by = 2)), 
                         type = "unsigned", 
                         plot_top = TRUE, 
                         plot_dir = "./", 
                         plot_file = "thresh.sft.pdf", 
                         nThreads = NULL, 
                         ...){
    require(WGCNA)
    n_genes <- ncol(x)
    n_samples <- nrow(x)
    message("Analyzing topology in ", n_genes, " genes and ", 
            n_samples, " samples")

    enableWGCNAThreads(nThreads = nThreads)

    sft <- pickSoftThreshold(x, 
                             powerVector = powers, 
                             networkType = type, 
                             verbose = 5, ...)

    if (!plot_top)
        return(sft)

    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    outfn <- paste0(plot_dir, "/", plot_file)
    
    pdf(outfn, width = 7, height = 4)
    par(mfrow = c(1,2))
    cex1 = 0.8

    # Left plot
    xc <- sft$fitIndices[,1]
    yc <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
    plot(xc, yc, 
         xlab="Soft Threshold (power)",
         ylab="Scale Free Topology Model Fit,signed R^2",
         type="n",
         main = paste("Scale independence"))
    text(xc, yc, 
         labels=powers, 
         cex=cex1, 
         col="red");
    abline(h=0.90,col="red")

    # Right plot
    xc <- sft$fitIndices[,1]
    yc <- sft$fitIndices[,5]
    plot(xc, yc, 
         xlab="Soft Threshold (power)", 
         ylab="Mean Connectivity", 
         type="n", 
         main = paste("Mean connectivity"))
    text(xc, yc, 
         labels=powers, 
         cex=cex1, 
         col="red")

    dev.off()

    return(sft)
}

#' Calculate adjacencies and topological overlap
#'
#' The adjacency function returns a large gene by gene matrix. The 
#' matrix is the correlation across genes, with the gene-gene coefficients 
#' raised to a power of soft_power. TOMsimilarity calculates the 
#' topological overlap matrix. The topological overlap measures  
#' similarity by how many neighbors are shared bewteen two genes. The 
#' neighbors are specified by the adjacency matrix.
#' 
#' @param x Normalized expression data. Rows correspond to samples and 
#'  columns correspond to genes.
#' @param soft_power Powers to test for scale-free topology.
#' @param type Network type to pass to adjacency and TOMsimilarity.
#'
#' @return a TOM dissimilarity matrix
get_adj_tom <- function(x, 
                        soft_power = 6, 
                        type = "unsigned", 
                        nThreads = NULL){
    require(WGCNA)
    n_genes <- ncol(x)
    n_samples <- nrow(x)
    message("Calculating adjacencies and TOM in ", n_genes, " genes and ", 
            n_samples, " samples")

    enableWGCNAThreads(nThreads = nThreads)

    adj <- adjacency(x, power = soft_power, type = type)
    TOM <- TOMsimilarity(adj, TOMType = type)
    dissTOM <- 1 - TOM

    return(list("adj" = adj, "dissTOM" = dissTOM))
}

#' Clustering
#'
#'
cluster_tom <- function(x, 
                        dissTOM, 
                        method = "average", 
                        minModuleSize = 30, 
                        deepSplit = 2, 
                        pamRespectsDendro = FALSE, 
                        MEDissThres = 0.25,
                        plots = TRUE, 
                        plot_dir = "./", 
                        dendro_plot = "dendro.pdf", 
                        w = 7, 
                        h = 5,
                        me_plot = "me_clust.pdf", 
                        nThreads = NULL){

    require(WGCNA)
    n_genes <- ncol(dissTOM)
    n_samples <- nrow(x)
    message("Clustering TOM using ", n_genes, " genes")

    enableWGCNAThreads(nThreads = nThreads)

    geneTree <- hclust(as.dist(dissTOM), method = method)

    # returns a vector of labels assigning each gene to a cluster
    dynamicMods <- cutreeDynamic(dendro = geneTree, 
                                 distM = dissTOM, 
                                 deepSplit = deepSplit, 
                                 pamRespectsDendro = pamRespectsDendro, 
                                 minClusterSize = minModuleSize)
    dynamicColors <- labels2colors(dynamicMods)

    MEList <- moduleEigengenes(x, colors = dynamicColors)
    MEs <- MEList$eigengenes
    MEDiss <- 1 - cor(MEs)
    METree <- hclust(as.dist(MEDiss), method = "average")
    merged <- mergeCloseModules(x, 
                                dynamicColors, 
                                cutHeight = MEDissThres, 
                                verbose = 3)

    if (plots){
        dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
        pfile <- paste0(plot_dir, me_plot)
        pdf(pfile, width = w, height = h)
        plot(METree, main = "Clustering of module eigengenes", 
             xlab = "", sub = "")
        abline(h = MEDissThres, col = "red")
        dev.off()
    }

    if (plots){
        dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
        pfile <- paste0(plot_dir, dendro_plot)
        pdf(pfile, width = w, height = h)
        plotDendroAndColors(geneTree, merged$colors, "Dynamic Tree Cut",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05,
                            main = "Gene dendrogram and module colors")
        dev.off()
    }

    names(merged$colors) <- colnames(x)

    return(merged)
}

#' Plot variance explained by each module eigengene
plot_me_ve <- function(x, gene_lab, order_ve = TRUE, ret = TRUE){
    require(ggplot2)
    mods <- unique(gene_lab)
    ve <- sapply(mods, function(m){
                   g <- names(gene_lab)[gene_lab == m]
                   pcs <- prcomp(x[,g,drop=FALSE], scale.=TRUE, center = TRUE)
                   vars <- pcs$sdev^2
                   vars <- vars / sum(vars)
                   return(vars[1])
                            })
    datf <- data.frame("ME" = mods, 
                       "VarExpl" = ve)
    if (order_ve){
        o <- order(datf[,"VarExpl"], decreasing = TRUE)
        datf <- datf[o,,drop=FALSE]
    }
    datf[,"ME"] <- factor(datf[,"ME"], levels = datf[,"ME"])

    p <- ggplot(datf, aes(x = ME, y = VarExpl)) + 
        geom_point() + 
        theme_classic() + 
        theme(axis.text.x = element_text(hjust = 1, angle = 90)) + 
        xlab("Module eigengene") + 
        ylab("Variance explained")

    if (ret){
        return(p)
    } else {
        print(p)
    }
}

#' Plot variance explained by each module eigengene
plot_m_ngene <- function(gene_lab, order_n = TRUE, ret = TRUE){
    require(ggplot2)
    ngene <- tapply(names(gene_lab), gene_lab, length)
    datf <- data.frame("Module" = names(ngene), 
                       "NumGene" = ngene)
    if (order_n){
        o <- order(datf[,"NumGene"], decreasing = TRUE)
        datf <- datf[o,,drop=FALSE]
    }
    datf[,"Module"] <- factor(datf[,"Module"], levels = datf[,"Module"])

    p <- ggplot(datf, aes(x = Module, y = NumGene)) + 
        geom_bar(stat = "identity") + 
        theme_classic() + 
        theme(axis.text.x = element_text(hjust = 1, angle = 90)) + 
        xlab("Module") + 
        ylab("Number of genes")

    if (ret){
        return(p)
    } else {
        print(p)
    }
}

cor_mat <- function(x, y){

    cor_r <- matrix(nrow = ncol(x), ncol = ncol(y))
    rownames(cor_r) <- colnames(x)
    colnames(cor_r) <- colnames(y)
    cor_p <- cor_r

    for (i1 in 1:ncol(x)){
        for (i2 in 1:ncol(y)){
            r <- cor.test(x[,i1], y[,i2], use = "p")
            cor_r[i1, i2] <- r$estimate
            cor_p[i1, i2] <- r$p.value
        }
    }

    return(list("est" = cor_r, "p" = cor_p))
}


#' Correlate module eigengenes with traits
#'
#' @param MEs Data frame of module eigengenes. Samples in rows, MEs 
#'  in columns
#' @param traits Data frame of traits. Samples in rows, traits in columns.
#'
cor_me_trait <- function(MEs, traits){
    si <- intersect(rownames(MEs), rownames(traits))
    MEs <- MEs[si,,drop=FALSE]
    traits <- traits[si,,drop=FALSE]
    cors <- cor_mat(MEs, traits)

    return(cors)
}


