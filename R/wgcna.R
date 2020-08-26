
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

    return(dissTOM)
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

    if (plots){
        dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
        pfile <- paste0(plot_dir, dendro_plot)
        pdf(pfile, width = 7, height = 5)
        plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05,
                            main = "Gene dendrogram and module colors")
        dev.off()
    }

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
        pdf(pfile, width = 7, height = 5)
        plot(METree, main = "Clustering of module eigengenes", 
             xlab = "", sub = "")
        abline(h = MEDissThres, col = "red")
        dev.off()
    }

    names(merged$colors) <- colnames(x)

    return(merged)
}

