
#'
#' @param p is the path to the STAR solo output that contains the matrix.mtx,
#'   features.tsv and barcodes.tsv files
#' @param gene_col is the column index to use for the gene names in the 
#'   barcodes.tsv file
#' @param sep is the separator if gene names are not unique, 
#'   the argument to make.unique function
#'
read_counts <- function(p, gene_col = 1, sep = "."){
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

#' Add feature percentages to SCE object.
#'
#' @param x SCE object.
#' @param ids Features.
#' @param colname Column name to place percent.
sce_feat_pct <- function(x, ids, colname){
    ids <- lapply(ids, function(i) grep(i, rownames(x@counts), value = TRUE))
    ids <- unlist(ids)
    ids <- intersect(ids, rownames(x@counts))
    x <- get_gene_pct(x = x, genes = ids, name = colname)
    return(x)
}

#' Remove features from SCE object.
#'
#' @param x SCE object.
#' @param ids Features.
sce_rm_feat <- function(x, ids){
    ids <- lapply(ids, function(i) grep(i, rownames(x@counts), value = TRUE))
    ids <- unlist(ids)
    keep_ids <- setdiff(rownames(x@counts), ids)
    x@counts <- x@counts[keep_ids,]
    x@gene_data <- x@gene_data[keep_ids,]
    x@droplet_data[,"total_counts"] <- colSums(x@counts[,rownames(x@droplet_data)])
    x@droplet_data[,"n_genes"] <- colSums(x@counts[,rownames(x@droplet_data)] > 0)
    d <- rownames(x@test_data)
    x@test_data[,"total_counts"] <- x@droplet_data[d,"total_counts"]
    x@test_data[,"n_genes"] <- x@droplet_data[d,"n_genes"]
    return(x)
}

#' Initialize DIEM
#'
#' Initialize DIEM SCE object from STARSolo output.
#'
#' @param sample_id Sample ID.
#' @param data_dir Directory prefix to place DIEM SCE object. Appends sample ID.
#' @param star_dir Base directory with STAR results. This looks for the 
#'  results in "solo_dir".
#' @param solo_dir Directory with STARsolo counts.
#' @param pct_feat A names character list. Each element is a character 
#'  vector that contains the feature IDs to calculate the percentage of. 
#'  The names contain the column IDs to place the percentages in.
#' @param genes_rm Character vector containing genes to remove.
#' @param drop_prepend String to prepend droplet IDs, for example the 
#'  sample ID
#' @param drop_delim Delimiter for drop_prepend (";" is the default).
#' @param min_counts Minimum counts for DIEM test set.
#' @param cpm_thresh Only include genes with at least this CPM.
#' @param sce_file File name for SCE object. Will be saved in a directory 
#'  specified by data_dir/sample_id/.
diem.init <- function(sample_id, 
                      data_dir, 
                      star_dir, 
                      solo_dir = "/Solo.out/GeneFull/raw/",
                      pct_feat = NULL, 
                      genes_rm = NULL, 
                      drop_prepend = NULL, 
                      drop_delim = ';', 
                      min_counts = 100, 
                      cpm_thresh = 0, 
                      sce_file = "sce.rds"){

    dir_counts <- file.path(star_dir, solo_dir)
    counts <- read_counts(dir_counts)

    if (!is.null(drop_prepend)){
        colnames(counts) <- paste0(drop_prepend, drop_delim, colnames(counts))
    }

    sce <- create_SCE(counts, name = sample_id)       

    if (!is.null(pct_feat)){
        for (i in names(pct_feat)){
            sce <- sce_feat_pct(sce, pct_feat[[i]], i)
        }
    }

    if (!is.null(genes_rm)){
        sce <- sce_rm_feat(sce, genes_rm)
    }

    sce <- set_debris_test_set(sce, min_counts = 100)
    sce <- filter_genes(sce, cpm_thresh = cpm_thresh)

    sce_dir <- file.path(data_dir, sample_id)
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    sce_f <- file.path(data_dir, sce_file)
    saveRDS(sce, sce_f)
}

#' Cluster with DIEM
#'
#' Cluster DIEM SCE object.
#'
#' @param sample_id Sample ID.
#' @param data_dir Directory prefix to place DIEM SCE object.
#' @param n_expect Expected number of droplets.
#' @param quant Quantile threshold.
#' @param k_init Initial number of clusters K for clustering.
#' @param aprior Add a non-informative prior by adding a count of
#'  ‘alpha_prior’ to all genes in each cluster. Only valid for
#'  the multinomial model.
#' @param pprior Add a non-informative prior by adding a count of ‘pi_prior’
#'  to the each cluster's membership.
#' @param max_iter Maximum number of iterations for the EM estimation of the 
#'  mixture model.
#' @param eps Tolerance threshold for convergence.
#' @param psc Pseudocount.
#' @param n_threads Number of threads.
#' @param sce_file File name for SCE object. Will be saved in a directory 
#'  specified by data_dir/sample_id/.
diem.clust <- function(sample_id, 
                       data_dir, 
                       n_expect = 10e3, 
                       quant = 0.99, 
                       k_init = 30, 
                       aprior = 1, 
                       pprior = 1, 
                       max_iter = 1000, 
                       eps = 1e-4, 
                       psc = 0, 
                       n_threads = 1, 
                       sce_file = "sce.rds"){

    sce_f <- file.path(data_dir, sce_file)
    sce <- readRDS(sce_f)

    drop_dat <- droplet_data(sce)
    n_expect <- min(n_expect, nrow(drop_dat))
    o <- order(drop_dat[,"n_genes"], decreasing = TRUE)
    ix_o <- round(n_expect * (1 - quant))
    ix <- o[ix_o]
    min_g <- drop_dat[ix, "n_genes"] / 10

    sce <- get_pcs(sce, min_genes = min_g)
    sce <- init(sce, threads = n_threads, k_init = k_init)

    sce <- run_em(sce, alpha_prior = aprior, pi_prior = pprior,
                  max_iter = max_iter, eps = eps, psc = psc,
                  threads = n_threads)

    sce <- assign_clusters(sce)

    sce <- estimate_dbr_score(sce, max_genes = 100, thresh_genes = 100)

    saveRDS(sce, sce_f)
}

#' Run DE of DIEM clusters
#' 
#' @param sample_id Sample ID.
#' @param data_dir Directory prefix to with DIEM SCE object.
#' @param marker_dir Directory prefix to place cluster markers.
#' @param marker_file Marker file name.
#' @param scale_factor Scaling factor for normalizing counts.
#' @param sce_file File name for SCE object. Will be saved in a directory 
#'  specified by data_dir/sample_id/.
#'
diem.markers <- function(sample_id,
                         data_dir, 
                         marker_dir, 
                         marker_file = "clust_markers.txt", 
                         scale_factor = 1e3, 
                         sce_file = "sce.rds"){

    sce_dir <- file.path(data_dir, sample_id)
    sce_f <- file.path(sce_dir, sce_file)
    sce <- readRDS(sce_f)

    dir.create(marker_dir, recursive = TRUE, showWarnings = FALSE)

    rc <- raw_counts(sce)
    markers <- de_ttest_all(rc, labels = sce@test_data$Cluster,
                            normalize = TRUE, sf = scale_factor)

    ofn <- paste0(marker_dir, marker_file)
    write.table(markers, ofn, 
                row.names = TRUE, col.names = NA,
                quote = FALSE, sep = "\t")
}

#' Plot DIEM clusters
#'
#' Plot results from DIEM clustering.
#'
#' @param sample_id Sample ID.
#' @param data_dir Directory prefix with DIEM SCE object.
#' @param plot_dir Directory prefix to place plots.
#' @param plot_feat List of character vectors, each of length 2 containing 
#'  x and y variables.
#' @param logt List of boolean vectors, each of length 2 specifying 
#'  whether to log transform the axes.
#' @param sce_file File name for SCE object. Will be saved in a directory 
#'  specified by data_dir/sample_id/.
diem.plot.clust <- function(sample_id, 
                            data_dir,
                            plot_dir, 
                            plot_feat, 
                            logt = NULL, 
                            sce_file = "sce.rds"){
    sce_dir <- file.path(data_dir, sample_id)
    sce_f <- file.path(sce_dir, sce_file)
    sce <- readRDS(sce_f)
    
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    for (i in 1:length(plot_feat)){
        x <- plot_feat[[i]][1]
        y <- plot_feat[[i]][2]
        p <- plot_clust(sce, 
                        feat_x = x, feat_y = y,
                        log_x = logt[[i]][1], log_y = logt[[i]][2]) + 
            ggtitle(sample)
        ofn <- paste0(plot_dir, "cluster.", x, ".", y, ".jpeg")
        ggsave(ofn, width = 5, height = 5)
    }
                   
    # Plot number droplets in clusters
    p <- ggplot(sce@test_data, aes(Cluster)) + 
        geom_bar() + 
        theme_minimal() + 
        xlab("Cluster") + ylab("Number of Droplets") + 
        ggtitle(sample)
    
    ofn <- paste0(plot_prefix, "cluster.n_droplets.jpeg")
    ggsave(ofn, width = 5, height = 5)
}

#' Run DE of DIEM clusters
#' 
#' @param sample_id Sample ID.
#' @param data_dir Directory prefix to with DIEM SCE object.
#' @param bg_clust Background DIEM clusters to remove.
#' @param min_genes Filter out droplets with less than these many genes 
#'  detected.
#' @param seur_dir Directory prefix to place Seurat object.
#' @param seur_file Seurat file name.
#' @param sce_file File name for SCE object. Will be saved in a directory 
#'  specified by data_dir/sample_id/.
#'
diem.filter <- function(sample_id, 
                        data_dir, 
                        bg_clust = c(1), 
                        min_genes = 200, 
                        seur_dir = "./", 
                        seur_file = "seur.rds", 
                        sce_file = "sce.rds"){

    sce_dir <- file.path(data_dir, sample_id)
    sce_f <- file.path(sce_dir, sce_file)
    sce <- readRDS(sce_f)

    md <- sce@test_data

    sdk <- (! td[,"Cluster"] %in% bg_clust ) & ( td[,"n_genes"] >= min_genes )

    k <- rownames(td)[sdk]
    counts <- sce@counts[,k]
    td <- td[k,]

    seur <- Seurat::CreateSeuratObject(counts = counts,
                                       meta.data = md,
                                       names.delim = "}",
                                       project = sample_id)

    dir.create(seur_dir, recursive = TRUE, showWarnings = FALSE)
    ofn <- paste0(seur_dir, seur_file)
    saveRDS(seur, ofn)
}

