library(data.table)
library(ggplot2)
library(ggrepel)

plot_assoc_gistic <- function(
    path_to_assoc_results = NULL,
    tissue = "pancancer",
    signature = "DHX15",
    gene_labels = NULL,
    gene_focus = NULL) {
    # load association study results
    assoc_res <- fread(path_to_assoc_results)
    # colors for chromosomes
    farben <- rep(c("black", "gray60"), 11)
    names(farben) <- unique(assoc_res$chr)
    assoc_res[, chr := factor(chr, levels = unique(assoc_res$chr))]
    # label genes of interest
    targets <- c(
        "AQR", "DDX46", "DHX15", "DHX16", "DHX38", "PRPF8", "SF3B1", "U2AF2"
    )
    assoc_res[, label := NA_character_]
    if (!is.null(gene_labels)) {
        assoc_res[gene_name %in% gene_labels, label := gene_name]
    } else {
        if (signature %in% targets) {
            assoc_res[gene_name == signature, label := gene_name]
        } else {
            break()
        }
    }
    # subset to chromosome of gene focus window provided by user (if any)
    if (!is.null(gene_focus)) {
        assoc_res <- assoc_res[chr == assoc_res[gene_name == gene_focus[1], seqid]]
    }
    # Manhattan plot
    set.seed(123)
    p <- ggplot(assoc_res, aes(
        x = rank,
        y = association_score,
        color = chr,
        label = label
    )) +
        labs(
            title = paste0(
                "TCGA ", tissue, " - ", signature, " Score Association Study"
            ),
            y = "Association Score",
            x = "Genes ranked by genomic position"
        ) +
        geom_hline(yintercept = 0) +
        geom_point(size = 0.5) +
        ggrepel::geom_label_repel(
            box.padding = 5, max.overlaps = Inf, color = "blue2",
            seed = 0, min.segment.length = 0, nudge_y = 1, na.rm = TRUE
        ) +
        scale_color_manual(values = farben) +
        scale_x_continuous(
            name = "Genes ranked by genomic position",
            breaks = match(1:22, assoc_res$chr),
            labels = 1:22
        ) +
        theme_classic() +
            theme(
                legend.position = "none",
                axis.text = element_text(size = 8),
                axis.title = element_text(size = 16),
                plot.title = element_text(size = 18, hjust = 0.5)
            )
    p
}

