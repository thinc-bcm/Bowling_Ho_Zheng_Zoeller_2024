# Data Description

Description of data files used in Bowling_Ho_Zheng_2022

    diff_misproc_tabs/*.csv - differential intron level misprocessing
    intron_misproc/*_class_counts.tsv - intron level misprocessing
    multigene_misproc/*_class_counts.tsv - multigene level misprocessing
    sj_out/*_SJ.out.classified.tab - STAR junction output with junctions classified

    CSJ_signature_score*score.tsv - cryptic splice junction DHX15 score
    *classified_junctions.tsv - STAR junction output with junctions classified
    exp_introns.tsv - expressed introns for a given set of experiments

    bpp_outs/*ecdf.tsv - branch point predictor outputs

    genomic_sequences/*_canon_seqs.tsv - maxentscan scores for canonical splice sites
    genomic_sequences/*_trt_seqs.tsv - maxentscan scores for cryptic splice sites
    genomic_sequences/*_plus.tsv - gene info
    genomic_sequences/*_do_*_500 - splice_ai probability score for donor
    genomic_sequences/*_acc_*_500 - splice_ai probability score for acceptor
    genomic_sequences/spliceai* - splice_ai predictions (donor and acceptor)

    depmap/ - crpytic splice junction DHX15 score, Copy Number Variation, Demeter2 Dependency score and accompanying metadata for DepMap cell-lines

    eclip_overlaps/ip_peaks_exp_genes/*K562_IDR.txt - set of IDR peaks in genes with TPM > 1 in SUM159 & K562
    eclip_overlaps/summary_tables/*summary.csv - binding frequency heatmaps
    eclip_overlaps/*/*overlaps/*/peak_summary.txt - binding frequency line plots

    tcga_brca/*sj_genome_wide_association.tsv - association score cnv and crpytic splice junction score
    tcga_brca/*cnv.Rdata - TCGA cnv assay and metadata
    tcga_brca/*.csv - crpytic splice junction DHX15 score, known intron count and metadata
