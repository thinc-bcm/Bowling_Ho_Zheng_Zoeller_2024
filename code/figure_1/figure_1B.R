# Load data

filelist = c(
  '/storage/westbrook/projects/splice_dtag/guac_diff_tabs/U2AF2_SUM159_C5_250nM_vs_0nM_dTagV1_9hr.csv',
  '/storage/thinc/data/helicase_misprocessing/DDX46/SUM159/H5/ABSW1-ABSW72_DDX46_dTag13/diff_analysis/intron/DDX46_SUM159_H5_200nM_vs_0nM_dTag13_9hr.csv',
  '/storage/thinc/data/helicase_misprocessing/SF3B1/SUM159/D12/sum159_FKBP-SF3B1/diff_analysis/intron/SF3B1_SUM159_D12_200nM_vs_0nM_dTag13_9.5hr.csv',
  '/storage/westbrook/projects/splice_dtag/guac_diff_tabs/PRPF8_SUM159_G6_250nM_vs_0nM_dTagV1_9hr.csv',
  '/storage/thinc/data/helicase_misprocessing/DHX16/SUM159/H6/FKBP_DHX16_DHX38/diff_analysis/intron/DHX16_SUM159_H6_400nM_vs_0nM_dTag13_9hr.csv',
  '/storage/westbrook/projects/splice_dtag/guac_diff_tabs/AQR_SUM159_B11_500nM_vs_0nM_dTagV1_9hr.csv',
  '/storage/thinc/data/helicase_misprocessing/DHX38/SUM159/D9/FKBP_DHX16_DHX38/diff_analysis/intron/DHX38_SUM159_D9_400nM_vs_0nM_dTag13_9hr.csv',
  '/storage/thinc/data/helicase_misprocessing/DHX15/SUM159/D12/DHX15_D5_D12_G10_dTAG13/diff_analysis/intron/DHX15_SUM159_D12_500nM_vs_0nM_dTag13_9hr.csv'
)

exp_introns = read.delim(
  '/storage/westbrook/projects/splice_dtag/all_target_analysis/exp_introns.tsv',
  header=FALSE
)

# Figure 1B plots
plots <- list()

for (i in 1:length(filelist)){
  
  filnam <- strsplit(filelist[i],'/')[[1]]
  filnam <- filnam[length(filnam)]
  filnam <- strsplit(filnam,split='_', fixed=TRUE)[[1]][1]
  
  test = read.delim(filelist[i],row.names = 1)
  
  test_filt = test[exp_introns$V1,c('binomial_test_adj_p','delta_ratio')]
  
  test_filt[test_filt$binomial_test_adj_p <= 1e-100,'binomial_test_adj_p'] <- 1e-100
  
  plots[[i]] <- ggplot(test_filt, aes(delta_ratio, -log10(binomial_test_adj_p))) +
    geom_hex(bins = 50) +
    xlim(-0.4,0.5) +
    scale_fill_viridis_c(
      trans = 'log', 
      breaks = c(1, 10, 100,1000),
      limits=c(1,1000),
      na.value = "#FDE725FF") +
    xlab("Delta ratio") +
    ylab("-log10 p-value") +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    ggtitle(filnam) +
    theme_classic()
}

final <- gridExtra::grid.arrange(
  plots[[1]],
  plots[[2]],
  plots[[3]],
  plots[[4]],
  plots[[5]],
  plots[[6]],
  plots[[7]],
  plots[[8]], 
  ncol = 4)

ggsave(final, filename = "~/dhx15_uj_paper_figures/all_target_misproc_volcano.pdf",
       width = 16, height = 7)
