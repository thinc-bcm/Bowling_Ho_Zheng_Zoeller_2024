library(pheatmap)

prop_v2 = read.table(
  '../../data/dtag_experiments/intron_misproc_control_norm.tsv', 
  sep='\t',
  header=T,
  row.names=1)

meta = read.table(
  "../../data/dtag_experiments/intron_heatmap_metadata.tsv", 
  sep='\t',
  header=T,
  row.names=1)

introns <- c(scan("../../data/dtag_experiments/u2_introns_top_100.txt", what="", sep="\n"), 
             scan("../../data/dtag_experiments/catalytic_spliceosome_introns_top_100.txt", what="", sep="\n"), 
             scan("../../data/dtag_experiments/dhx15_introns_top_100.txt", what="", sep="\n"))


ann_colors = list(
  Target = c(
    'AQR'='#ff7f0e',
    'DDX46'='#8c564b',
    'DHX15'='#1f77b4',
    'DHX16'='#17becf',
    'DHX38'='#e377c2',
    'PRPF8'='#2ca02c',
    'SF3B1'='#d62728',
    'U2AF2'='#9467bd'
  ),
  Dosage.value=c(
    'Veh' = '#edf8fb',
    'High dose' = '#2ca25f'
  ),
  Time.value=c(
    '6'= '#edf8fb',
    '9'='#a0b2d4',
    '12'='#8856a7'
  )
)

pheatmap(prop_v2[introns,],
    scale = "row", 
    breaks = seq(-1, 3, length = 100),
    gaps_col = c(54,126),
    gaps_row = c(100,200),
    cluster_rows = F,
    col = viridisLite::magma(100),
    cluster_cols = F,
    annotation_col = meta[colnames(prop_v2), c("Target", "Dosage.value", "Time.value")],
    show_rownames = F, show_colnames = F,
    annotation_colors = ann_colors,
    border_color = NA
)