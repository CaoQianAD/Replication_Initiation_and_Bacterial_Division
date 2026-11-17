# Baby-machine CCR clustering and enrichment
#
# Purpose: Cluster corrected CCR trajectories and calculate functional enrichment.
# Input: MATLAB outputs in results/baby_machine/BabyMachineDataAnalysis_output;
#        Table S4 and Analysis_annotation_files.
# Output: Cluster assignments, k diagnostics, heatmaps, COG/TF enrichment,
#         and dcw component profiles.
#
# K-means cluster numbers are arbitrary identifiers. For Data S2 display,
# original clusters 4, 2, 3, 1 correspond to display clusters 1, 2, 3, 4.
# Relabeling does not change cluster membership or the clustering calculation.



# 1. Input paths and shared analysis settings
# Paths follow this saved script, independently of the current working directory.
# ID_SCRIPT_DIR optionally overrides the script folder for console execution.
rstudio_source_path <- function() {
  if(!requireNamespace("rstudioapi",quietly=TRUE)) return("")
  if(!rstudioapi::isAvailable()) return("")
  tryCatch({
    p<-rstudioapi::getSourceEditorContext()$path
    if(is.character(p)&&length(p)==1L&&!is.na(p)) p else ""
  },error=function(e) "")
}
locate_analysis_script_dir <- function() {
  override<-Sys.getenv("ID_SCRIPT_DIR",unset="")
  if(nzchar(override)) {
    if(!dir.exists(override)) stop("Script directory does not exist: ",override)
    return(normalizePath(override,winslash="/",mustWork=TRUE))
  }
  sourced<-vapply(sys.frames(),function(f) {
    p<-f$ofile
    if(is.character(p)&&length(p)==1L&&!is.na(p)) p else ""
  },character(1))
  sourced<-sourced[nzchar(sourced)]
  cli<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
  p<-if(length(sourced)) tail(sourced,1) else if(!interactive()&&length(cli)) cli[1] else rstudio_source_path()
  if(nzchar(p)) return(dirname(normalizePath(p,winslash="/",mustWork=TRUE)))
  normalizePath(getwd(),winslash="/",mustWork=TRUE)
}
script_dir<-locate_analysis_script_dir()
# Search the script directory first, then its ancestors. Never choose among
# multiple workbooks silently. TABLE_S4 can select an explicit input file.
ancestors <- function(path) {
  dirs<-path
  for(i in seq_len(8)) { parent<-dirname(tail(dirs,1)); if(parent==tail(dirs,1)) break; dirs<-c(dirs,parent) }
  dirs
}
search_dirs <- ancestors(script_dir)
table_s4 <- Sys.getenv("TABLE_S4", unset="")
if(!nzchar(table_s4)) {
  for(folder in search_dirs) {
    candidates<-list.files(folder,"^TableS4.*[.]xlsx$",full.names=TRUE)
    if(length(candidates)>1L) stop("Multiple Table S4 workbooks in ",folder,". Set TABLE_S4 to the intended file: ",paste(candidates,collapse="; "))
    if(length(candidates)==1L) {table_s4<-candidates[1];break}
  }
  if(!nzchar(table_s4)) stop("Table S4 not found. Source the saved script, run it with Rscript, or set ID_SCRIPT_DIR to its folder. TABLE_S4 can select an explicit workbook. Searched: ",paste(search_dirs,collapse="; "))
}
# Support Windows extended-length paths for deeply nested workbook files.
if(.Platform$OS.type=="windows" && !file.exists(table_s4) && grepl("^[A-Za-z]:[/\\\\]",table_s4))
  table_s4<-paste0("\\\\?\\",gsub("/","\\\\",table_s4))
if(!file.exists(table_s4)) stop("TABLE_S4 does not exist: ",table_s4)
table_s4<-normalizePath(table_s4,winslash="/",mustWork=TRUE)
project_dir <- dirname(table_s4)
annotation_candidates<-file.path(unique(c(script_dir,project_dir,search_dirs)),"Analysis_annotation_files")
annotation_candidates<-annotation_candidates[dir.exists(annotation_candidates)]
annotation_dir<-if(length(annotation_candidates)) annotation_candidates[1] else file.path(project_dir,"Analysis_annotation_files")
results_dir <- file.path(script_dir,"results")
# Numerical settings below are the manuscript analysis defaults.
# Results are written beside this script. Existing inputs are read-only.
analysis_config <- function(standard_root) {
  project <- project_dir
  table <- table_s4
  if (length(table) != 1L) stop("Expected exactly one TableS4 workbook in ", project)
  list(
    standard_root=standard_root, project_root=project, table_s4=table,
    annotation_dir=annotation_dir,
    results_dir=results_dir,
    font="Arial", font_size=9, figure_width=5.2, figure_height=3.8, dpi=300,
    # Changing these paper defaults makes an exploratory run.
    dnaa=list(low_min_tpm=1,low_max_tpm=10,cv_theta_q=.90,
      fdr_threshold=.05,fc_threshold_log2=.585,k_final=4L,seed=123L,
      nstart=50L,iter_max=100L),
    run_original_extended_plots=TRUE,
    plot_all_proteins=tolower(Sys.getenv("ID_PLOT_ALL_PROTEINS", "true")) %in% c("true","1","yes"),
    dnaa_markers=c("appY","rpoH","mioC","iraD","mutH","nrdA","nrdD","fliC"),
    dcw_main=c("mraZ","rsmH","ftsL","ftsI","murE","murF","mraY","murD",
               "ftsW","murG","murC","ddlB"),
    dcw_region=c("mraZ","rsmH","ftsL","ftsI","murE","murF","mraY","murD",
                 "ftsW","murG","murC","ddlB","ftsQ","ftsA","ftsZ"),
    colors=c(hns="#B34A88",zapC="#79ACCD",sdiA="#AF9B69",dcw="#AB94C2"),
    seed=123L
  )
}
STD_ROOT<-script_dir
cfg<-analysis_config(script_dir)

# Shared I/O and plotting only. This file does not change statistical methods.
options(stringsAsFactors=FALSE)
user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib)) .libPaths(unique(c(user_lib,.libPaths())))
default_user_lib <- file.path(Sys.getenv("LOCALAPPDATA"),"R","win-library",
                              paste(R.version$major,strsplit(R.version$minor,".",fixed=TRUE)[[1]][1],sep="."))
if (dir.exists(default_user_lib)) .libPaths(unique(c(default_user_lib,.libPaths())))
profile_user_lib <- file.path(Sys.getenv("USERPROFILE"),"AppData","Local","R","win-library",
                              paste(R.version$major,strsplit(R.version$minor,".",fixed=TRUE)[[1]][1],sep="."))
if (dir.exists(profile_user_lib)) .libPaths(unique(c(profile_user_lib,.libPaths())))

need <- function(pkgs) {
  absent <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
  if(length(absent)) stop("Missing R packages: ",paste(absent,collapse=", "),
    ". Install them in your R library before running. Nothing is installed automatically.")
}
need(c("readxl","ggplot2"))
if(.Platform$OS.type=="windows") {
  grDevices::windowsFonts(Arial=grDevices::windowsFont("Arial"))
}
dir_make <- function(p) { dir.create(p,recursive=TRUE,showWarnings=FALSE); p }
write_csv <- function(x,p) {
  dir_make(dirname(p)); utils::write.csv(x,p,row.names=FALSE,na="",fileEncoding="UTF-8")
}
read_sheet <- function(file,sheet) {
  as.data.frame(readxl::read_excel(file,sheet=sheet,.name_repair="minimal"),
                check.names=FALSE,stringsAsFactors=FALSE)
}
numeric_matrix <- function(d,cols) {
  m <- as.matrix(data.frame(lapply(d[cols],function(v) suppressWarnings(as.numeric(v))),
                            check.names=FALSE))
  colnames(m) <- cols
  m
}
sample_design <- function(headers,protein=FALSE) {
  rx <- if(protein) "^(dnaA|hns)upshift_Rep([0-9]+)_(-?[0-9.]+) min$" else
    "^Replicate([0-9]+)_(-?[0-9.]+) min$"
  cols <- grep(rx,headers,value=TRUE)
  if(!length(cols)) stop("No recognized sample columns.")
  x <- regmatches(cols,regexec(rx,cols))
  if(protein) data.frame(sample=cols,experiment=vapply(x,"[","",2),
    replicate=as.integer(vapply(x,"[","",3)),time_min=as.numeric(vapply(x,"[","",4))) else
    data.frame(sample=cols,replicate=as.integer(vapply(x,"[","",2)),
      time_min=as.numeric(vapply(x,"[","",3)))
}
mean_ok <- function(x) { x<-x[is.finite(x)]; if(length(x)) mean(x) else NA_real_ }
sd_ok <- function(x) { x<-x[is.finite(x)]; if(length(x)>1L) stats::sd(x) else NA_real_ }
relative_by_pre <- function(m,design,require_complete=TRUE) {
  rel<-m; rel[,]<-NA_real_; audits<-list()
  for(rep in unique(design$replicate)) {
    dd<-design[design$replicate==rep,,drop=FALSE]; pre<-dd$sample[dd$time_min<0]
    if(!length(pre)) stop("No pre-upshift samples in replicate ",rep)
    pm<-m[,pre,drop=FALSE]; n<-rowSums(is.finite(pm)&pm>0)
    b<-apply(pm,1,mean_ok)
    good<-is.finite(b)&b>0
    if(require_complete) good<-good & n==length(pre)
    rel[good,dd$sample]<-m[good,dd$sample,drop=FALSE]/b[good]
    audits[[length(audits)+1L]]<-data.frame(row_index=seq_len(nrow(m)),replicate=rep,
      pre_times_min=paste(dd$time_min[dd$time_min<0],collapse=";"),
      n_pre_expected=length(pre),n_pre_valid=n,baseline_PPM=b,baseline_valid=good)
  }
  list(relative=rel,audit=do.call(rbind,audits))
}
paper_theme <- function(cfg) {
  ggplot2::theme_classic(base_size=cfg$font_size,base_family=cfg$font)+
    ggplot2::theme(axis.line=ggplot2::element_line(linewidth=.35),
      axis.ticks=ggplot2::element_line(linewidth=.35),
      legend.title=ggplot2::element_text(size=cfg$font_size),plot.title=ggplot2::element_text(size=cfg$font_size),
      strip.background=ggplot2::element_blank())
}
save_panel <- function(p,name,out,cfg,width=cfg$figure_width,height=cfg$figure_height) {
  dir_make(out); p<-p+paper_theme(cfg)
  ggplot2::ggsave(file.path(out,paste0(name,".pdf")),p,width=width,height=height,
                 device=grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(out,paste0(name,".png")),p,width=width,height=height,dpi=cfg$dpi)
}
safe_name <- function(s) substr(gsub("[^a-zA-Z0-9_.-]","_",s),1,95)
# 2. Clustering, heatmaps and functional enrichment
# Clustering uses corrected expression trajectories. The within-cluster row
# order uses the phase reference replicate_time. Peak-time summaries use the
# first fitted peak peak_time1, including the source model's endpoint fallback.


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(cluster)
  library(cowplot)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(scales)
  library(grid)
})

## ================================================================
## 1. Paths and user settings
DATA_DIR <- file.path(cfg$results_dir,"baby_machine")
PATH_CCR <- file.path(DATA_DIR,"BabyMachineDataAnalysis_output","CCR_genes.csv")
PATH_BG <- file.path(DATA_DIR,"BabyMachineDataAnalysis_output","corrected_TPM_total_M_Ala.csv")
PATH_EXPRESSION_FILTER <- file.path(DATA_DIR,"BabyMachineDataAnalysis_output","expression_filter_summary.csv")
ANNOT_DIR <- cfg$annotation_dir
PATH_COG_ANN <- file.path(ANNOT_DIR,"250715_Ecoli_K12_COG_UniqueGene_Merged.csv")
PATH_COG_FUN <- file.path(ANNOT_DIR,"250717_cog-24.fun.csv")
PATH_ECO_TXT <- file.path(ANNOT_DIR,"260113_Ecocyc_New-SmartTable_allColumn_RUsed.txt")
PATH_ECO_TF <- file.path(ANNOT_DIR,"20260113_ManualChecked_Ecocyc_TF.csv")
OUTDIR <- dir_make(file.path(DATA_DIR,"ccr"))
DIR_TABLES <- file.path(OUTDIR, "tables")
DIR_PLOTS <- file.path(OUTDIR, "plots")
DIR_PER_CLUSTER <- file.path(DIR_PLOTS, "per_cluster")
DIR_COG <- file.path(OUTDIR, "COG_enrichment", "by_cluster")
DIR_TF <- file.path(OUTDIR, "TF_enrichment")
dir.create(DIR_TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PLOTS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PER_CLUSTER, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_COG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TF, recursive = TRUE, showWarnings = FALSE)

# Remove obsolete duplicate clustering plots from earlier script versions.
# This keeps the clustering output folder unambiguous after rerunning the script.
obsolete_cluster_outputs <- c(
  file.path(DIR_PLOTS, paste0(c("02b_cluster_means_zscore",
                                "02c_class_average_lines_zscore",
                                "02_cluster_mean_trajectories_zscore",
                                "03_class_average_lines_zscore",
                                "03_clustered_heatmap_zscore",
                                "03_overview_heatmap_zscore",
                                "06_overview_heatmap_zscore"), ".pdf")),
  file.path(DIR_PLOTS, paste0(c("02b_cluster_means_zscore",
                                "02c_class_average_lines_zscore",
                                "02_cluster_mean_trajectories_zscore",
                                "03_class_average_lines_zscore",
                                "03_clustered_heatmap_zscore",
                                "03_overview_heatmap_zscore",
                                "06_overview_heatmap_zscore"), ".png")),
  file.path(DIR_PER_CLUSTER, paste0("cluster_zscore_", seq_len(20), ".pdf")),
  file.path(DIR_PER_CLUSTER, paste0("cluster_zscore_", seq_len(20), ".png")),
  file.path(DIR_PER_CLUSTER, paste0("cluster", seq_len(20), "_zscore_trajectories.pdf")),
  file.path(DIR_PER_CLUSTER, paste0("cluster", seq_len(20), "_zscore_trajectories.png"))
)
invisible(file.remove(obsolete_cluster_outputs[file.exists(obsolete_cluster_outputs)]))

TIME_POINTS <- seq(0, 220, by = 11)
TIME_COLS <- as.character(TIME_POINTS)

K_FINAL <- 4
K_SCAN <- 2:10
SEED <- 123
KMEANS_NSTART <- 50
KMEANS_ITERMAX <- 100

# Cluster labels returned by k-means are arbitrary. After rerunning the upstream
# low-expression filter or changing K_FINAL, inspect the cluster-average
# trajectories and update this vector. For the manuscript k = 4 analysis, the
# intended display order is early decline, early induction, mid-cycle induction,
# and late-cycle induction. Exploratory values of k are left without biological
# temporal labels until their trajectories have been inspected.
ALL_CLUSTERS <- paste0("cluster", seq_len(K_FINAL))
DISPLAY_CLUSTER_ORDER <- c("cluster4", "cluster2", "cluster3", "cluster1")
TRAJECTORY_CLUSTER_ORDER <- DISPLAY_CLUSTER_ORDER
HEATMAP_CLUSTER_ORDER <- DISPLAY_CLUSTER_ORDER
COG_CLUSTER_ORDER <- DISPLAY_CLUSTER_ORDER
TF_CLUSTER_ORDER <- DISPLAY_CLUSTER_ORDER
 

# Cluster colors used by trajectory plots and heatmap cluster annotations.
# Edit the hex codes to set manuscript-specific colors.
CLUSTER_COLORS <- c(
  cluster4 = rgb(119, 148, 229, maxColorValue = 255),#蓝 
  cluster2 = rgb(241, 102, 103, maxColorValue = 255),#红
  cluster3 = rgb(229, 168,  75, maxColorValue = 255),#黄 
  cluster1 = rgb(129, 164, 115, maxColorValue = 255)#绿 
  )

Q_CUTOFF <- 0.05
M_MIN_TF <- 5
TF_UNIVERSE <- "both"  # "U1", "U2", or "both"
TOP_N_PER_CLUSTER <- 5

FIG_DPI <- 300
FONT_FAMILY <- cfg$font
BASE_SIZE <- cfg$font_size
Z_RANGE <- c(-3, 3)

PALETTE_HEAT <- colorRampPalette(c("#3B4CC0", "#88A1D2", "#F7F7F7", "#F3A982", "#B40426"))(101)
PAL_CLUSTER <- CLUSTER_COLORS

validate_cluster_order <- function(order_vec, order_name) {
  if (!setequal(order_vec, ALL_CLUSTERS) || length(order_vec) != length(ALL_CLUSTERS)) {
    stop(order_name, " must contain each cluster exactly once: ",
         paste(ALL_CLUSTERS, collapse = ", "))
  }
}
validate_cluster_order(TRAJECTORY_CLUSTER_ORDER, "TRAJECTORY_CLUSTER_ORDER")
validate_cluster_order(HEATMAP_CLUSTER_ORDER, "HEATMAP_CLUSTER_ORDER")
validate_cluster_order(COG_CLUSTER_ORDER, "COG_CLUSTER_ORDER")
validate_cluster_order(TF_CLUSTER_ORDER, "TF_CLUSTER_ORDER")
if (!setequal(names(CLUSTER_COLORS), ALL_CLUSTERS) || length(CLUSTER_COLORS) != length(ALL_CLUSTERS)) {
  stop("CLUSTER_COLORS must define one named color for each cluster: ",
       paste(ALL_CLUSTERS, collapse = ", "))
}
COLORS_DIV <- c("#3B4CC0", "#88A1D2", "#F7F7F7", "#F3A982", "#B40426")
NS_PT_COL <- "grey70"
NS_TILE_COL <- "grey90"

theme_pub <- function(base_size = BASE_SIZE) {
  theme_bw(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, colour = "grey35"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      strip.background = element_rect(fill = "grey95", colour = "grey80"),
      legend.position = "right"
    )
}

theme_enrichment <- function(base_size = BASE_SIZE) {
  theme_bw(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, colour = "grey35"),
      panel.grid.major = element_line(colour = "grey86", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      legend.position = "right"
    )
}

save_plot <- function(path_base, plot, width, height) {
  ggsave(paste0(path_base, ".pdf"), plot, width = width, height = height, device = grDevices::cairo_pdf)
  ggsave(paste0(path_base, ".png"), plot, width = width, height = height, dpi = FIG_DPI)
}

set.seed(SEED)

## ================================================================
## 2. Read CCR genes and construct the clustering matrix
## ================================================================

message("[01] Reading baby-machine CCR input tables")
if (!file.exists(PATH_CCR)) stop("CCR table not found: ", PATH_CCR)
if (!file.exists(PATH_BG)) stop("Background corrected TPM table not found: ", PATH_BG)

ccr_df <- fread(PATH_CCR, check.names = FALSE) %>% as.data.frame()
bg_df <- fread(PATH_BG, check.names = FALSE) %>% as.data.frame()
expression_filter_summary <- if (file.exists(PATH_EXPRESSION_FILTER)) {
  fread(PATH_EXPRESSION_FILTER, check.names = FALSE) %>% as.data.frame()
} else {
  NULL
}
if (!is.null(expression_filter_summary)) {
  fwrite(expression_filter_summary,
         file.path(DIR_TABLES, "upstream_expression_filter_summary.csv"))
  message(sprintf(
    "  Upstream expression filter retained %s / %s genes",
    expression_filter_summary$n_retained_genes[1],
    expression_filter_summary$n_input_genes[1]
  ))
} else {
  warning("Upstream expression-filter summary was not found: ", PATH_EXPRESSION_FILTER)
}

required_cols <- c("ID", "Gene name", TIME_COLS)
missing_cols <- setdiff(required_cols, colnames(ccr_df))
if (length(missing_cols) > 0) {
  stop("CCR table is missing required columns: ", paste(missing_cols, collapse = ", "))
}

ccr_df <- ccr_df %>%
  mutate(
    ID = as.character(ID),
    gene = as.character(`Gene name`),
    gene = if_else(!is.na(gene) & gene != "", gene, ID),
    row_key = make.unique(paste0(ID, "|", gene))
  )

expr_corrected <- as.matrix(ccr_df[, TIME_COLS])
storage.mode(expr_corrected) <- "numeric"
rownames(expr_corrected) <- ccr_df$row_key
colnames(expr_corrected) <- paste0("T", TIME_POINTS)

M_zscore <- t(scale(t(expr_corrected), center = TRUE, scale = TRUE))
M_zscore[!is.finite(M_zscore)] <- 0

# Relative corrected TPM is used only for visualizing cluster-average
# fold changes. Clustering itself remains based on z-score(corrected TPM).
first_tp_col <- paste0("T", TIME_POINTS[1])
baseline_tpm <- expr_corrected[, first_tp_col, drop = TRUE]
M_rel_tpm <- sweep(expr_corrected, 1, baseline_tpm, "/")
M_rel_tpm[!is.finite(M_rel_tpm)] <- NA_real_

fwrite(data.frame(row_key = rownames(M_zscore), M_zscore, check.names = FALSE),
       file.path(DIR_TABLES, "input_matrix_zscore_corrected_TPM.csv"))
fwrite(data.frame(row_key = rownames(M_rel_tpm), M_rel_tpm, check.names = FALSE),
       file.path(DIR_TABLES, "input_matrix_relative_corrected_TPM_T0eq1.csv"))
message(sprintf("  CCR genes retained for clustering: %d", nrow(M_zscore)))

## ================================================================
## 3. k-scan and cluster-number assessment
## ================================================================

message("[02] Scanning k values")
scan_kmeans_metrics <- function(mat, k_range) {
  d_mat <- dist(mat)
  bind_rows(lapply(k_range, function(k) {
    set.seed(SEED)
    km <- kmeans(mat, centers = k, nstart = KMEANS_NSTART, iter.max = KMEANS_ITERMAX)
    sizes <- table(km$cluster)
    sil <- mean(silhouette(km$cluster, d_mat)[, "sil_width"])
    tibble(
      k = k,
      wss = km$tot.withinss,
      silhouette = sil,
      max_size_frac = max(sizes) / sum(sizes),
      min_size = min(sizes),
      n_clusters = length(sizes)
    )
  }))
}

estimate_elbow_k <- function(k_vals, wss_vals) {
  x <- (k_vals - min(k_vals)) / (max(k_vals) - min(k_vals))
  y <- (wss_vals - min(wss_vals)) / (max(wss_vals) - min(wss_vals))
  start <- c(x[1], y[1])
  end <- c(x[length(x)], y[length(y)])
  line_vec <- end - start
  dist_to_line <- abs(line_vec[1] * (start[2] - y) - (start[1] - x) * line_vec[2]) /
    sqrt(sum(line_vec^2))
  k_vals[which.max(dist_to_line)]
}

k_metrics <- scan_kmeans_metrics(M_zscore, K_SCAN)
k_metrics <- k_metrics %>%
  mutate(
    silhouette_rank = rank(-silhouette, ties.method = "min"),
    size_balanced = max_size_frac <= 0.60 & min_size >= 5
  )
elbow_k <- estimate_elbow_k(k_metrics$k, k_metrics$wss)
silhouette_best_k <- k_metrics$k[which.max(k_metrics$silhouette)]
current_row <- k_metrics %>% filter(k == K_FINAL)
current_supported <- nrow(current_row) == 1 &&
  (current_row$silhouette_rank <= 2 || K_FINAL == elbow_k) &&
  current_row$size_balanced

k_assessment <- tibble(
  selected_K = K_FINAL,
  silhouette_best_K = silhouette_best_k,
  elbow_K = elbow_k,
  selected_K_silhouette_rank = current_row$silhouette_rank,
  selected_K_silhouette = current_row$silhouette,
  selected_K_max_cluster_fraction = current_row$max_size_frac,
  selected_K_min_cluster_size = current_row$min_size,
  selected_K_size_balanced = current_row$size_balanced,
  selected_K_supported_by_current_rules = current_supported
)

fwrite(k_metrics, file.path(DIR_TABLES, "k_scan_metrics.csv"))
fwrite(k_assessment, file.path(DIR_TABLES, "k_selection_summary.csv"))
writeLines(c(
  "Cluster-number assessment for baby-machine CCR genes",
  sprintf("Main matrix: z-score(corrected TPM), %d CCR genes x %d time points", nrow(M_zscore), ncol(M_zscore)),
  sprintf("Selected K_FINAL: %d", K_FINAL),
  sprintf("Best k by silhouette: %d", silhouette_best_k),
  sprintf("Approximate elbow k: %d", elbow_k),
  sprintf("Selected k silhouette rank: %s", current_row$silhouette_rank),
  sprintf("Selected k silhouette: %.4f", current_row$silhouette),
  sprintf("Selected k largest cluster fraction: %.3f", current_row$max_size_frac),
  sprintf("Selected k smallest cluster size: %d", current_row$min_size),
  sprintf("Selected k supported by scripted rules: %s", current_supported),
  "",
  "Interpretation rule used here:",
  "The selected k is considered supported when it is either the best/second-best by silhouette or matches the WSS elbow, and when cluster sizes are not dominated by one cluster.",
  "This rule is a diagnostic aid; final biological interpretation should also inspect the cluster mean trajectories and functional enrichment."
), file.path(DIR_TABLES, "k_selection_summary.txt"))

make_metric_panel <- function(df, ycol, ylab, ttl, log_y = FALSE) {
  p <- ggplot(df, aes(x = k, y = .data[[ycol]])) +
    geom_line(colour = "#4393C3", linewidth = 0.7) +
    geom_point(size = 2, colour = "#2166AC") +
    geom_vline(xintercept = K_FINAL, linetype = "dashed", colour = "firebrick", linewidth = 0.5) +
    geom_vline(xintercept = elbow_k, linetype = "dotted", colour = "grey35", linewidth = 0.5) +
    scale_x_continuous(breaks = K_SCAN) +
    labs(title = ttl, x = "k (number of clusters)", y = ylab) +
    theme_pub()
  if (log_y) p <- p + scale_y_log10()
  p
}

p_kscan <- plot_grid(
  make_metric_panel(k_metrics, "wss", "Total within-cluster SS", "Elbow plot (WSS)", TRUE),
  make_metric_panel(k_metrics, "silhouette", "Average silhouette width", "Silhouette"),
  make_metric_panel(k_metrics, "max_size_frac", "Largest cluster / N", "Cluster size balance"),
  make_metric_panel(k_metrics, "min_size", "Smallest cluster size", "Smallest cluster"),
  ncol = 2, align = "hv"
)
save_plot(file.path(DIR_PLOTS, "01_k_scan_diagnostics"), p_kscan, width = 9, height = 8)

## ================================================================
## 4. Final k-means clustering
## ================================================================

message("[03] Running final k-means clustering")
set.seed(SEED)
km <- kmeans(M_zscore, centers = K_FINAL, nstart = KMEANS_NSTART, iter.max = KMEANS_ITERMAX)
row_keys <- rownames(M_zscore)
cluster_vec <- setNames(paste0("cluster", km$cluster), row_keys)

meta_cols <- intersect(
  c("ID", "Gene name", "Start", "End", "Relative position", "R^2_sin",
    "p2", "p3", "p6", "replicate_time", "peak_time1", "peak_time2", "FC1_Sim", "FC2_Sim", "SSR_sin/SSR_lin"),
  colnames(ccr_df)
)

cluster_assignments <- tibble(
  row_key = row_keys,
  cluster = cluster_vec[row_keys]
) %>%
  left_join(ccr_df %>% select(row_key, any_of(meta_cols)), by = "row_key") %>%
  mutate(
    phase_reference_time = suppressWarnings(as.numeric(replicate_time)),
    peak_time = suppressWarnings(as.numeric(peak_time1)),
    max_simulated_FC = pmax(FC1_Sim, FC2_Sim, na.rm = TRUE)
  ) %>%
  arrange(factor(cluster, levels = TRAJECTORY_CLUSTER_ORDER), phase_reference_time, desc(max_simulated_FC))

fwrite(cluster_assignments, file.path(DIR_TABLES, "cluster_assignments.csv"))
fwrite(data.frame(cluster = ALL_CLUSTERS, km$centers, check.names = FALSE),
       file.path(DIR_TABLES, "kmeans_centers_zscore.csv"))

cluster_summary <- cluster_assignments %>%
  group_by(cluster) %>%
  summarise(
    n_genes = n(),
    median_peak_time_min = median(peak_time, na.rm = TRUE),
    median_period_min = median(p2, na.rm = TRUE),
    median_abs_damping = median(abs(p6), na.rm = TRUE),
    median_max_simulated_FC = median(max_simulated_FC, na.rm = TRUE),
    representative_genes = paste(head(`Gene name`, 10), collapse = ";"),
    .groups = "drop"
  )
fwrite(cluster_summary, file.path(DIR_TABLES, "cluster_summary.csv"))

temporal_labels <- if (K_FINAL == 4) {
  c(
    "early-cycle decline",
    "early-cycle induction",
    "mid-cycle induction",
    "late-cycle induction"
  )
} else {
  rep(NA_character_, K_FINAL)
}

display_cluster_mapping <- tibble(
  display_class = paste0("C", seq_along(DISPLAY_CLUSTER_ORDER)),
  cluster = DISPLAY_CLUSTER_ORDER,
  temporal_label = temporal_labels
) %>%
  left_join(cluster_summary, by = "cluster")
fwrite(display_cluster_mapping,
       file.path(DIR_TABLES, "display_cluster_mapping.csv"))
message(
  "  Wrote display_cluster_mapping.csv. Verify the temporal order after ",
  "each rerun because k-means cluster IDs are arbitrary."
)

cluster_gene_lists <- cluster_assignments %>%
  group_by(cluster) %>%
  summarise(
    n_genes = n(),
    genes = paste(`Gene name`, collapse = ";"),
    IDs = paste(ID, collapse = ";"),
    .groups = "drop"
  )
fwrite(cluster_gene_lists, file.path(DIR_TABLES, "cluster_gene_lists.csv"))

## ================================================================
## 5. Cluster plots in the dnaA-upshift style
## ================================================================

message("[04] Plotting cluster trajectories and heatmaps")
make_long <- function(mat, clusters) {
  as.data.frame(mat) %>%
    rownames_to_column("row_key") %>%
    pivot_longer(-row_key, names_to = "time_label", values_to = "value") %>%
    mutate(
      Time = as.numeric(str_remove(time_label, "^T")),
      cluster = factor(clusters[row_key], levels = TRAJECTORY_CLUSTER_ORDER)
    )
}

long_z <- make_long(M_zscore, cluster_vec)
long_rel <- make_long(M_rel_tpm, cluster_vec) %>%
  filter(is.finite(value))

summarise_cluster_timecourse <- function(long_df) {
  agg_df <- long_df %>%
    group_by(cluster, Time) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      sd_val = sd(value, na.rm = TRUE),
      n_genes = dplyr::n_distinct(row_key),
      se_val = sd_val / sqrt(pmax(n_genes, 1)),
      t_crit = qt(0.975, df = pmax(n_genes - 1, 1)),
      ci_lo = mean_val - t_crit * se_val,
      ci_hi = mean_val + t_crit * se_val,
      .groups = "drop"
    )
  agg_df
}

rel_cluster_mean_df <- summarise_cluster_timecourse(long_rel)
rel_cluster_mean_df <- rel_cluster_mean_df %>%
  mutate(cluster = factor(cluster, levels = TRAJECTORY_CLUSTER_ORDER))
rel_y_limits <- range(c(rel_cluster_mean_df$ci_lo, rel_cluster_mean_df$ci_hi, 1), na.rm = TRUE)
if (!all(is.finite(rel_y_limits)) || diff(rel_y_limits) <= 0) rel_y_limits <- c(0.5, 1.5)
rel_y_breaks <- pretty(rel_y_limits, n = 4)
rel_x_limits <- range(TIME_POINTS)
rel_x_breaks <- seq(min(TIME_POINTS), max(TIME_POINTS), by = 55)

plot_cluster_means_rel <- function(mean_df) {
  plot_df <- mean_df %>%
    mutate(cluster_lab = sprintf("%s (n=%d)", cluster, n_genes))

  ggplot(plot_df, aes(x = Time, y = mean_val, colour = cluster, fill = cluster)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.4) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    scale_colour_manual(values = PAL_CLUSTER, drop = FALSE) +
    scale_fill_manual(values = PAL_CLUSTER, drop = FALSE) +
    scale_x_continuous(limits = rel_x_limits, breaks = rel_x_breaks, expand = c(0, 0)) +
    scale_y_continuous(limits = rel_y_limits, breaks = rel_y_breaks, expand = c(0, 0)) +
    facet_wrap(~ cluster_lab, ncol = 2, scales = "fixed") +
    labs(title = sprintf("CCR cluster mean trajectories - relative corrected TPM, k = %d", K_FINAL),
         subtitle = sprintf("Clusters were assigned using z-score(corrected TPM); first time point = 1"),
         x = "Time (min)", y = "Relative corrected TPM (T0 = 1)",
         colour = "Cluster", fill = "Cluster") +
    theme_pub() +
    theme(legend.position = "bottom",
          aspect.ratio = 1)
}

p_means <- plot_cluster_means_rel(rel_cluster_mean_df)
save_plot(file.path(DIR_PLOTS, "02_cluster_mean_trajectories_relTPM"), p_means, width = 8.2, height = 8.2)
fwrite(rel_cluster_mean_df, file.path(DIR_TABLES, "cluster_mean_trajectories_relative_corrected_TPM_T0eq1.csv"))

plot_class_average_lines <- function(mean_df) {
  ggplot(mean_df, aes(x = Time, y = mean_val, colour = cluster, fill = cluster)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    scale_colour_manual(values = PAL_CLUSTER, drop = FALSE) +
    scale_fill_manual(values = PAL_CLUSTER, drop = FALSE) +
    scale_x_continuous(limits = rel_x_limits, breaks = rel_x_breaks, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, by =0.5), expand = c(0, 0))+
    labs(title = sprintf("Class-average CCR trajectories - relative corrected TPM, k = %d", K_FINAL),
         subtitle = "Clusters were assigned using z-score(corrected TPM)",
         x = "Time (min)", y = "Mean relative corrected TPM (T0 = 1)",
         colour = "Cluster", fill = "Cluster") +
    theme_pub() +
    theme(legend.position = "right",
          aspect.ratio = 1)
}

p_class_average <- plot_class_average_lines(rel_cluster_mean_df)
save_plot(file.path(DIR_PLOTS, "03_class_average_lines_relTPM"), p_class_average, width = 7.0 / 2.54, height = 5.04 / 2.54)

plot_one_cluster_trajectory <- function(cluster_id) {
  keys <- names(cluster_vec)[cluster_vec == cluster_id]
  n_genes <- length(keys)
  long_c <- long_z %>% filter(row_key %in% keys)
  mean_c <- long_c %>% group_by(Time) %>% summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")

  ggplot() +
    geom_line(data = long_c, aes(Time, value, group = row_key),
              colour = "grey70", linewidth = 0.2, alpha = 0.35) +
    geom_line(data = mean_c, aes(Time, mean_val),
              colour = PAL_CLUSTER[cluster_id], linewidth = 1.0) +
    geom_point(data = mean_c, aes(Time, mean_val),
               colour = PAL_CLUSTER[cluster_id], size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    labs(title = sprintf("%s | n = %d", cluster_id, n_genes),
         x = "Time (min)", y = "z-score(corrected TPM)") +
    theme_pub()
}

plot_one_cluster_detail <- function(cluster_id) {
  keys <- names(cluster_vec)[cluster_vec == cluster_id]
  n_genes <- length(keys)
  long_c <- long_z %>% filter(row_key %in% keys)
  mean_c <- long_c %>% group_by(Time) %>% summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")
  meta_c <- cluster_assignments %>% filter(row_key %in% keys)

  p1 <- plot_one_cluster_trajectory(cluster_id) +
    labs(title = "Individual and mean trajectories")

  p2 <- ggplot(meta_c, aes(x = peak_time)) +
    geom_histogram(binwidth = 10, fill = PAL_CLUSTER[cluster_id], colour = "white", linewidth = 0.25) +
    geom_vline(xintercept = median(meta_c$peak_time, na.rm = TRUE),
               linetype = "dashed", colour = "grey30", linewidth = 0.5) +
    labs(title = "Fitted peak-time distribution", x = "Peak time (min)", y = "Number of genes") +
    theme_pub()

  p3 <- meta_c %>%
    select(row_key, FC1_Sim, FC2_Sim) %>%
    pivot_longer(c(FC1_Sim, FC2_Sim), names_to = "cycle", values_to = "simulated_FC") %>%
    ggplot(aes(x = cycle, y = simulated_FC)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, fill = "grey92", colour = "black", linewidth = 0.35) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.45, colour = PAL_CLUSTER[cluster_id]) +
    labs(title = "Simulated fold-change distribution", x = NULL, y = "Simulated FC") +
    theme_pub()

  p4 <- ggplot(meta_c, aes(x = abs(p6), y = `SSR_sin/SSR_lin`)) +
    geom_point(size = 1.4, alpha = 0.65, colour = PAL_CLUSTER[cluster_id]) +
    geom_hline(yintercept = 0.7, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    labs(title = "Model-filter parameters", x = "|p6|", y = "SSR_sin / SSR_lin") +
    theme_pub()

  title_g <- ggdraw() +
    draw_label(sprintf("%s | z-score cluster | n = %d", cluster_id, n_genes),
               fontface = "bold", size = 11, x = 0.5, hjust = 0.5)
  combined <- plot_grid(p1, p2, p3, p4, ncol = 2, align = "hv", labels = c("A", "B", "C", "D"))
  plot_grid(title_g, combined, ncol = 1, rel_heights = c(0.08, 0.92))
}

for (cid in TRAJECTORY_CLUSTER_ORDER) {
  cid_num <- str_remove(cid, "^cluster")
  save_plot(file.path(DIR_PER_CLUSTER, sprintf("cluster_detail_zscore_%s", cid_num)),
            plot_one_cluster_detail(cid), width = 10, height = 8.5)
}

plot_overview_heatmap <- function(mat, clusters_vec, value_range, title_text, out_base) {
  order_df <- cluster_assignments %>%
    mutate(cluster = factor(cluster, levels = HEATMAP_CLUSTER_ORDER)) %>%
    arrange(cluster, phase_reference_time, desc(max_simulated_FC)) %>%
    select(row_key, cluster, ID, `Gene name`, peak_time, max_simulated_FC)
  gene_order <- order_df$row_key
  mat_ord <- mat[gene_order, , drop = FALSE]
  mat_ord <- pmax(pmin(mat_ord, value_range[2]), value_range[1])
  cluster_factor <- factor(clusters_vec[gene_order], levels = HEATMAP_CLUSTER_ORDER)

  fwrite(order_df, file.path(DIR_TABLES, "clustered_heatmap_gene_order.csv"))
  fwrite(data.frame(row_key = rownames(mat_ord), mat_ord, check.names = FALSE),
         file.path(DIR_TABLES, "clustered_heatmap_zscore_matrix.csv"))

  col_fun <- colorRamp2(seq(value_range[1], value_range[2], length.out = length(PALETTE_HEAT)), PALETTE_HEAT)
  ha_left <- rowAnnotation(
    Cluster = cluster_factor,
    col = list(Cluster = PAL_CLUSTER),
    width = unit(5, "mm"),
    annotation_name_gp = gpar(fontsize = BASE_SIZE, fontfamily = FONT_FAMILY)
  )
  
  time_labels <- rep("", length(TIME_POINTS))
  time_labels[TIME_POINTS %in% c(0, 110, 220)] <- as.character(TIME_POINTS[TIME_POINTS %in% c(0, 110, 220)])
  
  sample_tick_anno <- HeatmapAnnotation(
    sampling = anno_points(
      rep(1, length(TIME_POINTS)),
      pch = 124,
      size = unit(2, "mm"),
      gp = gpar(col = "black"),
      axis = FALSE
    ),
    show_annotation_name = FALSE,
    height = unit(3, "mm")
  )
  
  ht <- Heatmap(
    mat_ord,
    name = "z-score",
    col = col_fun,
    na_col = "grey90",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_labels = time_labels,
    column_names_side = "bottom",
    column_names_rot = 0,
    bottom_annotation = sample_tick_anno,
    column_title = title_text,
    column_title_gp = gpar(fontsize = BASE_SIZE + 1, fontfamily = FONT_FAMILY, fontface = "bold"),
    row_split = cluster_factor,
    row_gap = unit(1, "mm"),
    border = TRUE,
    heatmap_legend_param = list(at = c(value_range[1], 0, value_range[2]))
  )

  grDevices::cairo_pdf(paste0(out_base, ".pdf"), width = 9, height = 10, family = FONT_FAMILY)
  draw(ha_left + ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  png(paste0(out_base, ".png"), width = 9, height = 10, units = "in", res = FIG_DPI)
  draw(ha_left + ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

plot_overview_heatmap(
  M_zscore, cluster_vec, Z_RANGE,
  sprintf("CCR genes - z-score(corrected TPM), k = %d", K_FINAL),
  file.path(DIR_PLOTS, "04_overview_heatmap_zscore")
)

## ================================================================
## 6. COG enrichment by cluster
## ================================================================

message("[05] Running COG enrichment by cluster")
required_annot <- c(PATH_COG_ANN, PATH_COG_FUN)
if (!all(file.exists(required_annot))) {
  warning("COG annotation files not found. Skipping COG enrichment.")
} else {
  cog_raw <- read_csv(PATH_COG_ANN, show_col_types = FALSE)
  cog_fun <- read_csv(PATH_COG_FUN, show_col_types = FALSE)
  colnames(cog_fun) <- str_trim(colnames(cog_fun))

  gene_cate_long <- cog_raw %>%
    select(Gene, COG_Category) %>%
    filter(!is.na(Gene), !is.na(COG_Category), COG_Category != "") %>%
    separate_rows(COG_Category, sep = ";") %>%
    mutate(Gene = str_trim(as.character(Gene)),
           COG_Category = str_trim(COG_Category)) %>%
    filter(Gene != "", COG_Category != "") %>%
    distinct(Gene, COG_Category)

  needed_cols <- c("Functional category ID", "Functional category description")
  if (!all(needed_cols %in% colnames(cog_fun))) {
    stop("COG function table is missing required columns: ", paste(setdiff(needed_cols, colnames(cog_fun)), collapse = ", "))
  }
  cate_order_df <- cog_fun %>%
    mutate(Order = row_number()) %>%
    transmute(
      COG_Category = str_trim(`Functional category ID`),
      COG_Desc = `Functional category description`,
      Order
    )

  bg_genes <- bg_df$`Gene name` %>% as.character() %>% unique() %>% intersect(unique(gene_cate_long$Gene))
  fisher_enrich_one <- function(fg_genes, set_genes, bg_genes) {
    N <- length(bg_genes)
    M <- length(intersect(set_genes, bg_genes))
    n <- length(intersect(fg_genes, bg_genes))
    k <- length(intersect(fg_genes, set_genes))
    a <- k; b <- n - k; c <- M - k; d <- N - M - b
    ft <- tryCatch(fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater"),
                   error = function(e) list(p.value = 1, estimate = NA_real_))
    OR <- suppressWarnings(unname(ft$estimate))
    tibble(k = k, n = n, M = M, N = N,
           GeneRatio = ifelse(n > 0, k / n, NA_real_),
           bg_rate = ifelse(N > 0, M / N, NA_real_),
           OR = OR,
           log2OR = ifelse(is.finite(OR) & OR > 0, log2(OR), NA_real_),
           pvalue = ft$p.value,
           hit_genes = paste(sort(intersect(fg_genes, set_genes)), collapse = "/"))
  }

  cog_enrich <- bind_rows(lapply(COG_CLUSTER_ORDER, function(cid) {
    fg <- cluster_assignments %>% filter(cluster == cid) %>% pull(`Gene name`) %>% as.character() %>% intersect(bg_genes)
    bind_rows(lapply(cate_order_df$COG_Category, function(cate) {
      cate_genes <- gene_cate_long %>% filter(COG_Category == cate) %>% pull(Gene) %>% intersect(bg_genes)
      fisher_enrich_one(fg, cate_genes, bg_genes) %>%
        mutate(cluster = cid, COG_Category = cate, .before = 1)
    }))
  })) %>%
    group_by(cluster) %>%
    mutate(qvalue = p.adjust(pvalue, method = "BH")) %>%
    ungroup() %>%
    left_join(cate_order_df, by = "COG_Category") %>%
    mutate(is_sig = !is.na(qvalue) & qvalue < Q_CUTOFF,
           neg_log10q = ifelse(qvalue > 0, -log10(qvalue), NA_real_))

  fwrite(cog_enrich, file.path(DIR_COG, "COG_enrichment_by_cluster_all.tsv"), sep = "\t")
  fwrite(cog_enrich %>% filter(is_sig) %>% arrange(qvalue),
         file.path(DIR_COG, "COG_enrichment_by_cluster_sig.tsv"), sep = "\t")

  label_order <- cate_order_df %>%
    arrange(Order) %>%
    mutate(Label = paste0(COG_Category, " - ", COG_Desc)) %>%
    pull(Label)

  cog_plot_df <- cog_enrich %>%
    mutate(
      Label = factor(paste0(COG_Category, " - ", COG_Desc), levels = rev(label_order)),
      cluster = factor(cluster, levels = COG_CLUSTER_ORDER),
      neg_log10q_plot = ifelse(is_sig, neg_log10q, NA_real_),
      log2OR_plot = ifelse(is_sig, log2OR, NA_real_)
    )

  max_nlq <- max(cog_plot_df$neg_log10q_plot, na.rm = TRUE)
  if (!is.finite(max_nlq) || max_nlq <= 0) max_nlq <- 1
  abs_lor <- max(abs(cog_plot_df$log2OR_plot), na.rm = TRUE)
  if (!is.finite(abs_lor) || abs_lor <= 0) abs_lor <- 1
  plot_h <- max(6, 0.30 * nrow(cate_order_df) + 1.6)
  plot_w <- 8.5

  p_cog_dot <- ggplot() +
    geom_point(data = cog_plot_df %>% filter(!is_sig),
               aes(x = cluster, y = Label, size = GeneRatio),
               color = NS_PT_COL, alpha = 0.65) +
    geom_point(data = cog_plot_df %>% filter(is_sig),
               aes(x = cluster, y = Label, size = GeneRatio, color = neg_log10q_plot),
               alpha = 0.95) +
    scale_size_continuous(name = "GeneRatio (k/n)", range = c(1, 8), labels = percent_format(accuracy = 1)) +
    scale_color_gradientn(name = expression(-log[10]~"(q)"), colours = COLORS_DIV,
                          limits = c(0, max_nlq), oob = squish) +
    labs(title = "COG enrichment dotplot by CCR cluster",
         subtitle = sprintf("Background: corrected RNA-detected genes with COG annotation (N=%d); FDR<%.2f",
                            length(bg_genes), Q_CUTOFF),
         x = NULL, y = NULL) +
    theme_enrichment() +
    theme(axis.text.y = element_text(size = 8))
  save_plot(file.path(DIR_COG, "COG_dotplot_by_cluster"), p_cog_dot, width = plot_w, height = plot_h)

  p_cog_heat <- ggplot() +
    geom_tile(data = cog_plot_df %>% filter(!is_sig),
              aes(x = cluster, y = Label), fill = NS_TILE_COL, color = "white", linewidth = 0.3) +
    geom_tile(data = cog_plot_df %>% filter(is_sig),
              aes(x = cluster, y = Label, fill = log2OR_plot), color = "white", linewidth = 0.3) +
    scale_fill_gradientn(name = "log2(OR)", colours = COLORS_DIV,
                         limits = c(-abs_lor, abs_lor), oob = squish, na.value = NS_TILE_COL) +
    labs(title = "COG enrichment heatmap by CCR cluster",
         subtitle = sprintf("Background: corrected RNA-detected genes with COG annotation (N=%d); FDR<%.2f",
                            length(bg_genes), Q_CUTOFF),
         x = NULL, y = NULL) +
    theme_enrichment() +
    theme(axis.text.y = element_text(size = 8))
  save_plot(file.path(DIR_COG, "COG_heatmap_by_cluster"), p_cog_heat, width = plot_w, height = plot_h)
}

## ================================================================
## 7. TF enrichment by cluster
## ================================================================

message("[06] Running TF enrichment by cluster")
required_tf <- c(PATH_ECO_TXT, PATH_ECO_TF)
if (!all(file.exists(required_tf))) {
  warning("TF annotation files not found. Skipping TF enrichment.")
} else {
  eco_raw <- read_tsv(PATH_ECO_TXT, col_types = cols(.default = "c"), locale = locale(encoding = "UTF-8"))
  colnames(eco_raw) <- colnames(eco_raw) %>%
    str_replace_all('"', '') %>% str_replace_all("\r", "") %>% str_trim()
  eco_raw <- eco_raw %>% mutate(Gene = str_trim(str_replace_all(Gene, "\r", "")))
  TF_COL <- "Regulation - transcription factor regulators of gene"
  if (!TF_COL %in% colnames(eco_raw)) stop("EcoCyc table is missing TF regulator column.")

  tf_csv <- read.csv(PATH_ECO_TF, fileEncoding = "UTF-8-BOM", check.names = FALSE)
  if (!all(c("tf_raw", "tf_short_manual", "genes_all_txt") %in% colnames(tf_csv))) {
    stop("TF mapping table must contain tf_raw, tf_short_manual, and genes_all_txt.")
  }

  gene2tf <- tf_csv %>%
    select(tf = tf_short_manual, genes_all_txt) %>%
    mutate(gene_list = str_split(genes_all_txt, "/")) %>%
    unnest(gene_list) %>%
    transmute(gene = str_trim(gene_list), tf = str_trim(tf)) %>%
    filter(gene != "", tf != "", !is.na(gene), !is.na(tf)) %>%
    distinct(gene, tf)

  rna_genes <- bg_df$`Gene name` %>% as.character() %>% unique()
  eco_genes_all <- eco_raw$Gene %>% na.omit() %>% unique()
  U1 <- intersect(rna_genes, eco_genes_all)
  eco_genes_tfannot <- eco_raw %>%
    filter(!is.na(.data[[TF_COL]]),
           str_trim(as.character(.data[[TF_COL]])) != "",
           str_trim(as.character(.data[[TF_COL]])) != "NA") %>%
    pull(Gene) %>% na.omit() %>% unique()
  U2 <- intersect(rna_genes, eco_genes_tfannot)

  tf_sets_full <- gene2tf %>%
    group_by(tf) %>%
    summarise(genes = list(sort(unique(gene))), M_all = n_distinct(gene), .groups = "drop")

  universe_list <- switch(
    TF_UNIVERSE,
    "U1" = list(U1 = U1),
    "U2" = list(U2 = U2),
    "both" = list(U1 = U1, U2 = U2),
    stop("TF_UNIVERSE must be 'U1', 'U2', or 'both'.")
  )

  enrich_one_tf <- function(fg_genes, tf_genes, bg_genes) {
    N <- length(bg_genes)
    M <- length(intersect(tf_genes, bg_genes))
    n <- length(intersect(fg_genes, bg_genes))
    k <- length(intersect(fg_genes, tf_genes))
    pval <- if (n > 0 && M > 0) phyper(k - 1, M, N - M, n, lower.tail = FALSE) else 1
    GeneRatio <- ifelse(n > 0, k / n, NA_real_)
    bg_rate <- ifelse(N > 0, M / N, NA_real_)
    FE <- ifelse(!is.na(bg_rate) & bg_rate > 0, GeneRatio / bg_rate, NA_real_)
    tibble(k = k, n = n, M = M, N = N,
           GeneRatio = GeneRatio,
           bg_rate = bg_rate,
           FE = FE,
           log2FE = ifelse(!is.na(FE) & FE > 0, log2(FE), NA_real_),
           pvalue = pval,
           geneID = paste(sort(intersect(fg_genes, tf_genes)), collapse = "/"))
  }

  for (u_name in names(universe_list)) {
    universe_genes <- universe_list[[u_name]]
    tf_outdir <- file.path(DIR_TF, u_name, "by_cluster")
    dir.create(tf_outdir, recursive = TRUE, showWarnings = FALSE)

    tf_sets <- tf_sets_full %>%
      mutate(genes_in_U = map(genes, ~ intersect(.x, universe_genes)),
             M = map_int(genes_in_U, length)) %>%
      filter(M >= M_MIN_TF)

    tf_enrich <- bind_rows(lapply(TF_CLUSTER_ORDER, function(cid) {
      fg <- cluster_assignments %>% filter(cluster == cid) %>% pull(`Gene name`) %>% as.character() %>% intersect(universe_genes)
      bind_rows(lapply(seq_len(nrow(tf_sets)), function(i) {
        enrich_one_tf(fg, tf_sets$genes_in_U[[i]], universe_genes) %>%
          mutate(cluster = cid, tf = tf_sets$tf[i], .before = 1)
      }))
    })) %>%
      group_by(cluster) %>%
      mutate(qvalue = p.adjust(pvalue, method = "BH")) %>%
      ungroup() %>%
      mutate(is_sig = !is.na(qvalue) & qvalue < Q_CUTOFF,
             neg_log10q = ifelse(qvalue > 0, -log10(qvalue), NA_real_),
             universe = u_name)

    fwrite(tf_enrich, file.path(tf_outdir, "TF_enrich_by_cluster_all.tsv"), sep = "\t")
    fwrite(tf_enrich %>% filter(is_sig) %>% arrange(qvalue, desc(GeneRatio)),
           file.path(tf_outdir, "TF_enrich_by_cluster_sig.tsv"), sep = "\t")

    main_sel <- tf_enrich %>%
      filter(is_sig) %>%
      group_by(cluster) %>%
      arrange(qvalue, desc(GeneRatio), .by_group = TRUE) %>%
      slice_head(n = TOP_N_PER_CLUSTER) %>%
      ungroup()
    sig_tbl <- tf_enrich %>%
      filter(is_sig) %>%
      arrange(factor(cluster, levels = TF_CLUSTER_ORDER), qvalue, desc(GeneRatio))
    fwrite(main_sel, file.path(tf_outdir, "TF_main_selected.tsv"), sep = "\t")

    plot_tf_set <- function(selected_tbl, dot_base, heatmap_base, title_prefix, subtitle_text) {
      if (nrow(selected_tbl) == 0) return(invisible(NULL))

      tf_levels <- selected_tbl %>%
        arrange(factor(cluster, levels = TF_CLUSTER_ORDER), qvalue, desc(GeneRatio)) %>%
        pull(tf) %>%
        unique()

      plot_df <- tf_enrich %>%
        filter(tf %in% tf_levels) %>%
        mutate(
          tf = factor(tf, levels = tf_levels),
          cluster = factor(cluster,  levels = rev(TF_CLUSTER_ORDER)),
          neg_log10q_plot = ifelse(is_sig, neg_log10q, NA_real_),
          log2FE_plot = ifelse(is_sig, log2FE, NA_real_)
        )

      max_nlq <- max(plot_df$neg_log10q_plot, na.rm = TRUE)
      if (!is.finite(max_nlq) || max_nlq <= 0) max_nlq <- 1
      abs_lfe <- max(abs(plot_df$log2FE_plot), na.rm = TRUE)
      if (!is.finite(abs_lfe) || abs_lfe <= 0) abs_lfe <- 1
      plot_w <- max(8, 0.35 * length(tf_levels) + 2)

      p_tf_dot <- ggplot() +
        geom_point(data = plot_df %>% filter(!is_sig),
                   aes(x = tf, y = cluster, size = GeneRatio),
                   color = NS_PT_COL, alpha = 0.65) +
        geom_point(data = plot_df %>% filter(is_sig),
                   aes(x = tf, y = cluster, size = GeneRatio, color = neg_log10q_plot),
                   alpha = 0.95) +
        scale_size_continuous(name = "GeneRatio (k/n)", range = c(1, 8), labels = percent_format(accuracy = 1)) +
        scale_color_gradientn(name = expression(-log[10]~"(q)"), colours = COLORS_DIV,
                              limits = c(0, max_nlq), oob = squish) +
        labs(title = paste0(title_prefix, " dotplot by CCR cluster (", u_name, ")"),
             subtitle = subtitle_text,
             x = NULL, y = NULL) +
        theme_enrichment() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      save_plot(file.path(tf_outdir, dot_base), p_tf_dot, width = plot_w, height = 2.4)

      p_tf_heat <- ggplot() +
        geom_tile(data = plot_df %>% filter(!is_sig),
                  aes(x = tf, y = cluster), fill = NS_TILE_COL, color = "white", linewidth = 0.3) +
        geom_tile(data = plot_df %>% filter(is_sig),
                  aes(x = tf, y = cluster, fill = log2FE_plot), color = "white", linewidth = 0.3) +
        scale_fill_gradientn(name = "log2(FE)", colours = COLORS_DIV,
                             limits = c(-abs_lfe, abs_lfe), oob = squish, na.value = NS_TILE_COL) +
        labs(title = paste0(title_prefix, " heatmap by CCR cluster (", u_name, ")"),
             subtitle = subtitle_text,
             x = NULL, y = NULL) +
        theme_enrichment() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      save_plot(file.path(tf_outdir, heatmap_base), p_tf_heat, width = plot_w, height = 4.8)
    }

    plot_tf_set(
      main_sel,
      "TF_dotplot_main",
      "TF_heatmap_main",
      "TF enrichment",
      sprintf("M_min=%d; top %d significant TFs per cluster; FDR<%.2f",
              M_MIN_TF, TOP_N_PER_CLUSTER, Q_CUTOFF)
    )

    plot_tf_set(
      sig_tbl,
      "SI_TF_dotplot_sig",
      "SI_TF_heatmap_sig",
      "SI TF enrichment",
      sprintf("M_min=%d; all significant TFs; FDR<%.2f",
              M_MIN_TF, Q_CUTOFF)
    )
  }
}

## ================================================================
# 3. Export dcw component profiles
out<-file.path(cfg$results_dir,"baby_machine")
calcDir<-file.path(out,"BabyMachineDataAnalysis_output")
corr<-read.csv(file.path(calcDir,"corrected_TPM_total_M_Ala.csv"),check.names=FALSE)
sim<-read.csv(file.path(calcDir,"corrected_TPM_Simulation.csv"),check.names=FALSE)
needgenes<-cfg$dcw_main

times<-as.character(seq(0,220,11))
keep<-corr[["Gene name"]]%in%needgenes
trace<-list()
for(i in which(keep)) {
  j<-match(corr$ID[i],sim$ID)
  y<-as.numeric(corr[i,times]);f<-as.numeric(sim[j,times])
  trace[[length(trace)+1L]]<-data.frame(ID=corr$ID[i],gene=corr[["Gene name"]][i],
    time_min=as.numeric(times),corrected_TPM=y,fitted_TPM=f,
    relative_corrected_TPM=y/mean(y),relative_fitted_TPM=f/mean(y))
}
z<-do.call(rbind,trace)

# Correction and fitting precede display-only aggregation.
dcw<-z[z$gene%in%cfg$dcw_main,,drop=FALSE]
write_csv(dcw,file.path(out,"ccr","tables","dcw_component_profiles.csv"))
if(!setequal(unique(dcw$gene),cfg$dcw_main)) stop("Missing dcw components")
if(anyDuplicated(dcw[c("gene","time_min")]) ||
   any(!is.finite(as.matrix(dcw[c("corrected_TPM","fitted_TPM")]))))
  stop("Invalid or duplicated dcw component values")
message("CCR clustering and enrichment completed. Output: ", OUTDIR)
