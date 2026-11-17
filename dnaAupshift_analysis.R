# dnaA upshift analysis
#
# Purpose: Recompute RNA response, clustering and functional enrichment.
# Input: Table S4 sheets dnaA_RNA and dnaA_DNA; Analysis_annotation_files.
# Output: results/dnaA_upshift: filtered and corrected matrices, statistics, cluster memberships, enrichment and figures.
# Scope: processed measurements reported in Table S4. 
## Cluster labels: For presentation in the main-text figures and Data S2,
# script clusters 2, 1, 4 and 3 were relabeled as clusters 1, 2, 3 and 4, respectively. 
# This relabeling does not change cluster membership.This script retains the original cluster labels in its outputs.


# 1. Input paths and shared analysis settings
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
  if(!nzchar(override)) override<-Sys.getenv("DNAA_SCRIPT_DIR",unset="")
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
# 2. In-memory workbook adapter and step-specific analysis functions
# Convert Table S4 measurements into the analysis matrix and sample metadata.
prepare_dnaa <- function(cfg) {
  rna<-read_sheet(cfg$table_s4,"dnaA_RNA")
  dna<-read_sheet(cfg$table_s4,"dnaA_DNA")
  rs<-sample_design(names(rna));ds<-sample_design(names(dna))
  groups<-c("Replicate1","Replicate2","Replicate3")
  rs$Group<-paste0(groups[rs$replicate],"RNAseq")
  rs$Colnames<-sprintf("RNA_R%d_T%02d_TPM",rs$replicate,ave(seq_len(nrow(rs)),rs$replicate,FUN=seq_along))
  ds$Group<-paste0(groups[ds$replicate],"DNA")
  ds$Colnames<-sprintf("DNA_R%d_T%02d",ds$replicate,ave(seq_len(nrow(ds)),ds$replicate,FUN=seq_along))
  metadata<-function(x)data.frame(Colnames=x$Colnames,Group=x$Group,
    Treatment=ifelse(x$time_min<0,"Pre_1aTc","Post_10aTc"),
    Time=x$time_min,Time_h=x$time_min/60,Real_Uchanged_Time=x$time_min)
  rr<-rna[c("gene","locus_tag","product","Rel.Pos",rs$sample)]
  names(rr)[5:ncol(rr)]<-rs$Colnames
  dd<-data.frame(genome_position=dna$Relative_position,
    Relative_position=dna$Relative_position,dna[ds$sample],check.names=FALSE)
  names(dd)[3:ncol(dd)]<-ds$Colnames
  inputs<-list(rna_TPM.csv=rr,rna_metadata.csv=metadata(rs),
    dna_DPM.csv=dd,dna_metadata.csv=metadata(ds))
  # Hold the workbook adapters as CSV text in memory. This preserves the
  # original readers' column names, numeric precision and missing-value rules
  # without writing intermediate input files into the output directory.
  for(n in names(inputs)) {
    buffer<-character();con<-textConnection("buffer","w",local=TRUE)
    utils::write.csv(inputs[[n]],con,row.names=FALSE,na="")
    close(con)
    assign(n,paste(buffer,collapse="\n"),envir=analysis_inputs)
  }
  invisible(list(rna_samples=rs,dna_samples=ds))
}

analysis_inputs<-new.env(parent=emptyenv())
analysis_result_paths<-c(
  filtered_genes.csv="Step2_RNA_filtering/genes_C2_step2_after_preFluct.csv",
  corrected_TPM.csv="Step4a_DNA_dosage_correction/corrTPM_ALL_groups_merged.csv",
  response_all_stats.csv="Step5_response_identification/tables/all_genes_stats.csv",
  responsive_genes.csv="Step5_response_identification/tables/responsive_genes.csv",
  cluster_assignments.csv="Step6_clustering/tables/cluster_assignments.csv")
analysis_path<-function(file) {
  if(length(file)!=1L || !is.character(file)) stop("Expected one input path")
  if(file %in% names(analysis_result_paths)) return(unname(analysis_result_paths[[file]]))
  candidate<-file.path(cfg$annotation_dir,file)
  if(!grepl("[/\\\\]",file) && file.exists(candidate)) return(candidate)
  file
}
analysis_file_exists<-function(...) {
  paths<-c(...)
  vapply(paths,function(p)exists(p,envir=analysis_inputs,inherits=FALSE)||file.exists(analysis_path(p)),logical(1),USE.NAMES=FALSE)
}
analysis_read_csv<-function(file,...) {
  if(exists(file,envir=analysis_inputs,inherits=FALSE))
    return(utils::read.csv(text=get(file,envir=analysis_inputs),...))
  utils::read.csv(analysis_path(file),...)
}
analysis_read_tibble<-function(file,...) {
  if(exists(file,envir=analysis_inputs,inherits=FALSE))
    return(readr::read_csv(I(get(file,envir=analysis_inputs)),...))
  readr::read_csv(analysis_path(file),...)
}
analysis_fread<-function(file,...) {
  if(exists(file,envir=analysis_inputs,inherits=FALSE))
    return(data.table::fread(text=get(file,envir=analysis_inputs),...))
  data.table::fread(analysis_path(file),...)
}
analysis_read_tsv<-function(file,...)readr::read_tsv(analysis_path(file),...)

# Analysis section: DNA_copy_number_heatmaps
DNA_copy_number_heatmaps <- function() {


df.DNA_DPM <- analysis_read_csv("dna_DPM.csv", header = TRUE, sep = ",")
df.DNA_ExpInfo <- analysis_read_csv("dna_metadata.csv", header = TRUE, sep = ",")

# ============================================================
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(scales)
})

# ============================================================
# [USER SETTINGS]
# ============================================================

DATA_CSV <- "dna_DPM.csv"
INFO_CSV <- "dna_metadata.csv"

COL_SAMPLEID <- "Colnames"
COL_TIME     <- "Time_h"       # ← 改为小时列
COL_GROUP    <- "Group"
COL_TREAT    <- "Treatment"
PRE_LABEL    <- "Pre_1aTc"

#
#
#
CALC_TIME_RANGE <- NULL          # 推荐 NULL，全部参与计算
PLOT_TIME_RANGE <- NULL #c(-0.5, 3.7)  # 热图显示范围（h），NULL = 全显示

VALUE_TYPE <- "DPM"

#   "log2"   : log2(value + PSEUDO)
SCALE_MODE <- "zscore"

PSEUDO <- 0

#   ZRANGE_MODE:
#
#
ZRANGE_MODE   <- "shared"
ZRANGE_SHARED <- NULL

PALETTE <- c("#3B4CC0", "#88A1D2", "#F7F7F7", "#F3A982", "#B40426")

N_XLABELS <- 8

OUT_DIR   <- "Step1_DNA_heatmaps"
HEAT_W_MM <- 150
HEAT_H_MM <- 130
TOPBAR_H  <- 0.04

SAMPLE_TICK_LEN <- 10       # 短线长度
SAMPLE_TICK_LWD <- 0.4      # 线宽
SAMPLE_TICK_COL <- "black"  # 颜色

# ============================================================
# ============================================================
mod_load <- function(data_csv, info_csv,
                     col_sampleid, col_time,
                     col_group, col_treat) {
  mat_raw  <- analysis_read_csv(data_csv, header = TRUE, check.names = FALSE)
  info_raw <- analysis_read_csv(info_csv, header = TRUE, check.names = FALSE)

  coord <- mat_raw[, 1:2]
  colnames(coord) <- c("genome_position", "Relative_position")  # ← 此处对调

  valid_cols <- info_raw[[col_sampleid]][
    info_raw[[col_sampleid]] %in% colnames(mat_raw)
  ]

  mat <- as.matrix(mat_raw[, valid_cols, drop = FALSE])
  rownames(mat) <- coord$genome_position

  info <- info_raw %>%
  mutate(
    Time  = suppressWarnings(as.numeric(.data[[col_time]])),
    Group = as.character(.data[[col_group]]),
    Treat = as.character(.data[[col_treat]])
  )
  list(mat = mat, info = info, coord = coord)
}

# ============================================================
# ============================================================
mod_rel_dpm_group <- function(mat_group, info_group,
                               col_sampleid, pre_label) {
  pre_cols <- info_group[[col_sampleid]][
    info_group$Treat == pre_label
  ]

  if (length(pre_cols) == 0) {
    message("  [注意] 本组无 pre 样品，使用原始 DPM。")
    return(list(mat = mat_group, value_type_used = "DPM"))
  }

  pre_mean <- rowMeans(mat_group[, pre_cols, drop = FALSE],
                       na.rm = TRUE)
  pre_mean[pre_mean == 0 | !is.finite(pre_mean)] <- NA

  mat_rel <- sweep(mat_group, 1, pre_mean, "/")
  mat_rel[!is.finite(mat_rel)] <- NA

  list(mat = mat_rel, value_type_used = "rel.DPM")
}

# ============================================================
# ============================================================
mod_scale <- function(mat, mode = "log2", pseudo = 0) {
  if (mode == "log2") {
    mat <- log2(mat + pseudo)
    mat[!is.finite(mat)] <- NA
  } else if (mode == "zscore") {
    mu  <- rowMeans(mat, na.rm = TRUE)
    sdv <- apply(mat, 1, sd, na.rm = TRUE)
    sdv[sdv == 0 | !is.finite(sdv)] <- 1
    mat <- sweep(mat, 1, mu,  "-")
    mat <- sweep(mat, 1, sdv, "/")
    mat[!is.finite(mat)] <- 0
  }
  mat
}

# ============================================================
# ============================================================
mod_zrange <- function(mat, value_type_used, scale_mode) {
  vals <- mat[is.finite(mat)]
  if (length(vals) == 0) return(c(-3, 3))

  if (scale_mode == "zscore") return(c(-3, 3))

  if (value_type_used == "rel.DPM") {
    lim <- quantile(abs(vals), 0.98, na.rm = TRUE)
    lim <- ceiling(lim * 4) / 4
    return(c(-lim, lim))
  }

  lo <- floor(quantile(vals, 0.01, na.rm = TRUE) * 4) / 4
  hi <- ceiling(quantile(vals, 0.99, na.rm = TRUE) * 4) / 4
  c(lo, hi)
}

# ============================================================
# ============================================================
mod_xspec <- function(time_vec) {
  x <- sort(unique(as.numeric(time_vec)))
  x <- x[is.finite(x)]
  if (length(x) == 1)
    return(tibble(Time = x, xmin = x - 0.5, xmax = x + 0.5))
  mid  <- (x[-1] + x[-length(x)]) / 2
  xmin <- c(x[1] - (mid[1] - x[1]), mid)
  xmax <- c(mid, x[length(x)] + (x[length(x)] - mid[length(mid)]))
  tibble(Time = x, xmin = xmin, xmax = xmax)
}

# ============================================================
# ============================================================
mod_yspec <- function(coord) {
  coord_ord <- coord %>% arrange(desc(Relative_position))
  n <- nrow(coord_ord)
  tibble(
    genome_position   = as.character(coord_ord$genome_position),
    Relative_position = coord_ord$Relative_position,
    y    = seq(n, 1),
    ymin = seq(n, 1) - 0.5,
    ymax = seq(n, 1) + 0.5
  )
}

# ============================================================
# ============================================================
mod_plotdata <- function(mat_scaled, info_plot,
                          col_sampleid, xspec, yspec) {
  cols    <- info_plot[[col_sampleid]]
  mat_sub <- mat_scaled[, cols, drop = FALSE]

  as.data.frame(mat_sub) %>%
    rownames_to_column("genome_position") %>%
    mutate(genome_position = as.character(genome_position)) %>%
    pivot_longer(-genome_position,
                 names_to  = col_sampleid,
                 values_to = "value") %>%
    left_join(
      info_plot %>% select(all_of(c(col_sampleid, "Time"))),
      by = col_sampleid
    ) %>%
    left_join(xspec, by = "Time") %>%
    left_join(
      yspec %>%
        mutate(genome_position = as.character(genome_position)) %>%
        select(genome_position, y, ymin, ymax),
      by = "genome_position"
    ) %>%
    filter(is.finite(xmin), is.finite(ymin), is.finite(value))
}

# ============================================================
# ============================================================
mod_plot_single <- function(plotdata, xspec, yspec,
                             z_range, palette, n_xlabels,
                             title       = NULL,
                             value_label = "value",
                             topbar_h    = 0.04,
                             sample_tick_len  = 10,    # 短线长度（行数单位，向上）
                             sample_tick_lwd  = 0.4,
                             sample_tick_col  = "black") {

  x_lim <- range(c(xspec$xmin, xspec$xmax), na.rm = TRUE)
  y_lim <- range(c(yspec$ymin, yspec$ymax), na.rm = TRUE)

  xb <- breaks_pretty(n = n_xlabels)(range(xspec$Time))
  xb <- xb[is.finite(xb) & xb >= x_lim[1] & xb <= x_lim[2]]

  sample_times <- sort(unique(xspec$Time))
  tick_df <- tibble(
    x    = sample_times,
    xend = sample_times,
    y    = y_lim[1],
    yend = y_lim[1] + sample_tick_len   # 向上延伸
  )

  theme_heat <- theme_bw(base_size = 9) +
    theme(
      panel.grid        = element_blank(),
      panel.border      = element_rect(fill = NA, colour = "black",
                                       linewidth = 0.4),
      axis.ticks        = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = unit(1.0, "mm"),
      axis.text         = element_text(colour = "black"),
      plot.title        = element_text(size = 9, hjust = 0.5,
                                       face = "bold"),
      plot.margin       = margin(0, 4, 4, 4)
    )

  p_top <- ggplot(xspec) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = 0, ymax = 1, fill = Time),
              colour = NA) +
    scale_fill_gradientn(colours = palette, guide = "none") +
    scale_x_continuous(limits = x_lim, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(
      panel.border = element_rect(fill = NA, colour = "black",
                                  linewidth = 0.4),
      plot.margin  = margin(4, 4, 0, 4)
    )

target_rp  <- c(-1, -0.5, 0, 0.5, 1)
yspec_lab  <- lapply(target_rp, function(rp) {
  yspec %>%
    mutate(dist = abs(Relative_position - rp)) %>%
    slice_min(dist, n = 1, with_ties = FALSE) %>%
    mutate(target = rp)
}) %>%
  bind_rows() %>%
  distinct(y, .keep_all = TRUE) %>%    # ← 改为按 y 去重
  arrange(desc(y))                     # ← 保证与 y 轴方向一致

  p_heat <- ggplot(plotdata) +

    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = value), colour = NA) +

    geom_segment(
      data        = tick_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      colour      = sample_tick_col,
      linewidth   = sample_tick_lwd,
      inherit.aes = FALSE
    ) +

    scale_fill_gradientn(
      colours = palette,
      limits  = z_range,
      oob     = scales::squish,
      name    = value_label,
      guide   = guide_colorbar(
        title.position = "top",
        title.hjust    = 0.5,
        barheight      = unit(40, "mm"),
        barwidth       = unit(4,  "mm"),
        frame.colour   = "black",
        ticks.colour   = "black"
      )
    ) +

    scale_x_continuous(
      breaks = xb,
      labels = round(xb, 2),
      limits = x_lim,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
  limits = y_lim,
  expand = c(0, 0),
  breaks = yspec_lab$y,
  labels = yspec_lab$target    # ← 用精确的目标值，而非实际Rel.Pos
) +
    labs(title = title,
         x     = "Time (h)",
         y     = "Relative position") +
    theme_heat +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.title    = element_text(size = 8),
      legend.text     = element_text(size = 7),
      legend.margin   = margin(2, 2, 2, 2)
    )

  cowplot::plot_grid(
    p_top, p_heat,
    ncol        = 1,
    rel_heights = c(topbar_h, 1),
    align       = "v",
    axis        = "lr"
  )
}

# ============================================================
# MAIN
# ============================================================
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat   <- mod_load(DATA_CSV, INFO_CSV,
                  COL_SAMPLEID, COL_TIME, COL_GROUP, COL_TREAT)
mat   <- dat$mat
info  <- dat$info
coord <- dat$coord

stopifnot(
  range(coord$Relative_position) == c(-1, 1) |
  (min(coord$Relative_position) >= -1 & max(coord$Relative_position) <= 1)
)
message("✓ Relative_position 范围正常: [",
        round(min(coord$Relative_position), 3), ", ",
        round(max(coord$Relative_position), 3), "]")

yspec  <- mod_yspec(coord)
groups <- unique(info[[COL_GROUP]])
message("检测到分组：", paste(groups, collapse = ", "))

group_results <- list()

for (g in groups) {
  message("\n── 处理组：", g, " ──")

  g_info_all <- info %>%
    filter(Group == g) %>%
    arrange(Time)

  if (!is.null(CALC_TIME_RANGE)) {
    g_info_calc <- g_info_all %>%
      filter(Time >= CALC_TIME_RANGE[1],
             Time <= CALC_TIME_RANGE[2])
    message("  计算范围：", CALC_TIME_RANGE[1], " ~ ",
            CALC_TIME_RANGE[2], " h，共 ",
            nrow(g_info_calc), " 个样品")
  } else {
    g_info_calc <- g_info_all
    message("  计算范围：全部 ", nrow(g_info_calc), " 个样品")
  }

  if (!is.null(PLOT_TIME_RANGE)) {
    g_info_plot <- g_info_calc %>%
      filter(Time >= PLOT_TIME_RANGE[1],
             Time <= PLOT_TIME_RANGE[2])
    message("  显示范围：", PLOT_TIME_RANGE[1], " ~ ",
            PLOT_TIME_RANGE[2], " h，共 ",
            nrow(g_info_plot), " 个样品")
  } else {
    g_info_plot <- g_info_calc
    message("  显示范围：全部计算样品")
  }

  if (nrow(g_info_plot) == 0) {
    warning("组 [", g, "] 在指定显示范围内无样品，跳过。")
    next
  }

  calc_cols <- g_info_calc[[COL_SAMPLEID]]
  mat_calc  <- mat[, calc_cols, drop = FALSE]

  if (VALUE_TYPE == "rel.DPM") {
    rel_res    <- mod_rel_dpm_group(mat_calc, g_info_calc,
                                    COL_SAMPLEID, PRE_LABEL)
    mat_val    <- rel_res$mat
    vtype_used <- rel_res$value_type_used
  } else {
    mat_val    <- mat_calc
    vtype_used <- "DPM"
  }

  pseudo_use <- ifelse(vtype_used == "DPM" && SCALE_MODE == "log2",
                       max(PSEUDO, 1), PSEUDO)
  mat_scaled <- mod_scale(mat_val, mode = SCALE_MODE,
                           pseudo = pseudo_use)

  plot_cols  <- g_info_plot[[COL_SAMPLEID]]
  mat_plot   <- mat_scaled[, plot_cols, drop = FALSE]

  xspec_g <- mod_xspec(g_info_plot[["Time"]])

  pd <- mod_plotdata(mat_plot, g_info_plot,
                     COL_SAMPLEID, xspec_g, yspec)

  val_label <- switch(
    paste(vtype_used, SCALE_MODE, sep = "_"),
    "DPM_log2"       = "log2(DPM+1)",
    "DPM_zscore"     = "z-score(DPM)",
    "rel.DPM_log2"   = "log2(rel.DPM)",
    "rel.DPM_zscore" = "z-score(rel.DPM)"
  )

  group_results[[g]] <- list(
    mat_scaled = mat_scaled,   # 全计算集合，用于共享 Z_RANGE 推断
    xspec      = xspec_g,
    pd         = pd,
    vtype_used = vtype_used,
    val_label  = val_label,
    n_plot     = nrow(g_info_plot)
  )
}

if (ZRANGE_MODE == "shared") {
  if (!is.null(ZRANGE_SHARED)) {
    zrange_shared <- ZRANGE_SHARED
    message("\n共享 Z_RANGE（手动）: [",
            zrange_shared[1], ", ", zrange_shared[2], "]")
  } else {
    match_groups <- names(group_results)[
      sapply(group_results, function(r) r$vtype_used == VALUE_TYPE)
    ]
    if (length(match_groups) == 0) match_groups <- names(group_results)
    message("\n共享 Z_RANGE 参与推断的组：",
            paste(match_groups, collapse = ", "))

    all_vals <- unlist(lapply(group_results[match_groups], function(r)
      r$mat_scaled[is.finite(r$mat_scaled)]))
    vtype_ref     <- group_results[[match_groups[1]]]$vtype_used
    zrange_shared <- mod_zrange(matrix(all_vals), vtype_ref, SCALE_MODE)
    message("共享 Z_RANGE（自动推断）: [",
            round(zrange_shared[1], 3), ", ",
            round(zrange_shared[2], 3), "]")
  }
}

for (g in names(group_results)) {
  r <- group_results[[g]]

  if (ZRANGE_MODE == "shared" && r$vtype_used == VALUE_TYPE) {
    zr <- zrange_shared
  } else {
    zr <- mod_zrange(r$mat_scaled, r$vtype_used, SCALE_MODE)
    message("组 [", g, "] Z_RANGE（独立推断）: [",
            round(zr[1], 3), ", ", round(zr[2], 3), "]")
  }

  time_tag <- if (!is.null(PLOT_TIME_RANGE))
    paste0("_t", PLOT_TIME_RANGE[1], "to", PLOT_TIME_RANGE[2], "h")
  else "_allTime"

  out_file <- file.path(
    OUT_DIR,
    paste0(g, "_", r$vtype_used, "_", SCALE_MODE,
           time_tag, ".pdf")
  )

  p <- mod_plot_single(
    plotdata         = r$pd,
    xspec            = r$xspec,
    yspec            = yspec,
    z_range          = zr,
    palette          = PALETTE,
    n_xlabels        = N_XLABELS,
    title            = g,
    value_label      = r$val_label,
    topbar_h         = TOPBAR_H,
    sample_tick_len  = SAMPLE_TICK_LEN,
    sample_tick_lwd  = SAMPLE_TICK_LWD,
    sample_tick_col  = SAMPLE_TICK_COL
  )

  ggsave(out_file, p,
         width       = HEAT_W_MM,
         height      = HEAT_H_MM,
         units       = "mm",
         useDingbats = FALSE,
         limitsize   = FALSE)

  message("✓ 已保存：", out_file,
          "  (显示 ", r$n_plot, " 个时间点)",
          "  Z_RANGE: [", round(zr[1], 3), ", ",
          round(zr[2], 3), "]")
}

message("\n全部完成，输出目录：",
        normalizePath(OUT_DIR, mustWork = FALSE))
}


# Analysis section: RNA_relative_abundance
RNA_relative_abundance <- function() {


df.RNA_raw <- analysis_read_csv("rna_TPM.csv", header = TRUE, sep = ",")
RNA_Exp <- analysis_read_csv("rna_metadata.csv", header = TRUE, sep = ",")


## ============================================================================
##  Combine per-experiment relative expression into ONE wide table
##  (same layout as the input matrix: ID columns + every sample column,
##   in the input's original column order)
##
##  For each experiment (rows of the sample sheet sharing a `group_col` value):
##      baseline(gene) = mean of the first `n_pre` pre-shift points
##                       (sheet rows with  time_col < pre_time_max)
##      rel(gene, t)   = value(gene, t) / baseline(gene)
##  An experiment whose baseline == 0 for a gene is written as `pre0_label`
##  for that experiment's columns only (so a gene can be numeric in one
##  experiment and pre0 in another).
##
##  REUSE: edit the CONFIG block only. No experiment/sample names are
##  hard-coded; the sample sheet drives everything.
## ============================================================================

## ------------------------------- CONFIG ------------------------------------
cfg <- list(
  data_file    = "rna_TPM.csv",        # wide matrix: ID cols + one column per sample
  sheet_file   = "rna_metadata.csv", # sample sheet: one row per sample
  out_file     = "Step2_RNA_filtering/relative_TPM.csv",                    # output (single wide table)

  name_col     = "Colnames",   # sheet column holding each sample's matrix column name
  group_col    = "Group",      # sheet column that splits samples into experiments
  time_col     = "Time",       # numeric sheet column used to locate pre-shift points

  pre_time_max = 0,            # pre-shift points = sheet rows with  time_col < pre_time_max
  n_pre        = 3,            # number of pre-shift points averaged into the baseline
  digits       = 6,            # significant digits for the written ratios
  pre0_label   = "pre0",       # text written when an experiment's baseline == 0

  id_cols        = NULL,       # ID/annotation columns to carry over.
                               #   NULL = auto-detect (data columns not named as samples)
  sheet_encoding = "UTF-8-BOM" # set to "" if the sample sheet has no byte-order mark
)

## ----------------------------- FUNCTIONS -----------------------------------

## Map a sample's recorded name to the real column in `data`.
## Copes with analysis_read_csv(check.names = TRUE), which prefixes "X" to names that

resolve_col <- function(nm, cols) {
  if (nm %in% cols)  return(nm)
  alt <- make.names(nm)
  if (alt %in% cols) return(alt)
  NA_character_
}

build_relative_table <- function(cfg) {

  data  <- analysis_read_csv(cfg$data_file,  check.names = FALSE, stringsAsFactors = FALSE)
  sheet <- analysis_read_csv(cfg$sheet_file, check.names = FALSE, stringsAsFactors = FALSE,
                    fileEncoding = cfg$sheet_encoding)

  ## resolve every sample to a real data column
  sheet$.col <- vapply(sheet[[cfg$name_col]], resolve_col, character(1),
                       cols = colnames(data))
  if (anyNA(sheet$.col))
    stop("Samples missing from the data matrix: ",
         paste(sheet[[cfg$name_col]][is.na(sheet$.col)], collapse = ", "))

  sample_cols <- sheet$.col
  id_cols <- if (is.null(cfg$id_cols)) setdiff(colnames(data), sample_cols) else cfg$id_cols

  ## character matrix to hold the result (genes x samples)
  rel <- matrix(NA_character_, nrow = nrow(data), ncol = length(sample_cols),
                dimnames = list(NULL, sample_cols))

  for (g in unique(sheet[[cfg$group_col]])) {
    s <- sheet[sheet[[cfg$group_col]] == g, ]
    s <- s[order(s[[cfg$time_col]]), ]

    pre_idx <- which(s[[cfg$time_col]] < cfg$pre_time_max)
    if (length(pre_idx) < cfg$n_pre)
      stop("Experiment '", g, "': fewer than ", cfg$n_pre, " pre-shift points (",
           cfg$time_col, " < ", cfg$pre_time_max, ").")
    pre_cols <- s$.col[head(pre_idx, cfg$n_pre)]

    baseline <- rowMeans(data[, pre_cols, drop = FALSE], na.rm = TRUE)
    block    <- sweep(as.matrix(data[, s$.col, drop = FALSE]), 1, baseline, FUN = "/")

    block_chr <- matrix(formatC(block, format = "g", digits = cfg$digits),
                        nrow = nrow(block))
    zero_base <- which(baseline == 0)
    if (length(zero_base) > 0) block_chr[zero_base, ] <- cfg$pre0_label

    rel[, s$.col] <- block_chr   # place this experiment's columns by name

    message(sprintf("[%s] %d samples | baseline = %s | %d gene(s) -> %s",
                    g, nrow(s), paste(pre_cols, collapse = ", "),
                    length(zero_base), cfg$pre0_label))
  }

  ## restore the input's original sample-column order, prepend the ID columns
  rel <- rel[, intersect(colnames(data), sample_cols), drop = FALSE]
  out <- data.frame(data[, id_cols, drop = FALSE], rel,
                    check.names = FALSE, stringsAsFactors = FALSE)

  write.csv(out, cfg$out_file, row.names = FALSE)
  message(sprintf("\nWrote %d genes x %d samples -> %s",
                  nrow(out), length(sample_cols), cfg$out_file))
  invisible(out)
}

## -------------------------------- RUN --------------------------------------
invisible(build_relative_table(cfg))
}


# Analysis section: RNA_expression_filters
RNA_expression_filters <- function() {


df.RNA_raw <- analysis_read_csv("rna_TPM.csv", header = TRUE, sep = ",")
RNA_Exp <- analysis_read_csv("rna_metadata.csv", header = TRUE, sep = ",")


dim(df.RNA_raw) #4276   76
df.RNA_raw[df.RNA_raw =="#VALUE!"]=NA #数据
which(!complete.cases(df.RNA_raw))


RNA_SampleInfo <- data.frame(colnames(df.RNA_raw[,5:76]),RNA_Exp$Time_h,RNA_Exp$Treatment,RNA_Exp$Group)
colnames(RNA_SampleInfo) <- c("Sample","Time","Treatment","Group")

## ============================================================
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
  library(grid)
  library(readr)
  library(scales)
})

## ============================================================
## ============================================================

data_file    <- "rna_TPM.csv"
expinfo_file <- "rna_metadata.csv"
outdir       <- "Step2_RNA_filtering"

low_min_tpm <- cfg$dnaa$low_min_tpm
low_max_tpm <- cfg$dnaa$low_max_tpm

## ratio = preCV / allCV
cv_theta_q <- cfg$dnaa$cv_theta_q

sensitivity_theta_q_vals <- c(0.80, 0.85, 0.90, 0.92, 0.94, 0.95, 0.96, 0.98)

groups_ordered <- c("Replicate1RNAseq", "Replicate2RNAseq", "Replicate3RNAseq")

time_unit <- "auto"

post_time_limit <- Inf

## ============================================================
## ============================================================
theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle    = element_text(hjust = 0.5, size = base_size - 0.5),
      legend.position  = "right",
      panel.grid       = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border     = element_rect(fill = NA, colour = "black", linewidth = 0.8),
      axis.text        = element_text(size = base_size),
      axis.title       = element_text(size = base_size + 1),
      axis.ticks       = element_line(colour = "black", linewidth = 0.6),
      axis.ticks.length = unit(1.0, "mm"),
      plot.margin      = margin(6, 6, 6, 6),
      strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 0.5),
      strip.text       = element_text(size = base_size, face = "bold")
    )
}

## ============================================================
## ============================================================
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## ============================================================
## ============================================================
message("== 1. 读入数据 ==")

df.RNA_raw <- analysis_read_csv(data_file, header = TRUE, sep = ",")
RNA_Exp    <- analysis_read_csv(expinfo_file, header = TRUE, sep = ",")

df.RNA_raw[df.RNA_raw == "#VALUE!"] <- NA

RNA_SampleInfo <- data.frame(
  Sample    = colnames(df.RNA_raw)[5:ncol(df.RNA_raw)],
  Time      = RNA_Exp$Time_h,
  Treatment = RNA_Exp$Treatment,
  Group     = RNA_Exp$Group
)

sample_cols <- intersect(colnames(df.RNA_raw), RNA_SampleInfo$Sample)
annot_cols  <- setdiff(colnames(df.RNA_raw), sample_cols)

df.RNA_raw[sample_cols] <- lapply(df.RNA_raw[sample_cols],
                                   function(x) suppressWarnings(as.numeric(x)))

gene_id <- if ("gene" %in% names(df.RNA_raw)) df.RNA_raw$gene else
  if ("locus_tag" %in% names(df.RNA_raw)) df.RNA_raw$locus_tag else
    sprintf("g%05d", seq_len(nrow(df.RNA_raw)))
rownames(df.RNA_raw) <- make.unique(as.character(gene_id))

anno_tbl <- df.RNA_raw[, annot_cols, drop = FALSE] %>%
  as_tibble(rownames = "Gene")

SampleInfo <- RNA_SampleInfo %>%
  filter(Sample %in% sample_cols) %>%
  mutate(Time = as.numeric(Time)) %>%
  arrange(Group, Time)

info_use <- SampleInfo %>%
  mutate(is_pre = (Treatment == "Pre_1aTc")) %>%
  filter(is_pre | Time <= post_time_limit) %>%
  arrange(Group, Time)

missing_groups <- setdiff(groups_ordered, unique(info_use$Group))
if (length(missing_groups) > 0)
  stop("以下组在数据中缺失：", paste(missing_groups, collapse = ", "))

sample_use <- info_use$Sample
expr_all0  <- as.matrix(df.RNA_raw[, sample_use, drop = FALSE])
orig_genes <- rownames(expr_all0)

time_vals    <- sort(unique(info_use$Time))
time_range   <- diff(range(time_vals, na.rm = TRUE))
time_min_gap <- min(diff(time_vals), na.rm = TRUE)

if (time_unit == "auto") {
  time_unit_detected <- if (time_range < 30) "hour" else "min"
  message(sprintf("  自动检测时间单位：%s（时间范围 %.2f，最小间隔 %.3f）",
                  time_unit_detected, time_range, time_min_gap))
  time_unit <- time_unit_detected
} else {
  message(sprintf("  手动指定时间单位：%s", time_unit))
}

grid_step_raw  <- time_min_gap / 3
n_target_points <- max(30, ceiling(time_range / grid_step_raw))
grid_step       <- time_range / n_target_points
time_unit_label <- if (time_unit == "hour") "h" else "min"

message(sprintf("  插值步长：%.4f %s，预期插值点数：%d",
                grid_step, time_unit_label, n_target_points))
message(sprintf("  初始基因数：%d，样本数：%d", nrow(expr_all0), ncol(expr_all0)))

## ============================================================
## ============================================================
message("== 2. 低表达过滤 ==")

calc_lowExpr_keep <- function(mat, info, min_th, max_th) {
  split(info, info$Group) %>%
    imap(function(df_g, g) {
      cols     <- df_g$Sample
      mm       <- mat[, cols, drop = FALSE]
      gene_min <- apply(mm, 1, min, na.rm = TRUE)
      gene_max <- apply(mm, 1, max, na.rm = TRUE)
      keep     <- !((gene_min < min_th) & (gene_max < max_th))
      setNames(keep, rownames(mat))
    })
}

keep_list_low <- calc_lowExpr_keep(expr_all0, info_use, low_min_tpm, low_max_tpm)
low_keep_ids  <- Reduce(intersect, lapply(keep_list_low, function(v) names(v)[v]))
low_drop_ids  <- setdiff(orig_genes, low_keep_ids)
expr_step1    <- expr_all0[low_keep_ids, , drop = FALSE]

message(sprintf("  低表达过滤后保留：%d 基因（剔除 %d）",
                length(low_keep_ids), length(low_drop_ids)))

## ============================================================
## ============================================================
message("== 3. 预波动过滤 ==")

calc_preFluct <- function(mat, info, theta_q, eps = 1e-8) {
  out <- split(info, info$Group) %>%
    imap(function(df_g, g) {
      pre_cols <- df_g %>% filter(Treatment == "Pre_1aTc") %>% pull(Sample)
      all_cols <- df_g$Sample

      if (length(pre_cols) < 2 || length(all_cols) < 2) {
        keep   <- rep(TRUE, nrow(mat))
        df_out <- tibble(Gene = rownames(mat), Group = g,
                         preCV = NA_real_, allCV = NA_real_,
                         ratio = NA_real_, theta_adapt = NA_real_, bad = FALSE)
        return(list(keep = setNames(keep, rownames(mat)), df = df_out))
      }

      mm_pre <- mat[, pre_cols, drop = FALSE]
      mm_all <- mat[, all_cols, drop = FALSE]

      preCV <- apply(mm_pre, 1, sd,   na.rm = TRUE) /
               pmax(apply(mm_pre, 1, mean, na.rm = TRUE), 1e-6)
      allCV <- apply(mm_all, 1, sd,   na.rm = TRUE) /
               pmax(apply(mm_all, 1, mean, na.rm = TRUE), 1e-6)
      ratio <- preCV / (allCV + eps)

      theta_adapt <- quantile(ratio[is.finite(ratio)],
                              probs = theta_q, na.rm = TRUE, names = FALSE)
      bad  <- is.finite(ratio) & (ratio > theta_adapt)

      df_out <- tibble(Gene = rownames(mat), Group = g,
                       preCV = as.numeric(preCV), allCV = as.numeric(allCV),
                       ratio = as.numeric(ratio),
                       theta_adapt = theta_adapt, bad = as.logical(bad))
      list(keep = setNames(!bad, rownames(mat)), df = df_out)
    })

  details_df <- map_dfr(out, "df")
  thresh_df  <- details_df %>%
    group_by(Group) %>%
    summarise(cv_theta_q  = theta_q,
              theta_adapt = first(theta_adapt),
              n_bad       = sum(bad, na.rm = TRUE),
              n_kept      = sum(!bad, na.rm = TRUE),
              pct_kept    = round(100 * n_kept / n(), 1),
              .groups     = "drop")

  keep_ids <- Reduce(intersect,
                     lapply(out, function(x) names(x$keep)[x$keep]))
  list(keep_ids = keep_ids, details = details_df, thresh = thresh_df)
}

pf_res     <- calc_preFluct(expr_step1, info_use, cv_theta_q)
pf_keep    <- pf_res$keep_ids
pf_drop    <- setdiff(rownames(expr_step1), pf_keep)
expr_step2 <- expr_step1[pf_keep, , drop = FALSE]

message(sprintf("  预波动过滤后保留：%d 基因（剔除 %d）",
                length(pf_keep), length(pf_drop)))

message("  → 敏感性分析扫描...")
sens_df <- map_dfr(sensitivity_theta_q_vals, function(tq) {
  res <- calc_preFluct(expr_step1, info_use, tq)
  tibble(theta_q   = tq,
         n_kept    = length(res$keep_ids),
         n_removed = nrow(expr_step1) - length(res$keep_ids),
         pct_kept  = round(100 * length(res$keep_ids) / nrow(expr_step1), 1))
})

## ============================================================
## ============================================================
gene_status <- tibble(Gene = orig_genes) %>%
  mutate(
    keep_after_lowExpr  = Gene %in% low_keep_ids,
    keep_after_preFluct = Gene %in% pf_keep,
    final_keep          = Gene %in% pf_keep,
    drop_step           = case_when(
      !keep_after_lowExpr  ~ "step1_lowExpr",
      !keep_after_preFluct ~ "step2_preFluct",
      TRUE                 ~ "kept"),
    drop_reason = case_when(
  !keep_after_lowExpr  ~
    sprintf("max < %g TPM AND min < %g TPM in >=1 replicate group",
            low_max_tpm, low_min_tpm),
  !keep_after_preFluct ~
    sprintf("ratio > P%.0f(theta_adapt) in >=1 replicate group",
            cv_theta_q * 100),
  TRUE ~ "kept")
  ) %>%
  left_join(anno_tbl, by = "Gene")

## ============================================================
## ============================================================
message("== 5. 导出 CSV 文件 ==")

write_expr_csv <- function(mat, label) {
  as_tibble(mat, rownames = "Gene") %>%
    left_join(anno_tbl, by = "Gene") %>%
    relocate(Gene, all_of(annot_cols), .before = 1) %>%
    write.csv(file.path(outdir, sprintf("genes_%s.csv", label)), row.names = FALSE)
}

filter_funnel <- tibble(
  step        = c("step0_raw", "step1_lowExpr", "step2_preFluct"),
  n_genes     = c(length(orig_genes), length(low_keep_ids), length(pf_keep)),
  n_removed   = c(0,
                  length(orig_genes) - length(low_keep_ids),
                  length(low_keep_ids) - length(pf_keep)),
  pct_kept    = round(100 * c(length(orig_genes), length(low_keep_ids),
                               length(pf_keep)) / length(orig_genes), 1),
 params_used = c("-",
  sprintf("max < %g AND min < %g (TPM); removed if condition met in >=1 replicate group",
          low_max_tpm, low_min_tpm),
  sprintf("ratio > quantile(ratio, %.2f) (per-group adaptive theta); removed if exceeded in >=1 replicate group",
          cv_theta_q))
)
write.csv(filter_funnel, file.path(outdir, "A_filter_funnel_summary.csv"), row.names = FALSE)

preCV_dist <- pf_res$details %>%
  select(Gene, Group, preCV) %>%
  group_by(Group) %>%
  mutate(preCV_rank_pct = rank(preCV, na.last = "keep") /
           sum(!is.na(preCV)) * 100) %>%
  ungroup()
write.csv(preCV_dist, file.path(outdir, "B1_preCV_distribution_perGroup.csv"), row.names = FALSE)

allCV_dist <- pf_res$details %>%
  select(Gene, Group, allCV) %>%
  group_by(Group) %>%
  mutate(allCV_rank_pct = rank(allCV, na.last = "keep") /
           sum(!is.na(allCV)) * 100) %>%
  ungroup()
write.csv(allCV_dist, file.path(outdir, "B2_allCV_distribution_perGroup.csv"), row.names = FALSE)

ratio_dist <- pf_res$details %>%
  select(Gene, Group, ratio, theta_adapt, bad) %>%
  group_by(Group) %>%
  mutate(ratio_rank_pct = rank(ratio, na.last = "keep") /
           sum(!is.na(ratio)) * 100) %>%
  ungroup()
write.csv(ratio_dist, file.path(outdir, "B3_ratio_distribution_perGroup.csv"), row.names = FALSE)

write.csv(pf_res$thresh, file.path(outdir, "B4_threshold_summary_perGroup.csv"), row.names = FALSE)

write_expr_csv(expr_all0,  "C0_step0_all")

as_tibble(expr_step1, rownames = "Gene") %>%
  left_join(anno_tbl, by = "Gene") %>%
  relocate(Gene, all_of(annot_cols), .before = 1) %>%
  write.csv(file.path(outdir, "genes_C1_step1_after_lowExpr.csv"), row.names = FALSE)

pf_wide <- pf_res$details %>%
  select(Gene, Group, preCV, allCV, ratio, theta_adapt, bad) %>%
  pivot_wider(names_from = Group,
              values_from = c(preCV, allCV, ratio, theta_adapt, bad),
              names_sep = "__")
as_tibble(expr_step2, rownames = "Gene") %>%
  left_join(anno_tbl, by = "Gene") %>%
  left_join(pf_wide,  by = "Gene") %>%
  relocate(Gene, all_of(annot_cols), .before = 1) %>%
  write.csv(file.path(outdir, "genes_C2_step2_after_preFluct.csv"), row.names = FALSE)

write.csv(gene_status, file.path(outdir, "genes_C3_all_status.csv"), row.names = FALSE)

write.csv(sens_df, file.path(outdir, "D1_sensitivity_theta_q.csv"), row.names = FALSE)

message("  CSV 导出完成。")

## ============================================================
## ============================================================
message("== 6. 过滤参数分布图 ==")

pal3 <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A"), groups_ordered)

save_plot <- function(p, name, w, h, dpi = 200) {
  ggsave(file.path(outdir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(outdir, paste0(name, ".png")), p, width = w, height = h, dpi = dpi)
}

p_funnel <- filter_funnel %>%
  mutate(step = factor(step,
    levels = c("step0_raw", "step1_lowExpr", "step2_preFluct"),
    labels = c("Raw", "After\nlow-expr filter", "After\npre-fluct filter"))) %>%
  ggplot(aes(step, n_genes)) +
  geom_col(fill = "#4292C6", width = 0.55, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n_genes, pct_kept)),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Fig F1. Gene filtering funnel", x = NULL, y = "Number of genes") +
  theme_pub(10) + theme(aspect.ratio = 1.2)
save_plot(p_funnel, "FigF1_filter_funnel", 4, 4.5)

p_preCV <- preCV_dist %>%
  filter(!is.na(preCV)) %>%
  ggplot(aes(preCV, color = Group, fill = Group)) +
  geom_density(alpha = 0.20, linewidth = 0.7) +
  scale_color_manual(values = pal3) +
  scale_fill_manual(values  = pal3) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Fig F2. preCV distribution per replicate",
       subtitle = "preCV = sd(pre-perturbation samples) / mean(pre-perturbation samples)",
       x = "preCV", y = "Density") +
  theme_pub(10) + theme(aspect.ratio = 0.7)
save_plot(p_preCV, "FigF2_preCV_distribution", 6, 4.5)

p_allCV <- allCV_dist %>%
  filter(!is.na(allCV)) %>%
  ggplot(aes(allCV, color = Group, fill = Group)) +
  geom_density(alpha = 0.20, linewidth = 0.7) +
  scale_color_manual(values = pal3) +
  scale_fill_manual(values  = pal3) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Fig F3. allCV distribution per replicate",
       subtitle = "allCV = sd(all samples) / mean(all samples)",
       x = "allCV", y = "Density") +
  theme_pub(10) + theme(aspect.ratio = 0.7)
save_plot(p_allCV, "FigF3_allCV_distribution", 6, 4.5)

ratio_xlim <- ratio_dist %>%
  filter(is.finite(ratio)) %>% pull(ratio) %>%
  quantile(0.99, na.rm = TRUE)

p_ratio <- ratio_dist %>%
  filter(is.finite(ratio)) %>%
  ggplot(aes(ratio, color = Group, fill = Group)) +
  geom_density(alpha = 0.20, linewidth = 0.7) +
  geom_vline(data = pf_res$thresh,
             aes(xintercept = theta_adapt, color = Group),
             linetype = "dashed", linewidth = 0.8) +
  geom_text(data = pf_res$thresh,
            aes(x = theta_adapt, y = Inf,
                label = sprintf("P%.0f=%.3f", cv_theta_q * 100, theta_adapt),
                color = Group),
            vjust = 1.6, hjust = -0.08, size = 2.8, show.legend = FALSE) +
  scale_color_manual(values = pal3) +
  scale_fill_manual(values  = pal3) +
  coord_cartesian(xlim = c(0, ratio_xlim)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Fig F4. ratio (preCV / allCV) distribution per replicate",
       subtitle = sprintf(
         "Dashed lines = theta_adapt = P%.0f of ratio distribution (per replicate)",
         cv_theta_q * 100),
       x = "ratio  (preCV / allCV)", y = "Density") +
  theme_pub(10) + theme(aspect.ratio = 0.7)
save_plot(p_ratio, "FigF4_ratio_distribution", 6, 4.5)

p_scatter <- pf_res$details %>%
  filter(is.finite(preCV), is.finite(allCV)) %>%
  mutate(status = ifelse(bad, "Removed", "Kept")) %>%
  ggplot(aes(allCV, preCV, color = status)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted", color = "grey40", linewidth = 0.5) +
  facet_wrap(~Group, nrow = 1) +
  scale_color_manual(values = c("Kept" = "#4292C6", "Removed" = "#E41A1C"),
                     name = NULL) +
  labs(title = "Fig F5. preCV vs allCV scatter (per replicate)",
       subtitle = "Dotted diagonal: ratio = 1 (preCV = allCV). Red = removed.",
       x = "allCV", y = "preCV") +
  theme_pub(9) + theme(aspect.ratio = 1, legend.position = "bottom")
save_plot(p_scatter, "FigF5_preCV_vs_allCV_scatter", 9, 4.2)

p_sens <- sens_df %>%
  ggplot(aes(theta_q, n_kept)) +
  geom_line(color = "#2171B5", linewidth = 0.9) +
  geom_point(color = "#2171B5", size = 2.2) +
  geom_point(data = filter(sens_df, theta_q == cv_theta_q),
             color = "#E41A1C", size = 4, shape = 21,
             fill = "#E41A1C", stroke = 1) +
  geom_label(aes(label = sprintf("n=%d\n%.1f%%", n_kept, pct_kept)),
             size = 2.3, label.padding = unit(0.15, "lines"),
             label.size = 0.2, vjust = -0.3, fill = "white") +
  scale_x_continuous(breaks = sensitivity_theta_q_vals,
                     labels = sprintf("P%.0f", sensitivity_theta_q_vals * 100)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(title = "Fig F6. Sensitivity analysis: cv_theta_q",
       subtitle = sprintf("Red dot = selected threshold (P%.0f); input n = %d genes",
                          cv_theta_q * 100, nrow(expr_step1)),
       x = "Threshold quantile (theta_q)",
       y = "Genes retained after pre-fluct filter") +
  theme_pub(10) + theme(aspect.ratio = 0.65,
                        axis.text.x = element_text(size = 8))
save_plot(p_sens, "FigF6_sensitivity_theta_q", 7, 4.5)

p_panel <- (p_funnel | p_ratio | p_sens) +
  plot_annotation(
    title = "RNA-seq gene filtering summary",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))
  )
save_plot(p_panel, "FigF7_filter_summary_panel", 15, 4.8)

message("  过滤图输出完成。")

## ============================================================
## ============================================================
message("== 7. S1 重复性分析 ==")

mats_s1 <- setNames(groups_ordered, groups_ordered) %>%
  lapply(function(g) {
    cols <- info_use %>% filter(Group == g) %>% pull(Sample)
    expr_step2[, cols, drop = FALSE]
  })

cor_within <- function(M, grp) {
  if (ncol(M) < 2) return(tibble())
  C   <- suppressWarnings(cor(M, use = "pairwise.complete.obs"))
  idx <- which(upper.tri(C), arr.ind = TRUE)
  tibble(Group = paste("Within", grp),
         SampleA = colnames(M)[idx[, 1]],
         SampleB = colnames(M)[idx[, 2]],
         r = C[idx])
}

cor_cross <- function(M1, M2, label) {
  common <- intersect(rownames(M1), rownames(M2))
  if (length(common) < 2) return(tibble())
  C <- suppressWarnings(
    cor(M1[common, , drop = FALSE], M2[common, , drop = FALSE],
        use = "pairwise.complete.obs"))
  expand.grid(SampleA = colnames(M1), SampleB = colnames(M2),
              stringsAsFactors = FALSE) %>%
    mutate(r = as.vector(C[cbind(match(SampleA, rownames(C)),
                                  match(SampleB, colnames(C)))]),
           Group = paste0("Cross (", label, ")")) %>%
    select(Group, SampleA, SampleB, r) %>% as_tibble()
}

desired_levels <- c(
  "Cross (Replicate1RNAseq vs Replicate2RNAseq)", "Cross (Replicate1RNAseq vs Replicate3RNAseq)",
  "Cross (Replicate2RNAseq vs Replicate3RNAseq)",
  "Within Replicate1RNAseq", "Within Replicate2RNAseq", "Within Replicate3RNAseq")

df_s1 <- bind_rows(
  cor_cross(mats_s1[["Replicate1RNAseq"]], mats_s1[["Replicate2RNAseq"]], "Replicate1RNAseq vs Replicate2RNAseq"),
  cor_cross(mats_s1[["Replicate1RNAseq"]], mats_s1[["Replicate3RNAseq"]], "Replicate1RNAseq vs Replicate3RNAseq"),
  cor_cross(mats_s1[["Replicate2RNAseq"]], mats_s1[["Replicate3RNAseq"]], "Replicate2RNAseq vs Replicate3RNAseq"),
  cor_within(mats_s1[["Replicate1RNAseq"]], "Replicate1RNAseq"),
  cor_within(mats_s1[["Replicate2RNAseq"]], "Replicate2RNAseq"),
  cor_within(mats_s1[["Replicate3RNAseq"]], "Replicate3RNAseq")
) %>% mutate(Group = factor(Group, levels = desired_levels))

df_s1_summary <- df_s1 %>%
  group_by(Group) %>%
  summarise(n = n(), median = median(r, na.rm = TRUE),
            IQR = IQR(r, na.rm = TRUE),
            min = min(r, na.rm = TRUE), max = max(r, na.rm = TRUE),
            .groups = "drop")

write.csv(df_s1,         file.path(outdir, "SI_S1_correlations_long.csv"),    row.names = FALSE)
write.csv(df_s1_summary, file.path(outdir, "SI_S1_correlations_summary.csv"), row.names = FALSE)

pal6 <- setNames(brewer.pal(6, "Set2"), desired_levels)

p_s1 <- ggplot(df_s1, aes(Group, r, fill = Group)) +
  geom_violin(trim = FALSE, width = 0.85, color = "black",
              alpha = 0.5, linewidth = 0.25) +
  geom_boxplot(width = 0.16, outlier.size = 0, alpha = 1,
               color = "black", linewidth = 0.25, coef = Inf) +
  stat_summary(fun = median, geom = "point",
               size = 1.8, color = "black", shape = 21, fill = "black") +
  scale_fill_manual(values = pal6) +
  scale_y_continuous(limits = c(0.7, 1), breaks = seq(0.7, 1, 0.1),
                     expand = c(0, 0)) +
  scale_x_discrete(drop = FALSE,
                   labels = function(x) gsub(" \\(", "\n(", x)) +
  labs(y = "Pearson r (genes as features)", x = NULL,
       title = "S1. Sample-level correlation distributions") +
  theme_pub(10) +
  theme(aspect.ratio = 0.75, legend.position = "none",
        axis.text.x = element_text(size = 7.5))

save_plot(p_s1, "SI_Fig_S1_correlation_distributions", 8.5, 4.5)
message("  S1 完成。")

## ============================================================
## ============================================================
message("== 8. S3 PCA + CCF ==")

mats_s3 <- setNames(groups_ordered, groups_ordered) %>%
  lapply(function(g) {
    df_g <- info_use %>% filter(Group == g) %>% arrange(Time)
    expr_step2[, df_g$Sample, drop = FALSE]
  })

X_pca  <- t(do.call(cbind, mats_s3))
keep_v <- apply(X_pca, 2, function(v) var(v, na.rm = TRUE) > 0)
X_pca  <- X_pca[, keep_v, drop = FALSE]
pca    <- prcomp(X_pca, center = TRUE, scale. = TRUE)

var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

# Match sample metadata to the PCA score rows by ID, independently of group order.
sample_index <- match(rownames(pca$x), info_use$Sample)
if(anyNA(sample_index) || anyDuplicated(info_use$Sample))
  stop("PCA sample IDs do not match the sample metadata uniquely")
meta_pca <- info_use[sample_index, ] %>% select(Sample, Group, Time)
pc_scores <- as_tibble(pca$x[, 1:2]) %>%
  setNames(c("PC1", "PC2")) %>%
  bind_cols(meta_pca) %>%
  arrange(Group, Time)

write.csv(pc_scores, file.path(outdir, "SI_S3_PC_scores_by_sample.csv"), row.names = FALSE)

pal3_s3 <- setNames(brewer.pal(max(3, length(groups_ordered)), "Set1"),
                    sort(groups_ordered))

pc1_label <- sprintf("PC1 (%.1f%% variance)", var_exp[1])

p_traj <- pc_scores %>%
  arrange(Group, Time) %>%
  ggplot(aes(Time, PC1, color = Group, group = Group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = pal3_s3) +
  scale_x_continuous(expand = c(0.02, 0)) +
  labs(title = "S3. PC1 trajectories over time",
       x = sprintf("Time (%s)", time_unit_label),
       y = pc1_label) +
  theme_pub(10) + theme(aspect.ratio = 0.6)

save_plot(p_traj, "SI_Fig_S3_PC1_trajectories", 7, 4)

grid_seq <- seq(min(meta_pca$Time, na.rm = TRUE),
                max(meta_pca$Time, na.rm = TRUE),
                by = grid_step)

interp_pc1 <- function(grp) {
  sub <- pc_scores %>% filter(Group == grp) %>% arrange(Time)
  approx(x = sub$Time, y = sub$PC1, xout = grid_seq, rule = 2)$y
}

grid_series <- tibble(Time = grid_seq) %>%
  bind_cols(setNames(lapply(groups_ordered, interp_pc1), groups_ordered))

write.csv(grid_series, file.path(outdir, "SI_S3_PC1_grid_series.csv"), row.names = FALSE)
message(sprintf("  插值序列长度：%d 点", length(grid_seq)))

vline_col <- "#D55E00"
ccf_pairs <- list(c("Replicate1RNAseq","Replicate2RNAseq"),
                  c("Replicate1RNAseq","Replicate3RNAseq"),
                  c("Replicate2RNAseq","Replicate3RNAseq"))

ccf_plots <- lapply(ccf_pairs, function(pair) {
  g1 <- pair[1]; g2 <- pair[2]
  x1 <- grid_series[[g1]]; x2 <- grid_series[[g2]]

  if (all(is.na(x1)) || all(is.na(x2))) {
    warning(sprintf("CCF: %s vs %s 数据全为 NA，跳过。", g1, g2))
    return(NULL)
  }
  if (length(x1) < 5) {
    warning(sprintf("CCF: %s vs %s 插值点数仅 %d，CCF 结果不可靠，请检查 grid_step 设置。",
                    g1, g2, length(x1)))
  }

  cc       <- ccf(x1, x2, plot = FALSE, na.action = na.pass)
  lag_best <- cc$lag[which.max(cc$acf)] * grid_step
  r_best   <- max(cc$acf, na.rm = TRUE)

  cc_df   <- tibble(lag_time = as.numeric(cc$lag) * grid_step,
                    acf      = as.numeric(cc$acf))
  best_pt <- cc_df[which.max(cc_df$acf), ]

  lag_label <- if (time_unit == "hour") {
    sprintf("%.2f h", lag_best)
  } else {
    sprintf("%d min", as.integer(round(lag_best)))
  }

  write.csv(cc_df,
    file.path(outdir, sprintf("SI_S3_CCF_curve_%s_vs_%s.csv", g1, g2)),
    row.names = FALSE)
  write.csv(tibble(g1 = g1, g2 = g2, max_r = r_best, lag = lag_best,
                   lag_unit = time_unit_label),
    file.path(outdir, sprintf("SI_S3_CCF_best_%s_vs_%s.csv", g1, g2)),
    row.names = FALSE)

  ggplot(cc_df, aes(lag_time, acf)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
    geom_line(color = "grey40", linewidth = 0.5, alpha = 0.7) +
    geom_point(size = 1.5, alpha = 0.85, color = "grey30") +
    geom_vline(xintercept = best_pt$lag_time,
               linetype = "dashed", color = vline_col, linewidth = 0.7) +
    geom_point(data = best_pt, shape = 21, size = 3, stroke = 1,
               fill = vline_col, color = vline_col) +
    annotate("text",
             x     = best_pt$lag_time,
             y     = best_pt$acf,
             label = sprintf("max r = %.2f\nlag = %s", r_best, lag_label),
             vjust = -0.5, hjust = 0.5, size = 3.2, color = vline_col) +
    scale_x_continuous(expand = c(0.03, 0)) +
    scale_y_continuous(limits = c(
        min(0, min(cc_df$acf, na.rm = TRUE)) - 0.05, 1.18),
      breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
    labs(title    = sprintf("S3. CCF: %s vs %s", g1, g2),
         subtitle = sprintf("n = %d interpolated points, step = %.4f %s",
                            length(grid_seq), grid_step, time_unit_label),
         x = sprintf("Lag (%s)  [+ = %s lags %s]", time_unit_label, g2, g1),
         y = "Cross-correlation") +
    theme_pub(10) + theme(aspect.ratio = 0.65)
})

for (i in seq_along(ccf_pairs)) {
  if (is.null(ccf_plots[[i]])) next
  g1 <- ccf_pairs[[i]][1]; g2 <- ccf_pairs[[i]][2]
  save_plot(ccf_plots[[i]],
            sprintf("SI_Fig_S3_CCF_%s_vs_%s", g1, g2), 6, 4)
}

valid_ccf <- Filter(Negate(is.null), ccf_plots)
if (length(valid_ccf) > 0) {
  p_s3_panel <- (p_traj /
    (wrap_plots(valid_ccf, nrow = 1))) +
    plot_annotation(
      title = "S3. PC1 dynamics and cross-correlation between replicates",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))
    )
  save_plot(p_s3_panel, "SI_Fig_S3_panel", 14, 8)
}

message("  S3 完成。")

## ============================================================
## ============================================================
message("\n========== 全部分析完成 ==========")
message("输出目录：", normalizePath(outdir))
message(sprintf("\n── 检测到时间单位：%s，插值步长：%.4f %s，插值点数：%d",
                time_unit, grid_step, time_unit_label, length(grid_seq)))
message("\n── 过滤参数分布图 ──")
message("  FigF1  过滤漏斗")
message("  FigF2  preCV 密度分布（纯分布，无多余标注）")
message("  FigF3  allCV 密度分布")
message("  FigF4  ratio 密度分布 + theta_adapt 截断线（各组独立）")
message("  FigF5  preCV vs allCV 散点（标注剔除基因）")
message("  FigF6  敏感性折线图（8个 theta_q 扫描点）")
message("  FigF7  综合面板（漏斗 + ratio + 敏感性）")
message("\n── GraphPad 绘图数据 ──")
message("  A_filter_funnel_summary.csv")
message("  B1_preCV_distribution_perGroup.csv")
message("  B2_allCV_distribution_perGroup.csv")
message("  B3_ratio_distribution_perGroup.csv")
message("  B4_threshold_summary_perGroup.csv")
message("  D1_sensitivity_theta_q.csv")
message("\n── 逐步过滤基因表 ──")
message("  genes_C0_step0_all.csv")
message("  genes_C1_step1_after_lowExpr.csv")
message("  genes_C2_step2_after_preFluct.csv")
message("  genes_C3_all_status.csv")
message("\n── S1/S3 质控图 ──")
message("  SI_Fig_S1_correlation_distributions.pdf/png")
message("  SI_Fig_S3_PC1_trajectories.pdf/png")
message("  SI_Fig_S3_CCF_*.pdf/png  （各对 + 拼合面板）")


}

# Analysis section: RNA_heatmaps
RNA_heatmaps <- function() {


# ============================================================
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(scales)
  library(readr)
})

# ============================================================
# ============================================================

RAW_TPM_CSV  <- "rna_TPM.csv"   # 原始TPM矩阵
EXPINFO_CSV  <- "rna_metadata.csv" # 样品信息
GENE_LIST_CSV <- "filtered_genes.csv" # 2918基因列表

COL_SAMPLEID <- "Colnames"    # 样品列名列
COL_TIME     <- "Time_h"      # 时间（小时）列
COL_GROUP    <- "Group"       # 实验组列
COL_TREAT    <- "Treatment"   # 处理类型列
PRE_LABEL    <- "Pre_1aTc"    # Pre时间点的Treatment标签

GROUPS_ORDERED <- c("Replicate3RNAseq", "Replicate1RNAseq", "Replicate2RNAseq")

VALUE_TYPES <- c("relTPM", "TPM")

SCALE_MODES <- c("log2", "zscore")

#   ZRANGE_MODE:
ZRANGE_MODE <- "shared"

#              TPM_log2 = NULL, TPM_zscore = c(-3, 3))
ZRANGE_MANUAL <- list(
  relTPM_log2   = NULL,     # 自动推断，通常约 ±1.5
  relTPM_zscore = c(-3, 3), # z-score 固定
  TPM_log2      = NULL,     # 自动推断
  TPM_zscore    = c(-3, 3)  # z-score 固定
)

PALETTE <- c("#3B4CC0", "#88A1D2", "#F7F7F7", "#F3A982", "#B40426")


GENE_LABEL_OFFSET_MM <- 2.5   # 引线从热图右边界向右延伸的距离（mm等效单位）
GENE_LABEL_SIZE      <- 2.2   # 基因名字体大小（ggplot单位，约6pt）
GENE_TICK_LWD        <- 0.3   # 引线线宽
GENE_TICK_COL        <- "black"

N_XLABELS <- 8

SAMPLE_TICK_LEN <- 50
SAMPLE_TICK_LWD <- 0.4
SAMPLE_TICK_COL <- "black"

HEAT_W_MM    <- 150    # 单组热图宽度（mm）
HEAT_H_MM    <- 130   # 单组热图高度（mm）
TOPBAR_H     <- 0.04  # 顶部时间色条占总高度比例
PANEL_GAP_MM <- 3     # 三组拼合图的组间距（mm）

OUT_DIR <- "Step3_RNA_heatmaps"

# ============================================================
# ============================================================
mod_load <- function(raw_csv, expinfo_csv, genelist_csv,
                     col_sampleid, col_time, col_group, col_treat) {

  message("== 读入数据 ==")

  mat_raw <- analysis_read_csv(raw_csv, header = TRUE, check.names = FALSE,
                      stringsAsFactors = FALSE)
  annot_cols <- c("gene", "locus_tag", "product", "Rel.Pos")
  stopifnot(all(annot_cols %in% colnames(mat_raw)))

  info_raw <- analysis_read_csv(expinfo_csv, header = TRUE, check.names = FALSE,
                       stringsAsFactors = FALSE)

  dup_cols <- duplicated(colnames(info_raw))
  if (any(dup_cols)) {
    message(sprintf("  [修复] ExpInfo 含重复列名：%s，已自动去除",
                    paste(colnames(info_raw)[dup_cols], collapse = ", ")))
    info_raw <- info_raw[, !dup_cols, drop = FALSE]
  }

  info_raw[[col_time]] <- suppressWarnings(as.numeric(info_raw[[col_time]]))

  gene_list <- analysis_read_csv(genelist_csv, header = TRUE, check.names = FALSE,
                        stringsAsFactors = FALSE)
  if (!"gene" %in% colnames(gene_list)) {
    if ("Gene" %in% colnames(gene_list)) gene_list$gene <- gene_list$Gene
    else stop("2918基因列表中找不到 'gene' 或 'Gene' 列")
  }
  genes_keep <- unique(gene_list$gene)
  message(sprintf("  2918基因列表：%d 个唯一基因名", length(genes_keep)))

  mat_filt <- mat_raw %>% filter(gene %in% genes_keep)
  message(sprintf("  原始矩阵过滤后：%d 行（目标 %d）",
                  nrow(mat_filt), length(genes_keep)))
  if (nrow(mat_filt) != length(genes_keep)) {
    miss <- setdiff(genes_keep, mat_filt$gene)
    warning(sprintf("  %d 个基因在原始矩阵中未找到：%s",
                    length(miss), paste(head(miss, 5), collapse = ", ")))
  }

  valid_cols <- info_raw[[col_sampleid]][
    info_raw[[col_sampleid]] %in% colnames(mat_filt)
  ]
  message(sprintf("  有效样本列：%d 个", length(valid_cols)))

  anno_tbl <- mat_filt %>%
    select(all_of(annot_cols)) %>%
    mutate(Rel.Pos = as.numeric(Rel.Pos))

  expr_mat <- as.matrix(mat_filt[, valid_cols, drop = FALSE])
  rownames(expr_mat) <- mat_filt$gene
  storage.mode(expr_mat) <- "numeric"

  info_use <- info_raw %>%
  filter(.data[[col_sampleid]] %in% valid_cols) %>%
  mutate(
    Time  = as.numeric(.data[[col_time]]),
    Group = as.character(.data[[col_group]]),
    Treat = as.character(.data[[col_treat]])
  ) %>%
  arrange(Group, Time)

  list(
    expr_mat = expr_mat,
    anno_tbl = anno_tbl,
    info_use = info_use,
    genes    = mat_filt$gene
  )
}

# ============================================================
# ============================================================
mod_yspec_rna <- function(anno_tbl, target_genes = NULL) {

  anno_ord <- anno_tbl %>%
    arrange(desc(Rel.Pos)) %>%
    mutate(
      gene_id = gene,
      y    = seq(n(), 1),
      ymin = y - 0.5,
      ymax = y + 0.5
    )

  n_genes <- nrow(anno_ord)
  message(sprintf("  Y轴：%d 基因，Rel.Pos +1 在顶部(y=%d)，-1 在底部(y=1)",
                  n_genes, n_genes))

  target_rp <- c(1.0, 0.5, 0.0, -0.5, -1.0)   # 从顶到底
  yspec_breaks <- lapply(seq_along(target_rp), function(i) {
    rp <- target_rp[i]
    anno_ord %>%
      mutate(dist = abs(Rel.Pos - rp)) %>%
      slice_min(dist, n = 1, with_ties = FALSE) %>%
      mutate(target_label = rp)        # 记录精确目标值
  }) %>%
    bind_rows() %>%
    distinct(target_label, .keep_all = TRUE) %>%   # 按目标值去重（非基因名）
    arrange(desc(y))                               # 从顶到底排列

  message("  Y轴刻度：")
  for (i in seq_len(nrow(yspec_breaks))) {
    r <- yspec_breaks[i, ]
    message(sprintf("    Rel.Pos=%+.1f → %s (实际=%.4f, y=%d)",
                    r$target_label, r$gene_id, r$Rel.Pos, r$y))
  }

  gene_annot <- NULL
  if (!is.null(target_genes) && length(target_genes) > 0) {
    gene_annot <- anno_ord %>%
      filter(gene_id %in% target_genes) %>%
      select(gene_id, Rel.Pos, y, ymin, ymax) %>%
      arrange(desc(y))
    found <- gene_annot$gene_id
    miss  <- setdiff(target_genes, found)
    message(sprintf("  目标基因：%d/%d 个找到", length(found), length(target_genes)))
    if (length(miss) > 0)
      message(sprintf("  未找到：%s", paste(miss, collapse = ", ")))
  }

  list(
    yspec       = anno_ord,         # 全部基因的y坐标
    ybreaks     = yspec_breaks,     # 刻度数据框
    gene_annot  = gene_annot,       # 目标基因坐标
    n_genes     = n_genes
  )
}

# ============================================================
# ============================================================
mod_xspec <- function(time_vec) {
  x <- sort(unique(as.numeric(time_vec)))
  x <- x[is.finite(x)]
  if (length(x) == 1)
    return(tibble(Time = x, xmin = x - 0.5, xmax = x + 0.5))
  mid  <- (x[-1] + x[-length(x)]) / 2
  xmin <- c(x[1] - (mid[1] - x[1]), mid)
  xmax <- c(mid, x[length(x)] + (x[length(x)] - mid[length(mid)]))
  tibble(Time = x, xmin = xmin, xmax = xmax)
}

# ============================================================
# ============================================================
mod_rel_tpm <- function(mat_group, info_group,
                        col_sampleid, pre_label,
                        eps = 1e-8) {

  pre_cols <- info_group[[col_sampleid]][
    info_group$Treat == pre_label
  ]

  if (length(pre_cols) == 0) {
    message("  [警告] 本组无 Pre 样品，直接使用原始 TPM")
    return(list(mat = mat_group, value_type = "TPM",
                n_degraded = 0))
  }

  pre_mean <- rowMeans(mat_group[, pre_cols, drop = FALSE], na.rm = TRUE)

  invalid <- !is.finite(pre_mean) | pre_mean <= eps
  n_degraded <- sum(invalid)
  if (n_degraded > 0) {
    message(sprintf("  [注意] %d 个基因 Pre 均值≤0，rel.TPM 设为 NA", n_degraded))
  }

  pre_mean[invalid] <- NA

  mat_rel <- sweep(mat_group, 1, pre_mean, "/")
  mat_rel[!is.finite(mat_rel)] <- NA

  list(mat = mat_rel, value_type = "relTPM",
       n_degraded = n_degraded)
}

# ============================================================
# ============================================================
mod_scale <- function(mat, mode, value_type) {

  if (mode == "log2") {
    pseudo <- if (value_type == "relTPM") 0 else 1
    mat <- log2(mat + pseudo)
    mat[!is.finite(mat)] <- NA

  } else if (mode == "zscore") {
    mu  <- rowMeans(mat, na.rm = TRUE)
    sdv <- apply(mat, 1, sd, na.rm = TRUE)
    sdv[!is.finite(sdv) | sdv < 1e-10] <- 1
    mat <- sweep(mat, 1, mu,  "-")
    mat <- sweep(mat, 1, sdv, "/")
    mat[!is.finite(mat)] <- 0
  }

  mat
}

# ============================================================
# ============================================================
mod_zrange <- function(mat, value_type, scale_mode) {

  if (scale_mode == "zscore") return(c(-3, 3))

  vals <- mat[is.finite(mat)]
  if (length(vals) == 0) return(c(-3, 3))

  if (value_type == "relTPM") {
    lim <- quantile(abs(vals), 0.99, na.rm = TRUE)
    lim <- ceiling(lim * 4) / 4
    return(c(-lim, lim))
  }

  lo <- floor(quantile(vals, 0.01, na.rm = TRUE) * 4) / 4
  hi <- ceiling(quantile(vals, 0.99, na.rm = TRUE) * 4) / 4
  c(lo, hi)
}

# ============================================================
# ============================================================
mod_plotdata <- function(mat_scaled, info_plot, col_sampleid,
                         xspec, yspec_obj) {

  yspec <- yspec_obj$yspec   # 含 gene_id, y, ymin, ymax

  cols    <- info_plot[[col_sampleid]]
  mat_sub <- mat_scaled[, cols, drop = FALSE]

  as.data.frame(mat_sub, check.names = FALSE) %>%
    rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id,
                 names_to  = col_sampleid,
                 values_to = "value") %>%
    left_join(
      info_plot %>% select(all_of(c(col_sampleid, "Time"))),
      by = col_sampleid
    ) %>%
    left_join(xspec, by = "Time") %>%
    left_join(
      yspec %>% select(gene_id, y, ymin, ymax),
      by = "gene_id"
    ) %>%
    filter(is.finite(xmin), is.finite(ymin))
}

# ============================================================
# ============================================================
mod_plot_single <- function(plotdata, xspec, yspec_obj,
                            z_range, palette, n_xlabels,
                            title       = NULL,
                            value_label = "value",
                            topbar_h    = 0.03,
                            sample_tick_len = 50,
                            sample_tick_lwd = 0.4,
                            sample_tick_col = "black",
                            show_gene_annot  = TRUE,
                            gene_label_size  = 2.2,
                            gene_tick_lwd    = 0.3,
                            gene_tick_col    = "black",
                            gene_label_xoffset = 0.05,  # x方向偏移（数据单位）
                            show_legend = TRUE) {

  yspec      <- yspec_obj$yspec
  ybreaks    <- yspec_obj$ybreaks
  gene_annot <- yspec_obj$gene_annot

  x_lim <- range(c(xspec$xmin, xspec$xmax), na.rm = TRUE)
  y_lim <- range(c(yspec$ymin, yspec$ymax), na.rm = TRUE)

  xb <- breaks_pretty(n = n_xlabels)(range(xspec$Time))
  xb <- xb[is.finite(xb) & xb >= x_lim[1] & xb <= x_lim[2]]

  sample_times <- sort(unique(xspec$Time))
  tick_df <- tibble(
    x    = sample_times,
    xend = sample_times,
    y    = y_lim[1],
    yend = y_lim[1] + sample_tick_len
  )

  theme_heat <- theme_bw(base_size = 8) +
    theme(
      panel.grid        = element_blank(),
      panel.border      = element_rect(fill = NA, colour = "black",
                                       linewidth = 0.4),
      axis.ticks        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(1.0, "mm"),
      axis.text         = element_text(colour = "black", size = 7),
      axis.title        = element_text(colour = "black", size = 8),
      plot.title        = element_text(size = 8, hjust = 0.5, face = "bold"),
      plot.margin       = margin(0, 2, 4, 2)
    )

  p_top <- ggplot(xspec) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = 0,    ymax = 1,
                  fill = Time), colour = NA) +
    scale_fill_gradientn(colours = palette, guide = "none") +
    scale_x_continuous(limits = x_lim, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(
      panel.border = element_rect(fill = NA, colour = "black",
                                  linewidth = 0.4),
      plot.margin  = margin(4, 2, 0, 2)
    )

  x_right <- x_lim[2]
  x_label <- x_right + gene_label_xoffset * diff(x_lim)

  gene_layers <- list()
  if (show_gene_annot && !is.null(gene_annot) && nrow(gene_annot) > 0) {

    gene_annot_plot <- gene_annot %>%
      arrange(desc(y)) %>%
      mutate(
        y_center = y,
        label    = gene_id
      )

    gene_seg <- gene_annot_plot %>%
      transmute(
        x    = x_right,
        xend = x_label - 0.01 * diff(x_lim),
        y    = y_center,
        yend = y_center
      )

    gene_layers <- list(
      geom_segment(
        data        = gene_seg,
        aes(x = x, xend = xend, y = y, yend = yend),
        colour      = gene_tick_col,
        linewidth   = gene_tick_lwd,
        inherit.aes = FALSE
      ),
      geom_text(
        data        = gene_annot_plot,
        aes(x = x_label, y = y_center, label = label),
        hjust       = 0,
        vjust       = 0.5,
        size        = gene_label_size,
        fontface    = "italic",
        colour      = gene_tick_col,
        inherit.aes = FALSE
      )
    )
  }

  legend_guide <- if (show_legend) {
    guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barheight      = unit(35, "mm"),
      barwidth       = unit(3,  "mm"),
      frame.colour   = "black",
      ticks.colour   = "black"
    )
  } else "none"

  p_heat <- ggplot(plotdata) +

    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = value),
              colour = NA) +

    geom_segment(
      data        = tick_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      colour      = sample_tick_col,
      linewidth   = sample_tick_lwd,
      inherit.aes = FALSE
    ) +

    gene_layers +

    scale_fill_gradientn(
      colours   = palette,
      limits    = z_range,
      oob       = scales::squish,
      na.value  = "grey80",       # NA基因显示灰色
      name      = value_label,
      guide     = legend_guide
    ) +

    scale_x_continuous(
      breaks = xb,
      labels = round(xb, 2),
      limits = x_lim,
      expand = c(0, 0)
    ) +

    scale_y_continuous(
      limits = y_lim,
      expand = c(0, 0),
      breaks = ybreaks$y,
      labels = ybreaks$target_label
    ) +

    labs(
      title = title,
      x     = "Time (h)",
      y     = "Relative position"
    ) +
    theme_heat +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y     = element_text(size = 7),
      legend.position = if (show_legend) "right" else "none",
      legend.title    = element_text(size = 7),
      legend.text     = element_text(size = 6),
      legend.margin   = margin(2, 2, 2, 2)
    )

  cowplot::plot_grid(
    p_top, p_heat,
    ncol        = 1,
    rel_heights = c(topbar_h, 1),
    align       = "v",
    axis        = "lr"
  )
}

# ============================================================
# ============================================================
mod_panel_3groups <- function(group_plots_list, groups,
                               shared_zrange, palette,
                               value_label, heat_h_mm, heat_w_mm,
                               gap_mm = 3) {


  cb_data <- data.frame(x = 0, y = seq(shared_zrange[1],
                                        shared_zrange[2],
                                        length.out = 100))
  p_cb <- ggplot(cb_data, aes(x = x, y = y, fill = y)) +
    geom_tile() +
    scale_fill_gradientn(
      colours  = palette,
      limits   = shared_zrange,
      name     = value_label,
      guide    = guide_colorbar(
        title.position = "top",
        title.hjust    = 0.5,
        barheight      = unit(40, "mm"),
        barwidth       = unit(4,  "mm"),
        frame.colour   = "black",
        ticks.colour   = "black"
      )
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title    = element_text(size = 7),
      legend.text     = element_text(size = 6)
    )

  cb_grob <- cowplot::get_legend(p_cb)

  panel <- cowplot::plot_grid(
    plotlist    = c(group_plots_list, list(cb_grob)),
    nrow        = 1,
    rel_widths  = c(rep(1, length(groups)), 0.25),
    align       = "h",
    axis        = "tb"
  )

  panel
}

# ============================================================
# ============================================================

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat <- mod_load(
  raw_csv      = RAW_TPM_CSV,
  expinfo_csv  = EXPINFO_CSV,
  genelist_csv = GENE_LIST_CSV,
  col_sampleid = COL_SAMPLEID,
  col_time     = COL_TIME,
  col_group    = COL_GROUP,
  col_treat    = COL_TREAT
)

expr_mat <- dat$expr_mat
anno_tbl <- dat$anno_tbl
info_use <- dat$info_use

message("\n== 构建Y轴坐标系 ==")
yspec_obj <- mod_yspec_rna(anno_tbl)

message(sprintf("  Y轴刻度验证：%d 个刻度点", nrow(yspec_obj$ybreaks)))

groups_avail <- intersect(GROUPS_ORDERED, unique(info_use$Group))
miss_groups  <- setdiff(GROUPS_ORDERED, groups_avail)
if (length(miss_groups) > 0)
  warning("以下组在数据中缺失：", paste(miss_groups, collapse = ", "))
message(sprintf("\n== 将处理 %d 个实验组：%s ==",
                length(groups_avail), paste(groups_avail, collapse = ", ")))

all_combos <- expand.grid(
  value_type = VALUE_TYPES,
  scale_mode = SCALE_MODES,
  stringsAsFactors = FALSE
)

for (ci in seq_len(nrow(all_combos))) {

  vtype <- all_combos$value_type[ci]
  smode <- all_combos$scale_mode[ci]
  combo_key <- paste(vtype, smode, sep = "_")

  message(sprintf("\n========== 模式：%s ==========", combo_key))

  group_results <- list()

  for (g in groups_avail) {
    message(sprintf("\n  ── 组：%s ──", g))

    g_info <- info_use %>%
      filter(Group == g) %>%
      arrange(Time)

    g_cols   <- g_info[[COL_SAMPLEID]]
    mat_g    <- expr_mat[, g_cols, drop = FALSE]

    if (vtype == "relTPM") {
      rel_res    <- mod_rel_tpm(mat_g, g_info, COL_SAMPLEID, PRE_LABEL)
      mat_val    <- rel_res$mat
      vtype_used <- rel_res$value_type
    } else {
      mat_val    <- mat_g
      vtype_used <- "TPM"
    }

    mat_scaled <- mod_scale(mat_val, mode = smode, value_type = vtype_used)

    xspec_g <- mod_xspec(g_info$Time)

    pd <- mod_plotdata(mat_scaled, g_info, COL_SAMPLEID,
                       xspec_g, yspec_obj)

    val_label <- switch(combo_key,
      "relTPM_log2"   = expression(log[2](rel.TPM)),
      "relTPM_zscore" = "z-score (rel.TPM)",
      "TPM_log2"      = expression(log[2](TPM + 1)),
      "TPM_zscore"    = "z-score (TPM)"
    )

    group_results[[g]] <- list(
      mat_scaled = mat_scaled,
      xspec      = xspec_g,
      pd         = pd,
      vtype_used = vtype_used,
      val_label  = val_label
    )
  }

  zrange_key <- combo_key
  manual_zr  <- ZRANGE_MANUAL[[zrange_key]]

  if (!is.null(manual_zr)) {
    zrange_shared <- manual_zr
    message(sprintf("  Colorbar范围（手动）: [%.2f, %.2f]",
                    zrange_shared[1], zrange_shared[2]))
  } else {
    if (ZRANGE_MODE == "shared") {
      all_vals <- unlist(lapply(group_results, function(r)
        r$mat_scaled[is.finite(r$mat_scaled)]))
      vtype_ref     <- group_results[[1]]$vtype_used
      zrange_shared <- mod_zrange(matrix(all_vals), vtype_ref, smode)
      message(sprintf("  Colorbar范围（共享自动）: [%.2f, %.2f]",
                      zrange_shared[1], zrange_shared[2]))
    }
  }

  plots_no_legend <- list()   # 供拼合面板用（无图例）

  for (g in groups_avail) {
    r  <- group_results[[g]]
    zr <- if (ZRANGE_MODE == "shared") zrange_shared else
      mod_zrange(r$mat_scaled, r$vtype_used, smode)

    p_single <- mod_plot_single(
      plotdata         = r$pd,
      xspec            = r$xspec,
      yspec_obj        = yspec_obj,
      z_range          = zr,
      palette          = PALETTE,
      n_xlabels        = N_XLABELS,
      title            = g,
      value_label      = r$val_label,
      topbar_h         = TOPBAR_H,
      sample_tick_len  = SAMPLE_TICK_LEN,
      sample_tick_lwd  = SAMPLE_TICK_LWD,
      sample_tick_col  = SAMPLE_TICK_COL,
      show_gene_annot  = TRUE,
      gene_label_size  = GENE_LABEL_SIZE,
      gene_tick_lwd    = GENE_TICK_LWD,
      gene_tick_col    = GENE_TICK_COL,
      gene_label_xoffset = GENE_LABEL_OFFSET_MM / HEAT_W_MM,
      show_legend      = TRUE
    )

    out_single <- file.path(OUT_DIR,
      sprintf("%s_%s.pdf", g, combo_key))
    ggsave(out_single, p_single,
           width = HEAT_W_MM, height = HEAT_H_MM,
           units = "mm", useDingbats = FALSE, limitsize = FALSE)
    message(sprintf("  ✓ 单组PDF: %s", basename(out_single)))

    plots_no_legend[[g]] <- mod_plot_single(
      plotdata         = r$pd,
      xspec            = r$xspec,
      yspec_obj        = yspec_obj,
      z_range          = zr,
      palette          = PALETTE,
      n_xlabels        = N_XLABELS,
      title            = g,
      value_label      = r$val_label,
      topbar_h         = TOPBAR_H,
      sample_tick_len  = SAMPLE_TICK_LEN,
      sample_tick_lwd  = SAMPLE_TICK_LWD,
      sample_tick_col  = SAMPLE_TICK_COL,
      show_gene_annot  = TRUE,
      gene_label_size  = GENE_LABEL_SIZE,
      gene_tick_lwd    = GENE_TICK_LWD,
      gene_tick_col    = GENE_TICK_COL,
      gene_label_xoffset = GENE_LABEL_OFFSET_MM / HEAT_W_MM,
      show_legend      = FALSE   # 拼合图不含图例
    )
  }

  val_label_str <- switch(combo_key,
    "relTPM_log2"   = "log2(rel.TPM)",
    "relTPM_zscore" = "z-score (rel.TPM)",
    "TPM_log2"      = "log2(TPM+1)",
    "TPM_zscore"    = "z-score (TPM)"
  )

  p_panel <- mod_panel_3groups(
    group_plots_list = plots_no_legend,
    groups           = groups_avail,
    shared_zrange    = zrange_shared,
    palette          = PALETTE,
    value_label      = val_label_str,
    heat_h_mm        = HEAT_H_MM,
    heat_w_mm        = HEAT_W_MM,
    gap_mm           = PANEL_GAP_MM
  )

  panel_w <- HEAT_W_MM * length(groups_avail) + 20   # +20mm给colorbar
  panel_h <- HEAT_H_MM

  out_panel <- file.path(OUT_DIR,
    sprintf("panel_3groups_%s.pdf", combo_key))
  ggsave(out_panel, p_panel,
         width = panel_w, height = panel_h,
         units = "mm", useDingbats = FALSE, limitsize = FALSE)
  message(sprintf("  ✓ 拼合面板: %s  (%dx%dmm)",
                  basename(out_panel), panel_w, panel_h))

  out_png <- sub("\\.pdf$", ".png", out_panel)
  ggsave(out_png, p_panel,
         width = panel_w, height = panel_h,
         units = "mm", dpi = 300, limitsize = FALSE)
}

message("\n== 输出参数日志 ==")

message("\n========== 全部完成 ==========")
message(sprintf("输出目录：%s", normalizePath(OUT_DIR, mustWork = FALSE)))
message(sprintf("单组PDF：%d 个（%d组 × %d种模式）",
                length(groups_avail) * nrow(all_combos),
                length(groups_avail), nrow(all_combos)))
message(sprintf("拼合面板：%d 个（每种模式1个）", nrow(all_combos)))
message("\n── 输出文件清单 ──")
message("  单组：{Group}_{value_type}_{scale_mode}.pdf")
message("  面板：panel_3groups_{value_type}_{scale_mode}.pdf/.png")



}


# Analysis section: DNA_dosage_correction
DNA_dosage_correction <- function() {


# =============================================================================
# dnaA Upshift Experiment: Genomic Copy Number Correction of Transcriptomics
# Version 2  (corrected from v1, 2026-06-25)
# =============================================================================
# Changes vs v1:
#            a clean, self-contained QC scatter using only qc_scatter2 logic.
#   [FIX 2] Column-name standardisation now stops with an informative error
#            instead of silently skipping if the rename map is incomplete.
#  [FIX 3]  Pre-shift baseline checks now use stop() instead of a
#            floor . Any gene with non-positive baseline
#            triggers an informative error, since the upstream low-expression
#            filter (max TPM < 10 AND min TPM < 1) should remove all such cases.
#   [FIX 4] Gene-count check after filtering now reports which genes are
#            missing before stopping, making debugging easier.
#   [FIX 5] Part IVb prepare_merge() now checks that at least one valid
#            data column survived the filter before proceeding.
#   [NOTE]  Biological logic unchanged: shared DNA copy-number profile
#           (averaged from Replicate1DNA + Replicate2DNA) applied uniformly to all three
#           across experiments .
# =============================================================================


library(tidyverse)
library(data.table)


path_dna_dpm   <- "dna_DPM.csv"
path_dna_info  <- "dna_metadata.csv"
path_rna_tpm   <- "rna_TPM.csv"
path_rna_info  <- "rna_metadata.csv"
path_gene_list <- "filtered_genes.csv"
dir_out        <- "Step4a_DNA_dosage_correction"
dir.create(dir_out, showWarnings = FALSE)


cat("Loading data...\n")
dna_dpm  <- as.data.frame(analysis_fread(path_dna_dpm))
dna_info <- analysis_read_tibble(path_dna_info, show_col_types = FALSE)
rna_tpm  <- as.data.frame(analysis_fread(path_rna_tpm))
rna_info <- analysis_read_tibble(path_rna_info, show_col_types = FALSE)
gene_list <- analysis_read_tibble(path_gene_list, show_col_types = FALSE)

# [FIX 2] Standardise rna_info first-column name
colnames(rna_info)[1] <- "Colnames"

# rna_tpm layout (confirmed from head(rna_tpm)):
#   col1 : gene      (character)
#   col2 : locus_tag (character)
#   col3 : product   (character)
#   col4 : Rel.Pos   (numeric, range -1 to 1)

#

# analysis_fread/analysis_read_csv adds "X" to column names that start with a digit.
# We only need to confirm that all sample columns exist in rna_info$Colnames.

META_COLS_TPM <- c("gene", "locus_tag", "product", "Rel.Pos")

# Guard: meta columns must actually be present with correct types
stopifnot(
  "rna_tpm must contain 'gene' column"      = "gene"      %in% colnames(rna_tpm),
  "rna_tpm must contain 'locus_tag' column" = "locus_tag" %in% colnames(rna_tpm),
  "rna_tpm must contain 'product' column"   = "product"   %in% colnames(rna_tpm),
  "rna_tpm must contain 'Rel.Pos' column"   = "Rel.Pos"   %in% colnames(rna_tpm)
)

sample_cols <- setdiff(colnames(rna_tpm), META_COLS_TPM)

# Check all sample columns are registered in rna_info
unregistered <- setdiff(sample_cols, rna_info$Colnames)
if (length(unregistered) > 0) {
  stop(sprintf(
    "%d rna_tpm sample column(s) not found in rna_info$Colnames:\n  %s\n",
    length(unregistered), paste(unregistered, collapse = ", ")),
    "Check that rna_info has one row per sample with a matching Colnames entry.")
}
cat(sprintf("  rna_tpm structure OK: %d genes, %d sample columns, all matched to rna_info.\n",
            nrow(rna_tpm), length(sample_cols)))


interp_safe <- function(t_obs, y_obs, t_target) {
  ord  <- order(t_obs); t_s <- t_obs[ord]; y_s <- y_obs[ord]
  t_lo <- t_s[1]; t_hi <- t_s[length(t_s)]
  result <- rep(NA_real_, length(t_target))
  for (i in seq_along(t_target)) {
    tt <- t_target[i]
    if (tt < t_lo || tt > t_hi) next
    if (tt %in% t_s) { result[i] <- y_s[t_s == tt][1]; next }
    j_hi <- which(t_s > tt)[1]; j_lo <- j_hi - 1L
    alpha <- (tt - t_s[j_lo]) / (t_s[j_hi] - t_s[j_lo])
    result[i] <- y_s[j_lo] + alpha * (y_s[j_hi] - y_s[j_lo])
  }
  result
}


# =============================================================================
# =============================================================================
cat("Processing DNA-seq data...\n")

dna_meta  <- dna_info %>%
  filter(Group %in% c("Replicate1DNA", "Replicate2DNA")) %>%
  select(Colnames, Group, time = Real_Uchanged_Time)
meta_Replicate1d <- dna_meta %>% filter(Group == "Replicate1DNA") %>% arrange(time)
meta_Replicate2d <- dna_meta %>% filter(Group == "Replicate2DNA") %>% arrange(time)
t_Replicate1d    <- meta_Replicate1d$time
t_Replicate2d    <- meta_Replicate2d$time

t_dna_union <- sort(union(t_Replicate1d, t_Replicate2d))
t_dna_safe  <- t_dna_union[
  t_dna_union >= max(min(t_Replicate1d), min(t_Replicate2d)) &
  t_dna_union <= min(max(t_Replicate1d), max(t_Replicate2d))]
t_dna_pre   <- t_dna_safe[t_dna_safe < 0]
t_dna_post  <- t_dna_safe[t_dna_safe >= 0]
cat(sprintf("  DNA safe common grid: %d points (%d pre, %d post)\n",
            length(t_dna_safe), length(t_dna_pre), length(t_dna_post)))

dpm_mat_Replicate1 <- as.matrix(dna_dpm[, meta_Replicate1d$Colnames])
dpm_mat_Replicate2 <- as.matrix(dna_dpm[, meta_Replicate2d$Colnames])
colnames(dpm_mat_Replicate1) <- as.character(t_Replicate1d)
colnames(dpm_mat_Replicate2) <- as.character(t_Replicate2d)

cat("  Interpolating Replicate1DNA...\n")
dpm_interp_Replicate1 <- t(apply(dpm_mat_Replicate1, 1, interp_safe,
                           t_obs = t_Replicate1d, t_target = t_dna_safe))
cat("  Interpolating Replicate2DNA...\n")
dpm_interp_Replicate2 <- t(apply(dpm_mat_Replicate2, 1, interp_safe,
                           t_obs = t_Replicate2d, t_target = t_dna_safe))
colnames(dpm_interp_Replicate1) <- colnames(dpm_interp_Replicate2) <- as.character(t_dna_safe)
stopifnot(!anyNA(dpm_interp_Replicate1), !anyNA(dpm_interp_Replicate2))

dpm_avg      <- (dpm_interp_Replicate1 + dpm_interp_Replicate2) / 2
pre_cols_dna <- as.character(t_dna_pre)
dpm_baseline <- rowMeans(dpm_avg[, pre_cols_dna, drop = FALSE])

if (any(dpm_baseline == 0))
  warning(sprintf("%d bins have zero DPM baseline — will produce Inf Rel.avgDPM",
                  sum(dpm_baseline == 0)))

rel_avg_dpm    <- sweep(dpm_avg, 1, dpm_baseline, FUN = "/")
rel_avg_dpm_df <- data.frame(
  genome_position   = dna_dpm$genome_position,
  Relative_position = dna_dpm$Relative_position,
  rel_avg_dpm)
colnames(rel_avg_dpm_df)[-(1:2)] <- paste0("t_", t_dna_safe)

fwrite(rel_avg_dpm_df, file.path(dir_out, "QC_rel_avg_DPM_500bins.csv"))
cat("  Saved: QC_rel_avg_DPM_500bins.csv\n")


# =============================================================================
# =============================================================================
cat("Mapping genes to genomic bins...\n")
gene_ids <- gene_list$gene

# [FIX 4] Informative error if gene IDs don't match
rna_filtered <- rna_tpm %>% filter(gene %in% gene_ids)
missing_genes <- setdiff(gene_ids, rna_filtered$gene)
if (length(missing_genes) > 0) {
  cat(sprintf("  Warning: %d gene(s) in gene_list not found in rna_tpm:\n",
              length(missing_genes)))
  cat("  ", paste(head(missing_genes, 10), collapse = ", "),
      if (length(missing_genes) > 10) "... (truncated)" else "", "\n")
}
stopifnot("Gene count mismatch after filtering" =
            nrow(rna_filtered) == length(gene_ids))

gene_relpos     <- rna_filtered$Rel.Pos
bin_relpos      <- rel_avg_dpm_df$Relative_position
nearest_bin_idx <- vapply(gene_relpos, function(rp)
  which.min(abs(bin_relpos - rp)), integer(1L))
max_map_err <- max(abs(gene_relpos - bin_relpos[nearest_bin_idx]))

gene_bin_mapping<-data.frame(gene=rna_filtered$gene,locus_tag=rna_filtered$locus_tag,
  gene_relative_position=gene_relpos,DNA_bin_index=nearest_bin_idx,
  DNA_bin_relative_position=bin_relpos[nearest_bin_idx])
fwrite(gene_bin_mapping,file.path(dir_out,"gene_to_DNA_bin_mapping.csv"))

cat(sprintf("  Max gene→bin mapping error: %.6f (bin width ≈ %.4f)\n",
            max_map_err, diff(range(bin_relpos)) / length(bin_relpos)))


# =============================================================================
# PART III: RNA-seq correction (per group)
# =============================================================================
cat("Correcting RNA-seq TPM values...\n")
rna_groups <- c("Replicate1RNAseq", "Replicate2RNAseq", "Replicate3RNAseq")
dna_t_min  <- min(t_dna_safe)
dna_t_max  <- max(t_dna_safe)
results    <- list()

for (grp in rna_groups) {
  cat(sprintf("  Processing %s...\n", grp))

  meta_grp <- rna_info %>%
    filter(Group == grp) %>% arrange(Real_Uchanged_Time) %>%
    select(Colnames, time = Real_Uchanged_Time)
  t_grp <- meta_grp$time
  t_pre <- t_grp[t_grp < 0]

  tpm_mat <- as.matrix(rna_filtered[, meta_grp$Colnames])
  colnames(tpm_mat) <- as.character(t_grp)
  rownames(tpm_mat) <- rna_filtered$gene

  t_correctable     <- t_grp[t_grp >= dna_t_min & t_grp <= dna_t_max]
  t_not_correctable <- t_grp[t_grp  < dna_t_min | t_grp  > dna_t_max]
  if (length(t_not_correctable) > 0)
    cat(sprintf("    Note: %d time point(s) outside DNA range, not corrected: %s\n",
                length(t_not_correctable),
                paste(t_not_correctable, collapse = ", ")))

  # Build CF matrix
  cf_mat <- matrix(NA_real_, nrow(rna_filtered), length(t_correctable),
                   dimnames = list(rna_filtered$gene, as.character(t_correctable)))
  for (g in seq_len(nrow(rna_filtered))) {
    dpm_ts <- as.numeric(rel_avg_dpm[nearest_bin_idx[g], ])
    cf_mat[g, ] <- interp_safe(t_dna_safe, dpm_ts, t_correctable)
  }
  n_bad_cf <- sum(!is.finite(cf_mat) | cf_mat <= 0, na.rm = TRUE)
  if (n_bad_cf > 0)
    warning(sprintf("  %d non-positive/non-finite CF values in %s", n_bad_cf, grp))

  # corrTPM
  tpm_raw_sub <- tpm_mat[, as.character(t_correctable), drop = FALSE]
  corrTPM_mat <- tpm_raw_sub / cf_mat
  correction_factors<-data.frame(gene=rna_filtered$gene,locus_tag=rna_filtered$locus_tag,
    cf_mat,check.names=FALSE)
  fwrite(correction_factors,file.path(dir_out,paste0("DNA_correction_factors_",sub("RNAseq","",grp),".csv")))

  if (length(t_not_correctable) > 0) {
    tpm_uncorr <- tpm_mat[, as.character(t_not_correctable), drop = FALSE]
    colnames(tpm_uncorr) <- paste0(t_not_correctable, "_rawTPM_no_CF")
    corrTPM_full <- cbind(corrTPM_mat, tpm_uncorr)
  } else {
    corrTPM_full <- corrTPM_mat
  }

  t_pre_corr <- t_pre[t_pre >= dna_t_min & t_pre <= dna_t_max]
  if (length(t_pre_corr) == 0) {
  stop(sprintf(
    "Group %s: no pre-shift time points fall within the DNA-seq time range [%g, %g]. ",
    grp, dna_t_min, dna_t_max),
    "All RNA-seq groups should have at least one pre-shift time point ",
    "covered by the DNA-seq grid. Check the experimental time alignment.")
}
   pre_baseline <- rowMeans(corrTPM_mat[, as.character(t_pre_corr), drop = FALSE])

  n_bad <- sum(!is.finite(pre_baseline) | pre_baseline <= 0, na.rm = TRUE)
  if (n_bad > 0) {
  bad_genes <- rna_filtered$gene[!is.finite(pre_baseline) | pre_baseline <= 0]
  stop(sprintf(
    "Group %s: %d gene(s) have non-positive or non-finite pre-shift baseline:\n  %s\n",
    grp, n_bad, paste(head(bad_genes, 10), collapse = ", ")),
    "Upstream low-expression filter (max TPM < 10 AND min TPM < 1) ",
    "should have removed these. Verify the gene list passed in.")
}

  rel_corrTPM_mat <- sweep(corrTPM_mat, 1, pre_baseline, FUN = "/")

  gene_meta_cols <- rna_filtered %>% select(gene, locus_tag, product, Rel.Pos)
  out_corrTPM    <- bind_cols(gene_meta_cols, as.data.frame(corrTPM_full))
  out_relcorrTPM <- bind_cols(gene_meta_cols, as.data.frame(rel_corrTPM_mat))

  grp_tag <- sub("RNAseq", "", grp)
  fwrite(out_corrTPM,    file.path(dir_out, sprintf("corrTPM_%s.csv",     grp_tag)))
  fwrite(out_relcorrTPM, file.path(dir_out, sprintf("Rel_corrTPM_%s.csv", grp_tag)))
  cat(sprintf("    Saved: corrTPM_%s.csv | Rel_corrTPM_%s.csv\n", grp_tag, grp_tag))

  results[[grp]] <- list(meta = meta_grp, t_correctable = t_correctable,
                          corrTPM = corrTPM_mat, rel_corrTPM = rel_corrTPM_mat,
                          pre_baseline = pre_baseline)
}


# =============================================================================
# PART IV: QC outputs
# =============================================================================
cat("Generating QC outputs...\n")

# IV.1 Rel.avgDPM genome profile
qc_times <- intersect(c("-5", "56", "126", "220"), colnames(rel_avg_dpm_df))
if (length(qc_times) > 0) {
  qc_dpm_long <- rel_avg_dpm_df %>%
    select(Relative_position, all_of(paste0("t_", qc_times))) %>%
    pivot_longer(-Relative_position, names_to = "time_label", values_to = "Rel_avgDPM") %>%
    mutate(time_min = as.integer(sub("t_", "", time_label)))
  p_dpm <- ggplot(qc_dpm_long,
                  aes(Relative_position, Rel_avgDPM, colour = factor(time_min))) +
    geom_line(linewidth = 0.6) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
    scale_colour_brewer(palette = "Set1", name = "Time (h)") +
    labs(title = "Rel.avgDPM genomic profile",
         x = "Relative genome position (oriC = 0)", y = "Rel.avgDPM") +
    theme_bw(base_size = 11)
  ggsave(file.path(dir_out, "QC_rel_avgDPM_profile.pdf"), p_dpm, width = 7, height = 4)
  cat("  Saved: QC_rel_avgDPM_profile.pdf\n")
}

# PART IVb: Merge group outputs
# =============================================================================
cat("Merging group outputs...\n")

prepare_merge <- function(file_path, group_prefix) {
  df <- as.data.frame(analysis_fread(file_path))
  meta_c <- c("gene", "locus_tag", "product", "Rel.Pos")
  keep   <- setdiff(colnames(df), meta_c)
  keep   <- keep[!grepl("_rawTPM_no_CF", keep)]
  keep   <- keep[!is.na(suppressWarnings(as.numeric(keep)))]
  # [FIX 5] Abort if no valid data columns remain
  if (length(keep) == 0)
    stop(sprintf("prepare_merge: no valid numeric time-point columns in %s", file_path))
  df_sub <- df[, c(meta_c, keep)]
  colnames(df_sub)[colnames(df_sub) %in% keep] <- paste0(group_prefix, "_", keep)
  df_sub
}

grp_tags <- c("Replicate1", "Replicate2", "Replicate3")
corr_files <- setNames(
  file.path(dir_out, paste0("corrTPM_",     grp_tags, ".csv")), grp_tags)
rel_files  <- setNames(
  file.path(dir_out, paste0("Rel_corrTPM_", grp_tags, ".csv")), grp_tags)

merge_list <- function(flist) {
  lst <- lapply(names(flist), function(g) {
    cat(sprintf("  Loading %s...\n", flist[g]))
    prepare_merge(flist[g], g)
  })
  Reduce(function(a, b)
    full_join(a, b, by = c("gene", "locus_tag", "product", "Rel.Pos")),
    lst)
}

merged_corr <- merge_list(corr_files)
merged_rel  <- merge_list(rel_files)

fwrite(merged_corr, file.path(dir_out, "corrTPM_ALL_groups_merged.csv"))
fwrite(merged_rel,  file.path(dir_out, "Rel_corrTPM_ALL_groups_merged.csv"))
cat(sprintf("  Saved: corrTPM_ALL_groups_merged.csv       (%d genes × %d cols)\n",
            nrow(merged_corr), ncol(merged_corr)))
cat(sprintf("  Saved: Rel_corrTPM_ALL_groups_merged.csv  (%d genes × %d cols)\n",
            nrow(merged_rel), ncol(merged_rel)))


# =============================================================================
# PART V: Session info
# =============================================================================

cat("Done. All outputs written to:", dir_out, "\n")

# =============================================================================
# Output files (unchanged from v1 except the QC scatter no longer crashes)
# =============================================================================
# output_genomic_correction/
# =============================================================================

}


# Analysis section: Corrected_RNA_heatmaps
Corrected_RNA_heatmaps <- function() {


# ============================================================
#
#
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(readr)
})

# ============================================================
# ============================================================

CORRTPM_CSV   <- "corrected_TPM.csv"  # corrTPM 合并表
GENE_LIST_CSV <- "filtered_genes.csv"  # 2918基因列表

GROUPS_ORDERED <- c("Replicate1", "Replicate2", "Replicate3")

VALUE_TYPES <- c("relTPM", "TPM")

SCALE_MODES <- c("log2", "zscore")

ZRANGE_MODE <- "shared"   # "shared" 或 "per_group"

ZRANGE_MANUAL <- list(
  relTPM_log2   = NULL,
  relTPM_zscore = c(-3, 3),
  TPM_log2      = NULL,
  TPM_zscore    = c(-3, 3)
)

PALETTE <- c("#3B4CC0", "#88A1D2", "#F7F7F7", "#F3A982", "#B40426")


GENE_LABEL_OFFSET_MM <- 2.5
GENE_LABEL_SIZE      <- 2.2
GENE_TICK_LWD        <- 0.3
GENE_TICK_COL        <- "black"

N_XLABELS <- 8

SAMPLE_TICK_LEN <- 50
SAMPLE_TICK_LWD <- 0.4
SAMPLE_TICK_COL <- "black"

HEAT_W_MM    <- 150
HEAT_H_MM    <- 130
TOPBAR_H     <- 0.04
PANEL_GAP_MM <- 3


OUT_DIR <- "Step4b_corrected_RNA_heatmaps"

# ============================================================
# ============================================================
mod_load_corrtpm <- function(corrtpm_csv, genelist_csv,
                             groups_ordered) {

  message("== 读入 corrTPM 数据 ==")

  mat_raw <- analysis_read_csv(corrtpm_csv, header = TRUE, check.names = FALSE,
                      stringsAsFactors = FALSE)

  annot_cols <- c("gene", "locus_tag", "product", "Rel.Pos")
  stopifnot(all(annot_cols %in% colnames(mat_raw)))

  gene_list <- analysis_read_csv(genelist_csv, header = TRUE, check.names = FALSE,
                        stringsAsFactors = FALSE)
  if (!"gene" %in% colnames(gene_list)) {
    if ("Gene" %in% colnames(gene_list)) gene_list$gene <- gene_list$Gene
    else stop("基因列表中找不到 'gene' 或 'Gene' 列")
  }
  genes_keep <- unique(gene_list$gene)
  message(sprintf("  基因列表：%d 个唯一基因名", length(genes_keep)))

  mat_filt <- mat_raw %>% filter(gene %in% genes_keep)
  message(sprintf("  corrTPM 过滤后：%d 行（目标 %d）",
                  nrow(mat_filt), length(genes_keep)))
  if (nrow(mat_filt) != length(genes_keep)) {
    miss <- setdiff(genes_keep, mat_filt$gene)
    warning(sprintf("  %d 个基因在 corrTPM 中未找到：%s",
                    length(miss), paste(head(miss, 5), collapse = ", ")))
  }

  all_data_cols <- setdiff(colnames(mat_filt), annot_cols)
  col_info <- data.frame(
    col   = all_data_cols,
    stringsAsFactors = FALSE
  )
  parsed <- regmatches(all_data_cols,
                       regexpr("^(Replicate1|Replicate2|Replicate3)_(-?\\d+)$",
                               all_data_cols, perl = TRUE))
  valid_mask <- nchar(parsed) > 0
  if (sum(!valid_mask) > 0)
    warning(sprintf("  %d 列无法解析为 {group}_{time}，已忽略：%s",
                    sum(!valid_mask),
                    paste(all_data_cols[!valid_mask], collapse = ", ")))

  valid_cols <- all_data_cols[valid_mask]
  parts <- strsplit(valid_cols, "_")
  col_info <- tibble(
  col   = valid_cols,
  group = sapply(parts, `[`, 1),
  time  = as.numeric(sapply(parts, `[`, 2)) / 60   # ← 转成小时
)

  for (g in sort(unique(col_info$group))) {
    ci_g <- col_info %>% filter(group == g) %>% arrange(time)
    pre_n  <- sum(ci_g$time < 0)
    post_n <- sum(ci_g$time >= 0)
    message(sprintf("  组 %s: %d 个时间点（pre<0：%d，post>=0：%d）  时间：%s",
                    g, nrow(ci_g), pre_n, post_n,
                    paste(ci_g$time, collapse = ", ")))
  }

  anno_tbl <- mat_filt %>%
    select(all_of(annot_cols)) %>%
    mutate(Rel.Pos = as.numeric(Rel.Pos))

  expr_mat <- as.matrix(mat_filt[, valid_cols, drop = FALSE])
  rownames(expr_mat) <- mat_filt$gene
  storage.mode(expr_mat) <- "numeric"

  list(
    expr_mat = expr_mat,
    anno_tbl = anno_tbl,
    col_info = col_info,          # 替代原 info_use
    genes    = mat_filt$gene
  )
}

# ============================================================
# ============================================================
mod_yspec_rna <- function(anno_tbl, target_genes = NULL) {

  anno_ord <- anno_tbl %>%
    arrange(desc(Rel.Pos)) %>%
    mutate(
      gene_id = gene,
      y    = seq(n(), 1),
      ymin = y - 0.5,
      ymax = y + 0.5
    )

  n_genes <- nrow(anno_ord)
  message(sprintf("  Y轴：%d 基因，Rel.Pos +1 在顶部(y=%d)，-1 在底部(y=1)",
                  n_genes, n_genes))

  target_rp <- c(1.0, 0.5, 0.0, -0.5, -1.0)
  yspec_breaks <- lapply(seq_along(target_rp), function(i) {
    rp <- target_rp[i]
    anno_ord %>%
      mutate(dist = abs(Rel.Pos - rp)) %>%
      slice_min(dist, n = 1, with_ties = FALSE) %>%
      mutate(target_label = rp)
  }) %>%
    bind_rows() %>%
    distinct(target_label, .keep_all = TRUE) %>%
    arrange(desc(y))

  gene_annot <- NULL
  if (!is.null(target_genes) && length(target_genes) > 0) {
    gene_annot <- anno_ord %>%
      filter(gene_id %in% target_genes) %>%
      select(gene_id, Rel.Pos, y, ymin, ymax) %>%
      arrange(desc(y))
    found <- gene_annot$gene_id
    miss  <- setdiff(target_genes, found)
    message(sprintf("  目标基因：%d/%d 个找到",
                    length(found), length(target_genes)))
    if (length(miss) > 0)
      message(sprintf("  未找到：%s", paste(miss, collapse = ", ")))
  }

  list(
    yspec      = anno_ord,
    ybreaks    = yspec_breaks,
    gene_annot = gene_annot,
    n_genes    = n_genes
  )
}

# ============================================================
# ============================================================
mod_xspec <- function(time_vec) {
  x <- sort(unique(as.numeric(time_vec)))
  x <- x[is.finite(x)]
  if (length(x) == 1)
    return(tibble(Time = x, xmin = x - 0.5, xmax = x + 0.5))
  mid  <- (x[-1] + x[-length(x)]) / 2
  xmin <- c(x[1] - (mid[1] - x[1]), mid)
  xmax <- c(mid, x[length(x)] + (x[length(x)] - mid[length(mid)]))
  tibble(Time = x, xmin = xmin, xmax = xmax)
}

# ============================================================
# ============================================================
mod_rel_tpm_corr <- function(mat_group, pre_cols, eps = 1e-8) {

  if (length(pre_cols) == 0) {
    message("  [警告] 本组无 Pre（time<0）样品，直接使用原始 corrTPM")
    return(list(mat = mat_group, value_type = "TPM", n_degraded = 0))
  }

  pre_mean <- rowMeans(mat_group[, pre_cols, drop = FALSE], na.rm = TRUE)

  invalid <- !is.finite(pre_mean) | pre_mean <= eps
  n_degraded <- sum(invalid)
  if (n_degraded > 0)
    message(sprintf("  [注意] %d 个基因 Pre 均值≤0，rel.TPM 设为 NA", n_degraded))

  pre_mean[invalid] <- NA
  mat_rel <- sweep(mat_group, 1, pre_mean, "/")
  mat_rel[!is.finite(mat_rel)] <- NA

  list(mat = mat_rel, value_type = "relTPM", n_degraded = n_degraded)
}

# ============================================================
# ============================================================
mod_scale <- function(mat, mode, value_type) {

  if (mode == "log2") {
    pseudo <- if (value_type == "relTPM") 0 else 1
    mat <- log2(mat + pseudo)
    mat[!is.finite(mat)] <- NA

  } else if (mode == "zscore") {
    mu  <- rowMeans(mat, na.rm = TRUE)
    sdv <- apply(mat, 1, sd, na.rm = TRUE)
    sdv[!is.finite(sdv) | sdv < 1e-10] <- 1
    mat <- sweep(mat, 1, mu,  "-")
    mat <- sweep(mat, 1, sdv, "/")
    mat[!is.finite(mat)] <- 0
  }

  mat
}

# ============================================================
# ============================================================
mod_zrange <- function(mat, value_type, scale_mode) {

  if (scale_mode == "zscore") return(c(-3, 3))

  vals <- mat[is.finite(mat)]
  if (length(vals) == 0) return(c(-3, 3))

  if (value_type == "relTPM") {
    lim <- quantile(abs(vals), 0.99, na.rm = TRUE)
    lim <- ceiling(lim * 4) / 4
    return(c(-lim, lim))
  }

  lo <- floor(quantile(vals, 0.01, na.rm = TRUE) * 4) / 4
  hi <- ceiling(quantile(vals, 0.99, na.rm = TRUE) * 4) / 4
  c(lo, hi)
}

# ============================================================
# ============================================================
mod_plotdata_corr <- function(mat_scaled, col_info_g, xspec, yspec_obj) {

  yspec <- yspec_obj$yspec
  cols  <- col_info_g$col

  mat_sub <- mat_scaled[, cols, drop = FALSE]

  as.data.frame(mat_sub, check.names = FALSE) %>%
    rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id,
                 names_to  = "col",
                 values_to = "value") %>%
    left_join(col_info_g %>% select(col, time) %>% rename(Time = time),
              by = "col") %>%
    left_join(xspec, by = "Time") %>%
    left_join(yspec %>% select(gene_id, y, ymin, ymax),
              by = "gene_id") %>%
    filter(is.finite(xmin), is.finite(ymin))
}

# ============================================================
# ============================================================
mod_plot_single <- function(plotdata, xspec, yspec_obj,
                            z_range, palette, n_xlabels,
                            title       = NULL,
                            value_label = "value",
                            topbar_h    = 0.03,
                            sample_tick_len = 50,
                            sample_tick_lwd = 0.4,
                            sample_tick_col = "black",
                            show_gene_annot  = TRUE,
                            gene_label_size  = 2.2,
                            gene_tick_lwd    = 0.3,
                            gene_tick_col    = "black",
                            gene_label_xoffset = 0.05,
                            show_legend = TRUE) {

  yspec      <- yspec_obj$yspec
  ybreaks    <- yspec_obj$ybreaks
  gene_annot <- yspec_obj$gene_annot

  x_lim <- range(c(xspec$xmin, xspec$xmax), na.rm = TRUE)
  y_lim <- range(c(yspec$ymin, yspec$ymax), na.rm = TRUE)

  xb <- breaks_pretty(n = n_xlabels)(range(xspec$Time))
  xb <- xb[is.finite(xb) & xb >= x_lim[1] & xb <= x_lim[2]]

  sample_times <- sort(unique(xspec$Time))
  tick_df <- tibble(
    x    = sample_times,
    xend = sample_times,
    y    = y_lim[1],
    yend = y_lim[1] + sample_tick_len
  )

  theme_heat <- theme_bw(base_size = 8) +
    theme(
      panel.grid        = element_blank(),
      panel.border      = element_rect(fill = NA, colour = "black",
                                       linewidth = 0.4),
      axis.ticks        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(1.0, "mm"),
      axis.text         = element_text(colour = "black", size = 7),
      axis.title        = element_text(colour = "black", size = 8),
      plot.title        = element_text(size = 8, hjust = 0.5, face = "bold"),
      plot.margin       = margin(0, 2, 4, 2)
    )

  p_top <- ggplot(xspec) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = 0,    ymax = 1,
                  fill = Time), colour = NA) +
    scale_fill_gradientn(colours = palette, guide = "none") +
    scale_x_continuous(limits = x_lim, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(
      panel.border = element_rect(fill = NA, colour = "black",
                                  linewidth = 0.4),
      plot.margin  = margin(4, 2, 0, 2)
    )

  x_right <- x_lim[2]
  x_label <- x_right + gene_label_xoffset * diff(x_lim)

  gene_layers <- list()
  if (show_gene_annot && !is.null(gene_annot) && nrow(gene_annot) > 0) {

    gene_annot_plot <- gene_annot %>%
      arrange(desc(y)) %>%
      mutate(y_center = y, label = gene_id)

    gene_seg <- gene_annot_plot %>%
      transmute(
        x    = x_right,
        xend = x_label - 0.01 * diff(x_lim),
        y    = y_center,
        yend = y_center
      )

    gene_layers <- list(
      geom_segment(
        data        = gene_seg,
        aes(x = x, xend = xend, y = y, yend = yend),
        colour      = gene_tick_col,
        linewidth   = gene_tick_lwd,
        inherit.aes = FALSE
      ),
      geom_text(
        data        = gene_annot_plot,
        aes(x = x_label, y = y_center, label = label),
        hjust       = 0,
        vjust       = 0.5,
        size        = gene_label_size,
        fontface    = "italic",
        colour      = gene_tick_col,
        inherit.aes = FALSE
      )
    )
  }

  legend_guide <- if (show_legend) {
    guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barheight      = unit(35, "mm"),
      barwidth       = unit(3,  "mm"),
      frame.colour   = "black",
      ticks.colour   = "black"
    )
  } else "none"

  p_heat <- ggplot(plotdata) +

    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = value),
              colour = NA) +

    geom_segment(
      data        = tick_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      colour      = sample_tick_col,
      linewidth   = sample_tick_lwd,
      inherit.aes = FALSE
    ) +

    gene_layers +

    scale_fill_gradientn(
      colours   = palette,
      limits    = z_range,
      oob       = scales::squish,
      na.value  = "grey80",
      name      = value_label,
      guide     = legend_guide
    ) +

    scale_x_continuous(
      breaks = xb,
      labels = round(xb, 2),
      limits = x_lim,
      expand = c(0, 0)
    ) +

    scale_y_continuous(
      limits = y_lim,
      expand = c(0, 0),
      breaks = ybreaks$y,
      labels = ybreaks$target_label
    ) +

    labs(
      title = title,
      x     = "Time (h)",   # corrTPM 时间单位为分钟
      y     = "Relative position"
    ) +
    theme_heat +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y     = element_text(size = 7),
      legend.position = if (show_legend) "right" else "none",
      legend.title    = element_text(size = 7),
      legend.text     = element_text(size = 6),
      legend.margin   = margin(2, 2, 2, 2)
    )

  cowplot::plot_grid(
    p_top, p_heat,
    ncol        = 1,
    rel_heights = c(topbar_h, 1),
    align       = "v",
    axis        = "lr"
  )
}

# ============================================================
# ============================================================
mod_panel_3groups <- function(group_plots_list, groups,
                               shared_zrange, palette,
                               value_label, heat_h_mm, heat_w_mm,
                               gap_mm = 3) {

  cb_data <- data.frame(x = 0, y = seq(shared_zrange[1],
                                        shared_zrange[2],
                                        length.out = 100))
  p_cb <- ggplot(cb_data, aes(x = x, y = y, fill = y)) +
    geom_tile() +
    scale_fill_gradientn(
      colours  = palette,
      limits   = shared_zrange,
      name     = value_label,
      guide    = guide_colorbar(
        title.position = "top",
        title.hjust    = 0.5,
        barheight      = unit(40, "mm"),
        barwidth       = unit(4,  "mm"),
        frame.colour   = "black",
        ticks.colour   = "black"
      )
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title    = element_text(size = 7),
      legend.text     = element_text(size = 6)
    )

  cb_grob <- cowplot::get_legend(p_cb)

  panel <- cowplot::plot_grid(
    plotlist    = c(group_plots_list, list(cb_grob)),
    nrow        = 1,
    rel_widths  = c(rep(1, length(groups)), 0.25),
    align       = "h",
    axis        = "tb"
  )

  panel
}

# ============================================================
# ============================================================

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat <- mod_load_corrtpm(
  corrtpm_csv   = CORRTPM_CSV,
  genelist_csv  = GENE_LIST_CSV,
  groups_ordered = GROUPS_ORDERED
)

expr_mat <- dat$expr_mat
anno_tbl <- dat$anno_tbl
col_info <- dat$col_info   # tibble: col, group, time

message("\n== 构建Y轴坐标系 ==")
yspec_obj <- mod_yspec_rna(anno_tbl)

groups_avail <- intersect(GROUPS_ORDERED, unique(col_info$group))
miss_groups  <- setdiff(GROUPS_ORDERED, groups_avail)
if (length(miss_groups) > 0)
  warning("以下组在数据中缺失：", paste(miss_groups, collapse = ", "))
message(sprintf("\n== 将处理 %d 个实验组：%s ==",
                length(groups_avail), paste(groups_avail, collapse = ", ")))

all_combos <- expand.grid(
  value_type = VALUE_TYPES,
  scale_mode = SCALE_MODES,
  stringsAsFactors = FALSE
)

for (ci_idx in seq_len(nrow(all_combos))) {

  vtype <- all_combos$value_type[ci_idx]
  smode <- all_combos$scale_mode[ci_idx]
  combo_key <- paste(vtype, smode, sep = "_")

  message(sprintf("\n========== 模式：%s ==========", combo_key))

  group_results <- list()

  for (g in groups_avail) {
    message(sprintf("\n  ── 组：%s ──", g))

    ci_g     <- col_info %>% filter(group == g) %>% arrange(time)
    g_cols   <- ci_g$col
    pre_cols <- ci_g %>% filter(time < 0) %>% pull(col)
    post_cols <- ci_g %>% filter(time >= 0) %>% pull(col)
    message(sprintf("  Pre 列（time<0）：%d 个，Post 列（time>=0）：%d 个",
                    length(pre_cols), length(post_cols)))

    mat_g <- expr_mat[, g_cols, drop = FALSE]

    if (vtype == "relTPM") {
      rel_res    <- mod_rel_tpm_corr(mat_g, pre_cols)
      mat_val    <- rel_res$mat
      vtype_used <- rel_res$value_type
    } else {
      mat_val    <- mat_g
      vtype_used <- "TPM"
    }

    mat_scaled <- mod_scale(mat_val, mode = smode, value_type = vtype_used)

    xspec_g <- mod_xspec(ci_g$time)

    pd <- mod_plotdata_corr(mat_scaled, ci_g, xspec_g, yspec_obj)

    val_label <- switch(combo_key,
      "relTPM_log2"   = expression(log[2](rel.corrTPM)),
      "relTPM_zscore" = "z-score (rel.corrTPM)",
      "TPM_log2"      = expression(log[2](corrTPM + 1)),
      "TPM_zscore"    = "z-score (corrTPM)"
    )

    group_results[[g]] <- list(
      mat_scaled = mat_scaled,
      xspec      = xspec_g,
      pd         = pd,
      vtype_used = vtype_used,
      val_label  = val_label
    )
  }

  zrange_key <- combo_key
  manual_zr  <- ZRANGE_MANUAL[[zrange_key]]

  if (!is.null(manual_zr)) {
    zrange_shared <- manual_zr
    message(sprintf("  Colorbar范围（手动）: [%.2f, %.2f]",
                    zrange_shared[1], zrange_shared[2]))
  } else {
    if (ZRANGE_MODE == "shared") {
      all_vals <- unlist(lapply(group_results, function(r)
        r$mat_scaled[is.finite(r$mat_scaled)]))
      vtype_ref     <- group_results[[1]]$vtype_used
      zrange_shared <- mod_zrange(matrix(all_vals), vtype_ref, smode)
      message(sprintf("  Colorbar范围（共享自动）: [%.2f, %.2f]",
                      zrange_shared[1], zrange_shared[2]))
    }
  }

  plots_no_legend <- list()

  for (g in groups_avail) {
    r  <- group_results[[g]]
    zr <- if (ZRANGE_MODE == "shared") zrange_shared else
      mod_zrange(r$mat_scaled, r$vtype_used, smode)

    p_single <- mod_plot_single(
      plotdata           = r$pd,
      xspec              = r$xspec,
      yspec_obj          = yspec_obj,
      z_range            = zr,
      palette            = PALETTE,
      n_xlabels          = N_XLABELS,
      title              = paste0(g, " (corrTPM)"),
      value_label        = r$val_label,
      topbar_h           = TOPBAR_H,
      sample_tick_len    = SAMPLE_TICK_LEN,
      sample_tick_lwd    = SAMPLE_TICK_LWD,
      sample_tick_col    = SAMPLE_TICK_COL,
      show_gene_annot    = TRUE,
      gene_label_size    = GENE_LABEL_SIZE,
      gene_tick_lwd      = GENE_TICK_LWD,
      gene_tick_col      = GENE_TICK_COL,
      gene_label_xoffset = GENE_LABEL_OFFSET_MM / HEAT_W_MM,
      show_legend        = TRUE
    )

    out_single <- file.path(OUT_DIR,
                            sprintf("%s_%s_corrTPM.pdf", g, combo_key))
    ggsave(out_single, p_single,
           width = HEAT_W_MM, height = HEAT_H_MM,
           units = "mm", useDingbats = FALSE, limitsize = FALSE)
    message(sprintf("  ✓ 单组PDF: %s", basename(out_single)))

    plots_no_legend[[g]] <- mod_plot_single(
      plotdata           = r$pd,
      xspec              = r$xspec,
      yspec_obj          = yspec_obj,
      z_range            = zr,
      palette            = PALETTE,
      n_xlabels          = N_XLABELS,
      title              = paste0(g, " (corrTPM)"),
      value_label        = r$val_label,
      topbar_h           = TOPBAR_H,
      sample_tick_len    = SAMPLE_TICK_LEN,
      sample_tick_lwd    = SAMPLE_TICK_LWD,
      sample_tick_col    = SAMPLE_TICK_COL,
      show_gene_annot    = TRUE,
      gene_label_size    = GENE_LABEL_SIZE,
      gene_tick_lwd      = GENE_TICK_LWD,
      gene_tick_col      = GENE_TICK_COL,
      gene_label_xoffset = GENE_LABEL_OFFSET_MM / HEAT_W_MM,
      show_legend        = FALSE
    )
  }

  val_label_str <- switch(combo_key,
    "relTPM_log2"   = "log2(rel.corrTPM)",
    "relTPM_zscore" = "z-score (rel.corrTPM)",
    "TPM_log2"      = "log2(corrTPM+1)",
    "TPM_zscore"    = "z-score (corrTPM)"
  )

  p_panel <- mod_panel_3groups(
    group_plots_list = plots_no_legend,
    groups           = groups_avail,
    shared_zrange    = zrange_shared,
    palette          = PALETTE,
    value_label      = val_label_str,
    heat_h_mm        = HEAT_H_MM,
    heat_w_mm        = HEAT_W_MM,
    gap_mm           = PANEL_GAP_MM
  )

  panel_w <- HEAT_W_MM * length(groups_avail) + 20
  panel_h <- HEAT_H_MM

  out_panel <- file.path(OUT_DIR,
    sprintf("panel_3groups_%s_corrTPM.pdf", combo_key))
  ggsave(out_panel, p_panel,
         width = panel_w, height = panel_h,
         units = "mm", useDingbats = FALSE, limitsize = FALSE)
  message(sprintf("  ✓ 拼合面板: %s  (%dx%dmm)",
                  basename(out_panel), panel_w, panel_h))

  out_png <- sub("\\.pdf$", ".png", out_panel)
  ggsave(out_png, p_panel,
         width = panel_w, height = panel_h,
         units = "mm", dpi = 300, limitsize = FALSE)
}

message("\n== 输出参数日志 ==")

message("\n========== 全部完成 ==========")
message(sprintf("输出目录：%s", normalizePath(OUT_DIR, mustWork = FALSE)))
message(sprintf("单组PDF：%d 个（%d组 × %d种模式）",
                length(groups_avail) * nrow(all_combos),
                length(groups_avail), nrow(all_combos)))
message(sprintf("拼合面板：%d 个（每种模式1个）", nrow(all_combos)))
message("\n── 输出文件命名规则 ──")
message("  单组：{Group}_{value_type}_{scale_mode}_corrTPM.pdf")
message("  面板：panel_3groups_{value_type}_{scale_mode}_corrTPM.pdf/.png")


}


# Analysis section: Response_identification
Response_identification <- function() {


# =============================================================================
#
#
# =============================================================================


suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(cowplot)
  library(ggrepel)
})



PATH_CORRTPM <- "corrected_TPM.csv"
DIR_OUT      <- "Step5_response_identification"
DIR_PLOTS    <- file.path(DIR_OUT, "plots")
DIR_TABLES   <- file.path(DIR_OUT, "tables")

FDR_THR <- cfg$dnaa$fdr_threshold    # 各组独立 BH-FDR 阈值
FC_THR  <- cfg$dnaa$fc_threshold_log2   # max|log2FC| 阈值（对应 1.5 倍，log2(1.5) = 0.585）


set.seed(42)
dir.create(DIR_PLOTS,  recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_TABLES, recursive = TRUE, showWarnings = FALSE)


theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid        = element_blank(),
      panel.border      = element_rect(fill = NA, colour = "black", linewidth = 0.8),
      axis.ticks        = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = unit(1.2, "mm"),
      plot.title        = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      plot.subtitle     = element_text(size = base_size - 1, hjust = 0.5, colour = "grey40"),
      legend.position   = "right",
      aspect.ratio      = 1
    )
}

PAL_GROUP <- c("Replicate1" = "#1BB09A", "Replicate2" = "#6A7FDB", "Replicate3" = "#F0A202")
PAL_DIR   <- c("Up" = "#D6604D", "Down" = "#4393C3",
               "Inconsistent_Up" = "#F4A582", "Inconsistent_Down" = "#92C5DE",
               "Inconsistent" = "grey70")


cat("Loading data...\n")
df <- as.data.frame(analysis_fread(PATH_CORRTPM))

META_COLS <- c("gene", "locus_tag", "product", "Rel.Pos")
data_cols <- setdiff(colnames(df), META_COLS)

parse_cols <- function(cols) {
  m     <- regmatches(cols, regexpr("^(Replicate1|Replicate2|Replicate3)_(-?\\d+)$", cols, perl = TRUE))
  valid <- nchar(m) > 0
  parts <- strsplit(m[valid], "_")
  tibble(
    col   = cols[valid],
    group = sapply(parts, `[`, 1),
    time  = as.integer(sapply(parts, `[`, 2))
  )
}

col_info <- parse_cols(data_cols)
groups   <- sort(unique(col_info$group))   # "Replicate3" "Replicate1" "Replicate2"

cat(sprintf("  %d genes  x  %d expression columns\n", nrow(df), length(data_cols)))
for (g in groups) {
  ci <- col_info %>% filter(group == g)
  pre_n  <- sum(ci$time < 0)
  post_n <- sum(ci$time >= 0)
  cat(sprintf("  Group %s: %d pre (%s)  |  %d post\n",
              g, pre_n,
              paste(sort(ci$time[ci$time < 0]), collapse = ","),
              post_n))
}


cat("\nStep 1-2: Computing log2FC series and per-group t-tests...\n")

compute_group_stats <- function(grp) {
  ci_g      <- col_info %>% filter(group == grp)
  pre_cols  <- ci_g %>% filter(time < 0)  %>% arrange(time) %>% pull(col)
  post_cols <- ci_g %>% filter(time >= 0) %>% arrange(time) %>% pull(col)
  n_post    <- length(post_cols)

  pre_mat  <- as.matrix(df[, pre_cols,  drop = FALSE])
  post_mat <- as.matrix(df[, post_cols, drop = FALSE])

  pre_mean <- rowMeans(pre_mat)    # 各基因的 pre 均值（已确保 > 0）

  l2fc_mat <- log2(post_mat / pre_mean)  # nGenes x n_post

  l2fc_mean    <- rowMeans(l2fc_mat)
  l2fc_sd      <- apply(l2fc_mat, 1, sd)
  l2fc_max_abs <- apply(abs(l2fc_mat), 1, max)

  max_idx        <- apply(abs(l2fc_mat), 1, which.max)
  l2fc_at_max    <- l2fc_mat[cbind(seq_len(nrow(l2fc_mat)), max_idx)]

  se    <- l2fc_sd / sqrt(n_post)
  t_val <- l2fc_mean / se
  p_val <- 2 * pt(-abs(t_val), df = n_post - 1)
  p_val <- pmax(p_val, .Machine$double.eps)   # 避免 p=0

  data.frame(
    gene        = df$gene,
    group       = grp,
    pre_mean    = pre_mean,
    l2fc_mean   = l2fc_mean,
    l2fc_sd     = l2fc_sd,
    l2fc_max_abs = l2fc_max_abs,
    l2fc_at_max  = l2fc_at_max,
    t_stat      = t_val,
    p_raw       = p_val,
    n_post      = n_post,
    stringsAsFactors = FALSE
  )
}

stats_list <- lapply(groups, compute_group_stats)
stats_df   <- bind_rows(stats_list)

cat(sprintf("  Computed stats for %d gene × group combinations\n", nrow(stats_df)))


cat("\nStep 3: Per-group BH-FDR correction...\n")

stats_df <- stats_df %>%
  group_by(group) %>%
  mutate(fdr = p.adjust(p_raw, method = "BH")) %>%
  ungroup()

for (g in groups) {
  n_sig <- stats_df %>% filter(group == g, fdr < FDR_THR) %>% nrow()
  cat(sprintf("  Group %s: %d genes with FDR < %.2f\n", g, n_sig, FDR_THR))
}


cat("\nStep 4: Dual filter (FDR + FC) per group...\n")

stats_df <- stats_df %>%
  mutate(PASS = fdr < FDR_THR & l2fc_max_abs >= FC_THR)

pass_per_group <- stats_df %>%
  group_by(group) %>%
  summarise(n_pass = sum(PASS), .groups = "drop")

for (i in seq_len(nrow(pass_per_group))) {
  cat(sprintf("  Group %s: %d genes pass FDR < %.2f AND max|log2FC| >= %.3f (%.1fx)\n",
              pass_per_group$group[i],
              pass_per_group$n_pass[i],
              FDR_THR, FC_THR, 2^FC_THR))
}


cat("\nStep 5: Three-group strict intersection...\n")

pass_sets <- stats_df %>%
  filter(PASS) %>%
  group_by(group) %>%
  summarise(genes = list(gene), .groups = "drop")

gene_sets  <- setNames(pass_sets$genes, pass_sets$group)
final_genes <- Reduce(intersect, gene_sets)

cat(sprintf("  Intersection of all 3 groups: %d genes\n", length(final_genes)))

n_12 <- length(intersect(gene_sets[["Replicate1"]], gene_sets[["Replicate2"]]))
n_13 <- length(intersect(gene_sets[["Replicate1"]], gene_sets[["Replicate3"]]))
n_23 <- length(intersect(gene_sets[["Replicate2"]], gene_sets[["Replicate3"]]))
cat(sprintf("  Replicate1 ∩ Replicate2: %d  |  Replicate1 ∩ Replicate3: %d  |  Replicate2 ∩ Replicate3: %d\n",
            n_12, n_13, n_23))


cat("\nStep 6: Direction annotation...\n")

direction_df_full <- stats_df %>%
  filter(gene %in% final_genes) %>%
  group_by(gene) %>%
  summarise(
    n_up   = sum(l2fc_mean > 0),
    n_down = sum(l2fc_mean < 0),
    .groups = "drop"
  ) %>%
  mutate(
    direction = case_when(
      n_up   == 3 ~ "Up",
      n_down == 3 ~ "Down",
      n_up   == 2 ~ "Inconsistent_Up",
      n_down == 2 ~ "Inconsistent_Down",
      TRUE        ~ "Inconsistent"
    )
  )

cat("  Direction breakdown (before sign-consistency filter):\n")
direction_df_full %>% count(direction) %>% print()

direction_df <- direction_df_full %>%
  filter(direction %in% c("Up", "Down"))

n_excluded_inconsistent <- nrow(direction_df_full) - nrow(direction_df)
cat(sprintf("  Excluded %d genes with inconsistent direction across replicates\n",
            n_excluded_inconsistent))

final_genes <- direction_df$gene
cat(sprintf("  Final responsive genes (direction-consistent): %d\n", length(final_genes)))


cat("\nStep 7: Assembling output tables...\n")

stats_wide <- stats_df %>%
  filter(gene %in% final_genes) %>%
  select(gene, group, pre_mean, l2fc_mean, l2fc_sd,
         l2fc_max_abs, l2fc_at_max, t_stat, p_raw, fdr, PASS) %>%
  pivot_wider(
    names_from  = group,
    values_from = c(pre_mean, l2fc_mean, l2fc_sd,
                    l2fc_max_abs, l2fc_at_max, t_stat, p_raw, fdr, PASS),
    names_glue  = "{group}_{.value}"
  )

meta_df <- df %>% select(gene, locus_tag, product, Rel.Pos)

responsive_genes <- stats_wide %>%
  left_join(direction_df %>% select(gene, direction, n_up, n_down),
            by = "gene") %>%
  left_join(meta_df, by = "gene") %>%
  mutate(
    median_l2fc_mean  = pmap_dbl(
      select(., ends_with("_l2fc_mean")),
      ~ median(c(...))),
    max_l2fc_max_abs  = pmap_dbl(
      select(., ends_with("_l2fc_max_abs")),
      ~ max(c(...))),
    max_fdr           = pmap_dbl(
      select(., ends_with("_fdr")),
      ~ max(c(...)))
  ) %>%
  arrange(direction, desc(max_l2fc_max_abs))

all_stats_wide <- stats_df %>%
  select(gene, group, l2fc_mean, l2fc_max_abs, p_raw, fdr, PASS) %>%
  pivot_wider(
    names_from  = group,
    values_from = c(l2fc_mean, l2fc_max_abs, p_raw, fdr, PASS),
    names_glue  = "{group}_{.value}"
  ) %>%
  mutate(
    PASS_all3 = `Replicate1_PASS` & `Replicate2_PASS` & `Replicate3_PASS`,
    median_l2fc_mean = pmap_dbl(
      select(., ends_with("_l2fc_mean")),
      ~ median(c(...))),
    max_l2fc_max_abs = pmap_dbl(
      select(., ends_with("_l2fc_max_abs")),
      ~ max(c(...)))
  ) %>%
 left_join(direction_df_full %>% select(gene, direction), by = "gene") %>%
  mutate(direction = replace_na(direction, "Not_responsive"))




fwrite(responsive_genes,  file.path(DIR_TABLES, "responsive_genes.csv"))
fwrite(all_stats_wide,    file.path(DIR_TABLES, "all_genes_stats.csv"))

param_summary <- tibble(
  parameter = c(
    "Input file", "Total genes (after pre-filter)", "FDR_THR", "FC_THR (log2)",
    "FC_THR (linear)", "Filter strategy",
    paste0("Group ", groups, " pass (FDR+FC)"),
    "Intersection of 3 groups (before direction filter)",
    "Direction — Up (3/3 concordant positive)",
    "Direction — Down (3/3 concordant negative)",
    "Direction — Inconsistent_Up (2 positive, 1 negative)",
    "Direction — Inconsistent_Down (2 negative, 1 positive)",
    "Direction — Inconsistent (other; excluded)",
    "Final responsive genes (Up + Down, retained)",
    "Excluded due to inconsistent direction"
  ),
  value = as.character(c(
    PATH_CORRTPM, nrow(df), FDR_THR, FC_THR,
    round(2^FC_THR, 3), "Per-group BH-FDR + FC, strict 3-group intersection",
    pass_per_group$n_pass,
    nrow(direction_df_full),
    sum(direction_df_full$direction == "Up"),
    sum(direction_df_full$direction == "Down"),
    sum(direction_df_full$direction == "Inconsistent_Up"),
    sum(direction_df_full$direction == "Inconsistent_Down"),
    sum(direction_df_full$direction == "Inconsistent"),
    length(final_genes),
    n_excluded_inconsistent
  ))
)
fwrite(param_summary, file.path(DIR_TABLES, "parameter_summary.csv"))

cat(sprintf("\n  Saved: responsive_genes.csv       (%d genes)\n", nrow(responsive_genes)))
cat(sprintf("  Saved: all_genes_stats.csv         (%d genes)\n", nrow(all_stats_wide)))
cat("  Saved: parameter_summary.csv\n")


# =============================================================================
# =============================================================================
cat("\nGenerating plots...\n")


p1_data <- stats_df %>%
  mutate(neg_log10_fdr = -log10(fdr),
         PASS_fc = l2fc_max_abs >= FC_THR)

p1 <- ggplot(p1_data, aes(x = neg_log10_fdr, colour = group, fill = group)) +
  geom_density(alpha = 0.15, linewidth = 0.7) +
  geom_vline(xintercept = -log10(FDR_THR),
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  annotate("text",
           x = -log10(FDR_THR) + 0.3, y = Inf,
           label = sprintf("FDR = %.2f", FDR_THR),
           hjust = 0, vjust = 1.5, size = 3, colour = "grey30") +
  scale_colour_manual(values = PAL_GROUP, name = "Experiment") +
  scale_fill_manual(values   = PAL_GROUP, name = "Experiment") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "Per-group FDR distribution",
    subtitle = sprintf("Single-sample t-test on log2FC series, BH correction per group"),
    x = expression(-log[10](FDR)),
    y = "Density"
  ) +
  theme_pub() +
  facet_wrap(~group, ncol = 1, scales = "free_y") +
  theme(aspect.ratio = NULL, strip.background = element_rect(fill = "grey95"))

ggsave(file.path(DIR_PLOTS, "01_fdr_distribution.pdf"),
       p1, width = 5, height = 7)
cat("  Saved: 01_fdr_distribution.pdf\n")


p2_data <- stats_df %>%
  filter(fdr < FDR_THR) %>%
  mutate(PASS_fc = l2fc_max_abs >= FC_THR)

p2 <- ggplot(p2_data, aes(x = l2fc_max_abs, colour = group, fill = group)) +
  geom_density(alpha = 0.15, linewidth = 0.7) +
  geom_vline(xintercept = FC_THR,
             linetype = "dashed", colour = "grey30", linewidth = 0.7) +
  annotate("text",
           x = FC_THR + 0.05, y = Inf,
           label = sprintf("%.3f (%.1fx)", FC_THR, 2^FC_THR),
           hjust = 0, vjust = 1.5, size = 3, colour = "grey30") +
  scale_colour_manual(values = PAL_GROUP, name = "Experiment") +
  scale_fill_manual(values   = PAL_GROUP, name = "Experiment") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "max|log2FC| distribution (FDR-passing genes)",
    x = "max |log2FC|",
    y = "Density"
  ) +
  theme_pub()

ggsave(file.path(DIR_PLOTS, "02_max_log2fc_distribution.pdf"),
       p2, width = 6, height = 5)
cat("  Saved: 02_max_log2fc_distribution.pdf\n")


volcano_data <- stats_df %>%
  mutate(
    neg_log10_fdr = pmin(-log10(fdr), 40),
    color_cat = case_when(
      PASS  ~ "Pass",
      TRUE  ~ "Background"
    ),
    color_cat = factor(color_cat, levels = c("Background", "Pass"))
  )

pal_volcano <- c(
  "Background"            = "grey80",
  "Pass" = "#4393C3"
)

p3 <- ggplot(volcano_data,
             aes(x = l2fc_mean, y = neg_log10_fdr, colour = color_cat)) +
  geom_point(
    data = volcano_data %>% filter(color_cat == "Background"),
    size = 0.4, alpha = 0.3
  ) +
  geom_point(
    data = volcano_data %>% filter(color_cat == "Pass"),
    size = 0.8, alpha = 0.6
  ) +
  geom_hline(yintercept = -log10(FDR_THR),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = c(-FC_THR, FC_THR),
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_colour_manual(values = pal_volcano, name = NULL) +
  facet_wrap(~group, ncol = 3) +
  labs(
    title    = "Volcano plot — per group",
    subtitle = sprintf("Dashed lines: FDR = %.2f, |log2FC| = %.3f (%.1fx)",
                       FDR_THR, FC_THR, 2^FC_THR),
    x        = "mean log2FC (post vs pre)",
    y        = expression(-log[10](FDR))
  ) +
  theme_pub() +
  theme(aspect.ratio = NULL,
        strip.background = element_rect(fill = "grey95"),
        legend.position = "bottom")

ggsave(file.path(DIR_PLOTS, "03_volcano_plot.pdf"),
       p3, width = 12, height = 5)
cat("  Saved: 03_volcano_plot.pdf\n")




dir_count <- direction_df_full %>%
  count(direction) %>%
  mutate(
    direction = factor(direction,
                       levels = c("Up", "Down",
                                  "Inconsistent_Up", "Inconsistent_Down",
                                  "Inconsistent")),
    pct = n / sum(n) * 100,
    label = sprintf("%d\n(%.0f%%)", n, pct),
    retained = ifelse(direction %in% c("Up", "Down"), "Retained", "Excluded")
  )

p5 <- ggplot(dir_count, aes(x = direction, y = n, fill = direction)) +
  geom_col(aes(alpha = retained), width = 0.6, colour = "black", linewidth = 0.3) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 3.5, colour = "white", fontface = "bold") +
  scale_fill_manual(values = PAL_DIR, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Direction composition of responsive genes",
    subtitle = sprintf("n = %d (intersection); %d retained after sign-consistency filter",
                   nrow(direction_df_full), length(final_genes)),
    x = "Direction", y = "Number of genes"
  ) +
  theme_pub() + theme(aspect.ratio = NULL)

ggsave(file.path(DIR_PLOTS, "05_direction_composition.pdf"),
       p5, width = 5, height = 4.5)
cat("  Saved: 05_direction_composition.pdf\n")


cat("  Running sensitivity analysis...\n")

calc_sensitivity <- function(fdr_thr, fc_thr) {
  pass_s <- stats_df %>%
    filter(fdr < fdr_thr, l2fc_max_abs >= fc_thr)
  sets_s <- split(pass_s$gene, pass_s$group)
  if (!all(groups %in% names(sets_s))) return(tibble(n_genes = 0L))
  inter  <- Reduce(intersect, sets_s[groups])
  tibble(n_genes = length(inter))
}

sens_fdr <- map_dfr(
  c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2),
  ~ calc_sensitivity(.x, FC_THR) %>%
    mutate(param = "FDR_threshold", value = .x)
)

sens_fc <- map_dfr(
  seq(0.263, 1.322, by = 0.1),
  ~ calc_sensitivity(FDR_THR, .x) %>%
    mutate(param = "FC_threshold_log2", value = .x)
)

sens_all <- bind_rows(sens_fdr, sens_fc)
fwrite(sens_all, file.path(DIR_TABLES, "sensitivity_analysis.csv"))

make_sens_plot <- function(df_s, param_name, x_lab,
                           vline_val, x_trans = "identity") {
  df_sub <- df_s %>% filter(param == param_name)

  p_n <- ggplot(df_sub, aes(x = value, y = n_genes)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    geom_vline(xintercept = vline_val,
               linetype = "dashed", colour = "firebrick", linewidth = 0.6) +
    scale_x_continuous(trans = x_trans) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = x_lab, y = "# responsive genes") +
    theme_pub() + theme(aspect.ratio = NULL)

  p_n
}

ps_fdr <- make_sens_plot(sens_all, "FDR_threshold",
                         "FDR threshold", FDR_THR, x_trans = "log10")
ps_fc  <- make_sens_plot(sens_all, "FC_threshold_log2",
                         "max|log2FC| threshold", FC_THR)

p6 <- plot_grid(
  ps_fdr + ggtitle("FDR scan"),
  ps_fc + ggtitle("FC scan"),
  ncol = 2, align = "hv"
)

ggsave(file.path(DIR_PLOTS, "06_sensitivity_analysis.pdf"),
       p6, width = 9, height = 4)
cat("  Saved: 06_sensitivity_analysis.pdf\n")


venn_data <- tibble(
  comparison              = c("Replicate1 only", "Replicate2 only", "Replicate3 only",
                               "Replicate1 ∩ Replicate2 only", "Replicate1 ∩ Replicate3 only", "Replicate2 ∩ Replicate3 only",
                               "Replicate1 ∩ Replicate2 ∩ Replicate3"),
  n_genes = c(
    length(setdiff(gene_sets[["Replicate1"]], union(gene_sets[["Replicate2"]], gene_sets[["Replicate3"]]))),
    length(setdiff(gene_sets[["Replicate2"]], union(gene_sets[["Replicate1"]], gene_sets[["Replicate3"]]))),
    length(setdiff(gene_sets[["Replicate3"]], union(gene_sets[["Replicate1"]], gene_sets[["Replicate2"]]))),
    length(setdiff(intersect(gene_sets[["Replicate1"]], gene_sets[["Replicate2"]]), gene_sets[["Replicate3"]])),
    length(setdiff(intersect(gene_sets[["Replicate1"]], gene_sets[["Replicate3"]]), gene_sets[["Replicate2"]])),
    length(setdiff(intersect(gene_sets[["Replicate2"]], gene_sets[["Replicate3"]]), gene_sets[["Replicate1"]])),
    length(final_genes)
  )
)

fwrite(venn_data, file.path(DIR_TABLES, "venn_diagram_data.csv"))

cat("  Saved: venn_diagram_data.csv (use for Venn diagram in GraphPad or ggvenn)\n")
cat("\n  Venn diagram numbers:\n")
print(venn_data, n = Inf)


cat(sprintf("\n%s\n", strrep("=", 52)))
cat(sprintf("  Input genes:                    %d\n", nrow(df)))
cat(sprintf("  FDR threshold:                  %.3f\n", FDR_THR))
cat(sprintf("  FC threshold:                   %.3f (%.2fx)\n", FC_THR, 2^FC_THR))
cat(sprintf("  Filter strategy:                Per-group, strict intersection\n"))
cat(sprintf("  Final responsive genes:         %d\n", length(final_genes)))
cat(sprintf("    -- Up (3/3 positive):              %d\n",
            sum(direction_df_full$direction == "Up")))
cat(sprintf("    -- Down (3/3 negative):            %d\n",
            sum(direction_df_full$direction == "Down")))
cat(sprintf("    -- Inconsistent_Up (excluded):     %d\n",
            sum(direction_df_full$direction == "Inconsistent_Up")))
cat(sprintf("    -- Inconsistent_Down (excluded):   %d\n",
            sum(direction_df_full$direction == "Inconsistent_Down")))
cat(sprintf("    -- Inconsistent (excluded):        %d\n",
            sum(direction_df_full$direction == "Inconsistent")))
cat(sprintf("%s\n", strrep("=", 52)))


cat(sprintf("\nDone. All outputs saved to: %s\n", normalizePath(DIR_OUT)))


# =============================================================================
# =============================================================================
# output_v3/


# =============================================================================


}


# Analysis section: Kmeans_clustering_and_diagnostics
Kmeans_clustering_and_diagnostics <- function() {


## ================================================================
##
##
##    Step6_clustering/
## ================================================================


## ================================================================
## ================================================================


PATH_TPM   <- "rna_TPM.csv"
PATH_INFO  <- "rna_metadata.csv"
PATH_RES   <- "response_all_stats.csv"

DIR_OUT   <- "Step6_clustering"
DIR_TBL   <- file.path(DIR_OUT, "tables")
DIR_DIAG  <- file.path(DIR_OUT, "diagnostics")
DIR_PERCK <- file.path(DIR_DIAG, "per_cluster")
DIR_HEAT  <- file.path(DIR_OUT, "heatmap")
DIR_LINE  <- file.path(DIR_OUT, "classavg")

for (d in c(DIR_TBL, DIR_DIAG, DIR_PERCK, DIR_HEAT, DIR_LINE)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}


## ================================================================
## ================================================================

K_FINAL        <- cfg$dnaa$k_final
K_SCAN         <- 2:10
KMEANS_NSTART  <- cfg$dnaa$nstart
KMEANS_ITERMAX <- cfg$dnaa$iter_max
SEED           <- cfg$dnaa$seed

PRE_LABEL  <- "Pre_1aTc"
TIME_COL   <- "Real_Uchanged_Time"  # 原始单位 min，读入后除以 60 转换为 h
EPS        <- 1e-8

PAL_CLUSTER <- c(
  "1" = rgb(241, 102, 103, maxColorValue = 255),
  "2" = rgb(119, 148, 229, maxColorValue = 255),
  "3" = rgb(129, 164, 115, maxColorValue = 255),
  "4" = rgb(229, 168,  75, maxColorValue = 255)
)

PAL_GROUP <- c(
  "Replicate1RNAseq" = "#1BB09A",
  "Replicate2RNAseq" = "#6A7FDB",
  "Replicate3RNAseq" = "#F0A202"
)

PAL_DIR <- c("Up" = "#D6604D", "Down" = "#4393C3")

PALETTE_HEAT <- colorRampPalette(
  c("#3B4CC0", "#88A1D2", "#F7F7F7", "#F3A982", "#B40426")
)(200)
Z_RANGE <- c(-3, 3)    # z-score 截断范围

CLUSTER_ORDER <- c(2, 1, 4, 3)   # 按需修改数字顺序

CLUSTER_GAP_ROWS <- 2      # cluster 间空白行数
CLUSTER_LW       <- 0.35   # cluster 边框线宽（pt）
RUG_NPC          <- 0.01   # rug 标记长度（NPC）
MARK_XPAD        <- 0.06   # 右侧标注留白（相对数据区宽度）
MARK_SEG         <- 0.02   # 连接线长度比例
MARK_TXT         <- 0.03   # 文字偏移比例
CH_ROW_GAP       <- 1.0    # ComplexHeatmap cluster 间距（mm）

BASE_SIZE    <- 7      # 基础字号（pt）
LINE_W       <- 0.6    # 折线宽度
PT_SIZE      <- 0.8    # 数据点大小
RIBBON_ALPHA <- 0.20   # 95% CI 色带透明度
AXIS_LW      <- 0.35   # 坐标轴线宽
PANEL_INCH   <- 1.2    # panel 物理宽高（英寸，正方形）
LEGEND_W     <- 1.1    # 图例列宽（英寸）
TOTAL_H      <- 2.2    # 含标题画布总高度（英寸）
PAD_W        <- 0.35   # PDF 左右 padding
PAD_H        <- 0.35   # PDF 上下 padding

SHADE_X1   <- 0
SHADE_X2   <- 0          # 根据实际实验窗口设置，单位 h
SHADE_FILL <- "grey90"

COMBINED_ORDER <- c("Replicate1RNAseq", "Replicate2RNAseq", "Replicate3RNAseq")
GROUP_GAP_X <- 30 / 60   # 单位 h  # 非等距合并热图组间空白（h）


## ================================================================
## ================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(cluster)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(cowplot)
  library(grid)
  library(scales)
})

set.seed(SEED)

col_fun_z <- colorRamp2(
  seq(Z_RANGE[1], Z_RANGE[2], length.out = length(PALETTE_HEAT)),
  PALETTE_HEAT
)


## ================================================================
## ================================================================

theme_diag <- function(base_size = BASE_SIZE) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(hjust = 0.5, colour = "grey40"),
      panel.grid       = element_blank(),
      panel.border     = element_rect(fill = NA, colour = "black",
                                      linewidth = 0.8),
      strip.background = element_rect(fill = "grey95", colour = "grey80"),
      legend.position  = "right"
    )
}

theme_graphpad <- function(base_size = BASE_SIZE) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text             = element_text(color = "black"),
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      hjust = 0),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size, color = "black"),
      axis.ticks       = element_line(color = "black", linewidth = AXIS_LW),
      axis.line        = element_line(color = "black", linewidth = AXIS_LW),
      panel.border     = element_rect(fill = NA, color = "black",
                                      linewidth = AXIS_LW),
      panel.grid       = element_blank(),
      legend.position  = "right",
      legend.title     = element_text(size = base_size),
      legend.text      = element_text(size = base_size),
      legend.key.height = unit(0.85, "lines"),
      legend.key.width  = unit(1.0,  "lines"),
      plot.margin      = margin(4, 6, 4, 4, unit = "pt")
    )
}

theme_heatmap <- function(base_size = BASE_SIZE) {
  theme_graphpad(base_size) +
    theme(
      panel.border = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin  = margin(5.5, 35, 5.5, 8)
    )
}


## ================================================================
## ================================================================

set_panel_size <- function(p, w = PANEL_INCH, h = PANEL_INCH) {
  gt  <- ggplotGrob(p)
  pid <- which(grepl("^panel", gt$layout$name))
  for (k in pid) {
    gt$widths[gt$layout$l[k]]  <- unit(w, "in")
    gt$heights[gt$layout$t[k]] <- unit(h, "in")
  }
  gt
}

no_clip <- function(gt) { gt$layout$clip <- "off"; gt }

save_pdf <- function(path, grob, w, h) {
  if (isTRUE(capabilities("cairo"))) {
    cairo_pdf(path, width = w + PAD_W, height = h + PAD_H, family = "sans")
  } else {
    pdf(path, width = w + PAD_W, height = h + PAD_H, useDingbats = FALSE)
  }
  grid.newpage(); grid.draw(grob); dev.off()
}

make_tilespec <- function(x) {
  x <- sort(x[is.finite(x)])
  if (anyDuplicated(x)) {
    for (i in seq_along(x)) {
      dup <- sum(x[seq_len(i - 1)] == x[i])
      if (dup > 0) x[i] <- x[i] + dup * 1e-6
    }
  }
  if (length(x) == 1)
    return(tibble(xmin = x - 0.5, xmax = x + 0.5))
  mid  <- (x[-1] + x[-length(x)]) / 2
  tibble(
    xmin = c(x[1] - (mid[1] - x[1]), mid),
    xmax = c(mid, x[length(x)] + (x[length(x)] - mid[length(mid)]))
  )
}

make_axis_params <- function(times, step_s = 0.5, step_l = 1.0,
                              thresh = 2.0, pad = 1) {
  span <- diff(range(times))
  step <- if (span <= thresh) step_s else step_l
  lo   <- floor(min(times)  / step) * step - pad * step
  hi   <- ceiling(max(times) / step) * step + pad * step
  brks <- seq(lo, hi, by = step)
  xs   <- make_tilespec(times)
  list(limits  = c(lo, hi),
       breaks  = brks,
       labels  = as.character(round(brks, 1)),
       x_left  = min(xs$xmin),
       x_right = max(xs$xmax))
}


scale_fill_z <- function() {
  scale_fill_gradientn(
    colours = PALETTE_HEAT,
    limits  = Z_RANGE,
    breaks  = seq(Z_RANGE[1], Z_RANGE[2], by = 1),
    oob     = scales::squish,
    name    = "z-score",
    guide   = guide_colorbar(frame.colour = "black",
                             barheight    = unit(30, "mm"))
  )
}


## ================================================================
## ================================================================

cat("\n========== 读入数据 ==========\n")

df_tpm <- analysis_fread(PATH_TPM) %>% as.data.frame()
colnames(df_tpm) <- gsub("^\ufeff", "", colnames(df_tpm))
GENE_COL    <- if ("gene" %in% colnames(df_tpm)) "gene" else "Gene"
META_COLS   <- intersect(c("gene","Gene","locus_tag","product","Rel.Pos"),
                         colnames(df_tpm))
sample_cols <- setdiff(colnames(df_tpm), META_COLS)
cat(sprintf("  TPM: %d 基因 × %d 样本\n", nrow(df_tpm), length(sample_cols)))

df_info <- analysis_fread(PATH_INFO) %>% as.data.frame()
colnames(df_info) <- gsub("^\ufeff", "", colnames(df_info))
df_info <- df_info[, !duplicated(colnames(df_info)), drop = FALSE]
stopifnot(all(c("Colnames","Group","Treatment",TIME_COL) %in% colnames(df_info)))

n_match <- length(intersect(df_info$Colnames, sample_cols))
if (n_match < nrow(df_info))
  stop("ExpInfo$Colnames 与 TPM 列名不完全匹配，请检查输入文件。")
cat(sprintf("  ExpInfo: %d 行（样本名全部匹配）\n", nrow(df_info)))

df_res <- analysis_fread(PATH_RES) %>% as.data.frame()
colnames(df_res) <- gsub("^\ufeff", "", colnames(df_res))
df_488 <- df_res %>% filter(direction %in% c("Up","Down"))
genes_488 <- intersect(df_488$gene, df_tpm[[GENE_COL]])
cat(sprintf("  响应基因: %d（Up %d / Down %d）\n",
            length(genes_488),
            sum(df_488$direction == "Up"),
            sum(df_488$direction == "Down")))


groups <- sort(unique(df_info$Group))
cat(sprintf("  实验组: %s\n", paste(groups, collapse = ", ")))
df_info[[TIME_COL]] <- df_info[[TIME_COL]] / 60

## ================================================================
##
## ================================================================

cat("\n========== 构建 z-score 矩阵 ==========\n")

build_zscore_one_group <- function(df_tpm, info_g, gene_set) {
  info_g     <- info_g %>% arrange(.data[[TIME_COL]])
  sample_use <- intersect(info_g$Colnames, colnames(df_tpm))
  info_g     <- info_g %>% filter(Colnames %in% sample_use)

  mat <- df_tpm %>%
    filter(.data[[GENE_COL]] %in% gene_set) %>%
    select(all_of(c(GENE_COL, info_g$Colnames)))
  genes   <- mat[[GENE_COL]]
  mat_num <- as.matrix(mat[, -1, drop = FALSE])
  rownames(mat_num) <- genes
  storage.mode(mat_num) <- "numeric"

  z_mat <- t(scale(t(mat_num), center = TRUE, scale = TRUE))
  z_mat[!is.finite(z_mat)] <- 0                        # 方差=0 的行置 0
  colnames(z_mat) <- paste0(info_g$Group[1], "|", info_g$Colnames)
  list(z = z_mat, info = info_g)
}

z_list    <- list()
info_list <- list()

for (g in groups) {
  info_g       <- df_info %>% filter(Group == g)
  out_g        <- build_zscore_one_group(df_tpm, info_g, genes_488)
  z_list[[g]]  <- out_g$z
  info_list[[g]] <- out_g$info
  cat(sprintf("  %s: %d 样本\n", g, ncol(z_list[[g]])))
}

all_genes <- sort(Reduce(union, lapply(z_list, rownames)))
Z_mat <- do.call(cbind, lapply(z_list, function(m) {
  mm <- matrix(0, nrow = length(all_genes), ncol = ncol(m),
               dimnames = list(all_genes, colnames(m)))
  mm[rownames(m), ] <- m
  mm
}))
cat(sprintf("\n  Z_mat: %d 基因 × %d 样本\n", nrow(Z_mat), ncol(Z_mat)))

fwrite(data.frame(Gene = rownames(Z_mat), Z_mat, check.names = FALSE),
       file.path(DIR_TBL, "input_matrix_zscore_tpm.csv"))


## ================================================================
## ================================================================

cat("\n========== k 扫描 ==========\n")

d_mat <- dist(Z_mat)
metrics_list <- lapply(K_SCAN, function(k) {
  set.seed(SEED)
  km    <- kmeans(Z_mat, centers = k,
                  nstart = KMEANS_NSTART, iter.max = KMEANS_ITERMAX)
  sizes <- table(km$cluster)
  sil   <- if (k >= 2) mean(silhouette(km$cluster, d_mat)[,"sil_width"])
           else NA_real_
  tibble(k = k, wss = km$tot.withinss, silhouette = sil,
         max_size_frac = max(sizes) / sum(sizes), min_size = min(sizes))
})
metrics_df <- bind_rows(metrics_list)
fwrite(metrics_df, file.path(DIR_TBL, "k_scan_metrics.csv"))

make_scan_panel <- function(ycol, ylab, log_y = FALSE, ttl = "") {
  p <- ggplot(metrics_df, aes(x = k, y = .data[[ycol]])) +
    geom_line(colour = "#4393C3", linewidth = 0.7) +
    geom_point(colour = "#4393C3", size = 2) +
    geom_vline(xintercept = K_FINAL, linetype = "dashed",
               colour = "firebrick", linewidth = 0.5) +
    scale_x_continuous(breaks = K_SCAN) +
    labs(title = ttl, x = "k", y = ylab) + theme_diag()
  if (log_y) p + scale_y_log10() else p
}

p_kscan <- plot_grid(
  make_scan_panel("wss",           "Total within-cluster SS",
                  log_y = TRUE, ttl = "Elbow (WSS)"),
  make_scan_panel("silhouette",    "Average silhouette width",
                  ttl = "Silhouette"),
  make_scan_panel("max_size_frac", "Largest cluster / N",
                  ttl = "Size balance"),
  make_scan_panel("min_size",      "Smallest cluster size",
                  ttl = "Smallest cluster"),
  ncol = 2, align = "hv"
)
ggsave(file.path(DIR_DIAG, "01_k_scan.pdf"), p_kscan, width = 8, height = 7)
cat("  Saved: 01_k_scan.pdf\n")


## ================================================================
## ================================================================

cat(sprintf("\n========== K-means 聚类 (k = %d) ==========\n", K_FINAL))

set.seed(SEED)
km_final <- kmeans(Z_mat, centers = K_FINAL,
                   nstart = KMEANS_NSTART, iter.max = KMEANS_ITERMAX)
cat(sprintf("  Cluster sizes: %s\n",
            paste(table(km_final$cluster), collapse = " / ")))


## ================================================================
## ================================================================

cat("\n========== 输出表格 ==========\n")

meta_df <- df_tpm %>%
  select(any_of(c("gene","Gene","locus_tag","product","Rel.Pos")))
if ("Gene" %in% colnames(meta_df) && !"gene" %in% colnames(meta_df))
  meta_df <- rename(meta_df, gene = Gene)

cluster_assignments <- tibble(
  gene    = rownames(Z_mat),
  cluster = as.integer(km_final$cluster[rownames(Z_mat)])
) %>%
  left_join(df_res %>% select(gene, direction,
                               ends_with("_l2fc_mean"),
                               ends_with("_l2fc_max_abs"),
                               median_l2fc_mean, max_l2fc_max_abs),
            by = "gene") %>%
  left_join(meta_df, by = "gene") %>%
  select(gene, locus_tag, product, Rel.Pos,
         direction, cluster,
         median_l2fc_mean, max_l2fc_max_abs, everything()) %>%
  arrange(cluster, direction, desc(max_l2fc_max_abs))

fwrite(cluster_assignments, file.path(DIR_TBL, "cluster_assignments.csv"))

cluster_summary <- cluster_assignments %>%
  group_by(cluster) %>%
  summarise(
    n_genes    = n(),
    n_up       = sum(direction == "Up"),
    n_down     = sum(direction == "Down"),
    top5_maxFC  = paste(gene[order(-max_l2fc_max_abs)][1:min(5,n())],
                        collapse = ";"),
    .groups = "drop"
  )
fwrite(cluster_summary, file.path(DIR_TBL, "cluster_summary.csv"))

fwrite(data.frame(cluster = seq_len(K_FINAL), km_final$centers,
                  check.names = FALSE),
       file.path(DIR_TBL, "kmeans_centers.csv"))


cat("  Saved: cluster_assignments.csv, cluster_summary.csv,\n")
cat("         kmeans_centers.csv
")


## ================================================================
## ================================================================

cat("\n========== 聚类诊断图 ==========\n")

make_long_z <- function(mat, clusters) {
  info_sub <- df_info %>%
    select(Colnames, Group, all_of(TIME_COL)) %>%
    rename(Time = !!TIME_COL)
  as.data.frame(mat) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "key", values_to = "z") %>%
    separate(key, into = c("Group","Colnames"), sep = "\\|", extra = "merge") %>%
    left_join(info_sub %>% rename(Group_check = Group), by = "Colnames") %>%
    filter(Group == Group_check) %>%
    select(-Group_check) %>%
    mutate(cluster = clusters[gene])
}

long_z <- make_long_z(Z_mat, km_final$cluster)

n_per_cluster <- cluster_assignments %>% dplyr::count(cluster, name = "n_genes")

agg_z <- long_z %>%
  group_by(cluster, Group, Time) %>%
  summarise(mean_z = mean(z, na.rm = TRUE),
            sd_z   = sd(z,   na.rm = TRUE), .groups = "drop") %>%
  left_join(n_per_cluster, by = "cluster") %>%
  mutate(cluster_lab = sprintf("Cluster %d  (n = %d)", cluster, n_genes))

p_means <- ggplot(agg_z, aes(x = Time, y = mean_z,
                              colour = Group, fill = Group, group = Group)) +
  geom_ribbon(aes(ymin = mean_z - sd_z, ymax = mean_z + sd_z),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dotted",
             colour = "grey40", linewidth = 0.35) +
  scale_colour_manual(values = PAL_GROUP) +
  scale_fill_manual(values   = PAL_GROUP) +
  facet_wrap(~ cluster_lab, ncol = 2, scales = "free_y") +
  labs(title = sprintf("Cluster mean trajectories  (k = %d)", K_FINAL),
       subtitle = "y = z-score(TPM);  shading = ±1 SD",
       x = "Time (h)", y = "z-score",
       colour = "Group", fill = "Group") +
  theme_diag() + theme(legend.position = "bottom")

ggsave(file.path(DIR_DIAG, "02_cluster_means_zscore.pdf"),
       p_means, width = 9, height = 7)
cat("  Saved: 02_cluster_means_zscore.pdf\n")

dir_df <- cluster_assignments %>%
  dplyr::count(cluster, direction) %>%
  left_join(n_per_cluster, by = "cluster") %>%
  mutate(cluster_lab = sprintf("Cluster %d\n(n=%d)", cluster, n_genes))

p_dir <- ggplot(dir_df,
                aes(x = factor(cluster), y = n, fill = direction)) +
  geom_col(position = "stack", colour = "black",
           linewidth = 0.3, width = 0.65) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 3, colour = "white", fontface = "bold") +
  scale_fill_manual(values = PAL_DIR) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Cluster × direction",
       x = "Cluster", y = "Number of genes", fill = "Direction") +
  theme_diag()

ggsave(file.path(DIR_DIAG, "03_direction_composition.pdf"),
       p_dir, width = 5, height = 4)
cat("  Saved: 03_direction_composition.pdf\n")


## ================================================================
##
## ================================================================

cat("\n========== Per-cluster 详细图 ==========\n")

plot_one_cluster <- function(cid) {
  genes_c <- names(km_final$cluster)[km_final$cluster == cid]
  meta_c  <- cluster_assignments %>% filter(cluster == cid)
  n_g     <- length(genes_c)

  long_c  <- long_z %>% filter(gene %in% genes_c)
  mean_c  <- long_c %>%
    group_by(Group, Time) %>%
    summarise(mean_z = mean(z, na.rm = TRUE), .groups = "drop")

  pA <- ggplot() +
    geom_line(data = long_c,
              aes(x = Time, y = z, group = interaction(gene, Group)),
              colour = "grey75", linewidth = 0.15, alpha = 0.35) +
    geom_line(data = mean_c,
              aes(x = Time, y = mean_z, colour = Group, group = Group),
              linewidth = 0.9) +
    geom_point(data = mean_c, aes(x = Time, y = mean_z, colour = Group),
               size = 1.4) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dotted",
               colour = "grey40", linewidth = 0.35) +
    scale_colour_manual(values = PAL_GROUP) +
    labs(title = "Individual genes & cluster mean",
         subtitle = sprintf("n = %d  |  Up = %d, Down = %d",
                            n_g, sum(meta_c$direction == "Up"),
                            sum(meta_c$direction == "Down")),
         x = "Time (h)", y = "z-score", colour = "Group") +
    theme_diag() + theme(legend.position = "bottom")

  pB <- ggplot(meta_c %>% dplyr::count(direction),
               aes(x = direction, y = n, fill = direction)) +
    geom_col(width = 0.6, colour = "black", linewidth = 0.3) +
    geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = PAL_DIR, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Direction", x = NULL, y = "n genes") +
    theme_diag()

  pC <- ggplot(meta_c,
               aes(x = direction, y = max_l2fc_max_abs, fill = direction)) +
    geom_violin(alpha = 0.4, linewidth = 0.4, trim = TRUE) +
    geom_jitter(aes(colour = direction),
                width = 0.12, size = 0.9, alpha = 0.65) +
    scale_fill_manual(values = PAL_DIR, guide = "none") +
    scale_colour_manual(values = c("Up"="#A0392B","Down"="#2B6A99"),
                        guide = "none") +
    labs(title = "Max |log2(Rel.TPM)|",
         x = NULL, y = "max |log2(Rel.TPM)|") +
    theme_diag()


  title_row <- ggdraw() +
    draw_label(sprintf("Cluster %d  |  n = %d  |  k = %d",
                       cid, n_g, K_FINAL),
               fontface = "bold", size = 11, x = 0.5, hjust = 0.5)

  plot_grid(title_row,
            pA + labs(tag = "A"),
            plot_grid(pB, pC, ncol = 2, align = "hv",
                      labels = c("B","C")),
            ncol = 1, rel_heights = c(0.10, 0.55, 0.35))
}

for (cid in sort(unique(km_final$cluster))) {
  p   <- plot_one_cluster(cid)
  fn  <- sprintf("cluster_%d.pdf", cid)
  ggsave(file.path(DIR_PERCK, fn), p, width = 9, height = 10)
  cat(sprintf("  Saved: per_cluster/%s\n", fn))
}


## ================================================================
##
## ================================================================

cat("\n========== 计算 Rel.TPM ==========\n")

reltpm_list <- lapply(groups, function(g) {
  info_g   <- info_list[[g]] %>% arrange(.data[[TIME_COL]])
  pre_cols <- info_g %>% filter(Treatment == PRE_LABEL) %>% pull(Colnames)

  mat <- df_tpm %>%
    filter(.data[[GENE_COL]] %in% genes_488) %>%
    select(all_of(c(GENE_COL, info_g$Colnames)))
  genes_g <- mat[[GENE_COL]]
  mat_num <- as.matrix(mat[,-1,drop=FALSE])
  rownames(mat_num) <- genes_g
  storage.mode(mat_num) <- "numeric"

  pre_mean <- rowMeans(mat_num[, pre_cols, drop=FALSE], na.rm = TRUE)
  rel_mat  <- mat_num / (pre_mean + EPS)
  rel_mat[!is.finite(rel_mat)] <- NA_real_

  as.data.frame(rel_mat) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "Colnames", values_to = "relTPM") %>%
    left_join(info_g %>% select(Colnames, !!TIME_COL, Treatment),
              by = "Colnames") %>%
    rename(Time = !!TIME_COL) %>%
    left_join(cluster_assignments %>% select(gene, cluster),
              by = "gene") %>%
    filter(!is.na(cluster), is.finite(Time), is.finite(relTPM)) %>%
    mutate(Group = g)
})

reltpm_long <- bind_rows(reltpm_list)

reltpm_summ <- reltpm_long %>%
  group_by(Group, Treatment, Time, cluster, gene) %>%
  summarise(relTPM = mean(relTPM, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group, Treatment, Time, cluster) %>%
  summarise(
    n_genes = n_distinct(gene),
    mean_val = mean(relTPM, na.rm = TRUE),
    sd_val   = sd(relTPM,   na.rm = TRUE),
    sem_val  = sd_val / sqrt(pmax(n_genes, 1)),
    t_crit   = qt(0.975, df = pmax(n_genes - 1, 1)),
    ci_lo    = mean_val - t_crit * sem_val,
    ci_hi    = mean_val + t_crit * sem_val,
    .groups  = "drop"
  )

fwrite(reltpm_summ, file.path(DIR_TBL, "reltpm_summary_allgroups.csv"))
cat("  Saved: reltpm_summary_allgroups.csv\n")


## ================================================================
## ================================================================

build_gene_order <- function() {
  gene_cor <- sapply(rownames(Z_mat), function(g) {
    cl <- km_final$cluster[g]
    if (is.na(cl)) return(NA)
    cor(Z_mat[g, ], km_final$centers[cl, ], method = "pearson")
  })

  cluster_assignments %>%
    mutate(cor_to_centroid = gene_cor[gene]) %>%
    mutate(cluster = factor(cluster, levels = CLUSTER_ORDER)) %>%
    arrange(cluster, desc(cor_to_centroid)) %>%
    pull(gene)
}

build_y_map <- function(gene_order) {
  y_map  <- numeric(0)
  blocks <- list()
  cur    <- 0
 for (k in CLUSTER_ORDER) {
    genes_k <- gene_order[cluster_assignments$cluster[
      match(gene_order, cluster_assignments$gene)] == k]
    genes_k <- genes_k[!is.na(genes_k)]
    if (!length(genes_k)) next
    y_seq   <- seq(cur + 1, cur + length(genes_k))
    y_map   <- c(y_map, setNames(y_seq, genes_k))
    blocks[[as.character(k)]] <- tibble(
      cluster = k,
      ymin    = min(y_seq) - 0.5,
      ymax    = max(y_seq) + 0.5
    )
    cur <- cur + length(genes_k) + CLUSTER_GAP_ROWS
  }
  list(y_map = y_map, blocks = bind_rows(blocks))
}

plot_irregular_heatmap <- function(g, y_map, cluster_blocks,
                                   file_pdf, width = 7.2, height = 9.5) {
  info_g  <- info_list[[g]] %>% arrange(.data[[TIME_COL]])
  times   <- info_g[[TIME_COL]]
  cols_g  <- intersect(paste0(g, "|", info_g$Colnames), colnames(Z_mat))

  mat_g   <- pmax(pmin(Z_mat[names(y_map), cols_g, drop=FALSE],
                        Z_RANGE[2]), Z_RANGE[1])

  xs    <- make_tilespec(times)
  ax    <- make_axis_params(times)
  xpad  <- diff(c(ax$x_left, ax$x_right)) * MARK_XPAD
  segL  <- diff(c(ax$x_left, ax$x_right)) * MARK_SEG
  txtO  <- diff(c(ax$x_left, ax$x_right)) * MARK_TXT

  tile_map <- tibble(col_key = cols_g, xmin = xs$xmin, xmax = xs$xmax)

  df_long <- as.data.frame(mat_g) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "col_key", values_to = "z") %>%
    left_join(tile_map, by = "col_key") %>%
    mutate(y = y_map[gene], ymin = y - 0.5, ymax = y + 0.5)

  y_max_all <- max(df_long$ymax, na.rm = TRUE)
  cb_plot   <- cluster_blocks %>%
    mutate(x_left = ax$x_left, x_right = ax$x_right)

  p <- ggplot(df_long) +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z),
              color = NA) +
    scale_fill_z() +
    scale_x_continuous(expand = c(0,0), limits = ax$limits,
                       breaks = ax$breaks, labels = ax$labels) +
    scale_y_reverse(expand = c(0,0), breaks = NULL) +
    coord_cartesian(xlim = c(ax$x_left, ax$x_right + xpad),
                    ylim = c(y_max_all, 0.5), clip = "off") +
    geom_rect(data = cb_plot,
              aes(xmin=x_left, xmax=x_right, ymin=ymin, ymax=ymax),
              inherit.aes = FALSE, fill = NA,
              color = "black", linewidth = CLUSTER_LW) +
    geom_rug(data = tibble(t = times), aes(x = t),
             inherit.aes = FALSE, sides = "b",
             length = unit(RUG_NPC, "npc"), linewidth = 0.3) +
    labs(title = g, x = "Time (h)", y = NULL) +
    theme_heatmap()


  ggsave(file_pdf, p, width = 8, height = 12, units = "in")
  invisible(p)
}

plot_equal_heatmap <- function(g_list, gene_order, file_pdf,
                               width = 9.0, height = 10.0) {
  clu_levels <- intersect(CLUSTER_ORDER, unique(km_final$cluster[gene_order]))
  clu_pal    <- PAL_CLUSTER[as.character(clu_levels)]
  clusters_o <- km_final$cluster[gene_order]
  df_row     <- data.frame(Cluster = factor(clusters_o, levels = clu_levels))
  rownames(df_row) <- gene_order

  ha_left <- rowAnnotation(
    Cluster = df_row$Cluster,
    col     = list(Cluster = clu_pal),
    width   = unit(5, "mm")
  )
  ht_list <- lapply(seq_along(g_list), function(i) {
    g       <- g_list[i]
    info_g  <- info_list[[g]] %>% arrange(.data[[TIME_COL]])
    cols_g  <- intersect(paste0(g,"|",info_g$Colnames), colnames(Z_mat))
    mat_g   <- pmax(pmin(Z_mat[gene_order, cols_g, drop=FALSE],
                          Z_RANGE[2]), Z_RANGE[1])
    colnames(mat_g) <- paste0(round(info_g[[TIME_COL]] * 60, 1), " min")

    Heatmap(
      mat_g, name = paste0("z_",g), col = col_fun_z, na_col = "grey90",
      cluster_rows = FALSE, cluster_columns = FALSE,
      show_row_names = FALSE, show_column_names = TRUE,
      column_title = g,
      column_title_gp = gpar(fontsize = BASE_SIZE+1, fontface = "bold"),
      row_split = df_row$Cluster, row_gap = unit(CH_ROW_GAP, "mm"),
      border = TRUE, column_names_rot = 45,
      column_names_gp = gpar(fontsize = BASE_SIZE-1),
      show_heatmap_legend = (i == length(g_list)),
      heatmap_legend_param = list(
        at = seq(Z_RANGE[1], Z_RANGE[2], by = 1), title = "z-score")
    )
  })

  ht_comb <- Reduce(`+`, ht_list)
  pdf(file_pdf, width = width, height = height)
  obj <- ha_left + ht_comb
  draw(obj, heatmap_legend_side = "right", annotation_legend_side = "right",
       ht_gap = unit(3, "mm"))
  dev.off()
}


## ================================================================
## ================================================================

build_line_plot <- function(g, show_legend = TRUE) {
  dfp <- reltpm_summ %>% filter(Group == g)
  if (!nrow(dfp)) return(NULL)

  t_min  <- min(dfp$Time, na.rm = TRUE)
  t_max  <- max(dfp$Time, na.rm = TRUE)

  brk_lo <- floor(t_min)      # 向下取整到整数小时
  brk_hi <- ceiling(t_max)    # 向上取整到整数小时

  x_ax <- list(
    limits = c(t_min, t_max),
    breaks = seq(brk_lo, brk_hi, by = 1)   # 整数小时间隔
  )
  ## ---------------------------------------------------------------

  y_hi <- max(ceiling(max(dfp$ci_hi, na.rm = TRUE)), 2)

  dfp <- dfp %>% mutate(cluster = factor(as.character(cluster),
                                          levels = names(PAL_CLUSTER)))
  n_tot <- n_distinct(
    reltpm_long$gene[reltpm_long$Group == g], na.rm = TRUE)

  ggplot(dfp, aes(x = Time, y = mean_val,
                  colour = cluster, fill = cluster, group = cluster)) +
    annotate("rect", xmin = SHADE_X1, xmax = SHADE_X2,
             ymin = -Inf, ymax = Inf, fill = SHADE_FILL, alpha = 0.20) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey55", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "black", linewidth = 0.3) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                alpha = RIBBON_ALPHA, colour = NA) +
    geom_line(linewidth = LINE_W) + geom_point(size = PT_SIZE) +
    scale_colour_manual(values = PAL_CLUSTER, drop = FALSE,
                        name = "Cluster") +
    scale_fill_manual(values = PAL_CLUSTER, drop = FALSE) +
    guides(fill = "none") +
    scale_x_continuous(
      limits = x_ax$limits,
      breaks = x_ax$breaks,
      labels = as.character(x_ax$breaks),   # 整数，无需as.integer截断
      expand = c(0, 0)
    ) +
    ## -------------------------------------------------------------
    scale_y_continuous(
      limits = c(0, y_hi),
      breaks = seq(0, y_hi, by = 1),
      labels = as.character(seq(0, y_hi, by = 1)),
      expand = c(0, 0)
    ) +
    labs(title    = paste0("Class-average Rel.TPM \u2014 ", g),
         subtitle = sprintf("mean \u00b1 95%% CI  (n = %d genes)", n_tot),
         x = "Time (h)", y = "Mean Rel.TPM (vs Pre baseline)") +
    theme_graphpad() +
    theme(legend.position = if (show_legend) "right" else "none")
}

export_line_plot <- function(g) {
  p_leg   <- build_line_plot(g, TRUE)
  p_noleg <- build_line_plot(g, FALSE)
  if (is.null(p_leg)) return(invisible(NULL))

  gt  <- no_clip(set_panel_size(p_noleg))
  leg <- cowplot::get_legend(p_leg)

  composed <- no_clip(ggplotGrob(plot_grid(
    ggdraw() + draw_grob(gt,  clip="off"),
    ggdraw() + draw_grob(leg, clip="off"),
    nrow = 1, rel_widths = c(PANEL_INCH, LEGEND_W), align = "h"
  )))
  save_pdf(file.path(DIR_LINE,
                     paste0("classavg_", gsub("[^A-Za-z0-9_]","_",g), ".pdf")),
           composed, PANEL_INCH + LEGEND_W, TOTAL_H)
  cat(sprintf("  Saved: classavg_%s.pdf\n", g))
}


## ================================================================
## ================================================================

cat("\n========== 逐组绘图 ==========\n")

gene_order    <- build_gene_order()
y_res         <- build_y_map(gene_order)
y_map         <- y_res$y_map
cluster_blocks <- y_res$blocks

for (g in groups) {
  cat(sprintf("\n--- %s ---\n", g))

  plot_irregular_heatmap(
    g, y_map, cluster_blocks,
    file.path(DIR_HEAT, paste0("heatmap_", g, "_irregular.pdf"))
  )
  cat(sprintf("  Saved: heatmap_%s_irregular.pdf\n", g))

  plot_equal_heatmap(
    g_list     = g,
    gene_order = gene_order,
    file_pdf   = file.path(DIR_HEAT, paste0("heatmap_", g, "_equal_CH.pdf"))
  )
  cat(sprintf("  Saved: heatmap_%s_equal_CH.pdf\n", g))

  export_line_plot(g)
}


## ================================================================
## ================================================================

cat("\n========== 三组合并图 ==========\n")

valid_g <- intersect(COMBINED_ORDER, groups)

if (length(valid_g) >= 2) {

  ax_list    <- lapply(valid_g, function(g) {
    make_axis_params(info_list[[g]][[TIME_COL]])
  }) %>% setNames(valid_g)

  xspec_list <- lapply(valid_g, function(g) {
    make_tilespec(info_list[[g]][[TIME_COL]])
  }) %>% setNames(valid_g)

  offsets <- numeric(length(valid_g))
  names(offsets) <- valid_g
  for (i in 2:length(valid_g)) {
    g_prev     <- valid_g[i-1]
    offsets[i] <- offsets[i-1] +
      (ax_list[[g_prev]]$x_right - ax_list[[g_prev]]$x_left) + GROUP_GAP_X
  }

  df_comb_list <- lapply(valid_g, function(g) {
    info_g  <- info_list[[g]] %>% arrange(.data[[TIME_COL]])
    times_g <- info_g[[TIME_COL]]
    cols_g  <- intersect(paste0(g,"|",info_g$Colnames), colnames(Z_mat))
    mat_g   <- pmax(pmin(Z_mat[gene_order, cols_g, drop=FALSE],
                          Z_RANGE[2]), Z_RANGE[1])
    tile_map <- tibble(col_key = cols_g,
                       xmin = xspec_list[[g]]$xmin + offsets[g],
                       xmax = xspec_list[[g]]$xmax + offsets[g])
    as.data.frame(mat_g) %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "col_key", values_to = "z") %>%
      left_join(tile_map, by = "col_key") %>%
      mutate(y = y_map[gene], ymin = y-0.5, ymax = y+0.5, Group = g)
  })
  df_comb <- bind_rows(df_comb_list)

  rug_comb <- lapply(valid_g, function(g) {
    tibble(x_rug = info_list[[g]][[TIME_COL]] + offsets[g], Group = g)
  }) %>% bind_rows()

  all_breaks <- all_labels <- c()
  for (g in valid_g) {
    ax <- ax_list[[g]]; offset <- offsets[g]
    in_rng <- ax$breaks >= ax$x_left & ax$breaks <= ax$x_right
    all_breaks <- c(all_breaks, ax$breaks[in_rng] + offset)
    all_labels <- c(all_labels, ax$labels[in_rng])
  }

  x_gl <- min(df_comb$xmin, na.rm=TRUE)
  x_gr <- max(df_comb$xmax, na.rm=TRUE)
  g_last  <- valid_g[length(valid_g)]
  ax_last <- ax_list[[g_last]]
  xpad_last <- (ax_last$x_right - ax_last$x_left) * MARK_XPAD
  segL_last <- (ax_last$x_right - ax_last$x_left) * MARK_SEG
  txtO_last <- (ax_last$x_right - ax_last$x_left) * MARK_TXT
  x_rl      <- ax_last$x_right + offsets[g_last]

  cb_comb <- lapply(valid_g, function(g) {
    cluster_blocks %>%
      mutate(x_left  = ax_list[[g]]$x_left  + offsets[g],
             x_right = ax_list[[g]]$x_right + offsets[g], Group = g)
  }) %>% bind_rows()

  title_df <- tibble(
    Group = valid_g,
    x_mid = sapply(valid_g, function(g)
      (ax_list[[g]]$x_left + ax_list[[g]]$x_right)/2 + offsets[g]),
    y_pos = -3
  )

  sep_df <- if (length(valid_g) >= 2) tibble(x = sapply(2:length(valid_g), function(i) {
    g_p <- valid_g[i-1]; g_c <- valid_g[i]
    (ax_list[[g_p]]$x_right + offsets[g_p] +
     ax_list[[g_c]]$x_left  + offsets[g_c]) / 2
  })) else NULL


  p_comb <- ggplot(df_comb) +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=z),
              color = NA) +
    scale_fill_z() +
    scale_x_continuous(expand=c(0,0), breaks=all_breaks, labels=all_labels) +
    scale_y_reverse(expand=c(0,0), breaks=NULL) +
    coord_cartesian(xlim=c(x_gl, x_gr+xpad_last),
                    ylim=c(max(df_comb$ymax,na.rm=TRUE), 0.5), clip="off") +
    geom_rect(data=cb_comb,
              aes(xmin=x_left, xmax=x_right, ymin=ymin, ymax=ymax,
                  group=interaction(cluster,Group)),
              inherit.aes=FALSE, fill=NA, color="black", linewidth=CLUSTER_LW) +
    geom_text(data=title_df, aes(x=x_mid, y=y_pos, label=Group),
              inherit.aes=FALSE, fontface="bold",
              size=BASE_SIZE/ggplot2::.pt) +
    geom_rug(data=rug_comb, aes(x=x_rug), inherit.aes=FALSE,
             sides="b", length=unit(RUG_NPC,"npc"), linewidth=0.3) +
    labs(x = "Time (h)", y = NULL) +
    theme_heatmap() + theme(plot.margin = margin(12,40,5.5,8))

  if (!is.null(sep_df))
    p_comb <- p_comb + geom_vline(data=sep_df, aes(xintercept=x),
                                  inherit.aes=FALSE, linetype="dashed",
                                  colour="grey60", linewidth=0.4)


  total_w <- sum(sapply(valid_g, function(g)
    ax_list[[g]]$x_right - ax_list[[g]]$x_left)) / 50 * 2.0 + 3.0

  ggsave(file.path(DIR_HEAT, "heatmap_combined_irregular.pdf"),
         p_comb, width=max(total_w,12), height=10, units="in")
  cat("  Saved: heatmap_combined_irregular.pdf\n")
}

if (length(valid_g) >= 2) {
  plot_equal_heatmap(
    g_list     = valid_g,
    gene_order = gene_order,
    file_pdf   = file.path(DIR_HEAT, "heatmap_combined_equal_CH.pdf"),
    width      = 10* length(valid_g) + 3.0,
    height     = 10.0
  )
  cat("  Saved: heatmap_combined_equal_CH.pdf\n")
}

if (length(valid_g) >= 2) {
  gt_list <- lapply(valid_g, function(g) {
    p <- build_line_plot(g, FALSE)
    if (is.null(p)) return(NULL)
    no_clip(set_panel_size(p + labs(title=g, subtitle=NULL)))
  }) %>% `[`(!sapply(., is.null))

  if (length(gt_list)) {
    shared_leg <- cowplot::get_legend(build_line_plot(valid_g[1], TRUE))
    composed <- no_clip(ggplotGrob(plot_grid(
      plot_grid(plotlist = lapply(gt_list, function(gt)
        ggdraw() + draw_grob(gt, clip="off")),
        nrow=1, rel_widths=rep(PANEL_INCH, length(gt_list))),
      ggdraw() + draw_grob(shared_leg, clip="off"),
      nrow=1, rel_widths=c(PANEL_INCH*length(gt_list), LEGEND_W)
    )))
    save_pdf(file.path(DIR_LINE, "classavg_combined.pdf"),
             composed, PANEL_INCH*length(gt_list)+LEGEND_W, TOTAL_H)
    cat("  Saved: classavg_combined.pdf\n")
  }
}


## ================================================================
## ================================================================

cat(sprintf("\n%s\n", strrep("=", 55)))
cat(sprintf("  响应基因: %d  |  k = %d\n",
            nrow(cluster_assignments), K_FINAL))
cat(sprintf("  Cluster sizes: %s\n",
            paste(table(cluster_assignments$cluster), collapse=" / ")))
cat(sprintf("  输出目录: %s\n", normalizePath(DIR_OUT)))
cat(sprintf("%s\n", strrep("=", 55)))


cat("\nDone.\n")


## ================================================================
## ================================================================
# Step6_clustering/

## ================================================================

}


# Analysis section: Cluster_COG_TF_enrichment
Cluster_COG_TF_enrichment <- function() {


# ============================================================================
# ============================================================================

# ============================================================
# ============================================================

PATH_RNA      <- "rna_TPM.csv"
PATH_EXPINFO  <- "rna_metadata.csv"
PATH_CLUSTER  <- "cluster_assignments.csv"
PATH_COG_ANN  <- "250715_Ecoli_K12_COG_UniqueGene_Merged.csv"
PATH_COG_FUN  <- "250717_cog-24.fun.csv"
PATH_ECO_TXT  <- "260113_Ecocyc_New-SmartTable_allColumn_RUsed.txt"
PATH_ECO_TF   <- "20260113_ManualChecked_Ecocyc_TF.csv"

DNAA_IN_CLUSTER <- "dnaA"
DNAA_IN_RNA     <- "Insert_dnaA_frontfrag"

RUN_COG <- TRUE
RUN_TF  <- TRUE

TF_UNIVERSE <- "both"

Q_CUTOFF          <- 0.05
M_MIN_TF          <- 5
TOP_N_PER_CLUSTER <- 6

X_CLUSTER_ORDER_COG <- c("cluster2","cluster1","cluster4","cluster3")
Y_CLUSTER_ORDER_COG_T <- c("cluster3","cluster4","cluster1","cluster2")

FORCE_TF_FIRST   <- "DnaA"
X_CLUSTER_BLOCKS <- c("cluster2","cluster1","cluster4","cluster3")
Y_CLUSTER_LEVELS <- c("cluster3","cluster4","cluster1","cluster2")

FIG_DPI       <- 300
PLOT_H_MAIN   <- 3.4
PLOT_H_SI     <- 3.4
PLOT_W_MIN    <- 9.0
PLOT_W_PER_TF <- 0.28

OUTDIR <- "Step7_cluster_enrichment"

# ============================================================
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(purrr)
})

message("[00] 开始数据读入与预处理...")

rna_df <- analysis_read_tibble(PATH_RNA, show_col_types = FALSE)
tpm_cols <- colnames(rna_df)[str_detect(colnames(rna_df), "TPM")]
stopifnot(length(tpm_cols) == 72)
rna_genes <- rna_df$gene %>% na.omit() %>% as.character() %>% unique()
message(sprintf("  RNA可测基因数: %d", length(rna_genes)))

exp_df <- analysis_read_tibble(PATH_EXPINFO, show_col_types = FALSE)
stopifnot(all(exp_df$Colnames %in% tpm_cols), all(tpm_cols %in% exp_df$Colnames))
message(sprintf("  样品数: %d", nrow(exp_df)))

cluster_raw <- analysis_read_tibble(PATH_CLUSTER, show_col_types = FALSE)
cluster_tbl <- cluster_raw %>%
  mutate(gene = if_else(gene == DNAA_IN_CLUSTER, DNAA_IN_RNA, gene)) %>%
  transmute(gene = as.character(gene), cluster = paste0("cluster", as.integer(cluster))) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct()
n_replaced <- sum(cluster_raw$gene == DNAA_IN_CLUSTER, na.rm = TRUE)
message(sprintf("  响应基因数: %d（%s→%s 替换%d处）",
                nrow(cluster_tbl), DNAA_IN_CLUSTER, DNAA_IN_RNA, n_replaced))

cog_raw <- analysis_read_tibble(PATH_COG_ANN, show_col_types = FALSE)
gene_cate_long <- cog_raw %>%
  select(Gene, COG_Category) %>%
  filter(!is.na(Gene), !is.na(COG_Category), COG_Category != "") %>%
  separate_rows(COG_Category, sep = ";") %>%
  mutate(Gene = str_trim(as.character(Gene)), COG_Category = str_trim(COG_Category)) %>%
  filter(COG_Category != "") %>%
  distinct(Gene, COG_Category)
bg_cog <- intersect(rna_genes, unique(gene_cate_long$Gene))
message(sprintf("  COG注释基因数: %d，COG富集背景 bg_cog: %d",
                n_distinct(gene_cate_long$Gene), length(bg_cog)))

cogfun_raw <- analysis_read_tibble(PATH_COG_FUN, show_col_types = FALSE)
colnames(cogfun_raw) <- str_trim(colnames(cogfun_raw))
needed_cogfun_cols <- c("Functional category ID", "Functional category description",
                        "Major_Group_ID", "Major_Group")
if (!all(needed_cogfun_cols %in% colnames(cogfun_raw))) {
  stop("COG功能表缺少必要列")
}
cate_order_df <- cogfun_raw %>%
  mutate(Order = row_number()) %>%
  transmute(
    COG_Category = str_trim(`Functional category ID`),
    COG_Desc     = `Functional category description`,
    Major_Group_ID = as.character(Major_Group_ID),
    Major_Group    = Major_Group,
    Order
  )
message(sprintf("  COG类别数（功能表）: %d", nrow(cate_order_df)))

eco_raw <- analysis_read_tsv(PATH_ECO_TXT, col_types = cols(.default = "c"), locale = locale(encoding = "UTF-8"))
colnames(eco_raw) <- colnames(eco_raw) %>%
  str_replace_all('"', '') %>% str_replace_all('\r', '') %>% str_trim()
eco_raw <- eco_raw %>% mutate(Gene = str_trim(str_replace_all(Gene, '\r', '')))
TF_COL <- "Regulation - transcription factor regulators of gene"
stopifnot(TF_COL %in% colnames(eco_raw))

eco_genes_all <- eco_raw$Gene %>% na.omit() %>% unique()
U1 <- intersect(rna_genes, eco_genes_all)
eco_genes_tfannot <- eco_raw %>%
  filter(!is.na(.data[[TF_COL]]), str_trim(as.character(.data[[TF_COL]])) != "",
         str_trim(as.character(.data[[TF_COL]])) != "NA") %>%
  pull(Gene) %>% na.omit() %>% unique()
U2 <- intersect(rna_genes, eco_genes_tfannot)
message(sprintf("  EcoCyc全收录: %d基因，U1: %d，U2: %d",
                length(eco_genes_all), length(U1), length(U2)))

tf_csv <- analysis_read_csv(PATH_ECO_TF, fileEncoding = "UTF-8-BOM")
stopifnot(all(c("tf_raw","tf_short_manual","genes_all_txt") %in% colnames(tf_csv)))
gene2tf <- tf_csv %>%
  select(tf = tf_short_manual, genes_all_txt) %>%
  mutate(gene_list = str_split(genes_all_txt, "/")) %>%
  unnest(gene_list) %>%
  transmute(gene = str_trim(gene_list), tf = str_trim(tf)) %>%
  filter(gene != "", tf != "", !is.na(gene), !is.na(tf)) %>%
  distinct(gene, tf)
message(sprintf("  gene→tf映射对数: %d", nrow(gene2tf)))

tf_sets_full <- gene2tf %>%
  group_by(tf) %>%
  summarise(genes = list(sort(unique(gene))), M_all = n_distinct(gene), .groups = "drop")
message(sprintf("  TF总数: %d", nrow(tf_sets_full)))
message("[00] 预处理完成 ✓\n")

# ============================================================
# ============================================================

if (RUN_COG) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(scales)
    library(tibble)
  })

  message("[01] 开始COG富集分析（by_cluster）...")
  cog_outdir <- file.path(OUTDIR, "COG_enrichment", "by_cluster")
  dir.create(cog_outdir, recursive = TRUE, showWarnings = FALSE)

  fisher_enrich_one <- function(fg_genes, cate_genes, bg_genes) {
    N <- length(bg_genes)
    M <- sum(bg_genes %in% cate_genes)
    n <- sum(bg_genes %in% fg_genes)
    k <- sum(fg_genes %in% cate_genes)
    if (n == 0 || M == 0) {
      return(tibble(k=k, n=n, M=M, N=N, GeneRatio=NA_real_, bg_rate=NA_real_,
                    OR=NA_real_, log2OR=NA_real_, pvalue=1, hit_genes=""))
    }
    a <- k; b <- n - k; c <- M - k; d <- N - M - (n - k)
    ft <- tryCatch(fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater"),
                   error = function(e) list(p.value=1, estimate=NA_real_))
    OR_val <- if (!is.na(ft$estimate) && ft$estimate > 0) unname(ft$estimate)
              else if (b == 0 || c == 0) Inf else (a*d)/(b*c)
    tibble(k=k, n=n, M=M, N=N, GeneRatio=k/n, bg_rate=M/N,
           OR=OR_val, log2OR=if(is.finite(OR_val) && OR_val>0) log2(OR_val) else NA_real_,
           pvalue=ft$p.value, hit_genes=paste(sort(intersect(fg_genes, cate_genes)), collapse="/"))
  }

  clusters_use <- sort(unique(cluster_tbl$cluster))
  fg_list <- setNames(
    lapply(clusters_use, function(cl) intersect(cluster_tbl$gene[cluster_tbl$cluster == cl], bg_cog)),
    clusters_use
  )
  for (nm in names(fg_list)) message(sprintf("  前景 [%s]: %d 基因", nm, length(fg_list[[nm]])))

  categories <- cate_order_df$COG_Category
  enrich_list <- list()
  for (grp_name in names(fg_list)) {
    fg <- fg_list[[grp_name]]
    if (length(fg) == 0) next
    for (cate in categories) {
      cate_genes <- gene_cate_long %>% filter(COG_Category == cate) %>% pull(Gene) %>% intersect(bg_cog)
      if (length(cate_genes) == 0) next
      res <- fisher_enrich_one(fg, cate_genes, bg_cog)
      enrich_list[[length(enrich_list)+1]] <- res %>% mutate(Group=grp_name, COG_Category=cate, .before=1)
    }
  }
  enrich_df <- bind_rows(enrich_list) %>%
    mutate(qvalue = p.adjust(pvalue, method="BH"),
           neg_log10q = if_else(qvalue > 0, -log10(qvalue), NA_real_),
           is_sig = !is.na(qvalue) & qvalue < Q_CUTOFF) %>%
    left_join(cate_order_df %>% select(COG_Category, COG_Desc, Order), by="COG_Category")

  write_tsv(enrich_df, file.path(cog_outdir, "COG_enrichment_all.tsv"))
  write_tsv(enrich_df %>% filter(is_sig) %>% arrange(Group, qvalue, desc(GeneRatio)),
            file.path(cog_outdir, "COG_enrichment_sig.tsv"))
  message(sprintf("  显著富集条目数（q<%.2f）: %d", Q_CUTOFF, sum(enrich_df$is_sig, na.rm=TRUE)))

  label_order <- cate_order_df %>% arrange(Order) %>%
    mutate(Label = paste0(COG_Category, " — ", COG_Desc)) %>% pull(Label)
  grp_levels <- X_CLUSTER_ORDER_COG[X_CLUSTER_ORDER_COG %in% unique(enrich_df$Group)]
  plot_df <- enrich_df %>%
    mutate(Label = factor(paste0(COG_Category, " — ", COG_Desc), levels = rev(label_order)),
           Group = factor(Group, levels = grp_levels),
           log2OR_plot = if_else(is_sig, log2OR, NA_real_),
           neg_log10q_plot = if_else(is_sig, neg_log10q, NA_real_))

  max_nlq <- max(plot_df$neg_log10q_plot, na.rm=TRUE); if (!is.finite(max_nlq) || max_nlq<=0) max_nlq <- 1
  abs_lor <- max(abs(plot_df$log2OR_plot), na.rm=TRUE); if (!is.finite(abs_lor) || abs_lor<=0) abs_lor <- 1
  COLORS_DIV <- c("#3B4CC0","#88A1D2","#F7F7F7","#F3A982","#B40426")
  NS_PT_COL <- "grey70"; NS_TILE_COL <- "grey90"

  n_grps <- length(grp_levels)
  plot_w <- max(4.5, 1.6*n_grps+3.0); plot_h <- max(5.5, 0.3*nrow(cate_order_df)+1.5)
  theme_cog <- function() theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(),
                                                         axis.title.x=element_blank(),
                                                         axis.text.x=element_text(angle=30,hjust=1),
                                                         axis.text.y=element_text(size=9),
                                                         plot.margin=margin(6,6,6,6))

  # Dotplot
  p_dot <- ggplot() +
    geom_point(data=plot_df %>% filter(!is_sig), aes(x=Group, y=Label, size=GeneRatio),
               color=NS_PT_COL, alpha=0.7) +
    geom_point(data=plot_df %>% filter(is_sig), aes(x=Group, y=Label, size=GeneRatio, color=neg_log10q_plot),
               alpha=0.95) +
    scale_size_continuous(name="GeneRatio (k/n)", range=c(1,9), labels=percent_format(accuracy=1)) +
    scale_color_gradientn(name=expression(-log[10]~"(q)"), colours=COLORS_DIV,
                          limits=c(0,max_nlq), oob=squish) +
    theme_cog() +
    labs(title="COG enrichment dotplot (by_cluster)",
         subtitle=sprintf("Background: RNA∩COG (N=%d); FDR<%.2f (BH global)", length(bg_cog), Q_CUTOFF))
  ggsave(file.path(cog_outdir, "COG_dotplot.png"), p_dot, width=plot_w, height=plot_h, dpi=FIG_DPI)
  ggsave(file.path(cog_outdir, "COG_dotplot.pdf"), p_dot, width=plot_w, height=plot_h, device=cairo_pdf)

  group_levels_T <- Y_CLUSTER_ORDER_COG_T
  plot_df_T <- plot_df %>% mutate(Group=factor(Group, levels=group_levels_T),
                                  Label=factor(Label, levels=label_order))
  p_dot_T <- ggplot() +
    geom_point(data=plot_df_T %>% filter(!is_sig), aes(x=Label, y=Group, size=GeneRatio),
               color=NS_PT_COL, alpha=0.7) +
    geom_point(data=plot_df_T %>% filter(is_sig), aes(x=Label, y=Group, size=GeneRatio, color=neg_log10q_plot),
               alpha=0.95) +
    scale_size_continuous(name="GeneRatio (k/n)", range=c(1,9), labels=percent_format(accuracy=1)) +
    scale_color_gradientn(name=expression(-log[10]~"(q)"), colours=COLORS_DIV,
                          limits=c(0,max_nlq), oob=squish) +
    theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(),
                                   axis.title=element_blank(),
                                   axis.text.x=element_text(angle=45,hjust=1,size=8),
                                   plot.margin=margin(6,6,6,6)) +
    labs(title="COG enrichment dotplot — transposed (by_cluster)",
         subtitle=sprintf("Background: RNA∩COG (N=%d); FDR<%.2f (BH global)", length(bg_cog), Q_CUTOFF))
  ggsave(file.path(cog_outdir, "COG_dotplot_transposed.png"), p_dot_T,
         width=plot_h, height=plot_w*0.6, dpi=FIG_DPI)
  ggsave(file.path(cog_outdir, "COG_dotplot_transposed.pdf"), p_dot_T,
         width=plot_h, height=plot_w*0.6, device=cairo_pdf)

  # Heatmap
  p_heat <- ggplot() +
    geom_tile(data=plot_df %>% filter(!is_sig), aes(x=Group, y=Label), fill=NS_TILE_COL, color="white", linewidth=0.3) +
    geom_tile(data=plot_df %>% filter(is_sig), aes(x=Group, y=Label, fill=log2OR_plot),
              color="white", linewidth=0.3) +
    scale_fill_gradientn(name="log2(OR)", colours=COLORS_DIV,
                         limits=c(-abs_lor, abs_lor), oob=squish, na.value=NS_TILE_COL) +
    scale_x_discrete(expand=c(0,0)) + theme_cog() +
    labs(title="COG enrichment heatmap (by_cluster)",
         subtitle=sprintf("Background: RNA∩COG (N=%d); FDR<%.2f (BH global)", length(bg_cog), Q_CUTOFF))
  ggsave(file.path(cog_outdir, "COG_heatmap.png"), p_heat, width=plot_w, height=plot_h, dpi=FIG_DPI)
  ggsave(file.path(cog_outdir, "COG_heatmap.pdf"), p_heat, width=plot_w, height=plot_h, device=cairo_pdf)

  message(sprintf("[01] COG富集完成 ✓ → %s\n", cog_outdir))
}

# ============================================================
# ============================================================
# Replace only the original: 02_tf_enrichment.R (by_cluster)
# Run after the original configuration, data preparation and COG section.
if (RUN_TF) {
  tf_bycluster_v62_results <- local({
    # absolute: filtered genes; relative: union of clustered response genes.
    PATH_TF_BACKGROUND <- "filtered_genes.csv"
    TF_ENRICHMENT_TYPES <- c("absolute", "relative")
    MAKE_TF_PLOTS <- TRUE

    analysis_read_tibble <- function(p) {
      if (!analysis_file_exists(p)) stop("File not found: ", p)
      d <- analysis_read_csv(p, check.names = FALSE, stringsAsFactors = FALSE,
                    fileEncoding = "UTF-8-BOM")
      names(d) <- sub("^\ufeff", "", names(d))
      d
    }
    clean <- function(x) {
      x <- trimws(as.character(x))
      x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
      x
    }
    key <- function(x) {
      x <- tolower(clean(x))
      x[!is.na(x) & x == tolower(DNAA_IN_RNA)] <- tolower(DNAA_IN_CLUSTER)
      x
    }
    safe_fraction <- function(a, b) if (b > 0) a / b else NA_real_

    v6_ora_from_map_complete <- function(foreground, background, term_map, min_term_size, gene_label_map = NULL) {
      foreground <- unique(intersect(foreground, background))
      background <- unique(background)
      map <- unique(term_map[term_map$gene_key %in% background & !is.na(term_map$term), , drop = FALSE])
      if (!length(background) || !nrow(map)) return(data.frame())
      terms <- sort(unique(map$term))
      N <- length(background)
      n <- length(foreground)
      ans <- lapply(terms, function(term) {
        targets <- unique(map$gene_key[map$term == term])
        K <- length(targets)
        if (K < min_term_size) return(NULL)
        hit_keys <- sort(unique(intersect(foreground, targets)))
        a <- length(hit_keys)
        b <- n - a
        c0 <- K - a
        d <- N - K - b
        if (any(c(a,b,c0,d) < 0)) return(NULL)
        hit_labels <- hit_keys
        if (!is.null(gene_label_map) && length(gene_label_map) && length(hit_keys)) {
          mapped <- unname(gene_label_map[hit_keys])
          use_mapped <- !is.na(mapped) & nzchar(mapped)
          hit_labels[use_mapped] <- mapped[use_mapped]
        }
        tab <- matrix(c(a,b,c0,d), nrow = 2, byrow = TRUE)
        ft_enrich <- stats::fisher.test(tab, alternative = "greater")
        ft_ci <- stats::fisher.test(tab, alternative = "two.sided", conf.level = 0.95)
        data.frame(
          term = term,
          foreground_n = n,
          background_n = N,
          overlap_n = a,
          background_term_n = K,
          hit_genes = paste(hit_labels, collapse = ";"),
          hit_gene_keys = paste(hit_keys, collapse = ";"),
          gene_ratio = safe_fraction(a, n),
          background_ratio = safe_fraction(K, N),
          odds_ratio = unname(ft_enrich$estimate),
          CI_low_95 = unname(ft_ci$conf.int[1]),
          CI_high_95 = unname(ft_ci$conf.int[2]),
          p_value = ft_enrich$p.value,
          stringsAsFactors = FALSE
        )
      })
      ans <- ans[!vapply(ans, is.null, logical(1))]
      if (!length(ans)) return(data.frame())
      do.call(rbind, ans)
    }

    bg_raw <- analysis_read_tibble(PATH_TF_BACKGROUND)
    gene_col <- intersect(c("gene", "Gene", "gene_key"), names(bg_raw))[1]
    if (is.na(gene_col)) stop("Background requires a gene column.")
    bg_ids <- clean(bg_raw[[gene_col]])
    bg_keys <- key(bg_ids)
    if (anyNA(bg_keys)) stop("Missing background gene IDs.")
    collisions <- vapply(split(bg_ids, bg_keys),
                         function(z) length(unique(z)) > 1L, logical(1))
    if (any(collisions)) stop("Background alias collision; review gene IDs: ",
                             paste(names(collisions)[collisions], collapse = ", "))
    bg <- sort(unique(bg_keys))
    labels <- setNames(bg_ids, bg_keys)

    if (!all(c("gene", "cluster") %in% names(cluster_tbl)))
      stop("cluster_tbl requires gene and cluster columns.")
    membership <- unique(data.frame(
      gene_key = key(cluster_tbl$gene),
      cluster = clean(cluster_tbl$cluster), stringsAsFactors = FALSE))
    if (anyNA(membership)) stop("Missing gene or cluster in cluster_tbl.")
    if (anyDuplicated(membership$gene_key))
      stop("A gene belongs to multiple clusters; resolve before enrichment.")
    outside <- setdiff(membership$gene_key, bg)
    if (length(outside)) stop("Cluster genes outside the filtered background: ",
                             paste(outside, collapse = ", "))
    responders <- sort(unique(membership$gene_key))
    if (!length(responders)) stop("No clustered response genes.")
    clusters <- sort(unique(membership$cluster))
    x_order <- unique(c(intersect(X_CLUSTER_BLOCKS, clusters), clusters))
    y_order <- unique(c(intersect(Y_CLUSTER_LEVELS, clusters), clusters))
    groups <- setNames(lapply(clusters, function(cl)
      membership$gene_key[membership$cluster == cl]), clusters)

    ann <- analysis_read_tibble(PATH_ECO_TF)
    if (!all(c("tf_short_manual", "genes_all_txt") %in% names(ann)))
      stop("TF annotation requires tf_short_manual and genes_all_txt.")
    rows <- lapply(seq_len(nrow(ann)), function(i) {
      tf <- clean(ann$tf_short_manual[i])
      txt <- clean(ann$genes_all_txt[i])
      if (is.na(tf) || is.na(txt)) return(NULL)
      data.frame(term = tf, gene_key = key(strsplit(txt, "/", fixed = TRUE)[[1]]))
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) stop("No TF-target annotation.")
    term_map <- unique(do.call(rbind, rows))
    term_map <- term_map[!is.na(term_map$gene_key), , drop = FALSE]
    if (length(M_MIN_TF) != 1L || !is.finite(M_MIN_TF) ||
        M_MIN_TF < 1 || M_MIN_TF != as.integer(M_MIN_TF))
      stop("M_MIN_TF must be a positive integer.")
    if (length(Q_CUTOFF) != 1L || !is.finite(Q_CUTOFF) ||
        Q_CUTOFF <= 0 || Q_CUTOFF >= 1) stop("Invalid Q_CUTOFF.")
    if (length(TOP_N_PER_CLUSTER) != 1L || !is.finite(TOP_N_PER_CLUSTER) ||
        TOP_N_PER_CLUSTER < 1 || TOP_N_PER_CLUSTER != as.integer(TOP_N_PER_CLUSTER))
      stop("TOP_N_PER_CLUSTER must be a positive integer.")
    if (!all(TF_ENRICHMENT_TYPES %in% c("absolute", "relative")))
      stop("Invalid enrichment type.")

    root <- file.path(OUTDIR, "TF_enrichment", "V62_by_cluster")
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    membership$original_gene <- unname(labels[membership$gene_key])
    write.csv(membership, file.path(root, "cluster_membership.csv"),
              row.names = FALSE)
    write.csv(data.frame(gene_key = bg, original_gene = unname(labels[bg]),
      clustered_response = bg %in% responders),
      file.path(root, "filtered_background.csv"), row.names = FALSE)
    term_map$in_filtered_background <- term_map$gene_key %in% bg
    term_map$in_response_union <- term_map$gene_key %in% responders
    write.csv(term_map, file.path(root, "TF_target_mapping.csv"),
              row.names = FALSE)

    can_plot <- MAKE_TF_PLOTS && requireNamespace("ggplot2", quietly = TRUE)
    if (MAKE_TF_PLOTS && !can_plot)
      warning("ggplot2 is unavailable; statistics are still exported.")
    results <- list()
    for (etype in TF_ENRICHMENT_TYPES) {
      universe <- if (etype == "absolute") bg else responders
      dir_tf <- file.path(root, etype)
      dir.create(dir_tf, recursive = TRUE, showWarnings = FALSE)
      tabs <- lapply(names(groups), function(cl) {
        z <- v6_ora_from_map_complete(groups[[cl]], universe, term_map,
                                      M_MIN_TF, labels)
        if (!nrow(z)) return(NULL)
        z$cluster <- cl
        z
      })
      tabs <- Filter(Negate(is.null), tabs)
      if (!length(tabs)) {
        write.csv(data.frame(status = "no_eligible_TFs", background_n = length(universe)),
          file.path(dir_tf, "analysis_status.csv"), row.names = FALSE)
        next
      }
      z <- do.call(rbind, tabs)
      # One BH family across ALL cluster x eligible-TF tests for this background.
      z$FDR <- p.adjust(z$p_value, method = "BH")
      z$FDR_family <- paste("dnaA_RNA", "cluster", "TF", etype, sep = "|")
      z$FDR_method <- "Benjamini-Hochberg"
      z$n_tests_in_family <- nrow(z)
      z$significant <- is.finite(z$FDR) & z$FDR < Q_CUTOFF
      z$FE <- z$gene_ratio / z$background_ratio
      z$log2FE <- log2(z$FE)
      z$neg_log10FDR <- -log10(pmax(z$FDR, .Machine$double.xmin))
      z$analysis_id <- "dnaA_RNA"
      z$group_type <- "cluster"
      z$database <- "TF"
      z$enrichment_type <- etype
      z$background_definition <- if (etype == "absolute")
        "All genes retained by expression/baseline filters" else
        "Union of all clustered responsive genes"
      z <- z[order(z$cluster, z$FDR, z$p_value, z$term), ]
      rownames(z) <- NULL
      results[[etype]] <- z
      write.csv(z, file.path(dir_tf, "TF_enrich_all.csv"), row.names = FALSE)
      write.csv(z[z$significant, ], file.path(dir_tf, "TF_enrich_sig.csv"),
                row.names = FALSE)

      # Legacy column names for the original downstream workflow.
      legacy <- data.frame(cluster = z$cluster, tf = z$term,
        k = z$overlap_n, n = z$foreground_n, M = z$background_term_n,
        N = z$background_n, GeneRatio = z$gene_ratio, bg_rate = z$background_ratio,
        FE = z$FE, log2FE = z$log2FE, OR = z$odds_ratio,
        pvalue = z$p_value, qvalue = z$FDR, neg_log10q = z$neg_log10FDR,
        is_sig = z$significant, geneID = z$hit_genes, universe = etype)
      write.table(legacy, file.path(dir_tf, "TF_enrich_all.tsv"),
                  sep = "\t", quote = TRUE, row.names = FALSE)
      write.table(legacy[legacy$is_sig, ], file.path(dir_tf, "TF_enrich_sig.tsv"),
                  sep = "\t", quote = TRUE, row.names = FALSE)
      selected <- do.call(rbind, lapply(x_order, function(cl) {
        d <- z[z$cluster == cl & z$significant, ]
        d <- d[order(d$FDR, -d$gene_ratio, d$term), ]
        head(d, TOP_N_PER_CLUSTER)
      }))
      write.csv(selected, file.path(dir_tf, "TF_main_selected.csv"),
                row.names = FALSE)
      summary <- do.call(rbind, lapply(clusters, function(cl) {
        d <- z[z$cluster == cl, ]
        data.frame(cluster = cl, foreground_n = length(groups[[cl]]),
          background_n = length(universe), eligible_TFs = nrow(d),
          zero_hit_TFs = sum(d$overlap_n == 0), significant_TFs = sum(d$significant),
          foreground_TF_annotated = sum(groups[[cl]] %in% term_map$gene_key),
          background_TF_annotated = sum(universe %in% term_map$gene_key),
          FDR_family = unique(z$FDR_family), n_tests_in_family = nrow(z))
      }))
      write.csv(summary, file.path(dir_tf, "cluster_summary.csv"), row.names = FALSE)

      message(etype, ": background=", length(universe), ", clusters=", length(clusters),
              ", tests=", nrow(z), ", significant TF-cluster pairs=", sum(z$significant))

      # Original-style TF dotplots / heatmaps; no publication module.
      if (!can_plot || !any(z$significant)) next
      for (view in c("main", "all_significant")) {
        src <- if (view == "main") selected else z[z$significant, ]
        terms <- unique(unlist(lapply(x_order, function(cl)
          src$term[src$cluster == cl]), use.names = FALSE))
        first <- terms[tolower(terms) == tolower(FORCE_TF_FIRST)]
        terms <- c(first, setdiff(terms, first))
        if (!length(terms)) next
        pd <- z[z$term %in% terms, ]
        pd$term <- factor(pd$term, levels = terms)
        pd$cluster <- factor(pd$cluster, levels = y_order)
        pd$colour_value <- ifelse(pd$significant, pd$neg_log10FDR, NA_real_)
        palette <- c("#3B4CC0", "#88A1D2", "#F7F7F7", "#F3A982", "#B40426")
        common <- list(ggplot2::theme_bw(base_size = 12),
          ggplot2::theme(axis.title = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)),
          ggplot2::labs(title = paste("TF enrichment:", etype, view),
            subtitle = paste("Background N =", length(universe),
                             "; BH across", nrow(z), "TF-cluster tests")))
        dot <- ggplot2::ggplot(pd, ggplot2::aes(term, cluster)) +
          ggplot2::geom_point(ggplot2::aes(size = gene_ratio, colour = colour_value)) +
          ggplot2::scale_colour_gradientn(colours = palette, na.value = "grey70",
                                          name = "-log10(FDR)") +
          ggplot2::scale_size_continuous(range = c(1.2, 10), name = "GeneRatio") + common
        heat <- ggplot2::ggplot(pd, ggplot2::aes(term, cluster, fill = colour_value)) +
          ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
          ggplot2::scale_fill_gradientn(colours = palette, na.value = "grey90",
                                        name = "-log10(FDR)") + common
        h <- if (view == "main") PLOT_H_MAIN else PLOT_H_SI
        w <- max(PLOT_W_MIN, PLOT_W_PER_TF * length(terms))
        for (kind in c("dotplot", "heatmap")) {
          p <- if (kind == "dotplot") dot else heat
          stem <- file.path(dir_tf, paste("TF", kind, view, sep = "_"))
          ggplot2::ggsave(paste0(stem, ".pdf"), p, width = w, height = h)
          ggplot2::ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = FIG_DPI)
        }
      }
    }
    results
  })
}


}


# Analysis section: Responsive_set_COG_TF_enrichment
Responsive_set_COG_TF_enrichment <- function() {


# Integrated All_responsive COG + V6.2 TF ORA + publication TF plots.
# Source text through the end of COG is preserved verbatim.
# Run this file to execute the complete step-specific analysis.
#      Run Responsive_set_COG_TF_enrichment() from this file.
# For a TF-only replacement, copy sections 4 and 5 over the old section 4.
# TF outputs: Step8_all_responsive_enrichment/TF_enrichment/tf/
# Publication plot script: point CONFIG$input_csv to TF__Response__absolute.csv.

# ============================================================================
# ============================================================================

# ============================================================
# ============================================================

PATH_RNA      <- "rna_TPM.csv"
PATH_EXPINFO  <- "rna_metadata.csv"
PATH_CLUSTER  <- "cluster_assignments.csv"
PATH_COG_ANN  <- "250715_Ecoli_K12_COG_UniqueGene_Merged.csv"
PATH_COG_FUN  <- "250717_cog-24.fun.csv"
PATH_ECO_TXT  <- "260113_Ecocyc_New-SmartTable_allColumn_RUsed.txt"
PATH_ECO_TF   <- "20260113_ManualChecked_Ecocyc_TF.csv"

DNAA_IN_CLUSTER <- "dnaA"
DNAA_IN_RNA     <- "Insert_dnaA_frontfrag"

RUN_COG <- TRUE
RUN_TF  <- TRUE

TF_UNIVERSE <- "both"

Q_CUTOFF   <- 0.05
M_MIN_TF   <- 5   # TF 富集时，Universe 中靶基因数的最低要求

FIG_DPI       <- 300
PLOT_H        <- 8      # 图高度（英寸）
PLOT_W        <- 10     # 图宽度（英寸）

OUTDIR <- "Step8_all_responsive_enrichment"

# ============================================================
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(purrr)
})

message("[0] 开始数据读入与预处理...")

rna_df <- analysis_read_tibble(PATH_RNA, show_col_types = FALSE)
tpm_cols <- colnames(rna_df)[str_detect(colnames(rna_df), "TPM")]
stopifnot(length(tpm_cols) == 72)
rna_genes <- rna_df$gene %>% na.omit() %>% as.character() %>% unique()
message(sprintf("  RNA可测基因数: %d", length(rna_genes)))

exp_df <- analysis_read_tibble(PATH_EXPINFO, show_col_types = FALSE)
stopifnot(all(exp_df$Colnames %in% tpm_cols), all(tpm_cols %in% exp_df$Colnames))
message(sprintf("  样品数: %d", nrow(exp_df)))

cluster_raw <- analysis_read_tibble(PATH_CLUSTER, show_col_types = FALSE)
cluster_tbl <- cluster_raw %>%
  mutate(gene = if_else(gene == DNAA_IN_CLUSTER, DNAA_IN_RNA, gene)) %>%
  transmute(gene = as.character(gene), cluster = paste0("cluster", as.integer(cluster))) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct()
n_replaced <- sum(cluster_raw$gene == DNAA_IN_CLUSTER, na.rm = TRUE)
message(sprintf("  响应基因总数: %d（%s→%s 替换%d处）",
                nrow(cluster_tbl), DNAA_IN_CLUSTER, DNAA_IN_RNA, n_replaced))

all_responsive_genes <- unique(cluster_tbl$gene)
message(sprintf("  所有响应基因（unique）: %d", length(all_responsive_genes)))

cog_raw <- analysis_read_tibble(PATH_COG_ANN, show_col_types = FALSE)
gene_cate_long <- cog_raw %>%
  select(Gene, COG_Category) %>%
  filter(!is.na(Gene), !is.na(COG_Category), COG_Category != "") %>%
  separate_rows(COG_Category, sep = ";") %>%
  mutate(Gene = str_trim(as.character(Gene)), COG_Category = str_trim(COG_Category)) %>%
  filter(COG_Category != "") %>%
  distinct(Gene, COG_Category)
bg_cog <- intersect(rna_genes, unique(gene_cate_long$Gene))
message(sprintf("  COG注释基因数: %d，COG富集背景 bg_cog: %d",
                n_distinct(gene_cate_long$Gene), length(bg_cog)))

cogfun_raw <- analysis_read_tibble(PATH_COG_FUN, show_col_types = FALSE)
colnames(cogfun_raw) <- str_trim(colnames(cogfun_raw))
needed_cols <- c("Functional category ID", "Functional category description",
                 "Major_Group_ID", "Major_Group")
if (!all(needed_cols %in% colnames(cogfun_raw))) stop("COG功能表缺少必要列")
cate_order_df <- cogfun_raw %>%
  mutate(Order = row_number()) %>%
  transmute(
    COG_Category = str_trim(`Functional category ID`),
    COG_Desc     = `Functional category description`,
    Major_Group_ID = as.character(Major_Group_ID),
    Major_Group    = Major_Group,
    Order
  )
message(sprintf("  COG类别数（功能表）: %d", nrow(cate_order_df)))

eco_raw <- analysis_read_tsv(PATH_ECO_TXT, col_types = cols(.default = "c"), locale = locale(encoding = "UTF-8"))
colnames(eco_raw) <- colnames(eco_raw) %>%
  str_replace_all('"', '') %>% str_replace_all('\r', '') %>% str_trim()
eco_raw <- eco_raw %>% mutate(Gene = str_trim(str_replace_all(Gene, '\r', '')))
TF_COL <- "Regulation - transcription factor regulators of gene"
stopifnot(TF_COL %in% colnames(eco_raw))

eco_genes_all <- eco_raw$Gene %>% na.omit() %>% unique()
U1 <- intersect(rna_genes, eco_genes_all)
eco_genes_tfannot <- eco_raw %>%
  filter(!is.na(.data[[TF_COL]]), str_trim(as.character(.data[[TF_COL]])) != "",
         str_trim(as.character(.data[[TF_COL]])) != "NA") %>%
  pull(Gene) %>% na.omit() %>% unique()
U2 <- intersect(rna_genes, eco_genes_tfannot)
message(sprintf("  EcoCyc全收录: %d基因，U1: %d，U2: %d",
                length(eco_genes_all), length(U1), length(U2)))

tf_csv <- analysis_read_csv(PATH_ECO_TF, fileEncoding = "UTF-8-BOM")
stopifnot(all(c("tf_raw","tf_short_manual","genes_all_txt") %in% colnames(tf_csv)))
gene2tf <- tf_csv %>%
  select(tf = tf_short_manual, genes_all_txt) %>%
  mutate(gene_list = str_split(genes_all_txt, "/")) %>%
  unnest(gene_list) %>%
  transmute(gene = str_trim(gene_list), tf = str_trim(tf)) %>%
  filter(gene != "", tf != "", !is.na(gene), !is.na(tf)) %>%
  distinct(gene, tf)
message(sprintf("  gene→tf映射对数: %d", nrow(gene2tf)))

tf_sets_full <- gene2tf %>%
  group_by(tf) %>%
  summarise(genes = list(sort(unique(gene))), M_all = n_distinct(gene), .groups = "drop")
message(sprintf("  TF总数: %d", nrow(tf_sets_full)))
message("[0] 预处理完成 ✓\n")

# ============================================================
# ============================================================

if (RUN_COG) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(scales)
    library(tibble)
  })

  message("[1] 开始COG富集分析（All_responsive）...")
  cog_outdir <- file.path(OUTDIR, "COG_enrichment")
  dir.create(cog_outdir, recursive = TRUE, showWarnings = FALSE)

  fg_genes <- intersect(all_responsive_genes, bg_cog)
  message(sprintf("  前景基因（响应∩背景）: %d", length(fg_genes)))
  if (length(fg_genes) == 0) {
    message("  警告：前景为空，跳过COG分析")
  } else {
    fisher_enrich_one <- function(fg_genes, cate_genes, bg_genes) {
      N <- length(bg_genes)
      M <- sum(bg_genes %in% cate_genes)
      n <- sum(bg_genes %in% fg_genes)
      k <- sum(fg_genes %in% cate_genes)
      if (n == 0 || M == 0) {
        return(tibble(k=k, n=n, M=M, N=N, GeneRatio=NA_real_, bg_rate=NA_real_,
                      OR=NA_real_, log2OR=NA_real_, pvalue=1, hit_genes=""))
      }
      a <- k; b <- n - k; c <- M - k; d <- N - M - (n - k)
      ft <- tryCatch(fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater"),
                     error = function(e) list(p.value=1, estimate=NA_real_))
      OR_val <- if (!is.na(ft$estimate) && ft$estimate > 0) unname(ft$estimate)
                else if (b == 0 || c == 0) Inf else (a*d)/(b*c)
      tibble(k=k, n=n, M=M, N=N, GeneRatio=k/n, bg_rate=M/N,
             OR=OR_val, log2OR=if(is.finite(OR_val) && OR_val>0) log2(OR_val) else NA_real_,
             pvalue=ft$p.value, hit_genes=paste(sort(intersect(fg_genes, cate_genes)), collapse="/"))
    }

    categories <- cate_order_df$COG_Category
    enrich_list <- list()
    for (cate in categories) {
      cate_genes <- gene_cate_long %>% filter(COG_Category == cate) %>% pull(Gene) %>% intersect(bg_cog)
      if (length(cate_genes) == 0) next
      res <- fisher_enrich_one(fg_genes, cate_genes, bg_cog)
      enrich_list[[length(enrich_list)+1]] <- res %>% mutate(COG_Category = cate, .before = 1)
    }
    enrich_df <- bind_rows(enrich_list) %>%
      mutate(qvalue = p.adjust(pvalue, method="BH"),
             neg_log10q = if_else(qvalue > 0, -log10(qvalue), NA_real_),
             is_sig = !is.na(qvalue) & qvalue < Q_CUTOFF) %>%
      left_join(cate_order_df %>% select(COG_Category, COG_Desc, Order), by="COG_Category")

    write_tsv(enrich_df, file.path(cog_outdir, "COG_enrichment_all.tsv"))
    write_tsv(enrich_df %>% filter(is_sig) %>% arrange(qvalue, desc(GeneRatio)),
              file.path(cog_outdir, "COG_enrichment_sig.tsv"))
    message(sprintf("  显著富集条目数（q<%.2f）: %d", Q_CUTOFF, sum(enrich_df$is_sig, na.rm=TRUE)))

    label_order <- cate_order_df %>% arrange(Order) %>%
      mutate(Label = paste0(COG_Category, " — ", COG_Desc)) %>% pull(Label)

    plot_df <- enrich_df %>%
      mutate(Label = factor(paste0(COG_Category, " — ", COG_Desc), levels = rev(label_order)),
             Group = "All_responsive",
             log2OR_plot = if_else(is_sig, log2OR, NA_real_),
             neg_log10q_plot = if_else(is_sig, neg_log10q, NA_real_))

    max_nlq <- max(plot_df$neg_log10q_plot, na.rm=TRUE)
    if (!is.finite(max_nlq) || max_nlq <= 0) max_nlq <- 1
    abs_lor <- max(abs(plot_df$log2OR_plot), na.rm=TRUE)
    if (!is.finite(abs_lor) || abs_lor <= 0) abs_lor <- 1
    COLORS_DIV <- c("#3B4CC0","#88A1D2","#F7F7F7","#F3A982","#B40426")
    NS_PT_COL <- "grey70"; NS_TILE_COL <- "grey90"

    plot_h <- max(5.5, 0.3 * nrow(cate_order_df) + 1.5)
    plot_w <- 6

    # Dotplot
    p_dot <- ggplot() +
      geom_point(data = plot_df %>% filter(!is_sig),
                 aes(x = Group, y = Label, size = GeneRatio),
                 color = NS_PT_COL, alpha = 0.7) +
      geom_point(data = plot_df %>% filter(is_sig),
                 aes(x = Group, y = Label, size = GeneRatio, color = neg_log10q_plot),
                 alpha = 0.95) +
      scale_size_continuous(name = "GeneRatio (k/n)", range = c(1, 9),
                            labels = percent_format(accuracy = 1)) +
      scale_color_gradientn(name = expression(-log[10]~"(q)"), colours = COLORS_DIV,
                            limits = c(0, max_nlq), oob = squish) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            axis.text.y = element_text(size = 9),
            plot.margin = margin(6,6,6,6)) +
      labs(title = "COG enrichment dotplot (All_responsive)",
           subtitle = sprintf("Background: RNA∩COG (N=%d); FDR<%.2f (BH global)",
                              length(bg_cog), Q_CUTOFF))
    ggsave(file.path(cog_outdir, "COG_dotplot.png"), p_dot, width = plot_w, height = plot_h, dpi = FIG_DPI)
    ggsave(file.path(cog_outdir, "COG_dotplot.pdf"), p_dot, width = plot_w, height = plot_h, device = cairo_pdf)

    # Heatmap
    p_heat <- ggplot() +
      geom_tile(data = plot_df %>% filter(!is_sig),
                aes(x = Group, y = Label), fill = NS_TILE_COL, color = "white", linewidth = 0.3) +
      geom_tile(data = plot_df %>% filter(is_sig),
                aes(x = Group, y = Label, fill = log2OR_plot),
                color = "white", linewidth = 0.3) +
      scale_fill_gradientn(name = "log2(OR)", colours = COLORS_DIV,
                           limits = c(-abs_lor, abs_lor), oob = squish, na.value = NS_TILE_COL) +
      scale_x_discrete(expand = c(0,0)) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            axis.text.y = element_text(size = 9),
            plot.margin = margin(6,6,6,6)) +
      labs(title = "COG enrichment heatmap (All_responsive)",
           subtitle = sprintf("Background: RNA∩COG (N=%d); FDR<%.2f (BH global)",
                              length(bg_cog), Q_CUTOFF))
    ggsave(file.path(cog_outdir, "COG_heatmap.png"), p_heat, width = plot_w, height = plot_h, dpi = FIG_DPI)
    ggsave(file.path(cog_outdir, "COG_heatmap.pdf"), p_heat, width = plot_w, height = plot_h, device = cairo_pdf)

    message(sprintf("[1] COG富集完成 ✓ → %s\n", cog_outdir))
  }
}

# ============================================================
# 4. TF enrichment (All_responsive): V6.2 complete ORA
# Replace the original section 4 with this entire if (RUN_TF) block.
# The original COG configuration, preparation and analysis above are retained.
# This module uses explicit filtered/response gene lists, not U1/U2 or clusters.
# ============================================================
if (RUN_TF) {
  tf_v62_results <- local({
    # Edit these two paths if your exported gene-list filenames differ.
    tf_cfg <- list(
      background_csv = "filtered_genes.csv",
      response_csv = "responsive_genes.csv",
      annotation_csv = PATH_ECO_TF,
      output_dir = file.path(OUTDIR, "TF_enrichment", "tf"),
      min_targets = M_MIN_TF,
      fdr_threshold = Q_CUTOFF,
      make_plots = TRUE
    )
    # Preserve the original explicitly configured engineered-gene alias.
    # Other engineered fragments retain their own identities.
    aliases <- setNames(tolower(trimws(DNAA_IN_CLUSTER)),
                        tolower(trimws(DNAA_IN_RNA)))
    clean_text <- function(x) {
      x <- trimws(as.character(x))
      x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
      x
    }
    gene_key <- function(x) {
      z <- tolower(clean_text(x))
      hit <- !is.na(z) & z %in% names(aliases)
      z[hit] <- unname(aliases[z[hit]])
      z
    }
    safe_fraction <- function(a, b) if (b > 0) a / b else NA_real_
    read_tf_csv <- function(path) {
      if (!analysis_file_exists(path)) stop("TF input file not found: ", path)
      d <- analysis_read_csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                    fileEncoding = "UTF-8-BOM")
      names(d) <- sub("^\ufeff", "", names(d))
      d
    }
    gene_table <- function(d, label) {
      # Some original tables contain both Gene and gene: prefer exact 'gene'.
      nm <- intersect(c("gene", "Gene", "gene_key"), names(d))
      if (!length(nm)) stop(label, " requires gene, Gene or gene_key.")
      out <- data.frame(original_gene = clean_text(d[[nm[1]]]),
                        gene_key = gene_key(d[[nm[1]]]),
                        stringsAsFactors = FALSE)
      if (anyNA(out$gene_key)) stop(label, " contains empty/missing gene identifiers.")
      by_key <- split(out$original_gene, out$gene_key)
      collisions <- names(Filter(function(z) length(unique(z)) > 1L, by_key))
      if (length(collisions)) {
        stop(label, ": distinct original IDs collapse to the same canonical ID: ",
             paste(collisions, collapse = ", "),
             ". Review aliases before analysis.")
      }
      unique(out)
    }
    require_cols <- function(d, cols, label) {
      absent <- setdiff(cols, names(d))
      if (length(absent)) stop(label, " missing columns: ", paste(absent, collapse = ", "))
    }
    # Verbatim V6.2 statistical helper; surrounding IO/metadata are adapted
    # to the original standalone All_responsive analysis.
    v6_ora_from_map_complete <- function(foreground, background, term_map, min_term_size, gene_label_map = NULL) {
      foreground <- unique(intersect(foreground, background))
      background <- unique(background)
      map <- unique(term_map[term_map$gene_key %in% background & !is.na(term_map$term), , drop = FALSE])
      if (!length(background) || !nrow(map)) return(data.frame())
      terms <- sort(unique(map$term))
      N <- length(background)
      n <- length(foreground)
      ans <- lapply(terms, function(term) {
        targets <- unique(map$gene_key[map$term == term])
        K <- length(targets)
        if (K < min_term_size) return(NULL)
        hit_keys <- sort(unique(intersect(foreground, targets)))
        a <- length(hit_keys)
        b <- n - a
        c0 <- K - a
        d <- N - K - b
        if (any(c(a,b,c0,d) < 0)) return(NULL)
        hit_labels <- hit_keys
        if (!is.null(gene_label_map) && length(gene_label_map) && length(hit_keys)) {
          mapped <- unname(gene_label_map[hit_keys])
          use_mapped <- !is.na(mapped) & nzchar(mapped)
          hit_labels[use_mapped] <- mapped[use_mapped]
        }
        tab <- matrix(c(a,b,c0,d), nrow = 2, byrow = TRUE)
        ft_enrich <- stats::fisher.test(tab, alternative = "greater")
        ft_ci <- stats::fisher.test(tab, alternative = "two.sided", conf.level = 0.95)
        data.frame(
          term = term,
          foreground_n = n,
          background_n = N,
          overlap_n = a,
          background_term_n = K,
          hit_genes = paste(hit_labels, collapse = ";"),
          hit_gene_keys = paste(hit_keys, collapse = ";"),
          gene_ratio = safe_fraction(a, n),
          background_ratio = safe_fraction(K, N),
          odds_ratio = unname(ft_enrich$estimate),
          CI_low_95 = unname(ft_ci$conf.int[1]),
          CI_high_95 = unname(ft_ci$conf.int[2]),
          p_value = ft_enrich$p.value,
          stringsAsFactors = FALSE
        )
      })
      ans <- ans[!vapply(ans, is.null, logical(1))]
      if (!length(ans)) return(data.frame())
      do.call(rbind, ans)
    }

    bg_df <- gene_table(read_tf_csv(tf_cfg$background_csv), "Filtered background")
    fg_raw <- read_tf_csv(tf_cfg$response_csv)
    # Accept a responsive-only table or an all-genes statistics table.
    if ("direction" %in% names(fg_raw)) {
      fg_raw <- fg_raw[!is.na(fg_raw$direction) &
                         fg_raw$direction %in% c("Up", "Down"), , drop = FALSE]
    }
    fg_df <- gene_table(fg_raw, "Responsive genes")
    bg <- sort(unique(bg_df$gene_key))
    fg <- sort(unique(fg_df$gene_key))
    if (!length(fg) || !length(bg) || length(fg) >= length(bg)) {
      stop("TF analysis needs a nonempty foreground and a larger background.")
    }
    outside <- setdiff(fg, bg)
    if (length(outside)) {
      stop("Responsive genes outside the filtered background: ",
           paste(outside, collapse = ", "), ". No genes were silently dropped.")
    }
    if (!is.finite(tf_cfg$min_targets) || tf_cfg$min_targets < 1 ||
        tf_cfg$min_targets != as.integer(tf_cfg$min_targets)) {
      stop("M_MIN_TF must be a positive integer.")
    }
    if (!is.finite(tf_cfg$fdr_threshold) ||
        tf_cfg$fdr_threshold <= 0 || tf_cfg$fdr_threshold >= 1) {
      stop("Q_CUTOFF must be between 0 and 1.")
    }

    ann <- read_tf_csv(tf_cfg$annotation_csv)
    require_cols(ann, c("tf_short_manual", "genes_all_txt"), "Manual EcoCyc TF table")
    raw_rows <- lapply(seq_len(nrow(ann)), function(i) {
      tf <- clean_text(ann$tf_short_manual[i])
      txt <- clean_text(ann$genes_all_txt[i])
      if (is.na(tf) || is.na(txt)) return(NULL)
      tokens <- unlist(strsplit(txt, "/", fixed = TRUE), use.names = FALSE)
      data.frame(term = tf, original_target = clean_text(tokens),
                 gene_key = gene_key(tokens), stringsAsFactors = FALSE)
    })
    raw_rows <- Filter(Negate(is.null), raw_rows)
    if (!length(raw_rows)) stop("No TF-target relationships found.")
    raw_map <- do.call(rbind, raw_rows)
    raw_map <- raw_map[!is.na(raw_map$gene_key), , drop = FALSE]
    term_map <- unique(raw_map[, c("gene_key", "term"), drop = FALSE])
    label_map <- setNames(bg_df$original_gene, bg_df$gene_key)

    out <- v6_ora_from_map_complete(fg, bg, term_map,
                                     tf_cfg$min_targets, label_map)
    if (!nrow(out)) stop("No TF has enough targets in the filtered background.")
    # Exactly one BH family: all eligible TFs, including zero foreground hits.
    out$FDR <- p.adjust(out$p_value, method = "BH")
    out$FDR_family <- "dnaA_RNA|Response|TF|absolute"
    out$FDR_method <- "Benjamini-Hochberg"
    out$n_tests_in_family <- nrow(out)
    out$significant <- is.finite(out$FDR) & out$FDR < tf_cfg$fdr_threshold
    out$FE <- out$gene_ratio / out$background_ratio
    out$log2FE <- log2(out$FE)
    out$neg_log10FDR <- -log10(pmax(out$FDR, .Machine$double.xmin))
    out$Group <- "Response"
    out$analysis_level <- "individual"
    out$analysis_id <- "dnaA_RNA"
    out$group_type <- "Response"
    out$database <- "TF"
    out$enrichment_type <- "absolute"
    out$foreground_definition <- "Direction-concordant responsive genes from original Step5"
    out$background_definition <- "All genes retained after low-expression and pre-fluctuation filters"
    out$comparison_definition <- "Responsive genes vs filtered background minus responsive genes"
    # TF membership enrichment does not depend on peak/mean effect selection.
    out$effect_metric <- "not_used_membership_only"
    out$effect_definition <- "Original response membership; no effect-size re-screening"
    out <- out[order(out$FDR, out$p_value, out$term), , drop = FALSE]
    rownames(out) <- NULL

    # Keep the complete table intact even for OR=0/Inf, infinite CI or zero hits.
    dir.create(tf_cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
    write_csv <- function(x, name) {
      write.csv(x, file.path(tf_cfg$output_dir, name),
                row.names = FALSE, na = "", fileEncoding = "UTF-8")
    }
    write_csv(out, "TF__Response__absolute.csv")
    write_csv(out[out$significant, , drop = FALSE],
              "TF__Response__absolute__significant.csv")

    # Legacy columns support downstream inspection using the original names.
    legacy <- data.frame(tf = out$term, k = out$overlap_n,
      n = out$foreground_n, M = out$background_term_n, N = out$background_n,
      GeneRatio = out$gene_ratio, bg_rate = out$background_ratio,
      FE = out$FE, log2FE = out$log2FE, OR = out$odds_ratio,
      CI_low_95 = out$CI_low_95, CI_high_95 = out$CI_high_95,
      pvalue = out$p_value, qvalue = out$FDR, neg_log10q = out$neg_log10FDR,
      is_sig = out$significant, geneID = out$hit_genes,
      universe = "Filtered_background", Group = "All_responsive")
    write.table(legacy, file.path(tf_cfg$output_dir, "TF_enrich_all.tsv"),
                sep = "\t", quote = TRUE, row.names = FALSE, na = "",
                fileEncoding = "UTF-8")
    write.table(legacy[legacy$is_sig, , drop = FALSE],
                file.path(tf_cfg$output_dir, "TF_enrich_sig.tsv"),
                sep = "\t", quote = TRUE, row.names = FALSE, na = "",
                fileEncoding = "UTF-8")

    bg_df$is_responsive <- bg_df$gene_key %in% fg
    bg_df$has_TF_annotation <- bg_df$gene_key %in% term_map$gene_key
    bg_df$has_eligible_TF_annotation <- bg_df$gene_key %in%
      term_map$gene_key[term_map$term %in% out$term]
    write_csv(bg_df, "TF_background_gene_list.csv")
    write_csv(fg_df, "TF_foreground_gene_list.csv")
    raw_map$in_background <- raw_map$gene_key %in% bg
    raw_map$eligible_TF <- raw_map$term %in% out$term
    write_csv(unique(raw_map), "TF_target_mapping.csv")
    write_csv(data.frame(original_key = names(aliases), canonical_key = unname(aliases)),
              "TF_gene_aliases.csv")
    annotated_bg <- intersect(bg, term_map$gene_key)
    annotated_fg <- intersect(fg, term_map$gene_key)
    write_csv(data.frame(foreground_total = length(fg),
      foreground_annotated = length(annotated_fg),
      foreground_annotation_fraction = length(annotated_fg) / length(fg),
      background_total = length(bg), background_annotated = length(annotated_bg),
      background_annotation_fraction = length(annotated_bg) / length(bg)),
      "TF_annotation_coverage.csv")


    if (isTRUE(tf_cfg$make_plots) && any(out$significant)) {
      if (!requireNamespace("ggplot2", quietly = TRUE) ||
          !requireNamespace("scales", quietly = TRUE)) {
        warning("TF statistics saved; install ggplot2/scales to generate TF plots.")
      } else {
        pd <- out[out$significant, , drop = FALSE]
        pd$term <- factor(pd$term, levels = pd$term)
        pd$Group <- "All_responsive"
        palette <- c("#3B4CC0", "#88A1D2", "#F7F7F7", "#F3A982", "#B40426")
        subtitle <- sprintf("Filtered background N=%d; response n=%d; BH across %d TFs",
                            length(bg), length(fg), nrow(out))
        p_dot <- ggplot2::ggplot(pd, ggplot2::aes(x = term, y = Group)) +
          ggplot2::geom_point(ggplot2::aes(size = gene_ratio, colour = neg_log10FDR)) +
          ggplot2::scale_size_continuous(name = "GeneRatio (k/n)", range = c(2, 10),
                                         labels = scales::percent_format(accuracy = 1)) +
          ggplot2::scale_colour_gradientn(colours = palette, name = "-log10(FDR)") +
          ggplot2::theme_bw(base_size = 12) +
          ggplot2::theme(axis.title = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::labs(title = "TF enrichment (All_responsive)", subtitle = subtitle)
        p_heat <- ggplot2::ggplot(pd, ggplot2::aes(x = term, y = Group,
                                                 fill = neg_log10FDR)) +
          ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
          ggplot2::scale_fill_gradientn(colours = palette, name = "-log10(FDR)") +
          ggplot2::theme_bw(base_size = 12) +
          ggplot2::theme(axis.title = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
          ggplot2::labs(title = "TF enrichment heatmap (All_responsive)",
                        subtitle = subtitle)
        plot_width <- max(8, 0.5 * nrow(pd))
        for (kind in c("dotplot", "heatmap")) {
          p <- if (kind == "dotplot") p_dot else p_heat
          ggplot2::ggsave(file.path(tf_cfg$output_dir, paste0("TF_", kind, ".pdf")),
                          p, width = plot_width, height = 4)
          ggplot2::ggsave(file.path(tf_cfg$output_dir, paste0("TF_", kind, ".png")),
                          p, width = plot_width, height = 4, dpi = FIG_DPI)
        }
      }
    }
    message("TF V6.2 complete ORA: ", sum(out$significant), " significant TFs; ",
            normalizePath(tf_cfg$output_dir, winslash = "/"))
    list(results = out, config = tf_cfg)
  })
}

# ============================================================
# 5. Publication TF figures integrated from publication_V1.R
# This is a display/export module; it never recomputes P values or BH FDR.
# Only this configuration needs editing to change the main figure selection.
# ============================================================
TF_PUBLICATION_CONFIG <- list(
  enabled = TRUE,
  fdr_threshold = Q_CUTOFF,
  top_n = 3L,
  always_add_prespecified_tf = TRUE,
  prespecified_tf = "DnaA",
  evidence_colour_limits = c(1, 15),
  panel_width_cm = 2.54,
  panel_height_cm = 2.54,
  font_family = "Arial",
  font_size_pt = 7,
  output_all_significant_pdfs = FALSE
)

run_tf_publication_figures <- function(dat, source_csv, output_dir, settings = list()) {
  CONFIG <- list(
    input_csv = source_csv,
    output_dir = output_dir,
    top_n = 3L,
    fdr_threshold = 0.05,
    prespecified_tf = "DnaA",
    always_add_prespecified_tf = TRUE,
    evidence_palette = c("#67A9CF", "#EF8A62", "#B40426"),
    evidence_colour_limits = c(1, 15),
    panel_width_cm = 2.54,
    panel_height_cm = 2.54,
    frame_line_pt = 0.25,
    font_family = "Arial",
    font_size_pt = 7,
    pdf_page_width_cm = 9.5,
    pdf_page_height_cm = 5.1,
    bar_metric = "log2FE", # alternative: "neg_log10FDR"
    dot_cex_range = c(0.85, 2.0),
    forest_point_cex = 1.05,
    embed_vertical_colorbar = TRUE,
    colorbar_width_cm = 2.2,
    colorbar_height_cm = 5.5,
    output_all_significant_pdfs = FALSE
  )
  CONFIG <- utils::modifyList(CONFIG, settings)

  # Font fallback applies only where Cairo is unavailable.
  if (!isTRUE(capabilities("cairo"))) CONFIG$font_family <- "Helvetica"
  open_publication_pdf <- function(path, width, height) {
    if (isTRUE(capabilities("cairo"))) {
      grDevices::cairo_pdf(filename = path, width = width, height = height,
        family = CONFIG$font_family, pointsize = CONFIG$font_size_pt, onefile = TRUE)
    } else {
      grDevices::pdf(file = path, width = width, height = height,
        family = CONFIG$font_family, pointsize = CONFIG$font_size_pt,
        onefile = TRUE, useDingbats = FALSE)
    }
  }
  dir.create(CONFIG$output_dir, recursive = TRUE, showWarnings = FALSE)

  required <- c(
    "term", "foreground_n", "background_n", "overlap_n",
    "background_term_n", "gene_ratio", "background_ratio", "odds_ratio",
    "CI_low_95", "CI_high_95", "p_value", "FDR", "significant",
    "hit_genes", "FDR_family", "n_tests_in_family"
  )
  # dat is passed directly from the freshly computed V6.2 statistics.
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("Missing required column(s): ", paste(missing, collapse = ", "))

  as_flag <- function(x) {
    if (is.logical(x)) return(x)
    tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
  }
  num_cols <- c(
    "foreground_n", "background_n", "overlap_n", "background_term_n",
    "gene_ratio", "background_ratio", "odds_ratio", "CI_low_95",
    "CI_high_95", "p_value", "FDR", "n_tests_in_family"
  )
  for (nm in num_cols) dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
  dat$significant <- as_flag(dat$significant)
  dat$FE <- dat$gene_ratio / dat$background_ratio
  dat$log2FE <- log2(dat$FE)
  dat$neg_log10FDR <- -log10(pmax(dat$FDR, .Machine$double.xmin))

  valid <- is.finite(dat$FDR) & dat$FDR >= 0 & dat$FDR <= 1 &
    is.finite(dat$FE) & dat$FE > 0 & is.finite(dat$log2FE) &
    is.finite(dat$odds_ratio) & dat$odds_ratio > 0 &
    is.finite(dat$CI_low_95) & dat$CI_low_95 > 0 &
    is.finite(dat$CI_high_95) & dat$CI_high_95 > 0

  stat_sig <- !is.na(dat$significant) & dat$significant &
    is.finite(dat$FDR) & dat$FDR < CONFIG$fdr_threshold
  dat$statistically_significant <- stat_sig
  dat$plot_eligible <- valid
  dat$plot_exclusion_reason <- vapply(seq_len(nrow(dat)), function(i) {
    reason <- character()
    if (!is.finite(dat$FDR[i]) || dat$FDR[i] < 0 || dat$FDR[i] > 1)
      reason <- c(reason, "invalid_FDR")
    if (!is.finite(dat$FE[i]) || dat$FE[i] <= 0 || !is.finite(dat$log2FE[i]))
      reason <- c(reason, "invalid_FE")
    if (!is.finite(dat$odds_ratio[i]) || dat$odds_ratio[i] <= 0)
      reason <- c(reason, "nonpositive_or_nonfinite_OR")
    if (!is.finite(dat$CI_low_95[i]) || dat$CI_low_95[i] <= 0 ||
        !is.finite(dat$CI_high_95[i]) || dat$CI_high_95[i] <= 0)
      reason <- c(reason, "nonpositive_or_nonfinite_CI")
    paste(reason, collapse = ";")
  }, character(1))
  write.csv(dat, file.path(CONFIG$output_dir, "TF_plot_data.csv"),
            row.names = FALSE, na = "")
  write.csv(dat[stat_sig, , drop = FALSE],
            file.path(CONFIG$output_dir, "ALL_statistically_significant_TFs.csv"),
            row.names = FALSE, na = "")
  sig <- dat[stat_sig & valid, , drop = FALSE]
  sig <- sig[order(sig$FDR, -sig$log2FE, sig$term), , drop = FALSE]
  if (!nrow(sig)) {
    status <- data.frame(status = "skipped_no_plot_eligible_significant_TFs",
      n_tests = nrow(dat), n_statistically_significant = sum(stat_sig),
      n_plot_eligible_significant = 0L)

    message("Publication figures skipped; full statistics and plot data are saved.")
    return(invisible(list(status = status$status, selected_terms = character())))
  }
  if (anyDuplicated(sig$term)) stop("Significant TF names are not unique.")
  if (length(CONFIG$top_n) != 1L || !is.finite(CONFIG$top_n) ||
      CONFIG$top_n < 1L || CONFIG$top_n != as.integer(CONFIG$top_n))
    stop("TF_PUBLICATION_CONFIG$top_n must be a positive integer.")
  top_n <- min(as.integer(CONFIG$top_n), nrow(sig))
  top_terms <- sig$term[seq_len(top_n)]
  selected_terms <- top_terms
  # DnaA must be significant and plot-eligible. Match its spelling case-insensitively.
  prespecified_match <- sig$term[
    tolower(sig$term) == tolower(CONFIG$prespecified_tf)]
  if (length(prespecified_match) > 1L)
    stop("Multiple TF labels match the prespecified TF.")
  if (isTRUE(CONFIG$always_add_prespecified_tf) && length(prespecified_match)) {
    selected_terms <- c(setdiff(selected_terms, prespecified_match), prespecified_match)
  }

  # Continuous CIE-Lab mapping with fixed configured limits, independent of top_n.
  colour_limits <- as.numeric(CONFIG$evidence_colour_limits)
  if (length(colour_limits) != 2L || any(!is.finite(colour_limits)) ||
      colour_limits[2] <= colour_limits[1]) {
    stop("evidence_colour_limits must contain two increasing finite values.")
  }
  colour_ramp <- grDevices::colorRamp(CONFIG$evidence_palette, space = "Lab")
  map_colour <- function(x) {
    u <- pmin(1, pmax(0, (x - colour_limits[1]) / diff(colour_limits)))
    m <- round(colour_ramp(u))
    grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
  }

  sig$plot_colour_HEX <- map_colour(sig$neg_log10FDR)
  rgb_values <- grDevices::col2rgb(sig$plot_colour_HEX)
  sig$plot_colour_R <- as.integer(rgb_values["red", ])
  sig$plot_colour_G <- as.integer(rgb_values["green", ])
  sig$plot_colour_B <- as.integer(rgb_values["blue", ])
  sig$plot_colour_RGB <- sprintf(
    "RGB(%d,%d,%d)", sig$plot_colour_R, sig$plot_colour_G, sig$plot_colour_B
  )
  sig$rank_by_FDR <- seq_len(nrow(sig))
  sig$selected_by_top_n <- sig$term %in% top_terms
  sig$is_prespecified_DnaA <- tolower(sig$term) == tolower(CONFIG$prespecified_tf)
  sig$selected_for_main_plot <- sig$term %in% selected_terms
  sig$selection_reason <- ifelse(
    sig$selected_by_top_n, paste0("Top_", top_n, "_by_FDR"),
    ifelse(sig$is_prespecified_DnaA, "Prespecified_perturbation_target", "Not_selected")
  )
  sig$display_order <- match(sig$term, selected_terms)
  selected <- sig[match(selected_terms, sig$term), , drop = FALSE]

  # R vector devices use approximately 0.75 pt per lwd unit.
  frame_lwd <- CONFIG$frame_line_pt / 0.75
  open_pdf <- function(path, page_height_cm = CONFIG$pdf_page_height_cm,
                       panel_height_cm = CONFIG$panel_height_cm) {
    open_publication_pdf(path, CONFIG$pdf_page_width_cm / 2.54,
                         page_height_cm / 2.54)
    graphics::par(
      mar = c(3.2, 5.2, 0.7, 8.1), mgp = c(1.75, 0.45, 0),
      tcl = -0.22, las = 1, family = CONFIG$font_family,
      ps = CONFIG$font_size_pt, xaxs = "i", yaxs = "i"
    )
    # Set pin after margins: changing mar recomputes pin in base R.
    graphics::par(pin = c(CONFIG$panel_width_cm, panel_height_cm) / 2.54)
  }
  close_pdf <- function() grDevices::dev.off()

  legend_breaks <- function() unique(round(c(
    colour_limits[1], mean(colour_limits), colour_limits[2]
  ), 2))
  draw_evidence_legend <- function(position = "right", inset_x = -0.88, inset_y = 0) {
    br <- legend_breaks()
    graphics::legend(
      position, inset = c(inset_x, inset_y), xpd = NA,
      title = expression(-log[10](FDR)), legend = format(br, trim = TRUE),
      pch = 21, pt.bg = map_colour(br), col = "black",
      pt.cex = 1.15, pt.lwd = frame_lwd, bty = "n", cex = 0.82,
      y.intersp = 1, x.intersp = 0.7
    )
  }
  dot_cex <- function(hits) {
    z <- sqrt(hits)
    CONFIG$dot_cex_range[1] + (z - min(z)) /
      max(diff(range(z)), .Machine$double.eps) * diff(CONFIG$dot_cex_range)
  }
  draw_hits_legend <- function(hits, position = "topright", inset_x = -0.88, inset_y = 0) {
    br <- unique(round(stats::quantile(hits, c(0, 0.5, 1), names = FALSE)))
    all_z <- sqrt(hits)
    cex <- CONFIG$dot_cex_range[1] + (sqrt(br) - min(all_z)) /
      max(diff(range(all_z)), .Machine$double.eps) * diff(CONFIG$dot_cex_range)
    graphics::legend(
      position, inset = c(inset_x, inset_y), xpd = NA,
      title = "Hit genes", legend = br, pch = 21,
      pt.bg = "grey90", col = "black", pt.cex = cex,
      pt.lwd = frame_lwd, bty = "n", cex = 0.82,
      y.intersp = 1.3, x.intersp = 0.7
    )
  }

  plot_forest <- function(d, path, page_height_cm = CONFIG$pdf_page_height_cm,
                          panel_height_cm = CONFIG$panel_height_cm) {
    open_pdf(path, page_height_cm, panel_height_cm)
    on.exit(close_pdf(), add = TRUE)
    n <- nrow(d); y <- rev(seq_len(n))
    xmin <- min(1, d$CI_low_95) / 1.18; xmax <- max(d$CI_high_95) * 1.18
    graphics::plot(
      d$odds_ratio, y, type = "n", log = "x", xlim = c(xmin, xmax),
      ylim = c(0.5, n + 0.5), xlab = "Odds ratio", ylab = "",
      yaxt = "n", bty = "n", axes = FALSE
    )
    graphics::abline(v = 1, lty = 2, lwd = frame_lwd, col = "grey45")
    graphics::segments(d$CI_low_95, y, d$CI_high_95, y, lwd = 0.65, col = "grey35")
    graphics::points(
      d$odds_ratio, y, pch = 21, bg = d$plot_colour_HEX,
      col = "black", cex = CONFIG$forest_point_cex, lwd = frame_lwd
    )
    pool <- c(0.25, 0.5, 1, 2, 5, 10, 20, 50, 100)
    ticks <- pool[pool >= xmin & pool <= xmax]
    graphics::axis(1, at = ticks, labels = ticks, lwd = 0, lwd.ticks = frame_lwd)
    graphics::axis(2, at = y, labels = d$term, tick = FALSE, lwd = 0, cex.axis = 0.95)
    graphics::box(lwd = frame_lwd)
    graphics::mtext("TF", side = 2, line = 3.55, las = 0, cex = 0.95)
    graphics::text(
      rep(xmax * 1.08, n), y, paste0("n=", d$overlap_n),
      xpd = NA, adj = 0, cex = 0.78
    )
    if (isTRUE(CONFIG$embed_vertical_colorbar)) draw_vertical_colorbar_in_margin()
  }

  plot_bar <- function(d, path, page_height_cm = CONFIG$pdf_page_height_cm,
                       panel_height_cm = CONFIG$panel_height_cm) {
    open_pdf(path, page_height_cm, panel_height_cm)
    on.exit(close_pdf(), add = TRUE)
    if (identical(CONFIG$bar_metric, "log2FE")) {
      value <- d$log2FE; xlabel <- expression(log[2](fold~enrichment))
    } else if (identical(CONFIG$bar_metric, "neg_log10FDR")) {
      value <- d$neg_log10FDR; xlabel <- expression(-log[10](FDR))
    } else stop("bar_metric must be 'log2FE' or 'neg_log10FDR'.")
    xmax <- max(value) * 1.30
    bp <- graphics::barplot(
      rev(value), names.arg = rev(d$term), horiz = TRUE,
      col = rev(d$plot_colour_HEX), border = "black", lwd = frame_lwd,
      space = 0.35, axes = FALSE, xlab = xlabel,
      xlim = c(0, xmax), cex.names = 0.95, las = 1
    )
    graphics::axis(1, lwd = 0, lwd.ticks = frame_lwd)
    graphics::box(lwd = frame_lwd)
    graphics::text(
      rev(value) + xmax * 0.025, bp, paste0("n=", rev(d$overlap_n)),
      xpd = NA, adj = 0, cex = 0.78
    )
    if (isTRUE(CONFIG$embed_vertical_colorbar)) draw_vertical_colorbar_in_margin()
  }

  plot_dot <- function(d, path, page_height_cm = CONFIG$pdf_page_height_cm,
                       panel_height_cm = CONFIG$panel_height_cm) {
    open_pdf(path, page_height_cm, panel_height_cm)
    on.exit(close_pdf(), add = TRUE)
    n <- nrow(d); y <- rev(seq_len(n)); cex <- dot_cex(d$overlap_n)
    xmin <- min(0, min(d$log2FE)); xmax <- max(d$log2FE) * 1.15
    graphics::plot(
      d$log2FE, y, type = "n", xlim = c(xmin, xmax), ylim = c(0.5, n + 0.5),
      xlab = expression(log[2](fold~enrichment)), ylab = "",
      yaxt = "n", bty = "n", axes = FALSE
    )
    graphics::abline(v = 0, lty = 2, lwd = frame_lwd, col = "grey45")
    graphics::segments(0, y, d$log2FE, y, col = "grey80", lwd = frame_lwd)
    graphics::points(
      d$log2FE, y, pch = 21, bg = d$plot_colour_HEX,
      col = "black", cex = cex, lwd = frame_lwd
    )
    graphics::axis(1, lwd = 0, lwd.ticks = frame_lwd)
    graphics::axis(2, at = y, labels = d$term, tick = FALSE, lwd = 0, cex.axis = 0.95)
    graphics::box(lwd = frame_lwd)
    graphics::mtext("TF", side = 2, line = 3.55, las = 0, cex = 0.95)
    draw_hits_legend(d$overlap_n, position = "topright", inset_x = -0.60)
    if (isTRUE(CONFIG$embed_vertical_colorbar)) draw_vertical_colorbar_in_margin()
  }

  # Vertical continuous evidence colourbar. The same core is used inside every
  # plot PDF and in the standalone legend PDF.
  draw_vertical_colorbar_core <- function() {
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = colour_limits, xaxs = "i", yaxs = "i")
    gradient <- grDevices::as.raster(matrix(
      rev(map_colour(seq(colour_limits[1], colour_limits[2], length.out = 512))),
      ncol = 1
    ))
    graphics::rasterImage(
      gradient, 0.10, colour_limits[1], 0.48, colour_limits[2],
      interpolate = TRUE
    )
    graphics::rect(
      0.10, colour_limits[1], 0.48, colour_limits[2],
      border = "black", lwd = frame_lwd
    )
    ticks <- c(1, 5, 10, 15)
    ticks <- ticks[ticks >= colour_limits[1] & ticks <= colour_limits[2]]
    graphics::axis(
      4, at = ticks, labels = ticks, pos = 0.48,
      las = 1, lwd = 0, lwd.ticks = frame_lwd, cex.axis = 0.82
    )
    graphics::mtext(
      expression(-log[10](BH-FDR)), side = 2, line = -0.2,
      las = 0, cex = 0.78
    )
  }

  draw_vertical_colorbar_in_margin <- function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(
      fig = c(0.86, 0.96, 0.20, 0.82), new = TRUE,
      mar = c(0.2, 0.2, 0.2, 1.9), family = CONFIG$font_family,
      ps = CONFIG$font_size_pt, xaxs = "i", yaxs = "i"
    )
    draw_vertical_colorbar_core()
  }

  plot_shared_colorbar <- function(path) {
    open_publication_pdf(path, CONFIG$colorbar_width_cm / 2.54,
                         CONFIG$colorbar_height_cm / 2.54)
    on.exit(grDevices::dev.off(), add = TRUE)
    graphics::par(
      mar = c(0.7, 1.5, 0.7, 2.7), family = CONFIG$font_family,
      ps = CONFIG$font_size_pt, xaxs = "i", yaxs = "i"
    )
    draw_vertical_colorbar_core()
  }

  selection_tag <- paste0("Top", top_n, "_plus_", CONFIG$prespecified_tf)
  forest_pdf <- file.path(CONFIG$output_dir, paste0("01_TF_forest_", selection_tag, ".pdf"))
  bar_pdf <- file.path(CONFIG$output_dir, paste0("02_TF_horizontal_bar_", selection_tag, ".pdf"))
  dot_pdf <- file.path(CONFIG$output_dir, paste0("03_TF_horizontal_dotplot_", selection_tag, ".pdf"))
  colorbar_pdf <- file.path(CONFIG$output_dir, "04_shared_continuous_FDR_colorbar.pdf")
  plot_forest(selected, forest_pdf)
  plot_bar(selected, bar_pdf)
  plot_dot(selected, dot_pdf)
  plot_shared_colorbar(colorbar_pdf)

  if (isTRUE(CONFIG$output_all_significant_pdfs)) {
    ph <- max(2.54, 0.32 * nrow(sig)); page_h <- ph + 3.0
    plot_forest(sig, file.path(CONFIG$output_dir, "05_TF_forest_ALL_significant.pdf"), page_h, ph)
    plot_bar(sig, file.path(CONFIG$output_dir, "06_TF_horizontal_bar_ALL_significant.pdf"), page_h, ph)
    plot_dot(sig, file.path(CONFIG$output_dir, "07_TF_horizontal_dotplot_ALL_significant.pdf"), page_h, ph)
  }

  all_columns <- c(
    "rank_by_FDR", "term", "foreground_n", "background_n", "overlap_n",
    "background_term_n", "gene_ratio", "background_ratio", "FE", "log2FE",
    "odds_ratio", "CI_low_95", "CI_high_95", "p_value", "FDR",
    "neg_log10FDR", "hit_genes", "FDR_family", "n_tests_in_family",
    "selected_by_top_n", "is_prespecified_DnaA", "selected_for_main_plot",
    "selection_reason", "display_order", "plot_colour_HEX", "plot_colour_R",
    "plot_colour_G", "plot_colour_B", "plot_colour_RGB"
  )
  write.csv(
    sig[, all_columns, drop = FALSE],
    file.path(CONFIG$output_dir, "ALL_significant_TFs_plot_parameters_and_RGB.csv"),
    row.names = FALSE, na = ""
  )
  selected_out <- sig[match(selected_terms, sig$term), all_columns, drop = FALSE]
  write.csv(
    selected_out,
    file.path(CONFIG$output_dir, paste0("SELECTED_", selection_tag, "_plot_parameters.csv")),
    row.names = FALSE, na = ""
  )

  make_forest_table <- function(x) data.frame(
    TF = x$term, OR = x$odds_ratio, CI_low_95 = x$CI_low_95,
    CI_high_95 = x$CI_high_95, NegLog10FDR = x$neg_log10FDR,
    Hits = x$overlap_n, Colour_HEX = x$plot_colour_HEX,
    Colour_R = x$plot_colour_R, Colour_G = x$plot_colour_G,
    Colour_B = x$plot_colour_B, Colour_RGB = x$plot_colour_RGB,
    Selection_reason = x$selection_reason, stringsAsFactors = FALSE
  )
  make_bar_table <- function(x) data.frame(
    TF = x$term, Bar_log2FE = x$log2FE, NegLog10FDR = x$neg_log10FDR,
    Hits = x$overlap_n, Colour_HEX = x$plot_colour_HEX,
    Colour_R = x$plot_colour_R, Colour_G = x$plot_colour_G,
    Colour_B = x$plot_colour_B, Colour_RGB = x$plot_colour_RGB,
    Selection_reason = x$selection_reason, stringsAsFactors = FALSE
  )
  make_dot_table <- function(x) data.frame(
    TF = x$term, X_log2FE = x$log2FE,
    Y_display_order = seq_len(nrow(x)), Bubble_size_Hits = x$overlap_n,
    Colour_value_NegLog10FDR = x$neg_log10FDR,
    Colour_HEX = x$plot_colour_HEX, Colour_R = x$plot_colour_R,
    Colour_G = x$plot_colour_G, Colour_B = x$plot_colour_B,
    Colour_RGB = x$plot_colour_RGB, Selection_reason = x$selection_reason,
    stringsAsFactors = FALSE
  )

  # GraphPad tables: selected main set and all significant TFs.
  write.csv(make_forest_table(selected), file.path(CONFIG$output_dir, paste0("GraphPad_forest_", selection_tag, ".csv")), row.names = FALSE)
  write.csv(make_bar_table(selected), file.path(CONFIG$output_dir, paste0("GraphPad_bar_", selection_tag, ".csv")), row.names = FALSE)
  write.csv(make_dot_table(selected), file.path(CONFIG$output_dir, paste0("GraphPad_dotplot_", selection_tag, ".csv")), row.names = FALSE)
  write.csv(make_forest_table(sig), file.path(CONFIG$output_dir, "GraphPad_forest_ALL_significant_TFs.csv"), row.names = FALSE)
  write.csv(make_bar_table(sig), file.path(CONFIG$output_dir, "GraphPad_bar_ALL_significant_TFs.csv"), row.names = FALSE)
  write.csv(make_dot_table(sig), file.path(CONFIG$output_dir, "GraphPad_dotplot_ALL_significant_TFs.csv"), row.names = FALSE)


  readme <- c(
    "dnaA RNA TF enrichment publication figures V1",
    "",
    paste0("Input: ", normalizePath(CONFIG$input_csv, winslash = "/", mustWork = TRUE)),
    paste0("Plot-eligible significant TFs: ", nrow(sig), " at BH-FDR < ", CONFIG$fdr_threshold, "."),
    paste0("Main display: Top ", top_n, " by FDR plus prespecified ", CONFIG$prespecified_tf, "."),
    paste0("Displayed TFs: ", paste(selected_terms, collapse = ", "), "."),
    "Forest: OR on log axis with 95% CI; fixed point size; colour=-log10(FDR); Hits as text.",
    "Bar: length=log2FE; colour=-log10(FDR); Hits as text.",
    "Dot: x=log2FE; marker size represents Hits; colour=-log10(FDR).",
    "FE = gene_ratio/background_ratio.",
    paste0("Continuous colour limits: ", signif(colour_limits[1], 7), " to ", signif(colour_limits[2], 7), "."),
    paste0("Colour palette: ", paste(CONFIG$evidence_palette, collapse = " -> "), "."),
    "ALL_significant_TFs_plot_parameters_and_RGB.csv contains all plot-eligible significant TFs.",
    "GraphPad_*_ALL_significant_TFs.csv contains all plot-eligible significant TFs.",
    "ALL_statistically_significant_TFs.csv preserves every significant TF, including infinite OR/CI.",
    "TF_plot_data.csv records exclusions from plotting; statistical FDR is unchanged.",
    "GraphPad_*_Top*.csv contains the main selected set.",
    "Each plot PDF contains the complete vertical continuous FDR colourbar.",
    "04_shared_continuous_FDR_colorbar.pdf is also provided as a standalone vertical legend.",
    "The dot plot retains a separate Hit genes size legend.",
    paste0("Physical panel: ", CONFIG$panel_width_cm, " x ", CONFIG$panel_height_cm, " cm."),
    paste0("Plain frame: ", CONFIG$frame_line_pt, " pt; font: ", CONFIG$font_family, " ", CONFIG$font_size_pt, " pt."),
    "PDF page is larger than the panel to preserve labels and legends.",
    "TF target enrichment does not establish TF activity or direct regulation."
  )
  writeLines(readme, file.path(CONFIG$output_dir, "README.txt"), useBytes = TRUE)

  message("Output directory: ", normalizePath(CONFIG$output_dir, winslash = "/"))
  message("Plot-eligible significant TFs exported: ", nrow(sig))
  message("Main TFs: ", paste(selected_terms, collapse = ", "))
  message("PDFs: ", paste(basename(c(forest_pdf, bar_pdf, dot_pdf, colorbar_pdf)), collapse = "; "))

  status <- data.frame(status = "completed", n_tests = nrow(dat),
    n_statistically_significant = sum(stat_sig),
    n_plot_eligible_significant = nrow(sig), n_selected = nrow(selected),
    selected_TFs = paste(selected_terms, collapse = ";"))

  invisible(list(status = "completed", selected_terms = selected_terms,
    plot_parameters = selected_out, files = c(forest_pdf, bar_pdf, dot_pdf, colorbar_pdf)))

}

if (RUN_TF && isTRUE(TF_PUBLICATION_CONFIG$enabled)) {
  tf_publication_results <- run_tf_publication_figures(
    dat = tf_v62_results$results,
    source_csv = file.path(tf_v62_results$config$output_dir,
                          "TF__Response__absolute.csv"),
    output_dir = file.path(tf_v62_results$config$output_dir,
                           "plots"),
    settings = TF_PUBLICATION_CONFIG
  )
}
message("COG, TF enrichment and publication outputs completed.")

}

# 3. Execute the analysis and export step-specific results
run_dnaa_analysis<-function() {
  need(c("tidyverse","data.table","cluster","cowplot","patchwork","ComplexHeatmap","circlize","ggrepel"))
  suppressPackageStartupMessages({library(tidyverse);library(data.table);library(cluster);library(grid)})
  out<-dir_make(file.path(cfg$results_dir,"dnaA_upshift"))
  prepare_dnaa(cfg)
  oldwd<-getwd();on.exit(setwd(oldwd),add=TRUE);setwd(out)
  # Suppress automatic Rplots.pdf files opened for plot-layout measurement.
  old_device<-getOption("device")
  options(device=function(...)grDevices::pdf(file=NULL))
  on.exit(options(device=old_device),add=TRUE)
  dir_make("Step2_RNA_filtering")
  if(cfg$run_original_extended_plots) DNA_copy_number_heatmaps()
  RNA_relative_abundance()
  RNA_expression_filters()
  if(cfg$run_original_extended_plots) RNA_heatmaps()
  DNA_dosage_correction()
  if(cfg$run_original_extended_plots) Corrected_RNA_heatmaps()
  Response_identification()
  Kmeans_clustering_and_diagnostics()
  Cluster_COG_TF_enrichment()
  Responsive_set_COG_TF_enrichment()
  message("dnaA upshift analysis completed. Output: ",out)
  invisible(out)
}
run_dnaa_analysis()
