# PPM reproduction modules from the original dnaA V5 and hns V1 scripts.
# The analytical function bodies are identical in both scripts.
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Expected EXPERIMENT, LFQ_CSV, PPM_CSV, METADATA_CSV, OUTPUT_DIR.")
}
experiment <- args[[1]]
if (!experiment %in% c("dnaA", "hns")) stop("Unknown experiment.")
is_hns <- experiment == "hns"
library_path <- Sys.getenv("BATCH2_R_LIBRARY")
if (nzchar(library_path)) .libPaths(c(library_path, .libPaths()))
CONFIG <- list(
  analysis_prefix = experiment,
  lfq_file = args[[2]],
  ppm_file = args[[3]],
  metadata_file = args[[4]],
  output_dir = args[[5]],
  protein_group_col = "Protein.Group",
  matrix_gene_col = "gene",
  sample_col = "Sample",
  time_col = "Time_h",
  phase_col = "Phase",
  trajectory_col = "Biological_replicates",
  batch_col = "Proteinomics_Batch",
  group_col = "Group",
  selected_batch = "2",
  selected_group = paste0(experiment, "_Upshift"),
  selected_trajectories = paste0(experiment, "_Upshift_replicate", c("A", "B", "C")),
  baseline_values = "Pre",
  time_round_digits = 4,
  drop_metadata_samples_missing_from_matrices = is_hns,
  minimum_model_timepoint_n = if (is_hns) 2L else 3L,
  require_exact_timepoint_n = !is_hns,
  multi_gene_delimiter_regex = "[;,/|]",
  ppm_low_filter = TRUE,
  low_min = 1,
  low_max = 10,
  low_mean = 10,
  baseline_stability_hard_filter = TRUE,
  baseline_ratio_quantile = .90,
  baseline_ratio_floor = 1.15,
  fdr_threshold = .05,
  fc_threshold_fold = 1.5,
  fc_threshold_log2 = log2(1.5),
  minimum_peak_n = if (is_hns) 2L else 3L
)

dir_make <- function(...) dir.create(..., recursive = TRUE, showWarnings = FALSE)
path_out <- function(...) file.path(CONFIG$output_dir, ...)
module_dir <- function(x) { d <- path_out(x); dir_make(d); d }
read_csv_safe <- function(x) read.csv(x,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8-BOM")
read_annotation <- function(x) read_csv_safe(x)
write_csv_safe <- function(x, y) write.csv(x, y, row.names = FALSE, na = "")
clean_string <- function(x) trimws(as.character(x))
normal_key <- function(x) tolower(clean_string(x))
is_baseline <- function(x) clean_string(x) %in% CONFIG$baseline_values
mean_or_na <- function(x) { x <- x[is.finite(x)]; if(length(x)) mean(x) else NA_real_ }
track_dir <- function(track) module_dir(if(track == "LFQ") "05_LFQ" else "06_PPM")
track_label <- function(track) if(track == "LFQ") "log2 MaxLFQ" else "log2 PPM"
track_prefix <- function(track) paste0(CONFIG$analysis_prefix, "_", track, "_track_")
track_module_dir <- function(track, code, label) {
  d <- file.path(track_dir(track), paste0(code, "_", label)); dir_make(d); d
}
track_path <- function(track, code, label, filename) {
  file.path(track_module_dir(track, code, label), paste0(track_prefix(track), filename))
}

shared_path <- function(code, label, filename) {
  d<-module_dir(paste0(code,"_",label));file.path(d,paste0(CONFIG$analysis_prefix,"_track_",filename))
}

assert_config <- function() {
  if (!requireNamespace("limma", quietly = TRUE)) stop("Package limma is required.")
  if (CONFIG$fc_threshold_fold <= 1) stop("fc_threshold_fold must be > 1.")
  if (!is.numeric(CONFIG$minimum_peak_n) || length(CONFIG$minimum_peak_n) != 1L ||
      CONFIG$minimum_peak_n < 2L || CONFIG$minimum_peak_n > length(CONFIG$selected_trajectories)) {
    stop("minimum_peak_n must be between 2 and the number of selected biological trajectories.")
  }
  if (CONFIG$minimum_model_timepoint_n < 2L || CONFIG$minimum_model_timepoint_n > 3L) stop("minimum_model_timepoint_n must be 2 or 3.")
  if (CONFIG$baseline_ratio_quantile <= 0 || CONFIG$baseline_ratio_quantile >= 1) stop("Invalid baseline ratio quantile.")
}

matrix_header <- function(path) names(read.csv(path, check.names=FALSE, stringsAsFactors=FALSE, nrows=0, fileEncoding="UTF-8-BOM"))

load_metadata <- function() {
  m <- read_csv_safe(CONFIG$metadata_file)
  req <- c(CONFIG$sample_col, CONFIG$time_col, CONFIG$phase_col, CONFIG$trajectory_col,
           CONFIG$batch_col, CONFIG$group_col)
  if (!all(req %in% names(m))) stop("Metadata lacks: ", paste(setdiff(req, names(m)), collapse=", "))
  keep <- clean_string(m[[CONFIG$batch_col]]) == CONFIG$selected_batch &
    clean_string(m[[CONFIG$group_col]]) == CONFIG$selected_group &
    clean_string(m[[CONFIG$trajectory_col]]) %in% CONFIG$selected_trajectories
  m <- m[keep,,drop=FALSE]; m$.sample <- clean_string(m[[CONFIG$sample_col]])
  m$.time <- round(suppressWarnings(as.numeric(m[[CONFIG$time_col]])), CONFIG$time_round_digits)
  m$.phase <- clean_string(m[[CONFIG$phase_col]]); m$.trajectory <- clean_string(m[[CONFIG$trajectory_col]])
  if (!nrow(m) || any(!is.finite(m$.time)) || anyDuplicated(m$.sample)) stop("Invalid selected metadata subset.")
  lfq_samples <- matrix_header(CONFIG$lfq_file); ppm_samples <- matrix_header(CONFIG$ppm_file)
  availability <- data.frame(sample=m$.sample, trajectory=m$.trajectory, time_h=m$.time,
    phase=m$.phase, in_LFQ=m$.sample%in%lfq_samples, in_PPM=m$.sample%in%ppm_samples,
    retained_in_model=(m$.sample%in%lfq_samples)&(m$.sample%in%ppm_samples), stringsAsFactors=FALSE)
  if (isTRUE(CONFIG$drop_metadata_samples_missing_from_matrices)) {
    m <- m[availability$retained_in_model,,drop=FALSE]
  } else if (any(!availability$retained_in_model)) {
    stop("Selected metadata contain samples absent from LFQ/PPM matrices: ",
         paste(availability$sample[!availability$retained_in_model], collapse=", "))
  }
  m <- m[order(match(m$.trajectory, CONFIG$selected_trajectories), m$.time, m$.sample),,drop=FALSE]
  attr(m,"sample_availability_audit") <- availability
  m
}

matrix_file <- function(track) if(track == "LFQ") CONFIG$lfq_file else CONFIG$ppm_file
load_mapped_matrix <- function(track, eligible = TRUE) {
  p <- if(!eligible) shared_path("01","mapping",paste0(tolower(track),"_one_to_one.csv")) else
    shared_path("04","baseline_stability",paste0(tolower(track),"_eligible.csv"))
  if(!file.exists(p)) stop("Missing upstream matrix: ", p)
  read_csv_safe(p)
}
numeric_matrix_from <- function(d, m) {
  miss <- setdiff(m$.sample, names(d)); if(length(miss)) stop("Matrix lacks samples: ", paste(miss, collapse=", "))
  x <- as.matrix(sapply(d[,m$.sample,drop=FALSE], function(z) suppressWarnings(as.numeric(z))))
  if(is.null(dim(x))) x <- matrix(x, ncol=1)
  rownames(x) <- d$gene; colnames(x) <- m$.sample; x
}

# =============================================================================
# 2. mapping and shared design audit
# =============================================================================
parse_mapping <- function(d, source) {
  pg <- clean_string(d[[CONFIG$protein_group_col]]); raw <- clean_string(d[[CONFIG$matrix_gene_col]])
  parts <- lapply(raw, function(x) { if(is.na(x) || !nzchar(x)) return(character()); unique(clean_string(unlist(strsplit(x, CONFIG$multi_gene_delimiter_regex, perl=TRUE)))) })
  parts <- lapply(parts, function(z) z[!is.na(z) & nzchar(z)])
  data.frame(source=source, row_id=seq_len(nrow(d)), Protein.Group=pg, matrix_gene=raw,
             parsed_gene=vapply(parts, function(z) if(length(z)==1) z else NA_character_, ""),
             n_parsed_genes=lengths(parts), stringsAsFactors=FALSE)
}

run_mapping <- function() {
  lfq <- read_csv_safe(CONFIG$lfq_file); ppm <- read_csv_safe(CONFIG$ppm_file)
  req <- c(CONFIG$protein_group_col, CONFIG$matrix_gene_col)
  if(!all(req %in% names(lfq)) || !all(req %in% names(ppm))) stop("Both matrices require Protein.Group and gene.")
  a <- rbind(parse_mapping(lfq,"LFQ"), parse_mapping(ppm,"PPM"))
  a$status <- ifelse(!nzchar(a$Protein.Group), "exclude_missing_protein_group",
    ifelse(a$n_parsed_genes==0, "exclude_missing_gene", ifelse(a$n_parsed_genes>1, "exclude_multi_gene_protein_group", "candidate")))
  cand <- a$status=="candidate"; pg_n <- table(a$Protein.Group[cand]); g_n <- table(a$parsed_gene[cand])
  a$protein_group_candidate_rows <- as.integer(pg_n[a$Protein.Group]); a$protein_group_candidate_rows[is.na(a$protein_group_candidate_rows)] <- 0L
  a$gene_candidate_rows <- as.integer(g_n[a$parsed_gene]); a$gene_candidate_rows[is.na(a$gene_candidate_rows)] <- 0L
  a$status[cand & a$protein_group_candidate_rows != 2L] <- "exclude_absent_or_duplicate_between_tracks"
  a$status[cand & a$gene_candidate_rows != 2L] <- "exclude_gene_maps_multiple_protein_groups"
  # Require an identical unique PG-gene pair in both tracks.
  pair <- paste(a$Protein.Group, a$parsed_gene, sep="\r")
  pair_n <- table(pair[a$status=="candidate"])
  a$status[a$status=="candidate" & as.integer(pair_n[pair]) != 2L] <- "exclude_inconsistent_between_tracks"
  a$retain <- a$status=="candidate"
  write_csv_safe(a,shared_path("01","mapping","protein_group_gene_mapping_audit.csv"))
  write_csv_safe(as.data.frame(table(a$source,a$status),stringsAsFactors=FALSE),shared_path("01","mapping","mapping_funnel_by_track.csv"))
  keep_pairs <- unique(pair[a$retain]);
  for(track in c("LFQ","PPM")) {
    raw <- if(track=="LFQ") lfq else ppm; z <- a[a$source==track,,drop=FALSE]
    # `pair` is indexed on the combined LFQ+PPM audit.  Rebuild it after
    # subsetting; using `pair[z$retain]` here would apply a shorter logical
    # index to the combined vector and can misalign Protein.Group and gene.
    pair_track <- paste(z$Protein.Group, z$parsed_gene, sep="\r")
    keep <- z$retain & pair_track %in% keep_pairs
    out <- raw[keep,,drop=FALSE]; out$gene <- z$parsed_gene[keep]
    if(anyDuplicated(out$gene) || anyDuplicated(out[[CONFIG$protein_group_col]])) stop("Strict mapping invariant failed.")
    write_csv_safe(out,shared_path("01","mapping",paste0(tolower(track),"_one_to_one.csv")))
  }
}

run_validate <- function() {
  m <- load_metadata()
  tab <- do.call(rbind,lapply(CONFIG$selected_trajectories,function(tr) {z<-m[m$.trajectory==tr,,drop=FALSE]; data.frame(trajectory=tr,n_samples=nrow(z),n_pre=sum(is_baseline(z$.phase)),n_post=sum(!is_baseline(z$.phase)),pre_samples=paste(z$.sample[is_baseline(z$.phase)],collapse=";"),stringsAsFactors=FALSE)}))
  if(any(tab$n_pre!=2L) || any(tab$n_post<1L)) stop("Each A/B/C trajectory needs exactly two retained pre samples and post samples.")
  time_n <- aggregate(.sample~.time, m, length)
  if(any(time_n$.sample < CONFIG$minimum_model_timepoint_n)) stop("A model time point has fewer samples than minimum_model_timepoint_n.")
  if(isTRUE(CONFIG$require_exact_timepoint_n) && any(time_n$.sample != length(CONFIG$selected_trajectories))) stop("Aligned design requires every selected trajectory at every time point.")
  availability <- attr(m,"sample_availability_audit")
  write_csv_safe(m,shared_path("02","validation","metadata_validated_Batch2_ABC.csv"));write_csv_safe(tab,shared_path("02","validation","trajectory_design_audit.csv"));write_csv_safe(time_n,shared_path("02","validation","timepoint_n_audit.csv"))
  if(!is.null(availability)) write_csv_safe(availability,shared_path("02","validation","sample_availability_audit.csv"))
}

# =============================================================================
# 3. Shared PPM eligibility audit (pre-specified low abundance and stability)
# =============================================================================
run_low_filter <- function() {
  m<-load_metadata(); ppm<-load_mapped_matrix("PPM",FALSE); x<-numeric_matrix_from(ppm,m); tr<-CONFIG$selected_trajectories
  fail<-matrix(FALSE,nrow(x),length(tr),dimnames=list(rownames(x),tr)); stat<-list()
  for(r in tr){ z<-x[,m$.sample[m$.trajectory==r],drop=FALSE]; mn<-apply(z,1,function(v) if(all(!is.finite(v))) NA_real_ else min(v,na.rm=TRUE)); mx<-apply(z,1,function(v) if(all(!is.finite(v))) NA_real_ else max(v,na.rm=TRUE)); av<-apply(z,1,mean_or_na); f<-!is.finite(mn)|!is.finite(mx)|!is.finite(av)|((mn<CONFIG$low_min)&(mx<CONFIG$low_max))|(av<CONFIG$low_mean); fail[,r]<-f; stat[[r]]<-data.frame(min=mn,max=mx,mean=av,fail=f) }
  keep<-if(CONFIG$ppm_low_filter) rowSums(fail)==0 else rep(TRUE,nrow(x)); o<-data.frame(gene=rownames(x),keep_low=keep, stringsAsFactors=FALSE)
  for(r in tr){k<-make.names(r);o[[paste0("min_PPM_",k)]]<-stat[[r]]$min;o[[paste0("max_PPM_",k)]]<-stat[[r]]$max;o[[paste0("mean_PPM_",k)]]<-stat[[r]]$mean;o[[paste0("failed_low_",k)]]<-stat[[r]]$fail}
  write_csv_safe(o,shared_path("03","low_filter","ppm_low_filter_status.csv"));write_csv_safe(ppm[keep,,drop=FALSE],shared_path("03","low_filter","ppm_after_low_filter.csv"));write_csv_safe(data.frame(stage=c("one_to_one_input","retained","removed"),n=c(nrow(ppm),sum(keep),sum(!keep))),shared_path("03","low_filter","low_filter_funnel.csv"))
}

run_baseline_stability <- function() {
  m<-load_metadata();ppm<-read_csv_safe(shared_path("03","low_filter","ppm_after_low_filter.csv"));x<-numeric_matrix_from(ppm,m);tr<-CONFIG$selected_trajectories;ratio<-matrix(NA_real_,nrow(x),length(tr),dimnames=list(rownames(x),tr));theta<-setNames(rep(NA_real_,length(tr)),tr)
  for(r in tr){allc<-m$.sample[m$.trajectory==r];prec<-m$.sample[m$.trajectory==r & is_baseline(m$.phase)];precv<-apply(x[,prec,drop=FALSE],1,function(v)sd(v,na.rm=TRUE)/max(mean(v,na.rm=TRUE),1e-6));allcv<-apply(x[,allc,drop=FALSE],1,function(v)sd(v,na.rm=TRUE)/max(mean(v,na.rm=TRUE),1e-6));ratio[,r]<-precv/(allcv+1e-8);theta[r]<-max(stats::quantile(ratio[,r],CONFIG$baseline_ratio_quantile,na.rm=TRUE),CONFIG$baseline_ratio_floor)}
  fail<-sweep(ratio,2,theta,">");fail[!is.finite(fail)]<-TRUE;keep<-if(CONFIG$baseline_stability_hard_filter)rowSums(fail)==0 else rep(TRUE,nrow(x));o<-data.frame(gene=rownames(x),keep_baseline=keep,stringsAsFactors=FALSE)
  for(r in tr){k<-make.names(r);o[[paste0("preCV_over_allCV_",k)]]<-ratio[,r];o[[paste0("theta_",k)]]<-theta[r];o[[paste0("failed_baseline_",k)]]<-fail[,r]}
  write_csv_safe(o,shared_path("04","baseline_stability","baseline_stability_status.csv"));write_csv_safe(data.frame(trajectory=tr,observed_q90=vapply(tr,function(r)stats::quantile(ratio[,r],CONFIG$baseline_ratio_quantile,na.rm=TRUE),0.0),floor=CONFIG$baseline_ratio_floor,theta=theta),shared_path("04","baseline_stability","baseline_thresholds.csv"));write_csv_safe(data.frame(stage=c("after_low_filter","retained","removed"),n=c(nrow(ppm),sum(keep),sum(!keep))),shared_path("04","baseline_stability","baseline_stability_funnel.csv"))
  # PPM receives its PPM-specific eligibility subset.  LFQ retains the strict
  # mapping subset and later has its own non-imputed complete-case audit.
  write_csv_safe(ppm[keep,,drop=FALSE],shared_path("04","baseline_stability","ppm_eligible.csv"))
  write_csv_safe(load_mapped_matrix("LFQ",FALSE),shared_path("04","baseline_stability","lfq_eligible.csv"))
}

build_response <- function(track) {
  m<-load_metadata(); d0<-load_mapped_matrix(track,TRUE); x0<-numeric_matrix_from(d0,m)
  # Direct log2 PPM only for positive PPM; MaxLFQ input is already log2 transformed.
  x<-if(track=="PPM") {z<-x0;z[!is.finite(z)|z<=0]<-NA_real_;log2(z)} else x0
  cc<-complete.cases(x); x<-x[cc,,drop=FALSE]; ann<-d0[match(rownames(x),d0$gene),,drop=FALSE]
  time_num<-m$.time; time_label<-paste0("t_",gsub("[.-]","_",time_num)); time<-factor(time_label,levels=unique(time_label[order(time_num)]));design<-stats::model.matrix(~0+time);colnames(design)<-levels(time)
  pre<-unique(time_label[time_num<0]);post<-unique(time_label[time_num>=0]);if(length(pre)!=2L||!length(post))stop("Need exactly two pre and at least one post time labels.")
  ctext<-vapply(post,function(p)paste0(p,"-(",paste(pre,collapse="+"),")/",length(pre)),"");names(ctext)<-paste0("post_",post,"_vs_pre_mean");contrast<-limma::makeContrasts(contrasts=ctext,levels=design)
  block<-factor(sub(".*replicate","",m$.trajectory));corfit<-limma::duplicateCorrelation(x,design,block=block);fit<-limma::lmFit(x,design,block=block,correlation=corfit$consensus.correlation);fit2<-limma::eBayes(limma::contrasts.fit(fit,contrast),trend=TRUE,robust=TRUE)
  p<-fit2$F.p.value;q<-stats::p.adjust(p,"BH"); co<-fit2$coefficients
  # makeContrasts() may replace the supplied contrast names with their algebraic
  # expressions (for example `post_t_0_0333-(...)/2`).  Do not parse time from
  # those expressions: the coefficient columns are in the same fixed order as
  # `post`, so restore the explicit post-time labels before peak selection.
  if(ncol(co) != length(post)) stop("Contrast/post-time dimension mismatch.")
  colnames(co) <- post
  # This is deliberately post-only: baseline samples cannot become a response peak.
  peak_i<-apply(abs(co),1,which.max);peak_effect<-co[cbind(seq_len(nrow(co)),peak_i)];peak_time<-suppressWarnings(as.numeric(gsub("_",".",sub("^t_","",colnames(co)[peak_i]))))
  rep_effect<-matrix(NA_real_,nrow(x),length(post),dimnames=list(rownames(x),post)); # expanded below for selected peak
  n_available<-rep(NA_integer_,nrow(x)); concordant<-rep(FALSE,nrow(x)); rep_text<-rep("",nrow(x))
  for(i in seq_len(nrow(x))){t<-peak_time[i];v<-vapply(CONFIG$selected_trajectories,function(tr){post_s<-m$.sample[m$.trajectory==tr & m$.time==t];pre_s<-m$.sample[m$.trajectory==tr & m$.time<0];if(length(post_s)==1L&&length(pre_s)==2L) x[i,post_s]-mean(x[i,pre_s]) else NA_real_},0.0);good<-is.finite(v);n_available[i]<-sum(good);concordant[i]<-sum(good)>=CONFIG$minimum_peak_n && (all(v[good]>0)||all(v[good]<0));rep_text[i]<-paste(sprintf("%s=%.6f",LETTERS[seq_along(v)],v),collapse=";")}
  out<-ann;out$tested<-TRUE;out$Omnibus_F<-fit2$F;out$Omnibus_P_Value<-p;out$BH_q<-q;out$peak_time_h<-peak_time;out$peak_log2_effect_signed<-peak_effect;out$peak_abs_log2_effect<-abs(peak_effect);out$peak_n_available<-n_available;out$peak_replicate_effects<-rep_text;out$peak_direction_concordant<-concordant;out$Direction<-ifelse(peak_effect>0,"Up",ifelse(peak_effect<0,"Down",NA));out$Significant<-is.finite(q)&q<CONFIG$fdr_threshold&abs(peak_effect)>=CONFIG$fc_threshold_log2&concordant
  for(j in seq_len(ncol(co)))out[[paste0("effect_post_time_",colnames(co)[j])]]<-co[,j]
  write_csv_safe(data.frame(Protein.Group=d0[[CONFIG$protein_group_col]],gene=d0$gene,tested=cc,missing_samples=apply(is.na(x0),1,function(z)paste(m$.sample[z],collapse=";"))),track_path(track,"04","response","complete_case_audit.csv"))
  write_csv_safe(out,track_path(track,"04","response","formal_response_all.csv"))
  write_csv_safe(out[out$Significant,,drop=FALSE],track_path(track,"04","response","response_molecules.csv"))
  write_csv_safe(design,track_path(track,"04","response","design_matrix.csv"))
  write_csv_safe(as.data.frame(contrast),track_path(track,"04","response","contrast_matrix.csv"))
  write_csv_safe(m,track_path(track,"04","response","model_sample_manifest.csv"))
  write_csv_safe(data.frame(stage=c("eligible_input","complete_case","valid_omnibus_p","BH_family","responders"),n=c(nrow(d0),nrow(x),sum(is.finite(p)),sum(is.finite(p)),sum(out$Significant))),track_path(track,"04","response","response_funnel.csv"))
  tp_stats <- do.call(rbind,lapply(seq_along(post),function(j){tt<-suppressWarnings(as.numeric(gsub("_",".",sub("^t_","",post[j]))));data.frame(gene=rownames(co),Protein.Group=ann[[CONFIG$protein_group_col]],time_h=tt,n_available=sum(m$.time==tt),estimated_log2_effect=co[,j],omnibus_p=p,BH_q=q,stringsAsFactors=FALSE)}))
  write_csv_safe(tp_stats,track_path(track,"04","response","timepoint_response_statistics.csv"))
  time_counts<-table(m$.time);diagnostics<-data.frame(track=track,scale=track_label(track),post_only_peak=TRUE,design_rank=qr(design)$rank,design_columns=ncol(design),contrasts=ncol(contrast),consensus_correlation=unname(corfit$consensus.correlation),residual_df_min=min(fit2$df.residual),model_timepoint_n_min=min(time_counts),model_timepoint_n_max=max(time_counts),status="PASS",stringsAsFactors=FALSE)
  write_csv_safe(diagnostics,track_path(track,"04","response","model_diagnostics.csv"))
  if(requireNamespace("jsonlite",quietly=TRUE))jsonlite::write_json(as.list(diagnostics[1,]),track_path(track,"04","response","model_diagnostics.json"),auto_unbox=TRUE,pretty=TRUE)
  # Formal inference remains on log2 scale. For RNA-style clustering/display,
  # use linear relative abundance divided by the arithmetic pre-sample mean.
  # This is not 2^(log2 value - mean(log2 pre)), which uses a geometric pre mean.
  long<-do.call(rbind,lapply(seq_len(nrow(m)),function(j){tr<-m$.trajectory[j];precols<-m$.sample[m$.trajectory==tr&m$.time<0];base_log<-rowMeans(x[,precols,drop=FALSE],na.rm=TRUE);base_linear<-rowMeans(2^x[,precols,drop=FALSE],na.rm=TRUE);linear_value<-2^x[,m$.sample[j]];data.frame(gene=rownames(x),Protein.Group=ann[[CONFIG$protein_group_col]],sample=m$.sample[j],trajectory=tr,replicate=sub(".*replicate","",tr),time_h=m$.time[j],phase=m$.phase[j],log2_value=x[,m$.sample[j]],log2_relative_to_pre=x[,m$.sample[j]]-base_log,relative_abundance=linear_value/base_linear,stringsAsFactors=FALSE)}))
  write_csv_safe(long,track_path(track,"04","response","relative_abundance_long.csv"))
}

assert_config()
dir_make(CONFIG$output_dir)
run_mapping()
run_validate()
run_low_filter()
run_baseline_stability()
build_response("PPM")
writeLines(capture.output(sessionInfo()), path_out("R_sessionInfo.txt"))
