# Protein abundance: mean and sample SD
#
# Purpose: Summarize the protein groups supplied in Table S4.
# Input: Table S4 sheet dnaA&hns_Proteins. No other worksheets are required.
# Output: Four CSV files in results/proteomics: dnaA/hns x PPM/RelPPM.
# Each row is a sampling time in minutes. Each protein has adjacent mean and
# SD columns, suitable for plotting directly. No figures are generated.
# Run with Rscript or source the saved file in RStudio.
# Required package: readxl. Packages are not installed automatically.
#
# PPM is already normalized in Table S4 and is not renormalized here.
# Relative PPM is calculated within each biological replicate by dividing
# each protein by its mean negative-time PPM. 


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

need("readxl")

# 2. Read the published protein groups and sample design.
d<-as.data.frame(readxl::read_excel(table_s4,sheet="dnaA&hns_Proteins",
                                  .name_repair="minimal"),check.names=FALSE)
if(!all(c("Protein.Group","gene") %in% names(d)))
  stop("Missing Protein.Group or gene column in dnaA&hns_Proteins")
if(anyNA(d$Protein.Group) || any(!nzchar(as.character(d$Protein.Group))) ||
   anyDuplicated(d$Protein.Group)) stop("Protein.Group must contain unique, nonempty IDs")
rx<-"^(dnaA|hns)upshift_Rep([0-9]+)_(-?[0-9.]+) min$"
samples<-grep(rx,names(d),value=TRUE)
if(!length(samples)) stop("No recognized protein sample columns")
parts<-regmatches(samples,regexec(rx,samples))
design<-data.frame(sample=samples,experiment=vapply(parts,"[","",2),
  replicate=as.integer(vapply(parts,"[","",3)),
  time_min=as.numeric(vapply(parts,"[","",4)))
if(anyDuplicated(design[c("experiment","replicate","time_min")]))
  stop("Duplicate biological-sample keys")

# Prefer readable gene labels; distinguish ambiguous labels by protein-group ID.
labels<-trimws(as.character(d$gene))
labels[is.na(labels)]<-""
labels[!nzchar(labels)]<-as.character(d$Protein.Group[!nzchar(labels)])
duplicated_labels<-duplicated(labels)|duplicated(labels,fromLast=TRUE)
labels[duplicated_labels]<-paste0(labels[duplicated_labels],"__",
                                d$Protein.Group[duplicated_labels])
labels<-make.unique(labels,sep="__")

mean_ok<-function(x) {
  x<-x[is.finite(x)]
  if(length(x)) mean(x) else NA_real_
}
sd_ok<-function(x) {
  x<-x[is.finite(x)]
  if(length(x)>1L) stats::sd(x) else NA_real_
}

# 3. Normalize within replicates, then summarize at each sampling time.
out<-file.path(script_dir,"results","proteomics")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(ex in unique(design$experiment)) {
  ds<-design[design$experiment==ex,,drop=FALSE]
  ppm<-as.matrix(data.frame(lapply(d[ds$sample],function(x)
    suppressWarnings(as.numeric(x))),check.names=FALSE))
  colnames(ppm)<-ds$sample
  relative<-ppm;relative[,]<-NA_real_
  for(rep in unique(ds$replicate)) {
    pre<-which(ds$replicate==rep & ds$time_min<0)
    if(!length(pre)) stop("No pre-upshift samples in ",ex," replicate ",rep)
    baseline_values<-ppm[,pre,drop=FALSE]
    baseline<-apply(baseline_values,1,mean_ok)
    valid<-rowSums(is.finite(baseline_values)&baseline_values>0)==length(pre)
    valid<-valid & is.finite(baseline) & baseline>0
    columns<-which(ds$replicate==rep)
    relative[valid,columns]<-ppm[valid,columns,drop=FALSE]/baseline[valid]
  }
  for(metric in c("PPM","RelPPM")) {
    values<-if(metric=="PPM") ppm else relative
    times<-sort(unique(ds$time_min))
    means<-sapply(times,function(t)apply(values[,ds$time_min==t,drop=FALSE],1,mean_ok))
    sds<-sapply(times,function(t)apply(values[,ds$time_min==t,drop=FALSE],1,sd_ok))
    means<-matrix(means,nrow=nrow(d),ncol=length(times))
    sds<-matrix(sds,nrow=nrow(d),ncol=length(times))
    wide<-data.frame(time_min=times)
    for(i in seq_len(nrow(d))) {
      wide[[paste0(labels[i],"_mean")]]<-means[i,]
      wide[[paste0(labels[i],"_SD")]]<-sds[i,]
    }
    utils::write.csv(wide,file.path(out,paste0(ex,"_",metric,"_mean_SD.csv")),
                     row.names=FALSE,na="",fileEncoding="UTF-8")
  }
}
message("Protein mean/SD tables completed for ",nrow(d)," published protein groups. Output: ",out)
