# Project orchestration helpers; upstream algorithms are imported, not vendored.
read_cfg <- function(path, mode) {
  stopifnot(mode %in% c('smoke', 'formal'))
  if (!requireNamespace('jsonlite', quietly=TRUE)) stop('Install jsonlite in the isolated R library')
  cfg <- jsonlite::fromJSON(path, simplifyVector=TRUE)
  cfg$config_file <- normalizePath(path, mustWork=TRUE)
  if (!grepl('^(/|[A-Za-z]:[/\\\\])', cfg$project_root)) cfg$project_root <- file.path(dirname(cfg$config_file),cfg$project_root)
  cfg$project_root <- normalizePath(cfg$project_root, mustWork=TRUE)
  resolve <- function(p) {
    if (is.null(p) || !nzchar(p)) return(NULL)
    if (!grepl('^(/|[A-Za-z]:[/\\\\])', p)) p <- file.path(cfg$project_root, p)
    normalizePath(p, mustWork=FALSE)
  }
  for (n in c('module_dir','input_dir','official_source','pairing_tsv','state_thresholds_tsv','pbmc_rds')) {
    cfg[[n]] <- resolve(cfg[[n]])
  }
  if (!is.null(cfg$r_library_paths)) {
    libs <- vapply(cfg$r_library_paths, resolve, character(1))
    .libPaths(c(libs, .libPaths()))
  }
  cfg$mode <- mode
  smoke_dir <- cfg$smoke$directory_name
  if (is.null(smoke_dir)) smoke_dir <- 'technical_smoke'
  if (length(smoke_dir)!=1L || !grepl('^[A-Za-z0-9][A-Za-z0-9_-]*$', smoke_dir) || smoke_dir=='formal_new') {
    stop('smoke.directory_name must be one safe directory name distinct from formal_new')
  }
  cfg$mode_dir <- file.path(cfg$module_dir, if(mode=='formal') 'formal_new' else smoke_dir)
  for (d in c('input','inputs','r_objects','hgnn_outputs','prs_outputs','tables','logs','audit')) {
    dir.create(file.path(cfg$mode_dir, d), recursive=TRUE, showWarnings=FALSE)
  }
  cfg
}
utc_now <- function() format(Sys.time(), tz='UTC', usetz=TRUE)
write_json <- function(x, p) jsonlite::write_json(x, p, auto_unbox=TRUE, pretty=TRUE, na='null', null='null')
write_tsv <- function(x, p) write.table(x, p, sep='\t', quote=FALSE, row.names=FALSE, na='NA')
sha256_file <- function(p) digest::digest(file=p, algo='sha256')
get_data <- function(obj, layer='data') SeuratObject::GetAssayData(obj, assay='RNA', layer=layer)
source_function <- function(file, name, ns) {
  ex <- parse(file)
  found <- Filter(function(e) is.call(e) && as.character(e[[1]]) %in% c('<-','=') &&
                  is.symbol(e[[2]]) && as.character(e[[2]]) == name, as.list(ex))
  if (length(found)!=1) stop('Cannot identify upstream function: ', name)
  eval(found[[1]][[3]], envir=ns)
}
audit_api <- function(cfg) {
  ns <- asNamespace('SigXTalkR')
  mapping <- list('MLNconstruct.R'=c('Infer_CCI','Prepare_Input','PRS_calc','Run_Regression'),
                  'Utilities.R'=c('Get_Exp_Clu','Filter_results'),
                  'Crosstalk_analysis.R'=c('Calculate_Fidelity_Matrix','Calculate_Specificity_Matrix'))
  rows <- list()
  for (f in names(mapping)) for (n in mapping[[f]]) {
    installed <- get(n, envir=ns)
    upstream <- source_function(file.path(cfg$official_source, 'R', f), n, ns)
    same <- identical(deparse(body(installed)), deparse(body(upstream))) &&
            identical(formals(installed), formals(upstream))
    rows[[length(rows)+1L]] <- data.frame(function_name=n,
      actual_formals=paste(names(formals(installed)), collapse=','),
      body_matches_pinned_checkout=same,
      installed_body_sha256=digest::digest(deparse(body(installed)),algo='sha256'),
      upstream_file_sha256=sha256_file(file.path(cfg$official_source,'R',f)))
  }
  d <- do.call(rbind, rows)
  write_tsv(d, file.path(cfg$mode_dir, 'audit/R_API_AND_BODY_AUDIT.tsv'))
  if (!all(d$body_matches_pinned_checkout)) stop('Installed algorithm differs from pinned checkout; inspect explicit compatibility patch before use')
  writeLines(capture.output(sessionInfo()), file.path(cfg$mode_dir,'audit/R_SESSION.txt'))
  d
}
load_runtime <- function(cfg) {
  options(Seurat.object.assay.version='v3', future.globals.maxSize=8*1024^3)
  for (p in c('digest','Matrix','Seurat','SeuratObject','SigXTalkR','CellChat','dplyr')) {
    if (!requireNamespace(p, quietly=TRUE)) stop('Missing isolated R dependency: ',p)
  }
  suppressPackageStartupMessages({library(Matrix); library(Seurat); library(SigXTalkR); library(CellChat); library(dplyr)})
  if (requireNamespace('future',quietly=TRUE)) future::plan('sequential')
  audit_api(cfg)
  set.seed(as.integer(cfg$seed))
}
run_checked <- function(cfg, stage, fun) {
  dest <- file.path(cfg$mode_dir,'audit',paste0(stage,'_STATUS.json'))
  start <- utc_now()
  tryCatch({
    result <- fun()
    write_json(c(list(status='DONE', stage=stage, mode=cfg$mode, started_utc=start,
                     finished_utc=utc_now(), new_computation=TRUE,
                     source_label=if(cfg$mode=='formal') 'NEW_OFFICIAL_IVDD_RUN' else
                       if(cfg$smoke$input_kind=='pbmc_rds') 'OFFICIAL_PBMC_REAL_DATA_SMOKE' else
                       'REAL_IVDD_API_TECHNICAL_SMOKE_NOT_PBMC_DEMO'),result),dest)
    invisible(result)
  }, error=function(e) {
    write_json(list(status='FAILED',stage=stage,mode=cfg$mode,started_utc=start,
                   failed_utc=utc_now(),error=conditionMessage(e),new_computation=TRUE),dest)
    stop(e)
  })
}
