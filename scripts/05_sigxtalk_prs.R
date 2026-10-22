#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly=TRUE)
if(length(argv)!=2) stop('Usage: Rscript 05_sigxtalk_prs.R config.json smoke|formal')
this <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(this)), 'sigxtalk_common.R'))
cfg <- read_cfg(argv[1],argv[2]);load_runtime(cfg)

bounded_prs <- function(cfg) {
  # Same forest model, response and importance definition as the imported upstream function.
  # Worker/task seeds specify one stochastic realization; they do not promise old numeric identity.
  ns <- asNamespace('SigXTalkR');f <- SigXTalkR::PRS_calc
  original <- deparse(f,width.cutoff=500); patched <- original
  workers <- as.integer(cfg$prs$workers)
  stopifnot(workers>=1)
  at <- grep('min(24, Ncores)',patched,fixed=TRUE)
  if(length(at)!=1) stop('Pinned PRS worker expression changed; refuse an unreviewed patch')
  patched[at] <- gsub('min(24, Ncores)',paste0('min(',workers,', Ncores)'),patched[at],fixed=TRUE)
  patched <- append(patched,c(
    paste0('parallel::clusterSetRNGStream(cl, iseed=',as.integer(cfg$seed),')'),
    'parallel::clusterEvalQ(cl, options(ranger.num.threads=1))'),after=at)
  task <- grep('tftg_name <- TFTG_All[i]',patched,fixed=TRUE)
  if(length(task)!=1) stop('Cannot locate fixed PRS task seed insertion point')
  patched <- append(patched,paste0('set.seed(',as.integer(cfg$seed),'L + i)'),after=task-1L)
  writeLines(original,file.path(cfg$mode_dir,'audit/PRS_calc_upstream_function.R'))
  writeLines(patched,file.path(cfg$mode_dir,'audit/PRS_calc_runtime_resource_function.R'))
  fun <- eval(parse(text=paste(patched,collapse='\n')),envir=ns)
  environment(fun) <- ns
  write_json(list(upstream='SigXTalkR::PRS_calc',license='MIT;Jiawen Hou and official authors',
    changes=c('SOCK worker cap','explicit worker and per-task random seeds','ranger thread option'),
    unchanged=c('TF-times-target normalized-expression response','ranger regression forest',
      '500 trees','mtry=3','permutation importance','coexpression cutoff','importance>1e-5 collection filter'),
    workers=workers,seed=cfg$seed,
    mathematical_definition_unchanged=TRUE,old_numerical_values_expected_identical=FALSE),
    file.path(cfg$mode_dir,'audit/PRS_RESOURCE_PATCH_AUDIT.json'))
  fun
}

official_fid_spe <- function(data) {
  # Both matrix functions are imported unchanged. Denominators are the supplied universe.
  if(!nrow(data)) return(data.frame())
  stopifnot(!anyDuplicated(data[,1:3]),all(is.finite(data$Weight)),all(data$Weight>0))
  fid <- setNames(lapply(unique(data$Target),function(g)
    SigXTalkR::Calculate_Fidelity_Matrix(data,g)$all),unique(data$Target))
  spe <- setNames(lapply(unique(data$Receptor),function(g)
    SigXTalkR::Calculate_Specificity_Matrix(data,g)$all),unique(data$Receptor))
  target_denom <- tapply(data$Weight,data$Target,sum)
  rec_denom <- tapply(data$Weight,data$Receptor,sum)
  rows <- lapply(seq_len(nrow(data)),function(i) {
    z <- data[i,]
    data.frame(Receptor=z$Receptor,SSC=z$SSC,Target=z$Target,PRS=z$Weight,
      fidelity_target_all=fid[[z$Target]][z$Receptor,z$SSC],
      specificity_receptor_all=spe[[z$Receptor]][z$Target,z$SSC],
      fidelity_denominator_PRS=unname(target_denom[z$Target]),
      specificity_denominator_PRS=unname(rec_denom[z$Receptor]),
      metric_origin='official Calculate_Fidelity_Matrix$all / Calculate_Specificity_Matrix$all',
      clinical_or_causal_interpretation=FALSE,stringsAsFactors=FALSE)
  })
  out <- do.call(rbind,rows)
  # Arithmetic identity check of official results only; never substitute a new metric.
  stopifnot(all(is.finite(out$fidelity_target_all)),all(is.finite(out$specificity_receptor_all)),
    max(abs(out$fidelity_target_all-out$PRS/out$fidelity_denominator_PRS))<1e-10,
    max(abs(out$specificity_receptor_all-out$PRS/out$specificity_denominator_PRS))<1e-10)
  out
}

run_checked(cfg,'R_PRS',function() {
  wd <- cfg$mode_dir
  input <- jsonlite::fromJSON(file.path(wd,'audit/NEW_INPUT_PROVENANCE.json'))
  hgstatus <- jsonlite::fromJSON(file.path(wd,'audit/HGNN_STATUS.json'))
  if(hgstatus$status!='DONE') stop('New official HGNN stage has not completed')
  predfile <- file.path(wd,'hgnn_outputs',input$project_id,paste0('pathways_',input$target_type,'.csv'))
  pp <- read.csv(predfile,check.names=FALSE,stringsAsFactors=FALSE)
  stopifnot(all(c('Receptor','TF','TG','pred_label') %in% names(pp)),nrow(pp)>0,
            all(is.finite(pp$pred_label)),!anyDuplicated(pp[,c('Receptor','TF','TG')]))
  active <- pp[pp$pred_label>cfg$hgnn$activation_threshold,c('Receptor','TF','TG'),drop=FALSE]
  write_tsv(active,file.path(wd,'prs_outputs/ALL_NEW_PA_GT_THRESHOLD_AXES.tsv'))
  if(!nrow(active)) stop('No new pathways exceed the fixed PA threshold; PRS is not estimable')
  prs_axes <- active
  if(cfg$mode=='smoke' && cfg$smoke$max_prs_TFTG_pairs>0) {
    # Predeclared technical-scope cap, not a biological pathway or candidate selection.
    pair_keys <- paste(active$TF,active$TG,sep='::')
    selected_pairs <- head(sort(unique(pair_keys)),cfg$smoke$max_prs_TFTG_pairs)
    prs_axes <- active[pair_keys %in% selected_pairs,,drop=FALSE]
  }
  write_tsv(prs_axes,file.path(wd,'prs_outputs/NEW_PRS_INPUT_AXES.tsv'))
  exp <- readRDS(file.path(wd,'r_objects/receiver_logCP10K_new.rds'))
  stopifnot(all(unlist(prs_axes) %in% rownames(exp)))
  fun <- if(isTRUE(cfg$prs$bounded_resource_wrapper)) bounded_prs(cfg) else SigXTalkR::PRS_calc
  set.seed(as.integer(cfg$seed))
  raw <- fun(exp,prs_axes,cutoff=cfg$prs$coexpression_cutoff)
  if(!is.data.frame(raw) || !nrow(raw)) stop('Official PRS returned no rows; do not borrow historical PRS')
  stopifnot(identical(names(raw),c('Receptor','SSC','Target','Weight')),
            all(is.finite(raw$Weight)),all(raw$Weight>0),!anyDuplicated(raw[,1:3]))
  saveRDS(raw,file.path(wd,'prs_outputs/PRS_raw_new.rds'))
  write_tsv(raw,file.path(wd,'prs_outputs/PRS_RAW_OFFICIAL_NEW.tsv'))
  # Actual function formal is PRS_thres, not the tutorial markdown's PRS_threshold.
  filtered <- SigXTalkR::Filter_results(raw,PRS_thres=cfg$prs$relative_max_threshold,
                                      remove_genes=cfg$prs$remove_mito_ribosomal_targets)
  write_tsv(filtered,file.path(wd,'prs_outputs/PRS_FILTERED_OFFICIAL_NEW.tsv'))
  stopifnot(nrow(filtered)>0)
  fid_status <- 'DONE';fid_error <- NULL
  tryCatch({
    raw_metrics <- official_fid_spe(raw);raw_metrics$normalization_universe <- 'all_new_returned_positive_PRS_rows'
    filtered_metrics <- official_fid_spe(filtered)
    filtered_metrics$normalization_universe <- 'new_Filter_results_relative_max_filtered_rows'
    write_tsv(raw_metrics,file.path(wd,'prs_outputs/FID_SPE_RAW_UNIVERSE_OFFICIAL_NEW.tsv'))
    write_tsv(filtered_metrics,file.path(wd,'prs_outputs/FID_SPE_FILTERED_UNIVERSE_OFFICIAL_NEW.tsv'))
  },error=function(e) {fid_status <<- 'PARTIAL_WITHOUT_FID_SPEC';fid_error <<- conditionMessage(e)})
  # Attach original new values; unreturned PRS remains NA, including activated paths.
  rr <- raw;names(rr) <- c('Receptor','TF','TG','PRS')
  full <- merge(pp,rr,by=c('Receptor','TF','TG'),all.x=TRUE,sort=FALSE)
  stopifnot(nrow(full)==nrow(pp),sum(!is.na(full$PRS))==nrow(raw))
  full$PA_threshold_pass <- full$pred_label>cfg$hgnn$activation_threshold
  full$positive_official_PRS_returned <- !is.na(full$PRS)
  full$PRS_relative_max_filter_pass <- do.call(paste,c(full[,c('Receptor','TF','TG')],sep='::')) %in%
     do.call(paste,c(filtered[,c('Receptor','SSC','Target')],sep='::'))
  full$method_origin <- 'NEW_OFFICIAL_SIGXTALK_RAW_COUNTS_R_HGNN_PRS'
  full$source_project_id <- input$project_id
  full$source_run_started_from_legacy_predictions <- FALSE
  write_tsv(full,file.path(wd,'tables/ALL_TESTED_NEW_AXES_WITH_PRS.tsv'))
  write_tsv(full[full$PA_threshold_pass,,drop=FALSE],file.path(wd,'tables/ACTIVATED_NEW_AXES_WITH_PRS.tsv'))
  outputs <- list.files(file.path(wd,'prs_outputs'),full.names=TRUE)
  write_tsv(data.frame(file=basename(outputs),bytes=file.info(outputs)$size,
      sha256=vapply(outputs,sha256_file,character(1))),file.path(wd,'audit/PRS_OUTPUT_FILE_HASHES.tsv'))
  result <- list(complete_R_python_R=nrow(raw)>0,all_tested_HGNN_rows=nrow(pp),
    active_PA_rows=nrow(active),actual_PRS_input_axes=nrow(prs_axes),positive_PRS_rows=nrow(raw),
    optional_relative_max_filter_threshold=cfg$prs$relative_max_threshold,
    relative_filter_is_not_pvalue=TRUE,filtered_PRS_rows=nrow(filtered),
    activated_rows_without_returned_PRS=sum(full$PA_threshold_pass & is.na(full$PRS)),
    receiver_expression_dimension=dim(exp),coexpression_cutoff=cfg$prs$coexpression_cutoff,
    fidelity_specificity_status=fid_status,fidelity_specificity_error=fid_error,
    new_input_provenance_sha256=sha256_file(file.path(wd,'audit/NEW_INPUT_PROVENANCE.json')),
    new_HGNN_predictions_sha256=sha256_file(predfile),
    new_raw_PRS_sha256=sha256_file(file.path(wd,'prs_outputs/PRS_RAW_OFFICIAL_NEW.tsv')),
    no_old_pathways_or_PRS_used=TRUE)
  write_json(result,file.path(wd,'audit/NEW_R_PYTHON_R_COMPLETION.json'))
  result
})
