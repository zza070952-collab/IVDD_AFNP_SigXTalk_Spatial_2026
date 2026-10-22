#!/usr/bin/env Rscript
argv <- commandArgs(trailingOnly=TRUE)
if (length(argv)!=2) stop('Usage: Rscript 03_prepare_sigxtalk.R config.json smoke|formal')
this <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(this)), 'sigxtalk_common.R'))
cfg <- read_cfg(argv[1],argv[2]); load_runtime(cfg)
run_checked(cfg,'R_PREPARE',function() {
  wd <- cfg$mode_dir
  use_pbmc <- cfg$mode=='smoke' && cfg$smoke$input_kind=='pbmc_rds'
  raw_hashes <- list()
  if (use_pbmc) {
    stopifnot(file.exists(cfg$pbmc_rds))
    original <- readRDS(cfg$pbmc_rds)
    x <- get_data(original,'counts')
    meta <- data.frame(cell_id=colnames(original), original_group=as.character(Seurat::Idents(original)),stringsAsFactors=FALSE)
    meta$bridge_group <- meta$original_group
    meta$donor_id <- 'official_PBMC_demo_donor_not_IVDD'
    rownames(meta) <- meta$cell_id
    receiver <- cfg$smoke$pbmc_receiver
    raw_hashes$pbmc_rds_sha256 <- sha256_file(cfg$pbmc_rds)
  } else {
    files <- file.path(cfg$input_dir,c('counts.mtx.gz','genes.tsv','meta.tsv'))
    if (!all(file.exists(files))) stop('Real selected raw counts / genes / metadata input is missing')
    x <- as(Matrix::readMM(gzfile(files[1])), 'CsparseMatrix')
    genes <- readLines(files[2],warn=FALSE)
    meta <- read.delim(files[3],check.names=FALSE,stringsAsFactors=FALSE)
    need <- c('cell_id','donor_id','gsm','geo_series','condition','bridge_group')
    if (!all(need %in% names(meta))) stop('Required original cell identity fields missing')
    stopifnot(nrow(x)==length(genes),ncol(x)==nrow(meta),!anyDuplicated(genes),
              !anyDuplicated(meta$cell_id),all(nzchar(genes)),all(nzchar(meta$cell_id)),
              all(is.finite(x@x)),all(x@x>=0),all(x@x==floor(x@x)),
              all(meta$geo_series=='GSE230808'),all(meta$condition=='diseased'))
    rownames(x) <- genes; colnames(x) <- meta$cell_id;rownames(meta) <- meta$cell_id
    stopifnot(all(c('AF_stromal','NP_GEP2_high','NP_GEP2_low') %in% meta$bridge_group))
    if (is.null(cfg$pairing_tsv) || !file.exists(cfg$pairing_tsv)) stop('Provide actual donor AF/NP GSM pairing table')
    pairs <- read.delim(cfg$pairing_tsv,stringsAsFactors=FALSE,check.names=FALSE)
    stopifnot(all(c('donor_id','AF_gsm','NP_gsm') %in% names(pairs)),!anyDuplicated(pairs$donor_id))
    pp <- pairs[match(meta$donor_id,pairs$donor_id),]
    expected <- ifelse(meta$bridge_group=='AF_stromal',pp$AF_gsm,pp$NP_gsm)
    if(anyNA(expected) || !all(expected==meta$gsm)) stop('Actual donor/tissue/GSM pair check failed')
    if (!is.null(cfg$state_thresholds_tsv) && file.exists(cfg$state_thresholds_tsv) && 'GEP2_top25_UCell' %in% names(meta)) {
      qs <- read.delim(cfg$state_thresholds_tsv,stringsAsFactors=FALSE)
      qs <- qs[match(meta$donor_id,qs$donor_id),]
      hi <- meta$bridge_group=='NP_GEP2_high';lo <- meta$bridge_group=='NP_GEP2_low'
      stopifnot(all(meta$GEP2_top25_UCell[hi]>=qs$q75[hi]),all(meta$GEP2_top25_UCell[lo]<=qs$q25[lo]))
    }
    raw_hashes <- as.list(vapply(c(files,cfg$pairing_tsv),sha256_file,character(1)))
    names(raw_hashes) <- c('raw_matrix_sha256','raw_genes_sha256','cell_identity_sha256','patient_pairing_sha256')
    receiver <- cfg$receiver
  }
  if (cfg$mode=='smoke') {
    # Predeclared technical input, real cells only; never interpreted as IVDD results.
    set.seed(as.integer(cfg$seed))
    ss <- split(meta$cell_id,paste(meta$donor_id,meta$bridge_group,sep='::'))
    keep <- unlist(lapply(ss,function(z) {
      z <- sort(z);if(length(z)>cfg$smoke$cap_cells_per_donor_group) sample(z,cfg$smoke$cap_cells_per_donor_group) else z
    }),use.names=FALSE)
    x <- x[,keep,drop=FALSE];meta <- meta[keep,,drop=FALSE]
  }
  stopifnot(receiver %in% meta$bridge_group,all(Matrix::colSums(x)>0))
  original_genes <- rownames(x)
  obj <- Seurat::CreateSeuratObject(counts=x,meta.data=meta,min.cells=0,min.features=0)
  gene_map <- data.frame(original_gene=original_genes,seurat_feature=rownames(obj))
  write_tsv(gene_map,file.path(wd,'tables/SEURAT_FEATURE_IDENTITY_MAP.tsv'))
  stopifnot(!anyDuplicated(rownames(obj)),identical(colnames(obj),meta$cell_id))
  obj <- Seurat::NormalizeData(obj,normalization.method='LogNormalize',scale.factor=10000,verbose=FALSE)
  obj <- Seurat::FindVariableFeatures(obj,selection.method='vst',nfeatures=2000,verbose=FALSE)
  Seurat::Idents(obj) <- 'bridge_group'
  datatype <- if(use_pbmc) 'scale.data' else cfg$prepare$datatype
  if (datatype=='scale.data') obj <- Seurat::ScaleData(obj,features=rownames(obj),verbose=FALSE)
  # Frozen target rule; do not insert nominated genes into the selected target list.
  if(use_pbmc) {
    markers <- Seurat::FindMarkers(obj,ident.1=receiver,min.pct=.25,only.pos=TRUE,logfc.threshold=.25)
    fc <- intersect(c('avg_log2FC','avg_logFC'),names(markers))[1]
    if(is.na(fc)) stop('Unrecognized FindMarkers fold-change output')
    markers <- markers[order(markers$p_val_adj,-markers[[fc]],rownames(markers),method='radix'),,drop=FALSE]
    targets <- rownames(markers)[markers$p_val_adj<1e-3]
    targets <- head(targets,as.integer(cfg$smoke$maximum_targets))
    selection <- data.frame(gene=rownames(markers),markers,
                            passes_official_DEG_rule=markers$p_val_adj<1e-3,
                            selected_target=rownames(markers) %in% targets,check.names=FALSE)
    target_rule <- paste0('PBMC technical subset of official receiver DEG rule;min.pct .25;positive;logFC .25;p_adj<1e-3;max',
                         cfg$smoke$maximum_targets,';sort p_adj then descending logFC then gene;no candidate panel')
  } else {
    z <- get_data(obj,'data')[,colnames(obj)[Seurat::Idents(obj)==receiver],drop=FALSE]
    zr <- get_data(obj,'counts')[,colnames(z),drop=FALSE]
    mu <- Matrix::rowMeans(z); z2 <- z;z2@x <- z2@x^2
    vv <- pmax(Matrix::rowMeans(z2)-mu^2,0)
    selection <- data.frame(gene=rownames(z),detection=Matrix::rowMeans(zr>0),mean_logCP10K=mu,variance_logCP10K=vv)
    selection$eligible_expression <- selection$detection>=cfg$target_selection$minimum_detection &
       !grepl('^(MT-|RPL|RPS)',selection$gene)
    selection <- selection[order(-selection$variance_logCP10K,selection$gene,method='radix'),]
    targets <- head(selection$gene[selection$eligible_expression],cfg$target_selection$top_n)
    selection$selected_target <- selection$gene %in% targets
    target_rule <- 'top receiver normalized-expression variance among detection>=.1;exclude MT-/RPL/RPS;no gene panel insertion'
  }
  if (!length(targets)) stop('No targets under fixed selection rule; no target/threshold rescue allowed')
  write_tsv(selection,file.path(wd,'tables/TARGET_GENE_SELECTION_ALL.tsv'))
  writeLines(targets,file.path(wd,'input/targets.txt'))
  write_tsv(meta,file.path(wd,'input/selected_cell_metadata.tsv'))
  write_tsv(as.data.frame(table(donor_id=meta$donor_id,group=meta$bridge_group)),file.path(wd,'tables/DONOR_GROUP_INPUT_COUNTS.tsv'))
  saveRDS(obj,file.path(wd,'r_objects/seurat_new_raw_count_input.rds'))
  # The actual pinned Infer_CCI reads Seurat Idents, not a separate cell_anno argument.
  set.seed(as.integer(cfg$seed))
  cci <- SigXTalkR::Infer_CCI(obj,use_spatial=FALSE,db_use='human')
  saveRDS(cci,file.path(wd,'r_objects/cellchat_full_new.rds'))
  lr_all <- CellChat::subsetCommunication(cci,thresh=1)
  write_tsv(lr_all,file.path(wd,'tables/CELLCHAT_ALL_SOURCE_TARGET_LR.tsv'))
  if(use_pbmc) {
    scoped <- cci
    lr_focus <- lr_all[lr_all$target==receiver,,drop=FALSE]
  } else {
    lr_focus <- lr_all[lr_all$source==cfg$sender & lr_all$target==receiver,,drop=FALSE]
    scoped <- cci
    for(nm in c('prob','pval')) {
      ar <- scoped@net[[nm]];ar[] <- if(nm=='prob') 0 else 1
      ar[cfg$sender,receiver,] <- cci@net[[nm]][cfg$sender,receiver,]
      scoped@net[[nm]] <- ar
    }
  }
  write_tsv(lr_focus,file.path(wd,'tables/FOCUS_LR_COMPLEX_ALL.tsv'))
  saveRDS(scoped,file.path(wd,'r_objects/cellchat_scoped_for_PrepareInput_new.rds'))
  # Preserve original complex identity and split only the explicit interaction identifier components.
  cc <- list()
  if(nrow(lr_focus)) for(i in seq_len(nrow(lr_focus))) {
    ids <- strsplit(as.character(lr_focus$interaction_name[i]),'_',fixed=TRUE)[[1]]
    ligand_parts <- strsplit(as.character(lr_focus$ligand[i]),'_',fixed=TRUE)[[1]]
    recs <- if(length(ids)>length(ligand_parts) && identical(ids[seq_along(ligand_parts)],ligand_parts)) ids[-seq_along(ligand_parts)] else as.character(lr_focus$receptor[i])
    cc[[length(cc)+1L]] <- data.frame(source=lr_focus$source[i],target=receiver,ligand=lr_focus$ligand[i],
      receptor_complex=lr_focus$receptor[i],interaction_name=lr_focus$interaction_name[i],
      receptor_component=recs,component_not_single_receptor_claim=TRUE)
  }
  write_tsv(if(length(cc)) do.call(rbind,cc) else data.frame(),file.path(wd,'tables/LR_COMPLEX_COMPONENTS.tsv'))
  env <- new.env(parent=globalenv())
  utils::data('RTF_human',package='SigXTalkR',envir=env);utils::data('TFT_human',package='SigXTalkR',envir=env)
  rtf <- dplyr::distinct(env$RTF_human,from,to,.keep_all=TRUE)
  tft <- dplyr::distinct(env$TFT_human,from,to,.keep_all=TRUE)
  rtf <- SigXTalkR::Filter_DB(rtf,rownames(obj)); tft <- SigXTalkR::Filter_DB(tft,rownames(obj))
  # Diagnose source limitations without changing upstream ligand/receiver filtering.
  expr <- SigXTalkR::Get_Exp_Clu(obj,clusterID=receiver,assay='RNA',datatype=datatype,
                                cutoff=cfg$prepare$exp_threshold)
  helper_lr <- SigXTalkR::Extract_LR_Prob(scoped,target_type=receiver)
  helper_lr$ligand_receiver_gene_present <- helper_lr$From %in% rownames(expr)
  helper_lr$receptor_receiver_gene_present <- helper_lr$To %in% rownames(expr)
  write_tsv(helper_lr,file.path(wd,'tables/PREPARE_LR_RECEIVER_FILTER_AUDIT.tsv'))
  args <- list(SeuratObj=obj,target_type=receiver,TGs=targets,CCC_results=scoped,RTFDB=rtf,TFTGDB=tft,
    data_dir=file.path(wd,'inputs'),assay='RNA',datatype=datatype,
    exp_threshold=cfg$prepare$exp_threshold,CCC_threshold=cfg$prepare$CCC_threshold,
    Fisher_threshold=cfg$prepare$Fisher_threshold,species='human')
  do.call(SigXTalkR::Prepare_Input,args)
  exp_prs <- SigXTalkR::Get_Exp_Clu(obj,clusterID=receiver,assay='RNA',datatype='data',cutoff=cfg$prs$expression_cutoff)
  saveRDS(exp_prs,file.path(wd,'r_objects/receiver_logCP10K_new.rds'))
  f <- file.path(wd,'inputs',c('ExpressionCount.csv','RecTFDB.csv','TFTGDB.csv'))
  stopifnot(all(file.exists(f)),all(file.info(f)$size>0))
  rp <- read.table(f[2],header=TRUE);tp <- read.table(f[3],header=TRUE)
  if(!nrow(rp) || !nrow(tp)) stop('Official prior filters returned empty inputs')
  candidate <- unique(c('FN1','ITGB1','PLOD2','CCN2','PRSS23','THBS1'))
  tested_tg <- unique(tp$TG)
  rx <- get_data(obj,'counts')[,colnames(obj)[Seurat::Idents(obj)==receiver],drop=FALSE]
  sx <- get_data(obj,'counts')[,colnames(obj)[if(use_pbmc) Seurat::Idents(obj)!=receiver else Seurat::Idents(obj)==cfg$sender],drop=FALSE]
  detection <- function(m,g) if(g %in% rownames(m)) as.numeric(Matrix::rowMeans(m[g,,drop=FALSE]>0)) else NA_real_
  write_tsv(data.frame(gene=candidate,
    query_role=c('ligand_context','receptor_complex_component','program_target_or_readout','program_target_or_readout','program_target_or_readout','comparison_ligand_or_target'),
    in_raw_gene_table=candidate %in% original_genes,
    receiver_counts_detection=vapply(candidate,function(g)detection(rx,g),numeric(1)),
    sender_counts_detection=vapply(candidate,function(g)detection(sx,g),numeric(1)),
    in_selected_targets=candidate %in% targets,in_filtered_TFTG_targets=candidate %in% tested_tg,
    reason=ifelse(!candidate %in% rownames(obj),'not_in_feature_table',
      ifelse(!candidate %in% targets,'not_selected_by_fixed_variance_rule',
        ifelse(!candidate %in% tested_tg,'removed_by_official_expression_or_prior_filter','prepared_for_HGNN')))),
    file.path(wd,'tables/CANDIDATE_TESTABILITY.tsv'))
  rc <- data.frame(file=basename(f),bytes=file.info(f)$size,sha256=vapply(f,sha256_file,character(1)))
  write_tsv(rc,file.path(wd,'audit/HGNN_INPUT_FILE_HASHES.tsv'))
  provenance <- list(receiver=receiver,input_kind=if(use_pbmc) 'official_PBMC_RDS' else 'real_GSE230808_raw_counts',
    n_cells=ncol(obj),n_genes=nrow(obj),n_receiver=sum(Seurat::Idents(obj)==receiver),
    independent_patient_count=length(unique(meta$donor_id)),target_rule=target_rule,n_selected_targets=length(targets),
    actual_datatype=datatype,prepare_parameters=cfg$prepare,RTF_rows_after_gene_filter=nrow(rtf),
    TFT_rows_after_gene_filter=nrow(tft),prepared_RTF_rows=nrow(rp),prepared_TFTG_rows=nrow(tp),
    raw_RTF_prior_object_sha256=digest::digest(env$RTF_human,algo='sha256'),
    raw_TFT_prior_object_sha256=digest::digest(env$TFT_human,algo='sha256'),
    prepared_target_genes=length(tested_tg),PRS_receiver_expression_dim=dim(exp_prs),
    new_input_file_hashes=raw_hashes,config_sha256=sha256_file(cfg$config_file),
    target_type=receiver,project_id=if(cfg$mode=='formal')cfg$project_id else paste0(cfg$project_id,'_SMOKE'),
    started_from_old_pathway_output=FALSE,prior_rows_unmodified=TRUE)
  write_json(provenance,file.path(wd,'audit/NEW_INPUT_PROVENANCE.json'))
  provenance
})
