# Integrated LUAD analysis: bulk, LASSO-Cox, scRNA-seq and Visium
# Configurable template; not the original code of a completed study.
# Requires user data, reviewed annotations and prespecified signatures.
# No nine-gene result or stronger glycolysis association is forced.
set.seed(2026)

# Install dependencies separately:
# install.packages(c('Seurat','glmnet','survival','survminer',
#                    'ggplot2','Matrix','BiocManager'))
# BiocManager::install('GSVA')
packages <- c('Seurat','glmnet','survival','survminer','ggplot2','Matrix','GSVA')
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install: ', paste(missing, collapse=', '))
suppressPackageStartupMessages({
  library(Seurat); library(glmnet); library(survival)
  library(survminer); library(ggplot2); library(Matrix); library(GSVA)
})

# INPUT FORMATS
# Expression CSV: gene symbols in first column, samples in other columns;
# comparable log2(TPM+1) values in training and validation.
# Clinical CSV: sample_id,time,status (positive time; death=1, censored=0).
# sc manifest CSV: sample_id,patient_id,path (10x matrix directory).
# annotations CSV: sample_id,barcode,cell_type,malignancy.
# malignancy labels: Malignant / Non-malignant / Unknown.
# Visium manifest CSV: section_id,patient_id,path (Space Ranger outs).
# GMT: signature-name, description, genes separated by tabs.
# immune_signature_names.txt: one GMT signature name per line.
cfg <- list(
  output='LUAD_results',
  bulk_train='data/bulk_train_log2TPM.csv',
  clinical_train='data/clinical_train.csv',
  bulk_validation='data/bulk_validation_log2TPM.csv',
  clinical_validation='data/clinical_validation.csv',
  signatures='data/signatures.gmt',
  glycolysis_set='HALLMARK_GLYCOLYSIS',
  immune_signature_names='data/immune_signature_names.txt',
  single_cell_manifest='data/single_cell_manifest.csv',
  single_cell_annotations='data/single_cell_annotations.csv',
  spatial_manifest='data/spatial_manifest.csv',
  min_features=200, max_features=7500, max_mt_percent=20,
  spatial_min_features=200, npcs=30, cluster_resolution=0.5,
  lambda_choice='lambda.1se', univariate_fdr=0.05,
  expected_sections=56, expected_patients=25,
  save_spatial_objects=FALSE
)
dir.create(cfg$output, recursive=TRUE, showWarnings=FALSE)
out <- function(name) file.path(cfg$output,name)

read_expression <- function(path) {
  x <- as.matrix(read.csv(path,row.names=1,check.names=FALSE))
  storage.mode(x) <- 'double'
  if(anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)))
    stop('Duplicate gene/sample IDs: ',path)
  if(any(!is.finite(x))) stop('Non-finite expression: ',path)
  x
}
read_clinical <- function(path) {
  d <- read.csv(path,stringsAsFactors=FALSE)
  required <- c('sample_id','time','status')
  stopifnot(all(required %in% names(d)))
  if(anyDuplicated(d$sample_id) || anyNA(d[,required]))
    stop('Duplicate IDs or missing clinical values: ',path)
  if(!is.numeric(d$time) || any(d$time<=0) || !all(d$status %in% c(0,1)))
    stop('Check survival time/status: ',path)
  d
}
align_cohort <- function(expression,clinical) {
  ids <- intersect(clinical$sample_id,colnames(expression))
  if(length(ids)<10) stop('Too few matched samples.')
  clinical <- clinical[match(ids,clinical$sample_id),,drop=FALSE]
  expression <- expression[,ids,drop=FALSE]
  stopifnot(identical(clinical$sample_id,colnames(expression)))
  list(expression=expression,clinical=clinical)
}
read_gmt <- function(path) {
  lines <- readLines(path,warn=FALSE)
  f <- strsplit(lines[nzchar(lines)],'\t',fixed=TRUE)
  if(any(lengths(f)<3)) stop('Malformed GMT.')
  sets <- lapply(f,function(z) unique(z[-c(1,2)]))
  names(sets) <- vapply(f,`[`,character(1),1)
  if(anyDuplicated(names(sets))) stop('Duplicate signature names.')
  sets
}
gene_sets <- read_gmt(cfg$signatures)
stopifnot(cfg$glycolysis_set %in% names(gene_sets))
immune_names <- trimws(readLines(cfg$immune_signature_names,warn=FALSE))
immune_names <- unique(immune_names[nzchar(immune_names)])
stopifnot(length(immune_names)>0,all(immune_names %in% names(gene_sets)),
          !cfg$glycolysis_set %in% immune_names)

# BULK: training-only association reporting and penalized model fitting
training <- align_cohort(read_expression(cfg$bulk_train),
                         read_clinical(cfg$clinical_train))
train_expr <- training$expression
train_clin <- training$clinical
genes <- intersect(gene_sets[[cfg$glycolysis_set]],rownames(train_expr))
sds <- apply(train_expr[genes,,drop=FALSE],1,sd)
genes <- genes[is.finite(sds) & sds>0]
if(length(genes)<2) stop('Too few variable glycolysis genes.')
univariate <- do.call(rbind,lapply(genes,function(g) {
  d <- data.frame(time=train_clin$time,status=train_clin$status,
                  expression=as.numeric(train_expr[g,]))
  tryCatch({
    s <- summary(coxph(Surv(time,status)~expression,data=d))
    data.frame(gene=g,HR=s$conf.int[1,'exp(coef)'],
               lower95=s$conf.int[1,'lower .95'],
               upper95=s$conf.int[1,'upper .95'],
               p=s$coefficients[1,'Pr(>|z|)'])
  },error=function(e) data.frame(gene=g,HR=NA_real_,lower95=NA_real_,
                                upper95=NA_real_,p=NA_real_))
}))
univariate$FDR <- p.adjust(univariate$p,'BH')
univariate <- univariate[order(univariate$FDR),]
write.csv(univariate,out('bulk_univariate_cox.csv'),row.names=FALSE)
write.csv(subset(univariate,!is.na(FDR)&FDR<cfg$univariate_fdr),
          out('bulk_survival_associated_genes.csv'),row.names=FALSE)
# All eligible genes enter LASSO; univariate screening is not used.
x_train <- t(train_expr[genes,,drop=FALSE])
y_train <- with(train_clin,Surv(time,status))
nfolds <- min(10L,floor(nrow(train_clin)/5),sum(train_clin$status==1))
if(nfolds<3) stop('Insufficient samples/events for CV.')
foldid <- integer(nrow(train_clin))
for(event_status in c(0,1)) {
  idx <- which(train_clin$status==event_status)
  if(length(idx)) foldid[idx] <- sample(rep(seq_len(nfolds),length.out=length(idx)))
}
cvfit <- cv.glmnet(x_train,y_train,family='cox',alpha=1,
                   foldid=foldid,standardize=TRUE,type.measure='deviance')
pdf(out('lasso_cross_validation.pdf')); plot(cvfit); dev.off()
cm <- as.matrix(coef(cvfit,s=cfg$lambda_choice))
selected <- data.frame(gene=rownames(cm),coefficient=as.numeric(cm[,1]))
selected <- selected[selected$coefficient!=0,]
write.csv(selected,out('lasso_selected_genes.csv'),row.names=FALSE)
if(!nrow(selected)) stop('No genes retained at prespecified lambda.')
message('Selected genes: ',nrow(selected))
train_clin$risk_score <- as.numeric(predict(cvfit,newx=x_train,
                                   s=cfg$lambda_choice,type='link'))
risk_cutoff <- median(train_clin$risk_score)
assign_group <- function(score) factor(ifelse(score>risk_cutoff,'High','Low'),
                                       levels=c('Low','High'))
train_clin$risk_group <- assign_group(train_clin$risk_score)

# EXTERNAL VALIDATION: fixed coefficients and training cutoff; no refitting
validation <- align_cohort(read_expression(cfg$bulk_validation),
                           read_clinical(cfg$clinical_validation))
val_expr <- validation$expression
val_clin <- validation$clinical
if(length(intersect(train_clin$sample_id,val_clin$sample_id)))
  stop('Training/validation sample IDs overlap.')
stopifnot(all(selected$gene %in% rownames(val_expr)))
val_clin$risk_score <- as.numeric(t(val_expr[selected$gene,,drop=FALSE]) %*%
                                   selected$coefficient)
val_clin$risk_group <- assign_group(val_clin$risk_score)
evaluate_survival <- function(d,label) {
  write.csv(d,out(paste0(label,'_risk_scores.csv')),row.names=FALSE)
  fit <- coxph(Surv(time,status)~risk_score,data=d)
  s <- summary(fit)
  write.csv(data.frame(cohort=label,
    HR_per_score_unit=s$conf.int[1,'exp(coef)'],
    lower95=s$conf.int[1,'lower .95'],upper95=s$conf.int[1,'upper .95'],
    p=s$coefficients[1,'Pr(>|z|)'],C_index=unname(s$concordance[1]),
    C_index_SE=unname(s$concordance[2])),
    out(paste0(label,'_survival_metrics.csv')),row.names=FALSE)
  capture.output(cox.zph(fit),file=out(paste0(label,'_PH_diagnostic.txt')))
  if(length(unique(d$risk_group))==2) {
    km <- survfit(Surv(time,status)~risk_group,data=d)
    p <- ggsurvplot(km,data=d,pval=TRUE,risk.table=TRUE,
                    legend.title='Risk group',title=paste(label,'cohort'))
    pdf(out(paste0(label,'_KM.pdf')),width=8,height=8); print(p); dev.off()
  }
}
evaluate_survival(train_clin,'training_apparent')
evaluate_survival(val_clin,'external_validation')
saveRDS(list(cvfit=cvfit,lambda_choice=cfg$lambda_choice,
             selected_genes=selected,cutoff=risk_cutoff,foldid=foldid),
         out('lasso_cox_model.rds'))

# SINGLE-CELL: baseline Seurat workflow; evaluate doublets/batches separately
sm <- read.csv(cfg$single_cell_manifest,stringsAsFactors=FALSE)
stopifnot(all(c('sample_id','patient_id','path') %in% names(sm)),
          !anyDuplicated(sm$sample_id),nrow(sm)>0)
objects <- lapply(seq_len(nrow(sm)),function(i) {
  counts <- Read10X(sm$path[i])
  if(is.list(counts)) {
    stopifnot('Gene Expression' %in% names(counts))
    counts <- counts[['Gene Expression']]
  }
  obj <- CreateSeuratObject(counts=counts,project=sm$sample_id[i])
  obj$sample_id <- sm$sample_id[i]; obj$patient_id <- sm$patient_id[i]
  obj$original_barcode <- colnames(obj)
  obj[['percent.mt']] <- PercentageFeatureSet(obj,pattern='^MT-')
  obj <- subset(obj,subset=nFeature_RNA>=cfg$min_features &
                 nFeature_RNA<=cfg$max_features & percent.mt<=cfg$max_mt_percent)
  if(ncol(obj)<20) stop('Too few cells after QC: ',sm$sample_id[i])
  RenameCells(obj,add.cell.id=sm$sample_id[i])
})
sc <- if(length(objects)==1) objects[[1]] else merge(objects[[1]],y=objects[-1])
sc <- JoinLayers(sc,assay='RNA')
sc <- NormalizeData(sc,verbose=FALSE)
sc <- FindVariableFeatures(sc,nfeatures=2000,verbose=FALSE)
sc <- ScaleData(sc,features=VariableFeatures(sc),verbose=FALSE)
npcs <- min(cfg$npcs,ncol(sc)-1,length(VariableFeatures(sc))-1)
if(npcs<2) stop('Too few cells/features for PCA.')
sc <- RunPCA(sc,npcs=npcs,verbose=FALSE)
sc <- FindNeighbors(sc,dims=seq_len(npcs),verbose=FALSE)
sc <- FindClusters(sc,resolution=cfg$cluster_resolution,verbose=FALSE)
sc <- RunUMAP(sc,dims=seq_len(npcs),verbose=FALSE)
ggsave(out('single_cell_UMAP_samples.pdf'),DimPlot(sc,group.by='sample_id'),
       width=8,height=6)
markers <- FindAllMarkers(sc,only.pos=TRUE,min.pct=0.25,
                         logfc.threshold=0.25,verbose=FALSE)
write.csv(markers,out('single_cell_cluster_markers.csv'),row.names=FALSE)

# Reviewed annotations are inputs; malignancy is not assigned from GPR87.
a <- read.csv(cfg$single_cell_annotations,stringsAsFactors=FALSE)
stopifnot(all(c('sample_id','barcode','cell_type','malignancy') %in% names(a)))
key <- paste(a$sample_id,a$barcode,sep='::')
stopifnot(!anyDuplicated(key))
idx <- match(paste(sc$sample_id,sc$original_barcode,sep='::'),key)
if(anyNA(idx)) stop('Complete reviewed annotations using exported markers.')
sc$cell_type <- a$cell_type[idx]; sc$malignancy <- a$malignancy[idx]
stopifnot(!anyNA(sc$cell_type),!anyNA(sc$malignancy),
          all(sc$malignancy %in% c('Malignant','Non-malignant','Unknown')),
          'GPR87' %in% rownames(sc))
ggsave(out('single_cell_UMAP_cell_types.pdf'),
       DimPlot(sc,group.by='cell_type',label=TRUE),width=9,height=7)
ggsave(out('single_cell_GPR87_UMAP.pdf'),FeaturePlot(sc,features='GPR87'),
       width=8,height=6)
ggsave(out('single_cell_GPR87_malignancy.pdf'),
       VlnPlot(sc,features='GPR87',group.by='malignancy',pt.size=0),width=8,height=5)
gpr <- FetchData(sc,vars='GPR87',layer='data')[,1]
d <- data.frame(patient_id=sc$patient_id,cell_type=sc$cell_type,
                malignancy=sc$malignancy,expression=gpr,detected=as.numeric(gpr>0))
means <- aggregate(expression~patient_id+cell_type+malignancy,d,mean)
detection <- aggregate(detected~patient_id+cell_type+malignancy,d,mean)
write.csv(merge(means,detection,by=c('patient_id','cell_type','malignancy')),
          out('single_cell_GPR87_patient_summaries.csv'),row.names=FALSE)
saveRDS(sc,out('single_cell_processed.rds'))

# VISIUM: normalized expression maps and descriptive within-section rho
spatial_sets <- lapply(gene_sets[c(cfg$glycolysis_set,immune_names)],
                       function(g) setdiff(g,'GPR87'))
safe_rho <- function(x,y) {
  k <- is.finite(x)&is.finite(y)
  if(sum(k)<10 || sd(x[k])==0 || sd(y[k])==0) return(NA_real_)
  cor(x[k],y[k],method='spearman')
}
vm <- read.csv(cfg$spatial_manifest,stringsAsFactors=FALSE)
stopifnot(all(c('section_id','patient_id','path') %in% names(vm)),
          !anyDuplicated(vm$section_id),nrow(vm)>0)
message('Visium: ',nrow(vm),' sections; ',length(unique(vm$patient_id)),' patients')
if(nrow(vm)!=cfg$expected_sections ||
   length(unique(vm$patient_id))!=cfg$expected_patients)
  warning('Manifest does not match expected 56 sections / 25 patients.')
correlations <- qc <- vector('list',nrow(vm))
for(i in seq_len(nrow(vm))) {
  section <- vm$section_id[i]; patient <- vm$patient_id[i]
  safe_section <- gsub('[^A-Za-z0-9_-]','_',section)
  message('Processing section: ',section)
  sp <- Load10X_Spatial(data.dir=vm$path[i],
    filename='filtered_feature_bc_matrix.h5',assay='Spatial',
    slice=safe_section,filter.matrix=TRUE)
  before <- ncol(sp)
  sp <- subset(sp,subset=nFeature_Spatial>=cfg$spatial_min_features)
  if(ncol(sp)<10 || !'GPR87' %in% rownames(sp))
    stop('Insufficient spots or GPR87 missing: ',section)
  sp <- NormalizeData(sp,assay='Spatial',verbose=FALSE)
  DefaultAssay(sp) <- 'Spatial'
  ggsave(out(paste0('Visium_',safe_section,'_GPR87.pdf')),
    SpatialFeaturePlot(sp,features='GPR87',min.cutoff=0,max.cutoff='q95'),
    width=7,height=7)
  expression <- as.matrix(GetAssayData(sp,assay='Spatial',layer='data'))
  sets <- lapply(spatial_sets,function(g) intersect(g,rownames(expression)))
  if(any(lengths(sets)<5)) stop('Too few signature genes: ',section)
  scores <- as.matrix(GSVA::gsva(
    GSVA::ssgseaParam(expression,sets,minSize=5,normalize=TRUE),verbose=FALSE))
  gpr87 <- as.numeric(expression['GPR87',])
  correlations[[i]] <- data.frame(section_id=section,patient_id=patient,
    signature=rownames(scores),rho=vapply(seq_len(nrow(scores)),
      function(j) safe_rho(gpr87,scores[j,]),numeric(1)))
  write.csv(data.frame(spot_id=colnames(expression),GPR87=gpr87,
                       t(scores),check.names=FALSE),
    out(paste0('Visium_',safe_section,'_ssGSEA.csv')),row.names=FALSE)
  qc[[i]] <- data.frame(section_id=section,patient_id=patient,
    spots_before_QC=before,spots_after_QC=ncol(sp),
    GPR87_detected_fraction=mean(gpr87>0))
  if(cfg$save_spatial_objects) saveRDS(sp,out(paste0('Visium_',safe_section,'.rds')))
  rm(sp,expression,scores); gc(verbose=FALSE)
}
section_correlations <- do.call(rbind,correlations)
write.csv(section_correlations,out('spatial_section_correlations.csv'),row.names=FALSE)
write.csv(do.call(rbind,qc),out('spatial_QC_summary.csv'),row.names=FALSE)

# PATIENT-LEVEL comparisons; no spot-level inferential p-values
valid <- subset(section_correlations,is.finite(rho))
if(!nrow(valid)) stop('No estimable spatial correlations.')
valid$z <- atanh(pmax(pmin(valid$rho,0.999999),-0.999999))
pc <- aggregate(z~patient_id+signature,valid,mean)
pc$rho <- tanh(pc$z)
write.csv(pc,out('spatial_patient_correlations.csv'),row.names=FALSE)
wide <- reshape(pc[,c('patient_id','signature','rho')],
                idvar='patient_id',timevar='signature',direction='wide')
comparisons <- do.call(rbind,lapply(immune_names,function(sig) {
  x <- wide[[paste0('rho.',cfg$glycolysis_set)]]
  y <- wide[[paste0('rho.',sig)]]
  k <- is.finite(x)&is.finite(y); delta <- x[k]-y[k]
  p <- if(length(delta)<5) NA_real_ else if(all(delta==0)) 1 else
    suppressWarnings(wilcox.test(delta,mu=0,alternative='two.sided',exact=FALSE)$p.value)
  data.frame(immune_signature=sig,n_patients=length(delta),
    median_glycolysis_rho=if(length(delta)) median(x[k]) else NA_real_,
    median_immune_rho=if(length(delta)) median(y[k]) else NA_real_,
    median_rho_difference=if(length(delta)) median(delta) else NA_real_,p=p)
}))
comparisons$FDR <- p.adjust(comparisons$p,'BH')
write.csv(comparisons,out('glycolysis_vs_immune_patient_comparisons.csv'),row.names=FALSE)
p <- ggplot(pc,aes(signature,rho))+geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.12,alpha=0.65)+geom_hline(yintercept=0,linetype=2)+
  coord_flip()+theme_bw()+labs(x=NULL,
    y='Patient-level summarized Spearman correlation',
    title='Spatial GPR87 associations')
ggsave(out('spatial_GPR87_signature_correlations.pdf'),p,width=9,height=6)
saveRDS(cfg,out('analysis_configuration.rds'))
capture.output(sessionInfo(),file=out('sessionInfo.txt'))
message('Completed. Outputs: ',normalizePath(cfg$output))

# LIMITATIONS
# Training metrics are apparent performance, not unbiased validation.
# Ensure independent patients and comparable bulk preprocessing.
# Review doublets, batch effects, annotations and malignancy evidence.
# Visium spots contain mixed cells; correlations are exploratory and
# unadjusted for tumor content/confounders. Patient summaries use equal
# section weighting in Fisher-z space; compare correlations directly.
# No supplied datasets were used to execute this template.
