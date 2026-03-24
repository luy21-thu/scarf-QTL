# ================================================================
# This script reformats the released OneK1K resources into per–cell-type /
# per–chromosome inputs for scarf-QTL analyses.
# All data (expression, metadata/covariates, and genotype dosages) 
# are read directly from OneK1K_Angli_SeuratObj.rds, *_peer_factors.tsv and
# chr*.dose.filtered.R2_0.8.vcf.gz, then split/aligned and saved.
# The only additional processing beyond reformatting is pseudotime inference
# within each cell type (PHATE + Slingshot), saved for downstream analyses.
# ================================================================

setwd("~/sceQTL")
options(max.print = 100)
library(Seurat)
library(phateR)
library(slingshot)
library(ggplot2)
library(patchwork)
library(stringr)

# Cell-type labels and file naming conventions
{
  CT_14_ <- c("B_IN", "B_MEM", "CD4_ET", "CD4_NC", "CD4_SOX4", 
              "CD8_ET", "CD8_NC", "CD8_S100B", "DC", 
              "Mono_C", "Mono_NC", "NK", "NK_R", "Plasma")
  CT_14 <- c("B IN", "B Mem", "CD4 ET", "CD4 NC", "CD4 SOX4", 
             "CD8 ET", "CD8 NC", "CD8 S100B", "DC", 
             "Mono C", "Mono NC", "NK", "NK R", "Plasma")
  CT_14_PFs <- c("BimmNaive", "Bmem", "CD4effCM", "CD4all", "CD4TGFbStim", 
                 "CD8eff", "CD8all", "CD8unknown", "DC", 
                 "MonoC", "MonoNC", "NKmat", "NKact", "Plasma")
}

# Build (gene, cis-SNP) lists
{
  for (CT in CT_14) {
    ct <- stringr::str_replace(tolower(CT), " ", "")
    eQTL_Yazar_CT <- readr::read_tsv(sprintf("data/eQTL_Yazar/%s_eqtl_table.tsv.gz", ct))
    {
      Gene_CT <- eQTL_Yazar_CT[!duplicated(eQTL_Yazar_CT$GENE), c("GENE", "CHR")]
      Gene_CT <- tapply(Gene_CT$GENE, Gene_CT$CHR, c, simplify = F)
      saveRDS(Gene_CT, file = sprintf(
        "data/Gene_cis_SNP/Gene_%s.rds", ct))
    }
    {
      eQTL_Yazar_CT <- eQTL_Yazar_CT[eQTL_Yazar_CT$ROUND == 1,]
      for (chr in seq(22)) {
        eQTL_Yazar_CT_chr <- eQTL_Yazar_CT[eQTL_Yazar_CT$CHR == chr,]
        Gene_SNP_list <- tapply(eQTL_Yazar_CT_chr$SNPID, eQTL_Yazar_CT_chr$GENE, c, simplify = F)
        saveRDS(Gene_SNP_list, file = sprintf(
          "data/Gene_cis_SNP/%s_chr%s.rds", ct, chr))
      }
    }
  }
}

# Split and save scRNA-seq data
{
  library(Seurat)
  RNA_all <- readRDS("/PublicData/Onek1k/angli/OneK1K_Angli_SeuratObj.rds")
  RNA_all <- UpdateSeuratObject(RNA_all)
  
  cell_meta <- RNA_all@meta.data
  scRNA_SCT <- RNA_all@assays$SCT@data
  scRNA_PCA <- RNA_all@reductions$pca@cell.embeddings
  saveRDS(cell_meta, file = "data/scRNA_SCT/Cell_meta.rds")
  saveRDS(scRNA_SCT, file = "data/scRNA_SCT/scRNA_SCT.rds")
  saveRDS(scRNA_PCA, file = "data/scRNA_SCT/PCA.rds")
}

# Prepare covariates
{
  cell_meta <- readRDS("data/scRNA_SCT/Cell_meta.rds")
  cell_meta$cell_label[cell_meta$cell_label == "B MEM"] <- "B Mem"
  covariates_Cell <- cell_meta[,c("log_nUMI", "percent.mt", "individual")]
  remove_2_individual <- covariates_Cell$individual %in% c("88_88","966_967")
  for (CT_idx in seq(14)) {
    CT <- CT_14[CT_idx]
    CT_PF <- CT_14_PFs[CT_idx]
    ct <- str_replace(tolower(CT), " ", "")
    covariates_cell <- covariates_Cell[(cell_meta$cell_label == CT) & (!remove_2_individual),]
    saveRDS(covariates_cell, file = sprintf(
      "data/covariate/cell_%s.rds", ct))
    
    indi_meta <- read.table(sprintf("data/covariate/OneK1K_covariates_PCs_PFs_14_cell_types-1/%s_peer_factors.tsv", CT_PF), header = T)
    indi_meta <- indi_meta[!indi_meta$sampleid %in% c("88_88","966_967"),]
    covariates_indi <- indi_meta[,-1]; rownames(covariates_indi) <- indi_meta[,1]
    saveRDS(covariates_indi, file = sprintf("data/covariate/individual_%s.rds", ct))
  }
}

# Save expression matrices by chromosome and cell type
{
  scRNA_SCT <- readRDS("data/scRNA_SCT/scRNA_SCT.rds")
  for (CT in CT_14) {
    scRNA_SCT_CT <- scRNA_SCT[,(cell_meta$cell_label == CT) & (!remove_2_individual), drop = F]
    cat(ncol(scRNA_SCT_CT), "")
    ct <- str_replace(tolower(CT), " ", "")
    Gene_CT <- readRDS(sprintf("data/Gene_cis_SNP/Gene_%s.rds", ct))
    for (chr in seq(22)) {
      saveRDS(scRNA_SCT_CT[Gene_CT[[chr]], , drop = F], file = sprintf(
        "data/scRNA_SCT/scRNA_SCT_%s_chr%s.rds", ct, chr))
    }
  }
}

# Infer pseudotime
{
  RNA_all <- readRDS("/PublicData/Onek1k/angli/OneK1K_Angli_SeuratObj.rds")
  RNA_all <- UpdateSeuratObject(RNA_all)
  remove_2_individual <- RNA_all$individual %in% c("88_88","966_967")
  RNA_all <- RNA_all[,!remove_2_individual]
  
  for (CT_ in CT_14_) {
    ct <- stringr::str_replace(tolower(CT_), "_", "")
    RNA_CT <- RNA_all[,RNA_all$new_names == CT_]
    system.time(RNA_CT <- SCTransform(RNA_CT, 
                                      vars.to.regress = c("pool", "percent.mt"), 
                                      do.correct.umi = FALSE, 
                                      variable.features.n = 500))
    
    em <- t(RNA_CT@assays$SCT@scale.data)
    system.time(ph <- phate(em, npca = 10))
    em_ph <- ph[["embedding"]]
    
    sds <- slingshot(em_ph)
    pt <- slingCurves(sds)[[1]]$lambda
    
    covariates_cell <- readRDS(sprintf("data/covariate/cell_%s.rds", ct))
    cat(CT_, identical(rownames(covariates_cell), names(pt)))
    covariates_cell$pseudotime <- pt
    saveRDS(sds, file = sprintf("data/pseudotime/sds_%s.rds", CT_))
    saveRDS(em, file = sprintf("data/pseudotime/em_%s.rds", CT_))
    saveRDS(covariates_cell, file = sprintf("data/covariate/cell_pt_%s.rds", ct))
  }
}

# Reformat genotype dosages
for (chr in 1:22) {
  VDF_file <- paste0("data/genotype/filter_vcf_r08/chr", 
                     chr, ".dose.filtered.R2_0.8.vcf.gz")
  VCF_annot <- readLines(VDF_file, 17)
  VCF_raw <- read.table(VDF_file)
  colnames(VCF_raw) <- strsplit(VCF_annot[17],"\t")[[1]]
  VCF_data <- VCF_raw[,10:ncol(VCF_raw)]
  VCF_genotype <- sapply(VCF_data, function(GT) {
    as.numeric(substr(GT, 1, 1)) + as.numeric(substr(GT, 3, 3))
  })
  VCF_SNP_info <- VCF_raw[,1:9]
  rownames(VCF_genotype) <- paste0(VCF_SNP_info$ID, "_", VCF_SNP_info$ALT)
  saveRDS(VCF_genotype, file =
            sprintf("data/genotype/chr%s_genotype_matrix_rename.rds", chr))
  saveRDS(VCF_SNP_info, file =
            sprintf("data/genotype/chr%s_SNP_info.rds", chr))
  cat(chr, ncol(VCF_genotype), nrow(VCF_genotype), "\n")
}

# Collect rsIDs
{
  SNP_info_list <- parallel::mclapply(CT_14, function(CT) {
    ct <- str_replace(tolower(CT), " ", "")
    eQTL_Yazar_CT <- readr::read_tsv(sprintf("data/eQTL_Yazar/%s_eqtl_table.tsv.gz", ct))
    df <- dplyr::select(eQTL_Yazar_CT, SNPID, RSID, A1, A2, A2_FREQ_ONEK1K, A2_FREQ_HRC)
    df <- unique(df)
  }, mc.cores = 14)
  SNP_df <- do.call(rbind, SNP_info_list)
  SNP_df <- SNP_df[!duplicated(SNP_df$SNPID), ]
  saveRDS(SNP_df, file = "data/genotype/rsid_AF.rds")
  saveRDS(SNP_df[,1:2], file = "data/genotype/rsid.rds")
}