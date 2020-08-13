


extension <- gsub("(.{5}).+", "\\1", digest::digest(rnorm(1)))
print(paste0("Extension: ", extension))

set.seed(1)


library(SingleCellExperiment)
library(DESeq2)
library(glmGamPoi)
library(edgeR)

register(BiocParallel::SerialParam())
setAutoBPPARAM(BiocParallel::SerialParam())
# register(BiocParallel::MulticoreParam(12))
# setAutoBPPARAM(BiocParallel::MulticoreParam(12))
setwd("/g/huber/users/ahlmanne/projects/glmGamPoi-Paper")

print("Read Dataset")
sce <- TENxPBMCData::TENxPBMCData("pbmc4k")
print(dim(sce))

# # To speed things up for now
# sce <- sce[sample(seq_len(nrow(sce)), size = 1000), ]

# Make sure matrix is in memory
assay(sce) <- as.matrix(assay(sce))

# Calculate size factors

size_factors <- scran::calculateSumFactors(assay(sce))


# cell number
cell_numbers <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000)

print("Testing different cell numbers")

for(cell_number in cell_numbers){
  
  print(paste0("Cell Number: ", cell_number))
  selection <- sample(seq_len(ncol(sce)), cell_number, replace = FALSE)
  working_sce <- sce[, selection]
  working_size_factors <- size_factors[selection]
  
  print("Fit glm_gp overdispersion = 0")
  
  print(system.time({
    fit_disp0 <- glm_gp(working_sce, design =  ~ 1,
                        overdispersion = 0,
                        size_factors = working_size_factors, verbose = TRUE)
  }))
  
  print("Fit glm_gp")
  
  print(system.time({
    fit <- glm_gp(working_sce, design =  ~ 1,
                  size_factors = working_size_factors, verbose = TRUE)
  }))
  
  print("Fit glm_gp on disk")
  
  print(system.time({
    fit_ondisk <- glm_gp(working_sce, design =  ~ 1,
                         size_factors = working_size_factors, verbose = TRUE,
                         on_disk = TRUE)
  }))
  
  print("Fit DESeq2 with glmGamPoi")
  
  print(system.time({
    dds_gp <- DESeq2::DESeqDataSet(working_sce, design =  ~ 1)
    sizeFactors(dds_gp) <- working_size_factors
    dds_gp <- DESeq2::estimateDispersions(dds_gp, minmu = 1e-6, fitType = "glmGamPoi")
    dds_gp <- DESeq2::nbinomWaldTest(dds_gp, minmu = 1e-6)
  }))
  
  print("Fit edgeR")
  
  print(system.time({
    model_matrix <- fit$model_matrix
    edgeR_data <- edgeR::DGEList(assay(working_sce), lib.size = working_size_factors)
    edgeR_data <- edgeR::estimateDisp(edgeR_data, design = model_matrix)
    edgeR_fit <- edgeR::glmQLFit(edgeR_data, design = model_matrix)
  }))
  
  print("Fit DESeq2")
  
  print(system.time({
    dds <- DESeq2::DESeqDataSet(working_sce, design =  ~ 1)
    sizeFactors(dds) <- working_size_factors
    dds <- DESeq2::estimateDispersions(dds, minmu = 1e-6)
    dds <- DESeq2::nbinomWaldTest(dds, minmu = 1e-6)
  }))
  
}

gene_numbers <- c(3000, 6000, 9000, 12000, 15000, 18000, 21000, 24000, 27000, 30000)
print("Testing different gene numbers")

for(gene_number in gene_numbers){
  
  print(paste0("Gene Number: ", gene_number))
  selection <- sample(seq_len(nrow(sce)), gene_number, replace = FALSE)
  working_sce <- sce[selection,]
  
  print("Fit glm_gp overdispersion = 0")
  
  print(system.time({
    fit_disp0 <- glm_gp(working_sce, design =  ~ 1,
                        overdispersion = 0,
                        size_factors = size_factors, verbose = TRUE)
  }))
  
  print("Fit glm_gp")
  
  print(system.time({
    fit <- glm_gp(working_sce, design =  ~ 1,
                  size_factors = size_factors, verbose = TRUE)
  }))
  
  print("Fit glm_gp on disk")
  
  print(system.time({
    fit_ondisk <- glm_gp(working_sce, design =  ~ 1,
                         size_factors = size_factors, verbose = TRUE,
                         on_disk = TRUE)
  }))
  
  print("Fit DESeq2 with glmGamPoi")
  
  print(system.time({
    dds_gp <- DESeq2::DESeqDataSet(working_sce, design =  ~ 1)
    sizeFactors(dds_gp) <- size_factors
    dds_gp <- DESeq2::estimateDispersions(dds_gp, minmu = 1e-6, fitType = "glmGamPoi")
    dds_gp <- DESeq2::nbinomWaldTest(dds_gp, minmu = 1e-6)
  }))

  print("Fit edgeR")
  
  print(system.time({
    model_matrix <- fit$model_matrix
    edgeR_data <- edgeR::DGEList(assay(working_sce), lib.size = size_factors)
    edgeR_data <- edgeR::estimateDisp(edgeR_data, design = model_matrix)
    edgeR_fit <- edgeR::glmQLFit(edgeR_data, design = model_matrix)
  }))
  
  print("Fit DESeq2")
  
  print(system.time({
    dds <- DESeq2::DESeqDataSet(working_sce, design =  ~ 1)
    sizeFactors(dds) <- size_factors
    dds <- DESeq2::estimateDispersions(dds, minmu = 1e-6)
    dds <- DESeq2::nbinomWaldTest(dds, minmu = 1e-6)
  }))

}

# Study the complexity wrt the model matrix complexity

helper_df <- data.frame(
  cat1 = sample(LETTERS[1:3], ncol(sce), replace = TRUE),
  cat2 = sample(LETTERS[3:4], ncol(sce), replace = TRUE),
  cat3 = sample(LETTERS[5:6], ncol(sce), replace = TRUE),
  cat4 = sample(LETTERS[7:8], ncol(sce), replace = TRUE),
  cat5 = sample(LETTERS[9:10], ncol(sce), replace = TRUE),
  cat6 = sample(LETTERS[11:12], ncol(sce), replace = TRUE),
  cont1 = rnorm(n = ncol(sce)),
  cont2 = rnorm(n = ncol(sce)),
  cont3 = rnorm(n = ncol(sce)),
  cont4 = rnorm(n = ncol(sce)),
  cont5 = rnorm(n = ncol(sce)),
  cont6 = rnorm(n = ncol(sce))
)

working_model_matrix_cat <- model.matrix(~ cat1 + cat2 + cat3 + cat4 + cat5 + cat6, data = helper_df)
print("Starting cat model subset")


for(model_matrix_subsetter in seq_len(8)){
  
  print(paste0("Model matrix subsetter ", model_matrix_subsetter))
  
  print("Fit glm_gp overdispersion = 0")
  
  working_model_matrix <- working_model_matrix_cat[,seq_len(model_matrix_subsetter),drop=FALSE]
  
  print(system.time({
    fit_disp0 <- glm_gp(sce, design =  working_model_matrix,
                        overdispersion = 0,
                        size_factors = size_factors, verbose = TRUE)
  }))
  
  print("Fit glm_gp")
  
  print(system.time({
    fit <- glm_gp(sce, design =  working_model_matrix,
                  size_factors = size_factors, verbose = TRUE)
  }))
  
  print("Fit glm_gp on disk")
  
  print(system.time({
    fit_ondisk <- glm_gp(sce, design =  working_model_matrix,
                         size_factors = size_factors, verbose = TRUE,
                         on_disk = TRUE)
  }))
  
  print("Fit DESeq2 with glmGamPoi")
  
  print(system.time({
    dds_gp <- DESeq2::DESeqDataSet(sce, design =  working_model_matrix)
    sizeFactors(dds_gp) <- size_factors
    dds_gp <- DESeq2::estimateDispersions(dds_gp, minmu = 1e-6, fitType = "glmGamPoi")
    dds_gp <- DESeq2::nbinomWaldTest(dds_gp, minmu = 1e-6, modelMatrix = working_model_matrix)
  }))
  
  print("Fit edgeR")
  
  print(system.time({
    model_matrix <- working_model_matrix
    edgeR_data <- edgeR::DGEList(assay(sce), lib.size = size_factors)
    edgeR_data <- edgeR::estimateDisp(edgeR_data, design = model_matrix)
    edgeR_fit <- edgeR::glmQLFit(edgeR_data, design = model_matrix)
  }))
  
  print("Fit DESeq2")
  
  print(system.time({
    dds <- DESeq2::DESeqDataSet(sce, design =  working_model_matrix)
    sizeFactors(dds) <- size_factors
    dds <- DESeq2::estimateDispersions(dds, minmu = 1e-6)
    dds <- DESeq2::nbinomWaldTest(dds, minmu = 1e-6, modelMatrix = working_model_matrix)
  }))
}






working_model_matrix_cont <- model.matrix(~ cont1 + cont2 + cont3 + cont4 + cont5 + cont6, data = helper_df)
print("Starting cont model subset")

for(model_matrix_subsetter in seq_len(7)){
  
  print(paste0("Model matrix subsetter ", model_matrix_subsetter))
  
  print("Fit glm_gp overdispersion = 0")
  
  working_model_matrix <- working_model_matrix_cont[,seq_len(model_matrix_subsetter),drop=FALSE]
  
  print(system.time({
    fit_disp0 <- glm_gp(sce, design =  working_model_matrix,
                        overdispersion = 0,
                        size_factors = size_factors, verbose = TRUE)
  }))
  
  print("Fit glm_gp")
  
  print(system.time({
    fit <- glm_gp(sce, design =  working_model_matrix,
                  size_factors = size_factors, verbose = TRUE)
  }))
  
  print("Fit glm_gp on disk")
  
  print(system.time({
    fit_ondisk <- glm_gp(sce, design =  working_model_matrix,
                         size_factors = size_factors, verbose = TRUE,
                         on_disk = TRUE)
  }))
  
  print("Fit DESeq2 with glmGamPoi")
  
  print(system.time({
    dds_gp <- DESeq2::DESeqDataSet(sce, design =  working_model_matrix)
    sizeFactors(dds_gp) <- size_factors
    dds_gp <- DESeq2::estimateDispersions(dds_gp, minmu = 1e-6, fitType = "glmGamPoi")
    dds_gp <- DESeq2::nbinomWaldTest(dds_gp, minmu = 1e-6, modelMatrix = working_model_matrix)
  }))
  
  print("Fit edgeR")
  
  print(system.time({
    model_matrix <- working_model_matrix
    edgeR_data <- edgeR::DGEList(assay(sce), lib.size = size_factors)
    edgeR_data <- edgeR::estimateDisp(edgeR_data, design = model_matrix)
    edgeR_fit <- edgeR::glmQLFit(edgeR_data, design = model_matrix)
  }))
  
  print("Fit DESeq2")
  
  print(system.time({
    dds <- DESeq2::DESeqDataSet(sce, design =  working_model_matrix)
    sizeFactors(dds) <- size_factors
    dds <- DESeq2::estimateDispersions(dds, minmu = 1e-6)
    dds <- DESeq2::nbinomWaldTest(dds, minmu = 1e-6, modelMatrix = working_model_matrix)
  }))
}


print("Finished")


# --------------
sessionInfo()

