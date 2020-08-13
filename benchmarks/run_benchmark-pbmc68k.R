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
sce <- TENxPBMCData::TENxPBMCData("pbmc68k")
print(dim(sce))

# # To speed things up for now
# sce <- sce[sample(seq_len(nrow(sce)), size = 1000), ]

# Make sure matrix is in memory
assay(sce) <- as.matrix(assay(sce))

# Calculate size factors

size_factors <- scran::calculateSumFactors(assay(sce))

print("Fit DESeq2 with glmGamPoi")

system.time({
  dds_gp <- DESeq2::DESeqDataSet(sce, design =  ~ 1)
  sizeFactors(dds_gp) <- size_factors
  dds_gp <- DESeq2::estimateDispersions(dds_gp, minmu = 1e-6, fitType = "glmGamPoi")
  dds_gp <- DESeq2::nbinomWaldTest(dds_gp, minmu = 1e-6)
})

print("Fit glm_gp")

system.time({
  fit <- glm_gp(sce, design =  ~ 1,
                size_factors = size_factors, verbose = TRUE)
})

print("Fit glm_gp overdispersion = 0")

system.time({
  fit_disp0 <- glm_gp(sce, design =  ~ 1,
                overdispersion = 0,
                size_factors = size_factors, verbose = TRUE)
})




print("Fit glm_gp on disk")

system.time({
  fit_disk <- glm_gp(sce, design =  ~ 1,
                size_factors = size_factors, verbose = TRUE,
                on_disk = TRUE)
})


# 250 GB of RAM are not enough...
# print("Fit DESeq2")
# 
# system.time({
#   dds <- DESeq2::DESeqDataSet(sce, design =  ~ 1)
#   sizeFactors(dds) <- size_factors
#   dds <- DESeq2::estimateDispersions(dds, minmu = 1e-6)
#   dds <- DESeq2::nbinomWaldTest(dds, minmu = 1e-6)
# })

print("Fit edgeR")
system.time({
  model_matrix <- fit$model_matrix
  edgeR_data <- edgeR::DGEList(assay(sce))
  edgeR_data$samples$norm.factors <- size_factors
  edgeR_data <- edgeR::estimateDisp(edgeR_data, design = model_matrix)
  edgeR_fit <- edgeR::glmQLFit(edgeR_data, design = model_matrix)
})





print("Don't Store fits")

saveRDS(list(fit_disp0 = fit_disp0, fit = fit, fit = fit_disk, edgeR_fit = edgeR_fit,
             # dds_fit = dds, 
             dds_gp_fit = dds_gp),
        paste0("run_benchmark/tmp/benchmark_objects_full_", extension, ".RDS"))

# remove(fit, res, edgeR_fit, edgeR_data, res_edgeR, dds, res_dds, dds_gp, res_dds_gp)

# print("Run benchmark")
# 
# bench_res <- bench::mark(
#   glmGamPoi_memory = {
#     fit <- glm_gp(sce, design =  ~ pool + stage + tomato,
#                 size_factors = size_factors, verbose = TRUE)
#     res <- test_de(fit, reduced_design = ~ 1)
#   }, glmGamPoi_disk = {
#     fit <- glm_gp(sce, design =  ~ pool + stage + tomato,
#                   size_factors = size_factors, on_disk = TRUE, verbose = TRUE)
#     res <- test_de(fit, reduced_design = ~ 1)
#   }, edgeR = {
#     model_matrix <- fit$model_matrix
#     edgeR_data <- edgeR::DGEList(assay(sce))
#     edgeR_data$samples$norm.factors <- size_factors
#     edgeR_data <- edgeR::estimateDisp(edgeR_data, design = model_matrix)
#     edgeR_fit <- edgeR::glmQLFit(edgeR_data, design = model_matrix)
#     edgeR_test <- edgeR::glmQLFTest(edgeR_fit, coef = 2:4)
#     res_edgeR <- edgeR::topTags(edgeR_test)
#   }, DESeq2 = {
#     dds <- DESeq2::DESeqDataSet(sce, design =  ~ pool + stage + tomato)
#     sizeFactors(dds) <- size_factors
#     dds <- DESeq2::estimateDispersions(dds, minmu = 1e-6)
#     dds <- DESeq2::nbinomLRT(dds, reduced = ~ 1, minmu = 1e-6)
#     res_dds <- results(dds)
#   }, DESeq2_glmGamPoi = {
#     dds_gp <- DESeq2::DESeqDataSet(sce, design =  ~ pool + stage + tomato)
#     sizeFactors(dds_gp) <- size_factors
#     dds_gp <- DESeq2::estimateDispersions(dds_gp, minmu = 1e-6, fitType = "glmGamPoi")
#     dds_gp <- DESeq2::nbinomLRT(dds_gp, reduced = ~ 1, minmu = 1e-6, type = "glmGamPoi")
#     res_dds_gp <- results(dds_gp)
#   }, check = FALSE, min_iterations = 3)
# 
# print(bench_res)
# 
# saveRDS(bench_res, "run_benchmark/tmp/benchmark_result_full_4.RDS")



# --------------
sessionInfo()