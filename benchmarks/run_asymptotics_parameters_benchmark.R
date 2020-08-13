


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

# # To speed things up
sce <- sce[sample(seq_len(nrow(sce)), size = 5000), ]

# Make sure matrix is in memory
assay(sce) <- as.matrix(assay(sce))

# Calculate size factors

size_factors <- scran::calculateSumFactors(assay(sce))



# Study the complexity wrt the model matrix complexity

helper_df <- data.frame(
  cont1 = rnorm(n = ncol(sce)),
  cont2 = rnorm(n = ncol(sce)),
  cont3 = rnorm(n = ncol(sce)),
  cont4 = rnorm(n = ncol(sce)),
  cont5 = rnorm(n = ncol(sce)),
  cont6 = rnorm(n = ncol(sce)),
  cont7 = rnorm(n = ncol(sce)),
  cont8 = rnorm(n = ncol(sce)),
  cont9 = rnorm(n = ncol(sce)),
  cont10 = rnorm(n = ncol(sce)),
  cont11 = rnorm(n = ncol(sce)),
  cont12 = rnorm(n = ncol(sce)),
  cont13 = rnorm(n = ncol(sce)),
  cont14 = rnorm(n = ncol(sce)),
  cont15 = rnorm(n = ncol(sce)),
  cont16 = rnorm(n = ncol(sce)),
  cont17 = rnorm(n = ncol(sce)),
  cont18 = rnorm(n = ncol(sce)),
  cont19 = rnorm(n = ncol(sce)),
  cont20 = rnorm(n = ncol(sce))
)






working_model_matrix_cont <- model.matrix(~ cont1 + cont2 + cont3 + cont4 + cont5 + cont6 + cont7 + cont8 + cont9 + cont10 +
                                            cont11 + cont12 + cont13 + cont14 + cont15 + cont16 + cont17 + cont18 + cont19 + cont20 , data = helper_df)
print("Starting cont model subset")

for(model_matrix_subsetter in seq_len(21)){
  
  print(paste0("Model matrix subsetter ", model_matrix_subsetter))
  
  working_model_matrix <- working_model_matrix_cont[,seq_len(model_matrix_subsetter),drop=FALSE]
  
  # print("Fit glm_gp overdispersion = 0")
  
  # print(system.time({
  #   fit_disp0 <- glm_gp(sce, design =  working_model_matrix,
  #                       overdispersion = 0,
  #                       size_factors = size_factors, verbose = TRUE)
  # }))
  
  print("Fit glm_gp")
  
  print(system.time({
    fit <- glm_gp(sce, design =  working_model_matrix,
                  size_factors = size_factors, verbose = TRUE)
  }))
  
  # print("Fit glm_gp on disk")
  # 
  # print(system.time({
  #   fit_ondisk <- glm_gp(sce, design =  working_model_matrix,
  #                        size_factors = size_factors, verbose = TRUE,
  #                        on_disk = TRUE)
  # }))
  # 
  # print("Fit DESeq2 with glmGamPoi")
  # 
  # print(system.time({
  #   dds_gp <- DESeq2::DESeqDataSet(sce, design =  working_model_matrix)
  #   sizeFactors(dds_gp) <- size_factors
  #   dds_gp <- DESeq2::estimateDispersions(dds_gp, minmu = 1e-6, fitType = "glmGamPoi")
  #   dds_gp <- DESeq2::nbinomWaldTest(dds_gp, minmu = 1e-6, modelMatrix = working_model_matrix)
  # }))
  # 
  # print("Fit edgeR")
  # 
  # print(system.time({
  #   model_matrix <- working_model_matrix
  #   edgeR_data <- edgeR::DGEList(assay(sce), lib.size = size_factors)
  #   edgeR_data <- edgeR::estimateDisp(edgeR_data, design = model_matrix)
  #   edgeR_fit <- edgeR::glmQLFit(edgeR_data, design = model_matrix)
  # }))
  # 
  # print("Fit DESeq2")
  # 
  # print(system.time({
  #   dds <- DESeq2::DESeqDataSet(sce, design =  working_model_matrix)
  #   sizeFactors(dds) <- size_factors
  #   dds <- DESeq2::estimateDispersions(dds, minmu = 1e-6)
  #   dds <- DESeq2::nbinomWaldTest(dds, minmu = 1e-6, modelMatrix = working_model_matrix)
  # }))
}


print("Finished")


# --------------
sessionInfo()

