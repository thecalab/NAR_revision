library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(dplyr)

# Loading count matrix
read.counts <- read.table("featurecount_final.txt", header = TRUE, row.names = 1)
read.counts <- read.counts[,-c(1:5)]	
read.counts <- read.counts[,-c(15:16)] # just comparing 301 and 741
names(read.counts) <- c("TypeII_PD100_1","TypeII_PD100_2", "TypeII_PD400_1","TypeII_PD400_2",
                        "TypeII_PD800_1","TypeII_PD800_2","TypeI_PD100_1","TypeI_PD100_2",
                        "TypeI_PD350_1","TypeI_PD350_2", "TypeI_PD450_1","TypeI_PD450_2",
                        "TypeI_PD800_1","TypeI_PD800_2","TypeII_PD450_1","TypeII_PD450_2")
#To change column order 
read.counts <- read.counts[c(1:4,15,16,5:14)]
str(read.counts)

#Prepare metadata data frame
sample_info <- read.table("meta_all.txt", header = TRUE, as.is = FALSE, row.names = 1)
View(sample_info)
str(sample_info)

# Match the metadata and counts data
colnames(read.counts) == rownames(sample_info)

# Loading DESeq2 library
library(DESeq2)

# Create DESEq2 object
DESeq.ds <- DESeqDataSetFromMatrix(countData = read.counts,
                                   colData = sample_info,
                                   design = ~ type)
DESeq.ds
colData(DESeq.ds)
head(assay(DESeq.ds, "counts"))
rowData(DESeq.ds)

# Transform normalized counts using the rlog transformation
rld <- rlog(DESeq.ds, blind=TRUE)
## The blind=TRUE argument results in a transformation unbiased to sample condition information

# Hierarchical Clustering
# Extract the rlog matrix from the object
rld_mat <- assay(rld)
# Compute pairwise correlation values
rld_cor <- cor(rld_mat)
head(rld_cor)

## Plot heatmap for Figure 3B
library(pheatmap)
library(RColorBrewer)
pheatmap(rld_cor)
pheatmap(rld_cor, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20)
heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20)

# Principal components analysis (PCA) for Figure 3A
plotPCA(rld,intgroup="type")

# PCA customization (Lanscape, 4 x 6 inch)
pcaData <- plotPCA(rld,intgroup="type", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
library(ggplot2)
ggplot(pcaData, aes(PC1, PC2, color=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  xlim(-20,30) + ylim(-25,20)

#Differential analyses start

# The standard differential expression analysis steps are wrapped into a single function, DESeq
table(rowSums(counts(DESeq.ds))> 0)
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds))> 0, ]
DESeq.ds <- DESeq(DESeq.ds)
plotDispEsts(DESeq.ds)

# Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values.
str(colData(DESeq.ds)$condition)
library(dplyr)
library(tibble)
res <- results(DESeq.ds, contrast=c("condition","3sene","3pre")) %>% 
  data.frame() %>% 
  rownames_to_column(var="gene") 


#Log fold change shrinkage for visualization and ranking
library(apeglm)
resultsNames(DESeq.ds)
resLFC <- lfcShrink(DESeq.ds, coef="condition_3sene_vs_3pre", type="apeglm")
resLFC
summary(resLFC)
sum(resLFC$padj < 0.1, na.rm=TRUE)
mcols(resLFC)$description

# MA-plot
plotMA(resLFC, ylim=c(-2,2))

# Examining single gene
plotCounts(DESeq.ds, gene=which.min(resLFC$padj), intgroup="condition")

resOrdered <- resLFC[order(resLFC$log2FoldChange),]
DEG_results <- as.data.frame(resOrdered)

#From now on, Pathway enrichment anaylsis

library("AnnotationDbi")
library("EnsDb.Mmusculus.v79")
keytypes(EnsDb.Mmusculus.v79)

# Making gene lists
ens.str <- substr(rownames(DEG_results), 1, 18)
head(ens.str)

resOrdered$symbol <- mapIds(EnsDb.Mmusculus.v79,
                            keys=ens.str,
                            column="SYMBOL",
                            keytype="GENEID",
                            multiVals="first")

resOrdered$biotype <- mapIds(EnsDb.Mmusculus.v79,
                             keys=ens.str,
                             column="GENEBIOTYPE",
                             keytype="GENEID",
                             multiVals="first")

length(which(is.na(resOrdered$symbol)))
length(which(is.na(resOrdered$biotype)))

library(tidyverse)
a <- as.data.frame(resOrdered) %>%
  rownames_to_column(var="gene") 
write.csv(a, file="a.csv")

# DESeq2 creates a matrix when you use the counts() function
## First convert normalized_counts to a data frame and transfer the row names to a new column called "gene"

normalized_counts <- counts(DESeq.ds, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene") 

# Making annotation matrix
## Return the Ensembl IDs for a set of genes
annotations_orgDb <- AnnotationDbi::select(EnsDb.Mmusculus.v79, # database
                                           keys = ens.str,  # data to use for retrieval
                                           columns = c("SYMBOL","ENTREZID","GENEBIOTYPE"), # information to retreive for given data
                                           keytype = "GENEID") # type of data given in 'keys' argument

## Check for NA entries
length(which(is.na(annotations_orgDb$SYMBOL)))
## Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_orgDb$SYMBOL) == FALSE)
## Return only the non-duplicated genes using indices
annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]

# Merging nomralized_counts with annotation matrix
normalized_counts <- merge(normalized_counts, annotations_orgDb, by.x="gene", by.y="GENEID")

# Merging normalized_counts with DEG results
DEG_final <- merge(normalized_counts, a, by.x="gene", by.y="gene")
write.csv(DEG_final, file="DEG_final.csv")


#DEG list clustering with all genes for Supplementary Figure S3F

rownames(DEG_final) <- DEG_final$SYMBOL

DEG_final.1 <-DEG_final[c(2:17)]
DEG_final.1.matrix <- as.matrix(DEG_final.1)
clustRows <- hclust(as.dist(1-cor(t(DEG_final.1.matrix), method="pearson")), method="complete") 
clustColumns <- hclust(as.dist(1-cor(DEG_final.1.matrix, method="spearman")), method="complete") #cluster columns by spearman correlation
module.assign <- cutree(clustRows, k=8)

module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
myheatcolors <- colorRampPalette(colors=c("blue","white","red"))(100)

# Produce a static heatmap of DEGs 
#plot the hclust results as a heatmap
dev.off()
heatmap.2(DEG_final.1.matrix, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20)) 


#rerun the heatmap script above using diffData.AVG as input instead of diffData
colnames(DEG_final.1.matrix) <- meta$type

#now an old function from the limma package to average your replicates 
diffSets.AVG <- avearrays(DEG_final.1.matrix)

clustRows <- hclust(as.dist(1-cor(t(diffSets.AVG), method="pearson")), method="complete") 
clustColumns <- hclust(as.dist(1-cor(diffSets.AVG, method="spearman")), method="complete") #cluster columns by spearman correlation
module.assign <- cutree(clustRows, k=8)

module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 

dev.off()
heatmap.2(diffSets.AVG, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,4)) 

#Cluster analysis for Figure 3C and D
library(DEGreport)
cluster_rlog <- rld_mat
clusters <- degPatterns(cluster_rlog, metadata = colData(rld), time = "time", col=NULL)

a <- degPlotCluster(clusters$normalized, time="time", points=FALSE, lines=FALSE, color="colored")
library(ggplot2)


#to find GO terms of each cluster
cluster_df <- clusters$df

for (i in c(1:8)){
  eval(parse(text=paste("group",i,"<-filter(cluster_df,cluster==",i,")",sep="")))}

#GSEA analysis
library("AnnotationDbi")
library("EnsDb.Mmusculus.v79")
keytypes(EnsDb.Mmusculus.v79)
library(clusterProfiler)
library(org.Mm.eg.db)

## Run GO enrichment analysis for Supplementary Figure S4

for (i in c(1:8)){
  eval(parse(text=paste("ego",i,"<- enrichGO(gene = group",i,
                        "$genes,keyType = \"ENSEMBL\",
                      OrgDb = org.Mm.eg.db,ont = \"BP\",
                      pAdjustMethod = \"BH\",qvalueCutoff = 0.05,
                      readable = TRUE)",sep="")))
}

dotplot(ego1,  orderBy="GeneRatio", showCategory=20)

#With newly acquired ALT survivors, overall analses were same. We used Combat-seq for batch correction:
#Defining Combat-seq---------
#' Adjust for batch effects using an empirical Bayes framework in RNA-seq raw counts
#' 
#' ComBat_seq is an extension to the ComBat method using Negative Binomial model.
#' 
#' @param counts Raw count matrix from genomic studies (dimensions gene x sample) 
#' @param batch Batch covariate (only one batch allowed)
#' @param group Vector / factor for condition of interest 
#' @param covar_mod Model matrix for other covariates to include in linear model besides batch and condition of interest
#' @param full_mod Boolean, if TRUE include condition of interest in model
#' @param shrink Boolean, whether to apply empirical Bayes estimation on parameters
#' @param shrink.disp Boolean, whether to apply empirical Bayes estimation on dispersion
#' @param gene.subset.n Number of genes to use in empirical Bayes estimation, only useful when shrink = TRUE
#' @return data A probe x sample count matrix, adjusted for batch effects.
#' 
#' @examples 
#' 
#' 
#' 
#' @export
#' 

ComBat_seq <- function(counts, batch, group=NULL, covar_mod=NULL, full_mod=TRUE, 
                       shrink=FALSE, shrink.disp=FALSE, gene.subset.n=NULL){  
  ########  Preparation  ######## 
  counts <- as.matrix(counts)
  
  ## Does not support 1 sample per batch yet
  batch <- as.factor(batch)
  if(any(table(batch)<=1)){
    stop("ComBat-seq doesn't support 1 sample per batch yet")
  }
  
  ## Remove genes with only 0 counts in any batch
  keep_lst <- lapply(levels(batch), function(b){
    which(apply(counts[, batch==b], 1, function(x){!all(x==0)}))
  })
  keep <- Reduce(intersect, keep_lst)
  rm <- setdiff(1:nrow(counts), keep)
  countsOri <- counts
  counts <- counts[keep, ]
  
  # require bioconductor 3.7, edgeR 3.22.1
  dge_obj <- DGEList(counts=counts)
  
  ## Prepare characteristics on batches
  n_batch <- nlevels(batch)  # number of batches
  batches_ind <- lapply(1:n_batch, function(i){which(batch==levels(batch)[i])}) # list of samples in each batch  
  n_batches <- sapply(batches_ind, length)
  #if(any(n_batches==1)){mean_only=TRUE; cat("Note: one batch has only one sample, setting mean.only=TRUE\n")}
  n_sample <- sum(n_batches)
  cat("Found",n_batch,'batches\n')
  
  ## Make design matrix 
  # batch
  batchmod <- model.matrix(~-1+batch)  # colnames: levels(batch)
  # covariate
  group <- as.factor(group)
  if(full_mod & nlevels(group)>1){
    cat("Using full model in ComBat-seq.\n")
    mod <- model.matrix(~group)
  }else{
    cat("Using null model in ComBat-seq.\n")
    mod <- model.matrix(~1, data=as.data.frame(t(counts)))
  }
  # drop intercept in covariate model
  if(!is.null(covar_mod)){
    if(is.data.frame(covar_mod)){
      covar_mod <- do.call(cbind, lapply(1:ncol(covar_mod), function(i){model.matrix(~covar_mod[,i])}))
    }
    covar_mod <- covar_mod[, !apply(covar_mod, 2, function(x){all(x==1)})]
  }
  # bind with biological condition of interest
  mod <- cbind(mod, covar_mod)
  # combine
  design <- cbind(batchmod, mod)
  
  ## Check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  #if(!is.null(ref)){check[ref]=FALSE} ## except don't throw away the reference batch indicator
  design <- as.matrix(design[,!check])
  cat("Adjusting for",ncol(design)-ncol(batchmod),'covariate(s) or covariate level(s)\n')
  
  ## Check if the design is confounded
  if(qr(design)$rank<ncol(design)){
    #if(ncol(design)<=(n_batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n_batch+1)){stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat-Seq")}
    if(ncol(design)>(n_batch+1)){
      if((qr(design[,-c(1:n_batch)])$rank<ncol(design[,-c(1:n_batch)]))){stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      }else{stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat-Seq")}}
  }
  
  ## Check for missing values in count matrix
  NAs = any(is.na(counts))
  if(NAs){cat(c('Found',sum(is.na(counts)),'Missing Data Values\n'),sep=' ')}
  
  
  ########  Estimate gene-wise dispersions within each batch  ########
  cat("Estimating dispersions\n")
  ## Estimate common dispersion within each batch as an initial value
  disp_common <- sapply(1:n_batch, function(i){
    if((n_batches[i] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[i]], ])$rank < ncol(mod)){ 
      # not enough residual degree of freedom
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=NULL, subset=nrow(counts)))
    }else{
      return(estimateGLMCommonDisp(counts[, batches_ind[[i]]], design=mod[batches_ind[[i]], ], subset=nrow(counts)))
    }
  })
  
  ## Estimate gene-wise dispersion within each batch 
  genewise_disp_lst <- lapply(1:n_batch, function(j){
    if((n_batches[j] <= ncol(design)-ncol(batchmod)+1) | qr(mod[batches_ind[[j]], ])$rank < ncol(mod)){
      # not enough residual degrees of freedom - use the common dispersion
      return(rep(disp_common[j], nrow(counts)))
    }else{
      return(estimateGLMTagwiseDisp(counts[, batches_ind[[j]]], design=mod[batches_ind[[j]], ], 
                                    dispersion=disp_common[j], prior.df=0))
    }
  })
  names(genewise_disp_lst) <- paste0('batch', levels(batch))
  
  ## construct dispersion matrix
  phi_matrix <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(k in 1:n_batch){
    phi_matrix[, batches_ind[[k]]] <- vec2mat(genewise_disp_lst[[k]], n_batches[k]) 
  }
  
  
  ########  Estimate parameters from NB GLM  ########
  cat("Fitting the GLM model\n")
  glm_f <- glmFit(dge_obj, design=design, dispersion=phi_matrix, prior.count=1e-4) #no intercept - nonEstimable; compute offset (library sizes) within function
  alpha_g <- glm_f$coefficients[, 1:n_batch] %*% as.matrix(n_batches/n_sample) #compute intercept as batch-size-weighted average from batches
  new_offset <- t(vec2mat(getOffset(dge_obj), nrow(counts))) +   # original offset - sample (library) size
    vec2mat(alpha_g, ncol(counts))  # new offset - gene background expression # getOffset(dge_obj) is the same as log(dge_obj$samples$lib.size)
  glm_f2 <- glmFit.default(dge_obj$counts, design=design, dispersion=phi_matrix, offset=new_offset, prior.count=1e-4) 
  
  gamma_hat <- glm_f2$coefficients[, 1:n_batch]
  mu_hat <- glm_f2$fitted.values
  phi_hat <- do.call(cbind, genewise_disp_lst)
  
  
  ########  In each batch, compute posterior estimation through Monte-Carlo integration  ########  
  if(shrink){
    cat("Apply shrinkage - computing posterior estimates for parameters\n")
    mcint_fun <- monte_carlo_int_NB
    monte_carlo_res <- lapply(1:n_batch, function(ii){
      if(ii==1){
        mcres <- mcint_fun(dat=counts[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]], 
                           gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)
      }else{
        invisible(capture.output(mcres <- mcint_fun(dat=counts[, batches_ind[[ii]]], mu=mu_hat[, batches_ind[[ii]]], 
                                                    gamma=gamma_hat[, ii], phi=phi_hat[, ii], gene.subset.n=gene.subset.n)))
      }
      return(mcres)
    })
    names(monte_carlo_res) <- paste0('batch', levels(batch))
    
    gamma_star_mat <- lapply(monte_carlo_res, function(res){res$gamma_star})
    gamma_star_mat <- do.call(cbind, gamma_star_mat)
    phi_star_mat <- lapply(monte_carlo_res, function(res){res$phi_star})
    phi_star_mat <- do.call(cbind, phi_star_mat)
    
    if(!shrink.disp){
      cat("Apply shrinkage to mean only\n")
      phi_star_mat <- phi_hat
    }
  }else{
    cat("Shrinkage off - using GLM estimates for parameters\n")
    gamma_star_mat <- gamma_hat
    phi_star_mat <- phi_hat
  }
  
  
  ########  Obtain adjusted batch-free distribution  ########
  mu_star <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(jj in 1:n_batch){
    mu_star[, batches_ind[[jj]]] <- exp(log(mu_hat[, batches_ind[[jj]]])-vec2mat(gamma_star_mat[, jj], n_batches[jj]))
  }
  phi_star <- rowMeans(phi_star_mat)
  
  
  ########  Adjust the data  ########  
  cat("Adjusting the data\n")
  adjust_counts <- matrix(NA, nrow=nrow(counts), ncol=ncol(counts))
  for(kk in 1:n_batch){
    counts_sub <- counts[, batches_ind[[kk]]]
    old_mu <- mu_hat[, batches_ind[[kk]]]
    old_phi <- phi_hat[, kk]
    new_mu <- mu_star[, batches_ind[[kk]]]
    new_phi <- phi_star
    adjust_counts[, batches_ind[[kk]]] <- match_quantiles(counts_sub=counts_sub, 
                                                          old_mu=old_mu, old_phi=old_phi, 
                                                          new_mu=new_mu, new_phi=new_phi)
  }
  
  #dimnames(adjust_counts) <- dimnames(counts)
  #return(adjust_counts)
  
  ## Add back genes with only 0 counts in any batch (so that dimensions won't change)
  adjust_counts_whole <- matrix(NA, nrow=nrow(countsOri), ncol=ncol(countsOri))
  dimnames(adjust_counts_whole) <- dimnames(countsOri)
  adjust_counts_whole[keep, ] <- adjust_counts
  adjust_counts_whole[rm, ] <- countsOri[rm, ]
  return(adjust_counts_whole)
}
#Defining vec2mat and other functions------
####  Expand a vector into matrix (columns as the original vector)
vec2mat <- function(vec, n_times){
  return(matrix(rep(vec, n_times), ncol=n_times, byrow=FALSE))
}


####  Monte Carlo integration functions
monte_carlo_int_NB <- function(dat, mu, gamma, phi, gene.subset.n){
  weights <- pos_res <- list()
  for(i in 1:nrow(dat)){
    m <- mu[-i,!is.na(dat[i,])]
    x <- dat[i,!is.na(dat[i,])]
    gamma_sub <- gamma[-i]
    phi_sub <- phi[-i]
    
    # take a subset of genes to do integration - save time
    if(!is.null(gene.subset.n) & is.numeric(gene.subset.n) & length(gene.subset.n)==1){
      if(i==1){cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))}
      mcint_ind <- sample(1:(nrow(dat)-1), gene.subset.n, replace=FALSE)
      m <- m[mcint_ind, ]; gamma_sub <- gamma_sub[mcint_ind]; phi_sub <- phi_sub[mcint_ind]
      G_sub <- gene.subset.n
    }else{
      if(i==1){cat("Using all genes for Monte Carlo integration; the function runs very slow for large number of genes\n")}
      G_sub <- nrow(dat)-1
    }
    
    #LH <- sapply(1:G_sub, function(j){sum(log2(dnbinom(x, mu=m[j,], size=1/phi_sub[j])+1))})  
    LH <- sapply(1:G_sub, function(j){prod(dnbinom(x, mu=m[j,], size=1/phi_sub[j]))})
    LH[is.nan(LH)]=0; 
    if(sum(LH)==0 | is.na(sum(LH))){
      pos_res[[i]] <- c(gamma.star=as.numeric(gamma[i]), phi.star=as.numeric(phi[i]))
    }else{
      pos_res[[i]] <- c(gamma.star=sum(gamma_sub*LH)/sum(LH), phi.star=sum(phi_sub*LH)/sum(LH))
    }
    
    weights[[i]] <- as.matrix(LH/sum(LH))
  }
  pos_res <- do.call(rbind, pos_res)
  weights <- do.call(cbind, weights)
  res <- list(gamma_star=pos_res[, "gamma.star"], phi_star=pos_res[, "phi.star"], weights=weights)	
  return(res)
} 


####  Match quantiles
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi){
  new_counts_sub <- matrix(NA, nrow=nrow(counts_sub), ncol=ncol(counts_sub))
  for(a in 1:nrow(counts_sub)){
    for(b in 1:ncol(counts_sub)){
      if(counts_sub[a, b] <= 1){
        new_counts_sub[a,b] <- counts_sub[a, b]
      }else{
        tmp_p <- pnbinom(counts_sub[a, b]-1, mu=old_mu[a, b], size=1/old_phi[a])
        if(abs(tmp_p-1)<1e-4){
          new_counts_sub[a,b] <- counts_sub[a, b]  
          # for outlier count, if p==1, will return Inf values -> use original count instead
        }else{
          new_counts_sub[a,b] <- 1+qnbinom(tmp_p, mu=new_mu[a, b], size=1/new_phi[a])
        }
      }
    }
  }
  return(new_counts_sub)
}



mapDisp <- function(old_mu, new_mu, old_phi, divider){
  new_phi <- matrix(NA, nrow=nrow(old_mu), ncol=ncol(old_mu))
  for(a in 1:nrow(old_mu)){
    for(b in 1:ncol(old_mu)){
      old_var <- old_mu[a, b] + old_mu[a, b]^2 * old_phi[a, b]
      new_var <- old_var / (divider[a, b]^2)
      new_phi[a, b] <- (new_var - new_mu[a, b]) / (new_mu[a, b]^2)
    }
  }
  return(new_phi)
}