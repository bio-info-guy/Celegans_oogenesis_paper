suppressPackageStartupMessages({
  library(grImport2)
  library(GenomicFeatures)
  library(stringr)
  library(biomaRt)
  library(BASiCS)
  library(umap)
  library(apeglm)
  library(heatmap3)
  library(ggplot2)
  library(ggfortify)
  library(stringr)
  library(RColorBrewer)
  library(MKmisc)
  library(DESeq2)
  library(Rtsne)
  library(MAST)
  library(reticulate)
  library(ComplexHeatmap)
  library(edgeR)
  library(GGally)
  library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(stringr)
  library(NMF)
  library(rsvd)
  library(RColorBrewer)
  library(MAST)
  library(pcaMethods)
  library(segmented)
  library(robust)
  library(MASS)
  library(umap)
  library(QoRTs)
  library(tximport)
  library(Seurat)
  require(ReactomePA)
  require(clusterProfiler)
  require(org.Ce.eg.db)
  require(meshes)
  library(msigdbr)
  library(Seurat)
  library(doParallel)
  library(grid)
  library(gridExtra)
  library(pathview)
  library(igraph)
  library(ggraph)
  library(WGCNA)
  library(flashClust)
  library(ggrepel)
  library(monocle)
})

DOT_COLOR <- c("#9ce6fd","#4e76f8","#293acc","#4d0294",
"#a8fd1f","#ffd966","#ffd966","#B20000","#744700",
"#008331","#8218f2","#8cb4f9", "#a08cf9","#f48cf9","#000000",
"#ea9999","#e69138","#b4a7d6","#c5156c"
)

names(DOT_COLOR) <- c('S1', 'S2', 'S3', 'S4', 
'-4','-3','-2-3','-2','-1',
'P0','-2S4','S2S3','-3S4','-2-3S4','Unknown',
"2-cell","4-cell","8-cell","16-cell"

)

STAGE_ORDER <-  c('S1', 'S2', 'S3', 'S4','-3', '-2', '-1', 'P0')

OLD_LABEL_2_NEW <- c('S1'='S1', 'S2'='S2', 'S3'='S3', 'S4'='S4', 
                     'F3'='-3','F2'='-2','F1'='-1', 'P0'='P0', 
                     'F2F3S4'='-2-3S4', 'F2S4'='-2S4', 'F3S4'='-3S4',
                     'F4'='-4')

read_expression <- function(dir, mode = 'salmon', tx2gene=NULL){
  if(mode %in% c('star', 'hisat')){
    files <- list.files(dir, pattern = paste(mode, 'htseq.ct', sep = '.'), full.names = T)
    sample_names <- as.character(data.frame(strsplit(list.files(dir, pattern = paste(mode, 'htseq.ct', sep = '.')), split = '_'))[1,])
    cts <- do.call(cbind, lapply(files, FUN = function(x){read.table(x, row.names = 1, header = F)}))
    colnames(cts) <- sample_names
    cts <- cts[1:(nrow(cts)-5),]
    #cts <- cts[rowSums(cts != 0) > 0,]
    return(cts)
  }else if(mode  == 'salmon'){
    tx2gene <- read.csv(tx2gene,sep='\t' ,header = T, col.names = c("TXNAME", "GENEID"))
    salmon_files <- list.files(dir, recursive = T, pattern='*quant.sf', full.names = T)
    names(salmon_files) <- as.character(as.data.frame(strsplit(salmon_files, '/'))[length(strsplit(salmon_files, '/')[[1]])-1,])
    gene.salmon <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene, dropInfReps = T, importer=read.delim)
    gene.salmon$abundance <- Matrix(gene.salmon$abundance, sparse = T)
    gene.salmon$counts <- Matrix(gene.salmon$counts, sparse = T)
    #gene.salmon$counts <- gene.salmon$counts[rowSums(gene.salmon$counts != 0) > 0,]
    #gene.salmon$abundance <- gene.salmon$abundance[rowSums(gene.salmon$abundance != 0) > 0,]
    transcript.salmon <- tximport(salmon_files, type = "salmon", txOut = TRUE, dropInfReps = T, importer=read.delim)
    transcript.salmon$abundance <- Matrix(transcript.salmon$abundance, sparse = T)
    transcript.salmon$counts <- Matrix(transcript.salmon$counts, sparse = T)
    #transcript.salmon$counts <- transcript.salmon$counts[rowSums(transcript.salmon$counts != 0) > 0,]
    #transcript.salmon$abundance <- transcript.salmon$abundance[rowSums(transcript.salmon$abundance != 0) > 0,]
    return(list(gene = gene.salmon, transcript = transcript.salmon))
  }else if(mode == 'rsem'){
    file_names <- list.files(dir, recursive = F, pattern='*genes.results', full.names = F)
    gene_files <- list.files(dir, recursive = F, pattern='*genes.results', full.names = T)
    names(gene_files) <- as.character(as.data.frame(strsplit(file_names, '\\.'))[1,])
    print(length(gene_files))
    gene.rsem <- tximport(gene_files, type = "rsem", txIn = FALSE, txOut = FALSE, importer=read.delim)
    gene.rsem$abundance <- Matrix(gene.rsem$abundance, sparse = T)
    gene.rsem$counts <- Matrix(gene.rsem$counts, sparse = T)
    #gene.rsem$counts <- gene.rsem$counts[rowSums(gene.rsem$counts != 0) > 0,]
    #gene.rsem$abundance <- gene.rsem$abundance[rowSums(gene.rsem$abundance != 0) > 0,]
    file_names <- list.files(dir, recursive = F, pattern='*isoforms.results', full.names = F)
    file_names <- file_names[!duplicated(file_names)]
    transcript_files <- list.files(dir, recursive = F, pattern='*isoforms.results', full.names = T)
    transcript_files <- transcript_files[!duplicated(transcript_files)]
    names(transcript_files) <- as.character(as.data.frame(strsplit(file_names, '\\.'))[1,])
    print(length(transcript_files))
    transcript.rsem <- tximport(transcript_files, type = "rsem", txIn = TRUE, txOut = TRUE, importer=read.delim)
    transcript.rsem$abundance <- Matrix(transcript.rsem$abundance, sparse = T)
    transcript.rsem$counts <- Matrix(transcript.rsem$counts, sparse = T)
    #transcript.rsem$counts <- transcript.rsem$counts[rowSums(transcript.rsem$counts != 0) > 0,]
    #transcript.rsem$abundance <- transcript.rsem$abundance[rowSums(transcript.rsem$abundance != 0) > 0,]
    return(list(gene = gene.rsem, transcript = transcript.rsem))
  }
  else{
    cat('only support salmon, rsem, hisat (htseq) and star (htseq)\n')
  }
}
prepareMonocle <- function(exprs, meta, features,tpm2abs=T, ct2tpm=F, monocle3 = F)
{
  ct <- exprs
  genes <- intersect(row.names(ct)[which(rowSums(ct) > 0)], row.names(features))
  ct <- ct[genes,]
  tpm <- ct
  features <- features[genes,]
  
  #colnames(features) <- c('gene_id','gene_short_name','biotype')
  #worms.gene.id <- data.frame(worms.gene.id, row.names = 3)
  fd <- new("AnnotatedDataFrame", data= features[row.names(ct),])
  pd <- new("AnnotatedDataFrame", data= meta)
  lower_detect = 5; min_expr = 1
  
  if(tpm2abs){
    if(ct2tpm){
      len <- obj$ct$raw$len[genes,]
      ct <- t(t(ct/(len/1e3))/(colSums(ct/(len/1e3))/1e6))
      cds <- newCellDataSet(as.matrix(ct), phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = tobit(Lower=0.1))
      ct  <- relative2abs(cds, method = "num_genes")
    }else{
      ct <- tpm
      cds <- newCellDataSet(as.matrix(ct), phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = tobit(Lower=0.1))
      ct  <- relative2abs(cds, method = "num_genes")
      #ct <- floor(as.matrix(ct))
    }
    lower_detect = 0.5
    min_expr = 0.1
    cds <- newCellDataSet(as.matrix(ct), phenoData = pd, featureData = fd, lowerDetectionLimit = lower_detect, expressionFamily = negbinomial.size())
    cds <- estimateSizeFactors(cds)
    cds <- tryCatch({
      cds <- estimateDispersions(cds)
    }, error=function(cond) {
      message(cond)
      # Choose a return value in case of error
      cat('running local\n')
      return(estimateDispersions(cds, fitType = 'pool'))
    })
  }else{
    if(ct2tpm){
      len <- obj$ct$raw$len[genes,]
      ct <- t(t(ct/(len/1e3))/(colSums(ct/(len/1e3))/1e6))
      cds <- newCellDataSet(as.matrix(ct), phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = tobit(Lower=0.1))
      cds <- estimateSizeFactors(cds)
    }else{
      cds <- newCellDataSet(as.matrix(ct), phenoData = pd, featureData = fd, lowerDetectionLimit = lower_detect, expressionFamily = negbinomial.size())
      cds <- estimateSizeFactors(cds)
      cds <- estimateDispersions(cds)
    }
    
  }
  #rpc_matrix <- relative2abs(cds, method = "num_genes")
  #cds <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"), phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
  
  cds <- detectGenes(cds, min_expr = min_expr)
  cds <- cds[which(!grepl('msp-', cds@featureData@data$gene_short_name) & !grepl('^spe-', cds@featureData@data$gene_short_name)),]
  if(monocle3)
  {
    cds <- new_cell_data_set(expression_data = cds@assayData$exprs,
                             cell_metadata = cds@phenoData@data,
                             gene_metadata = cds@featureData@data)
    return(cds)
  }
  return(cds)
}
  

prepareMonocle <- function(exprs, meta, features,tpm2abs=T, ct2tpm=F, monocle3 = F)
{
  ct <- exprs
  genes <- intersect(row.names(ct)[which(rowSums(ct) > 0)], row.names(features))
  ct <- ct[genes,]
  tpm <- ct
  features <- features[genes,]
  
  #colnames(features) <- c('gene_id','gene_short_name','biotype')
  #worms.gene.id <- data.frame(worms.gene.id, row.names = 3)
  fd <- new("AnnotatedDataFrame", data= features[row.names(ct),])
  pd <- new("AnnotatedDataFrame", data= meta)
  lower_detect = 5; min_expr = 1
  
  if(tpm2abs){
    if(ct2tpm){
      len <- obj$ct$raw$len[genes,]
      ct <- t(t(ct/(len/1e3))/(colSums(ct/(len/1e3))/1e6))
      cds <- newCellDataSet(as.matrix(ct), phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = tobit(Lower=0.1))
      ct  <- relative2abs(cds, method = "num_genes")
    }else{
      ct <- tpm
      cds <- newCellDataSet(as.matrix(ct), phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = tobit(Lower=0.1))
      ct  <- relative2abs(cds, method = "num_genes")
      #ct <- floor(as.matrix(ct))
    }
    lower_detect = 0.5
    min_expr = 0.1
    cds <- newCellDataSet(as.matrix(ct), phenoData = pd, featureData = fd, lowerDetectionLimit = lower_detect, expressionFamily = negbinomial.size())
    cds <- estimateSizeFactors(cds)
    cds <- tryCatch({
      cds <- estimateDispersions(cds)
    }, error=function(cond) {
      message(cond)
      # Choose a return value in case of error
      cat('running local\n')
      return(estimateDispersions(cds, fitType = 'pool'))
    })
  }else{
    if(ct2tpm){
      len <- obj$ct$raw$len[genes,]
      ct <- t(t(ct/(len/1e3))/(colSums(ct/(len/1e3))/1e6))
      cds <- newCellDataSet(as.matrix(ct), phenoData = pd, featureData = fd, lowerDetectionLimit = 0.1, expressionFamily = tobit(Lower=0.1))
      cds <- estimateSizeFactors(cds)
    }else{
      cds <- newCellDataSet(as.matrix(ct), phenoData = pd, featureData = fd, lowerDetectionLimit = lower_detect, expressionFamily = negbinomial.size())
      cds <- estimateSizeFactors(cds)
      cds <- estimateDispersions(cds)
    }
    
  }
  #rpc_matrix <- relative2abs(cds, method = "num_genes")
  #cds <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"), phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
  
  cds <- detectGenes(cds, min_expr = min_expr)
  cds <- cds[which(!grepl('msp-', cds@featureData@data$gene_short_name) & !grepl('^spe-', cds@featureData@data$gene_short_name)),]
  if(monocle3)
  {
    cds <- new_cell_data_set(expression_data = cds@assayData$exprs,
                             cell_metadata = cds@phenoData@data,
                             gene_metadata = cds@featureData@data)
    return(cds)
  }
  return(cds)
}


geneHeatmap <- function(cds, exprs_mat, tpm = F,
                        genes = NULL, clusterGenes = T,
                        dist_row = 'pearson', dist_col = 'euclidean',
                        main = "Gene Expression Heatmap",
                        fontsize = 5, legend = T, row_color = NULL,
                        annotations, clust_n = 30, mmethod = 'ward.D2') {
  

  stages <- c('S1', 'S2', 'S3', 'S4', 'F4', 'F3','F2','F2F3','F2S4','F3S4','F1', 'P0', 'Unknown')
  cluster_mat <- cds
  
  if(!is.null(genes)){mat = mat[genes,order(match(annotations$cellType, stages))]}
  DISTMETHOD_row = dist_row
  DISTMETHOD_col = dist_col
  dist_func_col = function(df){
    #df <- 2**df-1
    #print(dim(df))
    #mat <- DESeq2::varianceStabilizingTransformation(as.matrix(round(t(df))))
    if (DISTMETHOD_col %in% c("pearson", "spearman", "kandall", "cosine", "mcd", "ogk"))
    {
      
      #mat <- df
      require(MKmisc)
      
      cdist <- corDist(t(cluster_mat)[rownames(df),], method = DISTMETHOD_col)
      #cdist <- quickWGCNA_dist(cluster_mat[,row.names(df)])
      
      #print("calculating correlation distance")
    }
    else
    {
      cdist <- dist(t(cluster_mat)[rownames(df),], method = DISTMETHOD_col)
      #print("calculating non-correlation distance")
    }
    print(DISTMETHOD_col)
    print(dim(as.matrix(cdist)))
    cdist
  }
  dist_func_row = function(df){
    #df <- 2**df-1
    #mat <- DESeq2::varianceStabilizingTransformation(as.matrix(round(df)))
    if (DISTMETHOD_row %in% c("pearson", "spearman", "kandall", "cosine", "mcd", "ogk", "wgcna", 'bicor'))
    {
      #mat <- DESeq2::varianceStabilizingTransformation(as.matrix(df))
      #mat <- df
      require(MKmisc)
      if(DISTMETHOD_row == 'wgcna'){
        cdist <- quickWGCNA_dist(cluster_mat)
      }else{
        if(DISTMETHOD_row == 'bicor'){
          cdist <- as.dist(1-bicor(t(cluster_mat)))
        }else{
          cdist <- corDist(cluster_mat, method = DISTMETHOD_row)
        }
      }
      #cdist <- quickWGCNA_dist(cluster_mat)
      #print("calculating correlation distance")
    }
    else
    {
      cdist <- dist(cluster_mat, method = DISTMETHOD_row)
      #print("calculating non-correlation distance")
    }
    print(DISTMETHOD_row)
    print(dim(as.matrix(cdist)))
    cdist
  }
  #View(as.matrix(cdist))
  #print(color)
  anno_color <-
    list(
      cellType = DOT_COLOR[unique(as.character(annotations$cellType))],
      logfc = c("negative" = 'red', 'positive' = 'green')
      
    )
  if (nrow(cds) >= 0) {
    show_rown <- F
  } else{
    show_rown <- T
  }
  
  annotations = annotations[,'cellType', drop = F]
  if(!is.null(row_color)){
    #logfc <- c("FALSE" = 'negative', 'TRUE' = 'positive')[as.character(row_color > 0)]
    #names(logfc) <- names(row_color)
    logfc = data.frame(row.names = names(row_color), pval = row_color)
  }else{
    logfc = NULL
  }
  if(clust_n == 'auto'){
    clust_n = NbClust(data = cds, diss = cdist, distance = NULL, method = mmethod, index = 'all', max.nc = 50)$Best.partition
  }
  breaksList = seq(0, 13, by = 1)
  print(dim(exprs_mat))
  print(class(exprs_mat))
  #color_pal <- colorRampPalette(c('black','gold','white'))(12)
  color_pal <- colorRampPalette(c('blue','white','red'))(12)
  heat_obj <- Heatmap(matrix = exprs_mat, col = color_pal, name = 'Gene Heatmap',
                      na_col = "grey",
                      color_space = "LAB",
                      rect_gp = gpar(col = NA),
                      border = T,
                      border_gp = gpar(col = "black"),
                      cell_fun = NULL,
                      layer_fun = NULL,
                      row_title = character(0),
                      row_title_side = c("left"),
                      row_title_gp = gpar(fontsize = 8),
                      row_title_rot = 0,
                      column_title = character(0),
                      column_title_side = c("top"),
                      column_title_gp = gpar(fontsize = 13.2),
                      column_title_rot = 0,
                      
                      cluster_rows = TRUE,
                      cluster_row_slices = F,
                      clustering_distance_rows = dist_func_row,
                      clustering_method_rows = "ward.D2",
                      row_dend_side = c("left"),
                      row_dend_width = unit(10, "mm"),
                      show_row_dend = F,
                      row_dend_gp = gpar(),
                      
                      cluster_columns = T,
                      cluster_column_slices =  F,
                      clustering_distance_columns = dist_func_col,
                      clustering_method_columns = "ward.D2",
                      column_dend_side = c("top"),
                      column_dend_height = unit(10, "mm"),
                      show_column_dend = F,
                      column_dend_gp = gpar(),
                      
                      row_order = NULL,
                      column_order = NULL,
                      
                      row_names_side = c( "left"),
                      show_row_names = F,
                      row_names_max_width = unit(6, "cm"),
                      row_names_gp = gpar(fontsize = 8),
                      row_names_rot = 0,
                      row_names_centered = FALSE,
                      column_labels = colnames(cds),
                      column_names_side = c("bottom"),
                      show_column_names = F,
                      column_names_max_height = unit(8, "cm"),
                      column_names_gp = gpar(fontsize = 7),
                      column_names_rot = 270,
                      column_names_centered = FALSE,
                      
                      top_annotation = HeatmapAnnotation(cellType = annotations$cellType,
                                                         col = list(cellType = DOT_COLOR)),
                      bottom_annotation = HeatmapAnnotation( cellType= anno_text(colnames(cds), rot = 270, 
                                                                                 gp = gpar(fontsize = 4), location = 0.5, 
                                                                                 just = "center"), which = c('column')),
                      left_annotation = NULL,
                      right_annotation = NULL,
                      
                      km = 1,
                      split = NULL,
                      row_km = 1,
                      row_km_repeats = 1,
                      row_split = clust_n,
                      column_km = 1,
                      column_km_repeats = 1,
                      column_split = factor(annotations$cellType),
                      gap = unit(1, "mm"),
                      row_gap = unit(3, "mm"),
                      column_gap = unit(1, "mm"),
                      
                      heatmap_width = unit(1, "npc"),
                      width = NULL,
                      heatmap_height = unit(1, "npc"),
                      height = NULL,
                      
                      show_heatmap_legend = TRUE,
                      heatmap_legend_param = list(title = 'Exprs')
                      
  )
  return(heat_obj)
}





addTermsHeatmap <- function(cds, genes, meta, 
                            exprs_mat=NULL, 
                            dist_row = 'pearson', 
                            dist_col = 'euclidean', 
                            hmap_obj=NULL, 
                            org = 'mouse', 
                            k=NULL, 
                            universe = NULL, 
                            plot = F, no_log = F,
                            tpm = F, scale_minmax = T){
  ct <- cds[rowSums(cds > 0) > 1,]
  genes <- intersect(row.names(ct), genes)
  if(!is.null(exprs_mat)){
    genes <- intersect(row.names(exprs_mat), genes)
  }else{
    if(!tpm){
      if(scale_minmax){
        exprs_mat <- scale_down(log2(t(t(ct)/estimateSizeFactorsForMatrix(round(ct)))+1))
      }
      else{
        exprs_mat <- as.matrix(log10(t(t(ct)/estimateSizeFactorsForMatrix(round(ct)))+1))
      }
      #exprs_mat <- log10(t(t(ct)/estimateSizeFactorsForMatrix(round(ct)))+1)
      if(no_log){
        ct <- t(t(ct)/estimateSizeFactorsForMatrix(ct))
      }else{
        #ct <- log(t(t(ct)/estimateSizeFactorsForMatrix(ct))+1)
        ct <- DESeq2::varianceStabilizingTransformation(as.matrix(round(ct)))
      }
      
    }else{
      if(scale_minmax){
      exprs_mat <- scale_down(log2(ct+1))
      }else{
        exprs_mat <- as.matrix(log10(ct+1))
      }
      #exprs_mat <- log10(as.matrix(ct)+1)
      if(no_log){
        ct <- ct
      }else{
        ct <- ct <- log(ct+1)
      }
      
    }
  }
  
  mat <- ct[genes,]
  exprs_mat <- exprs_mat[genes,]
  meta$cellType <- factor(meta$cellType, levels = STAGE_ORDER)
  anno <- meta[,'cellType',drop =F]
  
  if(is.null(hmap_obj)){
    hmap_obj <- geneHeatmap(mat, exprs_mat = exprs_mat, annotations = anno, tpm = tpm, dist_row  = dist_row,  dist_col = dist_col, clust_n = 2)
    return(hmap_obj)
  }else{
    hmap_obj <- geneHeatmap(mat,exprs_mat  = exprs_mat,  annotations = anno, tpm = tpm, dist_col = dist_col, dist_row  = dist_row, clust_n = k)
  }
  
  hm_order = row_order(draw(hmap_obj))
  clust_l <- rep(0, length(genes))
  names(clust_l) <- 1:length(genes)
  for(i in 1:length(hm_order)){
    clust_l[hm_order[[i]]] <- i
  }
  genes_clust <- lapply(1:length(hm_order), FUN = function(x){genes[hm_order[[x]]]})
  names(genes_clust) <- 1:length(hm_order)

  gene_terms <- lapply(1:length(hm_order), FUN = function(x){
    sig_terms <- enrich_CP(genes[hm_order[[x]]], organisms = org,universe = universe, logFC = NULL, GSE = F, GO_BP_only = T)
  })
  test_gene_terms <- lapply(gene_terms, FUN= function(sig_terms){
    if('GO_BP_ora' %in% names(sig_terms))
    {
      if(sum(grepl('death', sig_terms$GO_BP_ora@result$Description[1:5]))){View(sig_terms$GO_BP_ora@result)}
      sig_terms$GO_BP_ora@qvalueCutoff <- 0.1
      sig_terms$GO_BP_ora@pvalueCutoff <- 0.1
      sig_terms <- subset(clusterProfiler::simplify(sig_terms$GO_BP_ora)@result, qvalue < 0.05 | pvalue < 0.001)$Description
    }else{
      return(NULL)
    }
    if(length(sig_terms) > 0){
      return(sig_terms[1:min(3, length(sig_terms))])
    }else{
      return(NULL)}
  })
  clust_l <- factor(clust_l)
  names(test_gene_terms) <-  levels(clust_l)
  test_gene_terms[sapply(test_gene_terms, is.null)] <- NULL
  if(length(intersect(names(clust_l), names(test_gene_terms))) == 0)
  {return(hmap_obj)}
  print(test_gene_terms)
  heatmap_obj <- hmap_obj+rowAnnotation(textbox = anno_textbox(clust_l, test_gene_terms, gp = gpar(fontsize = 11), word_wrap= T, text_space = unit(10, "pt"), max_width = unit(70, 'mm')))
  if(plot != F){
    png(file=plot, res = 300, width = 10,height = 10, units = 'in')
    draw(heatmap_obj)
    dev.off()
  }
  list(hmap = heatmap_obj, clust = genes_clust, enriched=gene_terms)
}

scale_down <- function(m, min_v = 0, max_v = 1){
  new_m <- t(apply(m, 1, FUN = function(x){
  x = (x-min(x))*(max_v-min_v)/(max(x)-min(x)) + min_v
  }))
  row.names(new_m) <- row.names(m)
  colnames(new_m) <- colnames(m)
  return(new_m)
}


sample_PCA <- function(cds, meta, 
                       reduced_dims = NULL, 
                       color_by = 'cellType', shape = NULL,
                       umap = F,
                       labeling = FALSE,
                       point_size = 4,
                       dimension = 2,
                       reduce_noise = FALSE,
                       return_matrix = FALSE,
                       highlight = NULL,
                       main = 'plot',
                       umap.config = umap_config)
{
  if (!umap) {
    test <- rsvd::rpca(t(cds), center = T, scale = T)
    percentVar <-
      c(
        100 * round(test$sdev[1] ^ 2 / sum(test$sdev ^ 2), 3),
        100 * round(test$sdev[2] ^ 2 / sum(test$sdev ^ 2), 3),
        100 * round(test$sdev[3] ^ 2 / sum(test$sdev ^ 2), 3)
      )
    axis_x = paste0('PC1: ', percentVar[1], "% variance")
    axis_y = paste0('PC2: ', percentVar[2], "% variance")
    axis_z = paste0('PC3: ', percentVar[3], "% variance")
    test <- test$x
  }
  else{
    umap_test <-
      umap(
        t(cds),
        n_components = dimension,
        n_neighbors = umap.config[['n_neighbors']],
        min_dist = umap.config[['min_dist']],
        metric = umap.config[['metric']],
        transform_state = umap.config[['transform_state']],
        random_state = umap.config[['random_state']]
      )
    axis_x = paste0('UMAP Comp 1')
    axis_y = paste0('UMAP Comp 2')
    axis_z = paste0('UMAP Comp 3')
    test <- umap_test$layout
  }
  alpha = 1
  if (dimension == 2)
  {
    theme0 <-
      theme_bw() + theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
      )
    scale_color_manual()
    if (color_by %in% c('cellType', 'cell_stage')) {
      dot_color <- DOT_COLOR
      colors_used <- dot_color[CONDITIONS_NAMES[as.character(unique(meta[[color_by]]))]]
      meta$condition <- CONDITIONS_NAMES[as.character(meta[[color_by]])]
      alpha = 1
    }
    if(!('experiment' %in% colnames(meta))){
      meta$experiment <- 'A'
    }
    if(!is.null(reduced_dims)){
      test <- reduced_dims
    }
    print(dim(test))
    test <- data.frame(test[row.names(meta), c(1, 2)], 
                       condition = meta[['condition']],  
                       shape = meta[['experiment']])
    #test <- data.frame(test[row.names(meta), c(1, 2)], 
    #condition = row.names(meta) %in% row.names(possible_somatic_only),  shape = 'experiment')
    if(!is.null(shape)){
      test$shape <- meta[[shape]]
    }
    colnames(test) <- c('PC1', 'PC2', 'Stage', 'Batch')
    #View(test)
    
    #return(test)
    if (labeling == TRUE)
    {
      plot0 <-
        ggplot(test,
               label = T,
               aes(PC1, PC2, color = Stage, label = row.names(test))) + geom_text_repel(size =
                                                                                          4) + geom_point(size = 0) + xlab(axis_x) +
        ylab(axis_y) + coord_fixed() + theme0 + ggtitle(main)
    }
    else
    {
      plot0 <-
        ggplot(test, label = T, aes(PC1, PC2, color = Stage, shape = Batch)) + geom_point( size = point_size, alpha = alpha) + xlab(axis_x) +
        ylab(axis_y) + coord_fixed() + theme0 + ggtitle(main)
    }
    if(color_by %in% c('cellType', 'cell_stage')) {
      plot0 <- plot0 + scale_color_manual(values = colors_used)
    } else{
      plot0 <- plot0
    }
    
  }
  else
  {
    test <-
      data.frame(test[, c(1, 2, 3)], condition = as.factor(meta[[color_by]]))
    colnames(test) <- c('PC1', 'PC2', 'PC3', 'condition')
    if (labeling == TRUE)
    {
      plot0 <- plot_ly(
        test,
        x = ~ PC1,
        y = ~ PC2,
        z = ~ PC3,
        color = ~ condition,
        colors = DOT_COLOR[as.character(unique(test$condition))],
        text = row.names(test)
      ) %>% add_text() %>% plotly::layout(scene = list(
        xaxis = list(title = axis_x),
        yaxis = list(title = axis_y),
        zaxis = list(title = axis_z)
      ))
    }
    else
    {
      plot0 <- plot_ly(
        test,
        x = ~ PC1,
        y = ~ PC2,
        z = ~ PC3,
        color = ~ condition,
        colors = DOT_COLOR[as.character(unique(test$condition))]
      ) %>% plotly::layout(scene = list(
        xaxis = list(title = axis_x),
        yaxis = list(title = axis_y),
        zaxis = list(title = axis_z)
      ))
    }
  }
  print(plot0)
  #return(test)
}


#### Revised version of how Monocle runs DE, to produce Log2FC values and also shrinkage of Log2FC using apeglm
monocle_diff <- function(cds, cells = NULL, subtype = NULL, num_cells = 3, min.pct = 0.1, full_formula, reduced_formula, foldChange = 'mle', apeglm_method = 'nbinomC'){
  cds0 <- cds
  
  if(!is.null(cells)){
    cds0 <- cds0[,cells]
  }
  
  
  
  
  genes <-  c()
  if(!is.null(subtype)){
    if(length(subtype) == 1){
      pData(cds0)[pData(cds0)$cellType != subtype,]$cellType <- 'other'
      subtype <- c(subtype, 'other')
    }else{
      cds0 <- cds0[,row.names(subset(pData(cds0), cellType %in% subtype))]
    }
    for(s in subtype){
      ct <- exprs(cds0[,row.names(subset(pData(cds0), cellType == s))])
      print(dim(ct))
      print(round(ncol(ct)*min.pct))
      genes <- union(genes, row.names(ct[rowSums(ct > 0)>= max(num_cells, round(ncol(ct)*min.pct)),]))
      print(length(genes))
    }
    
  }else{
    genes <- row.names(fData(cds0))
  }
  
  cds0 <- cds0[genes,]
  cds0 <- detectGenes(cds0, 0)
  cds0 <- cds0[fData(cds0)$num_cells_expressed >= num_cells*2,]
  if(grepl('num_genes_expressed', full_formula)){
    pData(cds0)$num_genes_expressed <- scale(pData(cds0)$num_genes_expressed)
  }
  cds0 <- estimateSizeFactors(cds0)
  cds0 <- estimateDispersions(cds0)
  ptm <- proc.time()
  #diff_test <- differentialGeneTest(cds = cds0, fullModelFormulaStr = full_formula, reducedModelFormulaStr = reduced_formula)
  diff_test_res <- smartEsApply(cds0,1,diff_test_helper,
                                convert_to_dense=TRUE,
                                fullModelFormulaStr=full_formula,
                                reducedModelFormulaStr=reduced_formula, 
                                expressionFamily=cds0@expressionFamily, 
                                relative_expr=T,
                                disp_func=cds0@dispFitInfo[["blind"]]$disp_func,
                                verbose=F)
  res <- do.call(rbind.data.frame, lapply(diff_test_res, FUN = function(x){x$res}))
  
  res$qval <- 1
  res$qval[which(res$status == 'OK')] <- p.adjust(subset(res, status == 'OK')[, 'pval'], method="BH")
  
  res <- merge(res, fData(cds0), by="row.names")
  row.names(res) <- res[, 1] #remove the first column and set the row names to the first column
  res[, 1] <- NULL 
  
  res <- res[row.names(cds0), ]
  #cell_models0 <- fit_models(cds0,
  #                model_formula_str = full_formula,
  #               expression_family="negbinomial")
  #null_models0 <- fit_models(cds0,
  #                      model_formula_str = reduced_formula,
  #                     expression_family="negbinomial")
  #diff_test <- compare_models(cell_models0, null_models0)
  #res <- subset(diff_test, q_value < 0.05)
  
  
  #res <- diff_test
  print('diff test done')
  print(proc.time() - ptm)
  if(!is.null(subtype) & length(subtype) <= 2){
    # most direct calculation of log2FC, add pseudocounts of 1 for stability and slight shrunkage of values
    mat <- t(t(exprs(cds0))/sizeFactors(cds0))
    MeanExp <- tapply(colnames(mat), as.character(pData(cds0)[colnames(mat),]$cellType), function(x){rowMeans(mat[,x])})
    MeanExp <- log2((MeanExp[[subtype[2]]]+1)/(MeanExp[[subtype[1]]]+1))
    res$Log2FC_direct <- MeanExp[row.names(res)]
    
    if(foldChange == 'vst'){# using vst values to calculate log2FC, not recommended by authors
      mat <- vstExprs(cds0)
      MeanExp <- tapply(colnames(mat), as.character(pData(cds0)[colnames(mat),]$cellType), function(x){rowMeans(mat[,x])})
      MeanExp <- MeanExp[[subtype[2]]]- MeanExp[[subtype[1]]]
      res$Log2FC <- MeanExp[row.names(res)]
    }else if(foldChange == 'mle'){
      for(l in c("BiocGenerics", "VGAM", "Matrix")){library(l, character.only = T)}
      #ptm <- proc.time()
      #diff_test_res <- smartEsApply(cds0,1,diff_test_helper,
      #                      convert_to_dense=TRUE,
      #                      fullModelFormulaStr=full_formula,
      #                      reducedModelFormulaStr=reduced_formula, 
      #                      expressionFamily=cds0@expressionFamily, 
      #                      relative_expr=T,
      #                      disp_func=cds0@dispFitInfo[["blind"]]$disp_func,
      #                      verbose=F)
      #print('models fitted for effect extraction')
      #print(proc.time() - ptm)
      return(logfc_shrink(list(cds=cds0, res = res, models=diff_test_res), subtype = subtype, apeglm_method = apeglm_method ))
      
    }
    #res <- subset(res, abs(Log2FC) >= 1)
  }
  
  #res <- res[order(res$q_value),]
  #res <- res[order(res$qval),]
  return(res)
  #res <- subset(diff_test, qval < 0.05)
  #return(res)
}




logfc_shrink <- function(obj, subtype, apeglm_method = 'nbinomC'){
  res <- obj$res
  cds0 <- obj$cds
  diff_test_res <- obj$models
  oppose_signs = ifelse(sum(grepl(subtype[2], names(coef(diff_test_res[[1]]$full)))) > 0, F, T)
  ptm = proc.time()
  Log2FC_effect_disp <- do.call(rbind, lapply(row.names(res), FUN = function(x){
    disp = diff_test_res[[x]]$disp
    #print('1')
    coefs <- coef(diff_test_res[[x]]$full)
    coefs_SE <- sqrt(diag(vcov(diff_test_res[[x]]$full)))
    #print('2')
    if(!oppose_signs){
      effect <- coefs[grepl(subtype[2], names(coefs))]
      effect_SE <- coefs_SE[grepl(subtype[2], names(coefs_SE))]
    }else{
      effect <- -coefs[grepl(subtype[1], names(coefs))]
      effect_SE <- coefs_SE[grepl(subtype[1], names(coefs_SE))]
    }
    #print('3')
    return(c(effect, effect_SE, disp))
  }))
  print('Effects and Errors calculated')
  print(proc.time() - ptm)
  res$Log2FC_effect <- Log2FC_effect_disp[,1]
  res$Log2FC_effect_SE <- Log2FC_effect_disp[,2]
  res$disp <- Log2FC_effect_disp[,3]
  design <- diff_test_res[[1]]$full@x
  mle_prior=res[,c('Log2FC_effect', 'Log2FC_effect_SE')]
  if(oppose_signs){
    mle_prior[,1] <- -mle_prior[,1]
  }
  print('starting apeglm for logfc shrinkage')
  apeglm_mat <- exprs(cds0)[row.names(res),]
  ptm = proc.time()
  
  fit <- apeglm(Y= apeglm_mat, x=design, log.lik=logLikNB, param=as.matrix(res[,'disp', drop = F]),
                coef=2, mle=mle_prior, offset =  matrix(log(sizeFactors(cds0)), nrow=nrow(apeglm_mat), ncol=ncol(apeglm_mat), byrow=TRUE), method = apeglm_method)
  
  print('effects shrunk with apeglm')
  print(proc.time() - ptm)
  
  res$Log2FC_shrunk <- fit$map[row.names(res),2]*log2(exp(1))
  res$lfc_sval <- fit$svalue[row.names(res),]
  res$lfc_fsr <- fit$fsr[row.names(res),]
  
  if(oppose_signs){
    res$Log2FC_shrunk <- -res$Log2FC_shrunk
  }
  res$Log2FC_unshrunk <- Log2FC_effect_disp[,1]*log2(exp(1))
  
  res$emp_disp <- data.frame(cds0@dispFitInfo$blind$disp_table, row.names = 1)[row.names(res),'disp']
  res$est_disp <- cds0@dispFitInfo$blind$disp_func(data.frame(cds0@dispFitInfo$blind$disp_table, row.names = 1)[row.names(res),'mu'])
  return(res)
}

sparseApply <- function(Sp_X, MARGIN, FUN, convert_to_dense, ...){
  
  if (convert_to_dense){
    if (MARGIN == 1){
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[,i]), ...) 
      }, FUN, ...)
    }else{
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(as.matrix(Sp_X[,i]), ...) 
      }, FUN, ...)
    }
  }else{
    if (MARGIN == 1){
      Sp_X <- Matrix::t(Sp_X)
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[,i], ...) 
      }, FUN, ...)
    }else{
      res <- lapply(colnames(Sp_X), function(i, FUN, ...) {
        FUN(Sp_X[,i], ...) 
      }, FUN, ...)
    }
  }
  
  return(res)
  
}

smartEsApply <- function(X, MARGIN, FUN, convert_to_dense, ...) {
  isSparseMatrix <- function(x){
    class(x) %in% c("dgCMatrix", "dgTMatrix")
  }
  parent <- environment(FUN)
  if (is.null(parent))
    parent <- emptyenv()
  e1 <- new.env(parent=parent)
  Biobase::multiassign(names(pData(X)), pData(X), envir=e1)
  environment(FUN) <- e1
  
  if (sum(isSparseMatrix(exprs(X)))){
    res <- sparseApply(exprs(X), MARGIN, FUN, convert_to_dense, ...)
  }else{
    res <- apply(exprs(X), MARGIN, FUN, ...)
  }
  
  if (MARGIN == 1)
  {
    names(res) <- row.names(X)
  }else{
    names(res) <- colnames(X)
  }
  
  res
}


diff_test_helper <- function(x, 
                             fullModelFormulaStr, 
                             reducedModelFormulaStr, 
                             expressionFamily, 
                             relative_expr,
                             weights,
                             disp_func=NULL,
                             verbose=FALSE
){ 
  calculate_NB_dispersion_hint <- function(disp_func, f_expression, expr_selection_func=mean)
  {
    expr_hint <- expr_selection_func(f_expression)
    if (expr_hint > 0 && is.null(expr_hint) == FALSE) {
      disp_guess_fit <- disp_func(expr_hint)
      
      # For NB: Var(Y)=mu*(1+mu/k)
      f_expression_var <- var(f_expression)
      f_expression_mean <- mean(f_expression)
      
      disp_guess_meth_moments <- f_expression_var - f_expression_mean
      disp_guess_meth_moments <- disp_guess_meth_moments / (f_expression_mean^2) #fix the calculation of k
      
      #return (max(disp_guess_fit, disp_guess_meth_moments))
      return (disp_guess_fit)
    }
    return (NULL)
  }
  
  reducedModelFormulaStr <- paste("f_expression", reducedModelFormulaStr, sep="")
  fullModelFormulaStr <- paste("f_expression", fullModelFormulaStr, sep="")
  
  x_orig <- x
  disp_guess <- 0
  
  if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")){
    if (relative_expr == TRUE)
    {
      x <- x / Size_Factor
    }
    f_expression <- round(x)
    if (is.null(disp_func) == FALSE){
      disp_guess <- calculate_NB_dispersion_hint(disp_func, round(x))
      if (is.null(disp_guess) == FALSE && disp_guess > 0 && is.na(disp_guess) == FALSE  ) {
        # FIXME: In theory, we could lose some user-provided parameters here
        # e.g. if users supply zero=NULL or something. 
        if (expressionFamily@vfamily == "negbinomial")
          expressionFamily <- negbinomial(isize=1/disp_guess)
        else
          expressionFamily <- negbinomial.size(size=1/disp_guess)
      }
    }
  }else if (expressionFamily@vfamily %in% c("uninormal")){
    f_expression <- x
  }else if (expressionFamily@vfamily %in% c("binomialff")){
    f_expression <- x
    #f_expression[f_expression > 1] <- 1
  }else{
    f_expression <- log10(x)
  }
  
  test_res <- tryCatch({
    if (expressionFamily@vfamily %in% c("binomialff")){
      if (verbose){
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily)                         
      }else{
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily))                    
      }
    }else{
      if (verbose){
        full_model_fit <- VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily)
        reduced_model_fit <- VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily)                         
      }else{
        full_model_fit <- suppressWarnings(VGAM::vglm(as.formula(fullModelFormulaStr), epsilon=1e-1, family=expressionFamily))
        reduced_model_fit <- suppressWarnings(VGAM::vglm(as.formula(reducedModelFormulaStr), epsilon=1e-1, family=expressionFamily))                    
      }
    }
    
    #list(full=full_model_fit, reduce = reduced_model_fit, disp = disp_guess, res = res)
    #print(coef(reduced_model_fit))
    res <- compareModels(list(full_model_fit), list(reduced_model_fit))
    list(full=full_model_fit, reduce = reduced_model_fit, disp = disp_guess, res = res)
  }, 
  #warning = function(w) { FM_fit },
  error = function(e) { 
    if(verbose)
      print (e);
    list(full=NULL, reduce = NULL, disp = disp_guess, res = data.frame(status = "FAIL", family=expressionFamily@vfamily, pval=1.0, qval=1.0))
    #data.frame(status = "FAIL", pval=1.0) 
  }
  )
  test_res
}

plot_top_diff_genes <- function(cds, 
         diff_res = NULL, 
         genes = NULL, 
         num_genes = 20, 
         compare = NULL, 
         color_by='cellType', 
         plot_type = c('trajectory', 'jitter', 'violin', 'bar', 'diagram'), 
         grouping='stage', ncols = 5, sep_plots = F, sep_plots_path = NULL){
  theme0 <- theme_bw() + theme(plot.title = element_text(hjust = 0.5, size=26), panel.border = element_blank(), panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),  axis.text.y = element_text(size=18), 
                               axis.text.x = element_text(size=20), legend.position = "none", strip.text.x = element_text(size = 15))
  if(!is.null(compare)){
    cds <- cds[,cds$cellType %in% compare]
  }
  if(is.null(genes)){
    qvals <- subset(diff_res, q_value < 0.05)$q_value
    names(qvals)<- subset(diff_res, q_value < 0.05)$gene_id
    top_genes <- names(sort(qvals)[1:num_genes])
    cds_subset <- cds[top_genes,]
  }else{
    cds_subset <- cds[fData(cds)$gene_short_name %in% genes,]
    
  }
  if(!is.null(diff_res)){
    diff_sig <- check_diff(row.names(fData(cds_subset)), diff_res)
  }else{
    diff_sig = NULL
  }
  if(plot_type == 'diagram'){
    q <- plot_genes_diagram(cds_subset, 'cellType', ncols = ncols, sep_plots = sep_plots, sep_plots_path = sep_plots_path, diff_sig = diff_sig)
    return(q)
  }
  if(plot_type == 'bar'){
    df <- data.frame(t(exprs(cds_subset)))
    df$celltype <- cds_subset$celltype
  }
  if(plot_type == 'trajectory'){
    q <- plot_genes_in_pseudotime(cds_subset, color_cells_by =color_by)
  }else if(plot_type == 'jitter'){
    q <- plot_genes_in_pseudotime(cds_subset,color_by= color_by, min_expr=0.5, cell_size = 2, ncol = 1 )+
      scale_x_discrete(name ="Stage", limits=c('S1', 'S2', 'S3', 'S4', 'F3', 'F2', 'F1', 'P0'))+
      scale_color_manual(values = DOT_COLOR[unique(pData(cds)$cellType)])
  }else if(plot_type == 'violin'){
    q <- plot_genes_violin(cds_subset, color_by =color_by, grouping = grouping, log_scale = T)+scale_x_continuous(breaks=1:length(unique(pData(cds_subset)[,grouping])))
  }
  q <- q+ theme0+theme(legend.position = "none")
  q
}


plot_genes_diagram <- function(cds, color_by = color_by, single_scale = F, ncols = 5, sep_plots = F, sep_plots_path = NULL, diff_sig = NULL){

  average_exprs_group <- do.call(rbind, lapply(row.names(fData(cds)), function(g){
    tapply(row.names(pData(cds)), pData(cds)[,color_by], function(c){mean(log10(t(t(exprs(cds))/sizeFactors(cds))[g,c]+1))})
  }))
  row.names(average_exprs_group) <- row.names(fData(cds))
    
  svg  <- grImport2::readPicture('./diagram_cairo.svg')
  #nam_content <- list('S1' = c(1), 'S2' = c(2), 'S3' = c(3), 'S4' = c(9, 10, 11, 8), '-3' = c(5), '-2'= c(6), '-1' = c(7), 'P0' = c(26))
  
  nam_content <- list('S1' = c(1), 'S2' = c(3), 'S3' = c(5), 'S4' = c(17, 19, 15, 21), '-3' = c(9), '-2'= c(11), '-1' = c(13), 'P0' = c(37))
  #nam_content <- list('S1' = c(1), 'S2' = c(2), 'S3' = c(3), 'S4' = c(8,9,10,11), 'F3' = c(5), 'F2'= c(6), 'F1' = c(7), 'P0' = c(19))
  grob_list <- list()
  for(g in row.names(average_exprs_group)){
    gene_name <- fData(cds)[g, 'gene_short_name']
    new_svg <- svg
    rang <- c(floor(min(average_exprs_group[g,])), max(max(average_exprs_group[g,]), floor(min(average_exprs_group[g,]))+2))
    #print(rang)
    ii <- cut(average_exprs_group[g,], breaks = seq(floor(rang[1]), ceiling(rang[2]), len = 100), include.lowest = TRUE)
    
    cols = colorRampPalette(brewer.pal(9,"Reds"))(100)[ii]
    #print(ii)
    names(cols) <- colnames(average_exprs_group)
    #print(cols)
    for(cond in colnames(average_exprs_group)){
      for(i in nam_content[[cond]]){
        new_svg@content[[1]]@content[[i]]@gp$fill <- cols[cond]
      }
      for(i in setdiff(1:length(new_svg@content[[1]]@content), do.call(c, nam_content))){
        if(i == 7){
          new_svg@content[[1]]@content[[i]]@gp$fill <- 'white'
        }
        else{
          new_svg@content[[1]]@content[[i]]@gp$fill <- 'black'
        }
      }
    }
    legd <- ggplotify::base2grob(~color.bar(colorRampPalette(brewer.pal(9,"Reds"))(100), floor(rang[1]), ceiling(rang[2]), title = expression('Log'[10]*'Exprs'), sub_title = gene_name))
    if(!is.null(diff_sig)){
    err.bar <- ggplotify::base2grob(~error.bar(diff_sig[[g]]))
    }
    #gg <- packGrob(packGrob(packGrob(frameGrob(), pictureGrob(new_svg, width = 0.8, height = 1.5)), err.bar, width = unit(0.2, 'null'), side = 'left'),
     #              legd, height=unit(1,"null"), width = unit(0.45, 'null'), side="right")
    
    gg <- packGrob(packGrob(packGrob(frameGrob(), pictureGrob(new_svg, width = 0.6, height = 1.2)), err.bar, width = unit(0.1, 'null'), side = 'right'),
                   legd, height=unit(1,"null"), side="right", width = unit(0.9, 'null'))
    grob_list[[g]] <- gg
  }
  cowplot::plot_grid(plotlist = grob_list,  ncol =ncols)
}

check_diff <- function( gene, res){
  chk <- lapply(gene, function(x){
    sapply(res, function(r){
      if(is.na(r[x, 'qval']) | abs(r[x,'Log2FC_shrunk']) < log2(1.5) | r[x, 'qval'] > 0.05){
        return(NA)
      }
      sig <- '*'
      if(r[x, 'qval'] < 0.01){
        sig <- '**'
      }
      if(r[x, 'qval'] < 0.001){
        sig <- '***'
      }
      return(sig)
    })
  })
  names(chk) <- gene
  return(chk)
}




color.bar <- function(lut, min, max=-min, nticks=5, ticks=seq(min, max, len=nticks), title='', sub_title = '') {
  scale = (length(lut)-1)/(max-min)
  par(mar = c(12, 14, 12, 1))
  #dev.new(width=1.75, height=5)
  p <- plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
  axis(2, ticks, las=1, font = 2)
  for (i in 1:(length(lut)-1)) {
    
    y = (i-1)/scale + min
    
    rect(0,y,4,y+1/scale, col=lut[i], border=NA)
  }
  mtext(bquote(bolditalic(.(sub_title))), cex = 1.7, side = 3, padj = -0.5, at = -2, adj = 0.5, font = 2 )
  mtext(title, cex = 1.3, adj = 0, side = 1, at =-12, padj = 1, font = 2)
}




error.bar <- function(all) {
  p <- plot(c(16,25), c(0, 10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
  brac_l = 3.5
  if(!is.na(all['F1_P0'])){
    x= -10
    y = c(-1.5, 0.4)
    segments(x,y[1], x, y[2], lwd = 2.3)
    segments(x-brac_l,y[1], x, y[1], lwd = 2.3)
    segments(x-brac_l,y[2], x, y[2], lwd = 2.3)
    text(x, y = mean(y), labels = all['F1_P0'], srt = 90, adj = c(0.5, 1.5))
  }
  
  if(!is.na(all['F2_F1'])){
    x= -10
    y = c(0.6, 1.4)
    segments(x,y[1], x, y[2], lwd = 2.3)
    segments(x-brac_l,y[1], x, y[1], lwd = 2.3)
    segments(x-brac_l,y[2], x, y[2], lwd = 2.3)
    text(x, y = mean(y), labels = all['F2_F1'], srt = 90, adj = c(0.5, 1.5))
  }
  if(!is.na(all['F3_F2'])){
    x= -10
    y = c(1.6, 2.3)
    segments(x,y[1], x, y[2], lwd = 2.3)
    segments(x-brac_l,y[1], x, y[1], lwd = 2.3)
    segments(x-brac_l,y[2], x, y[2], lwd = 2.3)
    text(x, y = mean(y), labels = all['S3_S4'], srt = 90, adj = c(0.5, 1.5))
  }
  if(!is.na(all['S4_F3'])){
    x= -10
    y = c(2.5, 3.8)
    segments(x,y[1], x, y[2], lwd = 2.3)
    segments(x-brac_l,y[1], x, y[1], lwd = 2.3)
    segments(x-brac_l,y[2], x, y[2], lwd = 2.3)
    text(x, y = mean(y), labels = all['S4_F3'], srt = 90, adj = c(0.5, 1.5))
  }
  if(!is.na(all['S3_S4'])){
    x= -10
    y = c(4, 5.7)
    segments(x,y[1], x, y[2], lwd = 2.3)
    segments(x-brac_l,y[1], x, y[1], lwd = 2.3)
    segments(x-brac_l,y[2], x, y[2], lwd = 2.3)
    text(x, y = mean(y), labels = all['S3_S4'], srt = 90, adj = c(0.5, 1.5))
  }
  if(!is.na(all['S2_S3'])){
    x= -10
    y = c(6, 7.9)
    segments(x,y[1], x, y[2], lwd = 2.3)
    segments(x-brac_l,y[1], x, y[1], lwd = 2.3)
    segments(x-brac_l,y[2], x, y[2], lwd = 2.3)
    text(x, y = mean(y), labels = all['S2_S3'], srt = 90, adj = c(0.5, 1.5))
  }
  if(!is.na(all['S1_S2'])){
    x= -10
    y = c(8.4, 10.6)
    segments(x,y[1], x-6, y[2], lwd = 2.3)
    segments(x-brac_l,y[1], x, y[1], lwd = 2.3)
    segments(x-brac_l-6,y[2], x-6, y[2], lwd = 2.3)
    text(x-1.5, y = mean(y), labels = all['S1_S2'], srt = 95, adj = c(0.5, 1.5))
  }
  
  
}



color.bar2 <- function(lut, min, max=-min, nticks=5, ticks=seq(min, max, len=nticks), title='', sub_title = '') {
  scale = (length(lut)-1)/(max-min)
  par(mar = c(14, 1, 10, 1))
  #dev.new(width=1.75, height=5)
  p <- plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
  axis(2, ticks, las=1, font = 2)
  for (i in 1:(length(lut)-1)) {
    
    y = (i-1)/scale + min
    
    rect(0,y,4,y+1/scale, col=lut[i], border=NA)
  }
  mtext(bquote(bolditalic(.(sub_title))), cex = 1.6, side = 3, padj = -1, at = -12, adj = 0, font = 2 )
  mtext(title, cex = 1.2, adj = 0, side = 1, at =-8, padj = 1, font = 2)
}
error.bar2 <- function() {
  p <- plot(c(0,5), c(0, 10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
  segments(5,-1, 5,11, lwd = 2.3)
}

test_plot <- function(){
  png('test_2.png', width = 2.7, height = 5, units = 'in', res= 300)
  plot_top_diff_genes(ce.ct, genes = c('ced-8'), plot_type = 'diagram', ncols = 1)
  dev.off()
}

# Install and load reshape2 package
install.packages("reshape2")
library(reshape2)

# creating correlation matrix
plot_cor_heatmap <- function(cor_mat, other_data_name, plot_name){
  corr_mat <- round(cor_mat,5)
  NA_mat = matrix(NA, nrow(cor_mat), ncol(cor_mat))
  diag(NA_mat) = diag(corr_mat)
  for(r in 1:nrow(NA_mat)){
    pos = which(corr_mat[r,] == max(corr_mat[r,]))
    NA_mat[r,pos] = max(corr_mat[r,pos])
  }
  
  for(c in 1:ncol(NA_mat)){
    pos = which(corr_mat[,c] == max(corr_mat[,c]))
    NA_mat[pos,c] = max(corr_mat[pos,c])
  }
  
  max_cor = round(max(corr_mat), 1)+0.11
  min_cor = round(min(corr_mat), 1)-0.11
  # reduce the size of correlation matrix
  melted_corr_mat <- reshape2::melt(corr_mat)
  colnames(melted_corr_mat)[3] = 'Correlation'
  
  melted_corr_mat$Sig_cor = reshape2::melt(NA_mat)$value
  # plotting the correlation heatmap
  print(head(melted_corr_mat))
  ggplot(data = melted_corr_mat, aes(x=Var1, y=Var2,fill=Correlation, label=round(Sig_cor,3))) +
    geom_tile() +labs(x = other_data_name, y = 'Our Data', fill = "Pearson's\nCorrelation", title=paste("Correlations with ", other_data_name)) +scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", midpoint = (min_cor+max_cor)/2, limits=c(min_cor,max_cor))+geom_text()+
    theme_classic() +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) 
  #theme(text=element_text(family="Roboto"))
  ggsave(plot_name, width = 5, height= 4, units = 'in', dpi = 300)
}



#
read_rmats_res <- function(dirc){
  res_list <- list()
  for(type in c('A3SS', 'A5SS', 'MXE', 'SE', 'RI')){
    res_list[[type]] <- subset(read.csv(paste(paste(dirc, type, sep = '/'),'.MATS.JC.txt', sep =''), sep = '\t', header = T), FDR < 0.05 & abs(IncLevelDifference) > 0.05)
  }
  res_list[['universe']] = lapply(c('A3SS', 'A5SS', 'MXE', 'SE', 'RI'), FUN = function(x){
    unique(read.csv(paste(paste(dirc, 'fromGTF.', sep = '/'),x,'.txt', sep =''), sep = '\t', header = T)[,'GeneID'])
  })
  names(res_list[['universe']]) <- c('A3SS', 'A5SS', 'MXE', 'SE', 'RI')
  return(res_list)
}



deg_utr <- function(file, ct, compare, meta, impute = F, method = 'fisher.test', alpha = 0.05, combine_p = NULL, diff_thresh = 0.05){
  require(betareg)
  require(tidyr)
  require(qvalue)
  require(pbapply)
  require(lmtest)
  ##read in the new dapars file that includes long. short and PDUI values
  dapars <- read.csv(file, sep = '\t', header = T, row.names = 1)
  dapars <- dapars[,!grepl(regex('^X'), colnames(dapars))]
  dapars <- cbind(dapars[,1:3],dapars[,grepl(regex(paste(utr_samples, collapse = '|')), colnames(dapars))])
  colnames(dapars) <- str_remove(colnames(dapars), '^X')
  samples <- row.names(subset(meta, cellType %in% compare))
  dapars <- cbind(dapars[,c(1,2,3)], dapars[,grepl(regex(paste(samples, collapse = '|')), colnames(dapars))])
  
  ct <- ct[,samples]
  ## Filter first based on coverage of the each gene's entire body
  ## Only account for UTR coverage in genes that are assigned at have mean expression of at least 5 uniquely mapped reads across all samples
  genes <- row.names(ct)[rowMeans(ct) > 10]
  #genes <- intersect(genes, row.names(ct)[rowSums(ct[,row.names(subset(meta, cellType == compare[2]))] > 10) > 0])
  dapars$gene_short_names <- sapply(row.names(dapars), FUN = function(x){strsplit(x,"\\|")[[1]][2]})
  dapars_orig <- dapars
  all_genes <- unique(dapars$gene_short_names)
  cat('Total genes assessed: ', length(unique((dapars$gene_short_names))), 'genes\n')
  dapars <- dapars[dapars$gene_short_names %in% genes,]
  cat('filtering based on gene expression: left with ', length(unique((dapars$gene_short_names))), 'genes\n')
  #dapars <- subset(dapars, fit_value >= 10) # Maybe filter also based on regression fit value
  cat('filtering based on fit value: left with ', length(unique((dapars$gene_short_names))), 'genes\n')
  
  # vector to store gene name and UTR region correspondence in case gene names is lost with imputation
  gene2region <- dapars$gene_short_names
  names(gene2region) <- sapply(row.names(dapars), FUN = function(x){strsplit(x,"\\|")[[1]][1]})
  
  #split file into long, short and pdui
  d_long <- dapars[,grepl('long_exp', colnames(dapars))]
  d_short <- dapars[,grepl('short_exp', colnames(dapars))]
  d_pdui <- dapars[,grepl('PDUI', colnames(dapars))]
  
  # change column names
  colnames(d_long) <- sapply(strsplit(colnames(d_long), '_'), FUN = function(x){strsplit(x[1], "\\.")[[1]][1]})
  colnames(d_short) <- sapply(strsplit(colnames(d_short), '_'), FUN = function(x){strsplit(x[1], "\\.")[[1]][1]})
  colnames(d_pdui) <- sapply(strsplit(colnames(d_pdui), '_'), FUN = function(x){strsplit(x[1], "\\.")[[1]][1]})
  # select filtering genes based on number of passes (Non NAs) and overall coverage (average > 2 either in long or short UTR in both conditions)
  grp = meta[colnames(d_pdui),'cellType']
  btch = meta[colnames(d_pdui), 'experiment']
  cond1_ind = which(grp == compare[1])
  cond2_ind = which(grp == compare[2])
  #print(d_pdui[dapars$gene_short_name == 'Cdk1',])
  ## NA FILTER
  # onyly keep genes that have at least than 3 non NAs in terms of coverage in both conditions
  na.filt.genes <- rowSums(!is.na(d_pdui[,cond1_ind])) >= 3 & rowSums(!is.na(d_pdui[,cond2_ind])) >= 3
  
  dapars <- dapars[na.filt.genes, ]
  
  cat('filtering based on number of non NAs in each condition ( > 3 in each condition): left with ', length(unique((dapars$gene_short_names))), 'genes\n')
  
  
  d_long <- d_long[na.filt.genes, ]
  d_short <- d_short[na.filt.genes, ]
  d_pdui <- d_pdui[na.filt.genes, ]
  
  # Filter based on coverage of UTR regions
  c1.filt.genes <- rowMeans(d_long[,cond1_ind], na.rm = T) > 1 | rowMeans(d_short[,cond1_ind], na.rm = T) > 1
  c2.filt.genes <- rowMeans(d_long[,cond2_ind], na.rm = T) > 1 | rowMeans(d_short[,cond2_ind], na.rm = T) > 1
  
  final.filt.genes <- c1.filt.genes & c2.filt.genes 
  
  # filtering matrices with genes selected prior
  dapars <- dapars[final.filt.genes, ] 
  d_long <- d_long[final.filt.genes, ] %>%  mutate(across(everything(), replace_na, 0))
  d_short <- d_short[final.filt.genes, ] %>%  mutate(across(everything(), replace_na, 0))
  d_pdui <- d_long/(d_long+d_short)
  cat('filtering based on number of mean coverage of long/short UTR in both condition: left with ', length(unique((dapars$gene_short_names))), 'genes\n')
  if(impute){
    dapars_out <- data.frame(Gene = row.names(dapars), dapars[,1:3], d_pdui)
    write.table(dapars_out, file = './temp.dp.tsv', sep = '\t', quote = F, row.names = F)
    d_pdui = scDaPars(raw_PDUI_file = './temp.dp.tsv',
                      out_dir = "apa/scDaPars_result",
                      filter_gene_thre = 0.2,
                      filter_cell_thre = 0.1)
    return(d_pdui)
    method = 'ks.test'
  }
  cat('performing tests\n')
  if(method == 'ks.test'){
    test <- apply(d_pdui, 1, FUN = function(x){
      if(sum(x[!is.na(x)]) == 0 | sum(!is.na(x[cond1_ind])) <= 3 | sum(!is.na(x[cond2_ind])) <= 3){
        c(1, 0)
      }else{
        c(ks.test(x[cond1_ind][!is.na(x[cond1_ind])], x[cond2_ind][!is.na(x[cond2_ind])])$p.value, mean(x[cond2_ind][!is.na(x[cond2_ind])]) - mean(x[cond1_ind][!is.na(x[cond1_ind])]))
      }
    })
    test <- t(test)
  }
  else if(method ==  'binom'){
    test <- do.call(rbind, pblapply(row.names(d_pdui),FUN = function(n){
      x <- d_pdui[n,]
      w <- d_long[n,]+d_short[n,]
      not_na <- !is.na(x)
      x0 <- c(x[not_na])
      w0 <- c(w[not_na])
      data <- data.frame(pdui = x0, cellType = grp[not_na], batch = btch[not_na])
      test_0 <- tryCatch({
        mylogit <- glm(pdui ~ cellType, data = data, family = "quasibinomial", weights =w0)
        mylogit0 <-glm(pdui ~ 1, data = data, family = "quasibinomial", weights = w0)
        pval <- anova(mylogit, mylogit0, test='F')[,"Pr(>F)"][2]
        diff <- predict(mylogit, data.frame(cellType = compare[2], batch = 'A'), type = 'response') - predict(mylogit, data.frame(cellType = compare[1], batch = 'A'), type = 'response')
        return(c(pval, diff))
      },error = function(cond) {
        return(c(1, 0))
      })
      return(test_0)
    }))
    test[is.na(test[,1]),] <- c(1,0)
  }else if(method ==  'betab'){
    test <-  do.call(rbind, pblapply(row.names(d_pdui),FUN = function(n){
      xx <- d_pdui[n,]
      x <- d_long[n,]
      w <- d_long[n,]+d_short[n,]
      not_na <- !is.na(xx)
      x0 <- round(c(x[not_na]))
      w0 <- round(c(w[not_na]))
      data <- data.frame(y = x0, n = w0, cellType = grp[not_na], batch = btch[not_na])
      if(sum(x0 == w0) == length(x0)){
        return(c(1,0))
      }
      test_0 <- tryCatch({
        mylogit <- betabin(formula=cbind(y, n - y) ~ cellType,random= ~ 1, data = data)
        mylogit0 <-betabin(formula =cbind(y, n - y) ~ 1, random=~ 1, data = data)
        pval <- anova(mylogit, mylogit0)@anova.table[,"P(> Chi2)"][2]
        preds = predict(mylogit, data.frame(cellType=c(compare[2],compare[1])))
        diff <- preds[2] - preds[1]
        return(c(pval, diff))
      },error = function(cond) {
        print(cond)
        print(data)
        return(c(1, 0))
      })
      return(test_0)
    }))
    test[is.na(test[,1]),] <- c(1,0)
  }
  else if(method == 'betareg'){
    test <- do.call(rbind, pblapply(row.names(d_pdui),FUN = function(n){
      x <- d_pdui[n,]
      not_na <- !is.na(x)
      x0 <- c(x[not_na])
      if(sum(x0 == 0) > 0 | sum(x0 == 1) > 0){
        x0 <- (x0*(length(x0)-1)+0.5)/length(x0)
      }
      if(length(unique(x0)) == 1){
        return(c(1,0))
      }
      data <- data.frame(pdui = x0, cellType = grp[not_na], batch = btch[not_na])
      data$cellType <- as.factor(data$cellType)
      test_0 <- tryCatch({
        mylogit <- betareg(pdui ~ cellType, data = data, link = 'log')
        mylogit0 <- betareg(pdui ~ 1, data = data, link = 'log')
        pval <- lrtest(mylogit, mylogit0)[,"Pr(>Chisq)"][2]
        diff <- predict(mylogit, data.frame(cellType = compare[2])) - predict(mylogit, data.frame(cellType = compare[1]))
        return(c(pval, diff))
      },error = function(cond) {
        print(cond)
        return(c(1, 0))
      })
      return(test_0)
    }))
  }
  else{
    l1 = length(cond1_ind)
    l2 = length(cond2_ind)
    utrl1_mean <- round(rowMeans(d_long[, cond1_ind], na.rm = T))
    utrl2_mean <- round(rowMeans(d_long[, cond2_ind], na.rm = T))
    utrs1_mean <- round(rowMeans(d_short[, cond1_ind], na.rm = T))
    utrs2_mean <- round(rowMeans(d_short[, cond2_ind], na.rm = T))
    test <- do.call(rbind, pblapply(row.names(d_long), FUN = function(x){
      #utr_l1 <- d_long[x,cond1_ind][!is.na(d_long[x,cond1_ind])] # long utr coverage in condition 1
      #utr_l2 <- d_long[x,cond2_ind][!is.na(d_long[x,cond2_ind])] # short utr coverage in condition 2
      #utr_s1 <- d_short[x,cond1_ind][!is.na(d_short[x,cond1_ind])] # long utr coverage in condition 1
      #utr_s2 <- d_short[x,cond2_ind][!is.na(d_short[x,cond2_ind])] # short utr coverage in condition 2
      pdui_1 <- mean(d_pdui[x,cond1_ind][!is.na(d_pdui[x,cond1_ind])])
      pdui_2 <- mean(d_pdui[x,cond2_ind][!is.na(d_pdui[x,cond2_ind])])
      #c(fisher.test(x = rbind(c(mean(utr_l1), mean(utr_s1)), c(mean(utr_l2), mean(utr_s2))))$p.value, pdui_2 - pdui_1)
      twobytwo <- rbind(c(utrl1_mean[x], utrs1_mean[x]), c(utrl2_mean[x], utrs2_mean[x]))
      c(fisher.test(x = twobytwo)$p.value, pdui_2 - pdui_1)
      
      
      #if(x == 'XM_039113041.1|Pou2f2|NC_051336.1|-'){
      # print(rbind(c(sum(utr_l1)/l1, sum(utr_s1)/l1), c(sum(utr_l2)/l1, sum(utr_s2)/l2)))
      #}
    }))
  }
  cat('finished\n')
  #print(dim(test))
  #print(dim(d_long))
  row.names(test) <- row.names(d_long)
  colnames(test) <- c('pval', 'mean.diff')
  test <- data.frame(test)
  test$pval[test$pval > 1] = 1
  test$padj <- p.adjust(test$pval, method = 'BH')
  test$fdr <- qvalue(test$pval)$qvalue
  
  if(impute){
    test$gene_short_names <- gene2region[row.names(test)]
  }else{
    test$gene_short_names <- sapply(row.names(test), FUN = function(x){strsplit(x,"\\|")[[1]][2]})
  }
  test$diff <- abs(test$mean.diff) > diff_thresh & test$fdr < alpha
  test$fit_value <- dapars[row.names(test),]$fit_value
  test$predicted_p_APA <- dapars_orig[row.names(test),]$Predicted_Proximal_APA
  test$loci <- dapars_orig[row.names(test),]$Loci
  test$strand = sapply(strsplit(row.names(test), '\\|'), FUN = function(x){x[4]})
  test$APA_dist = 0
  test[test$strand == '+',]$APA_dist <- abs(sapply(strsplit(test[test$strand == '+',]$loci, '-'), 
                                                   FUN = function(x){as.numeric(strsplit(x[1], ':')[[1]][2])}) - test[test$strand == '+',]$predicted_p_APA)-1
  test[test$strand == '-',]$APA_dist <- abs(sapply(strsplit(test[test$strand == '-',]$loci, '-'), 
                                                   FUN = function(x){as.numeric(x[2])}) - test[test$strand == '-',]$predicted_p_APA)-1
  cat('adjusting p values and combining p values\n')
  #test$APA_dist <- abs(sapply(strsplit(test$loci, '-'), FUN = function(x){as.numeric(x[2])}) - test$predicted_p_APA)-1
  
  gene_res <- data.frame(do.call(rbind, tapply(row.names(test), test$gene_short_names, function(x){
    df <- test[x,]
    min_pval <- min(df[,'pval'])
    min_pval_ind = which(df[,'pval'] == min_pval)
    min_pval_ind <- min_pval_ind[which(abs(df[min_pval_ind, 'mean.diff']) == max(abs(df[min_pval_ind, 'mean.diff'])))][1]
    if(!is.null(combine_p)){
      df[min_pval_ind,'pval'] <- metapod::combineParallelPValues(as.list(df[,'pval']), method = combine_p)$p.value
    }
    return(cbind(df[min_pval_ind, ], dapars[x[min_pval_ind],c(1,2,3)]))
  })))
  gene_res$padj <- p.adjust(gene_res$pval)
  gene_res$fdr <- qvalue(gene_res$pval)$qvalue
  gene_res$diff <- abs(gene_res$mean.diff) > diff_thresh & gene_res$fdr < alpha
  pdui <- dapars_orig[,grepl('PDUI', colnames(dapars_orig))]
  colnames(pdui) <- colnames(d_long)
  pdui <- pdui[row.names(d_long),]
  pdui_impute <- t(apply(pdui, 1, FUN = function(x){x[is.na(x)] = mean(x, na.rm =T); x}))
  
  return(list(deg= test, long = d_long, gene_res = gene_res, short = d_short, df = dapars_orig, pdui = pdui, pdui_imp = pdui_impute, gene_universe = all_genes))
}

dapars_impute <- function(file){
  dapars <- read.csv(file, sep = '\t', row.names = 1)
  d_pdui <- dapars[,grepl('PDUI', colnames(dapars))]
  dapars_out <- data.frame(Gene = row.names(dapars), dapars[,1:3], d_pdui)
  write.table(dapars_out, file = './temp.dp.tsv', sep = '\t', quote = F, row.names = F)
  d_pdui = scDaPars(raw_PDUI_file = './temp.dp.tsv',
                    out_dir = "apa/scDaPars_result",
                    filter_gene_thre = 0.2,
                    filter_cell_thre = 0.1)
  return(d_pdui)
}