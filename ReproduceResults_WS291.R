# Preparation With Mappings from WS291
require(GenomicFeatures)
worm_gtf <- makeTxDbFromGFF('./datasets/celegans_spike.gtf', 'gtf')
tx2gene_worm <- AnnotationDbi::select(worm_gtf, keys(worm_gtf, keytype = "TXNAME"), "GENEID", "TXNAME")
write.table(tx2gene_worm, './worm_tx2gene.tsv', sep = '\t', quote = F)


salmon_quant <- read_expression('./datasets/custom_param_quant/salmon/', 'salmon', tx2gene = './datasets/worm_tx2gene.tsv')
hisat_quant <- read_expression('./datasets/custom_param_quant/hisat/', 'hisat')
old_meta <- read.csv('./datasets/worms.meta.csv', header = T, row.names = 1)
WORM.GENES <- read.csv('./datasets/WormBase.WS291.genes.txt', header = T, sep = '\t', row.names = 1)


# Genes used for filtering and not analyzed: Mitochondrial genes, Vitellegenin genes, rRNA genes and Sperm specific genes
MITO_GENES <- row.names(subset(WORM.GENES, grepl('MtDNA',Chromosome.Coordinates)))



### obtain a list of intestine specifc genes for filtering
known_intestine_genes = c('WBGene00006926', 'WBGene00006925', 'WBGene00006927', 'WBGene00006928', 'WBGene00006929', 'WBGene00006930', 'WBGene00000781','WBGene00003473')
intestine_genes_2007 = read.csv('./datasets/Mcghee_et_al_2007.csv', header = T)
intestine_genes_2007 = subset(intestine_genes_2007, (Reliability == '++' | Reliability == '+++') & (In.situ == 'I' & GFP %in% c('_', 'I'))) 
intestine_genes = c('WBGene00006926', 'WBGene00006925', 'WBGene00006927', 'WBGene00006928', 
                    'WBGene00006929', 'WBGene00006930', 'WBGene00001500', 'WBGene00003090',
                    'WBGene00000218', 'WBGene00020662', 'WBGene00010769', 
                    'WBGene00003473', 'WBGene00000216', 'WBGene00000784', 'WBGene00000781',
                    'WBGene00000782', 'WBGene00000785', 'WBGene00000786', 'WBGene00003474')

genes_like_vit = t(cor(t(as.matrix(salmon_quant$gene$abundance[known_intestine_genes,])), t(as.matrix(salmon_quant$gene$abundance))))
genes_like_vit = genes_like_vit[rowMeans(salmon_quant$gene$abundance) > 0.5 | rowSums(salmon_quant$gene$abundance > 0.1) > 30,]
likely_intestine <- WORM.GENES[row.names(genes_like_vit[rowMaxs(genes_like_vit) > 0.5 & !is.na(genes_like_vit[,1]),]),]
likely_intestine$correlation = rowMax(genes_like_vit[row.names(genes_like_vit[rowMaxs(genes_like_vit) > 0.5,]), ])

all_possible_intestine_genes <- union(intestine_genes, row.names(subset(likely_intestine, Expr_pattern.Tissue == 'intestine' | correlation > 0.7)))


png('Intestine_correlation.png', width = 3, height = 3, units = 'in',res = 300)
par(mgp=c(2,1,0))
hist(rowMaxs(genes_like_vit), 100, col = 'black', main=NULL, xlab = 'Genes')
title('Max Correlation with Intestine Genes', line = 0, cex.main=0.7)
dev.off()
print(length(all_possible_intestine_genes))
print(length(all_possible_intestine_genes2))



#SPERM specific genes
micro_array_sperm <- read.csv('./datasets/Fig2a spermatogenesis-enriched gene set.txt', header = T,  sep = "\t")[,c(1,2,3,4, 5)]
micro_array_sperm$WormBase.Gene.ID <- WORM.GENES[match(micro_array_sperm$WormbaseID, WORM.GENES$Sequence.Name),]$WormBase.Gene.ID
micro_array_sperm$gene_name <- WORM.GENES[match(micro_array_sperm$WormbaseID, WORM.GENES$Sequence.Name),]$Public.Name


rnaseq_2014_sperm <- subset(read.csv('./datasets/2014_rnaseq.csv', header = T, row.names = 1), Gene.expression == 'Spermatogenic' )
rnaseq_2014_sperm[rnaseq_2014_sperm$log2.fold.change..of.normalized.reads. == 'Infinite',]$log2.fold.change..of.normalized.reads. = Inf
rnaseq_2014_sperm$log2FC = as.numeric(rnaseq_2014_sperm$log2.fold.change..of.normalized.reads.)
rnaseq_2014_sperm$gene_name = WORM.GENES[rnaseq_2014_sperm$WormBase.ID..WS240.,]$Public.Name
row.names(rnaseq_2014_sperm) <- rnaseq_2014_sperm$WormBase.ID..WS240.

CURATED_SPERM_GENES <- union(micro_array_sperm$WormBase.Gene.ID, row.names(subset(rnaseq_2014_sperm, log2.fold.change..of.normalized.reads. > 6)))


SPERM.SPECIFIC.GENES = intersect(row.names(salmon_quant$gene$abundance), subset(WORM.GENES, grepl('msp-', Public.Name) | grepl("^spe-", Public.Name) | grepl('ssp-', Public.Name))$WormBase.Gene.ID)
SPERM.SPECIFIC.GENES = setdiff(SPERM.SPECIFIC.GENES,c('WBGene00004975', 'WBGene00021771', 'WBGene00004960', 'WBGene00007732', 'WBGene00004959',
                                                      'WBGene00004969', 'WBGene00004976', 'WBGene00003433', 'WBGene00012211', 'WBGene00007117'))
MAJOR.SPERM.GENES = SPERM.SPECIFIC.GENES[salmon_quant$gene$abundance[SPERM.SPECIFIC.GENES, 'Ce497'] > 100]


genes_sperm = t(cor(t(as.matrix(salmon_quant$gene$abundance[SPERM.SPECIFIC.GENES,])), t(as.matrix(salmon_quant$gene$abundance))))
genes_sperm = genes_sperm[rowMeans(salmon_quant$gene$abundance) > 0.5  | (rowSums(salmon_quant$gene$abundance > 0.1) > 2 & rowSums(salmon_quant$gene$abundance > 0.1) < 50), ]
genes_sperm = genes_sperm[,!is.na(colMaxs(genes_sperm))]
genes_sperm_bycor = genes_sperm[rowMaxs(genes_sperm) > 0.9 & !is.na(rowSums(genes_sperm)),]
SPERM.ASSOCIATED.GENES = union(row.names(genes_sperm[rowMaxs(genes_sperm) > 0.99 & !is.na(rowSums(genes_sperm)),]), 
                               intersect(row.names(genes_sperm_bycor), CURATED_SPERM_GENES))


png('Sperm_correlation.png', width = 3, height = 3, units = 'in',res = 300)
par(mgp=c(2,1,0))
hist(rowMaxs(genes_sperm), 100, col = 'black', main=NULL, xlab = 'Genes')
title('Max Correlation with Sperm Genes', line = 0, cex.main=0.7)
dev.off()
print(length(SPERM.ASSOCIATED.GENES))
print(length(intersect(SPERM.ASSOCIATED.GENES, SPERM.ASSOCIATED.GENES2)))
#rRNA genes
rRNA_genes = intersect(row.names(salmon_quant$gene$abundance), subset(WORM.GENES, grepl('rrn-', Public.Name) | Biotype == 'rRNA_gene')$WormBase.Gene.ID)


spike_mito_msp_rRNA <- function(ct, cell_name){
  c('mito_rate'= as.numeric(sum(ct[MITO_GENES, cell_name])/colSums(ct[-c(1:92),])[cell_name]),
  'sperm_rate'=as.numeric(sum(ct[SPERM.SPECIFIC.GENES, cell_name])/colSums(ct[-c(1:92),])[cell_name]),
  'intestine_rate'=as.numeric(sum(ct[all_possible_intestine_genes, cell_name])/colSums(ct[-c(1:92),])[cell_name]),
  'spikein_rate'=as.numeric(sum(ct[c(1:92), cell_name])/colSums(ct)[cell_name]),
  'rRNA_rate'=as.numeric(sum(ct[rRNA_genes, cell_name])/colSums(ct[-c(1:92),])[cell_name]))
}



worms.mapping.info <- read.csv('./datasets/custom_param_quant/mapping.info', header = T, row.names = 1, sep = '\t')
worms.mapping.info$hisat.det.rate <- colSums(hisat_quant[,row.names(worms.mapping.info)] > 0)
worms.mapping.info$salmon.det.rate <- colSums(salmon_quant$gene$abundance[,row.names(worms.mapping.info)] > 0)
spike_mito_msp_rRNA_stat <- do.call(rbind, lapply(row.names(worms.mapping.info), FUN = function(x){spike_mito_msp_rRNA(salmon_quant$gene$abundance, x)}))
worms.mapping.info <- cbind(worms.mapping.info, spike_mito_msp_rRNA_stat)
worms.mapping.loc <- read.csv('./datasets/custom_param_quant/hisat_qorts_multiqc/multiqc_data/mqc_qorts_alignments_1.txt', header = T, row.names = 1, sep = '\t')
worms.mapping.info$hisat_uniq_cds <- worms.mapping.loc[row.names(worms.mapping.info), 'Unique.Gene..CDS']/worms.mapping.info$hisat_uniq
worms.mapping.info$hisat_uniq_utr <-  worms.mapping.loc[row.names(worms.mapping.info), 'Unique.Gene..UTR']/worms.mapping.info$hisat_uniq
worms.mapping.info$hisat_uniq_intron <-  worms.mapping.loc[row.names(worms.mapping.info), 'No.Gene..Intron']/worms.mapping.info$hisat_uniq
worms.mapping.loc <- read.csv('./datasets/custom_param_quant/hisat_qorts_multiqc/multiqc_data/mqc_qorts_splice_events_1.txt', header = T, row.names = 1, sep = '\t')
worms.mapping.info$reads_splice_known_loci <- rowSums(worms.mapping.loc[row.names(worms.mapping.info), c(1,2)])
worms.mapping.info$read.type <- 'se'; worms.mapping.info$Batch.ID <- 'Tintori_2016'
worms.mapping.info$orig_label <- worms.mapping.info$cellType <- worms.mapping.info$cellType1 <- 'P0'
worms.mapping.info$remarks <- worms.mapping.info$remarks1 <- worms.mapping.info$lane.ID <- worms.mapping.info$worm.id <- NA
worms.mapping.info$Analyze <- worms.mapping.info$Analyze1 <- T
worms.mapping.info[row.names(old_meta),]$orig_label <- old_meta$orig_label
worms.mapping.info[row.names(old_meta),]$worm.id <- old_meta$worm.id
worms.mapping.info[row.names(old_meta),]$cellType <- old_meta$alt.label
worms.mapping.info[row.names(old_meta),]$cellType1 <- old_meta$cellType2
worms.mapping.info[row.names(old_meta),]$read.type <- old_meta$read.type
worms.mapping.info[row.names(old_meta),]$lane.ID <- old_meta$lane.ID
worms.mapping.info[row.names(old_meta),]$Batch.ID <- old_meta$Batch.ID
worms.mapping.info[row.names(old_meta),]$remarks <- old_meta$Remarks
worms.mapping.info[row.names(old_meta),]$remarks1 <- old_meta$Remarks_A1
worms.mapping.info[row.names(old_meta),]$Analyze <- old_meta$Analyze
worms.mapping.info[row.names(old_meta),]$Analyze1 <- old_meta$Analyze1
worms.mapping.info$thresh_filt <- worms.mapping.info$orig_label != 'Unknown' & 
  worms.mapping.info$orig_label != '16-cell' & 
  worms.mapping.info$intestine_rate < 0.05 & 
  worms.mapping.info$sperm_rate < 0.01 & 
  worms.mapping.info$mito_rate < 0.05 & 
  worms.mapping.info$rRNA_rate < 0.05 & 
  worms.mapping.info$hisat_uniq_cds > 0.5 & 
  worms.mapping.info$hisat_uniq_rate > 70
worms.mapping.info$remarks1 <- str_replace_all(worms.mapping.info$remarks1, ',', ' ')
worms.mapping.info$remarks <- str_replace_all(worms.mapping.info$remarks, ',', ' ')
worms.mapping.info$tissue <- c('S1'='Mitotic Region and Transition Zone', 
                               'S2'='Pachytene Zone', 
                               'S3'='Diplotene Loop',
                               'S4'='Diakinesis', 
                               'F3'='-3 Oocyte', 
                               'F2'='-2 Oocyte', 
                               'F1'='-1 Oocyte', 
                               'P0'='Zygote',
                               'Unknown' = 'Zygote')[worms.mapping.info$cellType1]
worms.mapping.info$title <- paste('C. elegans gonad', worms.mapping.info$worm.id, worms.mapping.info$tissue, sep = ' ')
row.names(worms.mapping.info) <- str_replace_all(row.names(worms.mapping.info), 'SRR3170', '')
write.csv(worms.mapping.info, './worms.WS291.mapping.meta.csv', quote = F)

worms.meta <- read.csv('./worms.WS291.mapping.meta.csv', header = T, row.names = 1)




#Check Non-related gene expression content for a sample Function

# Samples to Remove
combined <- list()
combined$meta <- worms.mapping.info


combined$meta$group.id <- combined$meta$cellType1
combined$meta$unique.ids <- row.names(combined$meta)
combined$meta[combined$meta$group.id %in% c('Unknown', '1-cell'), ]$group.id <- 'P0'
combined$meta <- subset(combined$meta, cellType != 'tossed' & !grepl('cell', group.id))
combined$meta$cellType <- OLD_LABEL_2_NEW[combined$meta$group.id]
combined$meta$unique.ids <- str_replace(combined$meta$unique.ids, 'SRR3170', '')
row.names(combined$meta) <- combined$meta$unique.ids

combined$hisat <- hisat_quant[-c(1:92),]
colnames(combined$hisat) <- str_replace(colnames(combined$hisat), 'SRR3170', '')
combined$hisat <- combined$hisat[,row.names(combined$meta)]
combined$tpm <- t(t(salmon_quant$gene$abundance[-c(1:92),])/colSums(salmon_quant$gene$abundance[-c(1:92),]))*1e6
colnames(combined$tpm) <- str_replace(colnames(combined$tpm), 'SRR3170', '')
combined$tpm <- combined$tpm[,row.names(combined$meta)]
combined$salmon_ct <- salmon_quant$gene$counts[-c(1:92),]
colnames(combined$salmon_ct) <- str_replace(colnames(combined$salmon_ct), 'SRR3170', '')
combined$salmon_ct <- combined$salmon_ct[,row.names(combined$meta)]

#stress_genes

cengen_stress <- read.csv('./datasets/CenGen_stress_genes.csv', row.names = 3)
cengen_stress <- cengen_stress[setdiff(row.names(cengen_stress), c(MITO_GENES, rRNA_genes, SPERM.ASSOCIATED.GENES, all_possible_intestine_genes)),]
cengen_stress <- cengen_stress[intersect(row.names(cengen_stress), row.names(salmon_quant$gene$abundance[rowSums(salmon_quant$gene$abundance) > 0,])),]
stress <-  row.names(subset(WORM.GENES, grepl('hsp', Public.Name) | grepl('col', Public.Name) | grepl('clec', Public.Name)))
stress <- row.names(subset(cengen_stress, !is.na(Brunquell_source)))



#stress <- union(row.names(subset(cengen_stress, !is.na(Brunquell_source))),row.names(subset(WORM.GENES, grepl('hsp', Public.Name) | grepl('col', Public.Name) | grepl('clec', Public.Name))))
test_hmap_stress <- addTermsHeatmap(combined$tpm, stress, combined$meta,
                                    exprs_mat=NULL,
                                    dist_row = 'pearson',
                                    dist_col = 'euclidean',
                                    hmap_obj='sfd',
                                    org = 'celegans',
                                    k=5,
                                    universe = NULL, scale_minmax = F,
                                    plot = './test_hmap_stress.png',
                                    tpm = T)



#cengen_stress_cuticle <- intersect(row.names(combined$tpm), row.names(cengen_stress[grepl('col-',cengen_stress$Gene_name),]))
#cengen_stress_cuticle <- intersect(row.names(combined$tpm) , row.names(subset(WORM.GENES, grepl('hsp', Public.Name) | grepl('col', Public.Name) | grepl('clec', Public.Name))))
cengen_stress_cuticle <- test_hmap_stress$clust[[1]]


genes_stress = t(cor(t(as.matrix(combined$tpm)[cengen_stress_cuticle,]), t(as.matrix(combined$tpm))))
genes_stress = genes_stress[rowMeans(combined$tpm) > 0.5 | rowSums(combined$tpm > 0) > 2, ]
genes_stress = genes_stress[,!is.na(colMaxs(genes_stress))]
STRESS_GENES <- row.names(genes_stress[rowMaxs(genes_stress, na.rm = T) > 0.99,])

png('Stress_correlation.png', width = 3, height = 3, units = 'in',res = 300)
par(mgp=c(2,1,0))
hist(rowMaxs(genes_stress), 100, col = 'black', main=NULL, xlab = 'Genes')
title('Max Correlation with Stress Genes', line = 0, cex.main=0.7)
dev.off()

print(length(STRESS_GENES))
hist(rowMaxs(genes_stress), 100)
par(mfrow = c(1, 2))
hist(rowMeans(genes_sperm), 100, col = 'black', main = 'Pre-filtering Correlation with Sperm Genes', cex.main =0.7, xlab = 'Genes',)
hist(rowMeans(genes_sperm[setdiff(row.names(genes_sperm), SPERM.ASSOCIATED.GENES),]), 100, col = 'black', main = 'Post-Filter Correlation with Sperm Genes', cex.main = 0.7, xlab = 'Genes', ylab= '')


### COMPARISONS WITH PREVIOUS STUDIES
#Genes Found 
genes_by_type <- tapply(row.names(pData(ce.ct)), 
                        pData(ce.ct)$cellType, 
                        function(x){
                          mat = combined$hisat
                          row.names(mat)[rowMeans(mat) > 1 | rowSums(mat > 5) >= 3]
                        })
genes_by_type$F2F3 <- unique(c(genes_by_type$F2, genes_by_type$F3))
genes_by_type$OneCell <- genes_by_type$P0
genes_by_type$All_genes <- unique(do.call(c,genes_by_type[c('F1', 'F2', 'F3', 'S1', 'S2', 'S3', 'S4')]))
genes_by_type$Meiosis <- unique(c(genes_by_type$S3, genes_by_type$S2))
genes_by_type$Mitosis <- genes_by_type$S1
genes_by_type$Oocyte <- unique(c(genes_by_type$F3, genes_by_type$F2, genes_by_type$F1))

# Comparison with Tzur data
tzur_2018 <- read.csv('./datasets/Tzur_2018.csv', header = T, row.names = 1)
tzur_2018 <- subset(tzur_2018, WBGene != 'NAN')
row.names(tzur_2018) <- tzur_2018$WBGene
tzur_2018$WBGene <- NULL
tzur_2_su <- c('S1', 'S1', 'S2', 'S2', 'S2', 'S3', 'S3', 'S4', 'F2', 'F1', 'S1', 'S1', 'S2', 'S2', 'S2', 'S3', 'S3', 'S4', 'F2', 'F1')
names(tzur_2_su) <- colnames(tzur_2018)
tzur_genes_by_type <- tapply(colnames(tzur_2018), tzur_2_su, function(x){row.names(tzur_2018)[rowMeans(tzur_2018[,x]) > 1]})
tzur_genes_total <- unique(do.call(c, tzur_genes_by_type))
# Comparison with 
View(tzur_2018_mean)


# Microarray 2004
micro_array <- read.csv('./datasets/Fig2a oogenesis-enriched gene set.txt', header = T, row.names = 1, sep = "\t")
micro_array$WormBase.Gene.ID <- WORM.GENES[match(row.names(micro_array), WORM.GENES$Sequence.Name),]$WormBase.Gene.ID
micro_array$gene_name <- WORM.GENES[match(row.names(micro_array), WORM.GENES$Sequence.Name),]$Public.Name

marry_genes <- micro_array$WormBase.Gene.ID[!is.na(micro_array$WormBase.Gene.ID)]
micro_array <- read.csv('~/Downloads/dev00914-sup-data_s1/Fig1 IVa herm soma enriched genes.txt', header = T, sep = "\t")[,c(1,2,3,4)]
micro_array = data.frame(subset(micro_array, WormbaseID != '{{' & WormbaseID != ''), row.names = 1)
micro_array$WormBase.Gene.ID <- WORM.GENES[match(row.names(micro_array), WORM.GENES$Sequence.Name),]$WormBase.Gene.ID
micro_array$gene_name <- WORM.GENES[match(row.names(micro_array), WORM.GENES$Sequence.Name),]$Public.Name


#Gonad RNAseq 2014
rnaseq_2014 <- read.csv('./datasets/2014_rnaseq.csv', header = T, row.names = 1)
rnaseq_2014[rnaseq_2014$log2.fold.change..of.normalized.reads. == 'Infinite',]$log2.fold.change..of.normalized.reads. = Inf
rnaseq_2014$log2FC = as.numeric(rnaseq_2014$log2.fold.change..of.normalized.reads.)
rnaseq_2014
rnaseq_2014$gene_name = WORM.GENES[rnaseq_2014$WormBase.ID..WS240.,]$Public.Name
rnaseq_2014 <- subset(rnaseq_2014, Gene.expression == 'Oogenic' )

# OET RNAseq 2014
oet_2014 <- read.csv('./datasets/OET_2014.csv', header = T)
oet_2014 <- subset(oet_2014, !is.na(oet_2014$GENEWBID))[,c('GENEWBID','oocyte_rpkm','X1cell_rpkm')]
oet_2014 <- data.frame(do.call(cbind, lapply(colnames(oet_2014), function(n){
  tapply(oet_2014[,n], oet_2014$GENEWBID, function(x){if(class(x) == 'character'){unique(x)}else{sum(x)}})
})), row.names = 1)
oet_2014_genes <- list()
oet_2014_genes$Oocyte <- row.names(oet_2014)[oet_2014$X2 > 0.5]
oet_2014_genes$OneCell <- row.names(oet_2014)[oet_2014$X3 > 0.5]

# Liteseq 2018
liteseq_2018 <- read.csv('./datasets/liteseq_2018.csv', header = T, row.names = 1)
liteseq_2018_genes <- list()
liteseq_2018_genes$Mitosis <- row.names(liteseq_2018[liteseq_2018$Mitosis > 0,])
liteseq_2018_genes$Meiosis <- row.names(liteseq_2018[liteseq_2018$Meiosis > 0,])
liteseq_2018_genes$Oocyte <- row.names(liteseq_2018[liteseq_2018$Oocyte > 0,])

# Diag 2018
diag_2018 <- do.call(cbind, lapply(list.files("./datasets/Celseq_dataset/", pattern = ".tsv", full.names = T), FUN = function(x){
  f = read.csv(x, header = T, sep = "\t")
  counts = tapply(f$read_count, f$gene_id, sum)
}))
diag_2018 <- diag_2018[-c(1:86),]
diag_2018_meta <- t(sapply(list.files("./datasets/Celseq_dataset/", pattern = ".tsv", full.names = T), FUN = function(x){
  f = read.csv(x, header = T, sep = "\t")
  c(f$percent_distal_to_proximal[1],f$slice_index[1], f$barcode_sequence[1])
}))
diag_2018_meta <- data.frame(diag_2018_meta)
diag_2018_meta$replicate <- sapply(strsplit(row.names(diag_2018_meta), split = "_"), function(x){x[4]})
row.names(diag_2018_meta) <- paste(sapply(strsplit(row.names(diag_2018_meta), split = "_"), function(x){x[4]}), diag_2018_meta$X2, sep = "_")
diag_2018_meta$total_ct <- colSums(diag_2018)
colnames(diag_2018) <- row.names(diag_2018_meta)
diag_2018_meta$norm_pos <- do.call(c, tapply(row.names(diag_2018_meta), 
                                             diag_2018_meta$replicate, 
                                             FUN = function(x){
                                               min_x1 <- min(as.numeric(diag_2018_meta[x,'X1'])); max_x1 <- max(as.numeric(diag_2018_meta[x,'X1']))
                                               a <- as.numeric(diag_2018_meta[x,'X2'])*100/length(x)
                                               a <- (as.numeric(diag_2018_meta[x,'X1'])-min_x1)*100/(max_x1-min_x1)
                                               names(a) <- x;
                                               a}
))[paste(diag_2018_meta$replicate, row.names(diag_2018_meta), sep = ".")]
diag2su_by_pos <- function(x){
  if(x < 21)
  {cellType = 'S1'}
  else if(x < 52)
  {cellType = 'S2'}
  else if(x < 63)
  {cellType = 'S3'}
  else if(x < 84)
  {cellType = 'S4'}
  else if(x <(84+5))
  {cellType = 'F3'}
  else if(x < (84+5+5.2))
  {cellType = 'F2'}
  else if(x < (84+5+5.3+5.8))
  {cellType = 'F1'}
  return(cellType)
}


diag2su_by_slice <- function(x){
  if(x ==3 | x == 9){
    cellType = NA
  }
  else if(x < 3)
  {cellType = 'S1'}
  else if(x > 3 & x < 8)
  {cellType = 'S2'}
  else if(x == 8)
  {cellType = 'S3'}
  else if(x > 9 & x < 12)
  {cellType = 'S4'}
  else if(x < 13)
  {cellType = 'F3'}
  else if(x < 14)
  {cellType = 'F2'}
  else if(x < 16)
  {cellType = 'F1'}
  else{
    cellType = NA
  }
  return(cellType)
}


diag_2018_meta$su_label <- sapply(as.numeric(diag_2018_meta$X2), diag2su_by_slice)
diag_2018_meta$su_label_pos <- sapply(as.numeric(diag_2018_meta$norm_pos), diag2su_by_pos)
diag_2018_meta <- diag_2018_meta[diag_2018_meta$total_ct > 10000 & !is.na(diag_2018_meta$su_label),]
diag_2018 <- diag_2018[,row.names(diag_2018_meta)]
diag_genes_by_type <- tapply(colnames(diag_2018), diag_2018_meta$su_label, function(x){row.names(diag_2018)[rowMeans(diag_2018[,x]) > 1]})

comparePriorStudies <- function(dataset1, dataset2, prior_data_name){
  require(ggvenn)
  for(n in names(dataset2)){
    x <- list(a = dataset1[[n]],b = dataset2[[n]])
    names(x) <- c('Our data', prior_data_name)
    plt <- ggvenn(show_percentage = T, auto_scale = T,
                  x, 
                  fill_color = c("#0073C2FF", "#EFC000FF"),
                  stroke_size = 0.5, set_name_size = 4, 
    )+ggtitle(n)
    ggsave(
      filename = paste(prior_data_name, n,'venn_diagramm.png',sep="_"), plot = plt)
  }
}
comparePriorStudies(genes_by_type, diag_genes_by_type, 'Diag 2018')
comparePriorStudies(genes_by_type, tzur_genes_by_type, 'Tzur 2018')
comparePriorStudies(genes_by_type, liteseq_2018_genes, 'Liteseq 2018')
comparePriorStudies(genes_by_type, oet_2014_genes, 'OET 2014')

sapply(names(oet_2014_genes), FUN = function(x){c(length(intersect(genes_by_type[[x]], oet_2014_genes[[x]])), length(oet_2014_genes[[x]]))})
sapply(names(diag_genes_by_type), FUN = function(x){c(length(intersect(genes_by_type[[x]], diag_genes_by_type[[x]])), length(diag_genes_by_type[[x]]))})
sapply(names(tzur_genes_by_type), FUN = function(x){c(length(intersect(genes_by_type[[x]], tzur_genes_by_type[[x]])), length(tzur_genes_by_type[[x]]))})
sapply(names(liteseq_2018_genes), FUN = function(x){c(length(intersect(genes_by_type[[x]], liteseq_2018_genes[[x]])), length(liteseq_2018_genes[[x]]))})



length(intersect(genes_by_type$All_genes, marry_genes))/length(marry_genes)

length(intersect(genes_by_type$All_genes, rnaseq_2014$WormBase.ID..WS240.))/length(rnaseq_2014$Gene.ID)


# Correlation
lapply(unique(diag_2018_meta$su_label), function(x)
{
  dat1 = t(t(exprs(ce.ct))/sizeFactors(ce.ct))[,row.names(subset(pData(ce.ct), cellType == x))]
  dat1 <- vstExprs(ce.ct)[,row.names(subset(pData(ce.ct), cellType == x))]
  dat2 = diag_2018[, row.names(subset(diag_2018_meta, su_label == x))]
  dat1 = dat1[intersect(row.names(dat1), row.names(dat2)),]
  dat2 = dat2[intersect(row.names(dat1), row.names(dat2)),]
  print(dim(dat1))
  print(dim(dat2))
  mean(rowMax(cor(dat1, dat2)))
})


cor_tzur = do.call(cbind, lapply(unique(tzur_2_su), function(x)
{
  
  cors = do.call(c, lapply(unique(tzur_2_su), function(y){
    if(x == 'F2' & y == 'F2'){
      y=c('F2', 'F3')
    }
    dat1 = log(combined$hisat[,row.names(subset(pData(ce.ct), cellType %in% y))]+0.1)
    dat2 = log(tzur_2018[, names(tzur_2_su[tzur_2_su == x])]+0.1)
    dat1 = dat1[intersect(row.names(dat1), row.names(dat2)),]
    dat2 = dat2[intersect(row.names(dat1), row.names(dat2)),]
    #print(dim(dat1))
    #print(dim(dat2))
    mean(cor(dat1, dat2, method ='pearson'))
    #mean(rowMaxs(cor(dat1, dat2, method = 'pearson')))
  }))
  names(cors) = unique(tzur_2_su)
  cors
}))
colnames(cor_tzur) = unique(tzur_2_su)
plot_cor_heatmap(cor_tzur, 'Tzur et al, 2018', 'tzur_cor.png')


cor_diag = do.call(cbind, lapply(c('S1', 'S2', 'S3', 'S4', 'F3', 'F2', 'F1'), function(x)
{
  cors = do.call(c,lapply(c('S1', 'S2', 'S3', 'S4', 'F3', 'F2', 'F1'), function(y){
    dat1 = log(combined$hisat[,row.names(subset(pData(ce.ct), cellType == y))]+0.1)
    #dat1 <- vstExprs(ce.ct)[,row.names(subset(pData(ce.ct), cellType == y))]
    
    dat2 = diag_2018[, row.names(subset(diag_2018_meta, su_label == x))]
    dat1 = dat1[intersect(row.names(dat1), row.names(dat2)),]
    dat2 = log(dat2[intersect(row.names(dat1), row.names(dat2)),]+0.1)
    #dat2 = log(rowMeans(dat2)+0.1)
    #print(dim(dat1))
    #print(dim(dat2))
    #print(Hmisc::rcorr(as.matrix(dat1, dat2)))
    mean(cor(dat1, dat2))
    #mean(rowMaxs(cor(dat1, dat2, method = 'pearson')))
  }))
  names(cors) = c('S1', 'S2', 'S3', 'S4', 'F3', 'F2', 'F1')
  cors
  
}))
colnames(cor_diag) = c('S1', 'S2', 'S3', 'S4', 'F3', 'F2', 'F1')
plot_cor_heatmap(cor_diag, 'Diag et al, 2018', 'diag_cor.png')



#Run DEG with Monocle



filtered_genes <- c(MITO_GENES, all_possible_intestine_genes, SPERM.ASSOCIATED.GENES, STRESS_GENES,rRNA_genes )
filtered_hisat = combined$hisat[setdiff(row.names(combined$hisat), c(MITO_GENES, all_possible_intestine_genes, SPERM.ASSOCIATED.GENES, STRESS_GENES,rRNA_genes )),]
worms.genes.id <- WORM.GENES[row.names(filtered_hisat), c(1,2,4,7,16)]; worms.genes.id$gene_short_name <- worms.genes.id$Public.Name
ce.ct <- prepareMonocle(as.matrix(filtered_hisat)[,row.names(subset(combined$meta, thresh_filt))], meta = subset(combined$meta, thresh_filt), features = worms.genes.id, tpm2abs = F, monocle3 = F)
#pData(ce.ct )[pData(ce.ct )$cellType == 'Unknown',]$cellType <- 'P0'
pData(ce.ct )$experiment <- 'A'
pData(ce.ct )[pData(ce.ct )$Batch.ID == '150219HiSeq',]$experiment <- 'B'
pData(ce.ct )[c('Ce1x', 'Ce2x'),]$experiment <- 'B'
ce.ct <- ce.ct[,!row.names(pData(ce.ct)) %in% c('Ce533', 'Ce436', 'Ce435','Ce548','Ce460', 'Ce407', 'Ce415', 'Ce482')]
#ce.ct <- ce.ct[,which(pData(ce.ct )$experiment == 'A')]
sizeFactors(ce.ct) <- estimateSizeFactorsForMatrix(exprs(ce.ct))
#pData(ce.ct)$cellType <- pData(ce.ct)$orig_group
ce.ct <- dispersionSubset(ce.ct)



umap_config= list(n_neighbors = 10, min_dist = 0.5, metric='euclidean', random_state = 12345, transform_state = 12345)

sample_PCA(vstExprs(ce.ct), pData(ce.ct), umap = T, label = F, point_size = 2)+ggtitle('')
ggsave('./WS291_umap_110cells.png', width =4, height =4)


num_summary <- data.frame(num_of_Genes = colSums(exprs(ce.ct) > 0),
                          Batch = ce.ct$experiment,
                          cellType = c(ce.ct$cellType)) %>% group_by(cellType, Batch) %>% summarise(Num_of_Genes = mean(num_of_Genes), sd_num = sd(num_of_Genes), Batch = unique(Batch))

num_summary %>% mutate(cellType = factor(cellType, levels=c('S1','S2','S3','S4','F3','F2','F1','P0'))) %>%
  ggplot( aes(Batch, Num_of_Genes, fill = cellType)) + geom_bar(stat="identity", color="black",position=position_dodge()) + 
  geom_errorbar(aes(ymin = Num_of_Genes - sd_num, ymax = Num_of_Genes + sd_num), color = 'black',width = 0.2, position=position_dodge(.9))+
  theme_classic()+scale_fill_manual(values=DOT_COLOR)+ylab('Number of Expressed Genes')+theme(plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(0.4, 'cm'))
ggsave('number_gene_expressed.png', width = 3,height = 2.5)


num_summary %>% mutate(cellType = factor(cellType, levels=c('S1','S2','S3','S4','F3','F2','F1','P0'))) %>%
  ggplot( aes(x = cellType, group = Batch, y=Num_of_Genes, shape=Batch, linetype=Batch, color = Batch)) + geom_line(linewidth = 0.8, position=position_dodge(0.2)) + 
  geom_errorbar(aes(ymin = Num_of_Genes - sd_num, ymax = Num_of_Genes + sd_num), width = 0.6, position=position_dodge(0.2))+geom_point(position=position_dodge(0.2))+scale_linetype_manual(values=c("dashed", "solid"))+
  theme_classic()+scale_color_brewer(palette="Paired")+ylab('Number of Expressed Genes')+theme(plot.title = element_text(hjust = 0.5, size = 7), legend.key.size = unit(0.4, 'cm'))
ggsave('number_gene_expressed_line.png', width = 3,height = 2.5)

# Meta file for rMATs input
splicing_meta <- as.data.frame(pData(ce.ct))
row.names(splicing_meta)[!grepl('Ce', row.names(splicing_meta))] <- paste('SRR3170', row.names(splicing_meta)[!grepl('Ce', row.names(splicing_meta))], sep ='')
splicing_meta$reads_over_splice_rate <- splicing_meta$reads_splice_known_loci/splicing_meta$hisat_uniq
write.table(splicing_meta, 'WS291_filtered_splicing.meta.tsv', quote = F, sep ='\t')


test_01 <- ce.ct
sizeFactors(test_01) <- estimateSizeFactorsForMatrix(exprs(test_01))
test_01 <- estimateDispersions(test_01)

test01_fast <- list(S1_S2 = monocle_diff(test_01,subtype = c('S1', 'S2'),full_formula = '~cellType+num_genes_expressed+experiment', reduced_formula = '~num_genes_expressed+experiment', apeglm_method = 'nbinomC'),
                    S2_S3 = monocle_diff(test_01,subtype = c('S2', 'S3'),full_formula = '~cellType+num_genes_expressed+experiment', reduced_formula = '~num_genes_expressed+experiment', apeglm_method = 'nbinomC'),
                    S3_S4 = monocle_diff(test_01,subtype = c('S3', 'S4'),full_formula = '~cellType+num_genes_expressed+experiment', reduced_formula = '~num_genes_expressed+experiment', apeglm_method = 'nbinomC'),
                    S4_F3 = monocle_diff(test_01,subtype = c('S4', '-3'),full_formula = '~cellType+num_genes_expressed+experiment', reduced_formula = '~num_genes_expressed+experiment', apeglm_method = 'nbinomC'),
                    F3_F2 = monocle_diff(test_01,subtype = c('-3', '-2'),full_formula = '~cellType+num_genes_expressed+experiment',reduced_formula =  '~num_genes_expressed+experiment', apeglm_method = 'nbinomC'),
                    F2_F1 = monocle_diff(test_01,subtype = c('-2', '-1'),full_formula = '~cellType+num_genes_expressed+experiment',reduced_formula =  '~num_genes_expressed+experiment', apeglm_method = 'nbinomC'),
                    F1_P0 = monocle_diff(test_01,subtype = c('-1', 'P0'),full_formula = '~cellType+num_genes_expressed+experiment', reduced_formula = '~num_genes_expressed+experiment', apeglm_method = 'nbinomC')
)

sapply(test01_fast, function(x){c(nrow(subset(x, qval < 0.05 & Log2FC_shrunk > log2(1.5))), nrow(subset(x, qval < 0.05 & Log2FC_shrunk < -log2(1.5))))})


test01_fast_ora <- lapply(test01_fast, FUN = function(x){
  gse_list <- x$Log2FC_shrunk#effect*log(2))*(-log10(test01_S2_S3$S2_S3$pval))
  names(gse_list) <- row.names(x)
  gse_list <- sort(gse_list, T)
  gse_res <- enrich_CP(row.names(subset(x, qval < 0.05 & Log2FC_shrunk >= log2(1.5))), 'celegans', universe = row.names(x), logFC = gse_list, GSE = T, gsea_qval = 1.01)
})


aggregate_gse <- lapply(names(test01_fast_ora), FUN =function(nn){
  x <- test01_fast_ora[[nn]]
  df <- do.call(rbind, lapply(names(x), FUN = function(n){
    if(!grepl('gse', n) | grepl('MF', n) | grepl('CC', n)){return(NULL)}
    sub_df <- x[[n]]@result
    sub_df[,'Gene Set Type'] <- c('GO_BP_gse'='GO Biological Process', 'KEGG_gse'='KEGG', 'REACT_gse'='Reactome', 'WKP_gse'='Wikipathway')[n]
    sub_df$Description <- gsub(',', ' ', sub_df$Description)
    sub_df$leading_edge <- gsub(',', '|', sub_df$leading_edge)
    sub_df
  }))
  #write.csv(subset(df, qvalue < 0.05), paste(nn, '_gsea.csv', sep = ''), quote = F, row.names = F)
  #df$qvalue <- qvalue(df$pvalue)$qvalue
  #df$p.adjust <- p.adjust(df$pvalue)
  
  #obj <- x$GO_BP_gse
  #obj@result <- subset(df, qvalue < 0.05)
  #obj@params$pvalueCutoff <- 0.05
  #obj
})
names(aggregate_gse) <- names(test01_fast_ora)

test01_fast_gse <- lapply(test01_fast, FUN = function(x){
  gse_list <- x$Log2FC_shrunk#effect*log(2))*(-log10(test01_S2_S3$S2_S3$pval))
  names(gse_list) <- row.names(x)
  gse_list <- sort(gse_list, T)
  gse_res <- enrich_CP(row.names(subset(x, qval < 0.05 & Log2FC_shrunk >= log2(1.5))), 'celegans', universe = row.names(x), logFC = gse_list, GSE = T)
})

test01_fast_ora_up <- lapply(test01_fast, FUN = function(x){
  gse_list <- x$Log2FC_shrunk#effect*log(2))*(-log10(test01_S2_S3$S2_S3$pval))
  names(gse_list) <- row.names(x)
  gse_list <- sort(gse_list, T)
  gse_res <- enrich_CP(row.names(subset(x, qval < 0.05 & Log2FC_shrunk >= log2(1.5))), 'celegans', universe = deg_universe, logFC = gse_list, GSE = F)
})

test01_fast_ora_down <- lapply(test01_fast, FUN = function(x){
  gse_list <- x$Log2FC_shrunk#effect*log(2))*(-log10(test01_S2_S3$S2_S3$pval))
  names(gse_list) <- row.names(x)
  gse_list <- sort(gse_list, T)
  gse_res <- enrich_CP(row.names(subset(x, qval < 0.05 & Log2FC_shrunk <= -log2(1.5))), 'celegans', universe = deg_universe, logFC = gse_list, GSE = F)
})

#Gene Heatmap
# degs abs(Log2FC_shrunk) > log2(1.5), qval < 0.05; background all tested degs from all comparisons
# DESeq vst normalize -> pearson correlation distance -> wardD2 (WGCNA doesn't work as well, not a lot of enrichment) -> 20 clusters
# Enrichment GO biological process, qval < 0.05 or pvalue < 0.001

all_degs <- unique(do.call(c, lapply(test01_fast, FUN = function(x){row.names(subset(x, qval < 0.05 & abs(Log2FC_shrunk) >= log2(1.5)))})))

deg_universe <- unique(do.call(c, lapply(test01_fast, row.names)))

test_hmap <- addTermsHeatmap(exprs(ce.ct), all_degs, as.data.frame(pData(ce.ct)),
                             exprs_mat=NULL,
                             dist_row = 'pearson',
                             dist_col = 'euclidean',
                             hmap_obj='lkl',
                             org = 'celegans',
                             k=20,
                             universe = deg_universe, scale_minmax = T,
                             plot = './test_hmap_full.png',no_log = F,
                             tpm = F)
                 

genes_from_spermacae <- intersect(row.names(subset(test01_fast$F2_F1[row.names(subset(test01_fast$F1_P0, qval < 0.05 & Log2FC_shrunk < -7)),], qval < 0.05 & Log2FC_shrunk > 4)), row.names(ce.ct)[rowMeans(t(t(exprs(ce.ct))/sizeFactors(ce.ct))[, ce.ct$cellType == 'F1']) > 1000])
View(WORM.GENES[genes_from_spermacae,])




# F1_v_P0 histogram plot
logfc_in_other_comps <- data.frame(do.call(rbind, lapply(names(test01_fast)[1:6], function(n){
  x <- test01_fast[[n]]
  F1_P0 <-row.names(subset(test01_fast$F1_P0, Log2FC_shrunk < -log2(1.5) & qval < 0.05))
  logfc <- x[intersect(row.names(x), F1_P0),]$Log2FC_shrunk
  df <- data.frame(comparison=rep(paste(strsplit(n, '_')[[1]], collapse = '_v_'), length(logfc)), Log2FC=logfc)
  return(df)
})))
logfc_in_other_comps[logfc_in_other_comps$comparison != 'F2_v_F1', 'comparison'] <- 'Others'
ggplot(logfc_in_other_comps, mapping = aes(x=Log2FC, fill = comparison))+geom_density(position='identity', alpha = 0.6, adjust = 2)+theme_classic()+coord_trans(y = "log1p")+theme(legend.position = c(0.8, 0.8))+scale_fill_viridis_d()+geom_vline(xintercept = 4, color = 'red', linetype = 'dashed')
ggsave('./F1_P0_down_degs_in_F2_F1.png', height = 3, width = 3)
change_in_other_comps <- data.frame(do.call(rbind, lapply(names(test01_fast)[1:6], function(n){
  x <- test01_fast[[n]]
  F1_P0 <-row.names(subset(test01_fast$F1_P0, Log2FC_shrunk < -log2(1.5) & qval < 0.05))
  
  logfc_up <- sum(x[intersect(row.names(x), F1_P0),]$Log2FC_shrunk > log2(1.5))
  logfc_uup <- sum(x[intersect(row.names(x), F1_P0),]$Log2FC_shrunk > log2(4))
  logfc_down <- sum(x[intersect(row.names(x), F1_P0),]$Log2FC_shrunk < -log2(1.5))
  logfc_ddown <- sum(x[intersect(row.names(x), F1_P0),]$Log2FC_shrunk < -log2(4))
  df <- data.frame(comparison=rep(paste(strsplit(n, '_')[[1]], collapse = '_v_'), 4), Num_genes=c(logfc_up, logfc_uup, logfc_down, logfc_ddown), change=c('FC > 1.5', 'FC > 4', 'FC < 0.67', 'FC < 0.25'))
  return(df)
})))

# Write out gene lists to word or excel format
# Excel
write.csv(WORM.GENES[MITO_GENES, c(1,2)], 'mito_genes.csv', quote = F)
write.csv(WORM.GENES[known_intestine_genes, c(1,2)], 'intestine_specific_genes.csv', quote = F)
write.csv(WORM.GENES[intestine_genes, c(1,2)], 'curated_intestine_genes.csv', quote = F)
write.csv(WORM.GENES[SPERM.SPECIFIC.GENES, c(1,2)], 'sperm_specific.csv', quote = F)
write.csv(WORM.GENES[CURATED_SPERM_GENES, c(1,2)], 'curated_sperm_genes.csv', quote = F)
write.csv(WORM.GENES[SPERM.ASSOCIATED.GENES, c(1,2)], 'sperm_associated_genes.csv', quote = F)
write.csv(WORM.GENES[rRNA_genes, c(1,2)], 'rRNA_genes.csv', quote = F)
write.csv(WORM.GENES[all_possible_intestine_genes, c(1,2)], 'intestine_associated_genes.csv', quote = F)
write.csv(WORM.GENES[cengen_stress_cuticle, c(1,2)], 'curated_stress_genes.csv', quote = F)
write.csv(WORM.GENES[STRESS_GENES , c(1,2)], 'stress_associated_genes.csv', quote = F)

#Word paste
write(paste(WORM.GENES[SPERM.ASSOCIATED.GENES, 2], collapse = ', '), file = 'test_sperm.txt')
write(paste(WORM.GENES[all_possible_intestine_genes,2], collapse = ', '), file = 'test_intestine.txt')
write(paste(WORM.GENES[CURATED_SPERM_GENES,2], collapse = ', '), file = 'sperm_curate.txt')
write(paste(WORM.GENES[intestine_genes,2], collapse = ', '), file = 'intestine_curate.txt')
write(paste(WORM.GENES[known_intestine_genes,2], collapse = ', '), file = 'intestine_specific.txt')
write(paste(WORM.GENES[cengen_stress_cuticle,2], collapse = ', '), file = 'stress_curate.txt')
write(paste(WORM.GENES[STRESS_GENES,2], collapse = ', '), file = 'test_stress.txt')
write(paste(WORM.GENES[MITO_GENES,2], collapse = ', '), file = 'mito.txt')
write(paste(WORM.GENES[rRNA_genes,2], collapse = ', '), file = 'rRNA.txt')



# Writing out supplementary files to single excel tables
filter_genes_single_table <- rbind(data.frame('Gene Collection'=rep('Sperm Specific Genes', length(SPERM.SPECIFIC.GENES)), WORM.GENES[SPERM.SPECIFIC.GENES, c(1,2)]),
      data.frame('Gene Collection'=rep('Curated Sperm Genes', length(CURATED_SPERM_GENES)), WORM.GENES[CURATED_SPERM_GENES, c(1,2)]),
      data.frame('Gene Collection'=rep('SPERM ASSOCIATED GENES', length(SPERM.ASSOCIATED.GENES)), WORM.GENES[SPERM.ASSOCIATED.GENES, c(1,2)]),
      data.frame('Gene Collection'=rep('Intestine Specific Genes', length(known_intestine_genes)), WORM.GENES[known_intestine_genes, c(1,2)]),
      data.frame('Gene Collection'=rep('Curated Intestine Genes', length(intestine_genes)), WORM.GENES[intestine_genes, c(1,2)]),
      data.frame('Gene Collection'=rep('INTESTINE ASSOCIATED GENES', length(all_possible_intestine_genes)), WORM.GENES[all_possible_intestine_genes, c(1,2)]),
      data.frame('Gene Collection'=rep('Curated Stress Genes', length(cengen_stress_cuticle)), WORM.GENES[cengen_stress_cuticle, c(1,2)]),
      data.frame('Gene Collection'=rep('STRESS ASSOCIATED GENES', length(STRESS_GENES)), WORM.GENES[STRESS_GENES, c(1,2)]),
      data.frame('Gene Collection'=rep('rRNA Genes', length(rRNA_genes)), WORM.GENES[rRNA_genes, c(1,2)]),
      data.frame('Gene Collection'=rep('Mitochondrial Genes', length(MITO_GENES)), WORM.GENES[MITO_GENES, c(1,2)]))

write.csv(subset(filter_genes_single_table, !is.na(WormBase.Gene.ID)), file = 'Supplementary Table 2 v1.csv', quote = F, row.names = F)

DEG_single_table <- do.call(rbind, lapply(names(test01_fast), FUN = function(x){
  sub_df <- test01_fast[[x]][,c(5,6, 3,4,17)]
  df <- data.frame('Stages Tested'=rep(paste(strsplit(x, '_')[[1]], collapse = ' vs '), nrow(sub_df)), sub_df)
}))
write.csv(DEG_single_table, file = 'Supplementary Table S3.csv', quote = F, row.names = F)

hmap_single_table <- do.call(rbind, lapply(1:length(test_hmap$clust), FUN = function(x){
  sub_df <- test_hmap$clust[[x]]
  df <- data.frame('Wormbase.ID' = sub_df, 'Public.Name' = WORM.GENES[sub_df,2], 'Gene Cluster'=x )
}))
write.csv(hmap_single_table, file = 'Supplementary Table S4.csv', quote = F, row.names = F)


GSEA_single_table <- do.call(rbind, lapply(names(aggregate_gse), FUN = function(x){
  sub_df <- subset(aggregate_gse[[x]], qvalue < 0.05)
  df <- data.frame('Stages Tested'=rep(paste(strsplit(x, '_')[[1]], collapse = ' vs '), nrow(sub_df)), sub_df)
}))
write.csv(GSEA_single_table, file = 'Supplementary Table S5.csv', quote = F, row.names = F)


###Ribosomal protein expressions
rb_genes <- row.names(subset(WORM.GENES, !grepl(regex("trp|hrp|nrp"), Public.Name) & grepl(regex("rps|rpl|rla"), Public.Name)))
rb_genes <- intersect(row.names(exprs(ce.ct)), rb_genes)
rb_exprs <- data.frame(t(exprs(ce.ct)[rb_genes,])/sizeFactors(ce.ct))
rb_exprs$cellType <- factor(pData(ce.ct)$cellType, levels = c('S1', 'S2', 'S3', 'S4','F3', 'F2','F1','P0'))
rb_exprs$geneID <- subset(WORM.GENES, !grepl(regex("trp|hrp|nrp"), Public.Name) & grepl(regex("rps|rpl|rla"), Public.Name))[,'Public.Name']
rb_exprs %>% melt(id.vars = 'cellType') %>% ggplot(mapping = aes(fill = cellType, x = cellType, y = log10(value+1)))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_manual(values =DOT_COLOR)+ylab('log10 NormalizedExprs')+
  xlab('Cell Type')+
  theme(legend.position = 'none', axis.title=element_text(size=16))
ggsave('ribosome_protein_gene.png', width =6, height = 3)
### Draw Diagrams


png('test_2.png', width = 2.2, height = 4.8, units = 'in', res= 300, type = 'cairo')
plot_top_diff_genes(ce.ct, genes = c('ced-8'), plot_type = 'diagram', ncols = 1)
dev.off()




## APA
### UTR changes should be performed with samples that are processed similarly
## same batch, similar distribution of 3'UTR ends, similar levels of coverage
# extreme 3'prime bias can distort statistics without good batch correction methods
qorts_utr_stat <- read.csv('./datasets/custom_param_quant/multiqc_data/multiqc_qorts.txt', sep = '\t', row.names =1)
row.names(qorts_utr_stat) <- str_remove(row.names(qorts_utr_stat), 'SRR3170')
utr_file = './strict.txt'
test_dapars <- read.csv(utr_file, sep = '\t', row.names = 1); test_dapars <- test_dapars[,-c(1:3)];test_dapars <- test_dapars[,grepl('long', colnames(test_dapars))]; colnames(test_dapars) <- str_remove(colnames(test_dapars), '.depth_PDUI'); colnames(test_dapars) <- str_remove(colnames(test_dapars), '.depth_long_exp'); test_dapars_long <- test_dapars
test_dapars <- read.csv(utr_file, sep = '\t', row.names = 1); test_dapars <- test_dapars[,-c(1:3)];test_dapars <- test_dapars[,grepl('short', colnames(test_dapars))]; colnames(test_dapars) <- str_remove(colnames(test_dapars), '.depth_PDUI'); colnames(test_dapars) <- str_remove(colnames(test_dapars), '.depth_short_exp'); test_dapars_short <- test_dapars
test_dapars <- read.csv(utr_file, sep = '\t', row.names = 1); test_dapars <- test_dapars[,-c(1:3)];test_dapars <- test_dapars[,grepl('PDUI', colnames(test_dapars))]; colnames(test_dapars) <- str_remove(colnames(test_dapars), '.depth_PDUI'); test_dapars_pdui <- test_dapars
test_dapars_long <- test_dapars_long[rowSums(is.na(test_dapars_long)) == 0,]
test_dapars_short <- test_dapars_short[rowSums(is.na(test_dapars_short)) == 0,]
test_dapars_pdui <- test_dapars_pdui[rowSums(is.na(test_dapars_pdui)) == 0,]

colnames(test_dapars) <- str_remove(colnames(test_dapars), '^X')
test_dapars <- test_dapars[,!grepl(regex('X'), colnames(test_dapars))];
test_dapars <- test_dapars[sapply(row.names(test_dapars), FUN = function(x){strsplit(x,"\\|")[[1]][2]}) %in% row.names(exprs(ce.ct)[rowMeans(exprs(ce.ct)) > 10,]),]
test_dapars <- test_dapars[rowSums(is.na(test_dapars)) == 0,]

utr_cov <- data.frame(names(colSums(!is.na(test_dapars))), UTR_reads = qorts_utr_stat[colnames(test_dapars),'ReadPairs_UniqueGene_UTR'],
                      values=qorts_utr_stat[colnames(test_dapars),'ReadPairs_UniqueGene_UTR']/qorts_utr_stat[colnames(test_dapars),'ReadPairs_UniqueGene'], 
                      cellType = pData(ce.ct)[colnames(test_dapars), 'cellType'])
colnames(utr_cov)[1] <- 'samp'
ggplot(data = utr_cov, aes(x=values, fill = cellType))+geom_histogram(bins = 30)+facet_wrap(~cellType)
dapars_dist <- do.call(cbind, lapply(1:ncol(test_dapars), function(x){sapply(1:ncol(test_dapars), function(y){shared = intersect(which(!is.na(test_dapars[,x])), which(!is.na(test_dapars[,y]))); cor(test_dapars[shared,x], test_dapars[shared,y])})}));colnames(dapars_dist) <- colnames(test_dapars);row.names(dapars_dist) <- colnames(test_dapars)
utr_samples <- subset(utr_cov, values < 0.25 & values > 0.05 | cellType == 'P0')[,1]

sample_PCA(log(test_dapars_pdui*rowMeans(test_dapars_long+test_dapars_short)+1), pData(ce.ct)[utr_samples,], umap = T, color_by = 'cellType', labeling = F)+ggtitle('')+theme(legend.position = 'none')
ggsave('APA_strict_umap.png')
table(subset(utr_cov, values < 0.25 & values > 0.05 & !grepl(regex('Ce1|Ce2'), samp))$cellType)

utr_file = './strict.txt'
test_utr_F1_P0_all <- deg_utr(utr_file, exprs(ce.ct), compare = c('F1', 'P0'), meta = pData(ce.ct), combine_p = 'fisher')
test_utr_F1_P0_all$gene_res$SYMBOL <- WORM.GENES[test_utr_F1_P0_all$gene_res$gene_short_names, "Public.Name"]
test_utr_F1_P0 <- deg_utr(utr_file, exprs(ce.ct), compare = c('F1', 'P0'), meta = pData(ce.ct), combine_p = 'simes')
test_utr_F1_P0$gene_res$SYMBOL <- WORM.GENES[test_utr_F1_P0$gene_res$gene_short_names, "Public.Name"]
test_utr_F1_P0$gene_res$fdr[test_utr_F1_P0$gene_res$APA_dist < 40] <- 1

# seems like only F1 vs P0 produced significant changes in APA usage
utr_F1_P0_ora <- enrich_CP(subset(test_utr_F1_P0$gene_res, fdr < 0.05 & abs(mean.diff) > 0.05)$gene_short_names, universe = test_utr_F1_P0$deg$gene_short_names, organisms = 'celegans')
dot_test <- dotplot(clusterProfiler::simplify(utr_F1_P0_ora$GO_BP_ora))+theme_classic()
change_df <- do.call(rbind, lapply(lapply(strsplit(dot_test$data$geneID, '/'), function(x) row.names(WORM.GENES[match(x, WORM.GENES$Public.Name),])), function(x) c(sum(test_utr_F1_P0$gene_res[x, 'mean.diff'] >0.05), sum(test_utr_F1_P0$gene_res[x, 'mean.diff'] < -0.05))))
colnames(change_df) <- c('lengthened', 'shortened')
dot_test$data <- cbind(dot_test$data, change_df)


bar_test <- dot_test$data[,c('ID','Description','lengthened', 'shortened')] %>% melt(id.vars = c('ID', 'Description'), variable.name = 'change') %>% 
  ggplot( aes(fill=change, y=value, x=Description)) + coord_flip()+
  geom_bar(position="fill", stat="identity") +
  scale_fill_viridis(discrete = T) +
  ggtitle("") +
  theme_classic() +
  xlab("") + theme(axis.title.x = element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())+ylab('')

dap_summary <- data.frame(rbind(S1_v_S2 = c(nrow(subset(test_utr_S1_S2$gene_res, diff & mean.diff > 0)), nrow(subset(test_utr_S1_S2$gene_res, diff & mean.diff < 0))),
                                S2_v_S3 = c(nrow(subset(test_utr_S2_S3$gene_res, diff & mean.diff > 0)), nrow(subset(test_utr_S2_S3$gene_res, diff & mean.diff < 0))),
                                S3_v_S4 = c(nrow(subset(test_utr_S3_S4$gene_res, diff & mean.diff > 0)), nrow(subset(test_utr_S3_S4$gene_res, diff & mean.diff < 0))),
                                S4_v_F3 = c(nrow(subset(test_utr_S4_F3$gene_res, diff & mean.diff > 0)), nrow(subset(test_utr_S4_F3$gene_res, diff & mean.diff < 0))),
                                F3_v_F2 = c(nrow(subset(test_utr_F3_F2$gene_res, diff & mean.diff > 0)), nrow(subset(test_utr_F3_F2$gene_res, diff & mean.diff < 0))),
                                F2_v_F1 = c(nrow(subset(test_utr_F2_F1$gene_res, diff & mean.diff > 0)), nrow(subset(test_utr_F2_F1$gene_res, diff & mean.diff < 0))),
                                F1_v_P0 = c(nrow(subset(test_utr_F1_P0$gene_res, diff & mean.diff > 0)), nrow(subset(test_utr_F1_P0$gene_res, diff & mean.diff < 0)))))
colnames(dap_summary) <- c('lengthened', 'shortened')
dap_summary$Comparison <- factor(row.names(dap_summary), levels = row.names(dap_summary))

bar_test <- dap_summary %>% melt(id.vars = c('Comparison'), variable.name = 'change') %>% 
  ggplot( aes(fill=change, y=value, x=Comparison)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) +
  ggtitle("") +
  theme_classic() +
  xlab("Comparisons") +ylab('Number of Genes')

ggsave('./strict_DAP.png', width = 6, height = 3)

dot_test %>% insert_right(bar_test, width=.8) 

test_utr_S1_S2 <- deg_utr(utr_file, exprs(ce.ct), compare = c('S1', 'S2'), meta = pData(ce.ct), combine_p = 'simes')
test_utr_S1_S2$gene_res$SYMBOL <- WORM.GENES[test_utr_S1_S2$gene_res$gene_short_names, "Public.Name"]
test_utr_S2_S3 <- deg_utr(utr_file, exprs(ce.ct), compare = c('S2', 'S3'), meta = pData(ce.ct), combine_p = 'simes')
test_utr_S2_S3$gene_res$SYMBOL <- WORM.GENES[test_utr_S2_S3$gene_res$gene_short_names, "Public.Name"]
test_utr_S3_S4 <- deg_utr(utr_file, exprs(ce.ct), compare = c('S3', 'S4'), meta = pData(ce.ct), combine_p = 'simes')
test_utr_S3_S4$gene_res$SYMBOL <- WORM.GENES[test_utr_S3_S4$gene_res$gene_short_names, "Public.Name"]
test_utr_S4_F3 <- deg_utr(utr_file, exprs(ce.ct), compare = c('S4', 'F3'), meta = pData(ce.ct), combine_p = 'simes')
test_utr_S4_F3$gene_res$SYMBOL <- WORM.GENES[test_utr_S4_F3$gene_res$gene_short_names, "Public.Name"]
test_utr_F3_F2 <- deg_utr(utr_file, exprs(ce.ct), compare = c('F3', 'F2'), meta = pData(ce.ct), combine_p = 'simes')
test_utr_F3_F2$gene_res$SYMBOL <- WORM.GENES[test_utr_F3_F2$gene_res$gene_short_names, "Public.Name"]
test_utr_F2_F1 <- deg_utr(utr_file, exprs(ce.ct), compare = c('F2', 'F1'), meta = pData(ce.ct), combine_p = 'simes')
test_utr_F2_F1$gene_res$SYMBOL <- WORM.GENES[test_utr_F2_F1$gene_res$gene_short_names, "Public.Name"]

utr_F1_P0_ora <- enrich_CP(subset(test_utr_F1_P0$gene_res, fdr < 0.05 & abs(mean.diff) > 0.05)$gene_short_names, universe = test_utr_F1_P0$deg$gene_short_names, organisms = 'celegans')
utr_F1_P0_ora_all <- enrich_CP(subset(test_utr_F1_P0_all$gene_res, fdr < 0.05 & abs(mean.diff) > 0.05)$gene_short_names, universe = test_utr_F1_P0_all$deg$gene_short_names, organisms = 'celegans')

utr_F2_F1_ora <- enrich_CP(subset(test_utr_F2_F1$gene_res, fdr < 0.05 & abs(mean.diff) > 0.05)$gene_short_names, universe = test_utr_F2_F1$gene_res$gene_short_names, organisms = 'celegans')
utr_S1_S2_ora <- enrich_CP(subset(test_utr_S1_S2$gene_res, fdr < 0.05 & abs(mean.diff) > 0.05)$gene_short_names, universe = test_utr_S1_S2$gene_res$gene_short_names, organisms = 'celegans')
utr_S2_S3_ora <- enrich_CP(subset(test_utr_S2_S3$gene_res, fdr < 0.05 & abs(mean.diff) > 0.05)$gene_short_names, universe = test_utr_S2_S3$gene_res$gene_short_names, organisms = 'celegans')
utr_S3_S4_ora <- enrich_CP(subset(test_utr_S2_S3$gene_res, fdr < 0.05 & abs(mean.diff) > 0.05)$gene_short_names, universe = test_utr_S2_S3$gene_res$gene_short_names, organisms = 'celegans')
utr_S4_F3_ora <- enrich_CP(subset(test_utr_S4_F3$gene_res, fdr < 0.05 & abs(mean.diff) > 0.05)$gene_short_names, universe = test_utr_S4_F3$gene_res$gene_short_names, organisms = 'celegans')
utr_F3_F2_ora <- enrich_CP(subset(test_utr_F3_F2$gene_res, fdr < 0.05 & abs(mean.diff) > 0.05)$gene_short_names, universe = test_utr_F3_F2$gene_res$gene_short_names, organisms = 'celegans')



utr_F1_P0_ora <- enrich_CP(subset(test_utr_F1_P0$gene_res, fdr < 0.05 & abs(mean.diff) > 0.1)$gene_short_names, universe = test_utr_F1_P0$gene_res$gene_short_names, organisms = 'celegans')



# RMATS results


rmats_diff_res <- list()
rmats_diff_res[['S1_v_S2']] <- read_rmats_res('./datasets/custom_param_quant/rMATS/combined_out/S1_v_S2/')
rmats_diff_res[['S2_v_S3']] <- read_rmats_res('./datasets/custom_param_quant/rMATS/combined_out/S2_v_S3/')
rmats_diff_res[['S3_v_S4']] <- read_rmats_res('./datasets/custom_param_quant/rMATS/combined_out/S3_v_S4/')
rmats_diff_res[['S4_v_F3']] <- read_rmats_res('./datasets/custom_param_quant/rMATS/combined_out/S4_v_F3/')
rmats_diff_res[['F3_v_F2']] <- read_rmats_res('./datasets/custom_param_quant/rMATS/combined_out/F3_v_F2/')
rmats_diff_res[['F2_v_F1']] <- read_rmats_res('./datasets/custom_param_quant/rMATS/combined_out//F2_v_F1/')
rmats_diff_res[['F1_v_P0']] <- read_rmats_res('./datasets/custom_param_quant/rMATS/combined_out/F1_v_P0/')

rmats_summary <- do.call(cbind, lapply(rmats_diff_res, FUN = function(x){ sapply(x, FUN = function(m){length(m$GeneID)})}))
write.csv(rmats_summary, file = './rmats_num_genes.csv', quote = F)

rmats_ora_S1vS2 <- enrich_CP(do.call(c, lapply(rmats_diff_res$S1_v_S2, function(x){return(x$GeneID)})), 'celegans', universe = do.call(c, rmats_diff_res$S1_v_S2$universe))
rmats_ora_S1vS2$GO_BP_ora@result$qvalue <- round(rmats_ora_S1vS2$GO_BP_ora@result$qvalue, 3)
write.table(subset(rmats_ora_S1vS2$GO_BP_ora@result, qvalue < 0.05)[,c('Description', 'qvalue')], row.names = F,file = './results/diffsplice/rmats_S1_v_S2.bp.tsv', quote = F, sep ='\t')

rmats_ora_S2vS3 <- enrich_CP(do.call(c, lapply(rmats_diff_res$S2_v_S3, function(x){return(x$GeneID)})), 'celegans', universe = do.call(c, rmats_diff_res$S2_v_S3$universe))
rmats_ora_S2vS3$GO_BP_ora@result$qvalue <- round(rmats_ora_S2vS3$GO_BP_ora@result$qvalue, 3)
write.table(subset(rmats_ora_F1vP0$GO_BP_ora@result, qvalue < 0.05)[,c('Description', 'qvalue')], row.names = F,file = './results/diffsplice/rmats_F1_v_P0.bp.tsv', quote = F, sep ='\t')


rmats_ora_S3vS4 <- enrich_CP(do.call(c, lapply(rmats_diff_res$S3_v_S4, function(x){return(x$GeneID)})), 'celegans', universe = do.call(c, rmats_diff_res$S3_v_S4$universe))
rmats_ora_S3vS4$GO_BP_ora@result$qvalue <- round(rmats_ora_S3vS4$GO_BP_ora@result$qvalue, 3)
write.table(subset(rmats_ora_F1vP0$GO_BP_ora@result, qvalue < 0.05)[,c('Description', 'qvalue')], row.names = F,file = './results/diffsplice/rmats_F1_v_P0.bp.tsv', quote = F, sep ='\t')

rmats_ora_S4vF3 <- enrich_CP(do.call(c, lapply(rmats_diff_res$S4_v_F3, function(x){return(x$GeneID)})), 'celegans', universe = do.call(c, rmats_diff_res$S4_v_F3$universe))
rmats_ora_S4vF3$GO_BP_ora@result$qvalue <- round(rmats_ora_S4vF3$GO_BP_ora@result$qvalue, 3)
write.table(subset(rmats_ora_F1vP0$GO_BP_ora@result, qvalue < 0.05)[,c('Description', 'qvalue')], row.names = F,file = './results/diffsplice/rmats_F1_v_P0.bp.tsv', quote = F, sep ='\t')



rmats_ora_F1vP0 <- enrich_CP(do.call(c, lapply(rmats_diff_res$F1_v_P0, function(x){return(x$GeneID)})), 'celegans', universe = do.call(c, rmats_diff_res$F1_v_P0$universe))
rmats_ora_F1vP0$GO_BP_ora@result$qvalue <- round(rmats_ora_F1vP0$GO_BP_ora@result$qvalue, 3)
write.table(subset(rmats_ora_F1vP0$GO_BP_ora@result, qvalue < 0.05)[,c('Description', 'qvalue')], row.names = F,file = './results/diffsplice/rmats_F1_v_P0.bp.tsv', quote = F, sep ='\t')



rmats_ora_F1vP0 <- enrich_CP(do.call(c, lapply(rmats_diff_res$F1_v_P0, function(x){return(x$GeneID)})), 'celegans', universe = do.call(c, rmats_diff_res$F1_v_P0$universe))
rmats_ora_F1vP0$GO_BP_ora@result$qvalue <- round(rmats_ora_F1vP0$GO_BP_ora@result$qvalue, 3)
write.table(subset(rmats_ora_F1vP0$GO_BP_ora@result, qvalue < 0.05)[,c('Description', 'qvalue')], row.names = F,file = './results/diffsplice/rmats_F1_v_P0.bp.tsv', quote = F, sep ='\t')






