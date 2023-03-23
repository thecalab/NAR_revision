library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(UniprotR)
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(org.Mm.eg.db)
library(enrichplot) # great for making the standard GSEA enrichment plots
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(DT) #interactive and searchable tables of our GSEA results
library(ggrepel)
library(gplots)
library(limma)

#Input file was given as name of "SH_741_T_TEST" for type I ALT mESCs,
#and"SH_301_T_TEST" for type II ALT mESCs.

#For GSEA, start with gene name annotation: For Figure 2C and E
#https://www.uniprot.org/uploadlists/
#Annotate candidate list using website above and output was saved as :"typeI_ID_map" and "typeII_ID_map".

# GSEA with clusterprofiler
SH_741_T_TEST_2 <- merge(SH_741_T_TEST, typeI_ID_map, by.x="Accession", by.y="From")
length(which(is.na(SH_741_T_TEST_2$To)))
SH_741_T_TEST_2 <- filter(SH_741_T_TEST_2, t_test_Significant=="+") #only with significantly changed proteins

# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
TypeI<- dplyr::select(SH_741_T_TEST_2, To, Log_2_FC)
# construct a named vector
TypeI.gsea <- -(TypeI$Log_2_FC)
names(TypeI.gsea) <- as.character(toupper(TypeI$To))
TypeI.gsea <- sort(TypeI.gsea, decreasing = TRUE)
# run GSEA using the 'GSEA' function from clusterProfiler
c5.bp <- read.gmt("~/c5.go.bp.v7.5.1.symbols.gmt")

TypeI.res <- GSEA(TypeI.gsea, TERM2GENE=c5.bp, verbose=FALSE)
TypeI.df <- as_tibble(TypeI.res@result)


#bubble plot
TypeI.df.ordered <-TypeI.df[order(-TypeI.df$NES),] 
TypeI.df.ordered[1:20,] %>%
  arrange(NES) %>%
  mutate(ID=factor(ID, levels=ID)) %>%
  ggplot(aes(x = ID, y = NES)) + 
  geom_point(stat='identity', aes(size=setSize, col= -log10(p.adjust)))  +
  # ylim(2.0, 3.5) +
  coord_flip() + 
  theme_bw(base_size = 12) +
  theme(axis.title.y=element_blank(), legend.title = element_text(face = 2,color = 'black',
                                                                  
#TypeII analysis
SH_301_T_TEST_2 <- merge(SH_301_T_TEST, typeII_ID_map, by.x="Accession", by.y="From")
length(which(is.na(SH_741_T_TEST_2$To)))
SH_301_T_TEST_2 <- filter(SH_301_T_TEST_2, t_test_Significant=="+") #only with significantly changed proteins


# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
TypeII<- dplyr::select(SH_301_T_TEST_2, To, Log_2_FC)
# construct a named vector
TypeII.gsea <- -(TypeII$Log_2_FC)
names(TypeII.gsea) <- as.character(toupper(TypeII$To))
TypeII.gsea <- sort(TypeII.gsea, decreasing = TRUE)
# run GSEA using the 'GSEA' function from clusterProfiler
c5.bp <- read.gmt("~/c5.go.bp.v7.5.1.symbols.gmt")

#bubble plot
TypeII.res <- GSEA(TypeII.gsea, TERM2GENE=c5, verbose=FALSE)
TypeII.df <- as_tibble(TypeII.res@result)
TypeII.df.ordered <-TypeII.df[order(-TypeII.df$NES),] 
TypeII.df.ordered[1:20,] %>%
  arrange(NES) %>%
  mutate(ID=factor(ID, levels=ID)) %>%
  ggplot(aes(x = ID, y = NES)) + 
  geom_point(stat='identity', aes(size=setSize, col= -log10(p.adjust)))  +
  # ylim(2.0, 3.5) +
  coord_flip() + 
  theme_bw(base_size = 12) +
  theme(axis.title.y=element_blank(), legend.title = element_text(face = 2,color = 'black',
                                                                  size = 10))


## Volcano plot: Figure B and D
SH_741_T_TEST_3 <- SH_741_T_TEST_2 
SH_741_T_TEST_3 <- SH_741_T_TEST_3 %>% mutate(color = "")%>% 
  mutate(threshold =-log10(p.adjust) > 1 & abs(Log_2_FC) >= 0.58)

SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Hmgn1'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Srcap'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'znhit1'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Atrx'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Daxx'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Dnmt1'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Chd2'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'H3-5'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Kdm4b'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Kdm6b'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Vps72'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Dmap1'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Gatad2b'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Mta1'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Mbd2'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Rsf1'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Baz1b'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Parp1'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Zbtb1'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Kdm5b'),"color"] <- "check"

SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Smarcd3'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Satb1'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Tet2'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Brwd1'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Hmga1'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Dpf1'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'H1f0'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Hist2h2ab'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Kat2b'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Ing4'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Vrk1'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'H1fx'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Rcor1'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Sirt7'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Baz2b'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Kat6b'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Kdm3a'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Baz2a'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Baz2b'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Zbtb7a'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Actr6'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Cecr2'),"color"] <- "check"
SH_741_T_TEST_3[which(SH_741_T_TEST_3$To == 'Smarcc2'),"color"] <- "check"

ggplot(SH_741_T_TEST_3) +
  geom_point(aes(x = -(Log_2_FC), y = -log10(p.adjust),colour = threshold),  size = 1.5, alpha = 1) +
  ggtitle("Type I") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values = c("darkgrey", "#56B4E9")) +
  theme_bw(base_size = 14) + 
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_text_repel(data=subset(SH_741_T_TEST_3, color=="check"),aes(x=-(Log_2_FC), -log10(p.adjust), label=To),
                  size = 3, fontface="bold",  
                  box.padding = unit(1, "lines"),
                  point.padding = unit(0.3, "lines"),
                  max.overlaps = Inf)+
  lims(x=c(-3,3)) +
  lims(y=c(0,2.2)) +
  geom_hline(yintercept = 1, colour="#990000", linetype="dashed") + geom_vline(xintercept = 0.58, colour="#990000", linetype="dashed") + geom_vline(xintercept = -0.58, colour="#990000", linetype="dashed")


SH_301_T_TEST_3 <- SH_301_T_TEST_2

SH_301_T_TEST_3 <- SH_301_T_TEST_3 %>% mutate(color = "")%>%   mutate(threshold =-log10(p.adjust) > 1 & abs(Log_2_FC) >= 0.58)

SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Brca1'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Blm'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Atm'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Top3a'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Fanca'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Fancm'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Rad50'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Mre11a'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Nbs1'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'H2ax'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Ctip'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Topbp1'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Atr'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Rpa1'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Brip1'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Bard1'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Palb2'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Brca2'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Rad51'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Mus81'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Eme1'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Fancm'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Fancb'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Fancg'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Fance'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Fancd2'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Fanci'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Fance'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Fance4'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Fance2'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Faap100'),"color"] <- "check"
SH_301_T_TEST_3[which(SH_301_T_TEST_3$To == 'Rb1'),"color"] <- "check"

ggplot(SH_301_T_TEST_3) +
  geom_point(aes(x = -(Log_2_FC), y = -log10(p.adjust),colour = threshold), size = 1.5, alpha = 1) +
  ggtitle("Type II") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values = c("darkgrey", "#56B4E9")) +
  theme_bw(base_size = 14) + 
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  geom_text_repel(data=subset(SH_301_T_TEST_3, color=="check"),aes(x=-(Log_2_FC), -log10(p.adjust), label=To),
                  size = 3, fontface="bold",  
                  box.padding = unit(1, "lines"),
                  point.padding = unit(0.3, "lines"),
                  max.overlaps = Inf )+
  lims(x=c(-3,3)) +
  lims(y=c(0,2.2)) +
  geom_hline(yintercept = 1, colour="#990000", linetype="dashed") + geom_vline(xintercept = 0.58, colour="#990000", linetype="dashed") + geom_vline(xintercept = -0.58, colour="#990000", linetype="dashed")


#GO anaylsis for Figure 2F-I
typeI_sig <- SH_741_T_TEST %>% filter(t_test_Significant=="+")
typeI_sig_inc <- typeI_sig %>% filter(Log_2_FC <0)
typeI_sig_dec <- typeI_sig %>% filter(Log_2_FC >0)


typeII_sig <- SH_301_T_TEST %>% filter(t_test_Significant=="+")
typeII_sig_inc <- typeII_sig %>% filter(Log_2_FC <0)
typeII_sig_dec <- typeII_sig %>% filter(Log_2_FC >0)


both_inc <- merge(typeI_sig_inc, typeII_sig_inc, by="Accession")
both_dec <- merge(typeI_sig_dec, typeII_sig_dec, by="Accession")
typeI_inc_typeII_dec <- merge(typeI_sig_inc, typeII_sig_dec, by="Accession")
typeI_dec_typeII_inc <- merge(typeI_sig_dec, typeII_sig_inc, by="Accession")

both_inc_name <- merge(both_inc, typeI_ID_map, by.x="Accession", by.y="From")
both_inc_name <- both_inc_name$To

both_dec_name <- merge(both_dec, typeI_ID_map, by.x="Accession", by.y="From")
both_dec_name <- both_dec_name$To

both_ID <- rbind(typeI_ID_map,typeII_ID_map)
both_ID <- distinct(both_ID)

typeI_inc_typeII_dec_name <- merge(typeI_inc_typeII_dec, both_ID, by.x="Accession", by.y="From")
typeI_inc_typeII_dec_name <- typeI_inc_typeII_dec_name$To

typeI_dec_typeII_inc_name <- merge(typeI_dec_typeII_inc, both_ID, by.x="Accession", by.y="From")
typeI_dec_typeII_inc_name <- typeI_dec_typeII_inc_name$To


typeI_sig_inc_name <- merge(typeI_sig_inc, both_ID, by.x="Accession", by.y="From")
typeI_sig_inc_name <- typeI_sig_inc_name$To
typeI_sig_dec_name <- merge(typeI_sig_dec, both_ID, by.x="Accession", by.y="From")
typeI_sig_dec_name <- typeI_sig_dec_name$To

only_typeI_inc <- setdiff(typeI_sig_inc_name,both_inc_name)
#only_typeI_inc <- setdiff(only_typeI_inc,typeI_inc_typeII_dec_name)
only_typeI_dec <- setdiff(typeI_sig_dec_name,both_dec_name)
#only_typeI_dec <- setdiff(only_typeI_dec,typeI_dec_typeII_inc_name)

typeII_sig_inc_name <- merge(typeII_sig_inc, both_ID, by.x="Accession", by.y="From")
typeII_sig_inc_name <- typeII_sig_inc_name$To
typeII_sig_dec_name <- merge(typeII_sig_dec, both_ID, by.x="Accession", by.y="From")
typeII_sig_dec_name <- typeII_sig_dec_name$To

only_typeII_inc <- setdiff(typeII_sig_inc_name,both_inc_name)
#only_typeII_inc <- setdiff(only_typeII_inc,typeI_dec_typeII_inc_name)
only_typeII_dec <- setdiff(typeII_sig_dec_name,both_dec_name)
#only_typeII_dec <- setdiff(only_typeII_dec,typeI_inc_typeII_dec_name)

#GO analysis--------------
ego_both_inc <- enrichGO(gene = both_inc_name, 
                         universe = all_genes,
                         keyType = "SYMBOL",
                         OrgDb = org.Mm.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = F)


ego_both_dec <- enrichGO(gene = both_dec_name, 
                         universe = all_genes,
                         keyType = "SYMBOL",
                         OrgDb = org.Mm.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = F)


ego_I_inc_II_dec <- enrichGO(gene = typeI_inc_typeII_dec_name, 
                             universe = all_genes,
                             keyType = "SYMBOL",
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = F)


ego_I_dec_II_inc <- enrichGO(gene = typeI_dec_typeII_inc_name, 
                             universe = all_genes,
                             keyType = "SYMBOL",
                             OrgDb = org.Mm.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = F)

dotplot(ego_both_inc, showCategory=20)
dotplot(ego_both_dec, showCategory=20)
dotplot(ego_I_inc_II_dec, showCategory=20)
dotplot(ego_I_dec_II_inc, showCategory=20)


ego_only_typeI_inc <- enrichGO(gene = only_typeI_inc, 
                               universe = all_genes,
                               keyType = "SYMBOL",
                               OrgDb = org.Mm.eg.db, 
                               ont = "BP", 
                               pAdjustMethod = "BH", 
                               qvalueCutoff = 0.05, 
                               readable = F)
ego_only_typeI_dec <- enrichGO(gene = only_typeI_dec, 
                               universe = all_genes,
                               keyType = "SYMBOL",
                               OrgDb = org.Mm.eg.db, 
                               ont = "BP", 
                               pAdjustMethod = "BH", 
                               qvalueCutoff = 0.05, 
                               readable = F)
ego_only_typeII_inc <- enrichGO(gene = only_typeII_inc, 
                                universe = all_genes,
                                keyType = "SYMBOL",
                                OrgDb = org.Mm.eg.db, 
                                ont = "BP", 
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05, 
                                readable = F)
ego_only_typeII_dec <- enrichGO(gene = only_typeII_dec, 
                                universe = all_genes,
                                keyType = "SYMBOL",
                                OrgDb = org.Mm.eg.db, 
                                ont = "BP", 
                                pAdjustMethod = "BH", 
                                qvalueCutoff = 0.05, 
                                readable = F)



dotplot(ego_only_typeI_inc, showCategory=20)
dotplot(ego_only_typeI_dec, showCategory=20)
dotplot(ego_only_typeII_inc, showCategory=20)
dotplot(ego_only_typeII_dec, showCategory=20)