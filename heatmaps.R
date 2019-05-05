rm(list =ls())
options(stringsAsFactors = F)

# the ExpressionMatrix is calculated by "RPKM.py" to confirm the RPKM value of each gene
a = read.csv('ExpressionMatrix.csv',head =T,row.names = 1)
a = round(a , digits = 0)
# tmp = a[1:10,1:10]
# meta =a[,1:5]
exprSet = data.frame(a[,2:ncol(a)])
group_list = c("0hours","0hours","0hours","2hours","2hours","2hours","4hours","4hours","4hours","6hours","6hours","6hours","8hours","8hours","8hours")
exprSet = data.matrix(exprSet)
# condition = factor(c("0hours","0hours","0hours","2hours","2hours","2hours","4hours","4hours","4hours"))


suppressMessages(library(DESeq2)) 
(colData <- data.frame(row.names=colnames(exprSet), group_list=group_list) )
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)
dds <- DESeq(dds)
# png("qc_dispersions.png", 1000, 1000, pointsize=20)
# plotDispEsts(dds, main="Dispersion plot")
 
res <- results(dds, contrast=c("group_list","4hours","0hours"))
resOrdered <- res[order(res$padj),]   
head(resOrdered)
DEG_2hours =as.data.frame(resOrdered)
# dev.off()  
DEG_2hours = na.omit(DEG_2hours)

nrDEG=DEG_2hours
## heatmap
library(pheatmap)

choose_gene=c('AT3G57260','AT1G04580', 'AT1G05530', 'AT1G05560', 'AT1G07240', 'AT1G15520', 'AT1G16540', 'AT1G20780', 'AT1G30100', 'AT1G30400', 'AT1G42550', 'AT1G50030', 'AT1G52340', 'AT1G52400', 'AT1G64670', 'AT1G65690', 'AT1G67080', 'AT1G69850', 'AT1G71960', 'AT1G78390', 'AT2G04240', 'AT2G18790', 'AT2G23250', 'AT2G23260', 'AT2G27150', 'AT2G29090', 'AT2G34660', 'AT3G14440', 'AT3G19270', 'AT3G21780', 'AT3G24220', 'AT3G26790', 'AT3G43600', 'AT4G18350', 'AT4G18780', 'AT4G19230', 'AT4G34138', 'AT5G08520', 'AT5G11240', 'AT5G20960', 'AT5G45340', 'AT5G52050', 'AT5G67030')
# the gene IDs were selected by "obtain_ABA_related_Gene.py" and copy here
choose_matrix=exprSet[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
choose_matrix = na.omit(choose_matrix)
choose_matrix = data.frame(choose_matrix)
write.csv(choose_matrix,"choose_matrix.csv")
# then run the python file "gene_to_symbol.py" with the output file to convert the tairID to symbolID
choose_matrix = read.csv("choose_matrix2.csv",,head =T,row.names = 1)
choose_matrix = subset(choose_matrix )
choose_matrix = data.matrix(choose_matrix)
pheatmap(choose_matrix,filename = 'DEG_top_1_50_0hourvs8hours_10_heatmap.png',
         cluster_cols = F)




## volcano plot
colnames(nrDEG)
plot(nrDEG$log2FoldChange,-log10(nrDEG$pvalue))

DEG=DEG_2hours


logFC_cutoff <- 1
# logFC_cutoff=1

DEG$change = as.factor(ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
                              ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
)
library(ggplot2)
g = ggplot(data=DEG, 
           aes(x=log2FoldChange, y=-log10(pvalue), 
               color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p_value") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)

print(g)

ggsave(g,filename = 'volcano_0h_vs_8h.png')
