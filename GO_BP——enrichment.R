rm(list = ls())
library(topGO)
library(clusterProfiler)
library(Rgraphviz)
library(pathview)


library(org.At.tair.db)
# the outfile is from "select_differential.py" to find the what DEGs excatly are.
file= read.csv("outfile.csv")
gene_id = file$Unnamed..0
columns(org.At.tair.db)
enrich_go_bp= enrichGO(gene = gene_id,
                       OrgDb = org.At.tair.db,
                       keyType = "TAIR",
                       ont = "BP",
                       pvalueCutoff = 0.01,
                       qvalueCutoff = 0.05,
                       readable =T
                       )
write.csv(enrich_go_bp,"enrichment.csv")
barplot(enrich_go_bp)
png("enrich_GO_bp_dotplot.png",width = 1500,height = 3000,)
dotplot(enrich_go_bp)
dev.off()
pdf(file = "enrich_go_bp.pdf",width = 10,height = 15)
plotGOgraph(enrich_go_bp)
dev.off()



# the file is from David GO enrichment, download the files and and run them one by one
  library(ggplot2)
  ##读取刚才保存的富集分析结果文件
  goenrich <- read.csv("ABA_association_8hours.csv",header = T,encoding = "UTF-8")##header = T 第一行是表头，sep = "\t"表示以制表符（tab）分割
  ##利用levels来设置factor中的顺序，保证最后出图时按照我们之前排好的顺序排列
  
  x=goenrich$Fold.Enrichment
  y=factor(goenrich$Term,levels = goenrich$Term
           )
  ##先出一个框架
  p = ggplot(goenrich,aes(x,y))
  ##数据特征包括bubble大小来源为匹配到这个term上的基因数，颜色为qvalue，颜色变化为从低到高:"SpringGreen"到"DeepPink"而类型则是class，
  p1 = p + geom_point(aes(size=Pop.Hits,color=-1*log(FDR),shape=X.U.FEFF.Category,))+
    scale_color_gradient(low = "SpringGreen", high = "DeepPink")
  ##设置横纵坐标名字，标题，legend名字
  p2 = p1 + labs(color=expression(-log[10](Qvalue)),
                 size="Gene Number",
                 x="EnrichmentScore",
                 y="Go_term",
                 title="Go enrichment of ABA-related Genes")
  ##搞个主题，把边框画出来也可以通过ggplot Theme assistant修改
  p3 = p2 +theme_bw() 
  png("GO_BP_ABA-related_8hour.png",width = 1047,height = 900)
  p3
  dev.off()
  

# provide the data of tairID to geneID for python file
my_key <- keys(org.At.tair.db, keytype="ENTREZID")
# accession_list = c('AT3G14440', 'AT3G21780', 'AT3G24220', 'AT1G71960', 'AT1G67080', 'AT1G50030', 'AT1G30100', 'AT5G45340', 'AT4G19230', 'AT5G52050', 'AT2G34660', 'AT4G34138', 'AT5G20960', 'AT1G30400', 'AT5G11240', 'AT2G27150', 'AT1G69850', 'AT1G78390', 'AT2G23260', 'AT1G04580', 'AT3G43600')
my_col <- c('SYMBOL', 'TAIR','GO')
table = select(org.At.tair.db, keys=my_key, columns=my_col, keytype="ENTREZID")
write.csv(table,"table.csv")
# gene_id =  mapIds(org.At.tair.db,keys=accession_list,column="SYMBOL",keytype="TAIR",multiVals="first")
# print(gene_id)
