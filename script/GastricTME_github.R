library(ggplot2);library(data.table);
library(Seurat);library(dplyr);library(tidyr);library(ggpubr)
library(pheatmap);
library('ggsignif');
#### Fig 1D ####
load(file='meta.Rdata')
aa<-table(meta$sample);
sample<-names(aa)[aa>= 200]
data=meta[meta$sample %in% sample,c('tissue.level.8','celltype.big.L1')];
data$celltype.big.L1.v2<-factor(data$celltype.big.L1,levels = c('B cell','Plasma','CD4T','CD8T','T cell','NK','pDCs','Myeloid','Mast','Stromal'))
data<-data[,c(2,4)];
colnames(data)<-c('sample','celltype.big')
data1=as.data.frame(table(data))
data2 = data1 %>% group_by(sample) %>% summarise(count=sum(Freq))
data3<-inner_join(data1,data2,by='sample')
data3$percentage<- data3$Freq*100/data3$count;
data3$sample<-factor(as.character(data3$sample),levels = c('NGT','NAT','CAG','IM','GAC_Primary','Metastasis','GAC_PBMC','Healthy_PBMC'))
data3$tissue.level<-as.numeric(data3$sample)
data3<-data3[data3$sample != 'NGT',]
color.cell.type<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ae017e','#f768a1')
names(color.cell.type)<-c('B cell','granulocytes','Mast','Myeloid','NK','pDCs','Plasma','Proliferative plasma cell','Stromal','T cell','CD4T','CD8T')
ggplot(data3, aes(x=tissue.level, y=percentage, fill=celltype.big)) + 
  geom_area(aes(fill = celltype.big),stat="identity")+
  scale_fill_manual(values = color.cell.type)

#### Fig 1E ####
load(file='meta.Rdata')

meta<-meta[,c('sample','tissue.level.8','celltype.big.L1')]
aa<-table(meta$sample);
sample<-names(aa)[aa>= 200]
data=meta[meta$sample %in% sample,c('sample','celltype.big.L1')];
colnames(data)<-c('sample','celltype.big')
data1=as.data.frame(table(data))
data2 = data1 %>% group_by(sample) %>% summarise(count=sum(Freq))
data3<-inner_join(data1,data2,by='sample')
data3$percentage<- data3$Freq*100/data3$count;
data3$tissue.level<-meta[match(data3$sample,meta$sample),'tissue.level.8'];
data3<-data3[data3$tissue.level %in% c('NAT','CAG','IM','GAC_Primary','Metastasis'),]
data3$tissue.level<-factor(data3$tissue.level,levels = c('NAT','CAG','IM','GAC_Primary','Metastasis'));
color.cell.type<-c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf');
names(color.cell.type)<-c('CAG','GAC_PBMC','GAC_Primary','Healthy_PBMC','IM','Metastasis','NAT','NGT')
celltype<-c('Plasma','Mast','Myeloid')

gp.list<-list();
for(i in 1:length(celltype)){
  dat<-data3[data3$celltype.big %in% as.character(celltype[i]),]
  dat<-dat[dat$count>200,]
  plotx<-ggboxplot(dat, "tissue.level", "percentage",
                   color = 'tissue.level', 
                   palette =color.cell.type,
                   add = "jitter",
                   title=celltype[i]) +
    stat_compare_means()+
    theme(axis.line = element_line(colour = "black"),
          axis.text.x=element_text(angle = 90));
  gp.list[[length(gp.list) + 1]] = plotx
}
do.call(ggarrange,c(gp.list,ncol = 1,nrow = 1)) -> combined.gp
pdf('Fig1E.pdf',width =5,height=5)
print(combined.gp)
dev.off()
#### Fig 2D #####
# heatmap
load(file='/rsrch3/scratch/genomic_med/rwang12/Gastric.scs/merge/TME/meta.77392.Rdata')
meta<-meta[meta$celltype.big == 'CD8',c('tissue.level.8','celltype.harmony.2')];
data<-table(meta[,c('tissue.level.8','celltype.harmony.2')]);
data<-t(data)
res<-chisq.test(data)
ratio<-res$observed/res$expected;
ratio<-ratio[c('CD8_C1_NR4A2_SPRY1','CD8_C4_NR4A1_FOS_JUN','CD8_C8_TIM3_LAG3_TIGIT_PD-1',
               'CD8_C5_TRAV1-2_SLC4A10_KLRB1','CD8_C6_GNLY_KLRD1','CD8_C2_HLA-II_GZMH_LGALS1',
               'CD8_C0_CMC1_GZMK_HLA-II','CD8_C7_ISG15_IFI6_IFIT3/1','CD8_C3_CCR7_SELL_LEF1_TCF7','CD8_C9_FCMR_TIM3_CD27_GZMK'),
             c('NAT','CAG','IM','GAC_Primary','Metastasis')];

pheatmap(ratio,
         color = c(colorRampPalette(c('#313695','white'))(36),
                   colorRampPalette(c('white','#ff7f00'))(60)),
         cluster_rows = F,
         cluster_cols = F);

# barplot by tissuelevel
load(file='/rsrch3/scratch/genomic_med/rwang12/Gastric.scs/merge/TME/meta.77392.Rdata')
meta<-meta[meta$celltype.big == 'CD8',c('tissue.level.8','celltype.harmony.2')];
data<-as.data.frame(table(meta[,c('tissue.level.8','celltype.harmony.2')]));

color.cell.type<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')
names(color.cell.type)<-c('CD8_C0_CMC1_GZMK_HLA-II','CD8_C1_NR4A2_SPRY1','CD8_C2_HLA-II_GZMH_LGALS1','CD8_C3_CCR7_SELL_LEF1_TCF7','CD8_C4_NR4A1_FOS_JUN','CD8_C5_TRAV1-2_SLC4A10_KLRB1','CD8_C6_GNLY_KLRD1','CD8_C7_ISG15_IFI6_IFIT3/1','CD8_C8_TIM3_LAG3_TIGIT_PD-1','CD8_C9_FCMR_TIM3_CD27_GZMK')

data$celltype.harmony.2<-factor(as.character(data$celltype.harmony.2),
                                 levels = c('CD8_C1_NR4A2_SPRY1','CD8_C4_NR4A1_FOS_JUN',
                                            'CD8_C8_TIM3_LAG3_TIGIT_PD-1',
                                            'CD8_C5_TRAV1-2_SLC4A10_KLRB1',
                                            'CD8_C6_GNLY_KLRD1',
                                            'CD8_C2_HLA-II_GZMH_LGALS1',
                                            'CD8_C0_CMC1_GZMK_HLA-II',
                                            'CD8_C7_ISG15_IFI6_IFIT3/1',
                                            'CD8_C3_CCR7_SELL_LEF1_TCF7',
                                            'CD8_C9_FCMR_TIM3_CD27_GZMK'));
data<-data[data$tissue.level.8 %in% c('NAT','CAG','IM','GAC_Primary','Metastasis'),]
data$tissue.level.8<-factor(data$tissue.level.8,levels=c('NAT','CAG','IM','GAC_Primary','Metastasis'))
ggplot(data, aes(x=tissue.level.8,y=Freq)) + 
  geom_bar(aes(fill = celltype.harmony.2),stat="identity",position="fill") + 
  scale_fill_manual(values = color.cell.type)+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90));
# barplot by celltype
color.tissue.type<-c('#C0AB34','#0D4C7A','#C25757','#72ADE0','#F8A729','#9A253E','#EDAEAE','#5AAA46');
names(color.tissue.type)<-c('CAG','GAC_PBMC','GAC_Primary','Healthy_PBMC','IM','Metastasis','NAT','NGT')
ggplot(data, aes(x=celltype.harmony.2,y=Freq)) + 
  geom_bar(aes(fill = tissue.level.8),stat="identity",position="fill") + 
  scale_fill_manual(values = color.tissue.type)+
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90));

#### Fig 3E ####
library('BiocNeighbors')
obj<-readRDS('myeloid_Fig3E.RDS') #cell1, metastatic myeloid cells; cell2, non-metastatic myeloid cells
data<-Embeddings(object = obj, reduction = "harmony")[,1:50]
cell1=colnames(obj)[obj$tissue.level.8 == 'Metastasis']
query=data[cell1,];
cell2=colnames(obj)[obj$tissue.level.8 != 'Metastasis'];
data2=data[cell2,]
meta2<-data.frame(tissue.level=obj@meta.data[cell2,'tissue.level.8'])
qout <- queryKNN(data2, query, k=30);
index=qout$index;
res<-NULL;
for(i in 1:nrow(index)){
  aa=as.character(meta2[index[i,],])
  res<-rbind(res,aa)
}
knn=data.frame(cell=cell1,GAC_Primary=rowSums(res=='GAC_Primary'),GAC_PBMC=rowSums(res=='GAC_PBMC'))

#permutation
res.per<-NULL;
for(j in 1:1000){
  bb<-table(sample(meta2[,1],30));
  res.per<-rbind(res.per,bb)
}
sort(table(paste0(res.per[,2],'_',res.per[,1])))
#p < 0.05 GAC_Primary>23 or GAC_PBMC>16
#p < 0.01 GAC_Primary>25 or GAC_PBMC>18

label<-rep('other',dim(knn)[1]);
label[knn$GAC_Primary>25]<-'GAC_Primary_similar';
label[knn$GAC_PBMC>18]<-'GAC_PBMC_similar';
knn$label<-label;
knn$celltype<-obj@meta.data[cell1,'celltype']
celltype<-unique(knn$celltype);
dat<-NULL;
for(i in 1:length(celltype)){
  aa=knn;
  aa[aa$celltype != celltype[i],'celltype']<-'other';
  aa$celltype<-factor(aa$celltype,levels=c(celltype[i],'other'))
  aa[aa$label != 'GAC_Primary_similar','label']<-'other';
  bb=table(aa[,c('celltype','label')]);
  test.res1<-fisher.test(bb);
  
  aa=knn;
  aa[aa$celltype != celltype[i],'celltype']<-'other';
  aa$celltype<-factor(aa$celltype,levels=c(celltype[i],'other'))
  aa[aa$label != 'GAC_PBMC_similar','label']<-'other';
  bb=table(aa[,c('celltype','label')]);
  test.res2<-fisher.test(bb);
  
  cc<-knn[knn$celltype==celltype[i],]
  data=c(test.res1$p.value,test.res1$estimate,test.res2$p.value,test.res2$estimate,
         length(which(cc$label=='GAC_Primary_similar')),length(which(cc$label=='GAC_PBMC_similar')),length(which(cc$label=='other')));
  dat<-rbind(dat,data)
}
dat<-as.data.frame(dat)
dat$celltype<-celltype;
colnames(dat)<-c('GAC_Primary_similar.p','GAC_Primary_similar.odds.ratio','GAC_PBMC_similar.p','GAC_PBMC_similar.odds.ratio',
                 'GAC_Primary_similar','GAC_PBMC_similar','other','celltype')
dat1<-dat[dat$GAC_Primary_similar > dat$GAC_PBMC_similar ,];
dat1<-dat1[dat1$GAC_Primary_similar.p < 0.01,];
dat1$type<-'GAC_Primary_similar';
dat1$odds.ratio<- 0 - dat1$GAC_Primary_similar.odds.ratio;
dat1$p.value<-dat1$GAC_Primary_similar.p;
dat2<-dat[dat$GAC_Primary_similar < dat$GAC_PBMC_similar ,];
dat2<-dat2[dat2$GAC_PBMC_similar.p < 0.01,];
dat2$type<-'GAC_PBMC_similar';
dat2$odds.ratio<- dat2$GAC_PBMC_similar.odds.ratio;
dat2$p.value<-dat2$GAC_PBMC_similar.p;
dat3<-rbind(dat1,dat2)
dat3[dat3$odds.ratio> 100,'odds.ratio']<- 100;
dat3[dat3$odds.ratio< -100,'odds.ratio']<- -100;
dat3$p.value<- -log10(dat3$p.value)
dat3[dat3$p.value > 50,'p.value']<- 50;


ggdotchart(dat3, x = "celltype", y = "odds.ratio",
                  color = "p.value",                                # Color by groups
                  add = "segments", 
                  #palette = c('grey','darkblue') ,               # Add segments from y = 0 to dots
                  rotate = TRUE,                                # Rotate vertically
                  dot.size = 6,                                 # Large dot size
                  ggtheme = theme_pubr()                        # ggplot2 theme
)+gradient_color(c( "red"))

#### Fig 4A #######
# 4A 
load('data_Fig4A.RData'); #including Embeddings and meta
res<-NULL;
cluster<-unique(meta$cluster)
for(i in 1:length(cluster)){
  cells=rownames(meta)[meta$cluster==cluster[i]]
  data<-Embeddings[,cells]
  data<-apply(data,1,mean);
  res<-cbind(res,data)
}
colnames(res)<-cluster;

dd <- as.dist((1 - cor(res))/2);
#dd<-dist(t(res))
data.hclust <- hclust(dd);

library('ape')
treeB<-as.phylo(data.hclust)
meta.data<-unique(meta[,c('celltype.big.L1','cluster')])
groupInfo<-list();
celltype<-unique(meta.data$celltype.big.L1)
for(i in 1:length(celltype)){
  groupInfo[[i]]<-meta.data[meta.data$celltype==celltype[i],'cluster']
}
names(groupInfo)<-celltype;
treeB <- groupOTU(treeB, groupInfo)
color.cell.type<-c('#1F78B4','#B2DF8A','#FB9A99','#E31A1C','#FDBF6F','#FF7F00','#CAB2D6','#33A02C','#6A3D9A','grey')
names(color.cell.type)<-c('B cell','CD4T','Mast','Myeloid','NK','pDCs','Plasma','CD8T','Stromal','T cell')

p1<-ggtree( treeB, aes(color=group),
            ladderize =F,right=T) + 
  geom_tiplab(size=3,vjust=-0.3,hjust=1)+
  scale_color_manual(values=c(color.cell.type,'grey') )

#Ro/e
data<-table(meta[,c('tissue.level.8','cluster')]);
data<-t(data)
res<-chisq.test(data)
ratio<-res$observed/res$expected;
ratio[ratio>3]<- 3;
ratio<-ratio[c(
  'Proliferative plasma cell','Proliferative NK cell','Proliferative B cell','Proliferative CD8 T cell','Proliferative CD4 T cell',
  'Plasma_C4_IGHG1-4','Plasma_C1_IGHA1/2','Plasma_C0_IGHA1/2_BCMA',
  'pDCs','B cell_C2_CXCR4_LAPTM5','B cell_C3_CD83_CD69','Mast',
  'DNT','NK_C6_XCL1/2_KLRC1_FCER1G','NK_C7_CD39_KIR2DL4','NK_C3_XCL1/2_KLRC1_GZMK','NKT',
  'NK_C2_CD16_CD69_CCL3','NK_C4_CD16_KLRC2','NK_C1_CD16_KLRF1',
  'CD8_C5_TRAV1-2_SLC4A10_KLRB1',
  'gdT','CD8_C9_FCMR_TIM3_CD27_GZMK','CD8_C7_ISG15_IFI6_IFIT3/1',
  'CD8_C0_CMC1_GZMK_HLA-II','CD8_C2_HLA-II_GZMH_LGALS1','CD8_C4_NR4A1_FOS_JUN',
  'CD8_C8_TIM3_LAG3_TIGIT_PD-1','CD8_C1_NR4A2_SPRY1','CD8_C6_GNLY_KLRD1',
  'CD4_C3_Treg','CD4_C5_CCL20_KLRB1','CD4_C6_CXCL13_TIGIT_TOX2',
  'CD4_C4_CCL5_GZMA/K','CD8_C3_CCR7_SELL_LEF1_TCF7','CD4_C2_FOS_JUN',
  'CD4_C1_S100A4_LTB','CD4_C0_CCR7_SELL_TCF7',
  'DC_C13_LAMP3','DC2_C5_FCER1A_CD1C_CLEC10A','DC1_C14_CLEC9A','Proliferating DC_C9_AXL_CD1C',
  'TAM_C0_CD163_APOE_C1QA/B/C','TAM_C4_APOC1_C1QC_APOE','Neutrophil_C3_IL1B_CXCL8_CXCL2/3',
  'Classical Mono_C2_S100A12','Classical Mono_C1_S100A12','Non-classical Mono_C6_LILRA1','TAM_C7_CCL2_FN1_CD163',
  'Endothelial_C13_FABP4_RGS5','Endothelial_C4_CLDN5_ACKR1','Endothelial_C7_VWF_ENG_ESM1',
  'Endothelial_C0_FABP5_PLVAP','Mesothelial_C8_CALB2_UPK3B_SLPI','Mesothelial_C12_proliferating',
  'VSMC_C14_ACTG2_MYH11','VSMC_C3_RGS5_ACTA2','Fibroblast_C9_NRG1_AREG_CXCL14_APOE',
  'Fibroblast_C6_POSTN_CXCL14_VSTM2A','Fibroblast_C1_CFD_FBLN1_CCDC80_SFRPs','Fibroblast_C15_CXCLs_MMP1_G0S2',
  'Fibroblast_C2_COL1A1/2_MMP11_ASPN'),
  rev(c('Healthy_PBMC','GAC_PBMC','NAT','CAG','IM','GAC_Primary','Metastasis'))];
library(pheatmap);
pdf('Allcell_Roe.pdf')
pheatmap(ratio,
         color = c(colorRampPalette(c('#313695','white'))(29),
                   #colorRampPalette(c('#e0f3f8','#ffffbf'))(2),
                   colorRampPalette(c('white','#ff7f00'))(60)),
         cluster_rows = F,
         cluster_cols = F);
dev.off();

#### Fig 4B #####
load(file='meta.Rdata')
meta3<-meta
data=meta3[,c('sample','celltype.harmony.2')];
colnames(data)[2]<-'celltype.big'
data1=as.data.frame(table(data))
data2 = data1 %>% group_by(sample) %>% summarise(count=sum(Freq))
data3<-inner_join(data1,data2,by='sample')
data3$percentage<- data3$Freq*100/data3$count;
data3$tissue.level<-meta[match(data3$sample,meta$sample),'tissue.level.8'];
aa=table(meta3$celltype.harmony.2);
data3<-data3[data3$celltype.big %in% names(aa)[aa>200],]

celltype<-unique(data3$celltype.big);
res<-NULL;
for(i in 1:length(celltype)){
  dat1<-data3[data3$celltype.big==celltype[i],];
  dat2<-data3[data3$celltype.big=='Endothelial_C7_VWF_ENG_ESM1',];
  cor.spearman=round(cor(dat1$percentage,dat2$percentage,method='spearman'),4)
  cor.spearman.p=cor.test(dat1$percentage,dat2$percentage,method='spearman')$p.val
  cor.pearson=round(cor(dat1$percentage,dat2$percentage),4) 
  cor.pearson.p=cor.test(dat1$percentage,dat2$percentage)$p.val
  aa<-data.frame(celltype1=celltype[i],celltype2='Endothelial_C7_VWF_ENG_ESM1',cor.spearman,cor.spearman.p,cor.pearson,cor.pearson.p);
  res<-rbind(res,aa);
}
res$cor.spearman.p.adjust<-p.adjust(res$cor.spearman.p,'BH');

res<-res[res$cor.spearman.p.adjust < 0.05,];
celltype1=res$celltype1[order(res$cor.spearman)];
res$celltype1<-factor(res$celltype1,levels = celltype1);
cor.type<-NULL;
cor.type[res$cor.spearman > 0]<-'pos';
cor.type[res$cor.spearman < 0]<-'neg';
res$cor.type<-cor.type;

plotx<-ggplot(res, aes(x=celltype1,y=cor.spearman)) + 
  geom_bar(aes(fill = cor.type),stat="identity") + 
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle = 90));
ggsave('barplot_Fig4B.pdf',useDingbats=F,plotx,width = 12)


#### Fig 5A ######
library(corrplot)
library(reshape2);
library(pheatmap);
res2<-read.table('relative_abundance_Fig5A.txt',sep='\t',header=T)
colnames(res2)<-gsub('\\.','-',colnames(res2))
load(file='meta.Rdata')
color.tissue<-c('#C0AB34','#0D4C7A','#C25757','#72ADE0','#F8A729','#9A253E','#EDAEAE','#5AAA46');
names(color.tissue)<-c('CAG','GAC_PBMC','GAC_Primary','Healthy_PBMC','IM','Metastasis','NAT','NGT')

color.cell.type<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a')
names(color.cell.type)<-c('B cell','CD4T','Mast','Myeloid','NK','pDCs','Plasma','CD8T','Stromal','T cell')

annotation_colors<-list(tissue=color.tissue,cell.type=color.cell.type);
load(file='meta.Rdata')
annotation_col=unique(meta[,c('sample','tissue.level.8')]);
annotation_col2=data.frame(annotation_col[,2])
rownames(annotation_col2)<-annotation_col$sample;
colnames(annotation_col2)<-'tissue'

annotation_row=unique(meta[,c('celltype.harmony.2','celltype.big.L1')]);
annotation_row2=data.frame(annotation_row[,2])
rownames(annotation_row2)<-annotation_row$celltype.harmony.2;
colnames(annotation_row2)<-'cell.type'

pdf('heatmap_prop_sample.pdf',height = 10,width = 14)
pheatmap(res2,
         color = c(colorRampPalette(c('#313695','white'))(29),
                   #colorRampPalette(c('#e0f3f8','#ffffbf'))(2),
                   colorRampPalette(c('white','#ff7f00'))(60)),
         #scale ='column',
         cluster_rows = F,
         cluster_cols = T,
         annotation_row=annotation_row2,
         annotation_col = annotation_col2,
         annotation_colors = annotation_colors);
dev.off();

