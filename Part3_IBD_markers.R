#宿主host的注释============================================================================================
ckv = read.csv('data/data_generation/00.Data/votu.ckv',sep='\t',check.names = F,stringsAsFactors = F)
ckv$ID = sprintf("v%04x",1:nrow(ckv))
Host = read.csv('data/data_generation/05.host/votu.host.tax',sep='\t',stringsAsFactors = F,header = F)
vtax = read.csv('data/mapping.file/votu.tax_family20220426',sep='\t',header = F,stringsAsFactors = F)
Host$V1 = ckv$ID[match(Host$V1,ckv$contig_id)]
Host$V8[!grepl('_unknown',Host$V8)] = paste(Host$V8,Host$V9,sep='_')[!grepl('_unknown',Host$V8)]

host = unique(Host[,c(1,8)])

Host=multip_vir_host(Host)

host$V8 = sub('_[A-Z]_','_',host$V8)
host$V8 = sub('_[A-Z]$','',host$V8)
host$V8 = sub('_sp\\d+$','_sp.',host$V8)
host = unique(host)
aa = aggregate(host$V8,list(host$V1),function(x){paste(x,collapse = ',')})
#write.table(aa,'aa',sep = '\t',quote = F,row.names = F)


#共有marker================================================================================================
vtax = read.csv('data/mapping.file/votu.tax_family20220426',sep='\t',header = F,stringsAsFactors = F)
prof = read.csv('data/profile/votu.tpm.prof20240306',row.names = 1,check.names = F,sep='\t',stringsAsFactors = F)
map = read.csv('data/mapping.file/mapping.file',sep='\t',stringsAsFactors = F)
prof = prof[rowMeans(prof)>1e-4,]

ff=c('CON','IBD')
mapf = map
tt1 = wtest2(prof,mapf[mapf$Type=='VLP',],'Group2')
tt2 = wtest2(prof,mapf[mapf$Type=='WMS',],'Group2')

tt = cbind(tt1,tt2)
tt$tax = vtax$V4[match(tt$ID,vtax$V1)]


#write.table(tt1,'../../tt1',sep = '\t',quote = F,row.names = F)
#write.table(tt2,'../../tt2',sep = '\t',quote = F,row.names = F)

tt1$logq = -log10(tt1$q)
tt1$logq[tt1[,ff[1]]>tt1[,ff[2]]] = -tt1$logq[tt1[,ff[1]]>tt1[,ff[2]]] 
tt2$logq = -log10(tt2$q)
tt2$logq[tt2[,ff[1]]>tt2[,ff[2]]] = -tt2$logq[tt2[,ff[1]]>tt2[,ff[2]]]

tab = as.data.frame(rowMeans(prof))
tab$vlp = tt1$logq
tab$wms = tt2$logq
colnames(tab)[1] = 'RelAbun' 
tab[is.na(tab)] = 0

cutoff = -log10(0.05)
tab$Fill = 'PASS'
tab$Fill[tab$vlp< -cutoff & tab$wms < -cutoff] = ff[1]
tab$Fill[tab$vlp> cutoff & tab$wms > cutoff] = ff[2]
tab$tax = vtax$V4[match(row.names(tab),vtax$V1)]
tab$Iden = vtax$V3[match(row.names(tab),vtax$V1)]

tabf = tab#[which(tab$RelAbun>0 & tab$tax=='Anelloviridae'),]
#tabf = tab[which(tab$RelAbun>1e-4),]

fisher.test(
  matrix(c(
    sum(tabf$vlp < 0 & tabf$wms < 0),
    sum(tabf$vlp < 0 & tabf$wms > 0),
    sum(tabf$vlp > 0 & tabf$wms < 0),
    sum(tabf$vlp > 0 & tabf$wms > 0)
  ),nrow=2)
)

fisher.test(
  matrix(c(
    sum(tabf$vlp < log10(0.05) & tabf$wms < log10(0.05)),
    sum(tabf$vlp < log10(0.05) & tabf$wms > -log10(0.05)),
    sum(tabf$vlp > -log10(0.05) & tabf$wms < log10(0.05)),
    sum(tabf$vlp > -log10(0.05) & tabf$wms > -log10(0.05))
  ),nrow=2)
)


sum((tabf$vlp< log10(0.05) & tabf$wms < log10(0.05)) | (tabf$vlp> -log10(0.05) & tabf$wms> -log10(0.05)))
sum((tabf$vlp<log10(0.05) & tabf$wms > -log10(0.05)) | (tabf$vlp> -log10(0.05) & tabf$wms<log10(0.05)))

Marker = tabf[tabf$Fill!='PASS',]
Marker$ID = row.names(Marker)
Marker$oldID = ckv$contig_id[match(Marker$ID,ckv$ID)]
#write.table(Marker,'aa',row.names=F,sep = '\t',quote = F)


ggplot(tabf)+
  geom_point(aes(x=wms,y=vlp,
                 #size=RelAbun,
                 color=Fill),size=.8)+
  #stat_smooth(aes(x=wms,y=vlp),method = 'glm')+
  geom_vline(xintercept = c(-cutoff,cutoff),linetype=3)+
  geom_hline(yintercept = c(-cutoff,cutoff),linetype=3)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  scale_color_manual(values = c('steelblue','red','grey'))+
  scale_y_continuous(breaks = seq(-9,9,3))+
  scale_x_continuous(breaks = seq(-9,9,3))+
  labs(x='-log10(q value) [wms]',y='-log10(q value) [vlp]')+
  theme_test()

cor.test(tabf$vlp,tabf$wms)

#write.table(tab,'../../tab',sep = '\t',quote = F)



#marker的类型分布
aa=plyr::count(tabf[,c('Fill','tax')])
aa=aa[aa$Fill!='PASS',]
aa[is.na(aa)] = 'Unclassified'

ggplot(aa)+
  geom_bar(aes(y=tax,x=freq,fill=Fill),position = 'dodge',stat='identity',width = 0.7)+
  theme_classic()+
  scale_x_continuous(expand = c(0,0))


#marker的宿主
marker_host = tabf[tabf$Fill!='PASS',]
#marker_host$phylum = 'unclassified'
marker_host$ID= row.names(marker_host)
marker_host$phylum = Host$V3[match(marker_host$ID,Host$V1)]
marker_host$genus = Host$V7[match(marker_host$ID,Host$V1)]
marker_host[is.na(marker_host)] = 'unclassified'

p = marker_host[,c('Fill','tax','genus')]
colnames(p) = c('a','b','c')
aa = plyr::count(p)
aa = aa[order(-aa$freq),]
ff = split(aa,aa[,c('a','b')])
Top=100
n=c()
for(x in ff){
  if(nrow(x)>Top){
    n=c(n,x$c[1:Top])
  }else{
    n=c(n,x$c)
  }
}
p$c[!p$c%in%n] = 'Others'

dat = plyr::count(p[,c('a','b','c')])
ff = aggregate(dat$freq,list(dat$b),sum)
ff = ff[order(-ff$x),]
dat$b = factor(dat$b,ff$Group.1)
Col = data.frame(id=unique(dat$c),
                 #Col=mycolor[1:length(unique(dat$c))])#,
                 Col=c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2'),brewer.pal(10,'Set3'))[1:length(unique(dat$c))]
)

dat$d = marker_host$phylum[match(dat$c,marker_host$genus)]
dat$d[dat$c%in%c('unclassified','multiple','Others')] = 'ZZZ'
dat = dat[order(dat$d,dat$freq),]
dat$Col=Col$Col[match(dat$c,Col$id)]
dat$Fill = paste(dat$a,dat$c,dat$b,sep=':')
dat$Fill = factor(dat$Fill,dat$Fill)
ggplot(dat)+
  geom_bar(aes(x=b,y=freq,fill=Fill),position = 'stack',stat='identity',width = .6)+
  scale_fill_manual(values = dat$Col)+
  facet_grid(~a,scales = 'free',space = 'free')+
  guides(fill=F)+
  theme_light()+
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust =1
    )
  )

Col$id =factor(Col$id,Col$id)
ggplot(Col)+
  geom_bar(aes(x=1,y=id,fill=id),stat='identity')+
  scale_fill_manual(values = Col$Col)


#marker的溶源性
p = read.csv('data/data_generation/11.marker_vib/marker.vib/VIBRANT_marker/VIBRANT_results_marker/VIBRANT_genome_quality_marker.tsv',
             sep='\t',stringsAsFactors = F)
p$scaffold = sub('_fragment.*| .*','',p$scaffold)
p = unique(p[,1:2])
p$Group = Marker$Fill[match(p$scaffold,Marker$oldID)]


#marker菌的功能比较========================================================================================================
kegg = read.csv('data/data_generation/06.kegg/votu_kegg.m8',sep='\t',header = F)
kegg_map = read.csv('/share/data2/guorc/Database/KEGG/20201208/all_ko.mapping',sep='\t',stringsAsFactors = F)
kegg = kegg[,1:2]
pfam = read.csv('data/data_generation/08.pfam/vOTU.pfam',sep='',header = F)
pfam = pfam[!grepl('predicted_active_site',pfam$V1),]
pfam = pfam[!pfam$V1%in%kegg$V1,]
kegg = rbind(kegg,data.frame(V1=pfam$V1,V2=paste(pfam$V6,pfam$V7,sep=' '),stringsAsFactors = F))

kegg$V2 = sub('.*\\|','',kegg$V2)
kegg$contigid = sub('_\\d+$','',kegg$V1)
kegg$ID = ckv$ID[match(kegg$contigid,ckv$contig_id)]
kegg$tax = vtax$V4[match(kegg$ID,vtax$V1)]

keggf = kegg[kegg$ID%in%Marker$ID,]
keggf$Enrich = Marker$Fill[match(keggf$ID,Marker$ID)]
keggf = unique(keggf[,c('ID','Enrich','V2')])
aa = plyr::count(keggf$ID)
bb = ckv[ckv$ID%in%Marker$ID,]
bb$Enrich = Marker$Fill[match(bb$ID,Marker$ID)]
bb$ko_num = aa$freq[match(bb$ID,aa$x)]
bb$rate = bb$ko_num/bb$gene_count
bb$rate[is.na(bb$rate)] = 0

summary(bb$rate[bb$Enrich=='IBD'])
summary(bb$rate[bb$Enrich=='CON'])

#keggf[is.na(keggf)] = 'unclassified'
keggf = plyr::count(keggf[,c('Enrich','V2')])
aa = plyr::count(Marker[,c('Fill')])
keggf$rate = keggf$freq/(aa$freq[match(keggf$Enrich,aa$x)])
keggf=dcast(V2~Enrich,data=keggf,value.var = 'rate')
keggf[is.na(keggf)] = 0

keggf = keggf[keggf$CON>0.1 | keggf$IBD>0.1,]
#keggf = keggf[keggf$CON*92>3 | keggf$IBD*165 > 3,]


keggf$pval=NA
keggf$odd=NA
for(x in 1:nrow(keggf)){
  tt = fisher.test(matrix(c(keggf$CON[x]*aa$freq[aa$x=='CON'],(1-keggf$CON[x])*aa$freq[aa$x=='CON'],
                            keggf$IBD[x]*aa$freq[aa$x=='IBD'],(1-keggf$IBD[x])*aa$freq[aa$x=='IBD']),nrow = 2))
  keggf$pval[x]= tt$p.value
  keggf$odd[x]= tt$estimate
}
keggf$anno = kegg_map$Anno[match(keggf$V2,kegg_map$KO)]
keggf$anno = sub(' \\[.*','',keggf$anno)
keggf$anno = paste(keggf$V2,keggf$anno,sep=': ')
keggf$anno = sub(': NA$','',keggf$anno)
keggf$qval = p.adjust(keggf$pval,method = 'BH')
#write.table(keggf,'aa',quote = F,sep='\t',row.names = F)

keggf=keggf[keggf$qval<0.05,]
keggf$enrich='CASE'
keggf$enrich[keggf$odd>1]='CON'
Or = keggf$anno[order(keggf$enrich,-rowSums(keggf[,c('IBD','CON')]))]
Or = data.frame(ID=Or,Group=rep(c('a','b'),100)[1:length(Or)],val=max(dat$value))
Or$ID = factor(Or$ID,Or$ID) 

dat = melt(keggf[,c('anno','CON','IBD')],id.vars = 'anno')
dat$anno = factor(dat$anno,Or$ID)
ggplot(dat)+
  geom_bar(data=Or,aes(x=ID,y=val,fill=Group),stat='identity',width=1,alpha=.05)+
  geom_point(aes(x=anno,y=value,fill=variable,color=variable))+
  scale_fill_manual(values = c('white','#4292C9','#6464ff','#ff6464'))+
  scale_color_manual(values = c('#6464ff','#ff6464'))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

#K00986
aa = kegg
aa = aa[aa$V2=='K00986' & !is.na(aa$enrich),]
aa$sp = Host$V3[match(aa$ID,Host$V1)]

bb = ckv[ckv$ID%in%aa$ID,]
bb$enrich=aa$enrich[match(bb$ID,aa$ID)]
bb$sp = Host$V3[match(bb$ID,Host$V1)]

#pcoa

keggf = kegg[kegg$ID%in%Marker$ID,]
keggf = plyr::count(keggf[,c('ID','V2')])
keggf = dcast(keggf,V2~ID)
row.names(keggf) = keggf$V2
keggf = keggf[,-1]
keggf[is.na(keggf)] = 0

keggff = keggf[apply(keggf,1,function(x){sum(x!=0)})>1,]
keggff['Others',] = colSums(keggf[!row.names(keggf)%in%row.names(keggff),])
#keggff = keggf


mapf = data.frame(ID=Marker$ID,Group=Marker$Fill,stringsAsFactors = F)
mapf = mapf[mapf$ID%in%colnames(keggff),]
data.distan = vegdist(t(keggff[,mapf$ID]),method = 'jaccard')

ado = adonis2(data.distan~Group,data=mapf)
ado
RsquareAdj(ado$R2[1],n=nrow(mapf),m=ado$Df[1])


obs.pcoa=pcoa(data.distan,correction="none", rn=NULL)
dv <- as.data.frame(obs.pcoa$vectors)
dp <- round(100*obs.pcoa$values$Eigenvalues/sum(obs.pcoa$values$Eigenvalues),0)
dv <- dv[mapf[,1],]
dv[,'Group']<-as.character(mapf[,2])
dv$Group = factor(dv$Group,c('CON','IBD'))



library(ggrepel)
ggplot(data = dv,
       aes(x=Axis.1,
           y=Axis.2)) +
  geom_hline(yintercept = 0,linetype=2,color='grey')+
  geom_vline(xintercept = 0,linetype=2,color='grey')+
  geom_point(aes(color=Group),size=0.5)+
  stat_ellipse(aes(color=Group,fill=Group),size=0.5, geom="polygon", level=0.8, alpha=0.3) +
  #geom_point(data=dfen,aes(x=A1,y=A2),color='grey',size=0.5)+
  #geom_text_repel(data=dfen,aes(x=A1,y=A2,label=LAB),position=position_stack(0.9),color='black',size=2.5) +
  labs(title="", x=paste("PC1 ",dp[1],'%',sep=''), y=paste("PC2 ",dp[2],'%',sep=''))+
  theme_classic() +
  #scale_color_manual(values = c("#6464ff","#ffb03f","#ff6464"))+
  #scale_fill_manual(values = c("#6464ff","#ffb03f","#ff6464"))+
  #guides(fill=FALSE,color=F)+
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    #axis.title = element_blank(),
    panel.border = element_rect(
      fill=NA,
    )
  )

dat = as.matrix(data.distan)
dat = mtrx2cols(dat)
dat$g1 = Marker$Fill[match(dat$Var1,Marker$ID)]
dat$g2 = Marker$Fill[match(dat$Var2,Marker$ID)]
dat = dat[dat$g1==dat$g2,]
dat = dat[dat$value!=0,]

ggplot(dat,aes(x=g1,y=value))+
  geom_boxplot(outlier.size = .5)+
  stat_compare_means(method = 'wilcox')+
  theme_bw()
