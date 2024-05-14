source('/share/data2/guorc/script/Rmylib.R')
library(vegan)
library(ape)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(pROC)
library(parallel)


mycolor = c("#5893c8","#e27933","#a5a5a5","#f4ba1b","#416cb4","#6ea647","#245a8a","#984823","#636363","#947226","#254173","#426630","#79a8d3","#e89259",
            "#b8b7b7","#f5c63b","#6587c4","#88b966","#3476b5","#c95e17","#848485","#c49611","#325697","#588638","#9abdde","#ecad81","#c9c9c8","#f6d368",
            "#8ba5d3","#a4ca8b","#1e4c75","#803c1f","#525252","#7b5f26","#1e3661","#37552a","#89b2d9","#e9a06d","#c0c0bf","#f5cc50","#7896cb","#95c179",
            "#2a689f","#b05321","#737474","#ad831e","#2a4c85","#4d7634")


# 不同family的宿主分布===============================================================================================
p = read.csv('data/data_generation/05.host/votu.host.tax',sep='\t',header=F,check.names = F)
#mtax = read.csv('data/mapping.file/MAG.tax.mapping.20210526',sep='\t',stringsAsFactors = F)
vtax = read.csv('data//mapping.file/votu.tax_family20220426',header=F,sep='\t',stringsAsFactors = F)
p = p[,c('V1','V3','V6','V7')]
colnames(p)[2:4] = c('phylum','family','genus')
p$phylum = sub('_[A-Z]$','',p$phylum)
vtax = vtax[,4:5]
colnames(vtax) = c('Family','ID')
p = merge(vtax[,c('ID','Family')],p,by.x = 'ID',by.y='V1',all.x = T)
#p = merge(p,mtax[,c('Ref','phylum','family')],by.x = 'V2',by.y='Ref',all.x = T)
p[is.na(p)] = 'Unknown'

aa = plyr::count(vtax$Family)
aa = aa[head(order(-aa$freq),n=15),]
p$Family[!p$Family%in%aa$x] = 'Others'

# aggregate(vtax$ContigSize,list(Group=vtax$Family),mean)
# aggregate(vtax$ContigSize,list(Group=vtax$Family),function(x){sd(x)})
# dcast(Family~checkv_quality,data=count(merge(vtax,checkv[,c('contig_id','checkv_quality')],
#                                              by.x='ID',by.y='contig_id')[,c('Family','checkv_quality')]))
# count(unique(p[p$V2!='Unknown',c('ID','Family')])$Family)
# 
# 

lab = 'phylum'
dat = plyr::count(p[,c('ID','Family',lab)])
ff = plyr::count(dat$ID)
ff = ff$x[ff$freq>1]
a1 = unique(dat[dat$ID%in%ff,c('ID','Family')])
a1[,lab] = 'Multiple'
a2 = dat[!dat$ID%in%ff,c('ID','Family',lab)]
dat = rbind(a1,a2)
#write.table(dat,'aa',sep='\t',quote = F,row.names = F)
dat = plyr::count(dat[,2:3])

aa = dat
aa[aa[,lab]!='Unknown',lab] = 'Known'

Or = aggregate(dat$freq,list(dat$Family),sum)
Or = Or$Group.1[order(-Or$x)]
dat$Family = factor(dat$Family,Or)
Or = aggregate(dat$freq,list(dat[,lab]),sum)
Or = Or$Group.1[order(-Or$x)]
dat[,lab] = factor(dat[,lab],Or)
#write.table(dat,'phylum.host',sep='\t',quote = F,row.names = F)


aa=dcast(Family~phylum,data=dat)
aa[is.na(aa)]=0
#write.table(aa,'aa',sep='\t',quote = F,row.names = F)

ff = levels(dat$Family)
dat$Group = 'G1'
dat$Group[dat$Family%in%ff[6:length(ff)]] = 'G2'

f2 = ggplot(dat)+
  geom_bar(aes(x=Family,y=freq,fill=phylum),position = 'stack',
           stat = 'identity',width = 0.6,color=NA)+
  scale_y_continuous(limits =c(0,160))+
  scale_fill_manual(values = c(brewer.pal(12,'Paired'),brewer.pal(12,'Set3'),
                               brewer.pal(9,'Set1')))+
  facet_grid(~Group,scales = 'free',space = 'free')+
  theme_bw()+
  theme(
    axis.text.x = element_text(
      #  size=6,
      angle=90,
      hjust=1,
      vjust=0.5
    ),
    legend.key.size = unit(0.5, "cm"),
    axis.line = element_blank(),
    #axis.text.y = element_text(
    #  size=5,
    #  angle=45
    #)
  )
f2
ggarrange(f1,f2,nrow = 2,align = 'hv')

#门以下用上==========================

lab = 'family'
dat = plyr::count(p[,c('ID','Family',lab)])
ff = plyr::count(dat$ID)
ff = ff$x[ff$freq>1]
a1 = unique(dat[dat$ID%in%ff,c('ID','Family')])
a1[,lab] = 'Multiple'
a2 = dat[!dat$ID%in%ff,c('ID','Family',lab)]
dat = rbind(a1,a2)
#write.table(dat,'aa',sep='\t',quote = F,row.names = F)
dat = plyr::count(dat[,2:3])

dat$Family = as.character(dat$Family)
dat[,lab] = as.character(dat[,lab])
aa=split(dat,dat$Family)
tab = c()
for(x in aa){
  p1=x
  if(nrow(p1)>7){
    p1=p1[order(-p1$freq),]
    pf=data.frame(Family=p1$Family[1],group='Others',freq=sum(p1$freq[11:nrow(p1)]))
    colnames(pf)[2]=lab 
    p1=rbind(p1[1:7,],pf)
  }
  tab = rbind(tab,p1)
}
tab[which(!tab$freq>1 & tab[,lab]!='Unknown'),lab] = 'Others'
tab=aggregate(tab$freq,tab[,c('Family',lab)],sum)
unique(tab$family)
colnames(tab)[3] = 'freq'
dat = tab

aa=dcast(Family~family,data=dat)
aa[is.na(aa)]=0
#write.table(aa,'bb',sep='\t',quote = F,row.names = F)


#Col=data.frame(Group=unique(dat[,lab]),Col=c(brewer.pal(12,'Paired'),brewer.pal(12,'Set3'),
#                                             brewer.pal(9,'Set1'),brewer.pal(8,'Accent'))[1:length(unique(dat[,lab]))])
Col=data.frame(Group=unique(dat[,lab]),Col=colorRampPalette(brewer.pal(8,'Spectral'))(length(unique(dat[,lab]))))
dat = merge(dat,Col,by.x = lab,by.y='Group')

#write.table(dat,'family.host',sep='\t',quote = F,row.names = F)

Or = aggregate(dat$freq,list(dat$Family),sum)
Or = Or$Group.1[order(-Or$x)]
f1 =c()
#f2 =c()
#f3 =c()
for(x in Or){
  pF=dat[dat$Family==x,]
  pF=pF[order(-pF$freq),]
  pF[,lab]=factor(pF[,lab],pF[,lab])
  f1[[x]]= ggplot(pF)+
    geom_bar(aes(x=Family,y=freq,fill=family),position = 'fill',
             stat = 'identity',width = 0.6,color=NA)+
    #scale_y_continuous(limits =c(0,2500) )+
    scale_fill_manual(values = pF$Col)+
    guides(fill=F,color=F)+
    theme_bw()+
    theme(
      axis.text.x = element_text(
        size=5,
        #angle=45,
        #hjust=1
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank()
      #axis.text.y = element_text(
      #  size=5,
      #  angle=45
      #)
    )
}
ggarrange(plotlist = f1,align = 'hv')

aa=unique(dat[,c(lab,'Col')])
aa$phylum = p$phylum[match(aa$family,p$family)]
aa = aa[order(aa$phylum),]
aa[,lab] = factor(aa[,lab],aa[,lab])
ggplot(aa)+
  geom_tile(aes(x=1,y=family,fill=family))+
  scale_fill_manual(values = aa$Col)

#map rate====================================================================================================
fpath = data.frame(ID=sub('.rc','',list.files('data/data_generation/01.prof/','*.rc')),stringsAsFactors = F)
for(x in 1:nrow(fpath)){
  aa = readLines(paste('data/data_generation/01.prof/',fpath$ID[x],'.log',sep=''))
  bb = aa[grep(' reads; of these:',aa)]
  if(length(bb)!=0){
    fpath[x,'allReads']=as.numeric(sub(' .*','',bb))
    fpath[x,'SP']='single'
    if(length(aa[grepl('were paired; of these:',aa)])==1){
      fpath[x,'allReads']=as.numeric(sub(' .*','',bb))*2
      fpath[x,'SP']='paired'
    }
  }
}

prof = read.csv('data/data_generation/ibd.raw.prof',sep='\t',check.names = F, row.names = 1,stringsAsFactors = F)
fpath$vs = colSums(prof[,fpath$ID])

#区分wms/vlp
clstr = read.csv('data/data_generation/00.Data/all.virus.clstr',sep='\t',stringsAsFactors = F,header = F)

clstr_map = read.csv('/share/data1/lish/fungus/2021_ibd_mycobiome/06.virome_profile/catalog2/tot.fna.bt.cvg.uniq',
                     sep='\t',stringsAsFactors = F,header = F)
clstr_map_len = read.csv('/share/data1/lish/fungus/2021_ibd_mycobiome/06.virome_profile/catalog2/tot.fna.len',
                         sep='\t',stringsAsFactors = F,header = F)
clstr_map$V3 = clstr_map$V3/clstr_map_len$V2[match(clstr_map$V1,clstr_map_len$V1)]
clstr_map = clstr_map[clstr_map$V3>=0.8 & clstr_map$V4>=95,]
clstr_map$group = 'WMS'
clstr_map$group[grepl('-V',clstr_map$V1)] = 'VLP'

clstr = tidyr::separate_rows(clstr,'V3', sep = ",")
clstr[,c('WMS','VLP')]= NA
for(x in 1:nrow(clstr)){
  ff = clstr_map$group[clstr_map$V2%in%clstr$V3[x]]
  if(grepl('-V',clstr$V3[x]) & !grepl('[|]',clstr$V3[x])){ff = c(ff,'VLP')}
  if(!grepl('-V',clstr$V3[x]) & !grepl('[|]',clstr$V3[x])){ff = c(ff,'WMS')}
  clstr[x,unique(ff)] = 1
}
clstr[is.na(clstr)] = 0

fpath$WMS = colSums(prof[unique(clstr$V1[clstr$WMS==1]),fpath$ID])
fpath$VLP = colSums(prof[unique(clstr$V1[clstr$VLP==1]),fpath$ID])
fpath[,4:6] = sweep(fpath[,4:6],1,fpath$allReads,'/')

dat = data.frame(ID=fpath$ID,val=fpath$vs,stringsAsFactors = F)
dat$Group = 'all_WMS'
dat$Group[grepl('-V',dat$ID)] = 'all_VLP'
ff = fpath[!grepl('-V',fpath$ID),]
dat = rbind(dat,data.frame(ID=ff$ID,val=ff$WMS,Group='WMS',stringsAsFactors = F))
ff = fpath[grepl('-V',fpath$ID),]
dat = rbind(dat,data.frame(ID=ff$ID,val=ff$VLP,Group='VLP',stringsAsFactors = F))

dat$Group =  factor(dat$Group,c('all_WMS','WMS','all_VLP','VLP'))
aggregate(dat$val,list(dat$Group),mean)

ggplot(dat)+
  geom_boxplot(aes(x=Group,y=val))+
  scale_y_sqrt()

#ncbi
fpath = data.frame(ID=sub('.rc','',list.files('data/data_generation/14.NCBI_map/','*.rc')),stringsAsFactors = F)
for(x in 1:nrow(fpath)){
  aa = readLines(paste('data/data_generation/14.NCBI_map/',fpath$ID[x],'.log',sep=''))
  bb = aa[grep(' reads; of these:',aa)]
  if(length(bb)!=0){
    fpath[x,'allReads']=as.numeric(sub(' .*','',bb))
    fpath[x,'SP']='single'
    if(length(aa[grepl('were paired; of these:',aa)])==1){
      fpath[x,'allReads']=as.numeric(sub(' .*','',bb))*2
      fpath[x,'SP']='paired'
    }
  }
  aa = read.csv(paste('data/data_generation/14.NCBI_map/',fpath$ID[x],'.rc',sep=''),sep='\t',header = F)
  fpath[x,'map']=sum(aa$V2)
}

datf = dat[grep('all',dat$Group),]
datf = data.frame(ID=datf$ID,val=datf$val,group='this',g='VLP')
datf $g[!grepl('-V',datf$ID)] = 'bulk'
dd = data.frame(ID=fpath$ID,val=fpath$map/fpath$allReads,group='RefSeq',g='VLP')
dd$g[!grepl('-V',dd$ID)] = 'bulk'
dd = rbind(dd,datf)
dd$g =factor(dd$g,c('VLP','bulk'))

ggplot(dd)+
  geom_boxplot(aes(x=g,y=val,fill=group),outlier.size = .5)+
  scale_y_sqrt(breaks=c(0.1,0.4,0.7,1))+
  theme_bw()



#病毒丰度表计算========================================================================================================================
map = read.csv('data/mapping.file/mapping.file',sep='\t',stringsAsFactors = F)
ckv = read.csv('data/data_generation/00.Data/votu.ckv',sep='\t',check.names = F,stringsAsFactors = F)
ckv$ID = sprintf("v%04x",1:nrow(ckv))

func = function(x){
  p = read.csv(paste('data/data_generation/19.review1/01.prof/',x,'.depth',sep=''),header = F,sep='\t')
  
  #使用kraken评估结果，至少得有1条的read能够特异地比对到指定的votu
  kra = read.csv(paste('data/data_generation/19.review1/02.kprof/',x,'.report',sep=''),header = F,sep='\t')
  kra = gsub(' ','',kra$V8[kra$V6=='S' & kra$V3!=0])
  p = p[p$V1%in%ckv$contig_id[ckv$ID%in%kra],]
  
  #保留每个votu的测序深度在10%-90%之间的位点，并计算这些位点的平均测序深度
  p = split(p,p$V1)
  dat = c()
  for(y in names(p)){
    pf = p[[y]]
    length_data = ckv$contig_length[ckv$contig_id==y]
    lis = sort(c(pf$V3,rep(0,length_data-nrow(pf))))
    dat = rbind(dat,data.frame(Sample=x,votu=ckv$ID[ckv$contig_id==y],rb=mean(lis[ceiling(length_data * 0.1):floor(length_data * 0.9)])))
  }
  
  #归一化作为每个votu的相对丰度
  dat = dat[dat$rb!=0,]
  dat$rb = dat$rb/sum(dat$rb)
  return(dat)
}

library(parallel)
cl <- makeCluster(50) # 初始化线程数
clusterExport(cl, varlist = list("ckv"),envir=environment()) #设定全局变量, 多个直接加进list，如list("p","d","a")
results <- parLapply(cl,map$ID,func) # lapply的并行版本
pred <- do.call('rbind',results) # 整合结果
stopCluster(cl) # 关闭集群


dat = dcast(votu~Sample,data=pred)
dat[is.na(dat)] = 0
colnames(dat)[1] = 'ID'

#write.table(dat,'data/profile/votu.tpm.prof20240306',sep='\t',quote=F,row.names = F)


#pcoa=======================================================================================
map = read.csv('data/mapping.file/mapping.file',sep='\t',stringsAsFactors = F)
prof = read.csv('data/profile/votu.tpm.prof20240306',row.names = 1,check.names = F,sep='\t',stringsAsFactors = F)

#VLP vs WMS### 
mapf = data.frame(SampleID=map[,1],Group=map$Type,stringsAsFactors = F)
#CON vs IBD###
mapf = data.frame(SampleID=map[,1],Group=map$Group,stringsAsFactors = F)
mapf = mapf[map$Type=='WMS',]
####

#科水平====
vtax = read.csv('data/mapping.file/votu.tax_family20220426',sep='\t',header = F,stringsAsFactors = F)
vtax = vtax[vtax$V1%in%row.names(prof),]
gg = vtax$V4[match(row.names(prof),vtax$V1)]
gg[is.na(gg)] = 'Unclassified'
p = aggregate(prof,list(id=gg),sum)
row.names(p) = p$id
p = p[,-1]
prof = p
#===

prof = prof[,mapf$SampleID]
prof = sqrt(t(prof[rowSums(prof)!=0,]))
data.distan = vegdist(prof)

mm = mapf[mapf$Group!='CON',]
data.distanf = as.dist(as.matrix(data.distan)[mm$SampleID,mm$SampleID])
ado = adonis2(data.distanf~Group,data=mm)
ado
RsquareAdj(ado$R2[1],n=nrow(mapf),m=ado$Df[1])


obs.pcoa=pcoa(data.distan,correction="none", rn=NULL)
dv <- as.data.frame(obs.pcoa$vectors)
dp <- round(100*obs.pcoa$values$Eigenvalues/sum(obs.pcoa$values$Eigenvalues),0)
dv <- dv[mapf[,1],]
dv[,'Group']<-as.character(mapf[,2])

Down = dv[,c('Axis.1','Group')]
Left = dv[,c('Axis.2','Group')]

dv$Group = factor(dv$Group,c('CON','CD','UC'))
#dv$Group = factor(map$Stage[match(row.names(dv),map$ID)],c('CON','Remission','Moderately active','Severely active'))

ggplot(data = dv,
       aes(x=Axis.1,
           y=Axis.2)) +
  geom_hline(yintercept = 0,linetype=2,color='grey')+
  geom_vline(xintercept = 0,linetype=2,color='grey')+
  geom_point(aes(color=Group),size=0.5)+
  stat_ellipse(aes(color=Group,fill=Group),size=0.5, geom="polygon", level=0.8, alpha=0.3) +
  labs(title="", x=paste("PC1 ",dp[1],'%',sep=''), y=paste("PC2 ",dp[2],'%',sep=''))+
  theme_classic() +
  scale_color_manual(values = c("#ff6464","#6464ff","#ffb03f"))+
  scale_fill_manual(values = c("#ff6464","#6464ff","#ffb03f"))+
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

#科水平比较======================================================================

vtax = read.csv('data/data_generation/04.tax/votu.tax_family',sep='\t',stringsAsFactors = F,header = F)
vtax$V1 = ckv$ID[match(vtax$V1,ckv$contig_id)]
prof = read.csv('data/profile/votu.tpm.prof20240306',row.names = 1,check.names = F,sep='\t',stringsAsFactors = F)
map = read.csv('data/mapping.file/mapping.file',sep='\t',stringsAsFactors = F)

#VLP vs WMS### 
aa = plyr::count(map$Sample)
mf = map[map$Sample%in%aa$x[aa$freq>1],]
mf = mf[order(mf$ID),]

ff = intersect(vtax$V1,row.names(prof))
vtax = vtax[vtax$V1%in%ff,]
pf= prof[vtax$V1,mf$ID]
pf = aggregate(pf,list(tax=vtax$V4),sum)
pf = pf[head(order(-rowMeans(pf[,-1])),n=15),]
row.names(pf) = pf$tax; pf = pf[,-1]
pf['Others',] = 1 -colSums(pf)

tt=wtest2(pf,mf,'Type',T,'Sample')


Or1 = colnames(pf)[order(-pf[1,],-pf[2,])]
Or2 = row.names(pf)
pf = melt(as.matrix(pf))
pf$Group = map$Type[match(pf$Var2,map$ID)]
pf$Var2 =factor(pf$Var2,Or1)
pf$Var1 = factor(pf$Var1,Or2)

ggplot(pf)+
  geom_bar(aes(x=Var2,y=value,fill=Var1),position = 'stack',stat='identity')+
  facet_grid(~Group,scales = 'free',space = 'free')+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c(brewer.pal(7,'Set2'),brewer.pal(8,'Set3'),'#e5e5e5'))+
  theme_test()+
  theme(
    axis.text.x =element_blank(),
    strip.text.x = element_text(
      size = 10
    ),
    #strip.background = element_rect(
    #  fill= NA,
    #  color = NA
    #),
    axis.ticks.x = element_blank(),
    panel.spacing.x = unit(0.3, "lines")
  )

ggplot(pf)+
  geom_boxplot(aes(x=Var1,y=value,fill=Group),outlier.shape = 1,outlier.size = 0.5,width=0.6)+
  scale_y_continuous(trans='sqrt')+
  scale_fill_manual(values = c("#ffb03f","#6464ff"))+
  theme_classic()+
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

aa = aggregate(pf$value,pf[,c('Group','Var1')],mean)
aa$Var2=factor(aa$Var2,Or2)
ggplot(aa)+
  geom_bar(aes(x=1,y=x,fill=Var1),stat = 'identity',position = 'fill')+
  coord_polar(theta = 'y', start = 0, direction = 1)+
  facet_grid(Group~.)+
  theme_test()




#db可靠性比较=======================================================================================================
gg=c('votu','GVD','GPD','MGV','refseq')
myc = t(combn(gg,2))
myc = split(myc,1:nrow(myc))

dat = c()
for(x in gg){
  if(x == 'votu'){
    p = read.csv(paste('data/data_generation/00.Data/votu.ckv.d/quality_summary.tsv',sep=''),sep='\t',stringsAsFactors = F)
  }
  p = read.csv(paste('/share/data2/guorc/Project/oral_virome/2021.oral_virome/2021.Jul20.oral_virome/data/data_generation/02.quality/',x,'.ckv/quality_summary.tsv',sep=''),sep='\t',stringsAsFactors = F)
  p$rate = p$viral_genes/p$gene_count
  dat = rbind(dat,data.frame(ID=p$contig_id,value = p$rate,Group=x))
}


dat$Group = factor(dat$Group,gg) 
f1=ggplot(dat,aes(x=Group,y=value))+
  #stat_compare_means(comparisons = myc,method = 'wilcox')+
  geom_boxplot(aes(fill=Group),width=0.6,outlier.shape = 21,outlier.fill = NA,outlier.size = 0.1)+
  #geom_density(aes(x=value,group=Group,color=Group))+
  scale_fill_manual(values=brewer.pal(8,'Accent'))+
  guides(fill=F)+
  theme_light()+
  labs(x='',y='The rate of viral genes (Checkv)')+
  theme(
    axis.text.y = element_text(
      size=8,
      color = 'black'
    ),
    axis.text.x = element_text(
      size=8,
      color = 'black',
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )



dat = c()
for(x in gg){
  if(x == 'votu'){
    p = read.csv(paste('data/data_generation/00.Data/votu.vif/votu.fa_gt0bp_dvfpred.txt',sep=''),sep='\t',stringsAsFactors = F)
  }
  p = read.csv(paste('/share/data2/guorc/Project/oral_virome/2021.oral_virome/2021.Jul20.oral_virome/data/data_generation/02.quality/',x,'.vif/',x,'.fa_gt0bp_dvfpred.txt',sep=''),sep='\t',stringsAsFactors = F)
  p$Group = x
  dat = rbind(dat,p)
}

dat$Group = factor(dat$Group,gg) 
f2=ggplot(dat,aes(x=Group,y=score))+
  #stat_compare_means(comparisons = myc,method = 'wilcox')+
  geom_boxplot(aes(fill=Group),width=0.6,outlier.shape = 21,outlier.fill = NA,outlier.size = 0.1)+
  scale_fill_manual(values=brewer.pal(8,'Accent'))+
  guides(fill=F)+
  theme_light()+
  labs(x='',y='Score (DeepVirFinder)')+
  theme(
    axis.text.y = element_text(
      size=8,
      color = 'black'
    ),
    axis.text.x = element_text(
      size=8,
      color = 'black',
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

#geom_density(aes(x=value,group=Group,color=Group))


dat = c()
for(x in gg){
  if(x=='votu'){
    p = as.numeric(sub(' .*','',readLines('data/data_generation/00.Data/votu.vib/VIBRANT_votu/VIBRANT_log_run_votu.log')[10:11]))
  }else{
    p = as.numeric(sub(' .*','',readLines(paste("/share/data2/guorc/Project/oral_virome/2021.oral_virome/2021.Jul20.oral_virome/data/data_generation/02.quality/",x,".vib/VIBRANT_",x,"/VIBRANT_log_run_",x,".log",sep=''))[12:13]))
  }
  p = data.frame(Group=x,all=p[1],phage=p[2])
  dat = rbind(dat,p)
}
dat$value=dat$phage/dat$all
dat$LAB = paste(round(dat$value*100,2),'%',sep='')

dat$Group = factor(dat$Group,gg) 
f3=ggplot(dat,aes(x=Group,y=value))+
  #stat_compare_means(comparisons = myc,method = 'wilcox')+
  geom_bar(aes(fill=Group),width=0.6,stat='identity')+
  geom_text(aes(label=LAB))+
  scale_fill_manual(values=brewer.pal(8,'Accent'))+
  guides(fill=F)+
  theme_light()+
  labs(x='',y='The ratio of phages identified by VIBRANT')+
  theme(
    axis.text.y = element_text(
      size=8,
      color = 'black'
    ),
    axis.text.x = element_text(
      size=8,
      color = 'black',
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )



dat = c()
for(x in gg){
  if(x =='votu'){
    p = read.csv(paste("data/data_generation/00.Data/votu.visorter/final-viral-score.tsv",sep=''),sep='\t',stringsAsFactors = F)
    p2 = read.csv('data/data_generation/00.Data/votu.fa.len',header = F,sep='\t',stringsAsFactors = F)
  }
  p = read.csv(paste("/share/data2/guorc/Project/oral_virome/2021.oral_virome/2021.Jul20.oral_virome/data/data_generation/02.quality/",x,".vis/final-viral-score.tsv",sep=''),sep='\t',stringsAsFactors = F)
  p2 = read.csv(paste("/share/data2/guorc/Project/oral_virome/2021.oral_virome/2021.Jul20.oral_virome/data/data_generation/00.data/",x,".fa.sort.len",sep=''),header = F,sep='\t',stringsAsFactors = F)
  p = data.frame(Group=x,value=nrow(p)/nrow(p2))
  dat = rbind(dat,p)
}

dat$LAB = paste(round(dat$value*100,2),'%',sep='')

dat$Group = factor(dat$Group,gg) 
f4=ggplot(dat,aes(x=Group,y=value))+
  #stat_compare_means(comparisons = myc,method = 'wilcox')+
  geom_bar(aes(fill=Group),width=0.6,stat='identity')+
  geom_text(aes(label=LAB))+
  scale_fill_manual(values=brewer.pal(8,'Accent'))+
  guides(fill=F)+
  theme_light()+
  labs(x='',y='The ratio of phages identified by VirSorter')+
  theme(
    axis.text.y = element_text(
      size=8,
      color = 'black'
    ),
    axis.text.x = element_text(
      size=8,
      color = 'black',
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )


ggarrange(f1,f2,f3,f4,f1,f2,f3,f4,ncol = 4)