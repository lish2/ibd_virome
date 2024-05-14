
#alpha===========================================================================================
map = read.csv('data/mapping.file/mapping.file',sep='\t',stringsAsFactors = F)
prof = read.csv('data/profile/votu.tpm.prof20240306',row.names = 1,check.names = F,sep='\t',stringsAsFactors = F)


#分类型##
e_vir = c('Anelloviridae','Adenoviridae','Circoviridae','Genomoviridae','Geminiviridae',
          'Smacoviridae','Parvoviridae','Retroviridae','Iridoviridae','Hepadnaviridae','Redondoviridae',
          'Bornaviridae','Polyomaviridae','Metaviridae','Mimiviridae','Virgaviridae')
vtax = read.csv('data/data_generation/04.tax/votu.tax_family',sep='\t',header = F,stringsAsFactors = F)
vtax = vtax[which(vtax$V4%in%e_vir),]
vtax$V1 = ckv$ID[match(vtax$V1,ckv$contig_id)]
ff = intersect(vtax$V1,row.names(prof))
vtax = vtax[vtax$V1%in%ff,]
prof = prof[vtax$V1,]


#prof = read.csv('data/profile/mag.tpm.prof',row.names = 1,check.names = F,sep='\t',stringsAsFactors = F)

###


obs = apply(prof,2,function(x){sum(x!=0)})
shan = vegan::diversity(t(prof),index='shannon')
Alp = data.frame(ID=colnames(prof),obs=obs,shan=shan,Type=map$Type[match(colnames(prof),map$ID)],
                 Group=map$Group[match(colnames(prof),map$ID)])
Alp$Group = factor(Alp$Group,c('CON','CD','UC'))

fig=c()
for(j in c('VLP','WMS')){
  for(i in c('obs','shan')){
    if(i == 'obs'){
      fig[[paste(i,j)]]=ggplot(Alp[Alp$Type==j,],aes_string(x='Group',y=i,fill='Group'))+
        geom_boxplot(outlier.shape = 1,outlier.size = 0.5,width=.6)+
        stat_compare_means(comparisons = list(c('CON','UC'),c('CON','CD'),c('CD','UC')),method = 'wilcox',size=3)+
        scale_y_log10()+
        scale_color_manual(values = c("#6464ff","#ffb03f","#ff6464"))+
        scale_fill_manual(values = c("#6464ff","#ffb03f","#ff6464"))+
        theme_classic()
    }else{
      fig[[paste(i,j)]]=ggplot(Alp[Alp$Type==j,],aes_string(x='Group',y=i,fill='Group'))+
        geom_boxplot(outlier.shape = 1,outlier.size = 0.5,width=.6)+
        stat_compare_means(comparisons = list(c('CON','UC'),c('CON','CD'),c('CD','UC')),method = 'wilcox',size=3)+
        scale_color_manual(values = c("#6464ff","#ffb03f","#ff6464"))+
        scale_fill_manual(values = c("#6464ff","#ffb03f","#ff6464"))+
        theme_classic()
    }
  }
}
ggarrange(plotlist = fig,align = 'hv')

i='obs'
j='VLP'
wilcox.test(Alp[Alp$Type==j & Alp$Group=='CON',i],
            Alp[Alp$Type==j & Alp$Group!='CON',i])

##原核病毒与细菌丰度的多样性相关性

#vs_Alp=Alp
#bac_Alp=Alp

dat = vs_Alp
dat$Sample=sub('-V','',dat$ID)
dat$bac_obs = bac_Alp$obs[match(dat$Sample,bac_Alp$ID)]
dat$bac_shan = bac_Alp$shan[match(dat$Sample,bac_Alp$ID)]

fig=c()
for(i in c('obs','shan')){
  for(j in c('VLP','WMS')){
    datf = dat[dat$Type==j,]
    aa = cor.test(datf[,i],datf[,paste('bac',i,sep='_')],method = 'spearman')
    aa$p.value
    aa$estimate
    fig[[paste(i,j)]]=
      ggplot(datf,aes_string(x=i,y=paste('bac',i,sep='_')))+
      geom_point(alpha=.5,size=.5)+
      geom_smooth(method = 'lm')+
      labs(title = paste(j,' ',i,'\nspearman: rho=',round(aa$estimate,4),'; pval=',signif(aa$p.value,4),sep=''))+
      theme_test()+
      theme(
        title = element_text(
          size = 7
        )
      )
  }
}
ggarrange(plotlist = fig,align = 'hv')

#HC vs IBD科水平差异================================================================================================

vtax = read.csv('data/mapping.file/votu.tax_family20220426',sep='\t',stringsAsFactors = F,header = F)
prof = read.csv('data/profile/votu.tpm.prof20240306',row.names = 1,check.names = F,sep='\t',stringsAsFactors = F)
ff = intersect(vtax$V1,row.names(prof))
vtax = vtax[vtax$V1%in%ff,]
pf= prof[vtax$V1,]
pf = aggregate(pf,list(tax=vtax$V4),sum)
aa=colSums(pf[,-1])
pf = pf[head(order(-rowMeans(pf[,-1])),n=20),]
summary(colSums(pf[,-1])/aa)
row.names(pf) = pf$tax; pf = pf[,-1]
pf['Others',] = 1 -colSums(pf)

ff2=c('VLP','WMS')
tab = c()
for(y in 1:2){
  #x=1;y=1
  mapf = map[map$Group2%in%c('CON','IBD'),]
  tt1 = wtest2(pf,mapf[mapf$Type==ff2[y],],'Group2')
  tt1$logP = -log10(tt1$pvalue)
  tt1$logP[tt1$CON>tt1$IBD] = -tt1$logP[tt1$CON>tt1$IBD]
  tt1$g = ff2[y]
  tt1$Group = paste(ff2[y],ff2[y],sep=' ')
  tab=rbind(tab,tt1)
}

write.table(tab,'aa',sep='\t',quote = F,row.names = F)


tabf = tab#[grep('UC',tab$Group),]
#tabf = tabf[tabf$ID%in%unique(tabf$ID[which(abs(tabf$logP)>(-log10(0.05)))]),]
tabf = tabf[tabf$ID!='Others',]
aa = tabf[tabf$g == 'VLP',]
aa$logP[is.na(aa$logP)] = 0
tabf$ID = factor(tabf$ID,aa$ID[order(aa$logP)])
tabf$Col = 'a'
tabf$Col[which(tabf$logP<0)] = 'b'
tabf$logP[which(tabf$logP < -4)] = -4
tabf$logP[which(tabf$logP > 4)] = 4

# aa = as.data.frame(sort(rowMeans(pf)))
# aa$ID = row.names(aa)
# colnames(aa)[1] = 'val'
# aa = aa[order(aa$val), ]
# tabf$ID = factor(tabf$ID, aa$ID)

ggplot(tabf)+
  geom_bar(aes(y=ID,x=logP,fill=Col),stat='identity',width = .6)+
  #geom_point(aes(y=ID,x=logP,shape=Group))+
  #geom_segment(aes(x=0,xend=logP,y=ID,yend=ID))+
  geom_line(aes(x=logP,y=ID,group=ID))+
  geom_vline(xintercept = c(log10(0.05),0,-log10(0.05)))+
  facet_grid(~Group+g,space = 'free',scales = 'free')+
  xlim(-4,4)+
  theme_test()+
  theme(
    panel.spacing.y = unit(0, "lines")
  )

aa$ID = factor(aa$ID,aa$ID)
ggplot(aa)+
  geom_bar(aes(y=ID,x=val),stat='identity',position = 'dodge',width = .6)+
  #geom_point(aes(y=ID,x=val))+
  #geom_segment(aes(x=0,xend=logP,y=ID,yend=ID))+
  #geom_line(aes(x=logP,y=ID,group=ID))+
  #geom_vline(xintercept = c(log10(0.05),0,-log10(0.05)))+
  #facet_grid(~Group+g,space = 'free',scales = 'free')+
  scale_x_sqrt(breaks=c(0.01,0.2,0.4,0.6))+
  theme_classic()+
  theme(
    panel.spacing.y = unit(0, "lines")
  )

#Anelloviridae、Inoviridae=====================================================================================================
p = read.csv('data/profile/votu.tpm.prof20240306',sep='\t',row.names = 1,check.names = F)
vtax = read.csv('data/mapping.file/votu.tax_family20220426',sep='\t',header = F)
map  = read.csv('data/mapping.file/mapping.file',sep='\t',stringsAsFactors = F)

pf = p[row.names(p)%in%vtax$V1[vtax$V4=='Inoviridae'],]
tt1 = wtest2(pf,map[map$Type=='VLP',],'Group2')
tt2 = wtest2(pf,map[map$Type=='WMS',],'Group2')

ttf = tt2
ttf = ttf[!is.na(ttf$pvalue),]
ttf$gg = 'ibd_0'
ttf$gg[ttf$pvalue<0.05 & ttf$FC<1] = 'ibd_1'
ttf$gg[ttf$pvalue<0.05 & ttf$FC>1] = 'hc_1'
ttf$gg[ttf$pvalue>=0.05 & ttf$FC>1] = 'hc_0'
count(ttf$gg)


#Micro病毒特点=====================================================================================

p = read.csv('data/profile/votu.tpm.prof20240306',sep='\t',row.names = 1,check.names = F)
vtax = read.csv('data/mapping.file/votu.tax_family20220426',sep='\t',header = F)
map  = read.csv('data/mapping.file/mapping.file',sep='\t',stringsAsFactors = F)
vtax = vtax[vtax$V1%in%row.names(p),] 

pf = p[vtax$V1[vtax$V4=='Microviridae'],map$ID[map$Type=='WMS']] 
dat = data.frame(WMS=rowMeans(pf))
pf = p[vtax$V1[vtax$V4=='Microviridae'],map$ID[map$Type=='VLP']] 
dat$VLP = rowMeans(pf)

#write.table(dat,'aa',sep='\t',quote = F)

#Caudovirales、Petitvirales病毒特点=====================================================================================

p = read.csv('data/profile/votu.tpm.prof20240306',sep='\t',row.names = 1,check.names = F)
vtax = read.csv('data/mapping.file/votu.tax_family20220426',sep='\t',header = F)
map  = read.csv('data/mapping.file/mapping.file',sep='\t',stringsAsFactors = F)

db = read.csv('/share/data2/guorc/Database/virus/20220331/db.tax',stringsAsFactors = F,sep = '\t')
db = unique(db$tax)
tab = data.frame(ID=unique(vtax$V4),lab=NA)
for(x in 1:nrow(tab)){
  tab$lab[x]=db[grep(tab$ID[x],db)[1]]
}
tab$ord = NA
tab$ord[grep('Caudovirales',tab$lab)]='Caudovirales'
tab$ord[tab$ID%in%c('Flandersviridae','Gratiaviridae','Quimbyviridae')] = 'Caudovirales'
#tab$ord[grep('Petitvirales',tab$lab)] = 'Petitvirales'
tab$ord[is.na(tab$ord)] = 'Others'

p=p[vtax$V1[vtax$V4%in%tab$ID[tab$ord=='Caudovirales']],]
#p=p[vtax$V1[vtax$V4%in%tab$ID[tab$ord=='Petitvirales']],]
p = as.data.frame(t(as.data.frame(colSums(p,na.rm=T))))
p['others',] = 1-colSums(p) 

mapf = map
tt1 = wtest2(p,mapf[mapf$Type=='VLP',],'Group2')
tt2 = wtest2(p,mapf[mapf$Type=='WMS',],'Group2')

dat = data.frame(val=t(p[,mapf$ID])[,1],Group=mapf$Group,Group2=mapf$Group2,type=mapf$Type)
dat$Group = factor(dat$Group,c('CON','CD','UC'))
f1=ggplot(dat,aes(x=Group2,y=val,fill=Group2))+
  geom_boxplot()+
  facet_grid(~type)+
  scale_y_sqrt()+
  #stat_compare_means(comparisons = list(c('CON','CD'),c('CD','UC'),c('CON','UC')),method = 'wilcox')
  stat_compare_means(comparisons = list(c('CON','IBD')),method = 'wilcox')

ggarrange(f1,f2,align = 'hv')
