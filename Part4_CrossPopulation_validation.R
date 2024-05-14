
#公共数据的病毒丰度表计算========================================================================================================================
flist = list.files('data/data_generation/19.review1/03.pub/','*.depth')
map = data.frame(ID=sub('.depth','',flist))
ckv = read.csv('data/data_generation/00.Data/votu.ckv',sep='\t',check.names = F,stringsAsFactors = F)
ckv$ID = sprintf("v%04x",1:nrow(ckv))

func = function(x){
  p = read.csv(paste('data/data_generation/19.review1/03.pub/',x,'.depth',sep=''),header = F,sep='\t')
  
  #使用kraken评估结果，至少得有1条的read能够特异地比对到指定的votu
  kra = read.csv(paste('data/data_generation/19.review1/03.pub/',x,'.report',sep=''),header = F,sep='\t')
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

#write.table(dat,'data/profile/pub.votu.tpm.prof20240311',sep='\t',quote=F,row.names = F)



#所有数据集的多样性================================================================================================

map_pub = read.csv('data/mapping.file/pub.mapping.file20240311',sep='\t',check.names = F, stringsAsFactors = F)
map = read.csv('data/mapping.file/mapping.file',sep='\t',check.names = F, stringsAsFactors = F)
map$Group[map$Group=='CON'] = 'con'
map$Group[map$Group%in%c('UC','CD')] = 'case'
map = rbind(data.frame(SampleID=map$ID,Group=map$Group,Project=map$Type),map_pub[,1:3])


prof = read.csv('data/profile/votu.tpm.prof20240306',sep='\t',check.names = F, row.names = 1,stringsAsFactors = F)
pub = read.csv('data/profile/pub.votu.tpm.prof20240311',sep='\t',check.names = F, row.names = 1,stringsAsFactors = F)
pub = pub[row.names(prof),]
row.names(pub) = row.names(prof)
pub[is.na(pub)] = 0
Prof = cbind(prof,pub)
Prof = Prof[,map$SampleID]


#prof = Prof[,mapf$SampleID]
prof = Prof
prof = sqrt(t(prof[rowSums(prof)!=0,]))
Data.distan = vegdist(prof)


mapf = data.frame(SampleID=map[,1],Group=map$Group,Proj=map$Project,stringsAsFactors = F)
gg = c('VLP','Stockdale_ISR')
mapf = mapf#[mapf$Proj%in%gg,]

data.distan = as.dist(as.matrix(Data.distan)[mapf$SampleID,mapf$SampleID])

ado = adonis2(data.distan~Proj,data=mapf)
ado
RsquareAdj(ado$R2[1],n=nrow(mapf),m=ado$Df[1])


obs.pcoa=pcoa(data.distan,correction="none", rn=NULL)
dv <- as.data.frame(obs.pcoa$vectors)
dp <- round(100*obs.pcoa$values$Eigenvalues/sum(obs.pcoa$values$Eigenvalues),0)
dv <- dv[mapf[,1],]
dv[,'Group']<-as.character(mapf[,2])
dv[,'Proj']<-as.character(mapf$Proj)

Down = dv[,c('Axis.1','Group','Proj')]
Left = dv[,c('Axis.2','Group','Proj')]

#dv$Group = factor(dv$Group,c('CON','CD','UC'))
#dv$Group = factor(map$Stage[match(row.names(dv),map$ID)],c('CON','Remission','Moderately active','Severely active'))

p = ggplot(data = dv,
           aes(x=Axis.1,
               y=Axis.2)) +
  geom_hline(yintercept = 0,linetype=2,color='grey')+
  geom_vline(xintercept = 0,linetype=2,color='grey')+
  geom_point(aes(color=Proj,shape=Group),size=1)+
  #stat_ellipse(aes(color=Group,fill=Group),size=0.5, geom="polygon", level=0.8, alpha=0.3) +
  #labs(title="", x=paste("PC1 ",dp[1],'%',sep=''), y=paste("PC2 ",dp[2],'%',sep=''))+
  xlab("") +
  ylab("") +
  theme_classic() +
  #scale_color_manual(values = c("#ff6464","#6464ff"))+
  #scale_color_manual(values = mycolor)+
  scale_color_manual(values=brewer.pal(12,"Paired")[c(1:10,12)])+
  scale_shape_manual(values = c(1,16))+
  guides(fill=FALSE,color=F,shape=F)+
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
p<- ggplotGrob(p)

#right_data_box$Group = factor(right_data_box$Group, levels=f)
d <- ggplot(Left)+geom_boxplot(aes(x = Proj,y = Left[,1],fill = Proj),width = 0.5,outlier.size = .5)+
  theme_bw()+theme(panel.grid =element_blank())+
  #scale_fill_wsj("colors6", "Group")
  scale_fill_manual(values=brewer.pal(12,"Paired")[c(1:10,12)])+
  #scale_fill_manual(values = mycolor)+
  guides(fill=FALSE)+theme(axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(),
        text=element_text(
          family="Helvetica",
          face="bold",
          colour="black",
          size=12
        ))+
  #ylim(min(right_data_box[,2], 0)*yy, max(right_data_box[,2])*yy)+
  # ylim(min(right_data_box[,2], 0)*3, max(right_data_box[,2])*3)+
  xlab("")+ylab(paste("PC2 ",dp[2],'%',sep=''))
#lift_data_box$Group = factor(lift_data_box$Group, levels=f)
b<- ggplot(Down)+geom_boxplot(aes(x = Proj,y = Down[,1],fill = Proj),width = 0.5,outlier.size = .5)+
  theme_bw()+theme(panel.grid =element_blank())+coord_flip()+
  #scale_fill_manual(values = mycolor)+
  scale_fill_manual(values=brewer.pal(12,"Paired")[c(1:10,12)])+
  guides(fill=FALSE)+theme(axis.text.y = element_blank())+
  theme(axis.ticks = element_blank(),
        text=element_text(
          family="Helvetica",
          face="bold",
          colour="black",
          size=12
        ))+
  #scale_fill_manual(values=brewer.pal(9,"Set1"))+
  #ylim(min(lift_data_box[,2], 0)*xx, max(lift_data_box[,2])*xx)+
  xlab("")+ylab(paste("PC1 ",dp[1],'%',sep=''))
a<-ggplot()+theme_bw()+theme(panel.border = element_blank(),panel.grid =element_blank(),
                             axis.text = element_blank(),axis.title = element_blank(),
                             axis.ticks = element_blank())+
  annotate("text", x=1, y=40, label=paste("P.value =",round(ado[[1]][6][[1]][1],4),'\n',
                                          "R2      =",round(ado[[1]][5][[1]][1],4)), size=3.5)
a <- ggplotGrob(a)
d <- ggplotGrob(d)
b <- ggplotGrob(b)

gridExtra::grid.arrange(d,p,a,b,ncol=2,widths=c(1,3),heights = c(3,1))


#marker在不同cohort的分布==========================================================================================
Prof = read.csv('data/profile/pub.votu.tpm.prof20240311',sep='\t',stringsAsFactors = F,check.names = F)
map = read.csv('data/mapping.file/pub.mapping.file20240311',sep='\t',stringsAsFactors = F)
p = read.csv('data/profile/votu.tpm.prof20240306',sep='\t',stringsAsFactors = F,check.names = F)
m = read.csv('data/mapping.file/mapping.file',sep='\t',stringsAsFactors = F)

Prof = merge(Prof,p,all=T)
Prof[is.na(Prof)] = 0  
Prof = Prof[Prof$ID%in%Marker$ID,]
row.names(Prof) = Prof$ID
Prof = Prof[,-1]

map = data.frame(ID=c(map$SampleID,m$ID),Group=c(map$Group,m$Group2),Project=c(map$Project,m$Type))
map$Group[map$Group=='CON'] = 'con'
map$Group[map$Group=='IBD'] = 'case'

gg = c("VLP","Stockdale_ISR","WMS","Federic_FRA","Federic_GER","Federic_ISR","Federic_USA","Franzosa_NED",
       "Franzosa_USA","He_CHN","Weng_CHN")

dat = c()
for(x in gg){
  mf = map[map$Project==x,]
  tt = wtest2(Prof,mf,'Group')
  tt$Project = x
  dat = rbind(dat,tt)
  
}

dat$FC = log2(dat$case/dat$con)
dat$FC[is.na(dat$FC)] = 0
dat$FC[dat$FC > 10] = 10
dat$FC[dat$FC < -10] = -10
# ff = which(dat$con > dat$case)
# dat$FC = -dat$con/dat$case
# ff = which(dat$con < dat$case)
# dat$FC = dat$case/dat$con

Or = dcast(dat[,c('ID','Project','FC')],ID~Project)
row.names(Or) = Or$ID ; Or = Or[,-1]
Or=matrix_or(Or)

dat$Group = Marker$Fill[match(dat$ID,Marker$ID)]
dat$ID = factor(dat$ID,Or$r_order)
dat$Project = factor(dat$Project,rev(gg))

#dat$ID = factor(dat$ID,aa$ID)

ggplot(dat)+
  geom_tile(aes(y=Project,x=ID,fill=FC))+
  #geom_text(aes(x=test,y=train,label=auc))+
  #scale_fill_gradientn(colors = brewer.pal(11,'Spectral'),limit=c(-10,10))+
  scale_fill_gradient2(low='#6464FF',mid = 'white',high='#FF6464',midpoint = 0,limit=c(-10,10))+
  facet_grid(~Group,scales = 'free',space = 'free',switch = 'y')+
  theme_test()+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    # axis.text.x = element_text(
    #   angle = 45,
    #   hjust = 1,
    #   vjust = 1
    # )
  )

dat$host = marker_host$genus[match(dat$ID,marker_host$ID)]
aa = unique(dat[,c('ID','Group','host')])
aa = aa[order(aa$host),]
#aa$ID = factor(aa$ID,as.character(aa$ID))
ggplot(aa)+
  geom_tile(aes(x=ID,y=1,fill=host))+
  #geom_text(aes(x=test,y=train,label=auc))+
  #scale_fill_gradientn(colors = brewer.pal(11,'Spectral'),limit=c(-10,10))+
  scale_fill_manual(values = mycolor)+
  #scale_fill_gradient2(low='#6464FF',mid = 'white',high='#FF6464',midpoint = 0,limit=c(-10,10))+
  facet_grid(~Group,scales = 'free',space = 'free',switch = 'y')+
  theme_test()+
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )


dat$ly = ly$V4[match(dat$ID,ly$ID)]
aa = unique(dat[,c('ID','Group','ly')])

dat$tax = vtax$V4[match(dat$ID,vtax$V1)]
aa = unique(dat[,c('ID','Group','tax')])

ggplot(aa)+
  geom_tile(aes(x=ID,y=1,fill=ly))+
  #geom_text(aes(x=test,y=train,label=auc))+
  #scale_fill_gradientn(colors = brewer.pal(11,'Spectral'),limit=c(-10,10))+
  scale_fill_manual(values = mycolor)+
  #scale_fill_gradient2(low='#6464FF',mid = 'white',high='#FF6464',midpoint = 0,limit=c(-10,10))+
  facet_grid(~Group,scales = 'free',space = 'free',switch = 'y')+
  theme_test()+
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

#write.table(dat,'xx',sep='\t',row.names = F,quote = F)


#
d = dat
d$g[d$con>d$case] = 'CON'
d$g[d$con<d$case] = 'IBD'
d$gf = (d$Group==d$g) 
d = d[which(d$gf==T),]
d = plyr::count(d$Project)
d$rate = d$freq/nrow(Marker)
d = d[!d$x%in%c('VLP','WMS'),]
sd(d$rate)

#cross cohort预测================================================================================================

map_pub = read.csv('data/mapping.file/pub.mapping.file20240311',sep='\t',check.names = F, stringsAsFactors = F)
map = read.csv('data/mapping.file/mapping.file',sep='\t',check.names = F, stringsAsFactors = F)
map$Group[map$Group=='CON'] = 'con'
map$Group[map$Group%in%c('UC','CD')] = 'case'
map = rbind(data.frame(SampleID=map$ID,Group=map$Group,Project=map$Type),map_pub[,1:3])

pub = read.csv('data/profile/pub.votu.tpm.prof20240311',sep='\t',check.names = F, row.names = 1,stringsAsFactors = F)
pub = pub[Marker$ID,]
pub[is.na(pub)] = 0 
row.names(pub) = Marker$ID
prof = read.csv('data/profile/votu.tpm.prof20240306',sep='\t',check.names = F, row.names = 1,stringsAsFactors = F)
prof = prof[Marker$ID,]

Prof = cbind(prof,pub)


prof = as.data.frame(t(Prof[,map$SampleID]))
rowSums(prof)
prof$Group = as.factor(map$Group)

source('/share/data2/guorc/script/Rmylib/RandomF.R')
#gg = sort(unique(map$Project))
gg = c("VLP","Stockdale_ISR","WMS","Federic_FRA","Federic_GER","Federic_ISR",
       "Federic_USA","Franzosa_NED","Franzosa_USA","He_CHN","Weng_CHN")


#项目内预测===
impo=c()
Self=c()
for(x in gg){
  mf = map[map$Project==x,]
  pf = as.data.frame(t(prof[mf$SampleID,colnames(prof)!='Group']))
  mf = Sampling_kFold_dsc(mf,'Group',5)
  Ra = kFold_random(pf,mf,Threads = 25)
  #tab = rbind(tab,data.frame(train=x,test=x,auc=round(mean(kFold_assess(Ra$pred)),2)))
  Self = rbind(Self,data.frame(train=x,test=x,Ra$pred))
  aa = aggregate(Ra$Impo$MeanDecreaseAccuracy,list(Ra$Impo$ID),mean)
  aa = data.frame(train=x,test=x,ID=aa$Group.1,imp=aa$x)
  aa = aa[order(-aa$imp),]
  aa$or = 1:nrow(aa)
  impo = rbind(impo,aa)
  print(x)
}


lab=c()
for(x in 1:length(gg)){
  aa = Self[Self$train==gg[x],]
  bb = kFold_assess(aa)
  aa = aa[aa$seed==which.min(abs(bb-mean(bb))),]
  if(x==1){
    plot.roc(aa$true,aa$pred,col=brewer.pal(12,'Paired')[x]) 
  }else{
    lines.roc(aa$true,aa$pred,col=brewer.pal(12,'Paired')[x])
  }
  lab = c(lab,paste(gg[x],as.character(round(mean(bb),1)),sep=':'))
}
legend(x=.5,y=.8, legend=lab, col=brewer.pal(12,'Paired')[1:length(gg)],lwd=2,pt.cex = 1,cex = .6)


#重要marker挑选===
marker_impo = data.frame(Group.1=unique(impo$ID[impo$or%in%1:5]))
marker_impo$tax = Marker$tax[match(marker_impo$Group.1,Marker$ID)]
marker_impo$host = marker_host$genus[match(marker_impo$Group.1,marker_host$ID)]
marker_impo$enrich = Marker$Fill[match(marker_impo$Group.1,Marker$ID)]


d = impo[impo$ID%in%marker_impo$Group.1,]
d$enrich = Marker$Fill[match(d$ID,Marker$ID)]
d$val = d$or
d$val[d$val>60] = 60

dh = read.csv('data/data_generation/05.host/votu.host.tax',sep='\t',stringsAsFactors = F,header = F)
dh$V1 = ckv$ID[match(dh$V1,ckv$contig_id)]
dh = dh[dh$V1%in%d$ID,]
dh$sp = paste(dh$V8,dh$V9,sep=' ')
dh = aggregate(dh$sp,list(dh$V1),function(x){paste(sort(x),collapse = ', ')})
d$host = dh$x[match(d$ID,dh$Group.1)]
d$tax = Marker$tax[match(d$ID,Marker$ID)]
d$ID[!is.na(d$tax)] = paste(d$tax,d$ID,sep=' ')[!is.na(d$tax)]
d$lab=d$or
d$lab[d$lab>60]=NA
d$train=factor(d$train,gg)

or = aggregate(d$or,list(d$ID),sum)
or = or$Group.1[order(or$x)]
d$ID =  factor(d$ID,or)
#d1 = d[d$enrich=='CON',]
#d2 = d[d$enrich=='IBD',]
d$val = 60-d$val
d$val[d$enrich=='CON'] = -d$val[d$enrich=='CON']


f1=ggplot(d)+
  geom_tile(aes(x=train,y=ID,fill=val))+
  geom_text(aes(x=train,y=ID,label=lab),size=3)+
  scale_fill_gradient2(low='#6464FF',mid = 'white',high = '#FF6464',limit=c(-60,60))+
  #scale_fill_gradient(low='#FF6464',high = 'white',limit=c(0,60))+
  facet_grid(enrich~.,scales = 'free',space = 'free')+
  theme_test()+
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

f2=ggplot(unique(d[,c('ID','host','enrich')]))+
  #geom_tile(aes(x=train,y=ID,fill=val))+
  geom_text(aes(x=1,y=ID,label=host),size=3)+
  #scale_fill_gradient(low='#6464FF',high = 'white',limit=c(0,60))+
  #scale_fill_gradient(low='#FF6464',high = 'white',limit=c(0,60))+
  facet_grid(enrich~.,scales = 'free',space = 'free')+
  theme_classic()+
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

ggarrange(f1,f2,align = 'hv')



#高贡献marker项目间预测===

tab = c()
for(x in gg){
  Train = prof[map$SampleID[map$Project==x],c(marker_impo$Group.1,'Group')]
  for(y in gg){
    if(x!=y){
      Test = prof[map$SampleID[map$Project==y],]
      Pred = Random(Test,Train)
      Pred$pred$seed = 1
      tab = rbind(tab,data.frame(train=x,test=y,auc=round(mean(kFold_assess(Pred$pred)),2)))
    }
  }
}

#高贡献marker项目内预测===

Self_f=c()
for(x in gg){
  mf = map[map$Project==x,]
  pf = as.data.frame(t(prof[mf$SampleID,marker_impo$Group.1]))
  mf = Sampling_kFold_dsc(mf,'Group',5)
  Ra = kFold_random(pf,mf,Threads = 25)
  #tab = rbind(tab,data.frame(train=x,test=x,auc=round(mean(kFold_assess(Ra$pred)),2)))
  Self_f = rbind(Self_f,data.frame(train=x,test=x,Ra$pred))
}


lab=c()
for(x in 1:length(gg)){
  aa = Self_f[Self_f$train==gg[x],]
  bb = kFold_assess(aa)
  aa = aa[aa$seed==which.min(abs(bb-mean(bb))),]
  if(x==1){
    plot.roc(aa$true,aa$pred,col=brewer.pal(12,'Paired')[x]) 
  }else{
    lines.roc(aa$true,aa$pred,col=brewer.pal(12,'Paired')[x])
  }
  lab = c(lab,paste(gg[x],as.character(round(mean(bb),1)),sep=':'))
}
legend(x=.5,y=.8, legend=lab, col=brewer.pal(12,'Paired')[1:length(gg)],lwd=2,pt.cex = 1,cex = .6)


for(x in gg){
  aa = Self_f[Self_f$train==x,3:7]
  tab = rbind(tab,data.frame(train=x,test=x,auc=round(mean(kFold_assess(aa)),2)))
}


#LOCO预测===

tab2 = c()
for(x in gg){
  Train = prof[map$SampleID[map$Project!=x],c(marker_impo$Group.1,'Group')]
  Test = prof[map$SampleID[map$Project==x],c(marker_impo$Group.1,'Group')]
  Pred = Random(Test,Train)
  Pred$pred$seed = 1
  ass = kFold_assess(Pred$pred)
  tab2 = rbind(tab2,data.frame(train='others',test=x,auc=ass))
}
tab2$test = factor(tab2$test,gg)
tab2$lab = round(tab2$auc,1)


#绘图

tab$fill = tab$auc
tab$fill[tab$fill<50]=50
tab$lab = round(tab$auc,1)

tab$train = factor(tab$train,rev(gg))
tab$test = factor(tab$test,gg)

Right = aggregate(tab$auc[tab$train!=tab$test],list(tab$train[tab$train!=tab$test]),function(x){round(mean(x),1)})
Down = aggregate(tab$auc[tab$train!=tab$test],list(tab$test[tab$train!=tab$test]),function(x){round(mean(x),1)})


f1=ggplot(tab)+
  geom_tile(aes(x=test,y=train,fill=fill))+
  geom_text(aes(x=test,y=train,label=lab))+
  scale_fill_gradient2(low='black',mid = 'red',high='yellow',midpoint = 70,limit=c(50,100))+
  theme_test()
# theme(
#   axis.text.x = element_text(
#     angle = 45,
#     hjust = 1,
#     vjust = 1
#   )
# )

f2=ggplot(Right)+
  geom_tile(aes(x=1,y=Group.1,fill=x))+
  geom_text(aes(x=1,y=Group.1,label=x))+
  scale_fill_gradient2(low='black',mid = 'red',high='yellow',midpoint = 70,limit=c(50,100))+
  theme_test()

f3=ggplot(Down)+
  geom_tile(aes(x=Group.1,y=1,fill=x))+
  geom_text(aes(x=Group.1,y=1,label=x))+
  scale_fill_gradient2(low='black',mid = 'red',high='yellow',midpoint = 70,limit=c(50,100))+
  theme_test()

aa = round(mean(c(tab$auc[tab$train==tab$test])),1)
f4=ggplot()+
  geom_tile(aes(x=1,y=1,fill=aa))+
  geom_text(aes(x=1,y=1,label=aa))+
  scale_fill_gradient2(low='black',mid = 'red',high='yellow',midpoint = 70,limit=c(50,100))+
  theme_test()


f5=ggplot(tab2)+
  geom_tile(aes(x=test,y=train,fill=auc))+
  geom_text(aes(x=test,y=train,label=lab))+
  scale_fill_gradient2(low='black',mid = 'red',high='yellow',midpoint = 70,limit=c(50,100))+
  theme_test()

bb = round(mean(tab2$auc),1)
f6=ggplot()+
  geom_tile(aes(x=1,y=1,fill=bb))+
  geom_text(aes(x=1,y=1,label=bb))+
  scale_fill_gradient2(low='black',mid = 'red',high='yellow',midpoint = 70,limit=c(50,100))+
  theme_test()

ggarrange(f1,f2,f3,f3,f4,f4,f5,f6,f4,align = 'hv',widths = c(6,2,2),heights = c(6,1,1))

