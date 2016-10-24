#!/usr/bin/env Rscript

##Author: Samson Jacob
##Purpose : Read in the PSM matrix file and analyze samples; written in a verbose way to be very clear!

##requires FORMATED_PSM.txt to be in the same subfolder! ##the output of the python script
PSM_mat = read.table('./FORMATED_PSM.txt',sep='\t',header=TRUE,row.names = 1)

##Generate a matrix of proteins that are identified in each sample
PSM_mat$Missing = apply(PSM_mat[,2:4],1,function(x){sum(is.na(as.numeric(x)))})
COMPL = subset(PSM_mat,PSM_mat$Missing==0)  ##472 total

##remove the Missing column
COMPL$Missing=NULL

##Get the medians of the samples from the complete cases
PSM_4L_MED = median(as.numeric(COMPL$PSM_4L1)) ## 4L Median = 11
PSM_Skmel_MED = median(as.numeric(COMPL$PSM_Skmel2)) ## Mewo Median = 11
PSM_Mewo_MED = median(as.numeric(COMPL$PSM_Mewo1)) ## Skmel Median = 14


##graph the distributions of the complete cases
library(ggplot2)
library(reshape)

MeltedDF = melt(COMPL)
ggplot(MeltedDF,aes(x=value, fill=variable)) + geom_density(alpha=0.20) +xlim(-3,250)+ggtitle('Complete Cases')+scale_fill_brewer(palette="Set1") +xlab('PSM_Value')


##check if above Median
COMPL$AboveMed_4L=ifelse(COMPL$PSM_4L1>11,1,0)
COMPL$AboveMed_Skmel=ifelse(COMPL$PSM_Skmel2>11,1,0)
COMPL$AboveMed_Mewo=ifelse(COMPL$PSM_Mewo1>14,1,0)

##Sum the number
COMPL$SumofAboveMedians = apply(COMPL[,7:9],1,function(x){sum(x)})

##final df; contains above median values for all 3 samples
final_df = subset(COMPL,COMPL$SumofAboveMedians==3)


##total proteins per subset
print(nrow(subset(COMPL,COMPL$AboveMed_4L!=0))) #220
print(nrow(subset(COMPL,COMPL$AboveMed_Skmel!=0))) #275
print(nrow(subset(COMPL,COMPL$AboveMed_Mewo!=0))) #185

##printout
f = subset(COMPL,COMPL$AboveMed_4L!=0)
f = f[c(1,5,6)]
write.table(f,'AboveMed4L.txt',sep='\t')


sk = subset(COMPL,COMPL$AboveMed_Skmel!=0)
sk1 = sk[c(1,5,6)]

write.table(sk1,'AboveMed_SkMel.txt',sep='\t')

mew =subset(COMPL,COMPL$AboveMed_Mewo!=0)
mew1 = mew[c(1,5,6)]
write.table(mew1,'AboveMed_Mew.txt',sep='\t')




##Pairwise overlap
nrow(subset(COMPL,COMPL$AboveMed_4L==1 & COMPL$AboveMed_Skmel==1))  ##SkMel and 4L have 180 overlapping above median
nrow(subset(COMPL,COMPL$AboveMed_4L==1 & COMPL$AboveMed_Mewo==1))  ##MeWo and 4L have 153 overlapping
nrow(subset(COMPL,COMPL$AboveMed_Mewo==1 & COMPL$AboveMed_Skmel==1))##MeWo and SkMel have 174




##see the counts
barplot(table(COMPL$SumofAboveMedians))
p = ggplot(data=COMPL)+geom_density(aes(COMPL$PSM_4L1,y=..density..),color='black')


##write_out the file of the grand_overlap
write.table(final_df,'Final_Overlapping_Matrix.txt',sep='\t')


##Subset Graphs; 
library(grid)
library(gridBase)

## 4L1
f = table(COMPL$PSM_4L1)
names(f)=paste('PSM',names(f),sep='_') ## change the names of the table to include "_"
vps = baseViewports() ## list
pushViewport(vps$inner, vps$figure, vps$plot) ##access the subsets of vps
grid.text(names(f),x=unit(names(f),'native'),y=unit(-1,"lines"),just='right',rot=50) ##adjust barplot
midpts <- barplot(f, col='blue', names.arg="")
text(x=midpts, y=-2, names(f), cex=0.8, srt=45, xpd=TRUE)
grid.text(names(f),
          x = unit(midpts, "native"), y=unit(-1, "lines"),
          just="right", rot=50)
f2 = data.frame(f)
f2$val = COMPL$PSM_4L1
ggplot(f2,aes(Var1,Freq))+geom_bar(stat='identity',fill = "black", colour = "yellow")+
  theme(axis.text.y = element_text(size = 5),panel.background = element_rect(fill = "white"))+
  geom_rect(data=NULL,aes(xmin=11,xmax=72,ymin=-Inf,ymax=Inf),alpha=.01,fill="#ebb5e1")+
  geom_rect(data=NULL,aes(xmin = 29,  xmax = 72, ymin = -Inf, ymax = Inf),alpha=.01,fill = "light blue") +
  coord_flip()+xlab('PSM_4L')+geom_vline(xintercept = c(11,29),colour=c('red','blue'))

###Skmel
g = table(COMPL$PSM_Skmel2)
names(g)=paste('PSM',names(g),sep='_')
grid.text(names(g),x=unit(names(g),'native'),y=unit(-1,"lines"),just='right',rot=50)
midpts2 <- barplot(g, col='blue', names.arg="")
text(x=midpts2, y=-2, names(g), cex=0.8, srt=45, xpd=TRUE)
grid.text(names(g),
          x = unit(midpts2, "native"), y=unit(-1, "lines"),
          just="right", rot=50)
g2= data.frame(g)
ggplot(g2,aes(Var1,Freq))+geom_bar(stat='identity',fill = "black", colour = "yellow")+
  theme(axis.text.y = element_text(size = 5),panel.background = element_rect(fill = "white"))+
  geom_rect(data=NULL,aes(xmin=11,xmax=77,ymin=-Inf,ymax=Inf),alpha=.01,fill="#ebb5e1")+
  geom_rect(data=NULL,aes(xmin = 29,  xmax = 77, ymin = -Inf, ymax = Inf),alpha=.01,fill = "light blue") +
  coord_flip()+xlab('PSM_SkMel147')+geom_vline(xintercept = c(11,29),colour=c('red','blue'))

###Mewo
v = table(COMPL$PSM_Mewo1)
names(v)= paste('PSM',names(v),sep='_')
grid.text(names(v),x=unit(names(v),'native'),y=unit(-1,'lines'),just='right',rot=50)
midpts3 <- barplot(v, col='blue', names.arg="")
text(x=midpts3, y=-2, names(v), cex=0.8, srt=45, xpd=TRUE)
grid.text(names(v),
          x = unit(midpts2, "native"), y=unit(-1, "lines"),
          just="right", rot=50)
v2= data.frame(v)
ggplot(v2,aes(Var1,Freq))+geom_bar(stat='identity',fill = "black", colour = "yellow")+
  theme(axis.text.y = element_text(size = 5),panel.background = element_rect(fill = "white"))+
  geom_rect(data=NULL,aes(xmin=11,xmax=72,ymin=-Inf,ymax=Inf),alpha=.01,fill="#ebb5e1")+
  geom_rect(data=NULL,aes(xmin = 29,  xmax = 72, ymin = -Inf, ymax = Inf),alpha=.01,fill = "light blue") +
  coord_flip()+xlab('PSM_MeWo')+geom_vline(xintercept = c(11,29),colour=c('red','blue'))


##Check the tau estimate for each sample; and apply 
library(robustbase)

robust_4L = round(PSM_4L_MED + 2 * scaleTau2(as.numeric(COMPL$PSM_4L1)),0) ## 29
robust_Skmel = round(PSM_Skmel_MED + 2 *scaleTau2(as.numeric(COMPL$PSM_Skmel2)),0) ## 36
robust_Mewo = round(PSM_Skmel_MED + 2 * scaleTau2(as.numeric(COMPL$PSM_Mewo1)),0) ## 37

##Tau Estimate
COMPL$Robust4L = ifelse(COMPL$PSM_4L1>robust_4L,1,0)
COMPL$RobustSKMEL= ifelse(COMPL$PSM_Skmel2>robust_Skmel,1,0)
COMPL$RobustMEWO= ifelse(COMPL$PSM_Mewo1>robust_Mewo,1,0)


##Count
nrow(subset(COMPL,COMPL$Robust4L!=0)) #83
nrow(subset(COMPL,COMPL$RobustSKMEL!=0)) #65
nrow(subset(COMPL,COMPL$RobustMEWO!=0)) #51

##check the sum
COMPL$SumTau = apply(COMPL[,11:13],1,function(x){sum(x)})
