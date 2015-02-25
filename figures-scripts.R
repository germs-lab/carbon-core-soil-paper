library(phyloseq)
library(plyr)
library(ggplot2)
library(reshape)

setwd("/Users/adina/Google Drive/C_Metabolism_Agg_Paper_2013/reproducibility/")
##################################################################################
#SUPP FIGURE 1
counts <- read.delim(sep=',',file="core-counts.csv",header=FALSE)
colnames(counts) <- c('filecount','uniquecazy','bpcount','round')
counts$filecount <- factor(counts$filecount, levels=c("1","2","3","4"))
counts$bocount <- factor(counts$bpcount, levels=c("1","2","3","4"))
counts$bpcount <- as.numeric(counts$bpcount)
p = ggplot(counts, aes_string(x="filecount", y="uniquecazy", colour="bpcount"))
p + geom_point(stat="identity", size=4) + ylab("Number of core sequences") + xlab("Total combined metagenomes") + theme_bw()+theme(text=element_text(size=10, family="Helvetica"))+theme(axis.text.x=element_text(size=15, angle=90),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15))+scale_colour_gradientn(colours=rainbow(2), "Cumulative bp")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
ggsave(file="Supp_Fig_1.eps")
##################################################################################

##################################################################################
#FIGURE 1
abundance_data <- read.delim(sep='\t', file="./summary.unfiltered.txt",header=TRUE, strip.white=TRUE, row.names=1)
ann_data <- read.delim(sep='\t', file="./annotations_cazy.txt",header=TRUE, strip.white=TRUE, row.names=1)
ann_only_euk <- subset(ann_data, t1 == "Eukaryota")
summary(ann_only_euk$t2)
ann_data_filtered <- subset(ann_data, t1 != "Eukaryota" | t2 == "Fungi" )
metadata = read.delim(file="./cobs_metadata.txt", row.names = 1, sep="\t", header=TRUE)
ann_data_matrix <- as.matrix(ann_data_filtered)
data_norm <- read.delim(sep='\t', file="./cumulative-all-normreca.txt",header=TRUE, strip.white=TRUE, row.names=1)
data_norm_core <- read.delim(sep='\t', file="./core-ws-recanorm.txt",header=TRUE, strip.white=TRUE, row.names=1)
ann_data_matrix <- as.matrix(ann_data_filtered)
annotation <- tax_table(ann_data_matrix)
abundance_data_norm_matrix <- as.matrix(data_norm)
abundance_data_norm_core_matrix <- as.matrix(data_norm_core)
abundance_data <- as.matrix(data_norm)
abundance <- otu_table(abundance_data, taxa_are_rows=TRUE)
abundance <- otu_table(abundance_data_norm_matrix, taxa_are_rows=TRUE)
abundance_core <- otu_table(abundance_data_norm_core_matrix, taxa_are_rows=TRUE)
metadata <- sample_data(metadata)
all_agg <- phyloseq(metadata, annotation, abundance)
all_agg_core <- phyloseq(metadata, annotation, abundance_core)
all_agg <- subset_samples(all_agg, agg_frac == "WS")
all_agg_core <- subset_samples(all_agg_core, agg_frac == "WS")

write.table(rownames(tax_table(all_agg_core)), file="core-filtered.txt", quote=FALSE, sep="\n", row.names=FALSE, col.names=FALSE)

mdf = psmelt(all_agg_core)
mdf <- subset(mdf, Cazy_fam != "none")
#mdf <- subset(mdf, Cazy_fam2 == "CE")
#f <- ddply(mdf, .(Cazy_fam, Cazy_fam2, sample_name), summarise, SUM=sum(Abundance)) <- this is wrong
f <- ddply(mdf, .(Cazy_fam2, sample_name), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(Cazy_fam2), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f2$Cazy_fam2 <- factor(f2$Cazy_fam2, levels=f2$Cazy_fam2[order(-f2$MEAN)])
p = ggplot(f2, aes_string(x="Cazy_fam2", y="MEAN"))
p = p + geom_bar(stat = "identity")+ geom_errorbar(limits, width=0) 
p = p  + xlab("") + ylab("Abundance (per recA gene)")+ theme_bw()  + theme(axis.text.x = element_text(angle=90, hjust=1))
p+theme(text=element_text(size=10, family="Helvetica"))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15))
ggsave(file="Fig_1a.eps")
##################################################################################
#FIGURE 1b
mdf = psmelt(all_agg_core)
mdf <- subset(mdf, Cazy_fam != "none")
#mdf <- subset(mdf, Cazy_fam2 == "CE")
f <- ddply(mdf, .(Cazy_fam, Cazy_fam2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(Cazy_fam2, Cazy_fam, agg_frac), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f2 <- data.frame(f2)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
foo <- aggregate(f2$MEAN, by=list(f2$Cazy_fam), max)
foo <- data.frame(foo)
names(foo) <- c("Cazy_fam", "Agg_sort")
foo2 <- join(f2, foo)
f2 <- foo2
o <- factor(f2$Cazy_fam, levels=f2$Cazy_fam[order(-f2$MEAN)][1:10])
mdf$agg_frac<-factor(mdf$agg_frac, levels=c("micro","SM","MM","LM","WS"))
f2$Cazy_fam <- factor(f2$Cazy_fam, levels=c(levels(o)))
f2 <- subset(f2, Cazy_fam != "NA")
p = ggplot(f2, aes_string(x="Cazy_fam", y="MEAN"))
p = p + geom_point(stat = "identity", size=5, aes(color=Cazy_fam2, fill=Cazy_fam2))+ geom_errorbar(limits, width=0) 
p = p  + xlab("") + ylab("Abundance (per recA)")+ theme_bw()  + theme(axis.text.x = element_text(angle=90, hjust=1))
p+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+scale_fill_manual("CAZy Enzyme Class", values=c("red","blue","dark green"))+scale_colour_manual("CAZy Enzyme Class", values=c("red","blue","dark green"))
ggsave(file="Fig_1b.eps")
##################################################################################

##################################################################################
#FIGURE 2
f <- ddply(mdf, .(sample_name, t2), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(t2), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f2 <- data.frame(f2)
f2$t2 <- trim(f2$t2)
f2 <- subset(f2, t2 != "")
f2 <- subset(f2, t2 != "dsDNA")
f2 <- subset(f2, t2 != "environmental")
f2 <- subset(f2, t2 != "candidate")
f2 <- subset(f2, Cazy_fam2 == "GH")
temp <- subset(f2, Cazy_fam2 == "GT")
#o2 <- factor(f2$t2, levels=f$t2[order(-f2$MEAN)])
#f2$t2 <- factor(f2$t2, levels=c(levels(o2)))
f2$t2=factor(f2$t2, levels=f2$t2[with(f2, order(-MEAN))][1:15])
f2 <- subset(f2, t2 != "NA")
p = ggplot(f2, aes_string(x="t2", y="MEAN"))
p = p + geom_bar(stat = "identity")
p = p + theme_bw()+ theme(text=element_text(size=20))+ylab("Abundance (per recA)") +geom_errorbar(limits, width=0) + theme(strip.text.x = element_text(angle = 90)) + theme(strip.text.y = element_text(angle = 0)) + theme(legend.position="none")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + xlab("")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p
ggsave(file="Fig_2.eps")
##################################################################################

##################################################################################
#SUPP FIG 2A
all_agg_filter_empty <- all_agg_core
mdf = psmelt(all_agg_filter_empty)
mdf <- subset(mdf, Cazy_fam2 != "none")
mdf <- subset(mdf, t2 != " " | t2 != 'dsDNA ' | t2 != "environmental " | t2 != "candidate ")
f <- ddply(mdf, .(Cazy_fam2,Cazy_fam, sample_name, t2), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(Cazy_fam2, t2), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f2 <- data.frame(f2)
f2$t2 <- trim(f2$t2)
f2 <- subset(f2, t2 != "")
f2 <- subset(f2, t2 != "dsDNA")
f2 <- subset(f2, t2 != "environmental")
f2 <- subset(f2, t2 != "candidate")

f2 <- subset(f2, Cazy_fam2 == "GH")
temp <- subset(f2, Cazy_fam2 == "GH")
#o2 <- factor(f2$t2, levels=f$t2[order(-f2$MEAN)])
#f2$t2 <- factor(f2$t2, levels=c(levels(o2)))
f2$Cazy_fam2 = factor(f2$Cazy_fam2, levels = c("GT", "GH", "CE", "CB", "PL"))
f2$t2=factor(f2$t2, levels=f2$t2[with(f2, order(-MEAN))][1:10])
#f2$t2=factor(f2$t2, levels=c("Proteobacteria","Bacteroidetes","Firmicutes","Euryarchaeota","Actinobacteria","Chloroflexi","Viridiplantae","Verrucomicrobia","Cyanobacteria","Spirochaetes","Crenarchaeota","Acidobacteria","Planctomycetes","Korarchaeota","Fungi","Chlorobi","Metazoa","Aquificae","Mycetozoa","Chlamydiae","Thermotogae","Deinococcus-Thermus","Alveolata"))
f2 <- subset(f2, t2 != "NA")
p = ggplot(f2, aes_string(x="t2", y="MEAN"))
p = p + geom_bar(stat = "identity")
p = p + theme_bw()+ ylab("Abundance (per recA)") +geom_errorbar(limits, width=0) + facet_grid(~Cazy_fam2) +  theme(text=element_text(size=18), strip.text.x = element_text(angle = 90)) + theme(strip.text.y = element_text(angle = 0)) + theme(legend.position="none")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + xlab("")+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
p
ggsave(file="Supp_Fig_3a.eps")

#SUPP FIG 3B
f <- ddply(mdf, .(Cazy_fam2,Cazy_fam, sample_name, t2), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(Cazy_fam2, t2), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f2 <- data.frame(f2)
f2$t2 <- trim(f2$t2)
f2 <- subset(f2, t2 != "")
f2 <- subset(f2, t2 != "dsDNA")
f2 <- subset(f2, t2 != "environmental")
f2 <- subset(f2, t2 != "candidate")

f2 <- subset(f2, Cazy_fam2 == "GT")
temp <- subset(f2, Cazy_fam2 == "GT")
#o2 <- factor(f2$t2, levels=f$t2[order(-f2$MEAN)])
#f2$t2 <- factor(f2$t2, levels=c(levels(o2)))
f2$Cazy_fam2 = factor(f2$Cazy_fam2, levels = c("GT", "GH", "CE", "CB", "PL"))
f2$t2=factor(f2$t2, levels=f2$t2[with(f2, order(-MEAN))][1:10])
#f2$t2=factor(f2$t2, levels=c("Proteobacteria","Bacteroidetes","Firmicutes","Euryarchaeota","Actinobacteria","Chloroflexi","Viridiplantae","Verrucomicrobia","Cyanobacteria","Spirochaetes","Crenarchaeota","Acidobacteria","Planctomycetes","Korarchaeota","Fungi","Chlorobi","Metazoa","Aquificae","Mycetozoa","Chlamydiae","Thermotogae","Deinococcus-Thermus","Alveolata"))
f2 <- subset(f2, t2 != "NA")
p = ggplot(f2, aes_string(x="t2", y="MEAN"))
p = p + geom_bar(stat = "identity")
p = p + theme_bw()+ ylab("Abundance (per recA)") +geom_errorbar(limits, width=0) + facet_grid(~Cazy_fam2) +  theme(text=element_text(size=18), strip.text.x = element_text(angle = 90)) + theme(strip.text.y = element_text(angle = 0)) + theme(legend.position="none")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + xlab("")+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
p
ggsave(file="Supp_Fig_3b.eps")

#SUPP FIG 3C
f <- ddply(mdf, .(Cazy_fam2,Cazy_fam, sample_name, t2), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(Cazy_fam2, t2), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f2 <- data.frame(f2)
f2$t2 <- trim(f2$t2)
f2 <- subset(f2, t2 != "")
f2 <- subset(f2, t2 != "dsDNA")
f2 <- subset(f2, t2 != "environmental")
f2 <- subset(f2, t2 != "candidate")

f2 <- subset(f2, Cazy_fam2 == "CB")
temp <- subset(f2, Cazy_fam2 == "CB")
#o2 <- factor(f2$t2, levels=f$t2[order(-f2$MEAN)])
#f2$t2 <- factor(f2$t2, levels=c(levels(o2)))
f2$Cazy_fam2 = factor(f2$Cazy_fam2, levels = c("GT", "GH", "CE", "CB", "PL"))
f2$t2=factor(f2$t2, levels=f2$t2[with(f2, order(-MEAN))][1:10])
#f2$t2=factor(f2$t2, levels=c("Proteobacteria","Bacteroidetes","Firmicutes","Euryarchaeota","Actinobacteria","Chloroflexi","Viridiplantae","Verrucomicrobia","Cyanobacteria","Spirochaetes","Crenarchaeota","Acidobacteria","Planctomycetes","Korarchaeota","Fungi","Chlorobi","Metazoa","Aquificae","Mycetozoa","Chlamydiae","Thermotogae","Deinococcus-Thermus","Alveolata"))
f2 <- subset(f2, t2 != "NA")
p = ggplot(f2, aes_string(x="t2", y="MEAN"))
p = p + geom_bar(stat = "identity")
p = p + theme_bw()+ ylab("Abundance (per recA)") +geom_errorbar(limits, width=0) + facet_grid(~Cazy_fam2) +  theme(text=element_text(size=18), strip.text.x = element_text(angle = 90)) + theme(strip.text.y = element_text(angle = 0)) + theme(legend.position="none")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + xlab("")+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
p
ggsave(file="Supp_Fig_3c.eps")

#SUPP FIG 3D
f <- ddply(mdf, .(Cazy_fam2,Cazy_fam, sample_name, t2), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(Cazy_fam2, t2), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f2 <- data.frame(f2)
f2$t2 <- trim(f2$t2)
f2 <- subset(f2, t2 != "")
f2 <- subset(f2, t2 != "dsDNA")
f2 <- subset(f2, t2 != "environmental")
f2 <- subset(f2, t2 != "candidate")

f2 <- subset(f2, Cazy_fam2 == "CE")
temp <- subset(f2, Cazy_fam2 == "CE")
#o2 <- factor(f2$t2, levels=f$t2[order(-f2$MEAN)])
#f2$t2 <- factor(f2$t2, levels=c(levels(o2)))
f2$Cazy_fam2 = factor(f2$Cazy_fam2, levels = c("GT", "GH", "CE", "CB", "PL"))
f2$t2=factor(f2$t2, levels=f2$t2[with(f2, order(-MEAN))][1:10])
#f2$t2=factor(f2$t2, levels=c("Proteobacteria","Bacteroidetes","Firmicutes","Euryarchaeota","Actinobacteria","Chloroflexi","Viridiplantae","Verrucomicrobia","Cyanobacteria","Spirochaetes","Crenarchaeota","Acidobacteria","Planctomycetes","Korarchaeota","Fungi","Chlorobi","Metazoa","Aquificae","Mycetozoa","Chlamydiae","Thermotogae","Deinococcus-Thermus","Alveolata"))
f2 <- subset(f2, t2 != "NA")
p = ggplot(f2, aes_string(x="t2", y="MEAN"))
p = p + geom_bar(stat = "identity")
p = p + theme_bw()+ ylab("Abundance (per recA)") +geom_errorbar(limits, width=0) + facet_grid(~Cazy_fam2) +  theme(text=element_text(size=18), strip.text.x = element_text(angle = 90)) + theme(strip.text.y = element_text(angle = 0)) + theme(legend.position="none")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + xlab("")+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
p
ggsave(file="Supp_Fig_3d.eps")
##################################################################################




##################################################################################
#Supp FIGURE 4
mdf_core<-psmelt(all_agg_core)
mdf_noncore<-psmelt(all_agg)
f <- ddply(mdf_core, .(Cazy_fam2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(Cazy_fam2, agg_frac), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
g <- ddply(mdf_noncore, .(Cazy_fam2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
g2 <- ddply(g, .(Cazy_fam2, agg_frac), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
smple_wgt<-aggregate(f$SUM,by=list(f$sample_name),FUN=sum)

k<-1
f$SUMs<-0
while(k<=dim(smple_wgt)[1]) {
	smr<-which(f$sample_name==smple_wgt$Group.1[k])
	f$SUMs[smr]<-f$SUM[smr]/smple_wgt[k,2]
	k<-k+1
}

stderr <- function(x){
	result<-sd(x)/sqrt(length(x))
	return(result)}
	
fam_wgt<-aggregate(f$SUMs,by=list(f$Cazy_fam2),FUN=mean)
fam_wgt2<-aggregate(f$SUMs,by=list(f$Cazy_fam2),FUN=stderr)
names(fam_wgt)[2]<-"mean"
names(fam_wgt2)[2]<-"std_er"
fam_stat<-join(fam_wgt,fam_wgt2)
fam_stat$dataset <- "core"
smple2_wgt<-aggregate(g$SUM,by=list(g$sample_name),FUN=sum)
k<-1
g$SUMs<-0
while(k<=dim(smple_wgt)[1]) {
	smr<-which(g$sample_name==smple2_wgt$Group.1[k])
	g$SUMs[smr]<-g$SUM[smr]/smple2_wgt[k,2]
	k<-k+1
}

fam2_wgt<-aggregate(g$SUMs,by=list(g$Cazy_fam2),FUN=mean)
fam2_wgt2<-aggregate(g$SUMs,by=list(g$Cazy_fam2),FUN=stderr)
names(fam2_wgt)[2]<-"mean"
names(fam2_wgt2)[2]<-"std_er"
fam2_stat<-join(fam2_wgt,fam2_wgt2)
fam2_stat$dataset <- "noncore"

foo <- rbind(fam_stat, fam2_stat)
limits<-aes(ymin=mean-std_er, ymax=mean+std_er)
colnames(foo)[1] <- "Cazy_fam2"
foo <- subset(foo, Cazy_fam2 != "none")
foo <- subset(foo, Cazy_fam2 != "NA")
temp <- subset(foo, dataset=="noncore")
foo$Cazy_fam2 <- factor(foo$Cazy_fam2, levels=temp$Cazy_fam2[with(temp,order(-mean))][1:4])
foo <- subset(foo, Cazy_fam2 != "NA")
p = ggplot(foo, aes_string(x="Cazy_fam2", y="mean", colour="dataset", fill="dataset"))
p = p + geom_bar(stat="identity",position="dodge") + geom_errorbar(limits, width=0, position=position_dodge(.9)) 
p + theme_bw() + theme(text=element_text(size=20))+theme(axis.text.x=element_text(size=12, angle=90, hjust=1))+ylab("Relative Abundance")+xlab("")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+scale_fill_manual("Dataset",values=c("black","white"),labels=c("core","total"))+scale_colour_manual("Dataset",values=c("black","black"), labels=c("core","total"))+annotate("text",x=1,y=.47,label="***") + annotate("text",x=4, y=.11,label="***")
ggsave(file="Supp_Fig_5a.eps")

#Stats
f$dataset <- "core"
g$dataset <- "noncore"
combo <- rbind(f,g)
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam2 == "CB")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam2 == "GT")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam2 == "CE")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam2 == "GH")))

#Supp FIGURE 4B
f <- ddply(mdf_core, .(Cazy_fam, Cazy_fam2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(Cazy_fam, Cazy_fam2, agg_frac), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
g <- ddply(mdf_noncore, .(Cazy_fam, Cazy_fam2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
g2 <- ddply(g, .(Cazy_fam, Cazy_fam2, agg_frac), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
smple_wgt<-aggregate(f$SUM,by=list(f$sample_name),FUN=sum)

k<-1
f$SUMs<-0
while(k<=dim(smple_wgt)[1]) {
	smr<-which(f$sample_name==smple_wgt$Group.1[k])
	f$SUMs[smr]<-f$SUM[smr]/smple_wgt[k,2]
	k<-k+1
}

stderr <- function(x){
	result<-sd(x)/sqrt(length(x))
	return(result)}
	
fam_wgt<-aggregate(f$SUMs,by=list(f$Cazy_fam),FUN=mean)
fam_wgt2<-aggregate(f$SUMs,by=list(f$Cazy_fam),FUN=stderr)
names(fam_wgt)[2]<-"mean"
names(fam_wgt2)[2]<-"std_er"
fam_stat<-join(fam_wgt,fam_wgt2)
fam_stat$dataset <- "core"
smple2_wgt<-aggregate(g$SUM,by=list(g$sample_name),FUN=sum)
k<-1
g$SUMs<-0
while(k<=dim(smple_wgt)[1]) {
	smr<-which(g$sample_name==smple2_wgt$Group.1[k])
	g$SUMs[smr]<-g$SUM[smr]/smple2_wgt[k,2]
	k<-k+1
}

fam2_wgt<-aggregate(g$SUMs,by=list(g$Cazy_fam),FUN=mean)
fam2_wgt2<-aggregate(g$SUMs,by=list(g$Cazy_fam),FUN=stderr)
names(fam2_wgt)[2]<-"mean"
names(fam2_wgt2)[2]<-"std_er"
fam2_stat<-join(fam2_wgt,fam2_wgt2)
fam2_stat$dataset <- "noncore"

foo <- rbind(fam_stat, fam2_stat)
limits<-aes(ymin=mean-std_er, ymax=mean+std_er)
colnames(foo)[1] <- "Cazy_fam"
foo <- subset(foo, Cazy_fam != "none")
foo <- subset(foo, Cazy_fam != "NA")
temp <- subset(foo, dataset=="noncore")
foo$Cazy_fam <- factor(foo$Cazy_fam, levels=temp$Cazy_fam[with(temp,order(-mean))][1:10])
foo <- subset(foo, Cazy_fam != "NA")
p = ggplot(foo, aes_string(x="Cazy_fam", y="mean", colour="dataset", fill="dataset"))
p = p + geom_bar(stat="identity",position="dodge") + geom_errorbar(limits, width=0, position=position_dodge(.9)) 
p + theme_bw() + theme(text=element_text(size=20))+theme(axis.text.x=element_text(size=12, angle=90, hjust=1))+ylab("Relative Abundance")+xlab("")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+scale_fill_manual("Dataset",values=c("black","white"),labels=c("core","total"))+scale_colour_manual("Dataset",values=c("black","black"), labels=c("core","total"))+annotate("text",x=2,y=.13,label="***") + annotate("text",x=3,y=.05,label="***") + annotate("text",x=4,y=.06,label="***")+annotate("text",x=5,y=.04,label="**") + annotate("text",x=7,y=.05,label="***") + annotate("text",x=8,y=.03,label="***")+annotate("text",x=9,y=.03,label="***") +annotate("text",x=10,y=.03,label="**") 

ggsave(file="Supp_Fig_5b.eps")

f$dataset <- "core"
g$dataset <- "noncore"
combo <- rbind(f,g)
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam == "GT2")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam == "GT4")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam == "GT41")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam == "GH13")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam == "GH3")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam == "CE10")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam == "CE4")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam == "CE1")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam == "GT51")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam == "GH23")))





##################################################################################
#SUPP FIG 3A


mdf_core<-psmelt(all_agg_core)
mdf_noncore<-psmelt(all_agg)
f <- ddply(mdf_core, .(Cazy_fam,Cazy_fam2, t2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
f2 <- subset(f, Cazy_fam == "GT4")
g <- ddply(mdf_noncore, .(Cazy_fam, Cazy_fam2, t2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
g2 <- subset(g, Cazy_fam == "GT4")
smple_wgt<-aggregate(f$SUM,by=list(f$sample_name),FUN=sum)
smple2_wgt<-aggregate(f2$SUM,by=list(f2$sample_name),FUN=sum)

k<-1
f2$SUMs<-0
while(k<=dim(smple_wgt)[1]) {
	smr<-which(f2$sample_name==smple2_wgt$Group.1[k])
	f2$SUMs[smr]<-f2$SUM[smr]/smple2_wgt[k,2]
	k<-k+1
}

stderr <- function(x){
	result<-sd(x)/sqrt(length(x))
	return(result)}
	
fam_wgt<-aggregate(f2$SUMs,by=list(f2$t2),FUN=mean)
fam_wgt2<-aggregate(f2$SUMs,by=list(f2$t2),FUN=stderr)
names(fam_wgt)[2]<-"mean"
names(fam_wgt2)[2]<-"std_er"
fam_stat<-join(fam_wgt,fam_wgt2)
fam_stat$dataset <- "core"

smple2_wgt<-aggregate(g2$SUM,by=list(g2$sample_name),FUN=sum)

k<-1
g2$SUMs<-0
while(k<=dim(smple_wgt)[1]) {
	smr<-which(g2$sample_name==smple2_wgt$Group.1[k])
	g2$SUMs[smr]<-g2$SUM[smr]/smple2_wgt[k,2]
	k<-k+1
}

fam2_wgt<-aggregate(g2$SUMs,by=list(g2$t2),FUN=mean)
fam2_wgt2<-aggregate(g2$SUMs,by=list(g2$t2),FUN=stderr)
names(fam2_wgt)[2]<-"mean"
names(fam2_wgt2)[2]<-"std_er"
fam2_stat<-join(fam2_wgt,fam2_wgt2)
fam2_stat$dataset <- "noncore"


foo <- rbind(fam_stat, fam2_stat)
limits<-aes(ymin=mean-std_er, ymax=mean+std_er)
colnames(foo)[1] <- "t2"
foo <- subset(foo, t2 != "none")
temp <- subset(foo, dataset=="noncore")
foo$t2 <- factor(foo$t2, levels=temp$t2[with(temp,order(-mean))][1:9])
foo <- subset(foo, t2 != "NA")
foo$Cazy_fam2 <- "GT4"
p = ggplot(foo, aes_string(x="t2", y="mean", colour="dataset", fill="dataset"))
p = p + geom_bar(stat="identity",position="dodge") + geom_errorbar(limits, width=0, position=position_dodge(.9)) 
p + theme_bw() + theme(text=element_text(size=20))+theme(axis.text.x=element_text(size=12, angle=90, hjust=1))+ylab("Relative Abundance")+xlab("")+facet_grid(~Cazy_fam2)+theme(text=element_text(size=18), strip.text.x = element_text(angle = 90))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+scale_fill_manual("Dataset",values=c("black","white"),labels=c("core","total"))+scale_colour_manual("Dataset",values=c("black","black"), labels=c("core","total")) +annotate("text",x=2,y=.13,label="***") + annotate("text",x=4,y=.11,label="**") + annotate("text",x=5,y=.08,label="***")+annotate("text",x=6,y=.08,label="***") + annotate("text",x=7,y=.11,label="**") + annotate("text",x=8,y=.04,label="**")+annotate("text",x=9,y=.11,label="***") 
ggsave(file="Supp_Fig_6a.eps")

l = unique(foo$t2)
length_l = length(l)
for (i in 1:length_l){
	f_dat <- subset(f2, t2 == as.character(l[i]))
	if(dim(f_dat)[1] == 0){
		print(as.character(l[i]))
		print("error")} else{
	g_dat <- subset(g2, t2 == as.character(l[i]))
	f_dat$dataset <- "core"
	g_dat$dataset <- "noncore"
	combined <- rbind(f_dat, g_dat)
	print(as.character(l[i]))
	print(summary(aov.result <- aov(SUMs~dataset, data = combined)))
	}}
##################################################################################
#SUPP FIGURE 3B

f <- ddply(mdf_core, .(Cazy_fam,Cazy_fam2, t2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
f2 <- subset(f, Cazy_fam == "GH13")
g <- ddply(mdf_noncore, .(Cazy_fam, Cazy_fam2, t2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
g2 <- subset(g, Cazy_fam == "GH13")

smple2_wgt<-aggregate(f2$SUM,by=list(f2$sample_name),FUN=sum)

k<-1
f2$SUMs<-0
while(k<=dim(smple_wgt)[1]) {
	smr<-which(f2$sample_name==smple2_wgt$Group.1[k])
	f2$SUMs[smr]<-f2$SUM[smr]/smple2_wgt[k,2]
	k<-k+1
}

stderr <- function(x){
	result<-sd(x)/sqrt(length(x))
	return(result)}
	
fam_wgt<-aggregate(f2$SUMs,by=list(f2$t2),FUN=mean)
fam_wgt2<-aggregate(f2$SUMs,by=list(f2$t2),FUN=stderr)
names(fam_wgt)[2]<-"mean"
names(fam_wgt2)[2]<-"std_er"
fam_stat<-join(fam_wgt,fam_wgt2)
fam_stat$dataset <- "core"

smple2_wgt<-aggregate(g2$SUM,by=list(g2$sample_name),FUN=sum)

k<-1
g2$SUMs<-0
while(k<=dim(smple_wgt)[1]) {
	smr<-which(g2$sample_name==smple2_wgt$Group.1[k])
	g2$SUMs[smr]<-g2$SUM[smr]/smple2_wgt[k,2]
	k<-k+1
}

fam2_wgt<-aggregate(g2$SUMs,by=list(g2$t2),FUN=mean)
fam2_wgt2<-aggregate(g2$SUMs,by=list(g2$t2),FUN=stderr)
names(fam2_wgt)[2]<-"mean"
names(fam2_wgt2)[2]<-"std_er"
fam2_stat<-join(fam2_wgt,fam2_wgt2)
fam2_stat$dataset <- "noncore"

foo <- rbind(fam_stat, fam2_stat)
limits<-aes(ymin=mean-std_er, ymax=mean+std_er)
colnames(foo)[1] <- "t2"
foo <- subset(foo, t2 != "none")
temp <- subset(foo, dataset=="noncore")
foo$t2 <- factor(foo$t2, levels=temp$t2[with(temp,order(-mean))][1:9])
foo <- subset(foo, t2 != "NA")
foo_l<-dim(foo)[1]
foo[foo_l+1,]<-0
foo$t2[foo_l+1]<-"Bacteroidetes"
foo$dataset[foo_l+1]<-"core"
foo$Cazy_fam2 <- "GH13"
p = ggplot(foo, aes_string(x="t2", y="mean", colour="dataset", fill="dataset"))
p = p + geom_bar(stat="identity",position="dodge") + geom_errorbar(limits, width=0, position=position_dodge(.9)) 
p + theme_bw() + theme(text=element_text(size=20))+theme(axis.text.x=element_text(size=12, angle=90, hjust=1))+ylab("Relative Abundance")+xlab("")+facet_grid(~Cazy_fam2)+theme(text=element_text(size=18), strip.text.x = element_text(angle = 90))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+scale_fill_manual("Dataset",values=c("black","white"),labels=c("core","total"))+scale_colour_manual("Dataset",values=c("black","black"), labels=c("core","total"))+annotate("text",x=1,y=.43,label="***") + annotate("text",x=2,y=.33,label="***") + annotate("text",x=3,y=.1,label="***")+annotate("text",x=5,y=.1,label="*") + annotate("text",x=6,y=.1,label="***") + annotate("text",x=7,y=.1,label="***")+annotate("text",x=9,y=.05,label="*") 
ggsave(file="Supp_Fig_6b.eps")

l = unique(foo$t2)
length_l = length(l)
for (i in 1:length_l){
	f_dat <- subset(f2, t2 == as.character(l[i]))
	if(dim(f_dat)[1] == 0){
		print(as.character(l[i]))
		print("error")} else{
	g_dat <- subset(g2, t2 == as.character(l[i]))
	f_dat$dataset <- "core"
	g_dat$dataset <- "noncore"
	combined <- rbind(f_dat, g_dat)
	print(as.character(l[i]))
	print(summary(aov.result <- aov(SUMs~dataset, data = combined)))
	}}
##################################################################################
#SUPP FIGURE 3C

f <- ddply(mdf_core, .(Cazy_fam,Cazy_fam2, t2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
f2 <- subset(f, Cazy_fam == "CE10")
g <- ddply(mdf_noncore, .(Cazy_fam, Cazy_fam2, t2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
g2 <- subset(g, Cazy_fam == "CE10")

smple2_wgt<-aggregate(f2$SUM,by=list(f2$sample_name),FUN=sum)

k<-1
f2$SUMs<-0
while(k<=dim(smple_wgt)[1]) {
	smr<-which(f2$sample_name==smple2_wgt$Group.1[k])
	f2$SUMs[smr]<-f2$SUM[smr]/smple2_wgt[k,2]
	k<-k+1
}

stderr <- function(x){
	result<-sd(x)/sqrt(length(x))
	return(result)}
	
fam_wgt<-aggregate(f2$SUMs,by=list(f2$t2),FUN=mean)
fam_wgt2<-aggregate(f2$SUMs,by=list(f2$t2),FUN=stderr)
names(fam_wgt)[2]<-"mean"
names(fam_wgt2)[2]<-"std_er"
fam_stat<-join(fam_wgt,fam_wgt2)
fam_stat$dataset <- "core"

smple2_wgt<-aggregate(g2$SUM,by=list(g2$sample_name),FUN=sum)

k<-1
g2$SUMs<-0
while(k<=dim(smple_wgt)[1]) {
	smr<-which(g2$sample_name==smple2_wgt$Group.1[k])
	g2$SUMs[smr]<-g2$SUM[smr]/smple2_wgt[k,2]
	k<-k+1
}

fam2_wgt<-aggregate(g2$SUMs,by=list(g2$t2),FUN=mean)
fam2_wgt2<-aggregate(g2$SUMs,by=list(g2$t2),FUN=stderr)
names(fam2_wgt)[2]<-"mean"
names(fam2_wgt2)[2]<-"std_er"
fam2_stat<-join(fam2_wgt,fam2_wgt2)
fam2_stat$dataset <- "noncore"

foo <- rbind(fam_stat, fam2_stat)
limits<-aes(ymin=mean-std_er, ymax=mean+std_er)
colnames(foo)[1] <- "t2"
foo <- subset(foo, t2 != "none")
temp <- subset(foo, dataset=="noncore")
foo$t2 <- factor(foo$t2, levels=temp$t2[with(temp,order(-mean))][1:9])
foo <- subset(foo, t2 != "NA")
foo_l<-dim(foo)[1]
foo[foo_l+1,]<-0
foo$t2[foo_l+1]<-"Cyanobacteria"
foo$dataset[foo_l+1]<-"core"
foo_l<-dim(foo)[1]
foo[foo_l+1,]<-0
foo$t2[foo_l+1]<-"Firmicutes"
foo$dataset[foo_l+1]<-"core"
foo <- subset(foo, t2 != "")
foo$Cazy_fam2 <- "CE10"
p = ggplot(foo, aes_string(x="t2", y="mean", colour="dataset", fill="dataset"))
p = p + geom_bar(stat="identity",position="dodge") + geom_errorbar(limits, width=0, position=position_dodge(.9)) 
p + theme_bw() + theme(text=element_text(size=20))+theme(axis.text.x=element_text(size=12, angle=90, hjust=1))+ylab("Relative Abundance")+xlab("")+facet_grid(~Cazy_fam2)+theme(text=element_text(size=18), strip.text.x = element_text(angle = 90))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+scale_fill_manual("Dataset",values=c("black","white"),labels=c("core","total"))+scale_colour_manual("Dataset",values=c("black","black"), labels=c("core","total"))+annotate("text",x=1,y=.55,label="***") + annotate("text",x=4,y=.1,label="**") + annotate("text",x=8,y=.1,label="*") 
ggsave(file="Supp_Fig_6c.eps")

##################################################################################
#SUPP FIGURE 3D
mdf_core<-psmelt(all_agg_core)
mdf_noncore<-psmelt(all_agg)
f <- ddply(mdf_core, .(Cazy_fam,Cazy_fam2, t2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
f2 <- subset(f, Cazy_fam == "GT2")
g <- ddply(mdf_noncore, .(Cazy_fam, Cazy_fam2, t2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
g2 <- subset(g, Cazy_fam == "GT2")
smple_wgt<-aggregate(f$SUM,by=list(f$sample_name),FUN=sum)
smple2_wgt<-aggregate(f2$SUM,by=list(f2$sample_name),FUN=sum)


k<-1
f2$SUMs<-0
while(k<=dim(smple_wgt)[1]) {
	smr<-which(f2$sample_name==smple2_wgt$Group.1[k])
	f2$SUMs[smr]<-f2$SUM[smr]/smple2_wgt[k,2]
	k<-k+1
}

stderr <- function(x){
	result<-sd(x)/sqrt(length(x))
	return(result)}
	
fam_wgt<-aggregate(f2$SUMs,by=list(f2$t2),FUN=mean)
fam_wgt2<-aggregate(f2$SUMs,by=list(f2$t2),FUN=stderr)
names(fam_wgt)[2]<-"mean"
names(fam_wgt2)[2]<-"std_er"
fam_stat<-join(fam_wgt,fam_wgt2)
fam_stat$dataset <- "core"

smple2_wgt<-aggregate(g2$SUM,by=list(g2$sample_name),FUN=sum)

k<-1
g2$SUMs<-0
while(k<=dim(smple_wgt)[1]) {
	smr<-which(g2$sample_name==smple2_wgt$Group.1[k])
	g2$SUMs[smr]<-g2$SUM[smr]/smple2_wgt[k,2]
	k<-k+1
}

fam2_wgt<-aggregate(g2$SUMs,by=list(g2$t2),FUN=mean)
fam2_wgt2<-aggregate(g2$SUMs,by=list(g2$t2),FUN=stderr)
names(fam2_wgt)[2]<-"mean"
names(fam2_wgt2)[2]<-"std_er"
fam2_stat<-join(fam2_wgt,fam2_wgt2)
fam2_stat$dataset <- "noncore"

foo <- rbind(fam_stat, fam2_stat)
limits<-aes(ymin=mean-std_er, ymax=mean+std_er)
colnames(foo)[1] <- "t2"
foo <- subset(foo, foo$t2 == "Proteobacteria" | foo$t2 == "Bacteroidetes" | foo$t2 == "Actinobacteria" | foo$t2 == "Cyanobacteria" | foo$t2 == "Firmicutes"| foo$t2 == "Chloroflexi" | foo$t2 == "Acidobacteria" | foo$t2 == "Euryarchaeota" | foo$t2 == "Crenarchaeota")
temp <- subset(foo, dataset=="noncore")
foo$t2 <- factor(foo$t2, levels=temp$t2[with(temp,order(-mean))][1:9])
foo$Cazy_fam2 <- "GT2"
p = ggplot(foo, aes_string(x="t2", y="mean", colour="dataset", fill="dataset"))
p = p + geom_bar(stat="identity",position="dodge") + geom_errorbar(limits, width=0, position=position_dodge(.9)) 
p + theme_bw() + theme(text=element_text(size=20))+theme(axis.text.x=element_text(size=12, angle=90, hjust=1))+ylab("Relative Abundance")+xlab("")+facet_grid(~Cazy_fam2)+theme(text=element_text(size=18), strip.text.x = element_text(angle = 90))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+scale_fill_manual("Dataset",values=c("black","white"),labels=c("core","total"))+scale_colour_manual("Dataset",values=c("black","black"), labels=c("core","total"))+annotate("text",x=3,y=.2,label="**") + annotate("text",x=4,y=.12,label="***") + annotate("text",x=5,y=.08,label="***")+ annotate("text",x=6,y=.06,label="***") + annotate("text",x=8,y=.04,label="**")
ggsave(file="Supp_Fig_6d.eps")

l = unique(foo$t2)
length_l = length(l)
for (i in 1:length_l){
	f_dat <- subset(f2, t2 == as.character(l[i]))
	if(dim(f_dat)[1] == 0){
		print(as.character(l[i]))
		print("error")} else{
	g_dat <- subset(g2, t2 == as.character(l[i]))
	f_dat$dataset <- "core"
	g_dat$dataset <- "noncore"
	combined <- rbind(f_dat, g_dat)
	print(as.character(l[i]))
	print(summary(aov.result <- aov(SUMs~dataset, data = combined)))
	}}

#exploring differences in deeper
#In GT2, Actinobacteria and & Firmicutes are enriched in core
cazy_string = "GH13"
taxa_string = "Planctomycetes"
mdf_core_subset <- subset(mdf_core, mdf_core$Cazy_fam == cazy_string & mdf_core$t2 == taxa_string)
mdf_noncore_subset <- subset(mdf_noncore, mdf_noncore$Cazy_fam == cazy_string & mdf_noncore$t2 == taxa_string)
summary(mdf_core_subset)
summary(mdf_noncore_subset)
cat(intersect(mdf_core_subset$t3, mdf_noncore_subset$t3),sep='\n')
setdiff(mdf_core_subset$t3, mdf_noncore_subset$t3)
setdiff(mdf_noncore_subset$t3, mdf_core_subset$t3)
cat(intersect(mdf_core_subset$hit_id, mdf_noncore_subset$hit_id), sep="\n")
intersect(mdf_core_subset$description, mdf_noncore_subset$description)
setdiff(mdf_core_subset$hit_id, mdf_noncore_subset$hit_id)
length(setdiff(mdf_noncore_subset$hit_id, mdf_core_subset$hit_id))
##################################################################################
#Figure 5

library(reshape)
hist_data <- read.delim(sep='\t', file="./histogram_core.txt",header=FALSE, strip.white=TRUE)
colnames(hist_data)<-c('sample','soil_type', 'mgrastids','value')
mdf<-melt(hist_data)
mdf$sample<-factor(mdf$sample, levels=mdf$sample[order(-mdf$value)])
p = ggplot(mdf, aes_string(x="sample", y="value"))
p = p + geom_bar(stat = "identity")
p + theme_bw()+ ylab("Shared Core Sequences")  + theme(text=element_text(size=14)) + theme(legend.position="none")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), axis.text.y=element_text(angle=90)) + xlab("")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
ggsave(file="Fig_3.eps")
##################################################################################

##################################################################################
#Figure 8

data  <- read.delim(sep='\t', file="./blast-summary-against-cazy.txt",header=TRUE, strip.white=TRUE, row.names=1)
data$X.1 <- NULL
data$core.e.5.cazyfam <- NULL
data <- data[c(3,4,5,7,8),]
data$aggregates.combined.assembly1.x.cazy.blast8out.besthits.5.cazyfam <- NULL
d2 <- data/rep(colSums(data),each=length(data$X4477904.3.x.cazy.blastxout.best.e.5.cazyfam))
d2$Cazy_fam <- rownames(d2)

mdf = psmelt(all_agg_core)
mdf <- subset(mdf, Cazy_fam != "none")
#mdf <- subset(mdf, Cazy_fam2 == "CE")
f <- ddply(mdf, .(Cazy_fam2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(Cazy_fam2, agg_frac), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f2$agg_frac<-NULL
f2$SE<-NULL
f2$core <- f2$MEAN
f2$MEAN <- NULL
f2$core <- f2$core / sum(f2$core)
rownames(f2)<-f2$Cazy_fam2
f2$Cazy_fam2 <- NULL

mdf2 = psmelt(all_agg)
mdf2 <- subset(mdf2, Cazy_fam != "none")
#mdf <- subset(mdf, Cazy_fam2 == "CE")
f3 <- ddply(mdf2, .(Cazy_fam2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
f4 <- ddply(f3, .(Cazy_fam2, agg_frac), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f4$agg_frac<-NULL
f4$SE<-NULL
f4$ws_cum <- f4$MEAN
f4$MEAN <- NULL
f4$ws_cum <- f4$ws_cum / sum(f4$ws_cum)
f4 <- subset(f4, Cazy_fam2 != "LD")
f4 <- subset(f4, Cazy_fam2 != "LO")
rownames(f4)<- f4$Cazy_fam2
f4$Cazy_fam2 <- NULL
t2<-merge(f2, f4, by="row.names")
rownames(t2)<-t2$Row.names
t2$Row.names<-NULL
t2$Cazy_fam <- NULL

new_t <- merge(d2, t2, by="row.names")
rownames(new_t) <- new_t$Cazy_fam
new_t$Row.names <- NULL
new_t$Cazy_fam <- NULL

colnames(d2)
colnames(new_t)<- c("Polar desert EB026","Temp Grassland KP1", "Hot desert MD3", "Tropical Forest PE6","Hot desert SF2","Hot desert SV1","Arctic tundra TL1","Tropical Forest AR3","Boreal forest BZ1","Temp dec forest CL1","Temp con forest DF1","Polar desert EB017","Polar desert EB019","Polar desert EB020","Polar desert EB021","Polar desert EB024","Fert Prairie Core","Fert Prairie Cumul.")

new_t$Cazy_fam <- rownames(new_t)
mdf_melt <- melt.data.frame(new_t)
mdf_melt$variable <- factor(mdf_melt$variable, levels=c("Fert Prairie Core","Fert Prairie Cumul.","Temp Grassland KP1",  "Tropical Forest PE6","Tropical Forest AR3","Boreal forest BZ1","Temp dec forest CL1","Temp con forest DF1","Hot desert MD3","Hot desert SF2","Hot desert SV1","Polar desert EB026","Polar desert EB017","Polar desert EB019","Polar desert EB020","Polar desert EB021","Polar desert EB024","Arctic tundra TL1"))

p = ggplot(mdf_melt, aes_string(x="variable", y="value", color="Cazy_fam", fill="Cazy_fam"))
p + geom_bar(stat="identity")+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(axis.text.x=element_text(size=15, angle=90),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+xlab("Soil metagenome")+ylab("Relative Abundance")+scale_fill_manual("Cazy Enzyme Class",values=c("red","blue","dark green","orange","brown"))+scale_colour_manual("Cazy Enzyme Class",values=c("red","blue","dark green","orange","brown"))
ggsave(file="Fig_6.eps")







##################################################################################
#Supp Figure 4



f <- ddply(mdf_core, .(EC, sample_name, agg_frac), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(EC, agg_frac), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f2 <- data.frame(f2)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f3 <- subset(f2, EC != "n.d.")
f3 <- subset(f3, EC != "2.4.1.-")
f3 <- subset(f3, EC != "3.2.1.-")
f3$EC <- factor(f3$EC,  levels=f3$EC[order(-f3$MEAN)])
a1<-f3[which(f3$EC=="3.2.1.20"),]
#a2<-f3[which(f3$EC=="3.2.1.21"),]
#a3<-f3[which(f3$EC=="3.2.1.31"),] 
#a4<-f3[which(f3$EC=="3.2.1.23"),]
a5<-f3[which(f3$EC=="3.2.1.25"),]
#a6<-f3[which(f3$EC=="3.2.1.55"),]
#a7<-f3[which(f3$EC=="3.2.1.91"),] 
a8<-f3[which(f3$EC=="3.2.1.37"),]
#a9<-f3[which(f3$EC=="3.2.1.14"),]
#a10<-f3[which(f3$EC=="3.4.11.1"),]
#a11<-f3[which(f3$EC=="3.1.3.1"),]
alla<-rbind(a1,a5,a8)


p = ggplot(f3, aes_string(x="EC", y="MEAN"))
p = p + geom_point(stat = "identity")+ geom_errorbar(limits, width=0)+geom_point(data=alla, stat="identity", shape=21, size=5, colour='red')
p + theme_bw()+ theme(axis.text.x = element_text(angle=90, hjust=1))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+xlab("E.C. Number")+ylab("Abundance (per recA)")
ggsave(file="Supp_Fig_7.eps")


##################################################################################
mdf = psmelt(all_agg_core)
f <- ddply(mdf, .(OTU), summarise, MEAN=mean(Abundance), SE=sd(Abundance)/sqrt(4))
ann <- as.data.frame(tax_table(all_agg_core))
merged <- cbind(f, ann)
summary(merged$MEAN)
colnames(merged)<-c("Contig", "Mean Abundance", "SE","E.C.","CAZy Class","Family", "Hit ID", "Description","Taxonomy Lineage","","")
