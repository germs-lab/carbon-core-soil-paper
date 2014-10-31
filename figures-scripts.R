library(phyloseq)
library(plyr)
library(ggplot2)
library(reshape)

setwd("/Users/adina/Google Drive/C_Metabolism_Agg_Paper_2013/reproducibility/")
##################################################################################
#Figure 1
counts <- read.delim(sep=',',file="core-counts.csv",header=FALSE)
colnames(counts) <- c('filecount','uniquecazy','bpcount','round')
counts$filecount <- factor(counts$filecount, levels=c("1","2","3","4"))
counts$bocount <- factor(counts$bpcount, levels=c("1","2","3","4"))
counts$bpcount <- as.numeric(counts$bpcount)
p = ggplot(counts, aes_string(x="filecount", y="uniquecazy", colour="bpcount"))
p + geom_point(stat="identity", size=4) + ylab("Number of core sequences") + xlab("Total combined metagenomes") + theme_bw()+theme(text=element_text(size=10, family="Helvetica"))+theme(axis.text.x=element_text(size=15, angle=90),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15))+scale_colour_gradientn(colours=rainbow(2), "Cumulative bp")+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())
ggsave(file="Fig_1.eps")
##################################################################################

##################################################################################
#Figure 2a
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
p+theme(text=element_text(size=10, family="Helvetica"))+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15))
ggsave(file="Fig_2a.eps")
##################################################################################
#Figure 2b
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
#o <- factor(f2$Cazy_fam, levels=f2$Cazy_fam[order(-f2$Agg_sort)])
o <- factor(f2$Cazy_fam, levels=f2$Cazy_fam[order(-f2$MEAN)][1:10])
mdf$agg_frac<-factor(mdf$agg_frac, levels=c("micro","SM","MM","LM","WS"))
f2$Cazy_fam <- factor(f2$Cazy_fam, levels=c(levels(o)))
f2 <- subset(f2, Cazy_fam != "NA")
p = ggplot(f2, aes_string(x="Cazy_fam", y="MEAN"))
p = p + geom_point(stat = "identity", size=5, aes(color=Cazy_fam2, fill=Cazy_fam2))+ geom_errorbar(limits, width=0) 
p = p  + xlab("") + ylab("Abundance (per recA)")+ theme_bw()  + theme(axis.text.x = element_text(angle=90, hjust=1))
p+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+scale_fill_manual("CAZy Enzyme Class", values=c("red","blue","dark green"))+scale_colour_manual("CAZy Enzyme Class", values=c("red","blue","dark green"))
ggsave(file="Fig_2b.eps")
##################################################################################

##################################################################################
#Figure 3
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
p = p + theme_bw()+ theme(text=element_text(size=20))+ylab("Abundance (per recA)") +geom_errorbar(limits, width=0) + theme(strip.text.x = element_text(angle = 90)) + theme(strip.text.y = element_text(angle = 0)) + theme(legend.position="none")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + xlab("")+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())
p
ggsave(file="Fig_3.eps")
##################################################################################

##################################################################################
#Figure 4a-d

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
p = p + theme_bw()+ ylab("Abundance (per recA)") +geom_errorbar(limits, width=0) + facet_grid(~Cazy_fam2) +  theme(text=element_text(size=18), strip.text.x = element_text(angle = 90)) + theme(strip.text.y = element_text(angle = 0)) + theme(legend.position="none")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + xlab("")+opts(panel.grid.major=theme_blank(), panel.grid.minor=theme_blank())
p
ggsave(file="Supp_Fig_2a.eps")

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
p = p + theme_bw()+ ylab("Abundance (per recA)") +geom_errorbar(limits, width=0) + facet_grid(~Cazy_fam2) +  theme(text=element_text(size=18), strip.text.x = element_text(angle = 90)) + theme(strip.text.y = element_text(angle = 0)) + theme(legend.position="none")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + xlab("")+opts(panel.grid.major=theme_blank(), panel.grid.minor=theme_blank())
p
ggsave(file="Supp_Fig_2b.eps")
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
p = p + theme_bw()+ ylab("Abundance (per recA)") +geom_errorbar(limits, width=0) + facet_grid(~Cazy_fam2) +  theme(text=element_text(size=18), strip.text.x = element_text(angle = 90)) + theme(strip.text.y = element_text(angle = 0)) + theme(legend.position="none")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + xlab("")+opts(panel.grid.major=theme_blank(), panel.grid.minor=theme_blank())
p
ggsave(file="Supp_Fig_2c.eps")
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
p = p + theme_bw()+ ylab("Abundance (per recA)") +geom_errorbar(limits, width=0) + facet_grid(~Cazy_fam2) +  theme(text=element_text(size=18), strip.text.x = element_text(angle = 90)) + theme(strip.text.y = element_text(angle = 0)) + theme(legend.position="none")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + xlab("")+opts(panel.grid.major=theme_blank(), panel.grid.minor=theme_blank())
p
ggsave(file="Supp_Fig_2d.eps")
##################################################################################




##################################################################################
#Figure 5a
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
p = p + geom_bar(position="dodge") + geom_errorbar(limits, width=0, position=position_dodge(.9)) 
p + theme_bw() + theme(text=element_text(size=20))+theme(axis.text.x=element_text(size=12, angle=90, hjust=1))+ylab("Relative Abundance")+xlab("")+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+scale_fill_manual("Dataset",values=c("black","white"),labels=c("core","total"))+scale_colour_manual("Dataset",values=c("black","black"), labels=c("core","total"))#+theme
ggsave(file="Fig_4a.eps")

#Stats
f$dataset <- "core"
g$dataset <- "noncore"
combo <- rbind(f,g)
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam2 == "CB")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam2 == "GT")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam2 == "CE")))
summary(aov.result<-aov(SUMs ~ dataset, data = subset(combo, Cazy_fam2 == "GH")))

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
p = p + geom_bar(position="dodge") + geom_errorbar(limits, width=0, position=position_dodge(.9)) 
p + theme_bw() + theme(text=element_text(size=20))+theme(axis.text.x=element_text(size=12, angle=90, hjust=1))+ylab("Relative Abundance")+xlab("")+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+scale_fill_manual("Dataset",values=c("black","white"),labels=c("core","total"))+scale_colour_manual("Dataset",values=c("black","black"), labels=c("core","total"))#+theme

ggsave(file="Fig_4b.eps")


##################################################################################
#Figure 5

library(reshape)
hist_data <- read.delim(sep='\t', file="./histogram_core.txt",header=FALSE, strip.white=TRUE)
colnames(hist_data)<-c('sample','soil_type', 'mgrastids','value')
mdf<-melt(hist_data)
mdf$sample<-factor(mdf$sample, levels=mdf$sample[order(-mdf$value)])
p = ggplot(mdf, aes_string(x="sample", y="value"))
p = p + geom_bar(stat = "identity")
p + theme_bw()+ ylab("Shared Core Sequences")  + theme(text=element_text(size=14)) + theme(legend.position="none")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), axis.text.y=element_text(angle=90)) + xlab("")+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())
ggsave(file="Fig_5.eps")
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


# mdf_melt <- melt(t2)
p = ggplot(mdf_melt, aes_string(x="variable", y="value", color="Cazy_fam", fill="Cazy_fam"))
p + geom_bar(stat="identity")+theme_bw()+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+theme(axis.text.x=element_text(size=15, angle=90),axis.text.y=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+xlab("Soil metagenome")+ylab("Relative Abundance")+scale_fill_manual("Cazy Enzyme Class",values=c("red","blue","dark green","orange","brown"))+scale_colour_manual("Cazy Enzyme Class",values=c("red","blue","dark green","orange","brown"))
ggsave(file="Fig_6.eps")





##################################################################################
#Figure 6a
mdf_core<-psmelt(all_agg_core)
mdf_noncore<-psmelt(all_agg)
f <- ddply(mdf_core, .(Cazy_fam,Cazy_fam2, t2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
f2 <- subset(f, Cazy_fam == "GT4")
g <- ddply(mdf_noncore, .(Cazy_fam, Cazy_fam2, t2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
g2 <- subset(g, Cazy_fam == "GT4")

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
p = ggplot(foo, aes_string(x="t2", y="mean", colour="dataset", fill="dataset"))
p = p + geom_bar(position="dodge") + geom_errorbar(limits, width=0, position=position_dodge(.9)) 
p + theme_bw() + theme(text=element_text(size=20))+theme(axis.text.x=element_text(size=12, angle=90, hjust=1))+ylab("Relative Abundance")+xlab("")+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+scale_fill_manual("Dataset",values=c("black","white"),labels=c("core","total"))+scale_colour_manual("Dataset",values=c("black","black"), labels=c("core","total")) 
ggsave(file="Supp_Fig_3a.eps")
##################################################################################
#Figure 6b

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
p = ggplot(foo, aes_string(x="t2", y="mean", colour="dataset", fill="dataset"))
p = p + geom_bar(position="dodge") + geom_errorbar(limits, width=0, position=position_dodge(.9)) 
p + theme_bw() + theme(text=element_text(size=20))+theme(axis.text.x=element_text(size=12, angle=90, hjust=1))+ylab("Relative Abundance")+xlab("")+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+scale_fill_manual("Dataset",values=c("black","white"),labels=c("core","total"))+scale_colour_manual("Dataset",values=c("black","black"), labels=c("core","total")) 
ggsave(file="Supp_Fig_3b.eps")
##################################################################################
#Figure 6c

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

p = ggplot(foo, aes_string(x="t2", y="mean", colour="dataset", fill="dataset"))
p = p + geom_bar(position="dodge") + geom_errorbar(limits, width=0, position=position_dodge(.9)) 
p + theme_bw() + theme(text=element_text(size=20))+theme(axis.text.x=element_text(size=12, angle=90, hjust=1))+ylab("Relative Abundance")+xlab("")+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+scale_fill_manual("Dataset",values=c("black","white"),labels=c("core","total"))+scale_colour_manual("Dataset",values=c("black","black"), labels=c("core","total")) 
ggsave(file="Supp_Fig_3c.eps")






##################################################################################
#Figure 9



f <- ddply(mdf_core, .(EC, sample_name, agg_frac), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(EC, agg_frac), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f2 <- data.frame(f2)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f3 <- subset(f2, EC != "n.d.")
f3 <- subset(f3, EC != "2.4.1.-")
f3 <- subset(f3, EC != "3.2.1.-")
f3$EC <- factor(f3$EC,  levels=f3$EC[order(-f3$MEAN)])
a1<-f3[which(f3$EC=="3.2.1.20"),]#yes
#a2<-f3[which(f3$EC=="3.2.1.21"),]#yes
#a3<-f3[which(f3$EC=="3.2.1.31"),] 
#a4<-f3[which(f3$EC=="3.2.1.23"),]
a5<-f3[which(f3$EC=="3.2.1.25"),]
#a6<-f3[which(f3$EC=="3.2.1.55"),]
#a7<-f3[which(f3$EC=="3.2.1.91"),] #yes
a8<-f3[which(f3$EC=="3.2.1.37"),]#yes
#a9<-f3[which(f3$EC=="3.2.1.14"),]
#a10<-f3[which(f3$EC=="3.4.11.1"),]
#a11<-f3[which(f3$EC=="3.1.3.1"),]
alla<-rbind(a1,a5,a8)


p = ggplot(f3, aes_string(x="EC", y="MEAN"))
p = p + geom_point(stat = "identity")+ geom_errorbar(limits, width=0)+geom_point(data=alla, stat="identity", shape=21, size=5, colour='red')
p + theme_bw()+ theme(axis.text.x = element_text(angle=90, hjust=1))+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+xlab("E.C. Number")+ylab("Abundance (per recA)")
ggsave(file="Supp_Fig_6.eps")




##################################################################################
#Supp Figure X
CAZY_choice = "PL"
counts <- read.delim(sep='\t',file="cazy-annotations.info",header=FALSE)
colnames(counts)<-c("Cazy_fam2", "Cazy_fam", "hit_id", "t1","t2","t3")
f<-ddply(counts, .(Cazy_fam, t2), summarise, SUM=length(Cazy_fam))
f <- subset(f, Cazy_fam == CAZY_choice)
f2<-ddply(f, .(Cazy_fam, t2, SUM), summarise, NORM=SUM/sum(f$SUM))

mdf <- psmelt(all_agg_core)
f <- ddply(mdf, .(Cazy_fam2,t2, sample_name), summarise, SUM=sum(Abundance))
#f2 <- subset(f, Cazy_fam == "CE10")
#g <- ddply(mdf_noncore, .(Cazy_fam, Cazy_fam2, t2, sample_name, agg_frac), summarise, SUM=sum(Abundance))
#g2 <- subset(g, Cazy_fam == "CE10")
f <- subset(f, Cazy_fam2 == CAZY_choice)
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
	
fam_wgt<-aggregate(f$SUMs,by=list(f$t2),FUN=mean)
fam_wgt2<-aggregate(f$SUMs,by=list(f$t2),FUN=stderr)
names(fam_wgt)[2]<-"mean"
names(fam_wgt2)[2]<-"std_er"
fam_stat<-join(fam_wgt,fam_wgt2)
fam_stat$dataset <- "core"

f2$dataset<- "CAZy"
f2 <- subset(f2, f2 == CAZY_choice)
colnames(f2)<-c("Cazy_fam", "Group.1", "SUM","mean", "dataset")
f2$Cazy_fam <- NULL
f2$SUM <- NULL
f2$std_er <- 0

merge1 <- rbind(fam_stat, f2)
colnames(merge1) <- c("t2", "mean" ,"std_er", "dataset")
merge1 <- subset(merge1, t2 != "")
merge1 <- subset(merge1, t2 != "dsDNA")
merge1 <- subset(merge1, t2 != "environmental")
merge1 <- subset(merge1, t2 != "candidate")
#merge1$Cazy_fam2 <- factor(merge1$Cazy_fam2, levels=merge1$Cazy_fam2[with(merge1, order(-NORM))])
merge1$dataset <- factor(merge1$dataset,levels=c("core","CAZy"))
temp <- subset(merge1, dataset == "core")
merge1$t2 <- factor(merge1$t2, levels=merge1$t2[with(temp, order(-mean))])
merge1 <- subset(merge1, t2 != "NA")
#merge1$t2=factor(merge1$t2, levels=c("Proteobacteria","Bacteroidetes","Firmicutes","Euryarchaeota","Actinobacteria","Chloroflexi","Viridiplantae","Verrucomicrobia","Cyanobacteria","Spirochaetes","Crenarchaeota","Acidobacteria","Planctomycetes","Korarchaeota","Fungi","Chlorobi","Metazoa","Aquificae","Mycetozoa","Chlamydiae","Thermotogae","Deinococcus-Thermus","Alveolata"))

#merge1 <- subset(merge1, Cazy_fam2 != "NA")
#temp <- subset(merge1, dataset=="core")
#temp <- subset(temp, Cazy_fam2 == "GT")
#merge1$t2 <- factor(merge1$t2, levels=temp$t2[with(temp, order(-mean))][1:10])
#merge1$Cazy_fam2 <- factor(merge1$Cazy_fam2, levels=c("GT","GH","CE","PL"))
#merge1 <- subset(merge1, t2 != "<NA>")
limits<-aes(ymin=mean-std_er, ymax=mean+std_er)

p = ggplot(merge1, aes_string(x="t2", y="mean", colour="dataset", fill="dataset"))
p = p + geom_bar(position="dodge") + geom_errorbar(limits, width=0, position=position_dodge(.9)) 
p + theme_bw() + theme(text=element_text(size=20))+theme(axis.text.x=element_text(size=12, angle=90, hjust=1))+ylab("Relative Abundance")+xlab("")+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+scale_fill_manual("Dataset",values=c("black","white"),labels=c("core","total"))+scale_colour_manual("Dataset",values=c("black","black"), labels=c("core","total")) 
ggsave(file="Supp_fig_7_PL.eps")

#Counting core presence in other aggregates
min_count = 5
in_all_samples <- subset(abundance_data, abundance_data$PF_WS_H15 >=min_count &                           abundance_data$PF_WS_H07 >=min_count & abundance_data$PF_WS_H09 >=min_count  & abundance_data$PF_WS_H04 >=min_count)
abundance_data <- in_all_samples
in_all_samples_agg <- subset(abundance_data, abundance_data$PF_LM_H08 >=min_count |                           abundance_data$PF_LM_H14 >=min_count | abundance_data$PF_LM_H16 >=min_count |                           abundance_data$PF_LM_H03 >=min_count | abundance_data$PF_MI_H01 >=min_count |                           abundance_data$PF_MI_H06 >=min_count | abundance_data$PF_MI_H12 >=min_count |                           abundance_data$PF_MI_H13 >=min_count | abundance_data$PF_MM_H17 >=min_count |                           abundance_data$PF_MM_H19 >=min_count | abundance_data$PF_MM_H20 >=min_count |                           abundance_data$PF_SM_H02 >=min_count | abundance_data$PF_SM_H10 >=min_count |                           abundance_data$PF_SM_H11 >=min_count)
abundance_data <- as.matrix(in_all_samples_agg)
abundance <- otu_table(abundance_data, taxa_are_rows=TRUE)
all_agg_pf_agg <- phyloseq(metadata, annotation, abundance)
all_agg_pf_agg

in_all_samples <- subset(abundance_data, abundance_data$PF_WS_H15 >=min_count &                           abundance_data$PF_WS_H07 >=min_count & abundance_data$PF_WS_H09 >=min_count  & abundance_data$PF_WS_H04 >=min_count)
abundance_data <- in_all_samples
in_all_samples_up <- subset(abundance_data, abundance_data$UP_MI_H41 >=min_count |                           abundance_data$UP_MI_H43 >=min_count | abundance_data$UP_MI_H57 >=min_count |                           abundance_data$UP_MI_H59 >=min_count)
abundance_data <- as.matrix(in_all_samples_up)
abundance <- otu_table(abundance_data, taxa_are_rows=TRUE)
all_agg_up_agg <- phyloseq(metadata, annotation, abundance)
all_agg_up_agg

in_all_samples <- subset(abundance_data, abundance_data$PF_WS_H15 >=min_count &                           abundance_data$PF_WS_H07 >=min_count & abundance_data$PF_WS_H09 >=min_count  & abundance_data$PF_WS_H04 >=min_count)
abundance_data <- in_all_samples
in_all_samples_corn <- subset(abundance_data, abundance_data$CC_MI_H30  >=min_count |                           abundance_data$CC_MI_H34 >=min_count)
abundance_data <- as.matrix(in_all_samples_corn)
abundance <- otu_table(abundance_data, taxa_are_rows=TRUE)
all_agg_corn_agg <- phyloseq(metadata, annotation, abundance)
all_agg_corn_agg


####Venn diagram work###
core_taxa = data.frame(tax_table(all_agg_core))
cum_taxa = data.frame(tax_table(all_agg))
core_taxa_gh13 = subset(core_taxa, Cazy_fam == "GT2")
core_taxa_gh13 = subset(core_taxa_gh13, t2 == "Actinobacteria")
core_taxa_gh13 = subset(core_taxa_gh13, t2 == "Proteobacteria")
core_taxa_gh13 = subset(core_taxa_gh13, t2 == "Acidobacteria")
core_taxa_gh13 = subset(core_taxa_gh13, t2 == "Acidobacteria")
cum_taxa_gh13 = subset(cum_taxa, Cazy_fam == "GT2")
cum_taxa_gh13 = subset(cum_taxa_gh13, t2 == "Actinobacteria")
cum_taxa_gh13 = subset(cum_taxa_gh13, t2 == "Proteobacteria")
cum_taxa_gh13 = subset(cum_taxa_gh13, t2 == "Acidobacteria")
cum_taxa_gh13 = subset(cum_taxa_gh13, t2 == "Acidobacteria")
cat(as.character(unique(cum_taxa_gh13$hit_id)), sep="\n")
cat(as.character(unique(core_taxa_gh13$hit_id)), sep="\n")
cat(as.character(unique(cum_taxa_gh13$t3)), sep="\n")
cat(as.character(unique(core_taxa_gh13$t3)), sep="\n")
write.table(unique(core_taxa_gh13), file="proteo-GT2-core.tsv", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(unique(cum_taxa_gh13), file="proteo-GT2-cum.tsv", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)