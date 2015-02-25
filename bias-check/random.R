library(plyr)
library(ggplot2)
library(reshape)
setwd("/Users/adina/software/carbon-core-soil-paper/bias-check")

db_counts <- read.delim(sep=' ',file="cazy-db-1000.txt")
colnames(db_counts)<-c("rep",colnames(db_counts)[3:9])
db_counts$"NA" <- NULL
db_counts$rep <- NULL
db_counts_stats <- db_counts
db_counts_stats$dataset <- "CAZy database"
db_counts2 <- read.delim(sep=' ', file="cazy-contigs-1000.txt")
colnames(db_counts2)<-c("rep",colnames(db_counts2)[3:9])
db_counts2$"NA" <- NULL
db_counts2$rep <- NULL
db_counts2_stats <- db_counts2
db_counts_stats$dataset <- "CAZy database"
db_counts2_stats$dataset <- "Metagenome"
merged_stats <- rbind.fill(db_counts_stats, db_counts2_stats)
#x=subset(merged, dataset=="CAZy database")
#y=subset(merged, dataset=="Metagenome")
mdf <- melt(merged_stats)
mdf <- subset(mdf, mdf$variable != "LO")
mdf$variable <- factor(mdf$variable, levels=mdf$variable[order(-mdf$value)])
p<-ggplot(mdf, aes(variable,value))
p+geom_boxplot(aes(fill=dataset,colour=dataset))+theme_bw()+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+xlab('CAZy Class')+ylab('Counts')
ggsave(file="Supp_Fig_4A.eps")


summary(aov.results<-aov(CB~dataset,data=merged_stats))
summary(aov.results<-aov(GT~dataset,data=merged_stats))
summary(aov.results<-aov(CE~dataset,data=merged_stats))
summary(aov.results<-aov(GH~dataset,data=merged_stats))
summary(aov.results<-aov(LO~dataset,data=merged_stats))
TukeyHSD(aov.results)
summary(aov.results<-aov(PL~dataset,data=merged_stats))
TukeyHSD(aov.results)



db_counts <- read.delim(sep=' ',file="taxa-cazy-db-1000.txt")
colnames(db_counts)<-c("rep",colnames(db_counts)[3:57])
db_counts$"NA" <- NULL
db_counts$rep <- NULL
db_counts_stats <- db_counts
db_counts2 <- read.delim(sep=' ', file="taxa-contigs-db-1000.txt")
colnames(db_counts2)<-c("rep",colnames(db_counts2)[3:41])
db_counts2$"NA" <- NULL
db_counts2$rep <- NULL
db_counts2_stats <- db_counts2
db_counts_stats$dataset <- "CAZy database"
db_counts2_stats$dataset <- "Metagenome"
merged_stats <- rbind.fill(db_counts_stats, db_counts2_stats)
merged_stats[is.na(merged_stats)] <- 0
mdf <- melt(merged_stats)
mdf <- subset(mdf, mdf$variable == "Proteobacteria" | mdf$variable ==  "Firmicutes" | mdf$variable == "Actinobacteria" | mdf$variable == "Acidobacteria"| mdf$variable == "Fungi" | mdf$variable == "Bacterioidetes"| mdf$variable ==  "Cyanobacteria")
#mdf$variable <- factor(mdf$variable, levels=unique(mdf$variable[order(-mdf$value)])[1:10])
mdf$variable <- factor(mdf$variable, levels=c("Proteobacteria", "Firmicutes","Fungi","Actinobacteria", "Cyanobacteria","Acidobacteria"))
p<-ggplot(mdf, aes(variable,value))
p+geom_boxplot(aes(fill=dataset,colour=dataset))+theme_bw()+opts(panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+xlab('Phyla')+ylab('Counts')+theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave(file="Supp_Fig_4B.eps")



#summary(aov.results<-aov(value~dataset,data=subset(mdf, mdf$variable=="Acidobacteria")))
summary(aov.results<-aov(Acidobacteria~dataset,data=merged_stats))
summary(aov.results<-aov(Proteobacteria~dataset,data=merged_stats))
summary(aov.results<-aov(Fungi~dataset,data=merged_stats))
summary(aov.results<-aov(Cyanobacteria~dataset,data=merged_stats))
summary(aov.results<-aov(Actinobacteria~dataset,data=merged_stats))
TukeyHSD(aov.results)
summary(aov.results<-aov(Firmicutes~dataset,data=merged_stats))


