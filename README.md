# Genome-scans
Using genome scan programs and visualizing results

##################### Bayescan ########################3

###################################### REMOVE BALANCING SELECTION SNPs ###############################################

# Open the bayescan output file with the "_fst.txt" extension
bayescan=read.table("bayescan-24603snps-379ind_fst.txt")

# Download the list of SNPs
SNPb=read.table("24603snps-list.txt",header=FALSE)

# Merge the name of the outliers with the results from bayescan
bayescan=cbind(SNPb, bayescan)
colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST")
write.table(bayescan, "24603snps-bayescan-results.txt", quote=FALSE, sep="\t", row.names=FALSE)

# POST_PROB = 1 & Q_VALUE = 0 == 0.0001 (see comment below on LOG10_PO)
attach(bayescan)
class(bayescan$Q_VALUE)
bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE)
bayescan[bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001

bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4))
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4))
bayescan$ALPHA <- (round(bayescan$ALPHA, 4))
bayescan$FST <- (round(bayescan$FST, 6))

# ADD a column for selection grouping based on a Q-VALUE < 0.05
bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.05,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.05,"neutral","balancing"))
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION)

# Save the results of the SNPs potentially under selection (qvalue < 0.05)
positive <- bayescan[bayescan$SELECTION=="diversifying",]
neutral <- bayescan[bayescan$SELECTION=="neutral",]
balancing <- bayescan[bayescan$SELECTION=="balancing",]
xtabs(data=bayescan, ~SELECTION)

### Write the results of the SNPs potentially under selection (qvalue < 0.05)
write.table(neutral, "neutral.txt", row.names=F, quote=F)
write.table(balancing, "balancing.txt", row.names=F, quote=F)
write.table(positive, "positive.txt", row.names=F, quote=F)

#Transformation Log de la valeur Q pour le graph
range(bayescan$Q_VALUE)
bayescan$LOG10_Q <- -log10(bayescan$Q_VALUE)

##TODO change name to figure
#################################### SCATTER PLOT OF BAYESCAN RESULTS ###########################################
library(ggplot2)

x_title="Log(q-value)"
y_title="Fst"

graph_1<-ggplot(bayescan,aes(x=LOG10_Q,y=FST))
graph_1+geom_point(aes(fill=factor(bayescan$SELECTION)), pch=21, size=3)+
  scale_fill_manual(name="Selection",values=c("black","red","white"))+
  labs(x=x_title)+
  labs(y=y_title)+
  #geom_vline(xintercept=1.3,color="red")+
 #annotate("text", x= 0.7, y=0.02,label="10",colour="yellow")+
# annotate("text", x= 1.2, y=0.03,label="3",colour="orange")+
 #annotate("text", x= 1.7, y=0.05,label="3",colour="darkorange")+
  #annotate("text", x= 3.5, y=0.10,label="3",colour="red")+
  theme(axis.title=element_text(size=12, family="Helvetica",face="bold"), legend.position="none")+
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black",size=12))+
  theme(axis.text.x=element_text(colour="black",size=12))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold"))
  
# Save the file in a pdf format
ggsave("Bayescan-sebastes.pdf", dpi=600, width=5, height=5)
dev.off()

