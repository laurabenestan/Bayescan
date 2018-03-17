# Genome-scans
Using Bayescan progra and visualizing results in R

### Download libraries
```{r}
library(ggplot2) 
```

### Download data

* **Open the bayescan output file with the "_fst.txt" extension**. 
```{r}
bayescan=read.table("bayescan-13688snps-562ind.g_fst.txt") 
```

* **Download the list of SNPs**.
```{r}
SNPb=read.table("list-13688snps.txt",header=FALSE) 
```

* **Merge the name of the outliers with the results from bayescan**. 
```{r}
bayescan=cbind(SNPb, bayescan) 
colnames(bayescan)=c("SNP","PROB","LOG_PO","Q_VALUE","ALPHA","FST") 
write.table(bayescan, "24603snps-bayescan-results.txt", quote=FALSE, sep="\t", row.names=FALSE) 
```

### Edit the data

* **Change the value of the Q_VALUE column: 0 == 0.0001**.  
```{r}
attach(bayescan)
class(bayescan$Q_VALUE)  
bayescan$Q_VALUE <- as.numeric(bayescan$Q_VALUE) 
bayescan[bayescan$Q_VALUE<=0.0001,"Q_VALUE"]=0.0001 
```

* **Round the values**.  
```{r}
bayescan$LOG_PO <- (round(bayescan$LOG_PO, 4)) 
bayescan$Q_VALUE <- (round(bayescan$Q_VALUE, 4)) 
bayescan$ALPHA <- (round(bayescan$ALPHA, 4)) 
bayescan$FST <- (round(bayescan$FST, 6))
```

* **Add a column for the type of selection grouping based on a Q-VALUE < 0.05 (you can also choose a Q-VALUE < 0.01)**  
```{r}
bayescan$SELECTION <- ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE<=0.05,"diversifying",ifelse(bayescan$ALPHA>=0&bayescan$Q_VALUE>0.05,"neutral","balancing")) 
bayescan$SELECTION<- factor(bayescan$SELECTION)
levels(bayescan$SELECTION) 
```

### Save the results

* **Save the results of the SNPs potentially under positve (divergent) and balancing selection (qvalue < 0.05)**. 
```{r}
positive <- bayescan[bayescan$SELECTION=="diversifying",] 
neutral <- bayescan[bayescan$SELECTION=="neutral",] 
balancing <- bayescan[bayescan$SELECTION=="balancing",]  
```

* **Check the number of SNPs belonging to each category**. 
```{r}
xtabs(data=bayescan, ~SELECTION) 
```

* **Write the results of the SNPs potentially under selection (qvalue < 0.05)**. 
```{r}
write.table(neutral, "neutral.txt", row.names=F, quote=F)  
write.table(balancing, "balancing.txt", row.names=F, quote=F) 
write.table(positive, "positive.txt", row.names=F, quote=F) 
```

### Create a ggplot graph

* **Transformation Log of the Q value in order to create te ggplot graph**. 
```{r}
range(bayescan$Q_VALUE) 
bayescan$LOG10_Q <- -log10(bayescan$Q_VALUE) 
```

* **Create title for the ggplot graph**. 
```{r}
x_title="Log(q-value)" 
y_title="Fst" 
```

* **Make the ggplot graph**. 
```{r}
graph_1<-ggplot(bayescan,aes(x=LOG10_Q,y=FST)) 

graph_1<-ggplot(fst_loci2,aes(x=LOG10_Q,y=Fst_1, label=fst_loci2$SNP)) 
graph_1+geom_point(aes(fill=fst_loci2$SELECTION), pch=21, size=4)+ 
  #geom_text()+ 
  scale_fill_manual(name="Selection",values=c("white","darkgrey","white"))+ 
  labs(x=x_title)+ 
  labs(y=y_title)+ 
  theme(axis.title=element_text(size=12, family="Helvetica",face="bold"), legend.position="none")+ 
  theme(axis.text.x=element_text(colour="black"))+ 
  theme(axis.text.y=element_text(colour="black",size=12))+ 
  theme(axis.text.x=element_text(colour="black",size=12))+ 
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),  
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold")) 
```
   
### Save the file in a pdf format
```{r}
ggsave("bayescan_13688_562ind.pdf", dpi=600, width=5, height=5) 
dev.off()
```
