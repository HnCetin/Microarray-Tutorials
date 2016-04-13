#For DE analysis
library(Biobase)
library(limma)
library(lumiHumanIDMapping)
library(lumiHumanAll.db)
library(annotate)
library(lumi)
library(GEOquery)

#For batch correction and PEER
library(sva)
library(broom)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)

#Functional programming
library(magrittr)
library(purrr)

#String operations
library(stringr)

#Plotting
library(ggplot2)
library(Cairo)

#Reading and writing tables
library(readr)

library(WGCNA)

#Get GEO data
geo.object <- getGEO('GSE29378', destdir = './save')[[1]]
geo.supplemental <- getGEOSuppFiles('GSE29378', baseDir = './save')

#Format phenoData for GEO SOFT object
pdata.geo <- pData(geo.object) %>% select(title) %>% map(str_split_fixed, " - ", 5) %>% reduce(rbind) %>% data.frame
colnames(pdata.geo) <- c("Diagnosis", "Subject.ID", "Region", "Sex", "Age")
pdata.geo$Sample.ID <- sampleNames(geo.object)
pdata.geo$Age %<>% str_replace("yr", "") %>% as.character %>% as.numeric
pdata.geo$Subject.ID %<>% str_replace("#", "")
pdata.geo$Subject.ID <- str_replace(pdata.geo$Subject.ID, "#", "")
pdata.geo$Region %<>% str_replace("hippocampus ", "")
pdata.geo$Batch <- factor(geo.object$characteristics_ch1.10)

#Fix mesed up column names in raw probe profile
geo.reformat <- read_tsv('./save/GSE29378/GSE29378_non-normalized.txt')
pval.columns <- str_detect(colnames(geo.reformat), "Detection Pval") %>% which
colnames(geo.reformat)[pval.columns] <- paste(colnames(geo.reformat)[pval.columns - 1], "Detection Pval", sep = ".")
colnames(geo.reformat)[pval.columns - 1] %<>% paste("AVG_SIGNAL", sep = ".")
write.table(geo.reformat, "./save/geo_reformat.tsv", row.names = FALSE, sep = '\t')

#Read in fixed probe profile and fix pheno type data
lumi.raw <- lumiR('./save/geo_reformat.tsv', lib.mapping = 'lumiHumanIDMapping') 
pdata.lumi <- map(sampleNames(lumi.raw), str_split_fixed, "\\.", 7) %>% reduce(rbind) %>% data.frame
colnames(pdata.lumi) <- c("Diagnosis", "Subject.ID", "Region", "Sex", "Age", "X1", "X2")
pdata.lumi$sampleID <- lumi.raw$sampleID
pdata.lumi$Region %<>% str_replace("b", "")
pdata.lumi$Age %<>% as.character %>% as.integer
pdata.lumi$Diagnosis %<>% revalue(c(A = "Alzheimers", C = "Control"))
pData(lumi.raw) <- select(pdata.lumi, -(X1:X2))
sampleNames(lumi.raw) <- pdata.lumi$sampleID

#Create unique identifiers for each sample and remove any not found in SOFT object
lumi.names <- paste(pdata.lumi$Subject.ID, pdata.lumi$Region, sep = "_") 
geo.names <- paste(pdata.geo$Subject.ID, pdata.geo$Region, sep = "_")
known.names <- lumi.names %in% geo.names
rep.names <- str_detect(pdata.lumi$sampleID, "rep")
lumi.known <- lumi.raw[,known.names & !rep.names]
lumi.known$Batch <- str_replace(pdata.geo$Batch, "^.*: ", "") %>% factor %>% as.integer
sampleNames(lumi.known) <- paste(lumi.known$Subject.ID, lumi.known$Region, sep = "_")

lumi.ca3 <- lumi.known[,lumi.known$Region == "CA3"]
lumi.log2 <- lumiT(lumi.ca3, "log2")

#Boxplot
expr.log2 <- exprs(lumi.log2) %>% t %>% data.frame(Sample.Name = sampleNames(lumi.log2), Batch = lumi.log2$Batch, Diagnosis = lumi.log2$Diagnosis)
expr.log2.melt <- melt(expr.log2, id = c("Sample.Name", "Batch", "Diagnosis"))
colnames(expr.log2.melt)[4:5] <- c("Symbol", "Intensity")

p <- ggplot(expr.log2.melt, aes(x = Sample.Name, y = Intensity, fill = factor(Batch))) + geom_boxplot() + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0))     
p <- p + ggtitle("Log2 normalized signal intensity") #+ ylab("Intensity") + xlab("Sample") 
p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "boxplot_normalized.jpg", plot = p, family = "Oxygen", width = 15 , height = 8)

#Histogram
p <- ggplot(expr.log2.melt, aes(Intensity, group = Sample.Name, col = Diagnosis)) + geom_density() + theme_bw()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + ggtitle("Histogram of Log2 Expression") + ylab("Density") + xlab("Log2 Expression") 
CairoPDF("histogram_log2", height = 5, width = 9)
print(p)
dev.off()

#Principal components analysis
mds.log2 <- exprs(lumi.log2) %>% t %>% dist %>% cmdscale(eig = TRUE) 
mds.log2.plot <- data.frame(Sample.Name = rownames(mds.log2$points), Diagnosis = lumi.log2$Diagnosis, Batch = factor(lumi.log2$Batch), Component.1 = mds.log2$points[,1], Component.2 = mds.log2$points[,2])

p <- ggplot(mds.log2.plot, aes(x = Component.1, y = Component.2, col = Diagnosis)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle("MDS of Diagnosis")
CairoPDF(file = "mds_diagnosis", height = 6, width = 7)
print(p)
dev.off()

p <- ggplot(mds.log2.plot, aes(x = Component.1, y = Component.2, col = Batch)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle("MDS of Batch")
CairoPDF(file = "mds_batch", height = 6, width = 7)
print(p)
dev.off()

#Log2 normalize intensities, normalize between arrays, make QC object
lumi.norm <- lumiN(lumi.log2, "rsn")
lumi.cutoff <- detectionCall(lumi.norm) #Get the count of probes which passed the detection threshold per sample
lumi.expr <- lumi.norm[which(lumi.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi <- getSYMBOL(rownames(lumi.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.expr.annot <- lumi.expr[!symbols.lumi,] #Drop any probe which is not annotated

#Batch correction
model.combat <- model.matrix(~ Diagnosis + Sex + Age, data = pData(lumi.expr.annot)) %>% data.frame
expr.combat <- ComBat(dat = exprs(lumi.expr.annot), batch = factor(lumi.expr.annot$Batch), mod = model.combat) #Run ComBat
lumi.combat <- lumi.expr.annot 
exprs(lumi.combat) <- expr.combat 

mds.combat <- exprs(lumi.combat) %>% t %>% dist %>% cmdscale(eig = TRUE) 
mds.combat.plot <- data.frame(Sample.Name = rownames(mds.combat$points), Diagnosis = lumi.combat$Diagnosis, Batch = factor(lumi.combat$Batch), Component.1 = mds.combat$points[,1], Component.2 = mds.combat$points[,2])

p <- ggplot(mds.combat.plot, aes(x = Component.1, y = Component.2, col = Diagnosis)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle("MDS of Diagnosis")
CairoPDF(file = "mds_diagnosis_combat", height = 7, width = 7)
print(p)
dev.off()

p <- ggplot(mds.combat.plot, aes(x = Component.1, y = Component.2, col = Batch)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle("MDS of Batch")
CairoPDF(file = "mds_batch_combat", height = 6, width = 7)
print(p)
dev.off()

#Hierarchichal clustering
tree.combat <- exprs(lumi.combat) %>% t %>% dist %>% hclust(method = "average")
CairoPDF("clustering_combat", width = 13, height = 10)
plot(tree.combat, main = "Hierarchical Clustering Sammples")
dev.off()

#Connectivity outliers
normalized.adjacency <- (0.5 + 0.5 * bicor(exprs(lumi.combat))) ^ 2
network.summary <- fundamentalNetworkConcepts(normalized.adjacency)
connectivity <- network.summary$Connectivity
connectivity.zscore <- (connectivity - mean(connectivity)) / sqrt(var(connectivity))
connectivity.plot <- data.frame(Sample.Name = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))

p <- ggplot(connectivity.plot, aes(x = Sample.Num, y = Z.score, label = Sample.Name)) + geom_text(size = 4, colour = "red")
p <- p + geom_hline(aes(yintercept = -2)) 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle("Sample Connectivity")
CairoPDF("connectivity_combat", width = 10, height = 10)
print(p)
dev.off()

#Remove outliers
lumi.rmout <- lumi.combat[,abs(connectivity.zscore) < 2]

#New boxplot
expr.rmout <- exprs(lumi.rmout) %>% t %>% data.frame(Sample.Name = sampleNames(lumi.rmout), Batch = lumi.rmout$Batch, Diagnosis = lumi.rmout$Diagnosis)
expr.rmout.melt <- melt(expr.rmout, id = c("Sample.Name", "Batch", "Diagnosis"))
colnames(expr.rmout.melt)[4:5] <- c("Symbol", "Intensity")

p <- ggplot(expr.rmout.melt, aes(x = Sample.Name, y = Intensity, fill = factor(Batch))) + geom_boxplot() + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0))     
p <- p + ggtitle("Log2 normalized signal intensity") + ylab("Intensity") + xlab("Sample") 
p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "boxplot_rmout.jpg", plot = p, family = "Oxygen", width = 15 , height = 8)

#New histogram
p <- ggplot(expr.rmout.melt, aes(Intensity, group = Sample.Name, col = Diagnosis)) + geom_density() + theme_bw()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + ggtitle("Histogram of Log2 Expression") + ylab("Density") + xlab("Log2 Expression") 
CairoPDF("histogram_rmout", height = 5, width = 9)
print(p)
dev.off()

#New MDS
mds.rmout <- exprs(lumi.rmout) %>% t %>% dist %>% cmdscale(eig = TRUE) 
mds.rmout.plot <- data.frame(Sample.Name = rownames(mds.rmout$points), Diagnosis = lumi.rmout$Diagnosis, Batch = factor(lumi.rmout$Batch), Component.1 = mds.rmout$points[,1], Component.2 = mds.rmout$points[,2])

p <- ggplot(mds.rmout.plot, aes(x = Component.1, y = Component.2, col = Diagnosis)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle("MDS of Diagnosis")
CairoPDF(file = "mds_diagnosis_rmout", height = 6, width = 7)
print(p)
dev.off()

#Covariate analysis
lm.age <- lm(Age ~ Diagnosis, data = pData(lumi.rmout)) %>% anova %>% tidy
lm.sex <- lm(as.numeric(Sex) ~ Diagnosis, data = pData(lumi.rmout)) %>% anova %>% tidy
lm.batch <- lm(Batch ~ Diagnosis, data = pData(lumi.rmout)) %>% anova %>% tidy

#Covariate plots
p <- ggplot(pData(lumi.rmout), aes(x = Diagnosis, y = Age)) + geom_boxplot()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank()) + ggtitle(paste("Age (p < ", round(lm.age$p.value[1], 3), "*)", sep = ""))
CairoPDF("age_boxplot", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(lumi.rmout), aes(x = Diagnosis, fill = Sex)) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.title = element_blank()) 
p <- p + ggtitle(paste("Sex (p < ", round(lm.sex$p.value[1], 3), ")", sep = "")) 
CairoPDF("sex_barplot", height = 6, width = 6)
plot(p)
dev.off()

p <- ggplot(pData(lumi.rmout), aes(x = Diagnosis, fill = factor(Batch))) + geom_bar()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.title = element_blank(), legend.title = element_blank()) 
p <- p + ggtitle(paste("Batch (p < ", round(lm.batch$p.value[1], 3), ")", sep = "")) 
CairoPDF("batch_barplot", height = 6, width = 6)
plot(p)
dev.off()

#Collapse Rows
expr.symbols <- featureNames(lumi.rmout) %>% getSYMBOL('lumiHumanAll.db') %>% factor
expr.collapse <- collapseRows(exprs(lumi.rmout), rowGroup = expr.symbols, rowID = rownames(exprs(lumi.rmout)))
lumi.collapse <- ExpressionSet(assayData = expr.collapse$datETcollapsed, phenoData = phenoData(lumi.rmout))

#Differential Expression
model.limma <- model.matrix( ~ 0 + Diagnosis + Sex, data = pData(lumi.collapse))
fit.object <- lmFit(lumi.collapse, model.limma) 
contrasts.object <- makeContrasts(Alzheimers_vs_Control = DiagnosisAlzheimers - DiagnosisControl, levels = model.limma)
fit.contrasts <- contrasts.fit(fit.object, contrasts.object)
fit.ebayes <- eBayes(fit.contrasts)

top.table.contrasts <- topTable(fit.ebayes, number = Inf)
top.table <- topTable(eBayes(fit.object))

decide.table <- decideTests(fit.ebayes, adjust.method = "fdr", p = 0.05)
decide.df <- data.frame(Gene = rownames(decide.table), Direction = as.vector(decide.table))
