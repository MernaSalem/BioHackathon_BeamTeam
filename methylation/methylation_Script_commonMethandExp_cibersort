#install packages
install.packages("splitstackshape")


#load packages
library(ggfortify)#for pca
library(readr)
library(minfi)
library(matrixStats)
library(stringr)
library(splitstackshape) #divide col to n needed cols by separator
library(impute) # for imputing NAs
library(data.table) #for setDT
library(dplyr)
library(tibble)
library(limma)

#reading methylation data 
methylation=as.data.frame(read_tsv("../GSE224455_series_matrix.txt",comment="!"))
sum(is.na(methylation)) #6360
rownames(methylation) <- methylation$ID_REF
methylation= as.matrix(methylation[,-1])

#getting top 20k most variable features
sds=rowSds(methylation,na.rm = T)
methylation.sds=as.data.frame(cbind(methylation,sds))# because order (next) needs to be done on dataframe not simplt cbound vectors 
methylation.sds.sorted=methylation.sds[order(methylation.sds$sds, decreasing = T),]
top20kvariablemeth=rownames(methylation.sds.sorted)[1:20000]
top20kvariablemethmat=methylation[top20kvariablemeth,]
meth20kdf=as.data.frame(top20kvariablemethmat)
write.csv(meth20kdf,"posthackatho/meth20kdf.csv",) #copy pasted from GEO into csv file

#reading methylation metadata
meta=readxl::read_xlsx("meta_raw.xlsx", col_names =FALSE)
colnames(meta) <- c("ID", "comb")
meta= cSplit(meta, "comb", ",")
meta= cSplit(meta, "comb_1", " ")
meta$combo <- paste0(meta$comb_1_2,"_",toupper(str_sub(meta$comb_2, 1,1)))
colnames(meta)[2] <- "group"
colnames(meta)[4] <- "patno"
colnames(meta)[3] <- "patient"
meta$patients <- paste0(meta$patient,"_",meta$patno)


#making sure the methylation and metadata are ordered
identical(meta$ID,colnames(meth20kdf))#TRUE

#exploratory analysis of methylation data
sum(is.na(meth20kdf)) #2826
#make it numerical
Data_meth <- as.matrix(meth20kdf)
hist (Data_meth,col = "darkmagenta", main = "DNA_Methylation_Array",xlab = "Data") #not bimodal
Data_meth=Data_meth[rowMeans(is.na(Data_meth), na.rm=TRUE)<0.5,] #remove features with 50% or more missingness
#remove rows with zeros
Data_meth=Data_meth[rowSums(Data_meth,na.rm=TRUE)!=0,] #no changel no all zero features

#imputing NAs in methylation data
num.mat.imputed=impute.knn(as.matrix(Data_meth),k=10)$data
num.mat.df= as.data.frame(num.mat.imputed)
#sum(is.na(num.mat.df)) #0

write.csv(num.mat.df,"posthackatho/impmeth.csv")
################################
#Make GenomicRatioSet from Matrix

GRset <- makeGenomicRatioSetFromMatrix(
  mat = num.mat.imputed,
  array = "IlluminaHumanMethylationEPIC", 
  annotation = "ilm10b4.hg19"
)

Data_meth <- as.matrix(num.mat.df)
dim(Data_meth)
#grset <- makeGenomicRatioSetFromMatrix(Data_meth, array= "IlluminaHumanMethylationEPIC")anno.ilm10b4.hg19")
dim(grset)
IlluminaHumanMethylationEPICanno.ilm10b4.hg19
pal <- brewer.pal(8,"Dark2")

par(mfrow=c(1,3)) #Checking clustering by indv
plotMDS(getM(GRset), top=1000, gene.selection="common",col=pal[factor(meta$patno)],labels="*" ,dim=c(1,3))
plotMDS(getM(GRset), top=1000, gene.selection="common",col=pal[factor(meta$patno)],labels="*" ,dim=c(2,3))
plotMDS(getM(GRset), top=1000, gene.selection="common",col=pal[factor(meta$patno)],labels="*" ,dim=c(3,4))
#samples from the same individual cluster together, regardless of sample type

par(mfrow=c(1,3)) #Checking clustering by group
plotMDS(getM(GRset), top=1000, gene.selection="common",col=pal[factor(meta$group)],labels="*" ,dim=c(1,3))
plotMDS(getM(GRset), top=1000, gene.selection="common",col=pal[factor(meta$group)],labels="*" ,dim=c(2,3))
plotMDS(getM(GRset), top=1000, gene.selection="common",col=pal[factor(meta$group)],labels="*" ,dim=c(3,4))
#no prcomps show clustering by group at this global methylation pattern (no differential analysis done yet)

# remove probes with SNPs at CpG site
SNPremoved <- dropLociWithSnps(GRset)
dim(SNPremoved) #15215 18
dim(GRset) # 19893 18

## remove sex chromosome probes
keep <- !(featureNames(SNPremoved) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFlt <- SNPremoved [keep,]
dim(mSetSqFlt)

keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,]
dim(mSetSqFlt) #10136 18

# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,3))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",col=pal[factor(meta$group)],labels="*" ,dim=c(1,3))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",col=pal[factor(meta$group)],labels="*" ,dim=c(2,3))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",col=pal[factor(meta$group)],labels="*" ,dim=c(3,4))

par(mfrow=c(1,3))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",col=pal[factor(meta$patno)],labels="*" ,dim=c(1,3))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",col=pal[factor(meta$patno)],labels="*" ,dim=c(2,3))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",col=pal[factor(meta$patno)],labels="*" ,dim=c(3,4))

#still clustering is mostly showing samples from same individual together regardless of status lesio or nawm
# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)
write.csv(mVals, "posthackatho/mVals.csv")
write.csv(bVals, "posthackatho/bVals.csv")

par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=meta$group, main="Beta values",
            legend=FALSE, xlab="Beta values")
densityPlot(mVals, sampGroups=meta$group, main="M-values",
            legend=FALSE, xlab="M values")

####Probe-wise differential methylation analysis

individual <- factor(meta$patients)
type <- factor(meta$group)
design <- model.matrix(~0+group+individual, data=meta)
colnames(design) <-   c(levels(type),levels(individual)[-1])
colnames(design) <-   c(levels(type),levels(individual)[-1])

fit <- lmFit(mVals, design)

contMatrix <- makeContrasts(lesion-NAWM,levels=design)
contMatrix

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

# get the table of results for the specified contrast 

##annotation of probes
annEPIC <-as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
sub=annEPIC[match(rownames(mVals),rownames(annEPIC)),] #10136


ann_sub= d[rownames(num.mat.imputed),] #filtered out features with 50% NA
compann_sub <- ann_sub[complete.cases(cibersortdf), ]
compann_sub <- ann_sub[complete.cases(ann_sub), ]
annsubfinal= ann_sub %>% filter(chr != "NA")
length(rownames(annsubfinal))
length(annsubfinal$UCSC_RefGene_Name)

final=num.mut.imp.df[match(rownames(annsubfinal),rownames(num.mat.df)),]
head(rownames(final))
rownames(final) <- annsubfinal$UCSC_RefGene_Name
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=sub) #10136
topall <- topTable(fit2, coef=1, genelist=sub, number=Inf, adjust="fdr", p.value=0.05) #8343
#double checked that ater 8343 the adjpvalue is > 0.05

#topall <- topTable(fit2, coef=1, genelist=ann_sub, number=Inf, adjust="fdr", p.value=0.1, lfc=log2(1))
#in paper they didn't filter out by lfc, and i can see why, all lfc are low 
write.csv(topall,"posthackatho/topall.csv")

length(unique(topall$UCSC_RefGene_Name)) # 3978 unique genes are differentially methylated #actually 3622 r unique for reason below
#no this isn't correct for two reasons; first many annotated UCSC names are blank (probe not mapped to gene),second: as annotations are written as BRCA1;BRCA1;BRCA1 and in another line BRCA1;BRCA1 (so these r not couted as 1 unique )
#work on this

sum(is.na(newtopall$UCSC_RefGene_Name))
length(newtopall$UCSC_RefGene_Name[newtopall$UCSC_RefGene_Name==""])#1960 are blank ; no gene annotation so rem these if interested in genes (can keep and could be interestingly simply un-annotated yet yet regulate imp genes and r biologically imp)
newtopall= setDT(topall)[, uniqCnt := sapply(strsplit(UCSC_RefGene_Name, ";"), uniqueN)][]
newtopall= setDT(newtopall)[, uniq := sapply(strsplit(UCSC_RefGene_Name, ";"), unique)][]
table(newtopall$uniqCnt)
# 0    1    2    3    4   15   22 
# 1960 6007  354   18    2    1    1
#most probes map to 1 gene; but there are a few 100 instances where they map to multiple genes (so UCSC lists all those where probe matches perfectly)
head(newtopall$uniq)
#work on probes mapped to annotated genes (no blanks)
clean_newtopall= newtopall[newtopall$UCSC_RefGene_Name != "",] #6383
methgenes= unlist(clean_newtopall$uniq) #all gene names to which all probes are mapped
length(methgenes)#6814
length(unique(methgenes))#3622
unique(methgenes) -> uniquemethgenes

#counting how many times a probe was mapped to a gene (looking at all unique genes)
#'\\b' is that we count instances of "MBP" but not "CMBP" for example
data.frame(uniquemethgenes, t(sapply(paste0('\\b', uniquemethgenes, '\\b'), function(x) {
  tmp <- grepl(x, clean_newtopall$UCSC_RefGene_Name)
  c(perc = mean(tmp) * 100, 
    Freq = sum(tmp))
})), row.names = NULL) -> result
head(result)
#sorting results with gene counted most times
resultord=result[order(result$Freq, decreasing = TRUE),]

#double checking if grepl on UCSC_refgene_name col gives same results as the uniq col #it does
data.frame(uniquemethgenes, t(sapply(paste0('\\b', uniquemethgenes, '\\b'), function(x) {
  tmp <- grepl(x, clean_newtopall$uniq)
  c(perc = mean(tmp) * 100, 
    Freq = sum(tmp))
})), row.names = NULL) -> resultuniq
head(resultuniq)
resultuniqord=result[order(resultuniq$Freq, decreasing = TRUE),]
write_excel_csv(newtopall, "posthackatho/newtopall.csv")
#write.csv (didn't allow writing of lists in the uniq col)


dim(mVals)#8343 10136    18
length(rownames(topall)) #8343
#get matrix subetted to sig DMPs only from mVal and bVals
mVals_sig= mVals[topall$Name,] #Mvalues is better for analysis (use as matrix for MOFA, LDA etc)
bVals_sig= bVals[topall$Name,] #bvalues is better for visualization (as it has limits from 0-1; check distribution etc) 
write.csv(mVals_sig, file="posthackatho/mVals_sig.csv")
write.csv(bVals_sig, file="posthackatho/bVals_sig.csv")

##***************************************************************************************##

##reading in expression data
DEGS <- read.csv("deseq.deg.csv", row.names=1)# deseq2 toptable output
normexpcounts <- read.csv("normalized_counts.csv", row.names=1) #deseq2 normalized output
#View(normexpcounts[1:5,1:5]) #has geneID as rownas
#View(DEGS[1:5,1:5])##has geneID as rownas

####reading expression annotation file 
anngenes=as.data.frame(read_tsv("Human.GRCh38.p13.annot.tsv"))
#View(anngenes[1:5,1:5])
#colnames(anngenes)
##subset to degs only
degsann= anngenes[match(rownames(DEGS),anngenes$GeneID),]
degsannsymbol= anngenes[match(rownames(DEGS),anngenes$GeneID),]$Symbol
#search in DMPs for probes in genes that are DEGS from expression analysis; matching gene symbol using 'grepl'
#in all 20kmeth data
#the collapse "|" is to search for "geneA|geneB|geneC.." from vector of genes "geneA,B,C,.." which is degsannsymbol in (vctor of DMPs$USCS_REF_NAME)
ann_mVals= annEPIC[rownames(mVals),]
ann_mVals_sig= annEPIC[rownames(mVals_sig),]
mVals_keep <- ifelse(grepl(paste(degsannsymbol, collapse = "|"), ann_mVals$UCSC_RefGene_Name), "Keep","Discard")
#in sig. meth data
mVals_sig_keep <- ifelse(grepl(paste(degsannsymbol, collapse = "|"), newtopall$UCSC_RefGene_Name), "Keep","Discard")

mValskeepdf= as.data.frame(cbind(mVals, mVals_keep))
table(mValskeepdf$mVals_keep) # 
# Discard    Keep 
# 10109      27
probestokeep=mValskeepdf %>% filter(mVals_keep == "Keep") #matrix of the common dmps

#getting anno of these common DEG,DMPs writing probes as gene_probe to visually know its related to which gene
ann_mVals_common= ann_mVals[rownames(probestokeep),]
ann_mVals_common= setDT(ann_mVals_common)[, uniq := sapply(strsplit(UCSC_RefGene_Name, ";"), unique)][]
common_newprobenames= vector()

for (i in 1:nrow(ann_mVals_common)) {
  if (length(unlist(ann_mVals_common$uniq[i])) == 1) {
    common_newprobenames[i] <- paste0(ann_mVals_common$uniq[i][[1]],"_",rownames(ann_mVals_common)[i])
  } else if (length(unlist(ann_mVals_common$uniq[i])) > 1) {
    common_newprobenames[i] <- paste0(paste0(unlist(ann_mVals_common$uniq[i]),collapse = '_'),"_",rownames(ann_mVals_common)[i])
  } 
}

ann_mVals_common$geneprobeid <- common_newprobenames
#this above is a clear identifier for probes when we discuss meth and exp of common genes (common from differential analysis)

mVals_sig_common <- mVals[ann_mVals_common$Name,]
dim(mVals_sig_common)#27 *18
identical(ann_mVals_common$Name, rownames(renamed_mVals_sig_common)) #true

renamed_mVals_sig_common = mVals_sig_common
rownames(renamed_mVals_sig_common) <- ann_mVals_common$geneprobeid
#get common gene names
genesfromprobeuniq=unique(unlist(ann_mVals_common$uniq)) #as GeneIDs
sig_anngenes=anngenes[match(rownames(DEGS),anngenes$GeneID ),] #as Symbols
commongenes= intersect(sig_anngenes$Symbol,genesfromprobeuniq)#%in% rownames(DEGS)] # 13 genes
length(unique(anngenes$Symbol))
anngenes$Symbol[duplicated(anngenes$Symbol)]#  "TRNAV-CAC" "TRNAV-CAC"
#since that gene is not of interest (not part of the common genes; so instead of running aggregate just for this gene; remove from df)
dupl= c("107985753","107985615","107985614")
normexpwithoutdupl=as.data.frame(normexpcounts)[!(rownames(normexpcounts) %in% dupl),]
rownamesforexpwithdupl=anngenes[match(rownames(normexpwithoutdupl),anngenes$GeneID),]$Symbol
rownames(normexpwithoutdupl) <- rownamesforexpwithdupl 
subexp= normexpwithoutdupl[commongenes,]

#now ready for greplscript for fetching gene-probe instances and drawing boxplot figures
#need t of both matrices t(renamed_mVals_sig_common), t(subexp)
#metadata for both meta, exp_meta 

exp_meta <- read.csv("exp_meta.csv")
exp_meta$group <- factor(exp_meta$group, labels= c("lesion", "NAWM"))
meta$group <- as.factor(meta$group)


identical(exp_meta$ID,colnames(subexp)) #true
identical(meta$ID,colnames(renamed_mVals_sig_common)) #true

#for expression it's easy as features are gene names; not gene_probe which needs grepl by gene
boxplot(FGD6~group, data=  as.data.frame(cbind(t(subexp),group=exp_meta$group)), 
        main = "FGD6",
        xlab = "group",
        ylab = "exp", 
        names= levels(exp_meta$group))

#trying grepl script fetching meth and geneexp results from gene name and drawing boxplot
#for FGD6
rownames(renamed_mVals_sig_common)[grepl("FGD6", rownames(renamed_mVals_sig_common))]
#gives 2 probes

#plotting as script
par(mfrow=c(2,2))
boxplot(FGD6_cg23218272~group, data=  as.data.frame(cbind(t(renamed_mVals_sig_common),group=meta$group)), 
        main = "FGD6_cg23218272",
        xlab = "group",
        ylab = "Bvalues", 
        names= levels(exp_meta$group))
boxplot(FGD6_cg04874286~group, data=  as.data.frame(cbind(t(renamed_mVals_sig_common),group=meta$group)), 
        main = "FGD6_cg04874286",
        xlab = "group",
        ylab = "Bvalues", 
        names= levels(exp_meta$group))
boxplot(FGD6~group, data=  as.data.frame(cbind(t(subexp),group=exp_meta$group)), 
        main = "FGD6",
        xlab = "group",
        ylab = "exp", 
        names= levels(exp_meta$group))


####making it into an easier function # where we input methmatrix, methmetadata, gene name we want to grep+plot, expmatrix+expmetadata 
grepnplot=function(methmat,methmeta,gene,expmat,expmeta){
  greploutput=grepl(gene, rownames(methmat))
  probenames=rownames(methmat)[grepl(gene, rownames(methmat))]
  match(greploutput,rownames(methmat))
  indices=which(greploutput)
  rowexp=t(expmat[gene,])
  par(mfrow=c((length(indices)/2 +1),2))
  for (i in indices){
    skip_to_next=FALSE
    tryCatch({
      row1=t(methmat)[,i]
      boxplot(row1 ~ methmeta$group, main = colnames(t(methmat))[i],xlab = "group",ylab = "Bvalues", names= levels(methmeta$group))
      
    })}
  boxplot(rowexp ~ expmeta$group, main = gene,xlab = "group",ylab = "normexp", names= levels(expmeta$group))
}

#works nicely elhamdolelah :)

grepnplot(renamed_mVals_sig_common,meta,"FGD6",subexp,exp_meta)
grepnplot(renamed_mVals_sig_common,meta,"SLAIN1",subexp,exp_meta)# for plots which would have 10 mini plots; we'll get error margins too large; so better have that as seperate plots saved into one pdf

#another similar function but for saving all plots in 1 pdf instead of stacking mini-plots (for genes that have so many probes in meth)

grepnplot_pdf=function(methmat,methmeta,gene,expmat,expmeta){
  greploutput=grepl(gene, rownames(methmat))
  probenames=rownames(methmat)[grepl(gene, rownames(methmat))]
  match(greploutput,rownames(methmat))
  indices=which(greploutput)
  rowexp=t(expmat[gene,])
  #par(mfrow=c((length(indices)/2 +1),2))
  pdf("allfigures.pdf")
  for (i in indices){
    skip_to_next=FALSE
    tryCatch({
      row1=t(methmat)[,i]
      boxplot(row1 ~ methmeta$group, main = colnames(t(methmat))[i],xlab = "group",ylab = "Bvalues", names= levels(methmeta$group))
      
    })}
  boxplot(rowexp ~ expmeta$group, main = gene,xlab = "group",ylab = "normexp", names= levels(expmeta$group))
dev.off()
}

#test # works nicely
grepnplot_pdf(renamed_mVals_sig_common,meta,"SLAIN1",subexp,exp_meta) #works nicely; P.S need to delete/move the output pdf file from directory to run fcuntion on new gene

####################################################################################################

#CIBERSORT
#cibersort needs TPM normalized bulk RNA expression matrix; as its signature immune cells is in TPM normalized exp
#download TPM normalized counts of our experiment data from their GEO link
TPM_exp= read_tsv("GSE224377_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
anngenes=as.data.frame(read_tsv("Human.GRCh38.p13.annot.tsv"))
#TPM_exp$Symbol= anngenes$Symbol
#one issue is gene TRNAV-CAC is present in 3 rows (we can spare this 1 gene repreented via 3 probes on the array)
#first cbind the expTPMmatrix and the symbol col 
identical(anngenes$GeneID, TPM_exp$GeneID)#TRUE
TPM_exp_symbol=cbind(TPM_exp,Symbol=anngenes$Symbol)
anngenesa = anngenes %>% filter(Symbol != "TRNAV-CAC")
library(dplyr)
TPM_exp_sub = as.data.frame(TPM_exp_symbol) %>% filter(Symbol != "TRNAV-CAC") 
ncol(TPM_exp_sub) #we dont need the Symbol col anymore
TPM_exp_sub2= TPM_exp_sub %>% select(-Symbol)
rownames(TPM_exp_sub2) <- anngenesa$Symbol
forcibersortexp= TPM_exp_sub2[read_sig$`Gene symbol`,]
nrow(forcibersortexp)#547
#non_na_rows <- df[rowSums(is.na(df)) < ncol(df), ]
forcibersortexp_compcases <- forcibersortexp[rowSums(is.na(forcibersortexp)) < ncol(forcibersortexp), ]
nrow(forcibersortexp_compcases)#520 
notfound=setdiff(read_sig$`Gene symbol`, rownames(forcibersortexp_compcases))
foundornot=grepl(paste(notfound, collapse = "|"), anngenes$Synonyms)
sum(foundornot)#79 (but here grep gets gene abcd not just abc 
anngenes$Synonyms[foundornot] # these are the synonym elements which were matched using grepl on whole synonym col

#double check if theyre present as synonyms not symbols
anngenesa1= cSplit(anngenesa, "Synonyms", "|")
syn=colnames(anngenesa1)[grepl("Synonyms", colnames(anngenesa1))] #specifying the cols to search in 
indicessyn=grep("Synonyms", colnames(anngenesa1))#indices of the cols to search in 

#for loop to check on all synonym_x cols (which we extended from one synonym col which had genes as geneA;geneB;geneC as we couldnt grep exact word as "geneA;geneB;geneC" is not the exacy word as geneB for example while it does include geneB; so we had to divide each single word (single gene) into col then search into these cols )
anngenesa2= as.data.frame(anngenesa1) #need to make anngenesa1 into dataframe (even if is.data.frame gives true!), so we can access cols by i later in forloop

toad1=data.frame() 
toadrows=data.frame()
for (i in indicessyn) {
  toadrows= anngenesa2[(unlist(anngenesa2[,i]) %in% notfound),]
  toad1=rbind(toad1, toadrows)
}
#saves output as variable toad1
#toad1 is to-add rows to the easiy subsetted genes as their symbol name was directly matching cibersort gene name

#yes
vectorr=vector()
vectorrr= vector()
for (i in 1:length(notfound)) {
  vect= apply(toad1, 1, function(r) any(r == notfound[i]))
  if(sum(vect, na.rm=TRUE) > 0) {
  
    vectorr[i] <- which(vect==TRUE)
  } else if (sum(vect, na.rm=TRUE) == 0) {
    vectorrr[i] <- notfound[i]
  }
}
df= data.frame(notfound=notfound, indexrowtoad=vectorr)
complete_df <- na.omit(df)
toadcopy=toad1 
vectorrr
reallynotfound=vectorrr[!is.na(vectorrr)]
#intersect(notfound,complete_df$notfound)
#setdiff(notfound,complete_df$notfound)

complete_df[order(complete_df$indexrowtoad, decreasing = TRUE),]
toadcopytorename =toadcopy[complete_df$indexrowtoad,]
fromtoad1=TPM_exp_sub2[toadcopytorename$Symbol,]# we need to subset them using the symbol equivalent of the notfound names
rownames(fromtoad1) <- complete_df$notfound

forcibersortfinal= rbind(forcibersortexp_compcases, fromtoad1) #540/549 missing 9 features from cibersort signature 
forcibersortfinal1 = forcibersortfinal[,-1]
forcibersortfinal2 <- tibble::rownames_to_column(forcibersortfinal1, "GeneID")
write.table(forcibersortfinal2, "posthackatho/forcibersortfinal2.txt", sep= "\t", quote=F, row.names = F)

#############################################3
#cibersort command
#below is already previously run in script
#library(CIBERSORT)
#sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
#read_sig <- read_tsv(sig_matrix) #547 gene features, 23 cell types
input <- "posthackatho/forcibersortfinal2.txt"

res_ciber_MS_new <- cibersort(sig_matrix, input, perm = 100, QN = FALSE)#
res_ciber_MS_newdf <- as.data.frame(res_ciber_MS_new)
identical(rownames(res_ciber_MS_newdf), samplesheet_exp$ID)#true
res_ciber_MS_newdf$group= samplesheet_exp$group
#rename cols to remove whitespace to refer t these cols later in boxplt
colnames(res_ciber_MS_newdf) <- gsub(" ","",colnames(res_ciber_MS_newdf))

#for boxplot 2 subplots; one for lesion cell proportions one for NAWM
ciber_lesions=res_ciber_MS_newdf %>% filter(group=="lesion")
ciber_normal= res_ciber_MS_newdf %>% filter(group=="NAWM")
#cilcil
par(mar=c(15,2,1,1),mfrow=c(1,2))
boxplot(ciber_lesions[,1:22],main = "immune cells in lesions",ylab = "proportions", las=2)
#, group)
boxplot(ciber_normal[,1:22],main = "immune cells in NAWM",ylab = "proportions", las=2)


#carry out t-test between each cell type in both groups in a function:
ttest_results= data.frame()
ttesting=function(mat){
  for (i in 1:(length(colnames(mat))-4)){ #bcoz 1 col is group, 3 others are cibersort cols pval,corr,rmse
      ttest= t.test(mat[,i] ~ group, data= mat)        
      new_record = data.frame(name=colnames(mat)[i],estimate=paste0(ttest$estimate[[1]],ttest$estimate[[2]]),p_value=ttest$p.value, CI= paste0(ttest$conf.int[[1]],ttest$conf.int[[1]]))
      
      #     new_record = data.frame(name=ttest$data.name,estimate=paste0(ttest$estimate[[1]],ttest$estimate[[2]]),p_value=ttest$p.value, CI= paste0(ttest$conf.int[[1]],ttest$conf.int[[1]]))
      ttest_results =rbind(ttest_results ,new_record)
  }
  return(ttest_results)
}

ttesting(res_ciber_MS_newdf) -> ttest_results
#view the top sig. immune cell subtypes different between lesions and NAWM
head(ttest_results[order(ttest_results$p_value),])
sigcells=head(ttest_results[order(ttest_results$p_value),])$name

indexofsigdiffcells=which(colnames(res_ciber_MS_newdf) %in% sigcells)
par(mfrow=c(6,1))
pdf("posthackatho/boxplotstopcibersortcelltypes.pdf")
for (i in indexofsigdiffcells) {
  boxplot(res_ciber_MS_newdf[,i] ~ group, data = res_ciber_MS_newdf,
          main = colnames(res_ciber_MS_newdf)[i], #change to Celltype 
          xlab = "Group",
          ylab = "cell proportion")
}
dev.off()

#or one by one (cell types of interest)
boxplot(Bcellsmemory ~ group, data = res_ciber_MS_newdf,
        main = "Immune cell proportion in groups", #change to Celltype 
        xlab = "Group",
        ylab = "immune cell proportion")

################################################
#pca for cibersort results to see if immunecell proportions can cluster lesions vs. NAWM together
cell.pca=prcomp(res_ciber_MS_newdf[,1:22]) # - the last 4 cols as theyre pval,corr,rmse and group
df_out <- as.data.frame(cell.pca$x)
options(repr.plot.width=8,repr.plot.height=6)

ggplot(df_out,aes(x=PC1,y=PC2,shape=res_ciber_MS_newdf$group, color= as.factor(samplesheet_exp$ptno)))+
  geom_point()+ggtitle("")+labs(color='')+
  geom_point(size=8,alpha=0.5)+ #Size and alpha just for fun
  theme(  plot.title = element_text(hjust = 0.5,size=15,face = "bold"),
          axis.text.x = element_text( size = 15, angle = 45, hjust = .5, vjust = 0.5, face = "plain"),
          axis.text.y = element_text( size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
          axis.title.x = element_text( size = 15, angle = 0, hjust = .5, vjust = 0, face = "bold"),
          axis.title.y = element_text( size = 15, angle = 90, hjust = .5, vjust = .5, face = "bold"),
          #legend.title=element_text(size=20),
          legend.title=element_blank(), # remove legend title name
          legend.text = element_text(size=15,face="plain"),
          strip.text = element_text(size = 15,face="plain") ,
          legend.position="right",
          
          # for transparent background
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
          axis.line = element_line(colour = "black") # adding a black line for x and y axis
  ) +xlab(paste0("PC 1 (", round(cell.pca$sdev[1],1),"%)"))+
  ylab(paste0("PC 2 (", round(cell.pca$sdev[2],1),"%)"))

save.image(file='posthackatho/postSOLEimage.RData')
#load('posthackatho/postSOLEimage.RData')
