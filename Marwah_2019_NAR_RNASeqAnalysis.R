#================================================
# Created by - Dr. Ishita Marwah
# Created on - 10/10/2018
# This script reads the raw count files to produce
# differential expression of reads for a given set of 
# sample types (groups) using edgeR package for RNASeq analysis
#================================================
# Clear the workspace
WS = c(ls())
rm(WS, list = WS)
options("scipen"=100, "digits"=2)
# read transcript data
# convert ensembl transcript ids to hgnc_ids
library(org.Hs.eg.db)
library(stringi)
library('biomaRt')
library(edgeR)
library("corrplot")
library("psych")
library("corrplot")
library(RColorBrewer)



readDECountFile<-function(sampleSource,sampleType)
{
		df_dat<-read.csv(paste0("../",sampleSource,"_RawReadCount/transcript_count_matrix_", sampleType,".csv"),header = T, stringsAsFactors = F)
		return(df_dat)
}


getHGNCSymbols<-function(sampleSource,sampleType)
{
	outFileName<-paste0("Annot_", basename(paste0("../",sampleSource,"_RawReadCount/transcript_count_matrix_", sampleType,".csv")))
	
	df_dat<-readDECountFile(sampleSource,sampleType)

	G_list <- getBM(filters= "ensembl_transcript_id_version", attributes= c("ensembl_transcript_id_version","hgnc_symbol"),values= unique(df_dat$transcript_id),mart= ensembl)
  names(G_list)[1]<-"transcript_id"
  
  # merge with passed df_dat
  annotated_df_dat<- merge(df_dat, G_list, by = 'transcript_id', all = TRUE)	
  write.csv(annotated_df_dat,paste0("../",sampleSource,"_RawReadCount/",outFileName))
}

#----------------------------------------------------
# This function is to get HGNC symbols for the data 
# frame (df_dat) containing the complete set of DET
# as is created in the getDET function. 
# That means, it would have removed the transcripts 
# those were all zero across all the samples during 
# differential read count (expression) analysis. 
#----------------------------------------------------
writeDETWithHGNCSymbols <-function(df_dat,fileName)
{
	df_dat$transcript_id <-rownames(df_dat)
	G_list <- getBM(filters= "ensembl_transcript_id_version", attributes= c("ensembl_transcript_id_version","hgnc_symbol"),values= unique(df_dat$transcript_id),mart= ensembl)
  names(G_list)[1]<-"transcript_id"
# merge with passed df_dat
  annotated_df_dat<- merge(df_dat, G_list, by = 'transcript_id', all = TRUE)	
  write.csv(annotated_df_dat,paste0("../",fileName))
return(annotated_df_dat)
}

#-------------------------------------------------------------
# Differential expression will be alphabetically i.e. uses the
# alphabetically first group as baseline for example,
# comparison between Fix and UF,  will give transcripts 
# differentially expressed in UF Vs Fix. To reverse this
# comparison to be Fix Vs UF, set the reverseComparison = TRUE
#--------------------------------------------------------------
	
getDET<-function(sampleSource1,sampleSource2, sampleType1, sampleType2, reverseComparison=FALSE, removeTrizol=FALSE)
{
	# read transcript data
	transS1 <-readDECountFile(sampleSource1, sampleType1)
	transS2 <-readDECountFile(sampleSource2, sampleType2)
	
	names(transS1) <- paste0(names(transS1),paste0(".", sampleSource1,"_",sampleType1))
	names(transS2) <- paste0(names(transS2),paste0(".", sampleSource2,"_",sampleType2))

	trans_merge <- merge(transS1, transS2, 
	by.x = paste0("transcript_id.", sampleSource1,"_",sampleType1), 
	by.y = paste0("transcript_id.", sampleSource2,"_",sampleType2),all = TRUE)
	
	rownames(trans_merge)<-trans_merge[,1]
	trans_merge<-trans_merge[,-1]
	
	trans_merge[is.na(trans_merge)] <- 0
	
	sampleGp<-unlist(lapply(strsplit(paste(names(trans_merge)),".",fixed=T),function(l) l[2]))
	
	y <- DGEList(counts= trans_merge, group= sampleGp,genes = row.names(trans_merge), remove.zeros=T)
	#y$samples
	
	#Filter out lowly expressed genes
	keep <- rowSums(cpm(y)>1) >= 2
	y <- y[keep, , keep.lib.sizes=FALSE]
	#y$samples
	
	outFileName<-paste0(sampleSource1,"_",sampleType1,"_Vs_",sampleSource2,"_",sampleType2)

	if(reverseComparison)
	{
		sampleGp<-as.factor(sampleGp)
		revSampleGp<-rev(sampleGp)
		levels(revSampleGp)<-rev(levels(sampleGp))
		design <- model.matrix(~ revSampleGp)
	}
	else{
	design <- model.matrix(~ sampleGp)
	}
	
	y <- estimateDisp(y, design)
	
	fit <- glmFit(y, design)
	lrt.2vs1 <- glmLRT(fit)
	fitTreat<-glmTreat(fit,lfc = 1)

	# Save the R data objects with full calculations
	# about differential expression. 
	
	listDETRObjects<-list()
	
	listDETRObjects[["y"]]<- y
	listDETRObjects[["lrt.2vs1"]]<- lrt.2vs1
	listDETRObjects[["glmFit"]]<- fit
	listDETRObjects[["fitTreat"]]<- fitTreat
		
	saveRDS(listDETRObjects, paste0("../DET_RObjects/list_",outFileName,".rds"))
	
	df_Annot<-writeDETWithHGNCSymbols(lrt.2vs1@.Data[[15]], paste0("../DET_Files/DET_",outFileName,".csv"))

	return(df_Annot)
}


#------------------------------------
#
# Get the list of significantly diffrentially 
# expressed overlapping genes from
# the treat object following the glmFit
#
#------------------------------------
getOverlappingDET<-function(locRObjects,rObj1, rObj2)
{
  fitObject1<-readRDS(paste0(locRObj, rObj1))
  fitObject2<-readRDS(paste0(locRObj, rObj2))
  
  OverlappingDET=intersect(rownames(subset(fitObject1[[4]][[16]], PValue < 0.05)),
                           rownames(subset(fitObject2[[4]][[16]], PValue < 0.05)) )
  
  fitObject1_DET <- dim(subset(fitObject1[[4]][[16]], PValue < 0.05))[1]
  fitObject2_DET <- dim(subset(fitObject2[[4]][[16]], PValue < 0.05))[1]
  
  # We may need to get the original read count for the overlapping DETs
  outFileName=paste0(locRObjects,"OverlappingDET_",gsub(".rds", "", rObj1),"_Vs_",gsub(".rds", "", rObj2),".rds")
  
  if(length(OverlappingDET) > 0)
  {
    print (paste0("No of significant DET ", rObj1," :", fitObject1_DET, " Of ",dim(fitObject1[[4]][[16]])[1]))
    print (paste0("No of significant DET ", rObj2," :", fitObject2_DET, " Of ",dim(fitObject2[[4]][[16]])[1]))
    print (paste0("No of overlapping genes :",length(OverlappingDET)))
    print (paste0("Proportion of overlapping genes for all DET from ", rObj1," :", length(OverlappingDET)/fitObject1_DET))
    print (paste0("Proportion of overlapping genes for all DET from", rObj2," :", length(OverlappingDET)/fitObject2_DET))
    saveRDS(OverlappingDET, outFileName)
  }	
  else{
    print ("No overlap found! ")
  }
  #return (length(OverlappingDET))	
}


# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

#------------------------------------------------------
# Below function will calculate the correlation coefficients 
# along with confidence interval and pValue of significance for
# observed correlation and output to a tab delimited text file
# as CorrelationMatrix_'Comparison".txt. 
# The data (merged transcript 
# read count for the two groups of samples) used to calculate the 
# correlation will be written to a separate file. 
# All the zeros will be converted to NA
# All the transcript with NAs will be removed 
# i.e., only the complete cases will be considered
#------------------------------------------------------

prepareDataToDoReadCountCorrelation<-function(sampleSource1,sampleSource2, sampleType1, sampleType2, reverseComparison=FALSE, removeTrizol=FALSE)
{
  # read transcript data
  transS1 <-readDECountFile(sampleSource1, sampleType1)
  transS2 <-readDECountFile(sampleSource2, sampleType2)
  
  names(transS1) <- paste0(names(transS1),paste0(".", sampleSource1,"_",sampleType1))
  names(transS2) <- paste0(names(transS2),paste0(".", sampleSource2,"_",sampleType2))
  
  trans_merge <- merge(transS1, transS2, 
                       by.x = paste0("transcript_id.", sampleSource1,"_",sampleType1), 
                       by.y = paste0("transcript_id.", sampleSource2,"_",sampleType2),all = TRUE)
  
  rownames(trans_merge)<-trans_merge[,1]
  trans_merge<-trans_merge[,-1]
  
  trans_merge[trans_merge==0] <- NA
  
  df_aux<-trans_merge[which(complete.cases(trans_merge)),]
  
  trans_merge <- df_aux
  
  outFileName<-paste0(sampleSource1,"_",sampleType1,"_Vs_",sampleSource2,"_",sampleType2)
  options(digits = 3)
  
  results<-psych::corr.test(trans_merge, method = "p",ci=TRUE)
  
  corMatrix<-cor(trans_merge[which(complete.cases(trans_merge)),])
  
  # see the function above to calulate the pvalue
  # for the correlation coefficients
  p.mat <- cor.mtest(trans_merge[which(complete.cases(trans_merge)),])
  
  rownames(p.mat)<-c(paste0("$D1[",sampleSource1,"]^Fixed"),paste0("$D2[",sampleSource1,"]^Fixed"),paste0("$D1[",sampleSource1,"]^Unfixed"),paste0("$D2[",sampleSource1,"]^Unfixed"))
  rownames(corMatrix)<-c(paste0("$D1[",sampleSource1,"]^Fixed"),paste0("$D2[",sampleSource1,"]^Fixed"),paste0("$D1[",sampleSource1,"]^Unfixed"),paste0("$D2[",sampleSource1,"]^Unfixed"))
  
  colnames(corMatrix)<-c(paste0("$D1[",sampleSource1,"]^Fixed"),paste0("$D2[",sampleSource1,"]^Fixed"),paste0("$D1[",sampleSource1,"]^Unfixed"),paste0("$D2[",sampleSource1,"]^Unfixed"))
  
  colnames(p.mat)<-c(paste0("$D1[",sampleSource1,"]^Fixed"),paste0("$D2[",sampleSource1,"]^Fixed"),paste0("$D1[",sampleSource1,"]^Unfixed"),paste0("$D2[",sampleSource1,"]^Unfixed"))
  
  # Set the colour palette for the correlation plot
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  # use correlation matrix above to plot the heatmap
  pdf(paste0("../Correlate_ReadCounts/CorrelationPlot_",outFileName,".pdf"))
  corrplot(corMatrix
           , method="color", col=col(200),  
           type="upper", order="alphabet", 
           addCoef.col = "black", # Add coefficient of correlation
           tl.col="black", tl.srt=0, #Text label color and rotation
           # Combine with significance
           p.mat = p.mat, sig.level = 0.01, insig = "blank", 
           tl.offset = 1.25,
           # hide correlation coefficient on the principal diagonal
           diag=FALSE )
  dev.off()
  # Bootstrap ci - incase reviewers likes to see. For the 
  # first submitted version we will using the normal theory to
  # calculate cor and cis 
  #rci <-  psych::cor.ci(trans_merge ,overlap=TRUE, n.iter=10,main="Correct for overlap") 
  #show the confidence intervals
  #ci <- cor.plot.upperLowerCi(rci)
  
  write.table(results[[9]],paste0("../Correlate_ReadCounts/CorrelationMatrix_",outFileName,".txt"),sep = "\t")		
  write.csv(trans_merge, paste0("../Correlate_ReadCounts/ReadCount_",outFileName,".csv"))
}

# <- <- <- <- <- <- <- <- <- <- <- <- <- <-  <- <- <- 
#
# END OF FUNCTIONS USED FOR FOLLOWING ANALYSES
#
# <- <- <- <- <- <- <- <- <- <- <- <- <- <-  <- 



#for MAC=============
#mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice")
#ensembl <- useDataset("hsapiens_gene_ensembl",mart)
#for linux workstation=========
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org", path="/biomart/martservice")
#mart = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="grch37.ensembl.org",ensemblRedirect = FALSE)
ensembl <- useDataset("hsapiens_gene_ensembl",mart)

# <- <- <- <- <- <- <- <- <- <- <- <- <- <- 
#
# Read Input gene and transcript counts
#
# <- <- <- <- <- <- <- <- <- <- <- <- <- <- 

locRObj<-"../DET_RObjects/"
locCountFiles<- "../"

setwd(".")

dir.create("../DET_RObjects", showWarnings = TRUE, recursive = FALSE, mode = "0777")
dir.create("../Correlate_ReadCounts", showWarnings = TRUE, recursive = FALSE, mode = "0777")

#-------------------------------------------------------------
#
# GET DET FOR BAL SAMPLES 
#
#-------------------------------------------------------------

# BULK COMPARISON
DET_BAL_FixVsUF<-getDET("BAL", "BAL", "Fix", "UF",TRUE)

# CELL TYPE SPECIFIC COMPARISON
DET_BAL_T4MVsFix<-getDET("BAL", "BAL", "T4M", "Fix")
DET_BAL_T4MVsUF<-getDET("BAL", "BAL", "T4M", "UF",TRUE)
# To be done after calculating DET
getOverlappingDET(locRObj, "list_BAL_T4M_Vs_BAL_Fix.rds","list_BAL_T4M_Vs_BAL_UF.rds")

DET_BAL_T8MVsFix<-getDET("BAL", "BAL", "T8M", "Fix")
DET_BAL_T8MVsUF<-getDET("BAL", "BAL", "T8M", "UF",TRUE)
# To be done after calculating DET
getOverlappingDET(locRObj, "list_BAL_T8M_Vs_BAL_Fix.rds","list_BAL_T8M_Vs_BAL_UF.rds")

# GET ANNOTATIONS (GENE SYMBOLS)
getHGNCSymbols("BAL", "Fix")
getHGNCSymbols("BAL", "T4M")
getHGNCSymbols("BAL", "T8M")
getHGNCSymbols("BAL", "UF")


#-------------------------------------------------------------
#
# GET DET FOR BLOOD SAMPLES
#
#-------------------------------------------------------------

# BULK COMPARISON
DET_Blood_FixVsUF<-getDET("Blood", "Blood", "Fix", "UF",TRUE)

# CELL TYPE SPECIFIC COMPARISON
DET_Blood_T4MVsFix<-getDET("Blood", "Blood", "T4M", "Fix")
DET_Blood_T4MVsUF<-getDET("Blood", "Blood", "T4M", "UF",TRUE)
# to be done after calculating DET
getOverlappingDET(locRObj, "list_Blood_T4M_Vs_Blood_Fix.rds","list_Blood_T4M_Vs_Blood_UF.rds")

DET_Blood_T8MVsFix<-getDET("Blood", "Blood", "T8M", "Fix")
DET_Blood_T8MVsUF<-getDET("Blood", "Blood", "T8M", "UF",TRUE)
# to be done after calculating DET
getOverlappingDET(locRObj, "list_Blood_T8M_Vs_Blood_Fix.rds","list_Blood_T8M_Vs_Blood_UF.rds")

# GET ANNOTATIONS (GENE SYMBOLS)
getHGNCSymbols("Blood", "UF")
getHGNCSymbols("Blood", "Fix")
getHGNCSymbols("Blood", "T4M")
getHGNCSymbols("Blood", "T8M")

# to be done after calculating DET

#---------------

# GET Correlation Plots
prepareDataToDoReadCountCorrelation("BAL", "BAL", "Fix", "UF",TRUE)
prepareDataToDoReadCountCorrelation("Blood", "Blood", "Fix", "UF",TRUE)



