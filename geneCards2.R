# 
# # These functions grab the genes and the gene summaries from genecards.org
# # and some calculate the mean and median fold change values across 
# # samples of treatment/control or diseased/healthy etc.
# 
# # To extract the genes:
# # - find25genes(protein) will grab the 25 genes associated with the protein from web
# # - getProteinGenes(protein) will print the genes associated with the protein
# # - getSummaries(gene, protein) will grab the gene protein summaries from web
# # -getGeneSummaries(protein) will print the gene summary of protein gene
# # 


# This script can be used by calling
# source('geneCards.R') in other script inside same folder. The table are specific
# to the functions to get the fold change values and in this folder.

library(rvest)
library(lubridate)
library(dplyr)

Gene_Path <- './gene scrapes'

if (dir.exists(Gene_Path)){
  unlink(Gene_Path, recursive=TRUE)
  dir.create(Gene_Path)
} else {
  dir.create(Gene_Path)
}

find25genes <- function(protein){
  
  url <- 'https://www.genecards.org/Search/Keyword?queryString=protein'
  
  protein <- as.character(protein)
  protein <- tolower(protein)
  protein <- gsub(' ','%20',protein)
  
  url <- as.character(url)
  url <- gsub('protein',protein, url)
  
  #  webpage <- read_html(url,options=pedantic,encoding = "UTF-8")
  
  webpage <- read_html(url,encoding = "UTF-8")
  
  protein_html <- html_nodes(webpage,'.symbol-col a')
  protein1 <- html_text(protein_html)
  
  Protein <- as.data.frame(protein1)
  colnames(Protein) <- 'proteinType'
  Protein$proteinType <- as.character(paste(Protein$proteinType))
  Protein$proteinType <- gsub('\n','',Protein$proteinType)
  
  
  date <- as.data.frame(rep(date(),length(Protein$proteinType)))
  colnames(date) <- 'todaysDate'
  
  protein2 <- gsub('%20','-',protein)
  
  proteinName <- as.data.frame(rep(protein2,length(Protein$proteinType)))
  colnames(proteinName) <- 'proteinSearched'
  
  tableProtein <- cbind(Protein,proteinName,date)
  
  setwd(Gene_Path)
  
  
  write.table(tableProtein, 
              paste(protein2,".csv",sep=''), append=TRUE,
              col.names=FALSE, sep=",", quote=TRUE,qmethod="double",
              row.names=FALSE)
  names <- colnames(tableProtein)
  write.csv(names,paste('tableProteinHeader_',protein2,'.csv',sep=''),row.names=FALSE)
  
  setwd('../')
}


#find25genes('estrogen')


getProteinGenes <- function(protein){
  protein <- as.character(protein)
  protein <- tolower(protein)
  protein <- gsub(' ','-',protein)
  table <- read.csv(paste(Gene_Path,'/',protein,'.csv',sep=''),sep=',',
                    header=F,na.strings=c('',' ','NA'), stringsAsFactors = F)
  header <- read.csv(paste(Gene_Path,'/tableProteinHeader_',protein,'.csv',sep=''),
                     sep=',', header=T, na.strings=c('',' ','NA'), stringsAsFactors = F)
  names <- header$x
  colnames(table) <- names
  fileName <- paste('Top25',protein,'s.csv',sep='')
  write.csv(table, fileName, row.names=FALSE)
  return(list(table,fileName))
}


#getProteinGenes('estrogen')



getSummaries <- function(gene,protein){
  url <- 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE&keywords=protein'
  
  protein <- as.character(protein)
  protein <- tolower(protein)
  protein <- gsub(' ',',',protein)
  gene <- as.character(gene)
  gene <- tolower(gene)
  
  url <- as.character(url)
  url <- gsub('GENE',gene,url)
  url <- gsub('protein',protein, url)

  #webpage <- read_html(url,options=pedantic,encoding = "UTF-8")
  
  
  webpage <- read_html(url,encoding = "UTF-8")
  
  Entrez_html <- html_nodes(webpage, '.gc-section-header+ .gc-subsection p')
  Entrez <- html_text(Entrez_html) 
  
  GeneCards_html <- html_nodes(webpage, '.gc-subsection-header+ p')
  GeneCards <- html_text(GeneCards_html) 
  
  UniProt_html <- html_nodes(webpage, '#summaries li:nth-child(1) div')
  UniProtKB <- html_text(UniProt_html) 
  
  Entrez0 <- ifelse(length(Entrez)==0, 'no summary',as.character(paste(Entrez)))
  Entrez1 <- as.data.frame(Entrez0)
  colnames(Entrez1) <- 'EntrezSummary'
  
  GeneCards0 <- ifelse(length(GeneCards)==0,'no summary',
                       as.character(paste(GeneCards)))
  GeneCards1 <- as.data.frame(GeneCards0)
  colnames(GeneCards1) <- 'GeneCardsSummary'
  
  UniProtKB0 <- ifelse(length(UniProtKB)==0,'no summary',
                       as.character(paste(UniProtKB)))
  UniProtKB1 <- as.data.frame(UniProtKB0)
  colnames(UniProtKB1) <- 'UniProtKB_Summary'
  
  Entrez1$EntrezSummary <- as.character(paste(Entrez1$EntrezSummary))
  Entrez1$EntrezSummary <- gsub('\n','',Entrez1$EntrezSummary)
  
  GeneCards1$GeneCardsSummary <- as.character(paste(GeneCards1$GeneCardsSummary))
  GeneCards1$GeneCardsSummary <- gsub('\n','',GeneCards1$GeneCardsSummary)
  
  UniProtKB1$UniProtKB_Summary <- as.character(paste(UniProtKB1$UniProtKB_Summary))
  UniProtKB1$UniProtKB_Summary <- gsub('\n','',UniProtKB1$UniProtKB_Summary)
  
  date <- as.data.frame(rep(date(),length(Entrez1$EntrezSummary)))
  colnames(date) <- 'todaysDate'
  
  protein2 <- gsub(',','-',protein)
  
  proteinName <- as.data.frame(rep(protein2,length(Entrez1$EntrezSummary)))
  colnames(proteinName) <- 'proteinSearched'
  
  gene <- as.data.frame(rep(toupper(gene),length(Entrez1$EntrezSummary)))
  colnames(gene) <- 'gene'
  
  tableProtein <- cbind(proteinName,gene,Entrez1,GeneCards1,UniProtKB1,date)
  
  setwd(Gene_Path)
  
  
  write.table(tableProtein, 
              paste(protein2,"summary.csv",sep=''), append=TRUE,
              col.names=FALSE, sep=",", quote=TRUE,qmethod="double",
              row.names=FALSE)
  names <- colnames(tableProtein)
  write.csv(names,paste('geneHeader_summary_',protein2,'.csv',sep=''),row.names=FALSE)
  
  setwd('../')
  return(gene)
}


#getSummaries('TP53','estrogen')


getGeneSummaries <- function(protein){
  protein <- as.character(protein)
  protein <- tolower(protein)
  protein <- gsub(' ','-',protein)
  
  table <- read.csv(paste(Gene_Path,'/',protein,'summary.csv',sep=''),
                    sep=',',header=F,na.strings=c('',' ','NA'), stringsAsFactors = F)
  
  header <- read.csv(paste(Gene_Path,'/geneHeader_summary_',protein,'.csv',sep=''),
                     sep=',', header=T, na.strings=c('',' ','NA'), stringsAsFactors = F)
  names <- header$x
  colnames(table) <- names
  table$EnsemblID <- trimws(as.character(paste(table$EnsemblID)), 
                            which='right',  whitespace='\n')
  fileName <- paste('proteinGeneSummaries_',protein,'.csv',sep='')
  write.csv(table, fileName, row.names=FALSE)
  return(list(table,fileName))
}


#getGeneSummaries('estrogen')

#added the Ensembl ID to the data frame in the getSummaries() to be getSummaries2()

#get the Ensembl ID.

getSummaries2 <- function(gene,protein){
  url <- 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE&keywords=protein'
  
  protein <- as.character(protein)
  protein <- tolower(protein)
  protein <- gsub(' ',',',protein)
  gene <- as.character(gene)
  gene <- tolower(gene)
  
  url <- as.character(url)
  url <- gsub('GENE',gene,url)
  url <- gsub('protein',protein, url)
  
# webpage <- read_html(url,options=pedantic,encoding = "UTF-8")
  webpage <- read_html(url,encoding = "UTF-8")
  
  Entrez_html <- html_nodes(webpage, '.gc-section-header+ .gc-subsection p')
  Entrez <- html_text(Entrez_html) 
  
  GeneCards_html <- html_nodes(webpage, '.gc-subsection-header+ p')
  GeneCards <- html_text(GeneCards_html) 
  
  UniProt_html <- html_nodes(webpage, '#summaries li:nth-child(1) div')
  UniProtKB <- html_text(UniProt_html) 
  
  ensemblhtml <- html_nodes(webpage,'#aliases_descriptions .gc-subsection:nth-child(2) .gc-subsection-inner-wrap')
  ensemblID <- html_text(ensemblhtml)
  
  ensemble0 <- ifelse(length(ensemblID)==0,'no Ensembl_ID',
                      as.character(paste(ensemblID)))
  ensemble1 <- as.data.frame(ensemble0)
  colnames(ensemble1) <- 'EnsemblID'

  
  ##########
  ensemble1$EnsemblID <- trimws(ensemble1$EnsemblID, 
                                which='left',whitespace=' ')
  ensemble1$EnsemblID <- trimws(ensemble1$EnsemblID, 
                                which='right',whitespace=' ')
  
  split1 <- strsplit(ensemble1$EnsemblID,split='\\n',perl=T)
  df <- as.data.frame(split1)
  colnames(df) <- c('EnsemblID')
  df_ensemble <- df[4,]
  df_ensemble <- as.data.frame(df_ensemble)
  colnames(df_ensemble) <- 'EnsemblID'
  row.names(df_ensemble) <- NULL
  df_ensemble$EnsemblID <- trimws(df_ensemble$EnsemblID,which='left',
                                  whitespace=' ')
  df_ensemble$EnsemblID <- gsub('Ensemble: ','',df_ensemble$EnsemblID)
  ensemble1 <- df_ensemble
  ensemble0 <- ifelse(length(ensemble1$EnsemblID)==0,'no Ensembl_ID',
                      as.character(paste(ensemble1$EnsemblID)))
  ensemble1 <- as.data.frame(ensemble0)
  colnames(ensemble1) <- 'EnsemblID'
  ##########
  
  Entrez0 <- ifelse(length(Entrez)==0, 'no summary',as.character(paste(Entrez)))
  Entrez1 <- as.data.frame(Entrez0)
  colnames(Entrez1) <- 'EntrezSummary'
  
  GeneCards0 <- ifelse(length(GeneCards)==0,'no summary',
                       as.character(paste(GeneCards)))
  GeneCards1 <- as.data.frame(GeneCards0)
  colnames(GeneCards1) <- 'GeneCardsSummary'
  
  UniProtKB0 <- ifelse(length(UniProtKB)==0,'no summary',
                       as.character(paste(UniProtKB)))
  UniProtKB1 <- as.data.frame(UniProtKB0)
  colnames(UniProtKB1) <- 'UniProtKB_Summary'
  
  ensemble1$EnsemblID <- as.character(paste(ensemble1$EnsemblID))
  ensemble1$EnsemblID <- gsub('Ensembl: ','',
                              ensemble1$EnsemblID)
  ensemble1$EnsemblID <- trimws(as.character(paste(ensemble1$EnsemblID)),
                                which='right',
                                whitespace='\n')
  Entrez1$EntrezSummary <- as.character(paste(Entrez1$EntrezSummary))
  Entrez1$EntrezSummary <- gsub('\n','',Entrez1$EntrezSummary)
  
  GeneCards1$GeneCardsSummary <- as.character(paste(GeneCards1$GeneCardsSummary))
  GeneCards1$GeneCardsSummary <- gsub('\n','',GeneCards1$GeneCardsSummary)
  
  UniProtKB1$UniProtKB_Summary <- as.character(paste(UniProtKB1$UniProtKB_Summary))
  UniProtKB1$UniProtKB_Summary <- gsub('\n','',UniProtKB1$UniProtKB_Summary)
  
  date <- as.data.frame(rep(date(),length(Entrez1$EntrezSummary)))
  colnames(date) <- 'todaysDate'
  
  protein2 <- gsub(',','-',protein)
  
  proteinName <- as.data.frame(rep(protein2,length(Entrez1$EntrezSummary)))
  colnames(proteinName) <- 'proteinSearched'
  
  gene <- as.data.frame(rep(toupper(gene),length(Entrez1$EntrezSummary)))
  colnames(gene) <- 'gene'
  
  tableProtein <- cbind(proteinName,gene,ensemble1,Entrez1,GeneCards1,UniProtKB1,date)
  
  setwd(Gene_Path)
  
  
  write.table(tableProtein, 
              paste(protein2,"summary.csv",sep=''), append=TRUE,
              col.names=FALSE, sep=",", quote=TRUE,qmethod="double",
              row.names=FALSE)
  names <- colnames(tableProtein)
  write.csv(names,paste('geneHeader_summary_',protein2,'.csv',sep=''),row.names=FALSE)
  
  setwd('../')
  
  return(gene)
}




# from the monotonicGenes.Rmd of GSE152418 analysis on COVID-19
# tests four classees for accuracy returning precision and recall and accuracy
precisionRecallAccuracy <- function(df){
  
  colnames(df) <- c('pred','type')
  df$pred <- as.character(paste(df$pred))
  df$type <- as.character(paste(df$type))
  
  classes <- unique(df$type)
  
  class1a <- as.character(paste(classes[1]))
  class2a <- as.character(paste(classes[2]))
  class3a <- as.character(paste(classes[3]))
  class4a <- as.character(paste(classes[4]))
  
  #correct classes
  class1 <- subset(df, df$type==class1a)
  class2 <- subset(df, df$type==class2a)
  class3 <- subset(df, df$type==class3a)
  class4 <- subset(df, df$type==class4a)
  
  #incorrect classes
  notClass1 <- subset(df,df$type != class1a)
  notClass2 <- subset(df,df$type != class2a)
  notClass3 <- subset(df,df$type != class3a)
  notClass4 <- subset(df, df$type != class4a)
  
  #true positives (real positives predicted positive)
  tp_1 <- sum(class1$pred==class1$type)
  tp_2 <- sum(class2$pred==class2$type)
  tp_3 <- sum(class3$pred==class3$type)
  tp_4 <- sum(class4$pred==class4$type)
  
  #false positives (real negatives predicted positive)
  fp_1 <- sum(notClass1$pred==class1a)
  fp_2 <- sum(notClass2$pred==class2a)
  fp_3 <- sum(notClass3$pred==class3a)
  fp_4 <- sum(notClass4$pred==class4a)
  
  #false negatives (real positive predicted negative)
  fn_1 <- sum(class1$pred!=class1$type)
  fn_2 <- sum(class2$pred!=class2$type)
  fn_3 <- sum(class3$pred!=class3$type)
  fn_4 <- sum(class4$pred!=class4$type)
  
  #true negatives (real negatives predicted negative)
  tn_1 <- sum(notClass1$pred!=class1a)
  tn_2 <- sum(notClass2$pred!=class2a)
  tn_3 <- sum(notClass3$pred!=class3a)
  tn_4 <- sum(notClass4$pred!=class4a)
  
  
  #precision
  p1 <- tp_1/(tp_1+fp_1)
  p2 <- tp_2/(tp_2+fp_2)
  p3 <- tp_3/(tp_3+fp_3)
  p4 <- tp_4/(tp_4+fp_4)
  
  p1 <- ifelse(p1=='NaN',0,p1)
  p2 <- ifelse(p2=='NaN',0,p2)
  p3 <- ifelse(p3=='NaN',0,p3)
  p4 <- ifelse(p4=='NaN',0,p4)
  
  #recall
  r1 <- tp_1/(tp_1+fn_1)
  r2 <- tp_2/(tp_2+fn_2)
  r3 <- tp_3/(tp_3+fn_3)
  r4 <- tp_4/(tp_4+fn_4)
  
  r1 <- ifelse(r1=='NaN',0,r1)
  r2 <- ifelse(r2=='NaN',0,r2)
  r3 <- ifelse(r3=='NaN',0,r3)
  r4 <- ifelse(r4=='NaN',0,r4)
  
  #accuracy
  ac1 <- (tp_1+tn_1)/(tp_1+tn_1+fp_1+fn_1)
  ac2 <- (tp_2+tn_2)/(tp_2+tn_2+fp_2+fn_2)
  ac3 <- (tp_3+tn_3)/(tp_3+tn_3+fp_3+fn_3)
  ac4 <- (tp_4+tn_4)/(tp_4+tn_4+fp_4+fn_4)
  
  table <- as.data.frame(rbind(c(class1a,p1,r1,ac1),
                               c(class2a,p2,r2,ac2),
                               c(class3a,p3,r3,ac3),
                               c(class4a,p4,r4,ac4)))
  
  colnames(table) <- c('class','precision','recall','accuracy')
  acc <- (sum(df$pred==df$type)/length(df$type))*100
  cat('accuracy is: ',as.character(paste(acc)),'%')
  return(table)
  
  
}



transcriptList <- c('ENST00000000233','ENST00000000412','ENST00000000442','ENST00000001008')

Gene_Path <- './transcripts'

if (dir.exists(Gene_Path)){
  unlink(Gene_Path, recursive=TRUE)
  dir.create(Gene_Path)
} else {
  dir.create(Gene_Path)
}


transcript2Gene <- function(transcript){
  
  url <- 'https://www.genecards.org/Search/Keyword?queryString=transcript'
  
  #transcript are with the ENST prepended label ENST00000000233
  
  transcript <- as.character(transcript)
  #transcript <- tolower(transcript)
  
  url <- as.character(url)
  url <- gsub('transcript',transcript,url)
  
  webpage <- read_html(url)
  
  transcript_html <- html_nodes(webpage, '.symbol-col a')
  transcript1 <- html_text(transcript_html)
  
  transcript0 <- ifelse(length(transcript1)==0,'no ID',
                      as.character(paste(transcript1)))
  Transcript <- as.data.frame(transcript0)
  colnames(Transcript) <- 'transcriptGene'

  Transcript$transcriptGene <- as.character(paste(Transcript$transcriptGene))
  
  searchedTranscript <- as.data.frame(transcript)
  colnames(searchedTranscript) <- 'searchedTranscript'
  searchedTranscript$searchedTranscript <- as.character(paste(searchedTranscript$searchedTranscript))
  
  table <- cbind(searchedTranscript,Transcript)
  names <- colnames(table)
  
  setwd(Gene_Path)
  
  write.table(table, 
              'geneTranscriptTable.csv', append=TRUE,
              col.names=FALSE, sep=",", quote=TRUE,qmethod="double",
              row.names=FALSE)
  write.csv(names,'transcriptHeader.csv',row.names=FALSE)
  
  setwd('../')
  
}

for (i in transcriptList){
  transcript2Gene(i)
}


transcript2GeneTable <- function(){
  
  transcriptGenes <- read.csv('transcripts/geneTranscriptTable.csv', header=F)
  headerTranscripts <- read.csv('transcripts/transcriptHeader.csv')
  colnames(transcriptGenes) <- headerTranscripts$x
  
  write.csv(transcriptGenes,'TableOfTranscriptsAndGenes.csv',row.names=F)
  transcriptTable <- read.csv('TableOfTranscriptsAndGenes.csv', header=T)
  return(transcriptTable)
}

transcript2GeneTable()






