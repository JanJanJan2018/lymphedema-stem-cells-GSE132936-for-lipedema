

library(rvest)
library(lubridate)
library(dplyr)

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






