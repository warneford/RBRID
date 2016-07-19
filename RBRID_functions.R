# RBRID Functions 
# Author: Robert Warneford-Thomson
# Last Updated 18Jul2016

## Function to normalize and analyze RBR peptide intensities. MinusCols and plusCols are numeric vectors defining, respectively, columns containing -4sU and +4sU replicates.
# If cutoff = TRUE peptides are removed from final output that are not below p value cutoff or not depleted to specified factor
NormRBR <-function(Data = Prot_table, minusCols = c(8:9), plusCols = c(12:13), pvalcutoff = 0.1, depletionfactor = 0.95, cutoff = TRUE)
{
 
  MinusSamples <- names(Data[minusCols])
  PlusSamples <- names(Data[plusCols])
  RawSamples <- names(Data[c(minusCols, plusCols)])
  
  # Rename data columns
  colnames(Data)[grepl("contaminant", colnames(Data))] <- "contaminant"
  colnames(Data)[grepl("Proteins", colnames(Data))] <- "uniprotID"
  colnames(Data)[grepl("Protein.names", colnames(Data))] <- "name"
  colnames(Data)[grepl("Gene.names", colnames(Data))] <- "symbol"
  colnames(Data)[grepl("Unique..Groups.", colnames(Data))] <- "unique"
  
  # Convert Contaminants to logical vector
  Data$contaminant <- as.character(Data$contaminant)
  Data$contaminant[Data$contaminant == "+"] <- T
  Data$contaminant[Data$contaminant != "TRUE"] <- F
  Data$contaminant<-as.logical(Data$contaminant)
  Data$uniprotID <- as.character(Data$uniprotID)
  Data$symbol <- as.character(Data$symbol)
  Data$name <- as.character(Data$name)
  
  # Convert Unique to logical vector
  Data$unique <- as.character(Data$unique)
  Data$unique[Data$unique =="yes"]<-T
  Data$unique[Data$unique !="TRUE"]<-F
  Data$unique<-as.logical(Data$unique)
  
  # Assigns unique peptide ID to each peptide 
  Data[,RawSamples] <- as.numeric(unlist(Data[,c(minusCols, plusCols)]))
  Data <- cbind(Pept_ID = c(1:nrow(Data)), Data)
  
  # Normalize peptide intensity by Sum of total intensity (removing contaminants) in each replicate AKA set Euclidean length = 1
  
  Norm_data <- sweep(Data[,RawSamples], 2, colSums(Data[!Data$contaminant, RawSamples]), FUN = "/")
  names(Norm_data) <- paste0(RawSamples,"_Norm")
  Norm_data[Norm_data == 0] <- NA # Remove zero values
  Data <- cbind(Data, Norm_data)
  
  # Update sample names
  MinusSamples_Norm <- paste0(MinusSamples, "_Norm")
  PlusSamples_Norm <- paste0(PlusSamples, "_Norm")

  # Compute averages of biological replicates
  Data$minus_mean <- rowMeans(Data[,MinusSamples_Norm], na.rm = T)
  Data$plus_mean <- rowMeans(Data[,PlusSamples_Norm], na.rm = T)
  Data$minus_mean[is.nan(Data$minus_mean)] <- NA
  Data$plus_mean[is.nan(Data$plus_mean)] <- NA
  
  # Calculate +4su/-4sU log2-fold change
  Data$log2fold <- log2(Data$plus_mean/Data$minus_mean)
  Data$log2fold[is.infinite(Data$log2fold) | is.nan(Data$log2fold)] <- NA
  
  # Calculate p values from two-sided t-test
  Data$p.value <- apply(Data[,c(PlusSamples_Norm, MinusSamples_Norm)], 1, function(x) 
    

  { if (sum(!is.na(x[c(PlusSamples_Norm)])) < 2 | sum(!is.na(x[c(MinusSamples_Norm)])) < 2) {return(1)}
    
  else { t.test(x[PlusSamples_Norm], x[MinusSamples_Norm], paired = F, alternative = "two.sided")$p.value}})
  
  # Calculate RBR-ID score
  Data$RBRscore <- log(Data$p.value)*Data$log2fold
  
  # Apply Cut-offs 
  if (cutoff == T) {
  # Remove data with p value > 0.1
  Data <- Data[Data$p.value < pvalcutoff & !is.na(Data$p.value),]
  
  # Include data with depletion greater than 5%
  Data <- Data[Data$plus_mean/Data$minus_mean < depletionfactor ,]
  return(Data)}
  else {return(Data)}}
  
  

subsetID <- function(names = c("P19909"), category = col, data = df)
  # Function to extract rows that match desired value of string variable, e.g. specific uniprot IDs.
  # names is vector of strings listing desired values to subset, category gives column to search for desired names
{ match.list <-lapply(names, function(x) {grep(pattern = paste0(".*", x,".*$"), data[,category], value = FALSE)})
match.rows <- unlist(match.list)
return(match.rows)}


calculateCoverage.apply<-function(range,peps)
{
  hits<-subsetByOverlaps(peps,range)
  covered.regions<-reduce(hits)
  sum(width(covered.regions))/end(range)
}

residueScorePlot<-function(range,peps,smooth=T, show.coverage=T, no.negs=F,span=0.0001, return.values=F)
  ## Draws score plot, smoothing and white space on NA peptides can be controlled by options
  ## LAST UPDATED 4.1.2016
{
  score<-rep(0,end(range))
  hits<-subsetByOverlaps(peps,range)
  for(i in 1:length(hits)){
    pep.pos<-start(hits[i]):end(hits[i])
    score[pep.pos]=score[pep.pos]+hits$score[i]
  }
  gaps<-setdiff(range,reduce(hits))
  for(i in 1:length(gaps)){
    gap.pos<-start(gaps[i]):end(gaps[i])
    score[gap.pos]<-NA
  }
  score.mod<-score
  score.mod[is.na(score.mod)]<-0
  if(no.negs) score.mod[score.mod<0]<-0
  if(smooth){
    values<-supsmu(1:length(score.mod),score.mod,span=span)
  } else {
    values<-(list(x=1:length(score.mod),y=score.mod))
  }
  plot(values, type="l", main=paste0(range$symbol, " (",seqnames(range),")"), xlab="residue #", ylab="RBR-ID score")
  if(show.coverage) points(c(1:length(score))[is.na(score)],rep(0,sum(is.na(score))), pch=16, col="white", cex=500/end(range))
  if(return.values) return(values)
}

residueScore<-function(range,peps,na=T)
  ## Returns a vector with score per residue
{
  score<-rep(0,end(range))
  hits<-subsetByOverlaps(peps,range)
  for(i in 1:length(hits)){
    pep.pos<-start(hits[i]):end(hits[i])
    score[pep.pos]=score[pep.pos]+hits$score[i]
  }
  if(na){
    gaps<-setdiff(range,reduce(hits))
    for(i in 1:length(gaps)){
      gap.pos<-start(gaps[i]):end(gaps[i])
      score[gap.pos]<-NA
    }
  }
  score
}

binMean <- function(data, bins)
  ## Calculate mean in a bin
{
  total<-length(data)
  internal.bins<-bins-2
  binSize<-floor(total/internal.bins)
  extra<-floor((total-binSize*internal.bins)/2)
  binned<-rep(0,bins)
  binned[1]<-mean(data[1:extra],na.rm=T)
  for(i in 1:(internal.bins)){
    binStart<-(i*binSize-binSize+1)+extra
    binEnd<-(i*binSize)+extra
    binned[i+1]<-mean(data[binStart:binEnd], na.rm=T)
    # print(paste(binStart,binEnd,binned[i+1], sep="  "))
  }
  binned[bins]<-mean(data[binEnd:total], na.rm=T)
  return(binned)
}

listScoreplot <- function (Proteins, PDB, Range = prot.gr, peps = peps, Experiment = Experiment, Descriptor = Descriptor, Date = Date)
  # This is a wrapper function for residueScorePlot which will take a list of proteins, plot their RBRID scores, and output the scores
  # as an attribute file for use with the Chimera protein visualization software for a desired .pdb file
  # Experiment = experiment code
  # Descriptor = descriptive info about experiment
  # 
{
  lapply(Proteins, function(candidate, ...) {
    
    if (!any(grepl(candidate, Range$symbol))) {return(paste0("Candidate ", candidate, " is not identified"))} 
    else {
      pos<-which(Range$symbol==candidate)
      t<-residueScorePlot(Range[pos,],peps,smooth=T, no.negs=T, return.values=T, span=0.001)
      
      a<-as.data.frame(peps[seqnames(peps)==as.character(seqnames(Range)[pos]),])
      a<-a[order(a$score),]
      
      # to make score file that gets inserted in PDB for coloring
      t$label<-paste0(":",t$x)
      RBRscoretable<- cbind(rep("", length(t$x)), t$label, round(t$y*10)/10)
      
      fileName <- paste0("~/analyses/", Experiment, "/", Descriptor, Date, candidate, "scoreSmooth_", PDB, ".txt")
      write.table(RBRscoretable,sep="\t",quote=F, file= fileName, row.names = F, col.names = F) 
      
      header <- c(paste0("#  PDB entry ", PDB, " - " ,candidate, " RBR-ID RNA binding prediction score"),
                  "#  (4sU crosslinked peptide mass spectrometry depletion)", 
                  "#  Method credit to Chongsheng He & Roberto Bonasio",
                  paste0("#", "  Experiment: ", Experiment),
                  paste0("#", "  Details: ", Descriptor),
                  "#  Use this file to assign the attribute to crystal structures in in Chimera with the",
                  "#  Define Attribute tool or the command defattr.",
                  paste0("#  Make sure to specify attribute only to ", candidate, " chain"),
                  "#",
                  paste0("attribute: rbrscore", candidate),
                  "match mode: 1-to-1",
                  "recipient: residues")
      
      fConn <- file(fileName)
      Lines <- readLines(fConn) 
      writeLines(c(header, Lines), sep = "\n" , con = fConn)
      close(fConn)
    }})}