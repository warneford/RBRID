#************************************************************#
#----------------- BA-55 residue score   --------------------#
#************************************************************#

### Author: Roberto Bonasio
### This script is to look at residue-level score plots
### Created: 3.17.2016
### UPDATES:
### - 3.22.2016: replaced long calculations with intermediate files

##### INITIALIZATION ####

MACHINE="brainberg"
homed=Sys.getenv("HOME")
wd=paste0(homed,"/analyses/BA55")

#### LIBRARIES ####
require(GenomicRanges)
require(parallel)
require(gplots)

#### FUNCTIONS ####
color.palette <- function(steps, n.steps.between=NULL, ...)
  #This is a wrapper function for colorRampPalette. It allows for the
  #definition of the number of intermediate colors between the main colors.
  #Using this option one can stretch out colors that should predominate
  #the palette spectrum. Additional arguments of colorRampPalette can also
  #be added regarding the type and bias of the subsequent interpolation.
{ 
  if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
  if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")
  
  fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
  RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
  RGB[,fill.steps] <- col2rgb(steps)
  
  for(i in which(n.steps.between>0)){
    col.start=RGB[,fill.steps[i]]
    col.end=RGB[,fill.steps[i+1]]
    for(j in seq(3)){
      vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]  
      RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
    }
  }
  
  new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
  pal <- colorRampPalette(new.steps, ...)
  return(pal)
}

calculateScore<-function(values,test.string="312$")
  ## calculates RBR-ID score given a matrix with named columns
  ## and strings to match with grep for the samples of interest
{
  test.pval<-grepl(test.string, names(values)) & grepl("pval", names(values))
  test.logfold<-grepl(test.string, names(values)) & grepl("logFold", names(values))
  -values[,test.logfold]*(log10(values[,test.pval]))^2
}

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

generateDomainMatrix<-function(IPR.f,scoreList,aa=20,bins=30)
  ## Generate matrix with aveage bin score over a domain and neighboring regions
  ## Parameters: 
  ## - IPR.f is a data.frame with uniprotIDs start and end of the domain(s) of interest
  ## - scoreList is a residue-level score list of proteins with uniprotIDs as names
  ## - aa is the number of aa to extend before and after
  ## - bins is the number of bins to split the domain into (domains must be binned because they are different sizes)
{
  # create empty matrices
  before.matrix<-matrix(nrow=nrow(IPR.f),ncol=aa)
  domain.matrix<-matrix(nrow=nrow(IPR.f),ncol=bins)
  after.matrix<-matrix(nrow=nrow(IPR.f),ncol=aa)
  # populate matrices
  for(i in 1:nrow(IPR.f)){
    score<-scoreList[[IPR.f$uniprotID[i]]]
    start<-IPR.f$IPR_start[i]
    end<-IPR.f$IPR_end[i]
    # first deal with the pre-domain region
    if(start>1){
      if(start>aa) {
        before.matrix[i,]<-score[(start-aa):(start-1)]
      }
      else{
        before.matrix[i,(aa+2-start):aa]<-score[1:(start-1)]
      }
    }
    # then post-domain region
    if(end<length(score)){
      if(length(score)>end+aa){
        after.matrix[i,]<-score[(end+1):(end+aa)]
      }
      else{
        after.matrix[i,(aa+1-((length(score)-end))):aa]<-score[(end+1):length(score)]
      }
    }
    score.slice<-score[start:end]
    domain.matrix[i,]<-binMean(score.slice,bins)
  }
  
  cbind(before.matrix,domain.matrix,after.matrix)
}

metaDomainPlot<-function(combo.matrix,aa=20,bins=30, smooth=T)
{
  if(smooth){
    plot(supsmu(1:ncol(combo.matrix),colMeans(combo.matrix,na.rm=T)), type="l", lwd=3, xaxt="n", xlab="", ylab="RBR-ID score")
  }
  else{
    plot(1:ncol(combo.matrix),colMeans(combo.matrix,na.rm=T), type="l", lwd=3, xaxt="n", xlab="", ylab="RBR-ID score")
  }
  abline(v=c(aa,aa+bins), lty=2, col="gray")
  axis(1,at=c(0,aa,aa+bins,aa+bins+aa), labels=c(paste0("-",aa," aa"),"N","C",paste0("+",aa," aa")))
}

#### CONSTANTS ####

RRM="IPR000504"
KH="IPR004087"
dsRBD="IPR014720"
RBDs<-c(RRM,KH,dsRBD)

#### FILES ####
# annotations
## interpro
mouse_IPR<-read.delim("~/genomes/Mmus/uniprot/150906.protein2ipr.mouse", header=FALSE, stringsAsFactors=FALSE)
names(mouse_IPR)<-c("uniprotID","iprID","IPR_name","IPR_signature","IPR_start","IPR_end")
#uniprot
lengths<-read.delim(paste0(homed,"/genomes/Mmus/uniprot/160317.uniprot.mouse.lengths"))

# pval and logfold data
pl<-readRDS("output/pval/160317.2.BA55.pval.50runs.3.1.2.sum.0.paired.rds") # load pval & logFold matrix
pl$score.312<-calculateScore(pl,"312$")
pl$score.365<-calculateScore(pl,"365$")
pl$score.noUV<-calculateScore(pl,"noUV$")

# create unique list of uniprotIDs from peptide list
anno<-read.delim("output/peptideAnnotations/1.2.1.PC50.E14.depletion.50runs.peptideAnnotation.160304.tsv", fill=T, quote="", stringsAsFactors = F)		# loads annotation from peptideAnnotation table
prot.anno<-anno[!anno$contaminant,c("uniprotID","symbol","name")]
prot.anno<-prot.anno[order(prot.anno$name, decreasing=T),]
prot.anno<-prot.anno[!duplicated(prot.anno$uniprotID),]
prot.anno<-prot.anno[order(prot.anno$uniprotID),]
prot.anno<-merge(prot.anno, lengths, by="uniprotID")
# convert to grange
prot.gr<-GRanges(seqnames=Rle(prot.anno$uniprotID), ranges=IRanges(start=1, end=prot.anno$length))
prot.gr$symbol<-prot.anno$symbol
# remove problematic proteins (annotation changed between december 15 and march 16)
problem.proteins<-c("D3YU33")
problem.proteins.mask<-!seqnames(prot.gr) %in% problem.proteins
prot.gr<-prot.gr[problem.proteins.mask,]

# prepare peptide grange
pep.gr<-GRanges(seqnames=Rle(anno$uniprotID),ranges=IRanges(start=anno$start,end=anno$end))
pep.gr$pepID<-anno$pepID; pep.gr$unique<-anno$unique; pep.gr$contaminant<-anno$contaminant; pep.gr$sequence<-anno$sequence

#### MAIN SCRIPT ####

# calculate coverage per protein
{
  # rm(list=c("end","start")) # this is because the worker nodes behave strangely otherwise
  # cl<-makeCluster(15)
  # clusterExport(cl, list("subsetByOverlaps","reduce","width","end")) # must pass non-conventional functions to worker nodes
  # coverage<-parSapply(cl=cl,prot.gr,calculateCoverage.apply, peps=pep.gr) # this takes about 1'15" with 30 nodes
  # stopCluster(cl)
  # prot.gr$coverage<-coverage*100
  # saveRDS(prot.gr,"intermediate_files/prot.gr.coverage.160323.rds")
  prot.gr<-readRDS("intermediate_files/prot.gr.coverage.160323.rds")
}

# score plots
{
  ## calculate score
  pep.gr$score<-calculateScore(pl,"312$") #note that you have to use the $ sign to specify the column with 5th replicates (regex)
  saveRDS(pep.gr, "intermediate_files/pep.gr.score312.160418.rds")
  
  ## plot score
  candidate="Zmynd8"
  pos<-which(prot.gr$symbol==candidate)
  t<-residueScorePlot(prot.gr[pos,],pep.gr,smooth=T, no.negs=T, return.values=T, span=0.001)
  
  a<-as.data.frame(pep.gr[seqnames(pep.gr)==as.character(seqnames(prot.gr)[pos]),])
  a<-a[order(a$score),]
  
  # to make score file that gets inserted in PDB for coloring
  t$label<-paste0(":",t$x,".A")
  RBRscoretable<-cbind(t$label, round(t$y*10)/10)
  write.table(RBRscoretable,sep="\t",quote=F, file="~/analyses/BA57/160404_Polr2a_scoreSmooth_5FLM",row.names = F )
}

# make list of residue scores for all proteins
{
  # rm(list=c("end","start")) # this is because the worker nodes behave strangely otherwise
  # cl<-makeCluster(30)
  # clusterExport(cl, list("subsetByOverlaps","reduce","width","end","setdiff","start")) # must pass non-conventional functions to worker nodes
  # scoreList<-parLapply(cl=cl,prot.gr,residueScore,peps=pep.gr) # this takes ~3'20" with 30 cores
  # stopCluster(cl)
  # names(scoreList)<-seqnames(prot.gr)
  # saveRDS(scoreList,"output/score/E14.depletion.residue_score.160321.r")
  scoreList<-readRDS("output/score/E14.depletion.residue_score.160321.rds")
  allRes<-unlist(scoreList) # long vector with all residue scores
}

# metadomain plots, matrix, heatmaps
{
  ## select domain
  ipr.filt<-mouse_IPR[mouse_IPR$uniprotID %in% seqnames(prot.gr),] # no point keeping domain information for proteins not in the input list
  ## select IPR of interest
  IPR<-ipr.filt[ipr.filt$iprID %in% "IPR007275",]
  
  ## calculate mean score
  {
    for(i in 1:nrow(IPR)){
      IPR$score[i]<-mean(scoreList[[IPR$uniprotID[i]]][IPR$IPR_start[i]:IPR$IPR_end[i]], na.rm=T)
    }
    mean(IPR$score, na.rm=T)
    pval<-t.test(allRes,IPR$score)$p.value
    ## plot combo result
    boxplot(allRes,IPR$score, names=c("Total","Domain(s) of interest"), outline=F, ylab="RBR-ID score", main="Residue-level score") # plots score distribution if needed
    title(sub="dsRBD", line=2)
    legend("topleft",legend=paste0("P-value = ",sprintf("%0.1g", pval)))
  }
  
  ## generate matrix for metadomain
  IPR.f<-IPR[!is.na(IPR$score),] # remove rows with no signal
  aa=25
  bins=50
  combo.matrix<-generateDomainMatrix(IPR.f,scoreList,aa=aa,bins=bins)
  
  ## plot with function
  metaDomainPlot(combo.matrix,aa=25,bins=50,smooth=F)
  
  ## heatmap instead
  {
    require(gplots)
    
    mat1<-combo.matrix
    mat1[is.na(mat1)]<-0
    idx<-order(rowSums(mat1[,26:75]), decreasing=T)
    
    steps=256
    cols<-color.palette(c("white", "black"))(10)
    breaks<-c(min(mat1),2.5,5,7.5,10,12.5,15,17.5,20,22.5,25)
    heatmap.2(mat1[idx,],
              trace="none",Rowv=F,Colv=F,dendrogram="none",labCol=F,labRow=F, density.info="none",
              col=cols, symbreaks=F, breaks=breaks,
              colsep=c(0,aa,aa+bins,aa+bins+aa), sepcolor="black", scale="none",
              xlab=paste0("Domains of interest (binned) ± ",aa," amino acids"), main="IPR027417",
              key.xlab="RBR-ID score", symkey=F, keysize = 1)
  }
}

