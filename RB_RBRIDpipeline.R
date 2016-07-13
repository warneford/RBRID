#************************************************************#
#--- BA-59 Full reanalysis of E14 RBR-ID (5 reps 312 nm) ----#
#************************************************************#

### Author: Roberto Bonasio
### Created: 4.26.2016

# UPDATES:

##### INITIALIZATION ####

MACHINE="brainberg"
homed=Sys.getenv("HOME")
wd=paste0(homed,"/analyses/BA59")

####### LIBRARIES ########
require(GenomicRanges)
require(parallel)
require(gplots)
require(biomaRt)
require(Peptides)

####### FUNCTIONS ########

annotatePeptides<-function(index, data, IPR.info)
{
  ## Annotates overlapping IPRs for protein and peptide
  ## Computationally intensive, use parallelization if possible
  ## Takes 1 min for 5,000 proteins in mouse if no parallelization
  ## Takes 1.5 min for ~80,000 proteins with 30 nodes
  peptide<-data[index,]
  hits<-which(IPR.info$uniprotID==peptide$uniprotID)  # finds hits in IPR table
  prot_IPR_hits<-paste(unique(IPR.info$iprID[hits]),collapse="|") # annotates with all IPR hits for the protein (regardless of position)
  
  # the following decides whether the actual peptide overlaps the start, end, or is contained within the various IPRs
  overlap_start<-peptide$start<=IPR.info$IPR_start[hits] & peptide$end>=IPR.info$IPR_start[hits]
  overlap_end<-peptide$start<=IPR.info$IPR_end[hits] & peptide$end>=IPR.info$IPR_end[hits]
  inside<-peptide$start>=IPR.info$IPR_start[hits] & peptide$end<=IPR.info$IPR_end[hits]
  
  pept_IPR_hits<-paste(unique(IPR.info$iprID[hits][overlap_start | overlap_end | inside]),collapse="|") # annotates IPRs that overlap
  return(cbind(prot_IPR_hits,pept_IPR_hits))
}

logFold<-function(values, test.string, control.string, convert.zero.means=T)
{
  ## Calculates the log2 fold change for replicate columns in the matrix value
  ## The right columns are found by matching the test and control string to column names
  ## so make sure the column names are passed with the values as a data frame.
  
  test.columns<-grepl(test.string, names(values))  # find columns for the numerator ("test")
  control.columns<-grepl(control.string, names(values)) # find columns for the denominator ("control")
  test.means<-rowMeans(values[,test.columns], na.rm=T)
  control.means<-rowMeans(values[,control.columns], na.rm=T)
  if (convert.zero.means){
    test.means[test.means==0]<-min(test.means[test.means!=0])/2
    control.means[control.means==0]<-min(control.means[control.means!=0])/2
  }
  log2(test.means/control.means)
}

plotDomainDist<-function(domain.annotation, domains.of.interest, test.order, control.order=NULL, highlights=highlights)
{
  ## This function generates density plots from a vector with domain annotation
  ## an order vector (e.g. probability score, logFold, etc.) and a list of domain_of_interest to feed to grep
  
  domain.string<-paste0(domains.of.interest, collapse="|")  # prepares list of domains as a regex string
  test.dist<-grep(domain.string,domain.annotation[test.order]) # finds position of matching peptides after reordering
  # density plot 
  plot(density(test.dist, from=0, to=length(domain.annotation)),lwd=5,yaxt="n", xlab="",
       main=paste0("Density of ",paste0(domains.of.interest, collapse=", ")),
       sub=paste0(length(test.dist), " of ", length(domain.annotation), " peptides overlap this domain."))
  # plot position of single matching peptides at the bottom (rug)
  rug(test.dist,side=1, lwd=0.2)
  # plot control if present
  if(length(control.order)>0){
    control.dist<-grep(domain.string,domain.annotation[control.order])
    lines(density(control.dist, from=0, to=length(domain.annotation)),lwd=2,lty=2)
  }
  if(length(highlights)>0){
    for(i in highlights){
      abline(v=which(test.order==i), col="red", lwd=2)
      mtext(anno$symbol[i], side=3, at=which(test.order==i), las=3, cex=0.8)
    }
  }
}

plotDomainHist<-function(domain.annotation, domains.of.interest, test.order, control.order=NULL, highlights=NULL)
{
  ## This function generates density plots from a vector with domain annotation
  ## an order vector (e.g. probability score, logFold, etc.) and a list of domain_of_interest to feed to grep
  
  domain.string<-paste0(domains.of.interest, collapse="|")  # prepares list of domains as a regex string
  test.dist<-grep(domain.string,domain.annotation[test.order]) # finds position of matching peptides after reordering
  # histogram plot 
  hist(test.dist, breaks=50, freq=F, xlim=c(0,length(domain.annotation)), xlab="", 
       main=paste0("Density of ",paste0(domains.of.interest, collapse=", ")), 
       sub=paste0(length(test.dist), " of ", length(domain.annotation), " peptides overlap this domain."))
  # overlap density
  lines(density(test.dist), lwd=5)
  # plot position of single matching peptides at the bottom (rug)
  rug(test.dist,side=1, lwd=0.1)
  # plot control density if present
  if(length(control.order)>0){
    control.dist<-grep(domain.string,domain.annotation[control.order])
    lines(density(control.dist, from=0, to=length(domain.annotation)),lwd=2,lty=2)
  }
  if(length(highlights)>0){
    for(i in highlights){
      abline(v=which(test.order==i), col="red", lwd=2)
      mtext(anno$symbol[i], side=3, at=which(test.order==i), las=3, cex=0.8)
    }
  }
}

logFoldPlots<-function(values, domain.annotation, domains.of.interest=c("IPR000504"), hist=F, highlights=NULL, convert.zero.means=T)
{
  ## This function generates four standard plots for each comparison
  ## That is noUV, 312, and 365 ± 4SU as well as ± UV 254
  ## Different domains can be passed in as a character vector, default is RRM domain (IPR000504)
  ## The logical variable 'hist' directs the function toward the histogram plot or the density plot
  
  par(mfrow=c(2,2))
  if(hist){
    order<-order(logFold(values,"UV312.4SU","UV312.no4SU", convert.zero.means = convert.zero.means))
    control.order<-order(logFold(values,"UV312.4SU.b[123]","UV312.no4SU.b[123]"))
    plotDomainHist(domain.annotation, domains.of.interest, order, control.order, highlights=highlights)
    legend("topright",lwd=c(5,2), lty=c(1,2),legend=c("UV312 ± 4SU (5 reps)","UV312 ± 4SU (3 reps)"), cex=0.75)
    
    order<-order(logFold(values,"noUV.4SU","noUV.no4SU", convert.zero.means = convert.zero.means))
    plotDomainHist(domain.annotation, domains.of.interest, order, sample(order), highlights=highlights)
    legend("topright",lwd=5,"noUV ± 4SU", cex=0.75)
    
    order<-order(logFold(values,"UV365.4SU","UV365.no4SU", convert.zero.means = convert.zero.means))
    plotDomainHist(domain.annotation, domains.of.interest, order, sample(order), highlights=highlights)
    legend("topright",lwd=5,"UV365 ± 4SU", cex=0.75)
    
    order<-order(logFold(values,"UV254.no4SU","noUV.no4SU", convert.zero.means = convert.zero.means))
    plotDomainHist(domain.annotation, domains.of.interest, order, sample(order), highlights=highlights)
    legend("topright",lwd=5,"± UV254", cex=0.75)
  } 
  else {
    order<-order(logFold(values,"UV312.4SU","UV312.no4SU"))
    control.order<-order(logFold(values,"UV312.4SU.b[123]","UV312.no4SU.b[123]"))
    plotDomainDist(domain.annotation, domains.of.interest, order, control.order, highlights=highlights)
    legend("topright",lwd=c(5,2), lty=c(1,2),legend=c("UV312 ± 4SU (5 reps)","UV312 ± 4SU (3 reps)"), cex=0.75)
    
    order<-order(logFold(values,"noUV.4SU","noUV.no4SU"))
    plotDomainDist(domain.annotation, domains.of.interest, order, sample(order), highlights=highlights)
    legend("topright",lwd=5,"noUV ± 4SU", cex=0.75)
    
    order<-order(logFold(values,"UV365.4SU","UV365.no4SU"))
    plotDomainDist(domain.annotation, domains.of.interest, order, sample(order), highlights=highlights)
    legend("topright",lwd=5,"UV365 ± 4SU", cex=0.75)
    
    order<-order(logFold(values,"UV254.no4SU","noUV.no4SU"))
    plotDomainDist(domain.annotation, domains.of.interest, order, sample(order), highlights=highlights)
    legend("topright",lwd=5,"± UV254", cex=0.75)
  }
  par(mfrow=c(1,1))
}

pval.apply.varequal<-function(incoming, test.columns, control.columns, alternative="two.sided", paired=F)
{
  ## Function to calculate pvalue in the apply call
  
  x<-incoming[test.columns]
  y<-incoming[control.columns]
  if((sum(!is.na(x))<2 & sum(!is.na(y))<2) | sum(!is.na(x))==0 | sum(!is.na(y))==0)  pval=1
  else pval<-t.test(incoming[test.columns],incoming[control.columns], alternative=alternative, paired=paired, var.equal=T)$p.value
}

calculatePval<-function(values, test.string, control.string, paired=F, alternative="two.sided", var.equal=F)
{
  ## Calculates the t-test pvalue for replicate columns in the matrix value
  ## The right columns are found by matching the test and control string to column names
  ## so make sure the column names are passed in the same order as the column in the values matrix.
  ## If a paired test is desired all vectors must have same number of test and control values (i.e. no NA's).
  
  test.columns<-grepl(test.string,names(values))
  control.columns<-grepl(control.string,names(values))
  if(var.equal){
    p<-parApply(cl=cl,values, 1, pval.apply.varequal, test.columns=test.columns, control.columns=control.columns, alternative=alternative, paired=paired)
  }
  else {
    p<-parApply(cl=cl,values, 1, pval.apply, test.columns=test.columns, control.columns=control.columns, alternative=alternative, paired=paired)
  }
  p[is.na(p)]<-1 # replace all NA's with P=1
  p
}

significance.matrix<-function(values, paired=F, alternative="two.sided", var.equal=F, convert.zero.means=T)
{
  ## Calculates pvalues and logFold for all relevant comparisons
  
  pval<-calculatePval(values, "UV312.4SU", "UV312.no4SU", paired=paired, alternative=alternative, var.equal=var.equal)
  logFold<-logFold(values, "UV312.4SU", "UV312.no4SU", convert.zero.means = convert.zero.means)
  pval.ctrl<-calculatePval(values, "noUV.4SU", "noUV.no4SU", paired=paired, alternative=alternative, var.equal=var.equal)
  logFold.ctrl<-logFold(values, "noUV.4SU", "noUV.no4SU", convert.zero.means = convert.zero.means)
  
  as.data.frame(cbind(pval,logFold,pval.ctrl,logFold.ctrl))
}

calculateCoverage<-function(range,peps)
{
  ## Calculates protein-level coverage using genomeRanges
  
  hits<-subsetByOverlaps(peps,range)
  covered.regions<-reduce(hits)
  sum(width(covered.regions))/end(range)
}

###---- MAIN SCRIPT ----###

####### PEPTIDE ANALYSES ########

# Load and convert 1.2.PC50 (E14 50 runs, 5 replicates of 312)
{
  E14.50<-read.csv("input/1.2.1.PC50.E14.depletion.50runs.csv", stringsAsFactors=FALSE) # raw intensities for E14 depletion, 50 runs
  E14.50<-E14.50[E14.50$Proteins!="",] # clean up in case there are weird rows (happens sometimes)
  E14.50$pepID<-paste0("PEP",seq(1:nrow(E14.50))) # generate unique pepID (note: they are not same across experiments)
  
  ## reorganize annotation info
  E14.50<-E14.50[,c(59,1,3,2,5:8,4,9:58)] 
  names(E14.50)[1:9]<-c("pepID", "uniprotID", "symbol", "name", "sequence", "start", "end", "unique", "contaminant")
  
  ## reorganize data
  E14.50<-E14.50[,c(1:9,10:15,16:21,22:27,34:39,46:51,28:33,40:45,52:59)] 
  ### rename data columns
  sampleNames<-rep(c("noUV.no4SU", "UV254.no4SU", "UV312.no4SU", "UV365.no4SU", "noUV.4SU", "UV312.4SU", "UV365.4SU"), each=6)
  repNames<-rep(c("b1t1","b1t2","b2t1","b2t2","b3t1","b3t2"),7)
  runNames<-paste(sampleNames,repNames, sep=".")
  runNames<-c(runNames, "UV312.no4SU.b4t1", "UV312.no4SU.b4t2", "UV312.no4SU.b5t1", "UV312.no4SU.b5t2",  "UV312.4SU.b4t1", "UV312.4SU.b4t2", "UV312.4SU.b5t1", "UV312.4SU.b5t2")
  names(E14.50)[10:59]<-runNames
  ### move additional 312 replicates next to old ones
  E14.50<-E14.50[,c(1:27,52:55,28:45,56:59,46:51)] 
  
  ## change contaminant and unique columns to logical vectors
  ### contaminant
  E14.50$contaminant[E14.50$contaminant=="+"]<-T
  E14.50$contaminant[E14.50$contaminant!="TRUE"]<-F
  E14.50$contaminant<-as.logical(E14.50$contaminant)
  ### unique
  E14.50$unique[E14.50$unique=="yes"]<-T
  E14.50$unique[E14.50$unique!="TRUE"]<-F
  E14.50$unique<-as.logical(E14.50$unique)
  
  ## SAVE
  E14<-E14.50
  rm(E14.50)
  saveRDS(E14,file="intermediate/1.E14.50runs.raw.rds")
}
E14<-readRDS("intermediate/1.E14.50runs.raw.rds")

# Annotate peptide & protein overlap with Interpro domains
{
  ## Load and change column names for protein2ipr table (positions and ID's for all Interpro call on all mouse Uniprot)
  mouse_IPR<-read.delim("~/genomes/Mmus/uniprot/150906.protein2ipr.mouse", header=FALSE, stringsAsFactors=FALSE)
  names(mouse_IPR)<-c("uniprotID","iprID","IPR_name","IPR_signature","IPR_start","IPR_end")
  
  ## Use parallel sapply to find overlaps of peptides with annotated domains
  cl<-makeCluster(30)
  IPR.annotations<-t(parSapply(cl=cl,1:nrow(E14),annotatePeptides,data=E14,IPR.info=mouse_IPR))
  stopCluster(cl)
  
  ## Copy results to main data frame
  E14$prot_IPR_hits<-IPR.annotations[,1]
  E14$pept_IPR_hits<-IPR.annotations[,2]
  rm(IPR.annotations)
  
  ## Move new annotations to annotation section
  E14<-E14[,c(1:9,60:61,10:59)]
  
  ## SAVE
  saveRDS(E14,file="intermediate/2.E14.50runs.raw.IPR.rds")
}
E14<-readRDS("intermediate/2.E14.50runs.raw.IPR.rds")

# Calculate logFold, P-value, and score
{
  ## Load and split table
  anno<-E14[,1:11]
  raw<-E14[,12:61]
  row.names(raw)<-anno$pepID
  
  ## Some constants
  RRM="IPR000504"
  KH="IPR004087"
  dsRBD="IPR014720"
  RBDs<-c(RRM,KH,dsRBD)
  
  ## Find positions of candidate peptides
  L1td1<-which(anno$sequence=="FISDIPYLKDLLNNIH")
  Tet2<-which(anno$sequence=="HCLALWEAK")
  Nanog<-which(anno$sequence=="LSSPEADKGPEEEENKVLAR")
  Pcgf1<-which(anno$sequence=="MDPLRNEEEVR")
  Pcgf2<-which(anno$sequence=="SDKTLQDIVYK")
  Pou5f1<-which(anno$sequence=="ELEQFAK")
  Hnrnpc<-which(anno$sequence=="MIAGQVLDINLAAEPK")
  
  ## Presets
  var.equal=T
  highlights=c(L1td1,Tet2, Nanog, Pcgf1, Pcgf2, Pou5f1)
  highlights=c(Hnrnpc)
  
  ## Total normalization
  tot<-raw
  ### divide each column by its sum (without counting the 0's)
  for(i in 1:ncol(tot)) {
    tot[,i]<-tot[,i]/sum(tot[,i])
  }
  pep<-tot
  
  ## Draw a plot based on normalized counts
  logFoldPlots(pep, anno$pept_IPR_hits, hist=F, highlights=highlights, convert.zero.means = T)
  
  ## Calculate paired t-test P-values
  cl<-makeCluster(30)
  sigMat<-significance.matrix(pep, paired=T, var.equal=var.equal, convert.zero.means=T)
  stopCluster(cl)
  anno$logFold<-sigMat$logFold
  anno$pval<-sigMat$pval
  anno$score<-(-sigMat$logFold)*(log10(sigMat$pval))^2
  anno$score.ctrl<-(-sigMat$logFold.ctrl)*(log10(sigMat$pval.ctrl))^2
  
  ## Check results by plotting all RRM domain and some candidates
  mask<-grepl("IPR000504",anno$pept_IPR_hits)
  score.test.order<-order(anno$score, decreasing=T)
  plotDomainHist(anno$pept_IPR_hits, "IPR000504",score.test.order,highlights=highlights)
  
  ## SAVE
  saveRDS(anno,file="intermediate/3.1.E14.scores.rds")
  saveRDS(raw,file="intermediate/3.2.E14.raw.rds")
  rm(list=c("anno","raw","pep","tot","sigMat","E14"))
}
E14<-readRDS("intermediate/3.1.E14.scores.rds")

####### PROTEIN ANALYSES ########

# Clean up and extend annotation
require(biomaRt)
{
  ## load online annotations
  mENS84<-useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host="mar2016.archive.ensembl.org") # set up database using ENSEMBL v84
  
  mouse.ids<-getBM(c("ensembl_gene_id","external_gene_name","chromosome_name"),mart=mENS84) # get all genes and names in ENS84
  mouse.ids<-mouse.ids[mouse.ids$chromosome_name %in% unique(mouse.ids$chromosome_name)[1:22],]
  
  mouse.ids.uniprot<-getBM(c("ensembl_transcript_id","ensembl_peptide_id","uniprot_sptrembl"),mart=mENS84)
  mouse.ids.allprot<-getBM(c("ensembl_gene_id","external_gene_name","ensembl_transcript_id"),mart=mENS84)
  
  ## create unique list of uniprotIDs from peptide list
  anno<-E14
  prot<-anno[!anno$contaminant,c("uniprotID","symbol","name")] # remove contaminants
  prot<-prot[order(prot$uniprotID, decreasing=T),] # sort by uniprot ID
  prot<-prot[!duplicated(prot$uniprotID),] # make protein list (one row = one protein)
  
  ## Update uniprotIDs using newly downloaded file
  mUni<-read.delim("annotations/mouse_uniprot_GO.160406.tsv", stringsAsFactors=F)
  names(mUni)<-c("uniprotID","updated_uni_symbol","synonims","ensembl_transcript_id","subcellular_location","GO_BP","GO_MF","GO_CC")
  prot<-merge(prot,mUni,by="uniprotID", all.x=T)
  prot$updated_uni_symbol[is.na(prot$updated_uni_symbol)]<-prot$symbol[is.na(prot$updated_uni_symbol)] # for the 12 proteins that do not much uniprotID
  notSame<-prot$symbol != prot$updated_uni_symbol # mark proteins with changed symbol
  oldGood<-prot$symbol %in% mouse.ids$external_gene_name # symbols from Simone's table match ENS84
  newGood<-prot$updated_uni_symbol %in% mouse.ids$external_gene_name # symbols from new uniprot table match ENS84
  prot$consensus_symbol<-""
  prot$consensus_symbol[!notSame]<-prot$symbol[!notSame] # copy those that are the same in two versions
  prot$consensus_symbol[notSame & oldGood]<-prot$symbol[notSame & oldGood] # first copy those that were good already
  prot$consensus_symbol[notSame & newGood]<-prot$updated_uni_symbol[notSame & newGood] # then copy the newer ones, in the 4 cases when they both match ENS84 the newer dominates
  prot$consensus_symbol[notSame & !oldGood & !newGood & prot$updated_uni_symbol!=""]<-prot$updated_uni_symbol[notSame & !oldGood & !newGood & prot$updated_uni_symbol!=""]
  prot$consensus_symbol[notSame & !oldGood & !newGood & prot$updated_uni_symbol==""]<-prot$symbol[notSame & !oldGood & !newGood & prot$updated_uni_symbol==""]
  prot$symbol<-prot$consensus_symbol # consolidate 
  prot<-prot[,c(1:3,5:10)] # remove temp columns
  
  ## Add ENSG IDs using symbols alone (three genes have multiple ENSG ID and give warnings)
  prot$ensembl_gene_id.symb<-"" # empty column
  for(i in 1:nrow(prot)){
    if(!prot$symbol[i] %in% mouse.ids$external_gene_name) {
      prot$ensembl_gene_id.symb[i]<-""
    } else {
      prot$ensembl_gene_id.symb[i]<-mouse.ids$ensembl_gene_id[mouse.ids$external_gene_name==prot$symbol[i]]
    }
  }
  
  ## Add ENSG IDs matching the first transcript associated with uniprotID
  prot$ensembl_transcript_id[is.na(prot$ensembl_transcript_id)]<-""
  prot$ensembl_gene_id.enst<-""
  for(i in 1:nrow(prot)){
    if(prot$ensembl_transcript_id[i]=="") next
    trans<-unlist(strsplit(prot$ensembl_transcript_id[i],";"))
    trans.there<-trans %in% mouse.ids.allprot$ensembl_transcript_id
    if(sum(trans.there)==0) next
    this.trans<-trans[which(trans.there)[1]]
    prot$ensembl_gene_id.enst[i]<-mouse.ids.allprot$ensembl_gene_id[mouse.ids.allprot$ensembl_transcript_id==this.trans]
  }
  
  ## Consolidate ENSG
  prot$ensembl_gene_id<-prot$ensembl_gene_id.symb
  prot$ensembl_gene_id[prot$ensembl_gene_id.symb==""]<-prot$ensembl_gene_id.enst[prot$ensembl_gene_id.symb==""]
  E14.prot<-prot[,c(1,12,2,4,3)]
  saveRDS(E14.prot,"intermediate/4.E14.proteins.rds")
}
E14.prot<-readRDS("intermediate/4.E14.proteins.rds")

# Calculate coverage per protein
{
  ## Load and change column names for protein2ipr table & protein lengths
  mouse_IPR<-read.delim("~/genomes/Mmus/uniprot/150906.protein2ipr.mouse", header=FALSE, stringsAsFactors=FALSE)
  names(mouse_IPR)<-c("uniprotID","iprID","IPR_name","IPR_signature","IPR_start","IPR_end")
  lengths<-read.delim(paste0(homed,"/genomes/Mmus/uniprot/160317.uniprot.mouse.lengths"))
  
  ## mark proteins with annotation problems (changes between december 15 and march 16)
  problem.proteins<-c("D3YU33") 
  
  ## create unique list of uniprotIDs from peptide list
  anno<-E14.prot
  prot.anno<-anno[!prot.anno$uniprotID %in% problem.proteins,c("uniprotID","symbol","name")] # remove contaminants
  prot.anno<-merge(prot.anno, lengths, by="uniprotID") # add lengths (this drops 8 proteins due to changes in annotation)
  
  ## convert to grange
  prot.gr<-GRanges(seqnames=Rle(prot.anno$uniprotID), ranges=IRanges(start=1, end=prot.anno$length))
  prot.gr$symbol<-prot.anno$symbol
  
  ## prepare peptide grange
  pep.gr<-GRanges(seqnames=Rle(E14$uniprotID),ranges=IRanges(start=E14$start,end=E14$end))
  pep.gr$pepID<-E14$pepID; pep.gr$unique<-E14$unique; pep.gr$contaminant<-E14$contaminant; pep.gr$sequence<-E14$sequence
  
  ## calculate coverage per protein
  cl<-makeCluster(30)
  clusterExport(cl, list("subsetByOverlaps","reduce","width","end")) # must pass non-conventional functions to worker nodes
  coverage<-parSapply(cl=cl,prot.gr,calculateCoverage, peps=pep.gr) # this takes about 1'15" with 30 nodes
  stopCluster(cl)
  prot.gr$coverage<-coverage*100
  prot<-as.data.frame(prot.gr)
  prot<-prot[,c(1,4,7)]
  names(prot)<-c("uniprotID","length","coverage")
  E14.prot<-merge(E14.prot,prot,by="uniprotID", all.x=T)
  E14.prot$synonims[is.na(E14.prot$synonims)]<-""
  E14.prot<-E14.prot[,c(1:4,6:7,5)]
  saveRDS(E14.prot,"intermediate/5.E14.proteins.rds")
}
E14.prot<-readRDS("intermediate/5.E14.proteins.rds")

##### TEMP LIST GENERATION #####

a<-merge(E14,E14.prot,by="uniprotID", all.x=T)
a<-a[,c(1,16,17,3,19,20,10,9,4,5:8,11:15)]
names(a)[3:4]<-c("symbol","symbol_old")
names(a)[9]<-"name"
a$primary<-a$logFold<log2(0.95) & a$pval < 0.05
a$extended<-a$logFold<log2(0.95) & a$pval < 0.1
a<-a[order(a$score, decreasing=T),]
row.names(a)<-1:nrow(a)

saveRDS(a,"output/160517.RBRID.E14.list.rds")




























