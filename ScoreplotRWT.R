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