####################################################################################
#
# INTERNAL Function writes html source to a file. It builds a linkage list for
# locuslinkIDs. 
# Function based on annotates ll.htmlpage.
# Several modifications e.g. to bind in pictures etc.
#
####################################################################################
myhtml <- function (genelist, chrom, SLIDINGpic, CHROMpic, filename, mytitle, othernames,
                    table.head) {
    # Function based on annotates ll.htmlpage function

    require(annotate)
    outfile <- file(filename, "w")
    cat("<html>", "<head>", "<TITLE>Result Page for Chromosome</TITLE>", "</head>",
        "<body bgcolor= >", "<H2 ALIGN=CENTER> MACAT: MicroArray Chromosome Analysis Tool</H2>",
        file = outfile, sep = "\n")
    if (!missing(mytitle)) 
        cat("<CENTER><H2 ALIGN=\"CENTER\">", mytitle, " </H2></CENTER>\n", file = outfile,
            sep = "\n")


    cat("<br>", "<H3>Result of Kernel Smoothing</H3><p>Yellow dotted regions are considered significant.</p>", file=outfile)
    
    #---------------------------------------------
    # here the pictures
    
    cat("<img src=\"",SLIDINGpic,"\" alt=\"\" border=\"0\">\n", sep="",file = outfile)

    # adjust pic parameter for chromosom
    st="\"margin-left:70px; margin-top:0px\""
    w = "\"796\""
   
    
    cat("<p style=",st,"><img src=\"",CHROMpic,"\" width=",w,"height=\"80\" alt=\"\" border=\"0\"></p>" ,sep="" ,file = outfile)
    cat("<br><br>", file=outfile)
    cat("<H3>Genes within Significant Regions</H3>", file=outfile)


    #---------------------------------------------

    nrows <- length(genelist)
    if (nrows!=0) {
      cat("<TABLE BORDER=4>", file = outfile, sep = "\n")
      if (!missing(table.head)) {
        headout <- paste("<TH>", table.head, "</TH>")
        cat("<TR>", headout, "</TR>", file = outfile, sep = "\n")
      }
      
      # annotate function
      rows <- annotate:::getTDRows(genelist, "en")
      if (!missing(othernames)) {
        if (is.list(othernames)) {
          others <- ""
          for (nm in othernames) others <- paste(others, "<TD>", nm, "</TD>", sep = "")
        }
        else others <- paste("<TD>", othernames, "</TD>", sep = "")
        rows <- paste(rows, others)
      }
    
    
      for (i in 1:nrows) cat("<TR>", rows[i], "</TR>", file = outfile,sep = "\n")
      cat("</TABLE>", file = outfile)
    }
    else {
      cat("<H3> There are no genes in significant regions.</H3>", "\n", file=outfile)
    } 
    cat("</body>", "</html>", sep = "\n", file = outfile)
    close(outfile)
} #myhtml





####################################################################################
#
# INTERNAL Function finds all genes for given positions on chromosome.
# It builds a list of usefull things and assign it to myhtml for htmloutput.
# SLIDINGpic: Has to be 'png' or 'jpeg'
#
####################################################################################
getHtml <- function(MACATevalScoringOBJ, SLIDINGpic, HTMLfilename, mytitle){
   # SLIDINGpic: PNG of the scoring slidingAverage
   # HTMLfilename: output html filename
   # mytitle: Title of HTMLpage

   res=getResults(MACATevalScoringOBJ)

   # path for pics
   chrompic=paste("http://www.ebi.ac.uk/huber-srv/data/macatchrompics/ch",MACATevalScoringOBJ$chromosome,".png", sep="")
   
   if (is.null(res)){
     head=c()
     myhtml(c(), MACATevalScoringOBJ$chromosome, SLIDINGpic, chrompic , HTMLfilename,
            mytitle, table.head=head)
   }
   else {
     head=c("Entrezgene ID", "ProbeSet ID", "Cytoband", "Gene Symbol", "Gene Description", "Score", "p-Value")
     mylocusid <- unlist(res$locusid)
     myhtml(mylocusid, res$chromosome, SLIDINGpic, chrompic , HTMLfilename,
            mytitle, list(res$probeID, res$cytoband, res$geneSYM, res$genedescription, res$probeScore,
                        res$pvalue),table.head=head)
   }
   path=getwd()
   html = paste(path, "/", HTMLfilename, sep="" )
   slidepic <-  paste(path, "/", SLIDINGpic, sep="" )
   system(paste("chmod 0744",slidepic))
   system(paste("chmod 0744",html))
   browseURL(html)
} # getHtml

##########################################################
### auxiliary functions to work with newer .db packages:
############################################################
getFromDb <- function(ids, chip, envName){
  thisEnv <- paste(gsub("\\.db$","", chip), envName, sep="")
  as.character(get(thisEnv)[ids])
}

####################################################################################
# Functiopn returns a list with important information about interesting sides
# on chromosome
####################################################################################
getResults <- function(MACATevalScoringOBJ){

   stopifnot(inherits(MACATevalScoringOBJ,"MACATevalScoring")) # check class of object
   libname <- gsub("\\.db\\.db$", ".db",
           	paste(MACATevalScoringOBJ$chip,"db",sep="."))
   require(libname, character.only=TRUE)
      
   step.width <- MACATevalScoringOBJ$steps[2] - MACATevalScoringOBJ$steps[1]
   #----------------------------------------------
   # find all genes in significant regions
   # vec of significant Genes
   sig <- vector("logical", length(MACATevalScoringOBJ$original.geneid))
   # vec of significant sliding.values
   issig <- with(MACATevalScoringOBJ, (sliding.value>upper.permuted.border)|(sliding.value<lower.permuted.border))
   
   interpolate <- function(i, loc, values){
     y = values[i]+((values[i+1]-values[i])/(MACATevalScoringOBJ$steps[i+1]-MACATevalScoringOBJ$steps[i])) * (loc - MACATevalScoringOBJ$steps[i])
     return(y)
   }
  
   # find significant genes
   gene = 0
   for (loc in MACATevalScoringOBJ$original.loc) {
     gene = gene + 1     
     i = floor((loc-MACATevalScoringOBJ$steps[1])/ step.width) + 1

     if (i == length(MACATevalScoringOBJ$steps)) {
         sig[gene] = issig[i]   
     }
     else {
       if (issig[i] == FALSE & issig[i+1] == FALSE) {
         next 
       }
       else {
         sliding = interpolate(i, loc, MACATevalScoringOBJ$sliding.value)
         lower = interpolate(i, loc, MACATevalScoringOBJ$lower.permuted.border)
         upper = interpolate(i, loc, MACATevalScoringOBJ$upper.permuted.border)
         sig[gene] = ((sliding>upper)|(sliding<lower))
       }
     }
   }
  
   genes  <- MACATevalScoringOBJ$original.geneid[sig]
   if (length(genes)==0){
     return(NULL)
   }
   # now take care of genes annotated to more than one chromosomal position:
   #  if they are annotated to the same chromosomal band, discard copies
   areDuplicated <- duplicated(genes)
   genesU <- genes[!areDuplicated]
   #-------------------------------------------

   # collapse gene annotations for the same gene   
   contractMultiple <- function(charVec) {
     return(paste(charVec, collapse="; "))
   }
   genedescription <- getFromDb(genesU, MACATevalScoringOBJ$chip, "GENENAME")
   locusids <- getFromDb(genesU, MACATevalScoringOBJ$chip, "ENTREZID")
   symbol   <- getFromDb(genesU, MACATevalScoringOBJ$chip, "SYMBOL")
   cytoband <- getFromDb(genesU, MACATevalScoringOBJ$chip, "MAP")

   contract <- function(listX){
     return(sapply(listX, contractMultiple))
   }

   genedescription <- contract(genedescription)
   locusids        <- contract(locusids)
   symbol          <- contract(symbol)
   cytoband        <- contract(cytoband)

   # if one gene is annotated twice to the same chromosomal band remove the additional copies:
   pvalues <- MACATevalScoringOBJ$original.pvalue[sig][!areDuplicated]
   geneLoc <- MACATevalScoringOBJ$original.loc[sig][!areDuplicated]
   genesUScore <- round(MACATevalScoringOBJ$original.score[sig][!areDuplicated], digits=2)

   result <-  list(probeID=genesU, cytoband=cytoband, geneSYM=symbol,
                   pvalue=pvalues, locusid=locusids,
                   genedescription=genedescription,
                   probeScore= genesUScore, chromosome=MACATevalScoringOBJ$chromosome, class=class)

   #detach(MACATevalScoringOBJ)
   return(result)
} # getResults

