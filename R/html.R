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
        "<body bgcolor= >", "<H1 ALIGN=CENTER> MACAT: Microarray Chromosomal Abberation Tool</H1>",
        file = outfile, sep = "\n")
    if (!missing(mytitle)) 
        cat("<CENTER><H1 ALIGN=\"CENTER\">", mytitle, " </H1></CENTER>\n", file = outfile,
            sep = "\n")


    cat("<br>", "<H3>SlidingAverage Plot: <br >Yellow dotted regions seem to be interesting according to the statistic<\H3>", file=outfile)
    
    #---------------------------------------------
    # here the pictures
    
    cat("<img src=\"",SLIDINGpic,"\" width=\"\" height=\"\" alt=\"\" border=\"0\">\n", sep="",file = outfile)

    # adjust pic parameter for chromosom
    st="\"margin-left:70px; margin-top:0px\""
    w = "\"796\""
   
    
    cat("<p style=",st,"><img src=\"",CHROMpic,"\" width=",w,"height=\"80\" alt=\"\" border=\"0\"><\p>" ,sep="" ,file = outfile)
    cat("<br><br>", file=outfile)
    cat("<H3> List of interesting genes:<br>Yellow marked regions in the plot<H3>", file=outfile)


    #---------------------------------------------

    nrows <- length(genelist)
    if (nrows!=0) {
      cat("<TABLE BORDER=4>", file = outfile, sep = "\n")
      if (!missing(table.head)) {
        headout <- paste("<TH>", table.head, "</TH>")
        cat("<TR>", headout, "</TR>", file = outfile, sep = "\n")
      }
      
      # annotate function
      rows <- getTDRows(genelist, "ll")
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
}


####################################################################################
#
# INTERNAL Function finds all genes for given positions on chromosome.
# It builds a list of usefull things and assign it to myhtml for htmloutput.
# SLIDINGpic: Has to be 'png' or 'jpeg'
#
####################################################################################
get.html <- function(MACATevalScoringOBJ, SLIDINGpic, HTMLfilename, mytitle){
   # SLIDINGpic: PNG of the scoring slidingAverage
   # HTMLfilename: output html filename
   # mytitle: Title of HTMLpage

   res=get.results(MACATevalScoringOBJ)

   # path for pics
   #chrompic=paste("http://macat.sourceforge.net/pics/ch",MACATevalScoringOBJ$chromosome,".png", sep="")
   chrompic=paste("http://compdiag.molgen.mpg.de/software/macatchrompics/ch",MACATevalScoringOBJ$chromosome,".png", sep="")
   
   if (is.null(res)){
     head=c()
     myhtml(c(), MACATevalScoringOBJ$chromosome, SLIDINGpic, chrompic , HTMLfilename,
            mytitle, table.head=head)
   }
   else {
     head=c("LocusID", "ProbeID", "Cytoband", "GeneSYM", "GeneDescription", "p-Value")
     mylocusid <- unlist(res$locusid)
     myhtml(mylocusid, res$chromosome, SLIDINGpic, chrompic , HTMLfilename,
            mytitle, list(res$probeID, res$cytoband, res$geneSYM, res$genedescription,
                        res$pvalue),table.head=head)
   }
   path=getwd()
   html = paste(path, "/", HTMLfilename, sep="" )
   browseURL(html)
} # get.html


####################################################################################
# Functiopn returns a list with important information about interesting sides
# on chromosome
####################################################################################
get.results <- function(MACATevalScoringOBJ){

   attach(MACATevalScoringOBJ)
   step.width = steps[2] - steps[1]
   #----------------------------------------------
   # find all genes in significant regions
   # vec of significant Genes
   sig=rep(FALSE,length(original.geneid))
   # vec of significant sliding.values
   issig = ((sliding.value>upper.permuted.border)|(sliding.value<lower.permuted.border))
   
   interpolate <- function(i, loc, values){
    
     y = values[i]+((values[i+1]-values[i])/(steps[i+1]-steps[i])) * (loc - steps[i])
     return(y)
   }
  
   # find significant genes
   gene = 0
   for (loc in original.loc) {
     gene = gene + 1     
     i = floor((loc-steps[1])/ step.width) + 1

     if (i == length(steps)) {
         sig[gene] = issig[i]   
     }
     else {
       if (issig[i] == FALSE & issig[i+1] == FALSE) {
         next 
       }
       else {
         sliding = interpolate(i, loc, sliding.value)
         lower = interpolate(i, loc, lower.permuted.border)
         upper = interpolate(i, loc, upper.permuted.border)
         sig[gene] = ((sliding>upper)|(sliding<lower))
       }
     }
   }
  
   genes=original.geneid[sig]
   if (length(genes)==0){
     return(NULL)
   }
   pvalues= original.pvalue[sig]
   #-------------------------------------------

   # look for interesting things for genes on chip
   require(chip, character.only=TRUE)
   e1=paste(chip, "GENENAME", sep="")
   genedescription=mget(genes, env=eval(as.symbol(e1)))
   e2=paste(chip, "LOCUSID", sep="") 
   locusids=mget(genes, env=eval(as.symbol(e2)))
   e3=paste(chip, "SYMBOL", sep="")
   symbol=mget(genes, env=eval(as.symbol(e3)))
   e4=paste(chip, "MAP", sep="")
   cytoband=mget(genes, env=eval(as.symbol(e4)))

   result <-  list(probeID=genes,cytoband=cytoband,geneSYM=symbol,
                  pvalue=pvalues, locusid=locusids, genedescription=genedescription,
                  chromosome=chromosome, class=class)

   detach(MACATevalScoringOBJ)
   return(result)
} # get.results
