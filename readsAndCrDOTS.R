
#hella hella beta 

#you may need to get some packages from the bioconductor first if you no haz 
#source("https://bioconductor.org/biocLite.R")

library(GenomicRanges)
library(GenomicAlignments)
library(tidyverse)

getwd()
#setwd('Desktop/LCANAL/16s/')
list.files()

    ##These are the regions:
    data="180607.nL_hg38.bam" 
    target.bed="tertTSS_15kbPLusMinus.bed"
    
    #how file name is saved:
    nombre='180607.chr510Mb'
    
##########################    
    
    reads=readGAlignments(data, use.names=T)
    target=read_tsv(target.bed,col_names=c("chr","start","stop")) 

    target.gr=GRanges(target)
    reads.gr=GRanges(reads)

#subsetting reads to pos + neg stranded data overlapping the target   
    pos.reads <-  subset(reads, strand == '+' )
    neg.reads <-  subset(reads, strand == '-' )
    pos.reads.gr=GRanges(pos.reads)
    neg.reads.gr=GRanges(neg.reads)

    pos.targ=overlapsAny(pos.reads.gr, target.gr)
    neg.targ=overlapsAny(neg.reads.gr, target.gr)
    ontarg.pos=pos.reads[pos.targ]
    ontarg.neg=neg.reads[neg.targ] 

# make read data into tibbles 
    ontarg.tb.pos = tibble(start=start(ontarg.pos),end=end(ontarg.pos))
    # generate the y values for rectangular plot
    ontarg.tb.pos=mutate(ontarg.tb.pos,   y1=seq(length(ontarg.pos)), y0=y1-1)

#same for neg strand    
    ontarg.tb.neg = tibble(start=start(ontarg.neg),end=end(ontarg.neg))
    ontarg.tb.neg=mutate(ontarg.tb.neg,   y1=seq(length(ontarg.neg)), y0=y1-1)   
    #invert values for neg strand 
    ontarg.tb.neg=mutate(ontarg.tb.neg,   y1=-(y1), y0=-(y0)) 

#make plot object for rectangular read plot    
    g.rect=ggplot(data=ontarg.tb.pos)+
   # mapping=aes(xmin=start,xmax=end,ymin=y0,ymax=y1))
      geom_rect(mapping=aes(xmin=start,xmax=end,ymin=y0,ymax=y1), color='red')+ theme_bw()+
      geom_rect(mapping=aes(xmin=start,xmax=end,ymin=y0,ymax=y1), data=ontarg.tb.neg, color='blue')+
      labs(x="Coordinate",y="Read")+xlim(target$start,target$stop)

    #hardcoding in sites of crRNAs for rn (TERT)
    guide.df=data.frame(PAMs=c(1296412,1288556) , y0=c(0,0))
    g.rect <- g.rect +  geom_point( data=guide.df,  mapping=aes(x=PAMs, y=y0), size=0.1, color="purple"          )
    
    pdf(paste0(nombre, ".target_plot.pdf"))
    print(g.rect) 
#    print(g.guide)
#    print(g.cov) 
    dev.off()
