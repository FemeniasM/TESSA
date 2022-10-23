#!/usr/bin/Rscript

#outdir=as.character("/home/martin/work/misMS/TEESA/TESSA_20220914_out")
#summarize=as.character('repName')
#conditions=as.character('none')

outdir=as.character(commandArgs(TRUE)[1])
summarize=as.character(commandArgs(TRUE)[2])
conditions=as.character(commandArgs(TRUE)[3])

ref_tessa=paste0(outdir,"/references.csv")

split.vec <- function(vec, pattern, fragPos=1){
  lista <- stringr::str_split(vec, pattern)
  vectorFrags <- vector()

  for( i in 1:length(lista)){
    vectorFrags[i] <-  lista[[i]][fragPos]
  }
  vectorFrags
}

sal.ref.convert <- function(ref_tessa,rep_calssif ){

tb <- read.csv(file=ref_tessa,sep = ";", header = F)

df <- data.frame(seqID=tb$V1,
repName=split.vec(as.character(tb$V5),"/",1),
repFamily=split.vec(as.character(tb$V5),"/",3),
repClass=split.vec(as.character(tb$V5),"/",2)
)

df[,c("seqID",rep_calssif)]
}

import.RTEs <- function(quant_dir, reference_file, conditions, summarize = c('repName','repFamily','repClass','none')){
  dir <- list.dirs(quant_dir, recursive = F, full.names = F)
  files <- file.path(quant_dir, dir,"quant.sf")
  names(files) <- dir

if(all(file.exists(files))==F){
  stop("error loading files")
  }else{
  if(summarize=='none'){
  
    txi_raw <- tximport::tximport(files, type = "salmon",txOut=T)
    
    message("Writing raw tables to output directory")
    write.table(txi_raw$counts,paste0(outdir,"/counts.csv"),quote = F,sep = "\t")
    write.table(txi_raw$length,paste0(outdir,"/length.csv"),quote = F,sep = "\t")
    write.table(txi_raw$abundance,paste0(outdir,"/tpm.csv"),quote = F,sep = "\t")
    
  }else{
    
    ref <- sal.ref.convert(reference_file, summarize)

    group.0 <- factor(conditions)
    des <- cbind(names(files),conditions)
    message("The following experimental design has been created:\n")
    message(message=for(i in 1:nrow(des)){print(paste(des[i,1],des[i,2]))})
  
    txi_rte <- tximport::tximport(files, type = "salmon",tx2gene = ref)

    message("Writing tables to output directory")
    write.table(txi_rte$counts,paste0(outdir,"/counts.csv"),quote = F,sep = "\t")
    write.table(txi_rte$length,paste0(outdir,"/length.csv"),quote = F,sep = "\t")
    write.table(txi_rte$abundance,paste0(outdir,"/tpm.csv"),quote = F,sep = "\t")

      }
    }
}

import.RTEs(quant_dir=outdir, 
            reference_file=ref_tessa, 
            conditions=conditions, 
            summarize=summarize)

