library(data.table)
library(rpart)
library(lattice)
library(ggplot2)
library(caret)
#packages for preprocess
library(rpart.plot)
#packages for decision tree
library(mlbench)
#package with data sample
library(tidyr)
library(rattle)
library(RColorBrewer)
library(randomForest)
library(h2o)

source('/home/ev250/Cincinatti/Functions/various.R')



#'  check genotype concordance after prediction by prediction output category
#'
#' This function allows you to assess concordance in genotype based on a binary output
#' @param x data table input for making train and test datasets, names as in dna.imp.info.errors from rna.geno.trees.R
#' @param y data table with predictions and SNP IDs, names as in pred.train in rna.geno.trees.R
#' @param z vector with names of categories from prediction output
#' @param w data table with counts for number of errors across all snps and all samples (marginal), example t1 in rna.geno.trees.R, optional
#' @keywords concordance prediction
#' @export
#' @return data table with genotype concordance (and cumulative) when applying machine learning rules
#' 
pred.conc <- function(x,y,z,w=NULL){
    
   tmp <- lapply(seq_along(z), function(i) x[rsid_ref_alt %in% y[pred.out==z[i],SNP],grep("N_errors", names(x), value=T),with=F])
    t <- lapply(seq_along(tmp), function(i) { DT <- data.table(table(rowSums(tmp[[i]])))
    DT[,Cum:=cumsum(N)*100/sum(N)]
    setnames(DT, names(DT)[2:3],paste0(names(DT)[2:3],".",z[i]))
    return(DT)
    })
    if(is.null(w)){
    t2 <- Reduce(function(x,y) merge(x, y, by="V1", all=T), t)
    t2[,V1:=as.numeric(V1)]
    setkey(t2,V1)
    # replace NA in .err with 0 (means 0 snp) and recalculate cum accordingly
    for(i in seq_along(z)){
        t2[is.na(get(paste0("N.",z[i]))), paste0("N.", z[i]) :=0]
        t2[,paste0("Cum.",z[i]):=cumsum(get(paste0("N.",z[i])))*100/sum(get(paste0("N.",z[i])))]
        }
        return(t2)
    } else {
        t[[3]] <- w
        t2 <- Reduce(function(x,y) merge(x, y, by="V1", all=T), t)
        t2[,V1:=as.numeric(V1)]
        setkey(t2,V1)
    # replace NA in .err with 0 (means 0 snp) and recalculate cum accordingly
    for(i in seq_along(z)){
        t2[is.na(get(paste0("N.",z[i]))), paste0("N.", z[i]) :=0]
        t2[,paste0("Cum.",z[i]):=cumsum(get(paste0("N.",z[i])))*100/sum(get(paste0("N.",z[i])))]
    }
        return(t2)
        }
        
}

#'  check % of genotype concordance per sample after prediction by prediction output category
#'
#' This function allows you to assess concordance in genotype based on a binary output
#' @param x data table input for making train and test datasets, names as in dna.imp.info.errors from rna.geno.trees.R
#' @param y data table with predictions and SNP IDs, names as in pred.train in rna.geno.trees.R, if y=NULL gives number of errors without prediction
#' @keywords concordance prediction per sample
#' @export
#' @return list of data tables with % of genotype concordance per sample 
#'
pred.conc.per <- function(x,y=NULL){
    if(is.null(y)){
        DT <- x[,grep("N_errors", names(x)),with=F]
        e <- lapply(0:2, function (x) apply(DT,2,function(j) sum(j==x)))
    #calc % of 0,1 or 2 errors (1:3 elements in e)
        perc <- lapply(1:3, function(w) sapply(seq_along(e[[1]]), function(x) e[[w]][[x]]*100/sum(sapply(e,`[[`,x))))
        DT2 <- data.table(Errors=0:2,Mean=round(sapply(perc,mean),0),SD=round(sapply(perc,sd)))
        return(DT2)
    
        } else {
            
    conc <- lapply(unique(y$pred.out), function (i) {
    #select SNPs
    DT <- x[rsid_ref_alt %in% y[pred.out==i,SNP],grep("N_errors", names(x)),with=F]
    # sum number of errors (0,1,2) across samples
    e <- lapply(0:2, function (x) apply(DT,2,function(j) sum(j==x)))
    #calc % of 0,1 or 2 errors (1:3 elements in e)
    perc <- lapply(1:3, function(w) sapply(seq_along(e[[1]]), function(x) e[[w]][[x]]*100/sum(sapply(e,`[[`,x))))
    DT2 <- data.table(Errors=0:2,Mean=round(sapply(perc,mean),0),SD=round(sapply(perc,sd)))
    return(DT2)
    })

    return(conc)
}
}
