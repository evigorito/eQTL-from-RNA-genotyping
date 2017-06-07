#Functions

library(data.table)
library(biomaRt)
library(xtable)
library(dplyr)
library(tidyr)
library(parallel)
library(rasqualTools)

######################## DNA ########################################

#' Format txt converted vcf files: one file per chr
#'
#' This function allows you to get good format from vcf text.
#' @param file file name with full path to vcf
#' @keywords vcf format
#' @export
#' @return named list with info from vcf files
#' name()
name <- function(file){
    
    temp <- fread(file, sep='\t')
    temp2 <- fread(paste0(gsub(".txt","",file),".header.txt"))
    names(temp)<- gsub("^.*]","",names(temp2))
    names(temp)<-gsub(":","_", names(temp))

    return(temp)

}

#' extract info from vcf genome wide, one file per chr
#'
#' This function allows you to get good format from vcf text.
#' @param path Path to files, defaul current dir
#' @param pattern of files default=NULL
#' @param name for each element of list (chr)
#' @keywords vcf format
#' @export
#' @return named list with info from vcf files
#' name_list()

name_list <- function(path=NULL, pattern=NULL, name){
    if(is.null(path)){path <- "."}
    files <- list.files(path,pattern=pattern, full.names=T)
    temp <- lapply(files, name)
    names(temp) <- name
    return(temp)
}


#' format annotation field 
#'
#' This function allows you to get simplified annotation from vcf to text conversion.
#' @param DT data table with vcf txt, output from name
#' @keywords annotation format
#' @export
#' @return data table with simplified annotation information
#' ann_vcf ()

ann_vcf <- function(DT){
    s <- strsplit(DT$ANN,",")
    temp <- lapply(s, function(i) unique(do.call(rbind,strsplit(i, "|", fixed=T))[,2]))
    #create hierarchy of terms
    terms <-  data.table(ANN=unique(unlist(temp)))
    terms[,ann.sum:=c(rep("proximal",2),"intergenic","intronic","exonic","exonic","UTR","exonic","intronic","unknown","UTR","exonic","exonic","intronic","proximal","exonic","exonic","exonic","intronic","intronic",rep("exonic",6),"intronic")]
    #hierarchy
    h <- c("exonic","UTR","intronic","proximal","intergenic","unknown")
    temp2 <- lapply(temp, function(i) unique(terms[ANN %in% i,ann.sum])[order(h)][1])
    DT[,ANN:=unlist(temp2)]
   
        return(DT)
}

#' format annotation field in list
#'
#' This function allows you to get simplified annotation  from vcf to text conversion with input a list of data tables
#' @param list of DT data tables with vcf txt, output from name
#' @keywords annotation format
#' @export
#' @return list of data tables with simplified annotation information
#' ann_vcf_list ()

ann_vcf_list <- function(list){

    temp <- lapply(list,ann_vcf)

    return(temp)
}


#' Categorise R2 and maf
#'
#' This function allows you to split by low, medium and high imp quality and rare, low and common variants.
#' @param DT data table with vcf txt, output from annot_vcf.
#' @keywords imputation maf categories
#' @export
#' @return data table with maf and imputation quality categories
#' cats ()

cats <- function(DT){
    DT[,R2_cat:="low"][R2>=0.6 & R2<0.9, R2_cat:="medium"][R2>=0.9,R2_cat:="high"]
    phase3 <-  fread(paste0("zcat /scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr",unique(DT$CHROM),".legend.gz"))
    temp <- merge(DT,phase3[,.(position,a0,a1,EUR)], by.x=c("POS","REF","ALT"), by.y=c("position","a0","a1"), all.x=T)
    temp[,maf_cat:="rare"][EUR>=0.01 &EUR <0.05,maf_cat:="low"][EUR<=0.99 & EUR >0.95, maf_cat:="low"][EUR>=0.05 & EUR <=0.95, maf_cat:="common"][is.na(EUR),maf_cat:="N/A"]
    
        return(temp)
}

#' Number of SNPs called per chr, stratified by imp2 info and maf
#'
#' This function allows you to: count called SNPs per chr based on maf and imputation q
#' @param DT data table with DNA info, one chromosome, output from cats.
#' @param col variable to stratify
#' @keywords snps DNA maf imputation
#' @export
#' @return list
#' count_snps ()
count_snps <- function(DT,col){
    temp <- lapply(unique(unlist(DT[,col,with=F])), function(i) {
        if(is.na(i)){
            DT[is.na(get(col))][,.N,by=ANN][,PCT:=N*100/sum(N)]
            } else {
                DT[get(col)==i,][,.N,by=ANN][,PCT:=N*100/sum(N)]
            }
        })
        
    temp2 <- Reduce(function(...) merge(...,by="ANN", suffixes=paste0("_",unique(unlist(DT[,col,with=F])))),temp)
    setnames(temp2, c("Region", sapply(unique(unlist(DT[,col,with=F])), function(i) paste0(i,c("_N","_PCT")))))
    return(temp2)
}

#' Number of SNPs called genome-wide, stratified by imp2 info and maf
#'
#' This function allows you to: count called SNPs genome wide based on maf or imputation q
#' @param list list with DNA info, output from count_snps
#' @param xtab if xtable is required, indicate variable for columns, dafult NULL
#' @param cap caption for xtable
#' @keywords snps DNA maf imputation
#' @export
#' @return data table
#' snps_gw ()
snps_gw<- function(list, xtab=NULL, cap=NULL){
    temp <- rbindlist(list)
    temp <- temp[,c("Region",grep("_N", names(temp),value=T)), with=F]
    setkey(temp,Region)
    temp2 <-temp[,lapply(.SD,sum), by=Region]
    temp2[,total:=rowSums(temp2[,2:ncol(temp2),with=F])]
    temp2[,PCT:=round(total*100/sum(total),2)]
    temp2 <- rbindlist(list(temp2,as.list(c("Total",colSums(temp2[,2:ncol(temp2),with=F])))))
    setnames(temp2,gsub("_N","",names(temp2)))
    if(is.null(xtab)!=TRUE){
        addtorow <- list()
        addtorow$pos <- list(0,0)
        addtorow$command <- c(paste0("&\\multicolumn{",ncol(temp2)-1,"}{c}{",xtab,"}\\\\\n"), paste0(paste(names(temp2)," ", collapse="& ")," \\\\\n"))
        print(xtable(temp2,caption=cap), booktabs=TRUE, include.rownames=FALSE,include.colnames=FALSE ,add.to.row=addtorow)
        }
    return(temp2)
}


############################### RNA and DNA comparison ##########################

#' Select common samples from dna and rna
#'
#' This function allows you to: select variants that are in the same position and same REF and ALT in DNA and RNA; Add column Concordance for genotype; Add column counting # of genotyping discrepancies between RNA and DNA.
#' @param dna data table with DNA info, one chromosome
#' @param rna data table with RNA data (all chrs)
#' @param chr chromosome to subset
#' @keywords common variants RNA DNA
#' @export
#' @return data table containing variants genotyped commonly by DNA and RNA, indicating concordance based on same genotype.
#' convert ()
convert <- function(dna,rna, chr){
    DT <- merge(rna[CHROM==chr,][,CHROM:=as.numeric(CHROM)], dna, by=c("CHROM","POS","REF","ALT"), suffixes=c("_RNA","_DNA"))

    dna_GT<-grep("GT", names(dna), value=TRUE)
    rna_GT<-sub("[0-9]*_","",dna_GT)
    prefix <- sub("_GT","",rna_GT)
    #cols<-paste0(cols_pre, "_GT")
for(i in 1:length(prefix)){
	DT[,paste0(prefix[i],"_Concordance"):="Discordant"]
	DT[get(paste0(rna_GT)) == "./.", paste0(prefix[i],"_Concordance"):=NA]
	DT[get(dna_GT[i])==get(rna_GT[i]), paste0(prefix[i],"_Concordance"):="Concordant"]	
	DT[,paste0(prefix[i],"_N_errors"):=10][get(paste0(prefix[i],"_Concordance"))=="Concordant", paste0(prefix[i],"_N_errors"):=0]
        DT[ get(dna_GT[i]) == "0/0" & get(rna_GT[i]) == "0/1", paste0(prefix[i],"_N_errors"):=1]
	DT[ get(dna_GT[i]) == "0/0" & get(rna_GT[i]) == "1/1", paste0(prefix[i],"_N_errors"):=2]
	DT[ get(dna_GT[i]) == "1/1" & get(rna_GT[i]) == "0/1", paste0(prefix[i],"_N_errors"):=1]
	DT[get(rna_GT[i]) == "./.", paste0(prefix[i],"_N_errors"):=NA]
	DT[ get(dna_GT[i]) == "0/1" & get(rna_GT[i]) == "1/1", paste0(prefix[i],"_N_errors"):=1]
	DT[ get(dna_GT[i]) == "0/1" & get(rna_GT[i]) == "0/0", paste0(prefix[i],"_N_errors"):=1]
	DT[ get(dna_GT[i]) == "1/1" & get(rna_GT[i]) == "0/0", paste0(prefix[i],"_N_errors"):=1]
        

    }
return(DT)

}    

#' get N-errors per sample for a chr
#'
#' This function allows you to: summarize errors across samples per chr
#' @param dna_rna_10 data table with DNA and RNA info, one chromosome, output from convert function
#' @param format if not NULL gives output formatted for table
#' @param col variable to stratify results on (default null).
#' @keywords common variants RNA DNA
#' @export
#' @return data table containing variants genotyped commonly by DNA and RNA, indicating concordance based on same genotype indicating 0,1 or 2 errors of RNA relative to DNA
#' N_error ()
N_error <- function(dna_rna_10,format=NULL,col=NULL) {
        sum_errors <- N_error_sub(dna_rna_10)
        temp1 <- N_error_sub2(sum_errors)
      if(is.null(col)==TRUE){  
        return(temp1)
    } else {
        u <- unique(unlist(dna_rna_10[,col, with=F]))
        temp <- lapply(u, function(i) dna_rna_10[get(col)==i,])
        names(temp) <- u
        summ_errors_l <- lapply(temp,N_error_sub)
        temp2 <- lapply(sum_errors_l, N_error_sub2)
        temp2 <- Reduce(function(...) merge(...,by="Errors", suffixes=paste0("_",u)),temp2)
       temp2 <- merge(temp2,temp1, by="Errors", suffixes=c(paste0("_",u[length(u)]),""))
        return(temp2)
    }
}

#' subfunction 1 for get N-errors per sample for a chr
#'
#' This function allows you to: summarize errors across samples per chr
#' @param dna_rna_10 data table with DNA and RNA info, one chromosome
#' @keywords common variants RNA DNA
#' @export
#' @return data table containing variants genotyped commonly by DNA and RNA, indicating concordance based on same genotype indicating 0,1 or 2 errors of RNA relative to DNA
#' N_error_sub ()
    N_error_sub <- function(dna_rna_10){
 samples <- gsub("_Concordance","",grep("_Concordance",names(dna_rna_10), value=T))
    DT<-data.table(Errors=c(0,1,2,NA))
    
    for (i in seq_along(samples)) {
    cols<-paste0(samples, "_GT")
    if(length(unique(dna_rna_10[get(cols[i])!=".",get(paste0(samples[i],"_N_errors"))]))==4){
                  DT[,samples[i]:=as.numeric(data.table(table(dna_rna_10[get(cols[i])!=".",get(paste0(samples[i],"_N_errors"))], useNA="always"))$N)]
              } else{
                  w <- which(DT$Errors %in% unique(dna_rna_10[get(cols[i])!=".",get(paste0(samples[i],"_N_errors"))]))
                  DT[,samples[i]:=0]
                  DT[w,samples[i]:=as.numeric(data.table(table(dna_rna_10[get(cols[i])!=".",get(paste0(samples[i],"_N_errors"))], useNA="always"))$N)]
              }
    }
 return(DT)

}
        
#' subfunction 2 for get N-errors per sample for a chr
#'
#' This function allows you to: summarize errors across samples per chr
#' @param DT data table with RNA errors relative to DNA
#' @param format returns table with formatted Mean (SD)
#' @keywords common variants RNA DNA
#' @export
#' @return data table containing variants genotyped commonly by DNA and RNA, indicating concordance based on same genotype indicating 0,1 or 2 errors of RNA relative to DNA
#' N_error_sub2 ()

N_error_sub2 <- function(DT, format=NULL) {
sum_errors<-data.table(Errors=c(0,1,2,NA), Mean=round(apply(DT[,2:ncol(DT)],1,mean)), SD=round(apply(DT[,2:ncol(DT)],1,sd)))
sum_errors[,Pct:=round(Mean*100/sum(Mean[1:3]),2)]

    if(is.null(format)==TRUE){
        return(sum_errors[1:3,])
        } else {
    sum_errors[,`Mean (SD)`:=paste0(Mean," (",SD,")")]
sum_errors[,Mean:=NULL][,SD:=NULL]
    setcolorder(sum_errors, c("Errors","Mean (SD)","Pct"))
    return(sum_errors[1:3,])
    
        }
}    


#' get N-errors per sample genome wide with or without stratification
#'
#' This function allows you to: summarize errors across samples per genome wide
#' @param dna_rna list with DNA and RNA info, list output from convert
#' @param col variable to stratify by, defualt NULL
#' @keywords errors variants RNA DNA
#' @export
#' @return data table containing variants genotyped commonly by DNA and RNA, indicating concordance based on same genotype indicating 0,1 or 2 errors of RNA relative to DNA
#' N_error_gw ()

N_error_gw <- function(dna_rna, col=NULL){
    temp <- sub_error_gw(dna_rna)
     if(is.null(col)){
    return(temp)
     } else{
         u <- unique(unlist(dna_rna[[1]][,col, with=F]))
         v1 <- list()
         for(i in seq_along(dna_rna)){
             v <- lapply(u, function(j) dna_rna[[i]][get(col)==j,])
             names(v) <- u
             v1[[i]] <- v
         }
         v2 <- list()
         for (i in seq_along(u)){
             v2[[i]] <- rbindlist(lapply(seq_along(v1), function(j) v1[[j]][[i]]))
             names(v2)[[i]] <- u[i]
         }
         v3 <- lapply(v2, N_error)
         v4 <- Reduce(function(...) merge(...,by="Errors"),v3)
         setnames(v4, c("Errors",unlist(lapply(u, function(i) paste0(names(v3[[1]])[-1],"_",i)))))
       v5 <- merge(v4,temp, by="Errors")
    return(v5)
     }
}
    
#' subfunction to get N-errors per sample genome wide with or without stratification
#'
#' This function allows you to: summarize errors across samples per genome wide
#' @param x list with DNA and RNA info, each element is for a chr stratified by col
#' @keywords errors variants RNA DNA
#' @export
#' @return data table containing variants genotyped commonly by DNA and RNA, indicating concordance based on same genotype indicating 0,1 or 2 errors of RNA relative to DNA
#' sub_error_gw ()

sub_error_gw <- function(x){
    temp <- rbindlist(lapply(x,N_error_sub))
    temp <- temp[!is.na(Errors),]
    setkey(temp,Errors)
    temp2 <- temp[,lapply(.SD,sum), by=Errors]
    samples <- grep("Error",names(temp), invert=TRUE, value=TRUE)
    temp2[,Mean:=round(rowMeans(temp2[,samples, with=F]),0)]
    temp2[,SD:=round(apply(temp2[,samples,with=F], 1, sd),0)]
    temp2[,PCT:=round(Mean*100/sum(Mean),2)]
    temp2[,(samples):=NULL] 
    return(temp2)
}

#' Get concordance between DNA & RNA by genomic location stratified per N-errors per sample  with or without stratification
#'
#' This function allows you to: get overlap between RNA &  DNA by genomic location with or without further stratification
#' @param dna_rna list with DNA and RNA info, each element is for a chr, output from convert
#' @param col column to stratify by
#' @param Nerrors number of errors between RNA and DNA, default concordance
#' @param xtab whether the output also includes xtable. If so, take name for top title
#' @param cap caption for xtable
#' @keywords errors variants RNA DNA
#' @export
#' @return data table
#' conc_annot ()

conc_annot <- function(dna_rna,col=NULL,Nerrors=0, xtab=NULL, cap=NULL){

    a <-  rbindlist(lapply(dna_rna, function(i) sub_conc_annot(i,col=NULL,Nerrors)))
    a <- setkey(a,Region)
    b <- a[,lapply(.SD,sum,na.rm=TRUE),by=Region, .SDcols=names(a)[-1]]
    b[,Mean:=round(apply(b[,names(b)[-1],with=FALSE],1,mean,na.rm=TRUE),0)]
    b[,SD:=round(apply(b[,names(b)[-1],with=FALSE],1,sd,na.rm=TRUE),0)]
    b[,PCT:=round(Mean*100/sum(Mean,na.rm=TRUE),2)]
    b[, (names(a)[-1]):=NULL]
    if(is.null(col)==TRUE){
        return(b)
        } else {
    x <- rbindlist(lapply(dna_rna, function(i) sub_conc_annot(i,col,Nerrors)))
    g <- grep("N_errors",names(x),value=T)
    setkeyv(x,c("Region",col))
    temp <- x[,lapply(.SD,sum,na.rm=TRUE),by=.(Region,get(col)), .SDcols=g]
    temp[,Mean:=round(rowMeans(temp[,g,with=F],na.rm=TRUE),0)]
    temp[,SD:=round(apply(temp[,g,with=F],1,sd,na.rm=TRUE),0)]
    temp[,(g):=NULL]
    n <- unique(temp$get)
    temp2 <- dcast(temp,Region~get, value.var=list("Mean","SD"))
    setcolorder(temp2,c("Region",unlist(lapply(n, function(i) paste0(c("Mean","SD"), "_",i)))))
    temp3 <- merge(temp2,b, by="Region")
    temp3 <-  rbindlist(list(temp3,as.list(c("Total",colSums(temp3[,2:ncol(temp3),with=F], na.rm=TRUE)))))

    if(is.null(xtab)!=TRUE){
        temp4 <- copy(temp3)
        l <- lapply(n,function(i) grep(i,names(temp4),value=T))
        for(i in seq_along(l)){
            temp4[,l[[i]][1]:=paste0(get(l[[i]][1])," (",get(l[[i]][2]),")")]
        }
        temp4[,Mean:=paste0(Mean," (",SD,")")]
        temp4[,grep("SD", names(temp4),value=T):=NULL]
        addtorow <- list()
        addtorow$pos <- list(0,0)
        addtorow$command <- c(paste0("&\\multicolumn{",ncol(temp4)-2,"}{c}{",xtab," (Mean, SD)}\\\\\n"), paste0(paste(c(names(temp4)[1],n,"total","PCT")," ", collapse="& ")," \\\\\n"))
        print(xtable(temp4, caption=cap), booktabs=TRUE, include.rownames=FALSE,include.colnames=FALSE ,add.to.row=addtorow)
        }
    return(temp3)
     
        }
}

 
#' sub function to conc_annot get concordance between DNA & RNA by genomic location  per sample  with or without stratification
#'
#' This function allows you to: get overlap between RNA &  DNA by genomic location with or without further stratification
#' @param DT data table with DNA and RNA info for a chr, output from convert
#' @param col column to stratify by
#' @keywords errors variants RNA DNA
#' @export
#' @return data table
#' sub_conc_annot ()

sub_conc_annot <- function(DT,col=NULL,Nerrors=0){
        g <- grep("N_errors",names(DT),value=T)
        if(is.null(col)){
            temp <- lapply(g, function(i) {
                t <- DT[get(i)==Nerrors,]
                           setkey(t,ANN)
                           t2 <- t[,.N,by=ANN]
            setnames(t2,c("Region",i))
            })
            t3 <- Reduce(function(...) merge(...,by="Region", all=FALSE),temp)
            return(t3)
                           
            } else {
        temp <- lapply(g,function(i) {
            t <- DT[get(i)==Nerrors,]
            setkeyv(t,c("ANN",col))
            t2 <- t[,.N,by=.(ANN,get(col))]
            setnames(t2,c("Region",col,i))
            })
         t3 <- Reduce(function(...) merge(...,by=c("Region",col),all=FALSE),temp)
        
         return(t3)
            }
}

#' Coverage of RNA vs DNA by genomic location  with or without further stratification
#'
#' This function allows you to: get RNA coverage by genomic location with or without further stratification
#' @param dna_sum data table with summary DNA info w/wo stratification, output from snps_gw
#' @param dna_rna_sum data table with variants called in RNA concordant with DNA, output conc_annot
#' @keywords errors variants RNA DNA
#' @export
#' @return data table
#' cov_annot ()

    cov_annot <- function(dna_sum,dna_rna_sum){
        temp <- dna_rna_sum[,c("Region",grep("Mean",names(dna_rna_sum), value=T)), with=F]
        temp2 <- merge(dna_sum,temp, by="Region")
        temp2[,PCT:=NULL]
        setnames(temp2,"Mean","Mean_total")
        n <- names(dna_sum)[2:(ncol(dna_sum)-1)]
        l <- lapply(n, function(i) grep(i,names(temp2), value=T))

        for(i in seq_along(l)){
            temp2[,paste0(l[[i]][1],"_PCT"):=round(as.numeric(get(l[[i]][2]))*100/as.numeric(get(l[[i]][1])),1)]
        }
        temp2[,unlist(l):=NULL]
        return(temp2)
    }

#' Combine coverage of RNA before & after imputation vs DNA by genomic location  with or without further stratification
#'
#' This function allows you to: compare RNA coverage with/without imputation by genomic location with or without further stratification
#' @param noimp_ov data table comparing rna no imputation with DNA, output from cov_annot
#' @param imp_ov data table comparing rna imputed with DNA, output cov_annot
#' @param xtab whether to include xtable in output, provide title
#' @param cap caption for xtable
#' @keywords errors variants RNA DNA
#' @export
#' @return data table with or without xtable
#' ov_comp ()

ov_comp <- function(noimp_ov,imp_ov,xtab=NULL,cap=NULL){
    temp <- merge(noimp_ov,imp_ov,by="Region", suffixes=c("_NO","_YES"))
    n <- gsub("_PCT","",names(noimp_ov)[-1])
    l <- unlist(lapply(names(noimp_ov)[-1], function(i) grep(i,names(temp),value=T)))
    setcolorder(temp,c(names(temp)[1],l))
    if(is.null(xtab)==FALSE){
    addtorow <- list()
        addtorow$pos <- list(0,0,0,0)
    addtorow$command <- c(paste0("\\multicolumn{",ncol(temp)-1,"}{c}{Overlap with DNA by ",xtab,"} \\\\\n"),paste0(paste0("&\\multicolumn{2}{c}{",n," (\\%)}"," ", collapse=" ")," \\\\\n"),paste0("RNA imp &", paste0(rep(c(" NO "," YES "), length(n)), " ", collapse="&"), "\\\\\n"),
                          paste0(c("Region",rep("&",ncol(temp)-1), "\\\\\n"), collapse= " "))
                          
        print(xtable(temp,caption=cap), booktabs=TRUE, include.rownames=FALSE,include.colnames=FALSE ,add.to.row=addtorow)
    }
    return(temp)

}

######## Calculating probs for inferring genotypes ##############

#########  subfunction ############
##' Handy function, logsum
##'
##' This function calculates the log of the sum of the exponentiated
##' logs taking out the max, i.e. insuring that the sum is not Inf
##' @title logsum
##' @param x numeric vector
##' @return max(x) + log(sum(exp(x - max(x))))
##' @author Claudia Giambartolomei
logsum <- function(x) {
  my.max <- max(x)   ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max ))) 
  return(my.res)
}


#' open and selecting chrs to work after calling variants
#'
#' This function allows you to: open a file and select variants in a particular chormosome
#' @param file file with data, converted from vcf format
#' @param chr chromosome of interest to subset
#' @keywords select variants by chromosome
#' @export
#' @return data table
#' open_select ()

open_select <- function(file,chr){
    #calls name function defined above
    temp <- name(file)
    temp <- temp[CHROM %in% chr,]

    return(temp)

    }

#' subfunction to geno_chr to calculate genotype probabilities for samples with low read depth to inform imputation
#'
#' This function allows you to: calculate probs for genotyping using allele frequency
#' @param rna rna called variants converted from vcf format and opened with open_select
#' @param samples vector with the names of samples of interest
#' @param dp value for  depth of reads to use to compute probs
#' @param e value with probability of genotyping error
#' @param chr chromosome of interest to subset
#' @param phase3 data table with maf info for relevant chr
#' @keywords select variants by chromosome
#' @export
#' @return data table with calculated probabilities for genotypes and genotypes transformed to 0,1,2 scale.
#' geno_prob_e ()

geno_prob_e<-function(rna,samples, dp, e, chr, phase3) {
    all <- list()
   	for(i in seq_along(samples)){
	# for each sample select entries with DP=dp;
	temp<-rna[get(paste0(samples[i],"_DP"))==dp & CHROM==chr,c(1:6,grep(samples[i],names(rna))), with=F]
	# add maf to temp
       
	m<-merge(temp, phase3[,c("id","position","a0","a1","EUR"), with=F], by.x="POS", by.y="position", all.x=TRUE)
	# make sure I select the same snps and maf >0
        if(sum(names(m)=="REF_RNA")!=0){
            
            m1<-m[REF_RNA==a0  & ALT_RNA==a1 & EUR>0 & EUR<1,]
             } else {
            m1 <- m[REF==a0 & ALT==a1 &EUR>0 & EUR<1,]
        }
	#calculate p; AD is read in "ref,alt" alleles
	b<-as.numeric(unlist(strsplit(m1[,get(paste0(samples[i],"_AD"))], ","))[c(FALSE, TRUE)])
	a<-as.numeric(unlist(strsplit(m1[,get(paste0(samples[i],"_AD"))], ","))[c(TRUE,FALSE)])
	p<-m1[,EUR]
        # work in log scale to make it simpler
	ll_00<-a*log(1-e)+b*log(e)+2*log(1-p) 
	ll_01<-log(0.5)*a+log(0.5)*b+log(2*p*(1-p)) 
	ll_11<-a*log(e) + b*log(1-e) + 2*log(p)
         #l_00<-((1-e)^a)*(e^b)*(1-p)^2
	 #l_01<-(0.5^a)*0.5^b*2*p*(1-p)
	 #l_11<-(e^a)*(1-e)^b*p^2
         #rescale for prob (get rid of k, see formulae above), and transform to scale 0,2 (0.p(0)+1.p(1)+2.p(2)
         tot<-apply(rbind(ll_00,ll_01,ll_11),2, logsum)
         m1[,prob.00:=exp(ll_00)/exp(tot)]
        m1[,prob.01:=exp(ll_01)/exp(tot)]
        m1[,prob.11:=exp(ll_11)/exp(tot)]
        #m1[,test:=prob.00+ prob.01+ prob.11]
        #select info
        m2<-m1[,.( CHROM, POS, REF,ALT, prob.00,    prob.01,      prob.11)]
         all[[i]]<-m2
         names(all)[[i]]<-samples[i]
    }
    all_dt<-rbindlist(all, id=T)
    names(all_dt)[which(names(all_dt)==".id")]<-"sample"
    return(all_dt)
}


#' subfunction to geno_chr to calculate genotype probabilities for samples with low read depth to inform imputation
#'
#' This function allows you to: calculate probs for genotyping using allele frequency for a vector of "e"
#' @param rna rna called variants converted from vcf format and opened with open_select
#' @param samples vector with the names of samples of interest
#' @param dp single value with depth of reads to use to compute probs
#' @param e vector with probability of genotyping error
#' @param chr chromosome of interest to subset
#' @param phase3 data table with maf info for relevant chr
#' @keywords select variants by chromosome
#' @export
#' @return data table with calculated probabilities for genotypes and genotypes transformed to 0,1,2 scale for a vector of e.
#' geno_e ()

geno_e<-function(rna,samples, dp=1,e, chr,phase3){

   l<-lapply(e,function(i)  geno_prob_e(rna,samples,dp,i, chr,phase3))
   dt<-rbindlist(l)
   return(dt)

}

#' subfunction to geno_chr to calculate genotype probabilities for samples with low read depth to inform imputation
#'
#' This function allows you to: calculate probs for genotyping using allele frequency for a vector of "dp"
#' @param rna rna called variants converted from vcf format and opened with open_select
#' @param samples vector with the names of samples of interest
#' @param dp vector with  values with depth of reads to use to compute probs
#' @param e singel value probability for genotyping error
#' @param chr chromosome of interest to subset
#' @param phase3 data table with maf info for relevant chr
#' @keywords select variants by chromosome
#' @export
#' @return data table with calculated probabilities for genotypes and genotypes transformed to 0,1,2 scale for a vector of dp.
#' geno_dp ()

geno_dp <- function(rna,samples,dp,e, chr, phase3){

    l <- lapply(dp, function(i) geno_prob_e(rna,samples,dp=i,e, chr,phase3))
    names(l) <- dp
    dt <- rbindlist(l, id=T)
    return(dt)
}

#' Calculate genotype probabilities for samples with low read depth to inform imputation
#'
#' This function allows you to: get maf info for chrs of interest
#' @param rna rna called variants converted from vcf format and opened with open_select
#' @param samples vector with the names of samples of interest
#' @param dp vector of  values with depth of reads to use to compute probs
#' @param e vector of  values for the probability of genotyping error
#' @param chr vector with chromosomes of interest to subset
#' @keywords select variants by chromosome
#' @export
#' @return data table with calculated probabilities by dp and/or e, if dp or e are single values, a list with each element being a data table with probs for the range of dp per value of e
#' geno_chr ()

geno_chr <- function(rna,samples,dp,e, chr){
    if(length(dp>1) & length(e==1)){
   
                         l <- lapply(chr, function(i){
                             phase3 <-  fread(paste0("zcat /scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr",i,".legend.gz"))
                             geno_dp(rna,samples,dp,e,i,phase3)})
        return(rbindlist(l))
                         }
    if(length(dp)==1 & length(e)>1){
                           l <- lapply(chr, function(i){
                             phase3 <-  fread(paste0("zcat /scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr",i,".legend.gz"))
                             geno_e(rna,samples,dp,e,i,phase3)})
                            return(rbindlist(l))
    }
    if(length(dp)>1 & length(e)>1){
        temp <- list()
        for(j in e){
             l <- lapply(chr, function(i){
                             phase3 <-  fread(paste0("zcat /scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr",i,".legend.gz"))
            geno_dp(rna,samples,dp,j,i,phase3)
       
             })
             temp[[j]] <- rbindlist(l)
        }
        names(temp) <- e
        return(temp)
                     }
}
                                        
#' open rna called variants in bzip imp2 format (gen.gz) and  modify genotypes
#'
#' This function allows you to: Replace probabilities calculated by geno_chr in file with rna called variants to be used for imputation.
#' @param file path to file to gen.gz to open. Also requires ".samples" file in same location. This is the vcf file after calling variants converted to imp2 format for ease.
#' @export
#' @return data table with genotype replacements
#' new_genos ()

new_genos <- function(file, probs){
    gen.file <- fread(paste0("zcat ",file))
    sam <- paste0(gsub("gen.gz","",file),"samples")
    sam.file <- fread(sam)
#V6-V8 in gen.file correspond to first sample in vcf.file.

    #add column in prob.genos to match with gen.file
# may be quicker with set function for replacing values
probs[,chr.pos.ref.alt:=paste(paste(CHROM,POS,sep=":"),REF, ALT, sep="_")]
 for( i in 2:nrow(sam.file)){
    temp<-probs[sample==sam.file$ID_1[i],]
    setkey(temp,CHROM,POS)
    #gen.file has some positions out of order, which complicates matching
    mer <- merge(gen.file[,paste0("V",1:5), with=F],temp, by.x="V1", by.y="chr.pos.ref.alt", all.y=T, sort=F)
    w<-which(gen.file$V1 %in% mer$V1)
    #wt<-which(mer$chr.pos.ref.alt %in% gen.file$V1)
    if(identical(gen.file$V1[w], mer$V1)==TRUE) {
       #first sample is gen.file is in columns V6-V8.
        #samples in sam.file start in row 2
       n.sample=i-2
       cols=(n.sample*3+6):(n.sample*3+8)
       gen.file[w,paste0("V",cols):=mer[,.(prob.00,prob.01,prob.11)]]

   } else {

       print(paste("sample",i, "not identical"))

   }
    #print(paste("sample",i))

}
 return(gen.file)
}   


#################  RNA imputed files ###############################
# get imputed files
#' Reads and format imp2 files, each file one chr
#'
#' This function allows you to get good format for imp2 files with info score >=0.4
#' @param path path to imp2 file
#' @param file file name of imputed file
#' @keywords imp2 format
#' @export
#' @return data table 
#' read_imp
read_imp <- function(file, path='/scratch/ev250/Cincinatti/quant/imputation/imp_full_chr/') {
    
    temp <- fread(paste0(path, file))
    samples <- fread('/scratch/ev250/Cincinatti/quant/imputation/input/star_built37trimmed.imp2.samples')
    samps <- rep(samples$ID_1[2:nrow(samples)],3)
    samps <- samps[ order(match(samps,samples$ID_1[2:nrow(samples)]))]
    names(temp) <- c("imp","rsid_ref_alt","POS","ref","alt" ,paste(samps, 1:3, sep="_"))
    info <- fread(paste0(path,sub("\\.dosage","",file), '_info'), colClasses =list(character=3:10))
    info <- info[info>=0.4,]
    temp2 <- temp[rsid_ref_alt %in% info$rs_id,]
    return(temp2) 
}

#' Reads and format imp2 files, genome wide
#'
#' This function allows you to get good format for imp2 files, genome wide
#' @param path path to imp2 files
#' @param files file names of imputed file
#' @keywords imp2 format
#' @export
#' @return named list each element one chr
#' read_imp_gw

read_imp_gw <- function(path='/scratch/ev250/Cincinatti/quant/imputation/imp_full_chr/',files){
    imps <- lapply(files,function(i) read_imp(path, file=i))
    names(imps) <- paste0("chr",sub(".*\\.","",sub("\\.imp2","",files)))
    return(imps)
}


#' Compares RNA imputed file with DNA, one chr
#'
#' This function allows you to select common samples, select common alleles and recode imputed probabilities to 0,1,2 scale (0.p(0)+1.p(1)+2.p(2)) in both datasets to ease comparison: imp_convert function.
#' @param imp_10 DT with imputed data for a chr, output from read_imp
#' @param dna_chr10 DT with DNA info for a chr, output from convert
#' @keywords imp2 format
#' @export
#' @return data table comparing DNA with RNA
#' imp_convert
imp_convert <- function(imp_10,dna_chr10){
    dna_GT<-grep("GT", names(dna_chr10), value=TRUE)
    prefix <- sub("_GT","",sub("[0-9]*_","",dna_GT))
    cols <- unlist(lapply(prefix, function(i) grep(i,names(imp_10), value=T)))
    dna_rna <- merge(imp_10[,c(names(imp_10)[1:5],cols),with=F], dna_chr10, by.x=c("POS","ref","alt"),
                     by.y=c("POS","REF","ALT"))
    # transform to scale 0,2 (0.p(0)+1.p(1)+2.p(2))
    l <- lapply(prefix, function(i) unlist(dna_rna[,paste0(i,"_2"),with=F] + 2*dna_rna[,paste0(i, "_3"), with=F], use.names=F ))
    temp <- dna_rna[,(paste0(prefix, "_RNA")):=l]
    
    dna_gp <- grep("_GP", names(temp), value=T)
    l <- lapply(dna_gp, function(i) as.numeric(unlist(strsplit(unlist(temp[,i,with=F]), ","))[c(FALSE,TRUE,FALSE)]) + 2*as.numeric(unlist(strsplit(unlist(temp[,i,with=F]), ","))[c(FALSE,FALSE,TRUE)]))
    temp[,(paste0(gsub("_GP","",dna_gp), "_DNA")):=l]

    return(temp)

}

#' Compares RNA imputed file with DNA, genome wide
#'
#' This function allows you to select common samples, select common alleles and recode imputed probabilities to 0,1,2 scale (0.p(0)+1.p(1)+2.p(2)) in both datasets to ease comparison: imp_convert function.
#' @param imps list of DTs, each DT  with imputed data for a chr, output from read_imp
#' @param DNA_f list of DTs, each DT with DNA info for a chr, output from convert
#' @keywords imp2 format
#' @export
#' @return list comparing DNA with RNA, each element one chr
#' imp_convert_gw

imp_convert_gw <- function(imps,DNA_f) {
    s <- sapply(names(imps), function(i) which(i==names(DNA_f)))
    temp <- lapply(seq_along(imps), function(i) imp_convert(imps[[i]],DNA_f[[s[i]]]))
    return(temp)
}


#' Calculates number of errors DNA vs RNA, one chr
#'
#' This function allows you to compute Concordance and 1 or 2 discrepancies between DNA and RNA called variants
#' @param imp_rna_dna_10 DT merging DNA with RNA data for a chr, output from imp_convert
#' @keywords imp2 DNA RNA errors
#' @export
#' @return data table comparing DNA with RNA
#' N_error_imp

N_error_imp<- function(imp_rna_dna_10) {
    rna <- grep("_RNA",names(imp_rna_dna_10), value=T)
    dna <- grep("_DNA", names(imp_rna_dna_10), value=T)
    prefix <- gsub("_RNA","",rna)
    imp_rna_dna_10 <- N_error_imp_sample(imp_rna_dna_10)
     DT<-data.table(Errors=c(0,1,2)) 
    for (i in seq_along(rna)) {      
    if(length(unique(imp_rna_dna_10[,get(paste0(prefix[i],"_N_errors"))]))==3){
                  DT[,prefix[i]:=as.numeric(data.table(table(imp_rna_dna_10[,get(paste0(prefix[i],"_N_errors"))]))$N)]
              } else{
                  w <- which(DT$Errors %in% unique(imp_rna_dna_10[,get(paste0(prefix[i],"_N_errors"))]))
                  DT[,prefix[i]:=0]
                  DT[w,prefix[i]:=as.numeric(data.table(table(imp_rna_dna_10[,get(paste0(prefix[i],"_N_errors"))]))$N)]
              }
}

sum_errors<-data.table(Errors=c(0,1,2), Mean=round(apply(DT[,2:ncol(DT)],1,mean)), SD=round(apply(DT[,2:ncol(DT)],1,sd)))
sum_errors[,Pct:=round(Mean*100/sum(Mean[1:3]),2)]
    #sum_errors[,`Mean (SD)`:=paste0(Mean," (",SD,")")]
#sum_errors[,Mean:=NULL][,SD:=NULL]
    #setcolorder(sum_errors, c("Errors","Mean (SD)","Pct"))

 return(sum_errors)
    
}

#' Calculates number of errors DNA vs RNA, one chr, subfunction for N_error
#'
#' This function allows you to compute Concordance and 1 or 2 discrepancies between DNA and RNA called variants
#' @param imp_rna_dna_10 DT merging DNA with RNA data for a chr, output from imp_convert
#' @keywords imp2 DNA RNA errors
#' @export
#' @return original data table with columns indicating # of errors per sample
#' N_error_imp_sample

N_error_imp_sample <- function(imp_rna_dna_10){
     rna <- grep("_RNA",names(imp_rna_dna_10), value=T)
    dna <- grep("_DNA", names(imp_rna_dna_10), value=T)
    prefix <- gsub("_RNA","",rna)
    l <- lapply(seq_along(rna), function(i) round(abs(round(unlist(imp_rna_dna_10[,rna[i],with=F]))-round(unlist(imp_rna_dna_10[,dna[i],with=F])))))
     imp_rna_dna_10[,(paste0(prefix,"_N_errors")):=l]
     return(imp_rna_dna_10)
}


#' Calculates number of errors DNA vs RNA, genome wide
#'
#' This function allows you to compute Concordance and 1 or 2 discrepancies between DNA and RNA called variants
#' @param dna_imps list merging DNA with RNA data for a chr, output from imp_convert
#' @keywords imp2 DNA RNA errors
#' @export
#' @return data table comparing DNA with RNA
#' N_error_imp_gw

N_error_imp_gw <- function(dna_imps){
    temp <- rbindlist(lapply(dna_imps, N_error_imp))
    setkey(temp,Errors)
    temp2 <- temp[,.(Mean=round(mean(Mean),0), SD=round(sd(Mean),0)),by=Errors]
    temp2[,PCT:=round(Mean*100/sum(Mean),1)]
    
        return(temp2)

}

######################## Uncertanty in haplotypes ####################################

#' Calculates haplotype uncertanty from shapeit2 haplotype estimation
#'
#' This function allows you to average haplotype estimates by sampling from shapeit2 haplotype graphs and prepare format for RASQUAL (the custom allelic probability (AP) field consists of the allelic probabilities  of two haplotypes from one individual, separated by a comma:
#' @param file full path to shapit2 formatted output file
#' @keywords haplotype uncertanty
#' @export
#' @return data table with haplotype probabilities
#' hap_un

hap_un <- function(file){
    temp <- fread(file, header=F)
    setkey(temp,V1)
    DT <- temp[,lapply(.SD,mean),by=V1,.SDcols=paste0("V",7:(ncol(temp)))]
    temp2 <- temp[,.SD[1],by=V1, .SDcols=paste0("V",2:6)]
    temp3 <- merge(temp2, DT , by="V1")
    temp3[,V1:=NULL]
    for(i in seq(7,ncol(temp3)+1, 2)){
        temp3[,(paste0("V",i)):=list(paste0(get(paste0("V",i)),",",get(paste0("V",i+1))))]
    }
    temp3[,paste0("V",seq(8,ncol(temp3)+1, 2)):=NULL]
    return(temp3)
  
}




####################### RASQUAL input  #############################################
#' Function to prepare vcf from DNA files with GT:ASE and RSQ fields
#'
#' Format individual samples per chr with GT:AS and RSQ fields as required for vcf
#' @param path path to tab files with GT info and ASE output
#' @param pattern to match GT tab files
#' @param chr vector with chr of interest
#' @param path to file with R object containing a DNA list with each element corresponding to one chr, output from ann_vcf_list
#' @keywords vcf for RASQUAL
#' @export
#' @return saves tab delimited files for RASQUAL input per sample per chr 
#' DNAvcf4rasqual

DNAvcf4rasqual <- function(path,pattern,chr,DNA){
    DNA_f <- readRDS(DNA)
    lapply(chr, function(i) DNAvcf4rasqual_sub(path=path,pattern=pattern,chr=i,DNA_f=DNA_f[[which(names(DNA_f)==i)]]))

}

#' Subfunction to prepare vcf from DNA files with GT:ASE and RSQ fields
#'
#' Format individual samples per chr with GT:AS and RSQ fields as required for vcf
#' @param path path to tab files with GT info and ASE output
#' @param pattern to match GT tab files
#' @param chr chromosome of interest
#' @param DNA_f DNA data table with R2 information
#' @keywords DNA format RASQUAL
#' @export
#' @return saves tab delimited file with RASQUAL formatted input per sample for specified chr
#' DNAvcf4rasqual_sub

DNAvcf4rasqual_sub <- function(path,pattern,chr,DNA_f) {
    #open GT files
    setwd(path)
    t_files <- list.files(pattern=paste0(chr,".",pattern))
    tabs <- lapply(t_files,fread,header=F)
    names(tabs) <- gsub(paste0("\\.",pattern),"", t_files)
    #open phaser output with allele specific counts
    pat <- lapply(names(tabs), strsplit, split=".", fixed=TRUE)
    AS <- lapply(seq_along(pat), function(i) fread(paste0(pat[[i]][[1]][2],".",pat[[i]][[1]][1],".allelic_counts.txt")))
    AS <- lapply(seq_along(AS), function(i) setkey(AS[[i]], position))
    names(AS) <- names(tabs)
    tabs_AS <- lapply(seq_along(tabs), function(i) merge(tabs[[i]], AS[[i]], by.x=paste0("V",c(1:2,4:5)), by.y=c("contig","position","refAllele","altAllele"), all.x=T))  
     #format
    for(i in 1:length(tabs_AS)) {
        #replace NA with 0
        tabs_AS[[i]][is.na(refCount), refCount:=0][is.na(altCount),altCount:=0]
        # merge refCount with altCount
        tabs_AS[[i]][,GT:=rep("GT:AS",nrow(tabs_AS[[i]]))]
        tabs_AS[[i]][,V9:=paste0(V9,":",refCount,",",altCount)]
        #remove columns and reorder
        tabs_AS[[i]][,c('variantID', 'refCount', 'altCount', 'totalCount'):=NULL]
        setcolorder(tabs_AS[[i]], c(names(tabs[[i]])[1:8],"GT","V9"))
        setkey(tabs_AS[[i]], V2)
}

     # add info score
    
    for (i in seq_along(tabs_AS)) {
        #n <- which(names(DNA_f)==pat[[i]][[1]][2])
    tabs_AS[[i]] <- merge(tabs_AS[[i]], DNA_f[,.(CHROM,POS,REF,ALT,R2)], by.x=paste0("V",c(1:2,4:5)), by.y=c("CHROM","POS", "REF", "ALT" ), all.x=T)
    tabs_AS[[i]][,V8:=paste0("RSQ=",R2)][,R2:=NULL]
    setcolorder(tabs_AS[[i]],c(names(tabs[[i]])[1:8],"GT","V9"))
    }
    # save as tab delimed
    names(tabs_AS) <- names(tabs)
    for(i in seq_along(tabs_AS)){
        write.table(tabs_AS[[i]], file=paste0(names(tabs_AS)[[i]],".for.AS.tab"),row.names=F,col.names=F,quote=F,sep="\t")
    
    }

}

#' Function to prepare vcf from RNA genotyped files with GT:ASE and RSQ fields
#'
#' Format individual samples per chr with GT:AS and RSQ fields as required for vcf
#' @param path path to txt files with GT info and ASE output, files output from vcf4AS (bash)
#' @param pattern to match GT tab files
#' @param chr vector with chrs of interest, as indicated in teh input files
#' @param pat2 pattern match phaser output files, files in same dir as GT txt files
#' @param samples vector with the name of the samples to work with. should be part of the phaser output file name
#' @param info_path path to info score imp2 output file
#' @param info_pat pattern to match info file
#' @param prefix prefix for output file excluding chr (added automatically)
#' @keywords vcf for RASQUAL
#' @export
#' @return saves tab delimited files for RASQUAL input per sample per chr 
#' RNAvcf4rasqual
#'
RNAvcf4rasqual <- function(path,pattern,chr,pat2, samples, info_path, info_pat,prefix){
    setwd(path)
    info_files <- list.files(path=info_path,pattern=info_pat, full.names=TRUE)
    RNA <- lapply(chr,
                  function(i) RNAvcf4rasqual_sub(pattern,
                                                 i,
                                                pat2,
                                                samples,
                                                info=fread(grep(paste0("\\.",i, "\\."),info_files, value=T), colClasses =list(character=3:10)),
                                                prefix))

}

#' Subfunction to prepare vcf from RNA files with GT:ASE and RSQ fields
#'
#' Format individual samples per chr with GT:AS and RSQ fields as required for vcf
#' @param pattern to match GT txt files
#' @param chr a chromosomes to match
#' @param pat2 pattern match files with phaser output, files in same dir as GT txt files
#' @param samples vector with name of the samples to work with. should be part of the phaser output file name
#' @param info data table with info score for the chr of interest
#' @param prefix prefix for output file (excluding chr)
#' @keywords format GT:AS Rasqual
#' @export
#' @return saves tab delimited file with RASQUAL formatted input per sample for specified chr
#' RNAvcf4rasqual_sub

RNAvcf4rasqual_sub <- function(pattern,chr,pat2,samples,info,prefix) {
    #open GT files
    
    t_file <- name(list.files(pattern=paste0("^",chr,"\\.",pattern)))
    #open phaser output with allele specific counts in a list per chr
    files <- list.files(pattern=pat2)
    f_chr <- grep(pattern=paste0("^",chr,"\\."), x=files,value=T)
    AS <- lapply(f_chr,fread)
    AS <- lapply(seq_along(AS), function(i) setkey(AS[[i]], position))
    u <- unlist(strsplit(f_chr,"\\."))
    names(AS) <- u[which(u %in% samples)]
    #format AS

    AS <- lapply(seq_along(AS),function(i) AS[[i]][,GT:=paste0(refCount,",",altCount)][,c('refCount', 'altCount', 'totalCount'):=NULL])
    names(AS) <- u[which(u %in% samples)]
    AS <- rbindlist(AS, idcol="sample") 

    #wide
    setkey(AS,variantID)
    AS2 <- dcast(AS, contig + position + variantID + refAllele + altAllele ~ sample, value.var="GT")

    #checking each row has at least 1 value that !is.na
    fun <- function(x){sum(!is.na(x))}
    s <- apply(AS2[,6:ncol(AS2)], 1, fun)
    if(sum(s==0)!=0)
        stop("missing values for AS for all samples in at least one variant")

    #replace NA with 0,0
    AS2[is.na(AS2)] <- "0,0"

    #merge AS2 with GT
    AS_GT <- merge(t_file,AS2, by.x=c("CHROM","POS","REF","ALT"), by.y=c("contig","position","refAllele","altAllele"), all.x=T)
    AS_GT[is.na(AS_GT)] <- "0,0"

    # Join GT with AS:
    u <- unique(AS$sample)
    for(i in u){
        col <- paste0(i,"_GT")
        AS_GT[,(col):=paste0(get(col),":",get(i))]
    }
    #remove cols
    AS_GT[,c("variantID",u):=NULL]

    # add info score
    info <- info[info>=0.4][,REF:=sub("_.*","",(sub(":.*","",gsub(".*[0-9].","",rs_id))))][,ALT:=sub(".*_","",(sub(".*:","",gsub(".*[0-9].","",rs_id))))]
    info[,position:=as.integer(position)]

#merge AS_GT with info

    AS_GT_i <- merge(AS_GT,info[,.(position,REF,ALT,info)],by.x=c('POS','REF','ALT'), by.y=c('position','REF','ALT'), all.x=T)

    AS_GT_i[,GT:="GT:AS"]
    AS_GT_i[,info:=paste0("RSQ=",info)]

    setcolorder(AS_GT_i, c(names(t_file)[1:7],"info","GT",names(t_file)[9:ncol(t_file)]))

# add ID

    AS_GT_i[,ID:=paste0(CHROM,":",POS)]


#save as tab delim

    write.table(AS_GT_i, file=paste0(path,"/",prefix,".",chr,".for.AS.tab"),row.names=F, col.names=F, quote=F, sep="\t")

}

######## Input for RASQUAL using https://github.com/kauralasoo/rasqual/tree/master/rasqualTools

#' fsnp count per gene for rasqual input
#'
#' Takes snp coordinates to locate then withing genes
#' @param x path to input files with snp coordinates, output from bash functions vcf4rasqual or vcf4rasqualg. Default is current directory
#' @param pat pattern to match input files 
#' @param exons full path with file name containing exon bounderies for each gene, output from rasqualTools (script input_4_rasqual.R)
#' @param cis_window length of window to count snps, dafults to 5e5
#' @param out_path path to output files, defualt is current directory.
#' @param prefix vector with prefix names for output files (automatically added .txt)
#' @keywords format rasqual output
#' @export
#' @return saves file with snp_counts per gene
#' snp_coord

snp_coord <-  function(x=".",pat,exons, cis_window=5e5, out_path=".", prefix){

    files <- grep("header", list.files(path=x,pattern=pat, full.names=T), invert=T, value=T)
    coords <- lapply(files,fread, header=F, col.names=c("chr","pos","snp_id"))
    exons_all <- fread(file=exons)
    snp_counts <- lapply(coords, function(i) countSnpsOverlapingExons(exons_all,i,cis_window))
    snp_select <- lapply(snp_counts, function(i) dplyr::select(i, chromosome_name,gene_id, feature_snp_count, cis_snp_count))
    snp_count <- lapply(seq_along(snp_select), function(i) snp_select[[i]][snp_select[[i]]$chromosome_name==unique(coords[[i]]$chr),])

    lapply(seq_along(snp_count), function(i) write.table(snp_count[[i]],file=paste0(out_path,"/",prefix[i],".txt"), row.names=F, quote=F))
                         
}



#################### RASQUAL OUTPUT ######################
#' format rasqual output
#'
#' Selects and format rasqual output
#' @param x path to file with rasqual output
#' @param top.hits whether to limit results to strongest associated snp per gene or not
#' @keywords format rasqual output
#' @export
#' @return data table with rasqual formatted output
#' format_rasqual

format_rasqual <-  function(x,top.hits=c("yes","no")) {
    temp <- fread(x)
    names(temp) <- fread('/mrc-bsu/scratch/ev250/Cincinatti/quant/RASQUAL/output/names_rasqual.txt', header=F)$V1
    #convert effect size to fold change
    temp[,Effect_size:=Effect_size/(1-Effect_size)]
    setnames(temp,"Effect_size","Fold_change")
    #add p value
    temp[, p:=pchisq(chi_square, df = 1, lower = F)]
    # correct by MT (number of tested reg snps)
    temp[,p_adj:=p*r_SNPs][p_adj>1, p_adj:=1]
    #select for each gene the SNP with min p_adj, also some entries have 0 tested_SNPs, remove:
    temp <- temp[r_SNPs!=0, ]
    setkey(temp, gene_id,p_adj)
    top.hits <- match.arg(top.hits)
    if(top.hits=="yes") {

    DT <- temp[,.SD[1], by=gene_id]
    setkey(DT,p_adj)
    DT[, FDR:=p.adjust(p_adj,method="fdr")]
    return(DT)
    } else {
	return(temp)
}

}

#' format rasqual output genome wide
#'
#' Selects and format rasqual output for many files
#' @param x path to files with rasqual output
#' @param pattern pattern to match input files
#' @param top.hits whether to limit results to strongest associated snp per gene or not
#' @param failed whether to return sub-list with rasqual output that failed the run, dafults to TRUE
#' @keywords format rasqual output
#' @export
#' @return list with data tables with rasqual formatted output per input file, if failed=TRUE also a nested list with failed rasqual runs per chr (gene_id, Chrom, number of fSNPs and number of regSNPs
#' format_rasqual_gw

format_rasqual_gw <-  function(x,pattern ,top.hits=c("yes","no"), failed=TRUE){
    setwd(x)
    files <- list.files(pattern=pattern, full.names=T)
    tmp <- lapply(files, function(i) format_rasqual(i,top.hits))
    tmp2 <- lapply(tmp, function(i) i[rs_id !="SKIPPED",])
    if(failed==TRUE){
        f <- rbindlist(lapply(tmp, function(i) i[rs_id=="SKIPPED",.(gene_id,Chrom,f_SNPs,r_SNPs)]))
        tmp3 <- list(successful_snps=tmp2,failed_snps=f)
        return(tmp3)
    } else {
        return(tmp2)
  
}
}


#' rasqual merge dna and rna by FDR cut-off, add gene name and variant name from biomart
#'
#' Merges DNA and RNA genotyping results, select results based on FDR cut-off and adds gene name and rsid for variants according to biomart
#' @param rasq_top_hit_dna rasqual output for dna genotyped, output from format_rasqual
#' @param rasq_top_hit_rna rasqual output for rna genotyped, output from format_rasqual
#' @param fdr false discovery rate cut-off, defaults to 0.1
#' @param cols vector with columns to select from rasq_top_hit objects, defaults to all
#' @param shape whether long or wide format is preferred, long format is the input for ras.table function
#' @keywords compare_rasqual_DNA_RNA 
#' @export
#' @return data table 
#' merge_format

merge_format <- function(rasq_top_hit_dna,rasq_top_hit_rna,fdr=0.1,cols,shape=c("long","wide")){
    rasq_top_hit_dna[, Genotype:="DNA"]
    rasq_top_hit_rna[,Genotype:="RNA"]
    if(missing(cols)==FALSE){
    hits <- merge(rasq_top_hit_dna[FDR<=fdr,cols, with=FALSE],rasq_top_hit_rna[FDR<=fdr,cols,with=F], by="gene_id", suffixes=c(".DNA",".RNA"),all=TRUE)
    } else {
         hits <- merge(rasq_top_hit_dna[FDR<=fdr,],rasq_top_hit_rna[FDR<=fdr,], by="gene_id", suffixes=c(".DNA",".RNA"),all=TRUE)
    }
    if(dim(hits)[1]==0){
        return(paste0("No hits at ",fdr," FDR, neither in DNA nor in RNA datasets for CHR ",unique(rasq_top_hit_rna$Chrom)))
        }
        # Add ensembl data
    mart <- useEnsembl(biomart="ensembl", GRCh="37",dataset="hsapiens_gene_ensembl")
    snp <-useEnsembl(biomart="snp", GRCh="37", dataset="hsapiens_snp")

    gene_name <- data.table(getBM(attributes=c("ensembl_gene_id", "external_gene_name"), filters=c("ensembl_gene_id"), values=hits$gene_id, mart=mart))
    
    if(sum(!is.na(hits$SNP_pos.DNA))==0){
        
        snp_pos <-data.table(pos=c(hits[!is.na(SNP_pos.DNA),SNP_pos.DNA], hits[!is.na(SNP_pos.RNA),SNP_pos.RNA]), chr=sub(":.*","", c(hits[!is.na(rs_id.DNA),rs_id.DNA], hits[!is.na(rs_id.RNA),rs_id.RNA])), allele=paste0(hits[!is.na(Ref.RNA),Ref.RNA], "/" ,hits[!is.na(Alt.RNA),Alt.RNA]))
        
    } else if(sum(!is.na(hits$SNP_pos.RNA))==0){
        snp_pos <- data.table(pos=c(hits[!is.na(SNP_pos.DNA),SNP_pos.DNA], hits[!is.na(SNP_pos.RNA),SNP_pos.RNA]), chr=sub(":.*","", c(hits[!is.na(rs_id.DNA),rs_id.DNA], hits[!is.na(rs_id.RNA),rs_id.RNA])), allele=paste0(hits[!is.na(Ref.DNA),Ref.DNA], "/" ,hits[!is.na(Alt.DNA),Alt.DNA]))

        } else {
    
            snp_pos <- data.table(pos=c(hits[!is.na(SNP_pos.DNA),SNP_pos.DNA], hits[!is.na(SNP_pos.RNA),SNP_pos.RNA]), chr=sub(":.*","", c(hits[!is.na(rs_id.DNA),rs_id.DNA], hits[!is.na(rs_id.RNA),rs_id.RNA])), allele=c(paste0(hits[!is.na(Ref.DNA),Ref.DNA], "/" ,hits[!is.na(Alt.DNA),Alt.DNA]),paste0(hits[!is.na(Ref.RNA),Ref.RNA], "/" ,hits[!is.na(Alt.RNA),Alt.RNA])))

        }
 
    
    snp_name <- rbindlist(lapply(1:nrow(snp_pos), function(i) data.table(getBM(attributes=c("refsnp_id", "chrom_start","allele"), filters=c("start","end","chr_name"), values=list(snp_pos[i,pos],snp_pos[i,pos],snp_pos[i,chr]),  mart=snp))))
    #check if same alleles in biomart and my data
    snp_pos[,pos:=as.integer(pos)]
    snp_name[,chrom_start:=as.integer(chrom_start)]

    snp_check <- merge(snp_pos,snp_name, by.x="pos", by.y="chrom_start", all.x=T)
    snp_check[,same:=allele.x==allele.y]
    #same_check <- snp_check[pos %in% hits[!is.na(SNP_pos.DNA),SNP_pos.DNA],same]
    for(i in hits[!is.na(SNP_pos.DNA),SNP_pos.DNA]){
        if(length(snp_check[pos==i & same==TRUE,same])!=0){
            if(unique(snp_check[pos==i & same==TRUE,same])==TRUE){
                hits[SNP_pos.DNA==i, rs_id.DNA:=unique(snp_check[pos==i & same==TRUE,refsnp_id])]
            }
        }
        
    }
    for(i in hits[!is.na(SNP_pos.RNA),SNP_pos.RNA]){
         if(length(snp_check[pos==i & same==TRUE,same])!=0){
        if(unique(snp_check[pos==i & same==TRUE, same])==TRUE){
            hits[SNP_pos.RNA==i, rs_id.RNA:=unique(snp_check[pos==i & same==TRUE,refsnp_id])]
        }
         }
    }
    
    hits <- merge(gene_name,hits, by.x="ensembl_gene_id", by.y="gene_id")

    setkey(hits, rs_id.DNA, rs_id.RNA)

    if(shape=="long"){

        hits_long <- melt(hits, id.vars=c( "ensembl_gene_id", "external_gene_name"), measure.vars=patterns(paste0("^",names(rasq_top_hit_dna)[2:(ncol(rasq_top_hit_dna)-1) ], "\\.")), variable.name="Genotype", value.name=names(rasq_top_hit_dna)[2:(ncol(rasq_top_hit_dna)-1) ])

        setkey(hits_long,ensembl_gene_id)
        hits_long[Genotype==1, Genotype:="DNA"][Genotype==2, Genotype:="RNA"]
        return(hits_long)
    } else {
        return(hits)

    }
}

#' rasqual merge dna and rna GENOME WIDE by FDR cut-off, add gene name and variant name from biomart
#'
#' Merges DNA and RNA genotyping results, select results based on FDR cut-off and adds gene name and rsid for variants according to biomart
#' @param rasq_top_hit_dna list with rasqual output for dna genotyped, output from format_rasqual_gw
#' @param rasq_top_hit_rna list with rasqual output for rna genotyped, output from format_rasqual_gw
#' @param fdr false discovery rate cut-off, defaults to 0.1
#' @param cols vector with columns to select from rasq_top_hit objects, default to all
#' @param shape whether long or wide format is preferred, long format is the input for ras.table function
#' @keywords compare_rasqual_DNA_RNA 
#' @export
#' @return list
#' merge_format_gw

merge_format_gw<- function(DNA_top_hits,RNA_top_hits,fdr=0.1,cols=NULL,shape=c("long","wide")){
    if(is.null(cols)){
        cols <- names(DNA_top_hits[[1]])
        }
    tmp <- lapply(seq_along(DNA_top_hits), function(i) merge_format(DNA_top_hits[[i]],RNA_top_hits[[i]], fdr, cols, shape))

    #tmp <- list()
    #for(i in seq_along(DNA_top_hits)){
       # tmp[[i]] <-  merge_format(DNA_top_hits[[i]],RNA_top_hits[[i]], fdr, cols, shape)
        #print(i)
        #}   
    return(tmp)
}



#' rasqual: complete rna or dna hits with dna or rna data respectively
#'
#' info for snps not significant in dna but sig in rna and viceversa. 
#' @param rasq_hits_dna rasqual output for dna genotyped, output from format_rasqual
#' @param rasq_hits_rna rasqual output for rna genotyped, output from format_rasqual
#' @param hits output from merge_format with shape="long"
#' @param shape whether long or wide format is preferred, In wide format can only return the matches but not in the same table as the originals.
#' @keywords compare_rasqual_DNA_RNA 
#' @export
#' @return data table 
#' complete_matches

complete_matches <- function(rasq_hits_dna, rasq_hits_rna,hits, shape="long"){
    h <- copy(hits)
    if(is.null(dim(h))){
        return(h)
    }
    
    dna <- aux("DNA", shape, h,rasq_hits_dna,rasq_hits_rna)
    rna <- aux("RNA", shape, h,rasq_hits_dna,rasq_hits_rna)
    temp <- rbind(dna,rna)
   
    if(shape!="long"){
        return(temp)
    } else {
        setcolorder(temp,names(temp)[c(1,ncol(temp),2:(ncol(temp)-1))])
         temp[,FDR:=NA]
        for (i in 1:nrow(temp)) {
            if(is.na(h[ensembl_gene_id==temp$gene_id[i] & Genotype==temp$Genotype[i],rs_id])[1]){
                h[ensembl_gene_id==temp$gene_id[i] & Genotype==temp$Genotype[i],
                     names(h)[4:ncol(h)]:=as.list(temp[i,3:ncol(temp),with=F])]
            } else {
                h <- rbind(h,merge(unique(h[,1:2]), temp[i,], by.x="ensembl_gene_id", by.y="gene_id", all.y=T))
            
            }
            
    }
        setkey(h, ensembl_gene_id)
        h[,SNP:=paste(Chrom,SNP_pos,Ref,Alt, sep=":")]
        return(h)
    }
    
}

#' subfunction for rasqual: complete rna or dna hits with dna or rna data respectively
#'
#' info for snps not significant in dna but sig in rna and viceversa. 
#' @param m whether genotyping is DNA or RNA based, default DNA
#' @param shape input shape of h 
#' @param h, output from merge_format, dafault to  shape ="long"
#' @param rasq_hits_dna rasqual output for dna genotyped, output from format_rasqual
#' @param rasq_hits_rna rasqual output for rna genotyped, output from format_rasqual
#' @keywords compare_rasqual_DNA_RNA 
#' @export
#' @return data table 
#' aux

aux <- function(m="DNA", shape, h,rasq_hits_dna,rasq_hits_rna){
    o <- tolower(m)
    if(m=="DNA"){
        n="RNA"
    } else {
        n="DNA"
    }
     temp <- get(paste0("rasq_hits_",o))
    if(shape!="long"){
    setkeyv(h, paste0("rs_id.",m))
    
    to_match <- data.table(rs_id=paste0(unlist(h[which(is.na(h[,paste0("rs_id.",m), with=F])),paste0("Chrom.",n), with=F]), ":", unlist(h[which(is.na(h[,paste0("rs_id.",m), with=F])),paste0("SNP_pos.", n), with=F])),
                           gene_id=unlist(h[which(is.na(h[,paste0("rs_id.",m), with=F])), ensembl_gene_id]))
    
   
    l <- rbindlist(lapply(1:nrow(to_match), function(i) temp[gene_id==to_match$gene_id[i] & rs_id==to_match$rs_id[i], ]))
    l[,Genotype:=m]
    
    return(l)
    } else {
        #print(names(h))
        temp1 <- h[Genotype==n & !is.na(rs_id),]
        temp1[,rs.id:=paste0(Chrom,":",SNP_pos)]
        l <- rbindlist(lapply(1:nrow(temp1), function(i) temp[gene_id==temp1$ensembl_gene_id[i] & rs_id==temp1$rs.id[i],]))
        l[,Genotype:=m]
        return(l)
    }
  
  } 



#' rasqual: complete rna or dna hits with dna or rna data respectively, GENOME WIDE
#'
#' info for snps not significant in dna but sig in rna and viceversa. 
#' @param DNA_hits list of rasqual output for dna genotyped, output from format_rasqual_gw
#' @param RNA_hits rasqual output for rna genotyped, output from format_rasqual_gw
#' @param Hits list output from merge_format_gw with shape="long" (using DNA_top and RNA_top as inputs)
#' @param shape whether long or wide format is preferred, In wide format can only return the matches but not in the same table as the originals.
#' @keywords compare_rasqual_DNA_RNA 
#' @export
#' @return list of data tables 
#' complete_matches_gw

complete_matches_gw<- function(DNA, RNA,Hits, shape="long"){
    tmp <- lapply(seq_along(DNA[[1]]), function(i) complete_matches(DNA[[1]][[i]],RNA[[1]][[i]],Hits[[i]],shape))
    #exclude missing data tables from tmp
    tmp2 <- tmp[-which(sapply(tmp,is.character)==TRUE)]
    #identify rasqual failed runs
    tmp2 <- complete_matches_aux(DNA,RNA,tmp2)
    return(tmp)
}

#' rasqual: aux for complete rna or dna hits with dna or rna data respectively, GENOME WIDE
#'
#' Indicates whether the rasqual run failed, due to rasqual unable to process large number of SNPs
#' @param DNA list of rasqual output for dna genotypes, output from format_rasqual_gw
#' @param RNA rasqual output for rna genotypes, output from format_rasqual_gw
#' @param tmp2 object from complete_matches_gw to format
#' @keywords compare_rasqual_DNA_RNA 
#' @export
#' @return tmp2 list indicating whether rasqual failed in snp_id column 
#' complete_matches_aux

complete_matches_aux<- function(DNA, RNA,tmp2){
    geno <- c("DNA","RNA")
    for(i in seq_along(tmp2)){
        for(j in geno){
             f <- which(tmp2[[i]][Genotype==j ,ensembl_gene_id] %in% get(j)[[2]][,gene_id])
             m <- tmp2[[i]][Genotype==j,ensembl_gene_id][f]
             if(length(m)==0){
                 next}
             tmp2[[i]][ensembl_gene_id %in% m & Genotype==j,rs_id:="rasqual.failed"]
        }
    }
    return(tmp2)
}

        

#' table with rasqual output
#'
#' make xtable to present rasqual output
#' @param hits_f output from complete_matches with shape="long"
#' @param cols name of the columns of hits_f to include in the table
#' @param print whether to print the result or just return a formatted data table input
#' @param shape whether longtable or sidewaystable
#' @keywords compare_rasqual_DNA_RNA 
#' @export
#' @return formatted data table optional to print a xtable
#' ras.table

ras.table <- function(hits_f, cols,print=TRUE,shape){
    h <- copy(hits_f)
    h <- h[,cols,with=F]
    s <- data.table(class=sapply(h,class))
    s[class=="character" | class=="factor", type:="s"][class=="integer",type:="d"][class=="numeric",type:="g"]
    names(h)[which(names(h)=="external_gene_name")] <- "Gene"
    names(h) <- gsub("Allele_freq", "AF", names(h))
    names(h) <- gsub("_quality", "_q", names(h))
    if("Squared_correlation_prior_posterior_rSNP" %in% cols){
        names(h) <- gsub("Squared_correlation_prior_posterior_rSNP","r2.pr.pst.rSNP",names(h))
        }
   
    if(print==TRUE){
        #formatting numbers as text and duplicating columns to keep numeric value
        SDcols=names(h)[which(s$type=="g")]
        SDcols2=names(h)[which(s$type=="d")]
        h[,(paste0(SDcols,".f")):=lapply(.SD,function(i) formatC(i,digits=2,format="g")), .SDcols=SDcols]
        h[,(paste0(SDcols2,".f")):=lapply(.SD,function(i) formatC(i,digits=2,format="d")), .SDcols=SDcols2]

        tof <-  names(h)[-which(s$type %in% c("g","d"))]
        
        u <- h[,.N,by="Gene"]
        for(i in u$Gene){
            if (u[Gene==i,N]==2){
                d <- which(duplicated(h[Gene==i,SNP]))
                if(length(d)==1){
                    if(sum(is.na(h[Gene==i,FDR]))==0){
                        h[Gene==i , (tof):=lapply(.SD, function(j) paste0("\\textcolor{red}{\\textbf{",j,"}}")), .SDcols=tof]
                        } else {
                            h[Gene==i & FDR<0.1, (tof):=lapply(.SD, function(j) paste0("\\textbf{",j,"}")), .SDcols=tof]
                        }
                    }
                if(length(d)==0) {
                     if(sum(is.na(h[Gene==i,FDR]))==0){
                         h[Gene==i , (tof):=lapply(.SD, function(j) paste0("\\textcolor{blue}{\\textbf{",j,"}}")), .SDcols=tof]
                     } else {
                         h[Gene==i & FDR<0.1, (tof):=lapply(.SD, function(j) paste0("\\textbf{",j,"}")), .SDcols=tof]
                         
                    }
                }
            }
            
            
            if(u[Gene==i,N] %in% 3:4){
                d <- which(duplicated(h[Gene==i,SNP]))
                for(k in 1:length(d)){
                    s <- h[Gene==i,SNP][d[k]]
                    if(sum(is.na(h[Gene==i & SNP==s,FDR]))==0){
                        h[Gene==i & SNP==s, (tof):=lapply(.SD, function(j) paste0("\\textcolor{red}{\\textbf{",j,"}}")), .SDcols=tof]
                    } else {
                        h[Gene==i & FDR<0.1, (tof):=lapply(.SD, function(j) paste0("\\textcolor{blue}{\\textbf{",j,"}}")), .SDcols=tof]
                   

                    }
                }
            }
        }
        h[,(SDcols):=NULL][,(SDcols2):=NULL]
        names(h) <- gsub(".f","",names(h))
        names(h) <- gsub("_",".",names(h))
        #t <- xtable(h)
        if(shape=="side"){print(xtable(h, NA.string="-",floating = TRUE, floating.environment = "sidewaystable", booktabs=TRUE, include.rownames=F))
      
              }
        if(shape=="long"){
            print(xtable(h), tabular.environment="longtable",size="small",floating=FALSE, include.rownames=FALSE, rotate.colnames = TRUE, sanitize.text.function=function(x) x)
    } else {
        return(h)
    }
}
}
#print(xtable(h, NA.string="-",floating = TRUE, floating.environment = "sidewaystable", booktabs=TRUE, include.rownames=F)




    


#' table with rasqual output genome wide
#'
#' make xtable to present rasqual output genome wide
#' @param x list output from complete_matches_gw with shape="long"
#' @param cols name of the columns of x[[1]] to include in the table
#' @param print whether to get xtable output
#' @param shape whether to return long or sideways tables
#' @keywords compare_rasqual_DNA_RNA 
#' @export
#' @return xtable
#' ras.table_gw

ras.table_gw <- function(x,cols,print=TRUE,shape=c("long","side")) {
    #remove null elements (no hits at FDR for RNA or DNA)
    x <- x[-which(lapply(lapply(x,nrow),is.null)==TRUE)]
    lapply(x, function(i) ras.table(i,cols,print,shape))
    #return(tmp)
}

#' summary of rasqual output genome wide
#'
#' summarize rasqual output in list, first element is a table with # of eQTL uniquely called in RNA or DNA and second element list of genes called in both
#' @param x list output from complete_matches_gw with shape="long"
#' @keywords compare_rasqual_DNA_RNA 
#' @export
#' @return list with first element table of eQTL uniquely called in RNA or DNA and second element
#' ras.summ_gw

ras.summ_gw <- function(x){
     #remove null elements (no hits at FDR for RNA or DNA)
    x1 <- x[-which(lapply(lapply(x,nrow),is.null)==TRUE)]
    x1 <- rbindlist(x1)
    x2 <- x1[!is.na(FDR),]
    x3 <- x2[,.N, by=external_gene_name]
    x4 <- table(as.character(x2[external_gene_name %in% x3[N==1,external_gene_name],Genotype]))
   return(list(eQTL.unique.DNAorRNA=x4, eQTL.both.DNA.RNA=x3[N==2,external_gene_name]))
}

#' save rasqual output genome wide
#'
#' save
#' @param x list output from complete_matches_gw with shape="long"
#' @param file path to file
#' @keywords compare_rasqual_DNA_RNA 
#' @export
#' @return save csv
#' ras.save_gw

ras.save_gw <- function(x,file){
     #remove null elements (no hits at FDR for RNA or DNA)
    x1 <- x[-which(lapply(lapply(x,nrow),is.null)==TRUE)]
    x1 <- rbindlist(x1)
    write.csv(x1,file=file, row.names=F)
}






    
###### QC of quantification using Kallisto ##########

#' QC of quantification using Kallisto: process meta data with run variable
#'
#' Add run variable to meta-data Cincinatti files and format according to 
#' @param meta_data_file file with meta data
#' @keywords format meta data
#' @export
#' @return data table with formatted meta data
#' f_path

f_path <- function(meta_data_file,path_q_files){
  samp_info <- fread(meta_data_file)
  #'/scratch/wallace/FastQ_Cincinatti/samples-rnaseq.csv'
  # add run column, as indicator variable for batch effect
  samp_info[, run:=gsub("_.*","",`RNA Seq (Fast Q) Cincinatti File Name`)]
  # needs the column with the sample name to be named 'sample'
  setnames(samp_info,c("CHARMS P_ID","CHARMS cohort", "Gender M/F", "JIA ILAR Subtype", "ACR Response", "Caucasian?"),c("sample","CHARMS_cohort","GENDER", "JIA_ILAR_Subtype","ACR_response","Caucasian"))
  # add path to quant files
  samp_info[, path:=paste0(path_q_files,sample)]
  return(samp_info)
}




   



#################  QC STAR #######################

#' subfunction to picard or rseqc
#'
#' open files with output from picard or rseqc
#' @param files vector with filenames to open
#' @param tool whether picard or rseqc output is to be analysed
#' @param whether to include a column with sample names in the output 
#' @keywords format meta data
#' @export
#' @return data table with or without a column with sample id
#' read.files

read.files <- function(files, tool=c("picard","rseqc"), idcol=NULL){
    if(tool=="picard"){
    all <- lapply(files, function(i) fread(paste0('grep -A2 PF_ ', i), colClasses=rep("numeric",25)))
    names(all) <- sapply(files, function(i) gsub(".*/","", gsub(".RNA.*","",i)))
    } else {
        all <- lapply(files, function(i) fread(paste0('grep -A10 Group ', i))) #, colClasses=rep("numeric",25)))
        names(all) <- sapply(files, function(i) gsub(".*/","",gsub(".txt","",i)))
        }
    if(is.null(idcol)){
    all_noid <- rbindlist(all)
    return(all_noid)
    } else {
        all <- rbindlist(all, idcol="sample")
        return(all)
        }
    
}



#' Format picard output results
#'
#' open files with output from picard, merges into a datatable
#' @param files vector with filenames to open
#' @param xtable whether to print the data table as xtable
#' @keywords picard QC
#' @export
#' @return data table with or without xtable
#' picard

picard <- function(files,xtable=FALSE){
    tool="picard"
    all_noid <- read.files(files,tool)
    
    summ <- data.table(t(rbind(Mean <- all_noid[,lapply(.SD,mean)], stdev <- all_noid[,lapply(.SD,sd)])),keep.rownames=T)
    names(summ) <- c("Variable","Mean", "SD")
    summ <- summ[-c(8:10,18,23:25),]
    summ[Variable %in% grep("PCT_",Variable,value=T),':=' (Mean=100*Mean, SD=100*SD)]
    #summ[, ':=' (Mean=format(Mean,digits=2, width=2), SD=format(SD,digits=2, width=2))]

    if(xtable==TRUE) {
    print(xtable(summ, caption="BAM QC using PICARD RNA metrics. PF: passing filter, PCT: percentage."), booktabs=T,include.rownames = FALSE)
    }
    return(summ)
}

#' Format multiqc output
#'
#' open file with output from multiqc, merges into a datatable
#' @param file path with file name to multiqc output
#' @param xtable whether to print the data table as xtable
#' @keywords multiqc
#' @export
#' @return data table with or without xtable
#' multiqc

multiqc <- function(file, xtable=FALSE) {

    multi <- fread(file)
    multi_sum <- data.table(t(rbind(Mean <- multi[,lapply(.SD,mean),.SDcols=names(multi)[2:3]], stdev <-  multi[,lapply(.SD,sd),.SDcols=names(multi)[2:3]])), keep.rownames=T)
    names(multi_sum) <- c("Reads","Mean", "SD")
    multi_sum[, ':=' (Mean=format(Mean,digits=2, width=2), SD=format(SD,digits=2, width=2))]

    multi_stat <- fread('/scratch/ev250/Cincinatti/quant/STAR/QC/multiqc_data/multiqc_data/multiqc_star.txt')
    multi_stat_sum <- data.table(t(rbind(Mean <- multi_stat[,lapply(.SD,mean),.SDcols=names(multi_stat)[-1]], stdev <-  multi_stat[,lapply(.SD,sd),.SDcols=names(multi_stat)[-1]])), keep.rownames=T)
    names(multi_stat_sum) <- c("Reads","Mean", "SD")
    multi_stat_sum[, ':=' (Mean=format(Mean,digits=2, width=2), SD=format(SD,digits=2, width=2))]
    
    if(xtable==TRUE){
    print(xtable(multi_stat_sum, caption="BAM QC using multiqc."), booktabs=T,include.rownames = FALSE)
    }

    return(multi_stat_sum)

}

#' Format rseqc output results
#'
#' open files with output from picard, merges into a datatable
#' @param files vector with full path and filenames to open
#' @param xtable whether to print the data table as xtable
#' @keywords rseqc
#' @export
#' @return data table with or without xtable
#' rseqc

rseqc <- function(files, xtable=FALSE){
    tool <- "rseqc"
    all_noid <- read.files(files,tool)

    summ <- merge(all_noid[,lapply(.SD,mean), by=Group],  all_noid[,lapply(.SD,sd),by=Group] , by="Group", suffixes=c("_Mean", "_SD"))
    summ[,Total_bases_SD:=NULL]
    setnames(summ,"Total_bases_Mean", "Total_bases")
    setcolorder(summ,c("Group","Total_bases","Tag_count_Mean", "Tag_count_SD",  "Tags/Kb_Mean", "Tags/Kb_SD"))
    if(xtable==TRUE){

print(xtable(summ, caption="BAM QC using RSeQC. Total bases are the number of bases covered in the analysis. Tags are the number of reads except when reads are spliced that are counted as the number of splice events plus one."), booktabs=T, include.rownames=FALSE,size="small",floating = TRUE, floating.environment = "sidewaystable")
    }
    return(summ)

}


#' Combine PICARD with RSeQC to get better idea of intergenic coverage
#'
#' merges output from picard and rseqc
#' @param files.p vector with full names to picard output
#' @param files.r vector with full names to rseqc output
#' @param xtable whether to print the data table as xtable
#' @keywords rseqc picard
#' @export
#' @return data table with or without xtable
#' compbine_qc

combine_qc <- function(files.p,files.r,xtable=TRUE){
  
    # get raw data and intergenic/coding ratio
    pic <- read.files(files.p, tool="picard", idcol=TRUE)
    pic[,int2cod:= INTERGENIC_BASES/CODING_BASES]
    rseq <- read.files(files.r, tool="rseqc", idcol=TRUE)
    rseq <- rseq[sample %in% pic$sample,.(sample,Group,Tag_count)]
    rseq <- data.table(spread(rseq,Group,Tag_count))
    rseq[,tes102cod:=(TES_down_10kb + TSS_up_10kb)/CDS_Exons]
    # merge    
    temp <- merge(pic[,c(1:8,27), with=F], rseq[,c(1,12), with=F], by="sample")
    temp[,inter.far:=(int2cod-tes102cod)*CODING_BASES]
    temp[,inter.close:=tes102cod*CODING_BASES]
    temp[,c("sample","int2cod","tes102cod","INTERGENIC_BASES"):=NULL]
    temp[,paste0(names(temp)[1:ncol(temp)], "_PCT"):=lapply(names(temp)[1:ncol(temp)],function(i) get(i)*100/PF_ALIGNED_BASES)]

    summ <- data.table(t(rbind(Mean <- temp[,lapply(.SD,mean)], stdev <- temp[,lapply(.SD,sd)])),keep.rownames=T)
    summ_w <- cbind(summ[3:8,],summ[11:16,.(V1,V2)])
    names(summ_w) <- paste0("V",1:ncol(summ_w))
    summ_w[,V1:=c("Ribosomal","Coding", "UTR","Intronic","Intergenic 10KB+", "Intergenic 10KB")]
    summ_w <- summ_w[order(-V2),]
    summ_w[,paste0("V",4:5) := lapply(paste0("V",4:5), function(i) round(get(i),2))  ]
    
    summ_w[,paste0("V",2:3) := lapply(paste0("V",2:3), function(i) formatC(get(i),digits=2, format="e"))  ]
    names(summ_w) <- c("Region","Mean_bases","SD_bases", "Mean_%", "SD_%")
    if(xtable==TRUE){

        summ_w[,`Aligned bases (SD)`:=paste0(Mean_bases," (",SD_bases,")")]
        summ_w[,`Aligned bases % (SD)`:=paste0(`Mean_%`," (", `SD_%`,")")] 
        summ_w[, c("Mean_bases", "SD_bases", "Mean_%", "SD_%"):=NULL]
        addtorow <- list()
        addtorow$pos <- list(0,0)
        addtorow$command <- c(" &\\multicolumn{2}{c}{Aligned bases}\\\\\n", "Region & Number (SD) & Percentage (SD) \\\\\n")
        
        print(xtable(summ_w), booktabs=T, include.rownames=F,add.to.row=addtorow, include.colnames=FALSE )
    }
 return(summ_w)
}
