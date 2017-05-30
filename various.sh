#!/bin/bash


# split gen file per chr for impt2 input

input_chr () {
    # arg 1 path to gen file
    for i in `seq 22`; do
       out=$(echo $1 | sed "s/gen//")
	grep ^$i: $1  > $out$i.gen
    done
}



# get imputation range for each chromosome.
#arg 1, path to vcf file with all chrs
#arg2, paht to output file

imp_range () {
   for i in `seq 22`; do

       st=$(grep ^$i $1 | cut -d' ' -f3 | head -1)
       end=$(grep ^$i $1 | cut -d' ' -f3 | tail -1)
       inc=$(((end-st + 5000000-1)/5000000))
       echo $st $end $((st + inc*5000000))
      
   done > $2
   
 }

#concatenate files after imputation

#arg1 input path 
#arg2 output path
conc_imp2 () {
    cd $1
    for i in `seq 22`; do
	prefix=$(ls *.$i.phased*imp2 |sed "s/phased.*//" | uniq)
	ls *.$i.phased*imp2 | xargs cat  > $2/${prefix}imp2
	ls *.$i.phased*imp2_info | xargs cat  > $2/${prefix}imp2_info
	
done

}

## add genomic annotation to vcf file
# arg 1: path to files

genom_annot () {

    cd /home/ev250/bin/snpEff
    files=$(ls $1 | grep vcf | grep -v annot | grep -v gz)
    for file in $files; do
	out=$(echo $file | sed 's/\.vcf//')
	java -Xmx4g -jar snpEff.jar -v GRCh37.75 -canon -ud 10000 -t -noStats $file > $1/$out.annot.vcf

	done
	
    }

## extract fields out of vcf and convert to txt to read in R
#arg1: full path to files to be extracted
#arg2: prefix to distinguish files to extract
#arg3: list of fields to extract

vcf2R () {
    cd $1
    files=$(ls *$2*)
    for file in $files; do
	id=$(echo $file | sed 's/\.vcf//')
	bcftools query -f $3 $file > $1/${id}.txt
	bcftools query -f $3 $file -H | head -1 > $1/${id}.header.txt
	done
    }


## bgzip and tabix files
#arg1 dir containing files
#arg2 pattern to match files

vcr_index () {
    cd $1
    files=$(ls *$2)
    for file in $files; do
	bgzip -c $file > $file.gz
	tabix -p vcf $file.gz &
    done
    wait
}

#convert gen (imp2) format to vcf
#arg1 dir containing gen files
#arg2 pattern to match
#arg3 path to sample file

gen2vcf () {
    cd $1
    files=$(ls *$2)
    for file in $files; do
	out=$(echo $file | sed 's/.gen//')
	echo $file $out $1 $2 $3
	bcftools convert -G $file,$3 -Ov -o $out.vcf
    done
   
		
}


#format files to compute uncertanty for haplotype estimation by shapeit2

#arg1 input path 
#arg2 pattern to match
#arg3 prefix to save
hap_un () {
    cd $1
    files=$(ls | grep $2)
	ls $files | xargs cat  | awk '{print $1"_"$2"_"$3"_"$4"_"$5" "$0}' > $1/$3
}


# format vcf DNA or RNA file to add AS annotation

#arg1 path to files
#arg2 first chromosome to loop through
#arg3 last chromosome to loop through
#arg4 pattern in files to match, defaults to DNA
#arg5 path to output dir, defaults to DNA

vcf4AS () {
    cd $1
    chrs=$(seq $2 $3)
    for chr in $chrs; do
	if [ $# -eq 3 ]
	then
	      files=$(ls chr$chr.*.vcf.gz | grep -v ASE)
	else
	    files=$(ls $chr.$4)
	fi
	
	for file in $files; do
	    
	    id=$(echo $file | sed -e "s/chr$chr.//" -e 's/.vcf.gz//')
	    if [ $# -eq 3 ]
	    then
		
	    bcftools annotate -x ^FORMAT/GT $file | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t GT[\t%GT]\n' -o $id.chr$chr.less.tab

	    #header
	    bcftools annotate -x ^FORMAT/GT $file | bcftools view -h | sed -e '7i\
##FORMAT=<ID=AS,Number=2,Type=Integer,Description="Allele-specific expression counts from RNA-seq">' -e '/^##INFO=<ID/d' | sed '4i\
##INFO=<ID=RSQ,Number=1,Type=Float,Description="Genotype imputation quality from MaCH\/Thunder">' > $id.chr$chr.ASE.vcf

	    else
		 bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t GT[\t%GT]\n' -o $5/$id.rna.less.txt $file
		bcftools query -f  '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t GT[\t%GT]\n' $file -H | head -1 >$5/$id.rna.less.header.txt

		#header
		bcftools view -h $file | sed -e '7i\
##FORMAT=<ID=AS,Number=2,Type=Integer,Description="Allele-specific expression counts from RNA-seq">' -e '/^##INFO=<ID/d' | sed '5i\
##INFO=<ID=RSQ,Number=1,Type=Float,Description="Genotype imputation quality from MaCH\/Thunder">' > $5/$id.rna.ASE.vcf

	    fi	
	done
	
    done

}

#convert tab files with GT:AS and RSQ fields into vcf: one file per sample  (per chr)

#arg1 path to files
#arg2 pattern to match

tab2vcf () {
    cd $1
    files=$(ls | grep $2)
    for file in $files; do
	# needs to be converted into vcf: add header
	id=$(echo $file | sed "s/\\.$2//")
	cat $file >> $id.ASE.vcf
	#compress and tabix 
	bgzip $id.ASE.vcf
	bcftools index -t  $id.ASE.vcf.gz
    done

}


#alternative more general function

#arg1 path to files containing tab, vcf header
#arg2 full path with file name of a space separated file with column 1 listing tab files and column 2 the corresponding vcf header

tab2vcfg () {
    cd $1
    NL=$(cat $2 | wc -l)
    for i in `seq $NL`; do
	tab=$(awk "NR==$i" $2 | cut -d ' ' -f1)
	vcf=$(awk "NR==$i" $2 | cut -d ' ' -f2)
	cat $tab >> $vcf
	bgzip $vcf
	bcftools index -t $vcf.gz
    done
}

    


#merge vcfs into 1 for each chr

#arg1 path to files
#arg2 first chromosome to loop through
#arg3 last chromosome to loop through

mergevcf () {

    cd $1
    chrs=$(seq $2 $3)
    for chr in $chrs; do
	ls *chr$chr.ASE.vcf.gz > file.list
	bcftools merge -m none --threads 14 -l file.list -Oz -o chr$chr.ASE.allsamples.vcf.gz
	bcftools index -t  chr$chr.ASE.allsamples.vcf.gz
    done
}


#extract samples from multiple vcfs

#arg1 path to vcf files
#arg2 full path (if different from arg1) with file name for a space separated file with column 1 names of input vcf and col2 name of output vcf
#arg3 full path (if different from arg1) with file name to a list of samples to extract from vcf. One sample per line.

extSamplesVcf() {
    cd $1
    NL=$(cat $2 | wc -l)
    for i in `seq $NL`; do
	input=$(awk "NR==$i" $2 | cut -d ' ' -f1)
	output=$(awk "NR==$i" $2 | cut -d ' ' -f2)

	bcftools view -Oz -o $output.gz -S $3  $input.gz
	bcftools index -t $output.gz
    done
}


# extract info for rasqual
#arg1 path to files
#arg2 first chr to extract
#arg3 last chr to exract
#arg4 pattern in files to match
#arg5 pattern to extract, input for bcftools query
#arg6 path to output, optional if different to path to files
#arg7 suffix to file

vcf4rasqual () {
    cd $1
    chrs=$(seq $2 $3)
    if [ $# -eq 7 ]
    then
	var=$6
	var2=$7
    else
	var=$1
	var2=$6
    fi
    for chr in $chrs; do
	# "" to take argument ignoring blank spaces
	bcftools query -f "$5" chr$chr.$4 > $var/chr$chr.$var2
	bcftools query -f "$5" chr$chr.$4 -H | head -1 > $var/chr$chr.header.$var2
done
   
}

# extract info for rasqual general version
#arg1 path to files
#arg2 file with the list of vcf files to extract info from in column 1, names for corresponding output files in col2
#arg3 pattern to extract, input for bcftools query
#arg4 path to output, if different to path to files, optional

vcf4rasqualg () {
    cd $1
    if [ $# -eq 4 ]
    then
	var=$4
    else
	var=$1
    fi
    NL=$(cat $2 | wc -l)
    for i in `seq $NL`; do
	input=$(awk "NR==$i" $2 | cut -d ' ' -f1)
	output=$(awk "NR==$i" $2 | cut -d ' ' -f2)

	bcftools query -f "$3" $input > $var/$output
	bcftools query -f "$3" $input -H | head -1 > $var/header.$output
    done
}


