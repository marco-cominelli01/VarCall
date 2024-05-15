#!/bin/bash

###############################################################################################
# This script allows to do variant calling starting from the single-end short reads of a trio #
# child-mother-father.																		  #
# To do so, the type of the disease (either autosomic dominant -ad or recessive -ar) must     #
# be assumed and specified at the beginning.												  #
# The variant caller used is freebayes and is possible to specify the depth at which looking  #
# for variants: open "greppy" with a text editor to find out more.                            #
# At the end, all the necessary QCs are generated (aggregated by means of multiqc)            #
###############################################################################################
# How to run this script: nohup ./VarCall_pipeline.sh -ar case178 -ad case210 case221 case240 &       #
###############################################################################################
# This script was written by Marco Cominelli and Elena Sasso 								  #
# from Università degli Studi di Milano and Politecnico di Milano							  #
###############################################################################################

# You need:
# 1- An empty directory (see the first line of code)
# 2- This script (saved in the directory that contains the directory created at point 1)
# 3- The python script "greppy" (saved in the same directory of this script)

# Absolute path of the directory containing this script and "greppy" (must be the directory above the one in which all outputs will be saved)
script_dir=/home/BCG_2024_mcominelli

# Absolute path of already existing empty directory where to save all the outputs
my_wd=/home/BCG_2024_mcominelli/True_Project     
        
# Absolute path of the directory containing all .fq.gz, genome indexing, genome FASTA and BED files
file_dir=/home/BCG2024_genomics_exam

# Absolute path of the BED file for the intersect
bed_file=/home/BCG2024_genomics_exam/exons16Padded_sorted.bed

# Search depth for autosomic recessive cases (more info about this in 'greppy') 
ar_depth='basic'

# Search depth for autosomic dominant cases (more info about this in 'greppy') 
ad_depth='basic'     
          
create_directories() {
	mkdir "$my_wd/$1/"
	mkdir "$my_wd/$1/QC"
	mkdir "$my_wd/$1/bam/"
	mkdir "$my_wd/$1/bam/bamcov"
	mkdir "$my_wd/$1/vcf/"
}

quality_controls() {
	# FASTQC
	fastqc ${file_dir}/${1}* -o ${my_wd}/${1}/QC
	
	# BAMQC
	qualimap bamqc -bam ${1}_father.bam -gff ${bed_file} -outdir ./../QC/${1}_father
	qualimap bamqc -bam ${1}_child.bam -gff ${bed_file} -outdir ./../QC/${1}_child
	qualimap bamqc -bam ${1}_mother.bam -gff ${bed_file} -outdir ./../QC/${1}_mother

	# MULTIQC
	multiqc ./../QC/ --outdir ./../QC/multiqc
}

bowtie_alignment() {
	bowtie2 -U ${file_dir}/${1}_father.fq.gz -x ${file_dir}/uni --rg-id 'SF' --rg "SM:father" | samtools view -Sb | samtools sort -o ${1}_father.bam
	bowtie2 -U ${file_dir}/${1}_child.fq.gz -x ${file_dir}/uni --rg-id 'SC' --rg "SM:child" | samtools view -Sb | samtools sort -o ${1}_child.bam
    bowtie2 -U ${file_dir}/${1}_mother.fq.gz -x ${file_dir}/uni --rg-id 'SM' --rg "SM:mother" | samtools view -Sb | samtools sort -o ${1}_mother.bam
}

bam_indexing() {
	samtools index ${1}_father.bam 
    samtools index ${1}_child.bam
    samtools index ${1}_mother.bam
}

bam_coverage() {
	bedtools genomecov -ibam ../${1}_father.bam -bg -trackline -trackopts 'name="father"' -max 100 > ${1}_fatherCov.bg
	bedtools genomecov -ibam ../${1}_child.bam -bg -trackline -trackopts 'name="child"' -max 100 > ${1}_childCov.bg
	bedtools genomecov -ibam ../${1}_mother.bam -bg -trackline -trackopts 'name="mother"' -max 100 > ${1}_motherCov.bg
}

sorting() {
	bcftools query -l ${1}_filtered.vcf | sort > samples.txt     
    bcftools view -S samples.txt ${1}_filtered.vcf > ${1}_sorted.vcf 
}
		  
# Empty array to accomodate autosomic recessive cases
declare -a ar_cases=()

# Empty array to accomodate autosomic dominant cases  
declare -a ad_cases=()  

current_param=""

while [ $# -gt 0 ]; do
    case "$1" in
        -ad)
            current_param="-ad"
            shift
            ;;
        -ar)
            current_param="-ar"
            shift
            ;;
         -h|--help) 
	    echo ""
	    echo "                  VarCall_pipeline: Variant Calling Pipeline"
	    echo ""
    	    echo "   Built on 2024-04-22 05:59"
	    echo ""
	    echo "   Example of how to run it: ./VarCall_pipeline.sh -ar case510 case560 -ad case456 case438"
	    echo ""
	    echo "   -ar <arg(s)>                          Prefix of the autosomic recessive fastq file(s)"    
	    echo "   -ad <arg(s)>                          Prefix of the autosomic dominant fastq file(s)"
    	    echo "   -h --help                             Print this help file and exit"
	    exit 0
	    ;;
          *)
            if [ "$current_param" == "-ar" ]; then
                ar_cases+=("$1")
            elif [ "$current_param" == "-ad" ]; then
                ad_cases+=("$1") 
            fi
            shift
            ;;
    esac
done

echo "-ar vector contains the cases: ${ar_cases[@]}"
echo "-ad vector contains the cases: ${ad_cases[@]}"

# Complete analysis for autosomic recessive cases
for ar_trio in ${ar_cases[@]}
do
	# Creation of all the necessary directories
	create_directories $ar_trio
	
	cd ${my_wd}/${ar_trio}/bam

	# Alignment
	bowtie_alignment $ar_trio
	
	# FASTQC + BAMQC + MULTIQC
	quality_controls $ar_trio
	
	# BAM indexing
	bam_indexing $ar_trio
	
	cd ./bamcov

	# BAM coverage
	bam_coverage $ar_trio

	cd ./../../vcf

	# Variant calling
	freebayes -f ${file_dir}/universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 ../bam/${ar_trio}_mother.bam ../bam/${ar_trio}_child.bam ../bam/${ar_trio}_father.bam  > ${ar_trio}.vcf

 	# Variant prioritization
    grep "#" ${ar_trio}.vcf > ${ar_trio}_filtered.vcf
	cat ${ar_trio}.vcf | grep -v "#" | python ${script_dir}/greppy "ar" ${ar_depth} >> ${ar_trio}_filtered.vcf
	
	# Sorting columns by family members name
	sorting $ar_trio 

    # Intersection 
    bedtools intersect -a ${ar_trio}_sorted.vcf -b ${bed_file} -u > ${ar_trio}_final.vcf
	
done

# Complete analysis for autosomic dominant cases
for ad_trio in ${ad_cases[@]}
do
	# Creation of all the necessary directories
	create_directories $ad_trio

    cd ${my_wd}/${ad_trio}/bam
	
	# Alignment
    bowtie_alignment $ad_trio
	
	# FASTQC + BAMQC + MULTIQC
	quality_controls $ad_trio

	# BAM indexing
    bam_indexing $ad_trio

    cd ./bamcov
	
	# BAM coverage
    bam_coverage $ad_trio

    cd ./../../vcf

	# Variant calling
    freebayes -f ${file_dir}/universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 ../bam/${ad_trio}_mother.bam ../bam/${ad_trio}_child.bam ../bam/${ad_trio}_father.bam  > ${ad_trio}.vcf

    # Variant prioritization
    grep "#" ${ad_trio}.vcf > ${ad_trio}_filtered.vcf
    cat ${ad_trio}.vcf | grep -v "#" | python ${script_dir}/greppy "ad" ${ad_depth} >> ${ad_trio}_filtered.vcf

    # Sorting columns by family members name
    sorting $ad_trio 

    # Intersection
    bedtools intersect -a ${ad_trio}_sorted.vcf -b ${bed_file} -u > ${ad_trio}_final.vcf
	
done








