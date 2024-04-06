#!/bin/bash

###############################################################################################
# This script allows to do variant calling starting from the single-end short reads of a trio #
# child-mother-father.																		  #
# To do so, the type of the disease (either autosomic dominant -ad or recessive -ar) must     #
# be assumed and specified at the beginning.												  #
# The variant caller used is freebayes and is possible to specify the depth at which looking  #
# for variants: open "greppy.py" with a text editor to find out more.                         #
# At the end, all the necessary QCs are generated (fastqc, bamqc and multiqc)                 #
###############################################################################################
# How to run this script: nohup ./pipeline.sh -ar case178 -ad case210 case221 case240 &       #
###############################################################################################
# This script was written by Marco Cominelli and Elena Sasso 								  #
# from Università degli Studi di Milano and Politecnico di Milano							  #
###############################################################################################

# You need:
# 1- An empty directory (see the first line of code)
# 2- This script (saved in the directory that contains the directory created at point 1)
# 3- The python script "greppy.py" (saved in the same directory of this script)


# Absolute path of already existing empty directory where to save all the outputs
my_wd=/home/BCG_2024_mcominelli/project     
        
# Absolute path of the directory containing all .fq.gz, genome indexing, genome FASTA and BED files
file_dir=/home/BCG2022_genomics_exam      
          
create_directories() {
	mkdir "$my_wd/$1/"
	mkdir "$my_wd/$1/fastqc/"
	mkdir "$my_wd/$1/bam/"
	mkdir "$my_wd/$1/bam/bamcov"
	mkdir "$my_wd/$1/qualimap/"
	mkdir "$my_wd/$1/multiqc/"
	mkdir "$my_wd/$1/vcf/"
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
         -h) 
	    echo "Variant Calling
		
		Parameters are -ad and -ar for autosomic dominant and autosomic 
	    recessive cases, respectively"
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
	cat ${ar_trio}.vcf | grep -v "#" | python /home/BCG_2024_mcominelli/greppy "ar" "basic" >> ${ar_trio}_filtered.vcf
	
	# Sorting columns by family members name
	sorting $ar_trio 

    # Intersection 
    bedtools intersect -a ${ar_trio}_sorted.vcf -b ${file_dir}/targetsPad100.bed -u > ${ar_trio}_final.vcf
	
done

# Complete analysis for autosomic dominant cases
for ad_trio in ${ad_cases[@]}
do
	# Creation of all the necessary directories
	create_directories $ad_trio

    cd ${my_wd}/${ad_trio}/bam
	
	# Alignment
    bowtie_alignment $ad_trio

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
    cat ${ad_trio}.vcf | grep -v "#" | python /home/BCG_2024_mcominelli/greppy "ad" "basic" >> ${ad_trio}_filtered.vcf

    # Sorting columns by family members name
    sorting $ad_trio 

    # Intersection
    bedtools intersect -a ${ad_trio}_sorted.vcf -b ${file_dir}/targetsPad100.bed -u > ${ad_trio}_final.vcf
	
done







