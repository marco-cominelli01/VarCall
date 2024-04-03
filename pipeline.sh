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
	mkdir "$my_wd/$ar_trio/"
	mkdir "$my_wd/$ar_trio/fastqc/"
	mkdir "$my_wd/$ar_trio/bam/"
	mkdir "$my_wd/$ar_trio/bam/bamcov"
	mkdir "$my_wd/$ar_trio/qualimap/"
	mkdir "$my_wd/$ar_trio/multiqc/"
	mkdir "$my_wd/$ar_trio/vcf/"
	
	cd ${my_wd}/${ar_trio}/bam

	# Alignment
	bowtie2 -U ${file_dir}/${ar_trio}_father.fq.gz -x ${file_dir}/uni --rg-id 'SF' --rg "SM:father" | samtools view -Sb | samtools sort -o ${ar_trio}_father.bam
	bowtie2 -U ${file_dir}/${ar_trio}_child.fq.gz -x ${file_dir}/uni --rg-id 'SC' --rg "SM:child" | samtools view -Sb | samtools sort -o ${ar_trio}_child.bam
    bowtie2 -U ${file_dir}/${ar_trio}_mother.fq.gz -x ${file_dir}/uni --rg-id 'SM' --rg "SM:mother" | samtools view -Sb | samtools sort -o ${ar_trio}_mother.bam
	
	# BAM indexing
	samtools index ${ar_trio}_father.bam 
    samtools index ${ar_trio}_child.bam
    samtools index ${ar_trio}_mother.bam
	
	cd ./bamcov

	# BAM coverage
	bedtools genomecov -ibam ../${ar_trio}_father.bam -bg -trackline -trackopts 'name="father"' -max 100 > ${ar_trio}_fatherCov.bg
	bedtools genomecov -ibam ../${ar_trio}_child.bam -bg -trackline -trackopts 'name="child"' -max 100 > ${ar_trio}_childCov.bg
	bedtools genomecov -ibam ../${ar_trio}_mother.bam -bg -trackline -trackopts 'name="mother"' -max 100 > ${ar_trio}_motherCov.bg

	cd ./../../vcf

	# Variant calling
	freebayes -f ${file_dir}/universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 ../bam/${ar_trio}_mother.bam ../bam/${ar_trio}_child.bam ../bam/${ar_trio}_father.bam  > ${ar_trio}.vcf

 	# Variant prioritization
    grep "#" ${ar_trio}.vcf > ${ar_trio}_filtered.vcf
	cat ${ar_trio}.vcf | grep -v "#" | python /home/BCG_2024_mcominelli/try "ar" "standard" >> ${ar_trio}_filtered.vcf
	
	# Sorting columns by family members name
	bcftools query -l ${ar_trio}_filtered.vcf | sort > samples.txt     # extract samples, sort names, save to file
    bcftools view -S samples.txt ${ar_trio}_filtered.vcf > ${ar_trio}_sorted.vcf # print samples in “sorted” order to new file

    # Intersection 
    bedtools intersect -a ${ar_trio}_sorted.vcf -b ${file_dir}/targetsPad100.bed -u > ${ar_trio}_final.vcf
	
done

# Complete analysis for autosomic dominant cases
for ad_trio in ${ad_cases[@]}
do
	# Creation of all the necessary directories
    mkdir "$my_wd/$ad_trio/"
    mkdir "$my_wd/$ad_trio/fastqc/"
	mkdir "$my_wd/$ad_trio/bam/"
    mkdir "$my_wd/$ad_trio/bam/bamcov"
    mkdir "$my_wd/$ad_trio/qualimap/"
    mkdir "$my_wd/$ad_trio/multiqc/"
    mkdir "$my_wd/$ad_trio/vcf/"

    cd ${my_wd}/${ad_trio}/bam
	
	# Alignment
    bowtie2 -U ${file_dir}/${ad_trio}_father.fq.gz -x ${file_dir}/uni --rg-id 'SF' --rg "SM:father" | samtools view -Sb | samtools sort -o ${ad_trio}_father.bam
    bowtie2 -U ${file_dir}/${ad_trio}_child.fq.gz -x ${file_dir}/uni --rg-id 'SC' --rg "SM:child" | samtools view -Sb | samtools sort -o ${ad_trio}_child.bam
    bowtie2 -U ${file_dir}/${ad_trio}_mother.fq.gz -x ${file_dir}/uni --rg-id 'SM' --rg "SM:mother" | samtools view -Sb | samtools sort -o ${ad_trio}_mother.bam

	# BAM indexing
    samtools index ${ad_trio}_father.bam
    samtools index ${ad_trio}_child.bam
    samtools index ${ad_trio}_mother.bam

    cd ./bamcov
	
	# BAM coverage
    bedtools genomecov -ibam ../${ad_trio}_father.bam -bg -trackline -trackopts 'name="father"' -max 100 > ${ad_trio}_fatherCov.bg
    bedtools genomecov -ibam ../${ad_trio}_child.bam -bg -trackline -trackopts 'name="child"' -max 100 > ${ad_trio}_childCov.bg
    bedtools genomecov -ibam ../${ad_trio}_mother.bam -bg -trackline -trackopts 'name="mother"' -max 100 > ${ad_trio}_motherCov.bg

    cd ./../../vcf

	# Variant calling
    freebayes -f ${file_dir}/universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 ../bam/${ad_trio}_mother.bam ../bam/${ad_trio}_child.bam ../bam/${ad_trio}_father.bam  > ${ad_trio}.vcf

    # Variant prioritization
    grep "#" ${ad_trio}.vcf > ${ad_trio}_filtered.vcf
    cat ${ad_trio}.vcf | grep -v "#" | python /home/BCG_2024_mcominelli/try "ad" "standard" >> ${ad_trio}_filtered.vcf

    # Sorting columns by family members name
    bcftools query -l ${ad_trio}_filtered.vcf | sort > samples.txt     # extract samples, sort names, save to file
    bcftools view -S samples.txt ${ad_trio}_filtered.vcf > ${ad_trio}_sorted.vcf # print samples in “sorted” order to new file

    # Intersection
    bedtools intersect -a ${ad_trio}_sorted.vcf -b ${file_dir}/targetsPad100.bed -u > ${ad_trio}_final.vcf
	
done







