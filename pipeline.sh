#!/bin/bash

###############################################################################################
# This script allows to do variant calling starting from the single-end short reads of a trio
# child-mother-father 
# The implemented variant caller is freebayes
###############################################################################################
# This script is written by Marco Cominelli and Elena Sasso 
# from Università degli Studi di Milano and Politecnico di Milano


# You need 
# - empty dir
# - blah blah



my_wd=/home/BCG_2024_mcominelli/project             # Type here the directory in which all the output folders will be created (make sure                                                    it's empty)

file_dir=/home/BCG2022_genomics_exam                # Here insert the directory in which the fq.gz/index/fasta/BED files are in


declare -a ar_cases=()  # Array for autosomic recessive cases
declare -a ad_cases=()  # Array for autosomic dominant cases
current_param=""        # Parametro corrente

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
	    echo "Parameters are -ad and -ar for autosomic dominant and autosomic 
	    recessive cases, respectively"
	    exit 0
	    ;;
          *)
            if [ "$current_param" == "-ar" ]; then
                ar_cases+=("$1")  # Aggiungi il valore a ar_cases
            elif [ "$current_param" == "-ad" ]; then
                ad_cases+=("$1")  # Aggiungi il valore a ad_cases
            fi
            shift
            ;;
    esac
done

echo "AR CASES VECTOR: ${ar_cases[@]}"
echo "AD CASES VECTOR: ${ad_cases[@]}"

for ar_trio in ${ar_cases[@]}
do
	mkdir "$my_wd/$ar_trio/"
	mkdir "$my_wd/$ar_trio/fastqc/"
	mkdir "$my_wd/$ar_trio/bam/"
	mkdir "$my_wd/$ar_trio/bam/bamcov"
	mkdir "$my_wd/$ar_trio/qualimap/"
	mkdir "$my_wd/$ar_trio/multiqc/"
	mkdir "$my_wd/$ar_trio/vcf/"
	
	cd ${my_wd}/${ar_trio}/bam

	bowtie2 -U ${file_dir}/${ar_trio}_father.fq.gz -x ${file_dir}/uni --rg-id 'SF' --rg "SM:father" | samtools view -Sb | samtools sort -o ${ar_trio}_father.bam
	bowtie2 -U ${file_dir}/${ar_trio}_child.fq.gz -x ${file_dir}/uni --rg-id 'SC' --rg "SM:child" | samtools view -Sb | samtools sort -o ${ar_trio}_child.bam
        bowtie2 -U ${file_dir}/${ar_trio}_mother.fq.gz -x ${file_dir}/uni --rg-id 'SM' --rg "SM:mother" | samtools view -Sb | samtools sort -o ${ar_trio}_mother.bam
	
	samtools index ${ar_trio}_father.bam 
        samtools index ${ar_trio}_child.bam
        samtools index ${ar_trio}_mother.bam
	
	cd ./bamcov

	bedtools genomecov -ibam ../${ar_trio}_father.bam -bg -trackline -trackopts 'name="father"' -max 100 > ${ar_trio}_fatherCov.bg
	bedtools genomecov -ibam ../${ar_trio}_child.bam -bg -trackline -trackopts 'name="child"' -max 100 > ${ar_trio}_childCov.bg
	bedtools genomecov -ibam ../${ar_trio}_mother.bam -bg -trackline -trackopts 'name="mother"' -max 100 > ${ar_trio}_motherCov.bg

	cd ./../../vcf

	freebayes -f ${file_dir}/universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 ../bam/${ar_trio}_mother.bam ../bam/${ar_trio}_child.bam ../bam/${ar_trio}_father.bam  > ${ar_trio}.vcf


 	# filtration
        grep "#" ${ar_trio}.vcf > ${ar_trio}_filtered.vcf
	cat ${ar_trio}.vcf | grep -v "#" | python /home/BCG_2024_mcominelli/try "ar" "standard" >> ${ar_trio}_filtered.vcf
	
	# sorting of family member names
	bcftools query -l ${ar_trio}_filtered.vcf | sort > samples.txt     # extract samples, sort names, save to file
        bcftools view -S samples.txt ${ar_trio}_filtered.vcf > ${ar_trio}_sorted.vcf # print samples in “sorted” order to new file

        # bedtools
        bedtools intersect -a ${ar_trio}_sorted.vcf -b ${file_dir}/targetsPad100.bed -u > ${ar_trio}_final.vcf



	# fastqc /home/BCG2022_genomics_exam/$ar_trio* --outdir "$my_wd/$ar_trio/fasqtc"     #FASTQC for all members of TRIO
	

done




for ad_trio in ${ad_cases[@]}
do
        mkdir "$my_wd/$ad_trio/"
        mkdir "$my_wd/$ad_trio/fastqc/"
        mkdir "$my_wd/$ad_trio/bam/"
        mkdir "$my_wd/$ad_trio/bam/bamcov"
        mkdir "$my_wd/$ad_trio/qualimap/"
        mkdir "$my_wd/$ad_trio/multiqc/"
        mkdir "$my_wd/$ad_trio/vcf/"

        cd ${my_wd}/${ad_trio}/bam

        bowtie2 -U ${file_dir}/${ad_trio}_father.fq.gz -x ${file_dir}/uni --rg-id 'SF' --rg "SM:father" | samtools view -Sb | samtools sort -o ${ad_trio}_father.bam
        bowtie2 -U ${file_dir}/${ad_trio}_child.fq.gz -x ${file_dir}/uni --rg-id 'SC' --rg "SM:child" | samtools view -Sb | samtools sort -o ${ad_trio}_child.bam
        bowtie2 -U ${file_dir}/${ad_trio}_mother.fq.gz -x ${file_dir}/uni --rg-id 'SM' --rg "SM:mother" | samtools view -Sb | samtools sort -o ${ad_trio}_mother.bam

        samtools index ${ad_trio}_father.bam
        samtools index ${ad_trio}_child.bam
        samtools index ${ad_trio}_mother.bam

        cd ./bamcov

        bedtools genomecov -ibam ../${ad_trio}_father.bam -bg -trackline -trackopts 'name="father"' -max 100 > ${ad_trio}_fatherCov.bg
        bedtools genomecov -ibam ../${ad_trio}_child.bam -bg -trackline -trackopts 'name="child"' -max 100 > ${ad_trio}_childCov.bg
        bedtools genomecov -ibam ../${ad_trio}_mother.bam -bg -trackline -trackopts 'name="mother"' -max 100 > ${ad_trio}_motherCov.bg

        cd ./../../vcf

        freebayes -f ${file_dir}/universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 ../bam/${ad_trio}_mother.bam ../bam/${ad_trio}_child.bam ../bam/${ad_trio}_father.bam  > ${ad_trio}.vcf


        # filtration
        grep "#" ${ad_trio}.vcf > ${ad_trio}_filtered.vcf
        cat ${ad_trio}.vcf | grep -v "#" | python /home/BCG_2024_mcominelli/try "ad" "high" >> ${ad_trio}_filtered.vcf

        # sorting of family member names
        bcftools query -l ${ad_trio}_filtered.vcf | sort > samples.txt     # extract samples, sort names, save to file
        bcftools view -S samples.txt ${ad_trio}_filtered.vcf > ${ad_trio}_sorted.vcf # print samples in “sorted” order to new file

        # bedtools
        bedtools intersect -a ${ad_trio}_sorted.vcf -b ${file_dir}/targetsPad100.bed -u > ${ad_trio}_final.vcf



        # fastqc /home/BCG2022_genomics_exam/$ad_trio* --outdir "$my_wd/$ad_trio/fasqtc"     #FASTQC for all members of TRIO


done







