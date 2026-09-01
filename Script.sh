# Imputation pipeline

# Input data quality control
bcftools view -i '(REF~"A" || REF~"C" || REF~"G" || REF~"T") & (ALT~"A" || ALT~"C" || ALT~"G" || ALT~"T")' ${input} -o ${output}
bcftools plugin fill-tags --threads 8 ${input} -o ${output} -- -t AC,AF,HWE,MAF,F_MISSING
bcftools view -m2 -M2 -i 'F_MISSING < 0.05 & HWE > 1e-10 & AC > 0' -o ${output} ${input}

# Separate into chromosome level
for CHR in {1..22};do
bcftools view -t chr${CHR} ${input} -o output_chr${CHR}.vcf.gz
tabix -p vcf output_chr${CHR}.vcf.gz
done

# Imputation with Han-sv panel
for CHR in {1..22}; do
java -Xmx200g -jar beagle.28Jun21.220.jar \
ref=Han_SV_panel_chr${CHR}.vcf.gz \
gt=input_chr${CHR}.vcf.gz \
out=output_chr${CHR} \
chrom=chr${CHR} \
map=beagle_chr${CHR}.map \
impute=true
tabix -p vcf output_chr${CHR}.vcf.gz
done

# QC before GWAS
bcftools view -i "DR2 >= 0.7" ${input} -Oz -o ${output}

# GWAS pipeline
plink --vcf ${input} --maf 0.01 --geno 0.2 --recode 12 --output-missing-genotype 0 --transpose --out output --allow-extra-chr
emmax-kin-intel64 -v -s -d 10 ${tped_prefix}
emmax-intel64 -v -d 10 -t ${tped_prefix} -p ${pheno_file} -k ${kin_file} -c ${cov_file} -o ${output}

# Conditional analysis
gcta64 --bfile ${prefix} --maf 0.01 --cojo-p ${p_value} --cojo-file ${summary_statistics_of_phenotype}.ma --cojo-slct --cojo-actual-geno --out ${output}

# Fine mapping
Rscript run_susie_rss.R --z ${gwas_summary_statistics} --ld ${reference_ld.matrix} --n ${sample_size} --L ${max_causal_variants} --out ${output}

# Annotation
# annotation for SNPs and InDels
perl convert2annovar.pl -keepindelref -format vcf4 -allsample -withfreq ${input} -outfile ${output_prefix}
perl table_annovar.pl  ${input} /annovar/humandb -buildver hg38 -out ${output_prefix} -remove -protocol refGene -operation g -nastring . -polish
# annotation for SVs
AnnotSV -SVinputFile ${input} -outputFile ${output}









