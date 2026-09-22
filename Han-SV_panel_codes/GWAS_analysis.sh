# Input data quality control
bcftools view -i '(REF~"A" || REF~"C" || REF~"G" || REF~"T") & (ALT~"A" || ALT~"C" || ALT~"G" || ALT~"T")' ${raw_vcf}.gz -Oz -o ${filter_vcf}.gz
bcftools plugin fill-tags --threads 8 ${filter_vcf}.gz -Oz -o ${tagged_vcf}.gz -- -t AC,AF,HWE,MAF,F_MISSING
bcftools view -m2 -M2 -i 'F_MISSING < 0.05 & HWE > 1e-10 & AC > 0' -Oz -o ${tagged_filter_vcf}.gz ${tagged_vcf}.gz
tabix -p vcf ${tagged_filter_vcf}.gz

# imputation and QC before GWAS
for CHR in {1..22}; do
    # Split by chromosome
    bcftools view -t chr${CHR} -Oz -o seperate_chr${CHR}.vcf.gz ${tagged_filter_vcf}.gz
    tabix -p vcf seperate_chr${CHR}.vcf.gz

    # Imputation using Han-SV panel
    java -Xmx200g -jar beagle.28Jun21.220.jar \
        ref=Han_SV_panel_chr${CHR}.vcf.gz \
        gt=seperate_chr${CHR}.vcf.gz \
        out=HanSV_chr${CHR} \
        chrom=chr${CHR} \
        map=beagle_chr${CHR}.map \
        impute=true
    tabix -p vcf HanSV_chr${CHR}.vcf.gz

    # Post-imputation QC (imputation accuracy threshold)
    bcftools view -i "DR2 >= 0.7" -Oz -o HanSV_chr${CHR}.filter.vcf.gz HanSV_chr${CHR}.vcf.gz
    tabix -p vcf HanSV_chr${CHR}.filter.vcf.gz
done

# merge vcf file
bcftools concat HanSV_chr{1..22}.filter.vcf.gz -Oz -o multi_chrom_filter.vcf.gz
tabix -p vcf multi_chrom_filter.vcf.gz

# GWAS pipeline
plink --vcf multi_chrom_filter.vcf.gz --maf 0.01 --geno 0.2 --recode 12 --output-missing-genotype 0 --transpose --out ${tped_prefix} --allow-extra-chr
emmax-kin-intel64 -v -s -d 10 ${tped_prefix}
emmax-intel64 -v -d 10 -t ${tped_prefix} -p ${pheno_file} -k ${tped_prefix}.aIBS.kinf -c ${cov_file} -o ${GWAS_output}

# Conditional analysis
gcta64 --bfile ${bfile_prefix} --maf 0.01 --cojo-p ${p_value} --cojo-file ${summary_statistics_of_phenotype}.ma --cojo-slct --cojo-actual-geno --out ${cojo_output}

# Fine mapping
Rscript run_susie_rss.R --z ${gwas_summary_statistics} --ld ${reference_ld.matrix} --n ${sample_size} --L ${max_causal_variants} --out ${finemap_out}

# Annotation
# annotation for SNPs and InDels
java -Xmx16g -jar snpEff.jar -v -canon -stats ${output_prefix}.summary.html -csvStats ${output_prefix}.summary.csv GRCh38.105 ${snv_vcf} > ${snv_annoted}
# annotation for SVs
AnnotSV -SVinputFile ${sv_vcf} -genomeBuild GRCh38 -outputFile ${sv_annoted}









