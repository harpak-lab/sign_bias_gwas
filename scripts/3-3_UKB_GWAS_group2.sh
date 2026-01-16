#!/bin/bash
#SBATCH -J ukb_gwas_g2           
#SBATCH -A MCB25074       
#SBATCH -p skx              
#SBATCH -N 1                 
#SBATCH -n 1                 
#SBATCH -c 16                  
#SBATCH -t 5:00:00           
#SBATCH --array=1-22             
#SBATCH -o logs/ukb_gwas_g2_%A_%a.out
#SBATCH -e logs/ukb_gwas_g2_%A_%a.err

# Chromosome = array index
CHR=${SLURM_ARRAY_TASK_ID}

# Run PLINK2 GWAS for ALL phenotypes at once, for GROUP 2
/home1/06809/jmcole/software/plink2/plink2 --bgen ukb_imp_chr${CHR}_v3.bgen ref-first \
  --sample ukb61666_imp_chr${CHR}_v3_s487280.sample \
  --keep   /phenos/ukb_samples_group2.txt \
  --pheno  /phenos/ukb_phenotypes.txt \
  --pheno-name AD T1D T2D SCZ ASTHMA height bmi weight mono_pct neutro_pct baso_pct \
  --covar  /phenos/ukb.covariates.txt \
  --covar-name sex age age2 sex_age sex_age2 PC1-PC20 \
  --covar-variance-standardize \
  --maf 0.001 --max-maf 0.999 \
  --hwe 1e-10 midp keep-fewhet \
  --geno 0.05 \
  --glm hide-covar cols=+a1freq \
  --threads 16 \
  --memory 60000 \
  --out gwas_group2_chr${CHR}