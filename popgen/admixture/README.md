# Steps to run admixure incl filtering adn linkage pruning

mkdir admixture
cd admixture
mkdir lists log input_vcf bin

# get included individuals
cat ../../metadata/list_samples_ecu_lysimnia ../../metadata/list_samples_ecu_col_polymnia.txt > lists/samples.txt

#get list of vcf files
VCF_PATH=../ref_polymnia/variant_calling/output_vcf/vcfs/

ls $VCF_PATH/ | sed 's/\.vcf\.gz//' > lists/list_vcfs.txt

# Get depth for filtering
stat_file=/lustre/scratch125/tol/teams/meier/users/kn9/recombination/pop_gen/ref_polymnia/variant_calling/output_vcf/vcf_qc/tables/per_sample_summary.tsv

tail -n +2 $stat_file | cut -f7 | sort -n | awk '{a[NR]=$1; s+=$1}
  END{med=(NR%2)?a[(NR+1)/2]:(a[NR/2]+a[NR/2+1])/2
      printf "n=%d mean=%.2f median=%.2f\n", NR, s/NR, med}'

n=1760 mean=24.21 median=27.10

tail -n +2 $stat_file | cut -f1 |sort | uniq | wc -l

# or from random autosome
VCF=SUPER_11.call.MQ_BQ20
bcftools query -f '%INFO/DP\n' $VCF_PATH/${VCF}.vcf.gz |  grep -v '^\.$' | sort -n \
 | awk '{a[NR]=$1}
   END{
     med = (NR%2) ? a[(NR+1)/2] : (a[NR/2] + a[NR/2+1]) / 2
     printf "n=%d median=%.1f\n", NR, med
   }'
#n=14370243 median=2940.0

bcftools query -f '%INFO/DP\n' $VCF_PATH/${VCF}.vcf.gz | awk '{s+=$1; ss+=$1^2; n++} END{m=s/n; sd=sqrt(ss/n-m^2); printf "mean=%.1f sd=%.1f\n", m, sd}'
#mean=2837.3 sd=558.7


# Filter with bcftools
#add filters and make a script
cat > bin/bcftools_filter.sh << 'EOF'
#!/bin/bash
set -euo pipefail

VCF=$(sed -n -e ${LSB_JOBINDEX}p lists/list_vcfs.txt)

# filters
MIN=3
MAX=50
n_ind=110
MIN_AF=0.03
#3/110

MIN_FILTER=330
MAX_FILTER=5000

VCF_PATH=../ref_polymnia/variant_calling/output_vcf/vcfs/
OUT=input_vcf/${VCF}_filtered_DP${MIN_FILTER}_${MAX_FILTER}.vcf.gz

ml bcftools/1.20--h8b25389_0
bcftools view -S samples.txt $VCF_PATH/${VCF}.vcf.gz \
    | bcftools view -m2 -M2 -v snps \
        -i "QUAL>=30 && INFO/DP>=${MIN_FILTER} && INFO/DP<=${MAX_FILTER}"  \
    | bcftools +setGT -Oz -o $OUT -- -t q -n . -i "FMT/DP<${MIN}"
        
bcftools index $OUT
EOF

bsub -q normal -J bcftools_filter[1] -M 1000 -R "select[mem>1000] rusage[mem=1000] span[hosts=1]" -o log/bcftools.%J.%I.out -e log/bcftools.%J.%I.err -n2 "bash bin/bcftools_filter.sh"

bsub -q normal -J bcftools_filter[2-16] -M 1000 -R "select[mem>1000] rusage[mem=1000] span[hosts=1]" -o log/bcftools.%J.%I.out -e log/bcftools.%J.%I.err -n2 "bash bin/bcftools_filter.sh"

# Concatenate
#make list of filtered vcf and sort them , only autosomes
ls input_vcf/*.vcf.gz | sort -V | grep -E -v "SUPER_W|SUPER_Z" > lists/list_sorted_vcfs.txt 

bsub -q normal -M 1000 -R "select[mem>1000] rusage[mem=1000] span[hosts=1]" -o log/bcftools.%J.%I.out -e log/bcftools.%J.%I.err -n4 "
bcftools concat -f lists/list_sorted_vcfs.txt \
    --threads 4 -Oz -o input_vcf/filtered_concat.vcf.gz
"


# Step 1: Generate the input file in plink format

mkdir output

bsub -q normal -M 10000 -R "select[mem>10000] rusage[mem=10000] span[hosts=1]" -o log/plink.%J.out -e log/plink.%J.err -n8 "conda activate plink
plink2 --vcf input_vcf/filtered_concat.vcf.gz \
       --double-id --allow-extra-chr \
       --set-missing-var-ids @:# \
       --vcf-min-dp ${MIN} \
       --geno 0.1 --maf ${MIN_AF} \
       --make-bed \
       --threads 8 --out output/all
    conda deactivate"


Note: 14 nonstandard chromosome codes present.
27699118 variants loaded from output/all-temporary.pvar.zst.
Note: No phenotype data present.
Calculating allele frequencies... done.
--geno: 8569877 variants removed due to missing genotype data.
10271643 variants removed due to allele frequency threshold(s)
(--maf/--max-maf/--mac/--max-mac).
8857598 variants remaining after main filters.


# Stap 2: Find independent loci, prune, and make pca and bed-file for admixture

bsub -q normal -M 5000 -R "select[mem>5000] rusage[mem=5000] span[hosts=1]" -o log/plink.%J.out -e log/plink.%J.err -n4 "conda activate plink 
plink2 --bfile output/all --allow-extra-chr --indep-pairwise 50 10 0.2 --threads 4 --out output/all
plink2 --bfile output/all --allow-extra-chr --extract output/all.prune.in --make-bed --threads 4 --out output/pruned
plink2 --bfile output/pruned --allow-extra-chr --pca 10 --threads 4 --out output/pca
conda deactivate"



# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0

ADMIX_PATH=/software/team347/kn9/admixture_linux-1.4.0/

awk '{$1="0";print $0}' output/pruned.bim > output/pruned.bim.tmp;mv output/pruned.bim.tmp output/pruned.bim


cat > bin/admixture_K.sh << 'EOF'
#!/bin/bash
set -euo pipefail
ADMIX_PATH=/software/team347/kn9/admixture_linux-1.4.0/

K=${LSB_JOBINDEX}
SEED=42

$ADMIX_PATH/admixture --cv -j4 -s $SEED output/pruned.bed ${K} > output/log_K${K}.out
EOF

bsub -J "admix[1-8]" -q normal -M 4000 -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
     -o log/admixture.%J.%I.out -e log/admixture.%J.%I.err -n4 \
     "bash bin/admixture_K.sh"


grep -h "CV error" output/log*.out | sed 's/CV error (K=\([0-9]*\)): /\1 /' > output/cv.txt

# create a tsv metadata file with sample in col1 and pop in col2, no header
sed 's/.*_/M\. /' lists/samples.txt | paste lists/samples.txt - > metadata.tsv

Rscript bin/plot_admixure_pca.R


#check outlier

cut -f1,4 output/pca.eigenvec | sort -k2