
LIST=$1
#list with genomes in phylogenetic order

mkdir -p config_files

for species in $(cat $LIST)
do
ref=${species}
query=$(grep -A1 "${species}" $LIST |tail -n1)
echo "Reference is: $ref"
echo "Query is: $query"

echo "
# This is an example config file for odp/scripts/odp
#
# # To use this software first copy this config file to your analysis directory
# cp odp/example_configs/CONFIG_odp.yaml ./config.yaml
# # Then modify the config file to include your own data
# vim config.yaml
# # Then run the pipeline
# snakemake -r -p --snakefile odp/scripts/odp


ignore_autobreaks: True       # Skip steps to find breaks in synteny blocks
diamond_or_blastp: "diamond"  # "diamond" or "blastp"
duplicate_proteins: "fail"    # currently only "fail" or "best". Fail doesn't allow duplicate names or seqs
plot_LGs: False                # Plot the ALGs based on the installed databases
plot_sp_sp: True              # Plot the synteny between two species, if False just generates .rbh files

species:
  $ref:
    proteins: ../../input/peptide/${ref}.fasta
    chrom: ../../input/chrom_input/${ref}.chrom
    genome: ../../input/genome/${ref}_assembly.fasta
    minscafsize: 1000000  # Only plots scaffolds that are 1 Mbp or longer
  $query:
    proteins: ../../input/peptide/${query}.fasta
    chrom: ../../input/chrom_input/${query}.chrom
    genome: ../../input/genome/${query}_assembly.fasta
    minscafsize: 1000000  # Only plots scaffolds that are 8 Mbp or larger

" > config_files/${ref}_vs_${query}_odp.yaml
done


ls config_files/*_odp.yaml > list_configs.txt
