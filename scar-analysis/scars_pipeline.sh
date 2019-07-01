#/bin/bash

### Pipeline for scar mapping, Marloes

### Outline of pipeline:
# 1. Demultiplexing using demux.py
# 2. Trimming using TrimGalore
# 3. Mapping to full genome (mm93) using BWA
# 4. Change back digest tags
# 5. Make count table (with allele specificity for the IgH locus)

### This is run from a conda environment

### specify your references:
mmus_reference="/hpc/hub_oudenaarden/group_references/ensembl/95/mus_musculus/primary_assembly_NOMASK_ERCC92.fa"
allele_reference="/hpc/hub_oudenaarden/bdebarbanson/ref/gsnvMouse/129S1_SvImJ_subset_igh.vcf.gz"

### Make sure all data to process is in your folder of choice, where files from every sample are in their separate folders starting with 'MB' (in my case), and cd to this folder
# to submit: submission.py -s ./cluster -m 50 --nenv -N scars_pipe -y "source activate conda;/hpc/hub_oudenaarden/Marloes/bin/scars_pipeline.sh"

printf "\n[+] Scar pipeline  started at: $(date)\n"


### 1. DEMUX reads first with demux.py, with tag 'SCARC8R2':
for library in *;do echo "/hpc/hub_oudenaarden/group_scripts/modularDemultiplexer/demux.py $library/*fastq.gz --y -use SCARC8R2"; done | sh

# move to folder created by demux.py
cd raw_demultiplexed/


### 2. Trimming using TrimGalore

for fqfiles in $(find $PWD -name 'demultiplexedR1.fastq.gz')
do
  printf "\n[+] trimming: $fqfiles\n"
  fqdir=${fqfiles%demultiplexedR*}
  destdir="${fqdir}trimmed/"
  mkdir -p "$destdir"
  /hpc/hub_oudenaarden/group_binaries/TrimGalore-0.4.3/trim_galore --paired --path_to_cutadapt /hpc/hub_oudenaarden/Marloes/miniconda3/envs/conda/bin/cutadapt --output_dir ${destdir} "${fqdir}demultiplexedR1.fastq.gz" "${fqdir}demultiplexedR2.fastq.gz"
done


### 3. Mapping to full genome using BWA

for fqfiles in $(find $PWD -path '*trimmed/demultiplexedR2_val_2.fq.gz')
do
  printf "\n[+] mapping: $fqfiles\n"
  fqdir=${fqfiles%trimmed/demultiplexedR*}
  destdir="${fqdir}mapped/"
  mkdir -p "$destdir"
  /hpc/hub_oudenaarden/bdebarbanson/externalTools/bwa-0.7.16a/bwa mem -t 4 -M $mmus_reference "${fqdir}trimmed/demultiplexedR1_val_1.fq.gz" "${fqdir}trimmed/demultiplexedR2_val_2.fq.gz" | /hpc/hub_oudenaarden/bdebarbanson/bin/samtools view -Sb - > "${destdir}BWAmapped_mm95.bam" 
  /hpc/hub_oudenaarden/bdebarbanson/bin/samtools sort "${destdir}BWAmapped_mm95.bam" > "${destdir}BWAmapped_mm95_sorted.bam"
  /hpc/hub_oudenaarden/bdebarbanson/bin/samtools index "${destdir}BWAmapped_mm95_sorted.bam"
## rm "$destdir/BWAmapped_mm93.bam"
done


### 4. Change back digest tags
# TAGGING bam file, this script places all the data that was put in the readname for alignment back into the bam file as tags
# use addDigestTags.py --showtags to show all available tags in bam file, if in need of more tags contact Buys

for bamfiles in $(find $PWD -path '*mapped/BWAmapped_mm95_sorted.bam')
do
  printf "\n[+] tagging: $bamfiles\n"
  indir=${bamfiles%mapped/BWAmapped_mm95_sorted.*}
  destdir="${indir}tagged"
  universalBamTagger.py --ftag --scar --dedup -alleles $allele_reference -o $destdir ${bamfiles}
#  rm $bamfiles;rm $bamfiles.bai
#  rm ${indir}mapped/*.bam ##delete mapped bam
done


### 5. Make count table (with allele specificity for the IgH locus)
# make counts table, expects fastq directory to start with library name for annotation of count table
# use --showtags to show all available tags in bam file
for bamfiles in $(find $PWD -path '*tagged/BWAmapped_mm95_sorted.bam')
do
  printf "\n[+] making no dedup count table: $bamfiles\n"
  indir=${bamfiles%tagged*}
  destdir="${indir}countTable"
  sampledir=${bamfiles%/tagged/BWAmapped*}
  sample=${sampledir##*demultiplexed/}
  outfile="$destdir/${sample}_countTable_nodedup.csv"
  mkdir -p "$destdir"
  bamToCountTable.py -o $outfile -joinedFeatureTags chrom,DA,DS,SD -sampleTags SM ${bamfiles}
done

### 5b. Make dedup count table (optional)
for bamfiles in $(find $PWD -path '*tagged/*.dedup.bam')
do
  printf "\n[+] making dedup count table: $bamfiles\n"
  indir=${bamfiles%tagged*}
  destdir="${indir}countTable"
  sampledir=${bamfiles%/tagged/BWAmapped*}
  sample=${sampledir##*demultiplexed/}
  outfile="$destdir/${sample}_countTable_dedup.csv"
  bamToCountTable.py -o $outfile -joinedFeatureTags chrom,DA,DS,SD -sampleTags SM ${bamfiles}
done

printf "\n[+] Scar pipeline complete at: $(date)\n"
#
#
#
