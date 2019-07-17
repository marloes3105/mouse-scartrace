#!/bin/bash
# Marijn van Loenhout 2018

# Make a copy the script to your own directory and modify as needed, (indexes and GTF most likely, but aligner and trimming commands may also be changed examples are in the comments below, proper use of paired end reads will need modifications but your on your own there)

# create a conda environment (examples below use conda_env) and install:
# STAR, (or other aligner, for a BWA/BT2 direct to transcriptome implementation see commented code at the end of script)
# samtools
# umi-tools
# cutadapt (trimgalore/atropos)
# subread
# colorama

## preferably use a STAR index of readlenght -1 and specify index and GTF here:
STAR_index="/hpc/hub_oudenaarden/group_references/ensembl/95/mus_musculus/star_index_NOMASK_74/"
GTF="/hpc/hub_oudenaarden/group_references/ensembl/95/mus_musculus/annotations.gtf"

#DEMUX reads first with demux.py, for single end CS2 reads use the option -use CS2-CB8-U6:
# to check the barcodes present in the libary run the demux script without the -use and execute (--y) flag this will analyze the first 10000 reads and show the barcodes present:
## python hpc/hub_oudenaarden/bdebarbanson/internalTools/modularDemultiplexer/demux.py
##
## to execute the demultiuplexing then choose your libary type (e.g. -use CS2C8U6) and submit the job:
##submission.py -y --nenv 'source activate conda_env;/hpc/hub_oudenaarden/bdebarbanson/internalTools/modularDemultiplexer/demux.py -use CS2C8U6 */*.fastq.gz --y'

# After demultiplexing run this pipe from inside the directory containing demultiplexedR2.fastq.gz
## submission.py -y -m 50 --nenv 'source activate conda_env;"/hpc/hub_oudenaarden/group_scripts/CS2_pipe.sh"'
# the cluster/.stdout file will have a log documenting all steps of the pipe


### START OF PIPE, please read comments at each step to make sure it suits your needs:

printf "\n[+] CS2 pipe started at: $(date)\n"

### READ TRIMMIMG
##trim using cutadapt/trimgalore/atropos (all based on cutadapt)
#settings below for cutadapt. trim_galore may also be used with:
# trim_galore $fqfiles -o ${destdir} ## + your favourite options

#cutadapt will only trim one adapter from a read! if there is a need for 3' and 5' trimming the command needs to be executed twice, example below: 
#
# cutadapt -g ^G -j 6 -o "${destdir}demultiplexedR2.trimmed_1.fastq" $fqfiles > "${destdir}trimreport_1.txt"
# cutadapt -a "A{100};min_overlap=5" --nextseq-trim=20 --minimum-length 17 -j 6 -o "${destdir}demultiplexedR2.trimmed.fastq" "${destdir}demultiplexedR2.trimmed_1.fastq.gz" > "${destdir}trimreport.txt"
# rm "${destdir}demultiplexedR2.trimmed_1.fastq"
#

for fqfiles in $(find $PWD -name 'demultiplexedR2.fastq.gz')
do
  printf "\n[+] trimming: $fqfiles\n"
  fqdir=${fqfiles%demultiplexedR*}
  destdir="${fqdir}trimmed/"
  mkdir -p "$destdir"
  cutadapt -a "A{100};min_overlap=5" --nextseq-trim=20 -a "CTGTCTCTTATA;min_overlap=6" -a "GTTCAGAGTTCTACAGTCCGA;min_overlap=6" -a "TCGGACTGTAGAACTCTGAAC;min_overlap=6" -g "GTTCAGAGTTCTACAGTCCGA;min_overlap=6" --minimum-length 17 -j 6 -o "${destdir}demultiplexedR2.trimmed.fastq.gz" $fqfiles > "${destdir}trimreport.txt"
done
### check the trimreport file to asses if triming needs modification, generally STAR is pretty insensitive to trimming and will softclip the ends properly


###MAPPING
## read the well documented STAR manual to understand the options, to additionaly obtain a transcriptome aligned bam file use: 
# --quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend
for fqfiles in $(find $PWD -path '*trimmed/demultiplexedR2.trimmed.fastq.gz')
do
  printf "\n[+] mapping: $fqfiles\n"
  fqdir=${fqfiles%trimmed/demultiplexedR*}
  destdir="${fqdir}aligned/"
  mkdir -p "$destdir"
  STAR --runThreadN 8 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 20 --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat --readFilesIn ${fqfiles} --outFileNamePrefix ${destdir} --genomeDir $STAR_index
#  rm -r ${destdir}_STARtmp
#  rm -r $fqfiles
done

### TAGGING bam file, this script places all the data, that was put in the readname for alignment, back into the bam file as tags
### use addDigestTags.py --showtags to show all available tags in bam file, if in need of more tags contact Buys
for bamfiles in $(find $PWD -path '*aligned/Aligned.sortedByCoord.out.bam')
do
  printf "\n[+] tagging: $bamfiles\n"
  indir=${bamfiles%aligned/Aligned.sortedByCoord.*}
  destdir="${indir}tagged"
  samtools index $bamfiles
  python /hpc/hub_oudenaarden/group_scripts/addDigestTags.py --ftag -o $destdir ${bamfiles}
#  rm $bamfiles;rm $bamfiles.bai
#  rm ${indir}aligned/*.bam ##delete transcriptome bam
done

#### FEATURE COUNTING
## for mutimapper counting, use -M and set -Q 0 (map quality= 0 or a litte more to determine cutoff for number of allowed multimap locations, see STAR manual)
## featureCounts -s 1 -T 1 -R BAM -Q 2 -M -O -a ${GTF} -o "$destdir/fCounts.txt" ${fCinputfiles}
for fCinputfiles in $(find $PWD -path '*tagged/Aligned.sortedByCoord.out.bam')
do
  printf "\n[+] FeatureCounting: $fCinputfiles\n"
  indir=${fCinputfiles%tagged*}
  destdir="${indir}FCounts"
  mkdir -p "$destdir"
  featureCounts -s 1 -T 1 -R BAM -a $GTF -o "$destdir/fCounts.txt" ${fCinputfiles}
  samtools sort -T $destdir/temp "$destdir/Aligned.sortedByCoord.out.bam.featureCounts.bam" -o "$destdir/SortedByCoord.featureCounts.bam"
  samtools index "$destdir/SortedByCoord.featureCounts.bam"
#  rm "$destdir/Aligned.sortedByCoord.out.bam.featureCounts.bam"
#  rm ${fCinputfiles};rm ${fCinputfiles}.bai
done

###DEDUPLICATION, umitools directional deduplication at the gene level
##(gene tags  added by featureCounts are used to deduplicate) multiple transcripts can make up a gene so this is a more strict dedup accounting for the variable lenght of reads from the same molecule after IVT

for UTinputfiles in $(find $PWD -path '*FCounts/SortedByCoord.featureCounts.bam')
do
  printf "\n[+] deduplicating: $UTinputfiles\n"
  indir=${UTinputfiles%FCounts*}
  destdir="${indir}dedup"
  mkdir -p "$destdir"
  umi_tools dedup --extract-umi-method "tag" --umi-tag=MI:Z: --per-gene --gene-tag=XT:Z: -L "$destdir/dedup.log" -I $UTinputfiles -S "$destdir/dedup.bam"
done

### COUNT TABLES
### make dedup counts table, expects fastq directory to end in *_S* for sample/RPI index number
### use --showtags to show all available tag in bam file
for bamfiles in $(find $PWD -path '*dedup/dedup.bam')
do
  printf "\n[+] making dedup counts table: $bamfiles\n"
  sampledir=${bamfiles%/dedup/dedup*}
  sample=${sampledir##*_S}
  destdir="${sampledir}/counts"
  outfile="$destdir/S${sample}_dedup_counts.csv"
  mkdir -p "$destdir"
  python /hpc/hub_oudenaarden/group_scripts/modularDemultiplexer/featureCountsTaggedBamFileToCountTable.py -o $outfile -featureTags XS,XT -sampleTags aI,BI ${bamfiles} 
done

### make no dedup counts table, expects fastq directory to end in *_S* for sample/RPI index number
for bamfiles in $(find $PWD -path '*FCounts/SortedByCoord.featureCounts.bam')
do
  printf "\n[+] making no dedup counts table: $bamfiles\n"
  sampledir=${bamfiles%/FCounts/Sorted*}
  sample=${sampledir##*_S}
  destdir="${sampledir}/counts"
  outfile="$destdir/S${sample}_no_dedup_counts.csv"
  mkdir -p "$destdir"
  python /hpc/hub_oudenaarden/group_scripts/modularDemultiplexer/featureCountsTaggedBamFileToCountTable.py -o $outfile -featureTags XS,XT -sampleTags aI,BI ${bamfiles} 
done

printf "\n[+] CS2 pipe complete at: $(date)\n"



####mapping directly to the transcriptome with BWA mem (i.e. build and index from the transcriptome fasta), bowtie is similar just substitute bowtie in the mapping command
##
##replace the STAR pipeline from the mapping comand on down
#(if however aligning to the genome only the mapping command should be subsituted and featureCounts should be used to assign genes/transctipts as in the STAR implementation)
#
##bwa_index="path to bwa index"
#
#
#for fqfiles in $(find $PWD -path '*trimmed/demultiplexedR2.trimmed.fastq')
#do
#  printf "\n[+] mapping: $fqfiles\n"
#  fqdir=${fqfiles%trimmed/demultiplexedR*}
#  destdir="${fqdir}aligned/"
#  mkdir -p "$destdir"
#  bwa mem -t 8 $bwa_index ${fqfiles} | samtools sort -@8 -o "$destdir/demultiplexedR2.bam" -
#  #rm $fqfiles
#done
#
#### Filter stranded reads, Important for rnaseq!
# in the STAR implementation featureCounts takes care of the strandedness here it needs to be seperately implemented

#for bamfiles in $(find $PWD -path '*aligned/demultiplexedR2.bam')
#do
#  printf "\n[+] filtering sense strand: $bamfiles\n"
#  samtools view -F 16 -b $bamfiles > $bamfiles.forward
#  rm $bamfiles
#  mv $bamfiles.forward $bamfiles
#done
#
#### tag
#for bamfiles in $(find $PWD -path '*aligned/demultiplexedR2.bam')
#do
#  printf "\n[+] tagging: $bamfiles\n"
#  outdir=${bamfiles%aligned/demultiplexedR2.*}
#  destdir="${outdir}tagged"
#  samtools index $bamfiles
#  python /hpc/hub_oudenaarden/group_scripts/addDigestTags.py --ftag -o $destdir ${bamfiles}
#  #rm $bamfiles
#done
#
####dedup
#for UTinputfiles in $(find $PWD -path '*tagged/demultiplexedR2.bam')
#do
#  printf "\n[+] deduplicating: $UTinputfiles\n"
#  UTdir=${UTinputfiles%tagged*}
#  destdir="${UTdir}dedup"
#  mkdir -p "$destdir"
#  umi_tools dedup --extract-umi-method "tag" --umi-tag=MI:Z: -L "$destdir/dedup.log" -I $UTinputfiles -S "$destdir/dedup.bam"
#done
#
####count tables
#for bamfiles in $(find $PWD -path '*dedup/dedup.bam')
#do
#  printf "\n[+] making count table: $bamfiles\n"
#  outdir=${bamfiles%dedup/dedup.bam}
#  destdir="${outdir}count_table"
#  outfile="$destdir/S${sample}_dedup_counts.csv"
#  mkdir -p "$destdir"
#  python /hpc/hub_oudenaarden/group_scripts/modularDemultiplexer/taggedBamFileToCountTable.py -o $outfile -featureTags chrom -sampleTags aI,BI ${bamfiles} 
#done
#
#
#
