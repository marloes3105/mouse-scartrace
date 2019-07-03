# Mouse scartrace
Here you will find all information and code for processing of mouse scartrace data. This includes transcriptome and scar analysis.

# Usage
## Scar pipeline
The scar pipeline consists of two parts: 
1. scars_pipeline.sh - this processes raw fastq files to a count table. 
2. 20190627_ScarsPipeline_VAN2988.ipynb - this takes the count table as input and transforms raw counts to allele counts, to eventually cluster cells based on their scar pattern.

You can find these scripts in the folder scar-analysis.

### 1. scars_pipeline.sh
This pipeline runs on modules from the [SingleCellMultiOmics package](https://github.com/BuysDB/SingleCellMultiOmics). You need to download this first into your environment to be able to use the scripts. The script expects you to organise your fastqs in such a way that fastq files of every separate library are in a separate folder. To execute the script, simply move to the folder where you stored all folders with fastqs, activate your python or conda environment and execute the script.  

Step-by-step overview of the pipeline:
1. Raw fastq files are demultiplexed using demux.py from the SingleCellMultiOmics package. This will create a raw_demultiplexed folder with folders for every separate library, storing the demultiplexed and rejected reads for R1 and R2. The script will now move into this raw_demultiplexed folder to execute the next steps.
2. Demultiplexed reads are trimmed using trim_galore/cutadapt.
3. Reads are mapped to the **whole unmasked** genome using bwa mem.
4. The mapped bam file is tagged using addDigestTags.py from the SingleCellMultiOmics package. This places all information that was placed in the readname for alignment back into the bamfile as tags. The option '--scar' adds a tag describing the misalignment (scar) to the read, or states 'WT' when the read maps perfectly. This simultaneously also gives a mean read quality based on the average of phred scores for each base. The option '-alleles' requires a .vcf input file and adds allelic information in a tag. The following tags are added with this script:
- 'DA' - stores allele information
- 'SD' - stores scar information
- 'DS' - stores read start position (site)
- 'SQ' - stores mean base quality
5. The bam file is filtered based on the mean base quality ('SQ') tag - only reads with a mean base quality of >0.98 will be stored and used for generating the count table. Everything below this threshold is considered to be noise.
6. Two count tables are made for the **unfiltered** bam:
a) The mapped and tagged bam file is converted to a count table using bamToCountTable.py from the SingleCellMultiOmics package. It takes '-joinedFeatureTags' as columns and '-sampleTags' as rows for this count table. For columns, we use the 'SM' tag (samplename). For rows, we use the following tags: 'chrom' (chromosome), 'DA' (allele), 'DS' (site), 'SD' (scar). 
b) The second count table that is generated stores cells ('SM') as rows and the mean base quality score ('SQ') as columns. This count table can be used to check the SQ values for all reads in the unfiltered bam, and can be used to set a sensible threshold for filtering or to check if the threshold used makes sense.
7. Two count tables are made for the **filtered** bam - these are the same as specified in 6a and 6b. 
8. Optional: A count table of the deduplicated bam file is made, the same way as specified in 6a.

### 2. 20190627_ScarsPipeline_VAN2988.ipynb
This jupyter notebook takes the count table as input and processes raw counts into percentages and after that in allele counts. Cells are filtered based on their amount of counts, and cells that pass the threshold will be used for clustering using [IWSS](https://github.com/BuysDB/IWSS). All steps are annotated in the notebook itself, but if we go over them in short, the notebook does the following:
- The count table/dataframe is read in as df.
- First, some quality checks are done. Important to note here is that our target sites fall into the IgH locus on chromosome 12, so informative counts should be located on chromosome 12. It is essential to check whether there is noise (amplification on other chromosomes) since we know that our primers also amplify several regions on other chromosomes. We can consider to include these 'off targets' later as an internal control. Several figures are made for quality control, to begin with the total amount of counts per cell are plotted, and the total amount of counts for chromosome 12 only per cell are plotted. This will give an indication how many counts there are on average. Total counts per chromosome are also checked. Next, total counts for all gRNA target sites are checked. These are also checked in an allele-specific manner, so we can see exactly which target sites are not covered on either the B6 or the 129 genome. Lastly, individual gastruloids are checked by plotting the total amount of counts per chromosome for each gastruloid. 
- The following computing steps should only done once, and dataframes are saved so they can easily be loaded again next time for fast continuation with downstream analyses.
- After confirming that we have enough informative reads to work with, raw reads are transformed to percentages. This is done separately for the two alleles (129 and B6). Note that we also take along 'nonallelic', these reads/counts do not carry any information on whether they belong to the 129 or B6 allele. For now they are taken along in the analysis, since this way hopefully we can gather more information on these counts and figure out whether we can use them (and how). Note that we only save scars if they are present for >5% on that site + allele + cell, to reduce the size of the output matrix. For this reason, we later only save the 10 most common (highest percentage) scars for a specific cell + site + allele. We can do this because in theory this should already be black and white 100% or 0% - one cell on one site and one allele can in theory only have one sequence. However, we have to deal with a high background level of noise here so this is unfortunately not true in the data. In the nonallelic counts, we theoretically have 2 alleles (we cannot separate them) so here we can come across 2 sequences for a specific site within a cell.
- Next, we convert these percentages to allele counts. We set an arbitrary threshold of 70%: an allele count is only saved if the percentage of that sequence in a cell on one site and one allele is higher than 70%. For the nonallelic counts, we can in theory have 2 so here we say that a sequence that is present >70% is 2 counts, and between 30% and 70% is one count. Here we lose a lot of information, because of the high level of noise. I have done a short test with lowering the percentage threshold of 70% to 55% and this appears to keep much more information while not changing the outcome.
- After computing allele counts, we can transform the datasets to be able to implement them in downstream code. First we make the dataframe suitable for sparse distance matrix clustering using IWSS. We only save information of the gRNA target sites for this, since we want the clustering to be based on information within these sites.
- We also transform the dataset for heatmap plotting, where rows are cells and columns are sites-alleles. Here, black is always WT and all scars have different colours.
- If all these dataframes are computed once, there is no need to compute them again. They are saved and under 'load computed datasets', they can be easily loaded into the notebook again.
- Next, we can plot the results of computing percentages and allele counts. These plots are separated per gastruloid and per site.
- We can also plot the 'heatmap' dataframe as shown underneath the header 'heatmap' plot - to plot all scars for all sites. For this, we first define all scars present in the dataframe. Then every scar is coupled to a unique colour to visualize the variety of scars within the dataset.
- Lastly, we can cluster the cells based on their scar pattern using information weighted sparse sample distance (IWSS). First, we check the amount of counts per cell in a histogram. We can check the amount of counts in total (scars + wt counts) or scar counts only. We can also plot this for the whole dataframe or separated per gastruloid. Before clustering, we filter out cells with less than *n* counts. In this example, we only keep cells with >=4 counts (scars + wt) for clustering. We can still optimize this step to retain the most information (amount of cells). We can also test here if the clustering changes when we only use scar counts as input, or when we use both scar and wt counts. For now, both scar and wt counts are used. This script produces a distance matrix which can be displayed as a heatmap. We can also use this distance matrix to separate all cells of a gastruloid in *n* amount of clusters based on their similarity in counts. For now, I tested clustering cells for each gastruloid in 20 or in 40 clusters. This number still needs to be optimised as well.

### data
All computed dataframes can be found in the folder scars_dataframes. All computed plots can be found in the folder scars_plots.


## Transcriptome pipeline
