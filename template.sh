#!/bin/bash
# -*- coding: utf-8 -*-
# Description: Opentools pipeline
# creator: CZZ
# modified: 2025-Jun.29th
# version: 0.1.0

export PATH="/Users/server3/Documents/data/projects/opentools/target/release:$PATH"

project_id="SR250612M4"
sample_id="ICC"
raw_data_dir="/Volumes/B/data/spatial/raw_data/$project_id"
project_dir="/Volumes/B/openst/opentools/$project_id"
barcode_list_file="/Volumes/B/data/spatial/flowcell/241106/barcodes.txt.gz"
star_index_dir="/Volumes/B/openst/spacemake/species_data/human/genome/star_index"

mkdir -p $project_dir/raw $project_dir/report $project_dir/tile

# trimmed Adapter and PolyA tail
if [[ ! -f "$project_dir/report/cutadapt_trimmed_report.txt" ]]; then
	cutadapt \
		-A "PolyG=G{12};min_overlap=6" \
		-A "PolyA=A{20};min_overlap=10" \
		-A "Read1_Adapter_rev=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT;min_overlap=10" \
		-o $project_dir/raw/Read1_trimmed_adapter_polya.fastq.gz \
		-p $project_dir/raw/Read2_trimmed_adapter_polya.fastq.gz \
		--match-read-wildcards \
		-Q 20 \
		-m 32 \
		-j 10 \
		$raw_data_dir/*R1_001.fastq.gz $raw_data_dir/*R2_001.fastq.gz \
		> $project_dir/report/cutadapt_trimmed_report.txt || \
		rm -f $project_dir/report/cutadapt_trimmed_report.txt
fi

# Convert fastq into bam
if [[ ! -f "$project_dir/unaligned_tagged_trimmed_adapter_polya.bam" ]]; then
	opentools fq2bam \
		-1 $project_dir/raw/Read1_trimmed_adapter_polya.fastq.gz \
		-2 $project_dir/raw/Read2_trimmed_adapter_polya.fastq.gz \
		--sample-id $sample_id > $project_dir/unaligned_tagged_trimmed_adapter_polya.bam || \
		rm -f $project_dir/unaligned_tagged_trimmed_adapter_polya.bam
fi

if [[ ! -f "$project_dir/tile/barcode_whitelist.txt" ]]; then
	opentools dedupbarcode \
		-I $barcode_list_file \
		-o $project_dir/tile \
		--tile-list 22{1..3}0{2..6} || \
		rm -rf $project_dir/tile
fi

# Add file descriptors due to STAR
ulimit -n 5000

# STAR aligned read into genome and aligned read into gene
STAR \
	--runMode alignReads \
	--runThreadN 20 \
	--limitBAMsortRAM 5000000000 \
	--limitOutSJcollapsed 5000000 \
	--readFilesType SAM SE \
	--readFilesCommand samtools view -h \
	--outSAMattributes NH HI AS nM CR CY UR UY CB UB sS \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMprimaryFlag OneBestScore \
	--outMultimapperOrder Random \
	--outFilterType BySJout \
	--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
	--outSAMstrandField intronMotif \
	--soloType CB_UMI_Simple \
	--soloCellFilter None \
	--soloInputSAMattrBarcodeSeq CR UR \
	--soloInputSAMattrBarcodeQual CY UY \
	--soloBarcodeReadLength 0 \
	--soloUMIlen 9 \
	--soloCBlen 28 \
	--soloStrand Forward \
	--soloUMIdedup 1MM_CR \
	--soloFeatures GeneFull_Ex50pAS Velocyto \
	--soloMultiMappers EM \
	--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
	--soloUMIfiltering MultiGeneUMI_CR \
	--soloCellReadStats Standard \
	--genomeDir $star_index_dir \
	--soloCBwhitelist $project_dir/tile/barcode_whitelist.txt \
	--readFilesIn $project_dir/unaligned_tagged_trimmed_adapter_polya.bam \
	--outFileNamePrefix $project_dir/star.genome.