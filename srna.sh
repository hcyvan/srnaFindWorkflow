#!/usr/bin/env bash

# 1) ./srna.sh fastqc
# 1.1) ./srna.sh multiqc
# 2) ./srna.sh trimmoatic
# 3) ./srna.sh fastqc --trim
# 3.1) ./srna.sh multiqc
# 4) ./srna.sh gindex
# 4.1) ./srna.sh alignment
# 4.2) ./srna.sh sam2bam
# 5) ./srna.sh bam2bedgraph
########################################
# 6) ./srna.sh bam_size > bam_size.txt
# 7) create bedgraph.txt
# #./srna.py  scale -i ./bedgraph -o ./bedgraph -f ./bedgraph.txt
# ./srna.py find -g ./genome/GCF_000195955.2_ASM19595v2_genomic.gff -i ./bedgraph -o ./bedgraph -b ./sample.txt -c 10
# ./srna.py merge -co sample.txt -i ./bedgraph > ./gff/srna.gff
# ./srna.py merge2anno -f ./gff/srna.gff -a ./gff/GCF_000195955.2_ASM19595v2_genomic.gff |grep SRNA > ./gff/srna.nctype.gff

trimmed=false

case "$#:$1" in
    1:test)
        mode=run_test
        ;;
    *:fastqc)
        if [[ "$2" == "--trim" ]]; then
            trimmed=true
        fi
        mode=run_fastqc
        ;;
    *:multiqc)
        if [[ "$2" == "--trim" ]]; then
            trimmed=true
        fi
        mode=run_multiqc
        ;;
    1:trimmoatic)
        mode=run_trimmoatic
        ;;
    1:alignment)
        mode=run_alignment
        ;;
    1:gindex)
        mode=run_gindex
        ;;
    1:sam2bam)
        mode=run_sam2bam
        ;;
    1:bam2bedgraph)
        mode=run_bam2bedgraph
        ;;
    1:bam_size)
        mode=run_bam_size
        ;;

    *)
        echo >&2 "Usage: $0 [--verify|--list|--no-cipd]"
        exit 1
        ;;
esac

readonly SAMPLES_FILE=./sample.txt
readonly GENOME_DIR=./genome
readonly FASTQ_DIR=./fastq
readonly TRIM_DIR=./fastq_trim
readonly QC_DIR=./qc
readonly SAM_DIR=./sam
readonly BAM_DIR=./bam
readonly BEDGRAPH_DIR=./bedgraph
readonly TRIMMOMATIC_PATH=/home/c509/bin/Trimmomatic-0.38

.bam_size() {
    local sample="$1"
    size=$(samtools view -c -F4 ./bam/${sample}.bam)
    echo -e "$1\t$size"
}

.bam2bedgraph() {
    local sample="$1"
    bedtools genomecov  -strand + -d -ibam ./bam/${sample}.bam > ${BEDGRAPH_DIR}/${sample}.1.bedgraph
    bedtools genomecov  -strand - -d -ibam ./bam/${sample}.bam > ${BEDGRAPH_DIR}/${sample}.2.bedgraph
}

.sam2bam() {
    local sample="$1"
    samtools sort -@ 8 -o ${BAM_DIR}/${sample}.bam ${SAM_DIR}/${sample}.sam
}

.alignment() {
    local sample="$1"
    hisat2 -p 4 --new-summary --summary-file ${SAM_DIR}/${sample}.summary \
           -x ${GENOME_DIR}/GCF_000195955.2_ASM19595v2_genomic \
           -1 ${TRIM_DIR}/${sample}_1_paired.fq.gz \
           -2 ${TRIM_DIR}/${sample}_2_paired.fq.gz \
           -S ${SAM_DIR}/${sample}.sam
}

.trimmomatic() {
    local sample="$1"
    java -jar ${TRIMMOMATIC_PATH}/trimmomatic-0.38.jar PE \
         ${FASTQ_DIR}/${sample}_1.fq.gz ${FASTQ_DIR}/${sample}_2.fq.gz \
         ${TRIM_DIR}/${sample}_1_paired.fq.gz ${TRIM_DIR}/${sample}_1_unpaired.fq.gz \
         ${TRIM_DIR}/${sample}_2_paired.fq.gz ${TRIM_DIR}/${sample}_2_unpaired.fq.gz \
         ILLUMINACLIP:${TRIMMOMATIC_PATH}/adapters/TruSeq3-PE.fa:2:30:10 \
         LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
}

.fastqc() {
    local sample="$1"
    if $trimmed; then
        fastqc -t 8 -o ${QC_DIR}/fastqc_trim/ \
               ${TRIM_DIR}/${sample}_1_paired.fq.gz \
               ${TRIM_DIR}/${sample}_2_paired.fq.gz
    else
        fastqc -t 8 -o ${QC_DIR}/fastqc/ ${FASTQ_DIR}/${sample}_1.fq.gz ${FASTQ_DIR}/${sample}_2.fq.gz
    fi
}

.test() {
    echo $1
}

run_bam_size() {
    for_each_sample .bam_size
}

run_bam2bedgraph() {
    for_each_sample .bam2bedgraph
}

run_sam2bam() {
    for_each_sample .sam2bam
}

run_alignment() {
    for_each_sample .alignment
}

run_gindex() {
    hisat2-build -p 4 ${GENOME_DIR}/GCF_000195955.2_ASM19595v2_genomic.fna \
                 ${GENOME_DIR}/GCF_000195955.2_ASM19595v2_genomic
}

run_trimmoatic() {
    for_each_sample .trimmomatic
}

run_multiqc() {
    if $trimmed; then
        multiqc -f ${QC_DIR}/fastqc_trim/*_fastqc.zip -o ${QC_DIR} -n multiqc_report_trimmed.html
    else
        multiqc -f ${QC_DIR}/fastqc/*_fastqc.zip -o ${QC_DIR} -n multiqc_report.html
    fi
}

run_fastqc() {
    for_each_sample .fastqc
}

run_test() {
    for_each_sample .test
}

for_each_sample() {
    while read line; do
        echo "$line" >&2
        "$@" "$line"
    done < "$SAMPLES_FILE"
}

$mode
