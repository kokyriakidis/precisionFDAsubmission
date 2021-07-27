## MANUALLY SET THESE VARIABLES ##
WORKDIR="/home/kokyriakidis/Desktop/Submissions"
SAMPLENAME="HG002"
SUBMISSION="TZMTX"
FQ1FILEPATH="path/to/HG002.novaseq.pcr-free.35x.R1.fastq.gz"
FQ2FILEPATH="path/to/HG002.novaseq.pcr-free.35x.R2.fastq.gz"
NTHREADS=18


## SUBMISSION PIPELINE ##

mkdir -p ${WORKDIR}

cd ${WORKDIR}

mkdir -p ${SUBMISSION}

cd ${SUBMISSION}

mkdir -p reference

cd reference


# Download the full hg38/GRCh38 reference genome distributed by 1000 genomes
# Derived from NCBI set with HLA and decoy alternative alleles
# http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/url=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome
base=GRCh38_full_analysis_set_plus_decoy_hla
new=hg38
mkdir -p seq
for suffix in .dict .fa.fai
do
    [[ -f seq/$new$suffix ]] || wget -c -O seq/$new$suffix $url/$base$suffix
done
# gzipped references for use with ensembl-vep HGVS
urlgz=https://s3.amazonaws.com/biodata/hg38_bundle
for suffix in .fa.gz .fa.gz.fai .fa.gz.gzi
do
    [[ -f seq/$new$suffix ]] || wget --no-check-certificate -c -O seq/$new$suffix $urlgz/$new$suffix
done
gunzip -c seq/hg38.fa.gz > seq/hg38.fa
touch seq/hg38.fa.fai
touch seq/hg38.dict


## Download pre-build bwa indices ##
url=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome
base=GRCh38_full_analysis_set_plus_decoy_hla.fa
new=hg38.fa
mkdir -p bwa
for suffix in .bwt .amb .ann .pac .sa .alt
do
    [[ -f bwa/$new$suffix ]] || wget -c -O bwa/$new$suffix $url/$base$suffix
done


cd ${WORKDIR}

## ADAPTER TRIMMING WITH ATROPOS ##
atropos trim \
-a AACACTCTTTCCCT -A AGATCGGAAGAGCG \
--quality-base 64 --format fastq --overlap 8 \
--aligner adapter --quality-cutoff 5 \
-o ${WORKDIR}/${SUBMISSION}/ATROPOS/${SAMPLENAME}-${SUBMISSION}-trimmed.1.fq.gz -p ${WORKDIR}/${SUBMISSION}/ATROPOS/${SAMPLENAME}-${SUBMISSION}-trimmed.2.fq.gz \
-pe1 ${FQ1FILEPATH} -pe2 ${FQ2FILEPATH}


mkdir -p ${WORKDIR}/${SUBMISSION}/BWA
mkdir -p ${WORKDIR}/${SUBMISSION}/BWA/tmpdir


## READ MAPPING WITH BWA ##
bwa mem \
-c 250 \
-M \
-t ${NTHREADS} \
-R '@RG\tID:${SAMPLENAME}\tPL:illumina\tPU:${SAMPLENAME}\tSM:${SAMPLENAME}' \
-v 1 \
${WORKDIR}/${SUBMISSION}/reference/bwa/hg38.fa \
${WORKDIR}/${SUBMISSION}/ATROPOS/${SAMPLENAME}-${SUBMISSION}-trimmed.1.fq.gz \
${WORKDIR}/${SUBMISSION}/ATROPOS/${SAMPLENAME}-${SUBMISSION}-trimmed.2.fq.gz > ${WORKDIR}/${SUBMISSION}/BWA/${SAMPLENAME}-${SUBMISSION}-sorted.sam

## SORT BAM FILE BY CHROMOSOMAL COORDINATES ##
samtools sort -@ ${NTHREADS} -m 3G -O bam \
-T ${WORKDIR}/${SUBMISSION}/BWA/tmpdir \
-o ${WORKDIR}/${SUBMISSION}/BWA/${SAMPLENAME}-${SUBMISSION}-sorted.bam \
${WORKDIR}/${SUBMISSION}/BWA/${SAMPLENAME}-${SUBMISSION}-sorted.sam


## INDEX BAM FILE WITH SAMTOOLS ##
samtools index -@ ${NTHREADS} ${WORKDIR}/${SUBMISSION}/BWA/${SAMPLENAME}-${SUBMISSION}-sorted.bam

## REMOVE INTERMEDIATE FILES ##
rm ${WORKDIR}/${SUBMISSION}/BWA/${SAMPLENAME}-${SUBMISSION}-sorted.sam

rm -rf ${WORKDIR}/${SUBMISSION}/BWA/tmpdir

mkdir ${WORKDIR}/${SUBMISSION}/PICARD

## MARK DUPLICATES WITH PICARD ##
gatk MarkDuplicates \
    --INPUT ${WORKDIR}/${SUBMISSION}/BWA/${SAMPLENAME}-${SUBMISSION}-sorted.bam \
    --OUTPUT ${WORKDIR}/${SUBMISSION}/PICARD/${SAMPLENAME}-${SUBMISSION}-sorted-markdup.bam \
    --METRICS_FILE ${WORKDIR}/${SUBMISSION}/PICARD/${SAMPLENAME}-${SUBMISSION}-markdup-metrics.txt \
    --ASSUME_SORT_ORDER coordinate \
    --VALIDATION_STRINGENCY STRICT

## INDEX BAM FILE WITH SAMTOOLS ##
samtools index -@ ${NTHREADS} ${WORKDIR}/${SUBMISSION}/PICARD/${SAMPLENAME}-${SUBMISSION}-sorted-markdup.bam


mkdir ${WORKDIR}/${SUBMISSION}/DEEPVARIANT
mkdir ${WORKDIR}/${SUBMISSION}/DEEPVARIANT/tmpdir
mkdir ${WORKDIR}/${SUBMISSION}/DEEPVARIANT/logdir


## VARIANT CALLING WITH DEEPVARIANT ##
dv_make_examples.py \
--ref ${WORKDIR}/${SUBMISSION}/reference/seq/hg38.fa \
--reads ${WORKDIR}/${SUBMISSION}/PICARD/${SAMPLENAME}-${SUBMISSION}-sorted-markdup.bam \
--examples ${WORKDIR}/${SUBMISSION}/DEEPVARIANT/tmpdir \
--sample ${SAMPLENAME} \
--cores ${NTHREADS} \
--logdir ${WORKDIR}/${SUBMISSION}/DEEPVARIANT/logdir


dv_call_variants.py \
--cores ${NTHREADS} \
--outfile ${WORKDIR}/${SUBMISSION}/DEEPVARIANT/tmpdir/${SAMPLENAME}.tmp \
--sample ${SAMPLENAME} \
--examples ${WORKDIR}/${SUBMISSION}/DEEPVARIANT/tmpdir \
--model wgs

mkdir ${WORKDIR}/${SUBMISSION}/DEEPVARIANT/results

dv_postprocess_variants.py \
--ref ${WORKDIR}/${SUBMISSION}/reference/seq/hg38.fa \
--infile ${WORKDIR}/${SUBMISSION}/DEEPVARIANT/tmpdir/${SAMPLENAME}.tmp \
--outfile ${WORKDIR}/${SUBMISSION}/DEEPVARIANT/results/${SAMPLENAME}.vcf














