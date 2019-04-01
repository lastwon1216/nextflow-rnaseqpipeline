#!/usr/bin/env nextflow
/*
#$ -v PATH=/common/genomics-core/anaconda2/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin
#$ -cwd
*/

//single_end
//nextflow run rnaseqpipe.nf --genometype 'Mouse' --reads 'test*.fastq.gz' --genomeDir '/common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/reference_genome/' --singleEnd

//pair-end
//nextflow run rnaseqpipe.nf --genometype 'Mouse' --reads 'test{3,4}.fastq.gz' --genomeDir '/common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/reference_genome/'

/*
 * pipeline input parameters
 */
params.genomeDir="$PWD/reference_genome/"
params.outdir = "results"
//params.BAM="/common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/${prefix}_Aligned.toTranscriptome.out.bam"
//params.COR="$PWD/results/${prefix}_Aligned.sortedByCoord.out.bam"
params.reads="$PWD/*.fastq.gz"
params.singleEnd = false
params.genometype=Channel.from( 'Mouse', 'Human')
params.genomeTypeoutdir="Genome"

//BAM=params.BAM
//COR=params.COR
genome_file=params.genomeDir
//i=file(params.reads).simpleName

/*
 * Make Directories for each file type
 */

process mkDIR {
        executor 'sge'
        cpus 8
        penv 'smp'
        memory '6G'
        clusterOptions '-cwd -S /bin/bash -l h=csclprd1-c18v'
        publishDir "${params.outdir}", mode: 'copy'
input:
       // file html from multiqc
        // file "*.out" from start_mkdir
 //       file ('*.sam') from sam_files
output:
file starter into start
script:
"""
touch starter
mkdir /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/node_log
mkdir /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/bam
mkdir /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/fastq
mkdir /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/log
mkdir /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/genes_isoforms_results
mkdir /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/final_results
mkdir /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/others
mkdir /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/fastqc

"""
}

/*
 * Genome Loading
 */

if ( params.genometype == 'Mouse' ){

process mouseGenome{
executor 'sge'
        cpus 8
        penv 'smp'
        memory '4G'
        clusterOptions '-cwd -S /bin/bash -l h=csclprd1-c18v'
        publishDir "${params.outdir}", mode: 'symlink'

input:

file starter from start
output:
file ('*.sam') into sam_files
file "*.out"
file "*STARtmp" into STAR_load
file "QC_reference_files.txt" into RPATH

script:

"""
/common/genomics-core/anaconda2/bin/STAR --genomeDir '/home/choiw1/181227_NEXTFLOW_PIPELINE_WON/nextflow/reference_genome/' --genomeLoad LoadAndExit --outFileNamePrefix Mouse

cp /common/genomics-core/reference/RseQC/GRCm38_m8/QC_reference_files.txt .

"""
        }
}
if (params.genometype == 'Human' ){
process humanGenome{
executor 'sge'
        cpus 2
        penv 'smp'
        memory '4G'
        clusterOptions '-cwd -S /bin/bash -l h=csclprd1-c18v'
        publishDir "${params.outdir}", mode: 'symlink'

input:

output:
file ('*.sam') into sam_files
file "*.out" into start_mkdir
file "*STARtmp" into tmp

script:
"""
/common/genomics-core/anaconda2/bin/STAR --genomeDir '/common/genomics-core/reference/STAR/primary_GRCh38_23' --genomeLoad LoadAndExit --outFileNamePrefix Human
"""
        }
}


/*
 * Alignment
 */


if (params.singleEnd){

Channel
.fromPath( params.reads )
.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
.into { dataset; read_files_fastqc }

process star_SE {
       tag "$prefix"
         executor 'sge'
        cpus 8
        penv 'smp'
        memory '6G'
        clusterOptions '-cwd -S /bin/bash -l h=csclprd1-c18v'
        publishDir "${params.outdir}", mode: 'symlink'

input:
file genome from genome_file
file (reads) from dataset.collect()

file ('*.sam') from sam_files
output:
file ("*Aligned.toTranscriptome.out.bam") into bam_files
file ("*Aligned.sortedByCoord.out.bam") into sorted_bam
file ("*.out") into alignment_logs
file ("*SJ.out.tab") into SJ_tabs
file ("*Log.out") into star_log
script:
//prefix = reads.toString() - ~/(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
'''
for prefix in $(ls *.fastq.gz | rev | cut -c 10- | rev | uniq)
do
/common/genomics-core/anaconda2/bin/STAR --outSAMmode Full --readFilesCommand zcat --outSAMunmapped Within  --outFilterType Normal  --outSAMattributes All --outFilterMultimapNmax 20  --outFilterMismatchNmax 999  --outFilterMismatchNoverLmax 0.04  --alignIntronMin 20  --alignIntronMax 1000000  --alignMatesGapMax 1000000  --alignSJoverhangMin 8  --alignSJDBoverhangMin 1  --sjdbScore 1  --runThreadN 10 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000002 --outSAMtype BAM SortedByCoordinate  --quantMode TranscriptomeSAM --readFilesIn ${prefix}.fastq.gz --genomeDir /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/reference_genome/ --outFileNamePrefix ${prefix}_
done
'''
        }
}
else {
Channel
.fromFilePairs('./*R{1,2}*.fastq.gz')
.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
.into { dataset; read_files_fastqc }


process star_PE {
tag "$prefix"
         executor 'sge'
        cpus 8
        penv 'smp'
        memory '6G'
        clusterOptions '-cwd -S /bin/bash -l h=csclprd1-c18v'
        publishDir "${params.outdir}", mode: 'symlink'

input:
file genome from genome_file
file (reads) from dataset.collect()

file ('*.sam') from sam_files
output:
file ("*Aligned.toTranscriptome.out.bam") into bam_files
file ("*Aligned.sortedByCoord.out.bam") into sorted_bam
file ("*.out") into alignment_logs
file ("*SJ.out.tab") into SJ_tabs
file ("*Log.out") into star_log
script:
//prefix = reads.toString() - ~/(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
'''
for prefix in $(ls *.fastq.gz | rev | cut -c 28- | rev | uniq)
do
/common/genomics-core/anaconda2/bin/STAR --outSAMmode Full --outSAMunmapped Within --readFilesCommand zcat --outFilterType Normal  --outSAMattributes All --outFilterMultimapNmax 20  --outFilterMismatchNmax 999  --outFilterMismatchNoverLmax 0.04  --alignIntronMin 20  --alignIntronMax 1000000  --alignMatesGapMax 1000000  --alignSJoverhangMin 8  --alignSJDBoverhangMin 1  --sjdbScore 1  --runThreadN 3 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000002 --outSAMtype BAM SortedByCoordinate  --quantMode TranscriptomeSAM --readFilesIn ${prefix}_R{1,2}*.fastq.gz --genomeDir $genome_file --outFileNamePrefix ${prefix}_
done
'''
}
}


/*
 * Genome Unloading
 */


if ( params.genometype == 'Mouse' ){

process removeMouseGenome{
        executor 'sge'
        cpus 2
        penv 'smp'
        memory '4G'
        clusterOptions '-cwd -S /bin/bash -l h=csclprd1-c18v'

input:
file ("*.out") from alignment_logs
output:
file "*.out"

script:

"""
/common/genomics-core/anaconda2/bin/STAR --genomeDir '/home/choiw1/181227_NEXTFLOW_PIPELINE_WON/nextflow/reference_genome/' --genomeLoad Remove --outFileNamePrefix Mouse

echo "ref genome has been removed"
"""
        }
}
if (params.genometype == 'Human' ){
process removeHumanGenome{
executor 'sge'
        cpus 2
        penv 'smp'
        memory '4G'
        clusterOptions '-cwd -S /bin/bash'

input:
file ("*.out") from alignment_logs
output:
file "*.out"
file "*STARtmp"

script:
"""
/common/genomics-core/anaconda2/bin/STAR --genomeDir '/common/genomics-core/reference/STAR/primary_GRCh38_23' --genomeLoad Remove --outFileNamePrefix Mouse

echo "ref genome has been removed"
"""
        }
}


/*
 * FastQC
 */


process FastQC {
        tag "$prefix"
        executor 'sge'
        cpus 8
        penv 'smp'
        memory '6G'
        clusterOptions '-cwd -S /bin/bash -l h=csclprd1-c18v'
        publishDir "${params.outdir}", mode: 'symlink'
        beforeScript 'module load R samtools/1.6'
input:
        file(reads) from read_files_fastqc
file ("*SJ.out.tab") from SJ_tabs
output:
      file "*fastqc.html"
       file "*fastqc.zip" into fastqc
script:
prefix = reads.toString() - ~/(_R1)?(_R1)?(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
"""
/common/genomics-core/anaconda2/bin/fastqc ${prefix}.fastq.gz
"""
}


/*
 * RSEQC
 */


process rseqc {
        executor 'sge'
        cpus 8
        penv 'smp'
        memory '6G'
        clusterOptions '-cwd -S /bin/bash -l h=csclprd1-c18v'
        publishDir "${params.outdir}", mode: 'symlink'
        beforeScript 'module load R samtools/1.6'
input:
  //      file txt from rsem_summary
file zips from fastqc    
    file REFPATH from RPATH
       // file ("*Aligned.sortedByCoord.out.bam") from sorted_bam
        file bam from sorted_bam
output:
        file("${prefix}") into rseQc
script:
  if( params.singleEnd ) {
//prefix = txt[0].toString() - '.COUNTS.txt' - '.FPKM.txt' - '.TPM.txt'
prefix = zips.toString() - "_fastqc.zip"
COR = "$PWD/results/${prefix}_Aligned.sortedByCoord.out.bam"
         """
        /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/QC/Rse_QC_pipeline.pl -i ${COR} -t SE -d . -p ${prefix} -r $REFPATH

        echo "Done for RseQC!"

        """
    }
    else {
       """

        /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/QC/Rse_QC_pipeline.pl -i $COR -t PE -d ./QC/RseQC_results -p test -r $REFPATH

        echo "Done for RseQC!"

        """
    }
}


/*
 * RSEM
 */


process rsem {
        executor 'sge'
        cpus 8
        penv 'smp'
        memory '6G'
        clusterOptions '-cwd -S /bin/bash -l h=csclprd1-c18v'
        publishDir "${params.outdir}", mode: 'symlink'
input:
//file "*.out" from UnloadGen
//file "*STARtmp" from UnloadGen
 //   file("*.zip") from fastqc
file ("${prefix}_Aligned.toTranscriptome.out.bam") from bam_files
file txt from rseQc
//file zips from fastqc
//file "*fastqc.zip" from fastqc
output:
       file("*.time") into star_time
       file("*.results") into star_results
       file("*.stat") into star_stat
//       file( "${prefix}.genes.results" ) into rsem_results
       file( "*.genes.results" ) into rsem_results

script:
prefix = txt.toString() - '.final_results_summary.txt'
//prefix = zips.toString() - "_fastqc.zip" 
BAM = "/common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/${prefix}_Aligned.toTranscriptome.out.bam" 
  if( params.singleEnd ) {
"""
        rsem-calculate-expression --no-bam-output --append-names --bam -p 5 --time ${BAM} $genome_file ${prefix}
"""   
 }
    else {
        """
        rsem-calculate-expression --no-bam-output --paired-end --append-names --bam -p 5 --time $BAM $genome_file ${prefix}
"""
    }
}


/*
 * RSEM Summary
 */


process rsemSummary {
        executor 'sge'
        cpus 8
        penv 'smp'
        memory '6G'
        clusterOptions '-cwd -S /bin/bash -l h=csclprd1-c18v'
        publishDir "${params.outdir}", mode: 'symlink'
        beforeScript 'module load R samtools/1.6'
input:
//file( "*.genes.results" ) from rsem_results
file results from rsem_results
output:
file("*.txt") into rsem_summary
file ("*.csv")
script:
prefix = results.toString() - ".genes.results"
"""
/common/genomics-core/anaconda2/bin/rsem-generate-data-matrix /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/${prefix}.genes.results > /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/${prefix}.COUNTS.txt

/common/genomics-core/anaconda2/bin/rsem-generate-data-tpm /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/${prefix}.genes.results > /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/${prefix}.TPM.txt

/common/genomics-core/anaconda2/bin/rsem-generate-data-fpkm /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/${prefix}.genes.results > /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/${prefix}.FPKM.txt

/hpc/apps/R/3.4.1/bin/Rscript /common/genomics-core/apps/sequencing_data_distri/format.R /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/${prefix}.TPM.txt /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/${prefix}.COUNTS.txt /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/${prefix}.FPKM.txt ${prefix}

cp /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/${prefix}.*.txt .
"""
}


/*
 * RseQC Summary
 */


process rseQcSummary {
        executor 'sge'
        cpus 8
        penv 'smp'
        memory '6G'
        clusterOptions '-cwd -S /bin/bash -l h=csclprd1-c18v'
        publishDir "${params.outdir}", mode: 'symlink'
        beforeScript 'module load R samtools/1.6'
input:
       // file txt from rseQc
	file txt from rsem_summary
output:
        file("*txt") into rseqc_summary_txt
        file("*pdf") into rseqc_summary
//      file("test")
script:
//prefix = txt.toString() - '.final_results_summary.txt'
prefix = txt[0].toString() - '.COUNTS.txt' - '.FPKM.txt' - '.TPM.txt'
"""
Rscript /common/genomics-core/bin/Rse_QC_summary.R /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*/*.final_results_summary.txt all
"""
    }


/*
 * MultiQC
 */


process MultiQC {
        executor 'sge'
        cpus 8
        penv 'smp'
        memory '6G'
        clusterOptions '-cwd -S /bin/bash -l h=csclprd1-c18v'
        publishDir "${params.outdir}", mode: 'symlink'
input:
        //file txt from rseQc       
 file pdf from rseqc_summary
output:
        file("*.html") into multiqc
        file("*.multiqc_data")
	file starter into organize
script:
prefix = pdf.toString() - '_QC.pdf'
//prefix = txt.toString() - '.final_results_summary.txt'
"""
/common/genomics-core/anaconda2/bin/multiqc -n ${prefix}.multiqc -f -c /common/genomics-core/apps/multiqc/multiqc_config_RNAseq_HPC.yaml -e rsem /common/genomics-core/bin/Rse_QC_summary.R /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/

cp -f /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/starter .
"""
 }


/*
 * Directory Organizer
 */


process organizeDIR {
        executor 'sge'
        cpus 8
        penv 'smp'
        memory '6G'
        clusterOptions '-cwd -S /bin/bash -l h=csclprd1-c18v'
        publishDir "${params.outdir}", mode: 'symlink'
input:
file starter from organize
// file html from multiqc 
//       val "*results" from mkdir
output:
file starter into email
script:
"""
cp /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*_fastqc.* /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/fastqc
cp /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*.csv /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*QC.txt /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*QC.pdf /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*.multiqc.html /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/final_results/
cp /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*.genes.results /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/genes_isoforms_results
cp /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*.isoforms.results /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/genes_isoforms_results
cp /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*.bam /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/bam
cp /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*Log* /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*.time /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/log
cp /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*tab /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*.COUNTS.txt /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*.FPKM.txt /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*.TPM.txt /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/others
cp -rf /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*_STARtmp /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/*.stat /common/genomics-core/data/Internal_Tests/181227_NEXTFLOW_PIPELINE_WON/nextflow/results/others
"""
}


/*
 * Send Email
 */


process sendemail {
        executor 'sge'
        cpus 8
        penv 'smp'
        memory '6G'
        clusterOptions '-cwd -S /bin/bash -l h=csclprd1-c18v'
        publishDir "${params.outdir}", mode: 'symlink'
input:
file starter from email
// file html from multiqc
//       val "*results" from mkdir
output:
script:
""" 
echo "Subject: Mapping/QC is done" | sendmail -v wonje.choi@cshs.org

"""
}