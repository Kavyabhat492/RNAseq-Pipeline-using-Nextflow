
nextflow.enable.dsl=2

process FASTP_TRIM {

    publishDir "TRIMMED", mode:'copy'
    conda 'bioconda::fastp=1.0.1' 

    input:
        tuple val(sampleid), path(reads) 

    output:
        path "${sampleid}_R1.trimmed.fq.gz", emit: trimmed_R1
        path "${sampleid}_R2.trimmed.fq.gz", emit: trimmed_R2
        path "*" 

    script:
    """
    fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 ${sampleid}_R1.trimmed.fq.gz --out2 ${sampleid}_R2.trimmed.fq.gz --detect_adapter_for_pe --qualified_quality_phred 15 --unqualified_percent_limit 20 --trim_front1 9 --trim_front2 9 --json ${sampleid}.fastp.json --html ${sampleid}.fastp.html

    """
}

process FASTQC_TRIMMED {

    publishDir "QC_TRIMMED", mode: 'copy'
    conda 'bioconda::fastqc=0.12.1'
    
    cpus 2 
    memory '4 GB'

    input:
        tuple val(sampleid), path(r1), path(r2) 

    output:
        path "*_fastqc.html", emit: html_reports
        path "*_fastqc.zip", emit: zip_files

    script:
    """
    fastqc -o . -f fastq -t ${task.cpus} ${r1} ${r2}
    """
}

process HISAT2_INDEX {

    publishDir "HISAT2_INDEX", mode: 'copy'
    conda 'bioconda::hisat2=2.2.1'

    cpus 1 
    memory '12 GB'

    input:
        path(fasta) 

    output:
        path "genome.*.ht2", emit: hisat2_index_base

    script:
    """
    hisat2-build -p 1 ${fasta} genome
    """
}


process ALIGN_AND_COUNT {

    publishDir "COUNTS", mode: 'copy'
    
    conda 'bioconda::hisat2=2.2.1 bioconda::samtools=1.18 bioconda::subread=2.0.4'
    cpus 4
    memory '16 GB'

    input:
        tuple val(sampleid), 
              path(trimmed_R1), 
              path(trimmed_R2), 
              path(hisat2_index_base), 
              path(ref_gtf), 
              val(strand)

    output:
        path "${sampleid}.Aligned.sorted.bam", emit: aligned_bam
        path "${sampleid}.counts.txt", emit: gene_counts

    script:
    """
    

    # 1. ALIGNMENT (HISAT2)
    hisat2 -x hisat2_index_base \\
        -1 ${trimmed_R1} \\
        -2 ${trimmed_R2} \\
        --dta \\
        -p 4 | samtools view -bS
        
    # 2. SORT BAM
    samtools sort -@ 4 aligned.unsorted.bam -o ${sampleid}.Aligned.sorted.bam

    # 3. FEATURE_COUNTS (Quantification)
    featureCounts \\
        -p \\
        -t exon \\
        -g gene_name \\
        -a ${ref_gtf} \\
        -o ${sampleid}.counts.txt \\
        -s ${strand} \\
        -T ${task.cpus} \\
        ${sampleid}.Aligned.sorted.bam
    """
}


workflow {

    // Define input channels
    ref_fasta = Channel.fromPath(params.ref_fasta)
    ref_gtf = Channel.fromPath(params.ref_gtf)
    fastq_ch = Channel.fromFilePairs(params.reads)
    strand = Channel.of(params.strand)

    // 1. Trimming (Structural call - will be SKIPPED)
    FASTP_TRIM(fastq_ch).set{ trimmed }

    // 2. QC (Structural call - will be SKIPPED)
    fastqc_input_ch = trimmed.trimmed_R1.join(trimmed.trimmed_R2)
    FASTQC_TRIMMED(fastqc_input_ch).set { fastqc_reports }

    // 3. Indexing (Structural call - will be SKIPPED)
    HISAT2_INDEX(ref_fasta).set { hisat2_index_base }

    // 4. Alignment and Counting (The first incomplete job, where execution will start)
    final_alignment_channel = fastqc_input_ch 
        .combine(hisat2_index_base)
        .combine(ref_gtf)
        .combine(strand)

    ALIGN_AND_COUNT(final_alignment_channel).set { aligned_and_counted }

}
 
