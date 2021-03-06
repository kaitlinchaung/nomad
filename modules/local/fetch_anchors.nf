
process FETCH_ANCHORS {

    tag "${fastq_id}"
    label 'process_low'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1 numpy=1.22.3" : null)

    input:
    tuple val(fastq_id), path(fastq), val(group_id)
    val run_type
    val num_reads_first_pass
    val kmer_size
    val lookahead
    val anchor_mode
    val window_slide

    output:
    tuple val(fastq_id), path(outfile), val(group_id), emit: seqs

    script:
    outfile         = "sequences_${fastq_id}.txt.gz"
    """
    fetch_anchors.py \\
        --infile ${fastq} \\
        --run_type ${run_type} \\
        --fastq_id ${fastq_id} \\
        --num_lines ${num_reads_first_pass} \\
        --kmer_size ${kmer_size} \\
        --lookahead ${lookahead} \\
        --anchor_mode ${anchor_mode} \\
        --window_slide ${window_slide} \\
        --outfile ${outfile}
    """
}
