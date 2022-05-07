
process PARSE_ANCHORS {

    tag "${fastq_id}"
    label 'process_medium'
    conda (params.enable_conda ? "conda-forge::python=3.9.5 pandas=1.4.1" : null)

    input:
    tuple val(fastq_id), path(fastq), val(group_id)
    path anchors
    val num_parse_anchors_reads
    val consensus_length
    val kmer_size
    val direction
    val lookahead

    output:
    path "*_target_counts.tsv"  , emit: targets               , optional: true
    path "*.fasta"              , emit: consensus_fasta
    path "*.tab"                , emit: consensus_stats
    path "*log"                 , emit: log

    script:
    out_consensus_fasta_file    = "${fastq_id}_${direction}.fasta"
    out_counts_file             = "${fastq_id}_${direction}_counts.tab"
    out_fractions_file          = "${fastq_id}_${direction}_fractions.tab"
    out_target_file             = "${fastq_id}_target_counts.tsv"
    """
    parse_anchors.py \\
        --num_parse_anchors_reads ${num_parse_anchors_reads} \\
        --anchors_file ${anchors} \\
        --fastq_id ${fastq_id} \\
        --fastq_file ${fastq} \\
        --out_consensus_fasta_file ${out_consensus_fasta_file} \\
        --out_counts_file ${out_counts_file} \\
        --out_fractions_file ${out_fractions_file} \\
        --out_target_file ${out_target_file} \\
        --consensus_length ${consensus_length} \\
        --kmer_size ${kmer_size} \\
        --direction ${direction} \\
        --lookahead ${lookahead}
    """
}
