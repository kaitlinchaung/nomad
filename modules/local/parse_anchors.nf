
process PARSE_ANCHORS {

    input:
    tuple val(fastq_id), path(fastq)
    path anchors
    val num_lines
    val consensus_length
    val kmer_size
    val direction
    val read_len

    output:
    path "*_target_counts.tsv"  , emit: targets
    path "*.fasta"              , emit: consensus_fasta
    path "*.tab"                , emit: consensus_stats

    script:
    out_consensus_fasta_file    = "${fastq_id}_${direction}.fasta"
    out_counts_file             = "${fastq_id}_${direction}_counts.tab"
    out_fractions_file          = "${fastq_id}_${direction}_fractions.tab"
    out_adj_kmer_file           = "${fastq_id}_target_counts.tsv"
    """
    parse_anchors.py \\
        --num_input_lines ${num_lines} \\
        --anchors_file ${anchors} \\
        --fastq_id ${fastq_id} \\
        --fastq_file ${fastq} \\
        --out_consensus_fasta_file ${out_consensus_fasta_file} \\
        --out_counts_file ${out_counts_file} \\
        --out_fractions_file ${out_fractions_file} \\
        --out_adj_kmer_file ${out_adj_kmer_file} \\
        --consensus_length ${consensus_length} \\
        --kmer_size ${kmer_size} \\
        --direction ${direction} \\
        --read_len ${read_len}
    """
}
