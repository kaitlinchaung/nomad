
process COUNT_ANCHORS {

    label 'process_medium'

    input:
    tuple val(fastq_id), path(fastq), val(group_id)

    output:
    tuple val(fastq_id), path(counts)

    script:
    counts = "counted_${fastq_id}.txt"
    """
    zcat ${fastq} \\
        | sort \\
        | uniq -c > ${counts}
    """
}
