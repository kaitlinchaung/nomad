
process STRATIFY_ANCHORS {

    label 'process_high'
    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)

    input:
    tuple val(fastq_id), path(counts)
    val stratify_level

    output:
    path("stratified_*"), emit: seqs

    script:
    """
    stratify_anchors.py \\
        --infile \${counts} \\
        --stratify_level ${stratify_level}
    """
}
