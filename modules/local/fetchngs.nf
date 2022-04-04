
process fetchngs {

    label 'process_low'
    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)


    input:
    path ch_ids
    path outdir

    output:
    path outfile_counts_distances   , emit: anchor_target_counts
    path outfile_anchor_scores      , emit: anchor_scores

    script:
    outfile_counts_distances        = "anchor_targets_counts.tsv"
    outfile_anchor_scores           = "anchor_scores.tsv"
    """
    nextflow run kaitlinchaung/fetchngs --input ${ch_ids} --outdir ${outdir}
    """
}
