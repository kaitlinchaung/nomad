include { MERGE_ABUNDANT_ANCHORS    } from '../../modules/local/merge_abundant_anchors'
include { GET_ANCHORS_AND_SCORES    } from '../../modules/local/get_anchors_and_scores'
include { PARSE_ANCHORS             } from '../../modules/local/parse_anchors'
include { MERGE_TARGET_COUNTS       } from '../../modules/local/merge_target_counts'

workflow ANALYZE_BULK {

    take:
    ch_fastqs
    lookahead
    abundant_anchors

    main:

    /*
    // Process to merge abundant anchors
    */
    MERGE_ABUNDANT_ANCHORS(
        GET_ABUNDANT_ANCHORS.out.seqs.collect()
    )

    /*
    // Process to get significant anchors and their scores
    */
    GET_ANCHORS_AND_SCORES(
        MERGE_ABUNDANT_ANCHORS.out.seqs,
        file(params.input),
        params.distance_type,
        params.max_targets,
        params.max_dist,
        params.bonfer,
        params.pval_threshold,
        params.run_type
    )

    ch_anchors  = GET_ANCHORS_AND_SCORES.out.anchors.filter{ it.size() > 0 }
    ch_scores   = GET_ANCHORS_AND_SCORES.out.scores

    /*
    // Process to get consensus sequences and target counts for annchors
    */
    PARSE_ANCHORS(
        ch_fastqs,
        ch_anchors,
        params.num_parse_anchors_reads,
        params.consensus_length,
        params.kmer_size,
        params.direction,
        lookahead
    )

    // Create samplesheet of target counts files
    PARSE_ANCHORS.out.targets
        .collectFile() { file ->
            def X=file; X.toString() + '\n'
        }
        .set{ targets_samplesheet }

    /*
    // Process to get anchor scores and anchor-target counts
    */
    MERGE_TARGET_COUNTS(
        targets_samplesheet
    )

    emit:
    anchor_target_counts = MERGE_TARGET_COUNTS.out.anchor_target_counts.first()
    anchor_scores        = ch_scores

}
