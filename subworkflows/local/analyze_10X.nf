include { GET_ANCHORS_AND_SCORES    } from '../../modules/local/get_anchors_and_scores'
include { PARSE_ANCHORS             } from '../../modules/local/parse_anchors'
include { MERGE_TARGET_COUNTS       } from '../../modules/local/merge_target_counts'

workflow ANALYZE_FASTQS_10X {

    take:
    ch_fastqs
    lookahead
    abundant_anchors

    main:

    /*
    // Process to get significant anchors and their scores
    */
    GET_ANCHORS_AND_SCORES(
        abundant_anchors,
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
        lookahead,
        
    )

    

    emit:
    anchor_target_counts = MERGE_TARGET_COUNTS.out.anchor_target_counts.first()
    anchor_scores        = ch_scores

}
