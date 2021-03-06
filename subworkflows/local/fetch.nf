include { GET_UNMAPPED              } from '../../modules/local/get_unmapped'
include { TRIMGALORE                } from '../../modules/nf-core/modules/trimgalore/main'
include { FETCH_ANCHORS             } from '../../modules/local/fetch_anchors'
include { COUNT_ANCHORS             } from '../../modules/local/count_anchors'
include { STRATIFY_ANCHORS          } from '../../modules/local/stratify_anchors'
include { GET_ABUNDANT_ANCHORS      } from '../../modules/local/get_abundant_anchors'
include { MERGE_ABUNDANT_ANCHORS    } from '../../modules/local/merge_abundant_anchors'
include { COMPUTE_PVALS             } from '../../modules/local/compute_pvals'
include { SIGNIFICANT_ANCHORS       } from '../../modules/local/significant_anchors'


workflow FETCH {

    take:
    ch_fastqs
    lookahead

    main:

    if (params.run_trimming) {
        /*
        // Trim fastqs
        */
        TRIMGALORE(
            ch_fastqs
        )
        ch_fastqs = TRIMGALORE.out.fastq
    }

    /*
    // Check if we are only using unmapped reads
    */
    if (params.run_get_unmapped) {
        /*
        // Get unmapped reads
        */
        GET_UNMAPPED(
            ch_fastqs,
            params.index_bowtie
        )

        ch_fastqs = GET_UNMAPPED.out.fastq

    }

    /*
    // Process to get all candidate anchors and targets
    */
    FETCH_ANCHORS(
        ch_fastqs,
        params.run_type,
        params.num_reads_first_pass,
        params.kmer_size,
        lookahead,
        params.anchor_mode,
        params.window_slide

    )

    /*
    // Process to count anchors and targets
    */
    COUNT_ANCHORS(
        FETCH_ANCHORS.out.seqs
    )

    /*
    // Process to stratify counts files by 3mers
    */
    STRATIFY_ANCHORS(
        COUNT_ANCHORS.out.seqs.collect(),
        params.stratify_level
    )

    /*
    // Process to filter kmer counts for abundant anchors
    */
    GET_ABUNDANT_ANCHORS(
        STRATIFY_ANCHORS.out.seqs.flatten(),
        params.anchor_count_threshold,
        params.kmer_size
    )

    if (params.run_decoy) {
        MERGE_ABUNDANT_ANCHORS(
            GET_ABUNDANT_ANCHORS.out.anchor_counts.collect(),
            params.num_decoy_anchors
        )

        anchors_scores = MERGE_ABUNDANT_ANCHORS.out.seqs


    } else {
        /*
        // Process to get significant anchors and their scores
        */
        COMPUTE_PVALS(
            GET_ABUNDANT_ANCHORS.out.seqs,
            params.kmer_size,
            file(params.input),
            params.K_num_hashes,
            params.L_num_random_Cj
        )

        /*
        // Process to output top 5000 anchors as sorted by pvalue
        */
        SIGNIFICANT_ANCHORS(
            COMPUTE_PVALS.out.scores.collect(),
            params.fdr_threshold,
            params.all_anchors_pvals_file
        )

        anchors_scores = SIGNIFICANT_ANCHORS.out.scores
            .filter{
                it.countLines() > 1
            }

    }

    emit:
    anchors_scores = anchors_scores

}
