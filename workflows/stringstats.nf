/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowStringstats.initialise(params, log)


/*
========================================================================================
    CONFIG FILES
========================================================================================
*/


/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS     } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { GET_READ_LENGTH           } from '../modules/local/get_read_length'
include { GET_UNMAPPED              } from '../modules/local/get_unmapped'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { ANCHORS           } from '../subworkflows/local/anchors'
include { ANALYZE_BULK      } from '../subworkflows/local/analyze_bulk'
include { ANALYZE_10X       } from '../subworkflows/local/analyze_10X'
include { ANNOTATE          } from '../subworkflows/local/annotate'
include { SKIP_GET_ANCHORS  } from '../subworkflows/local/skip_get_anchors'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { TRIMGALORE        } from '../modules/nf-core/modules/trimgalore/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/



workflow STRINGSTATS {

     // Parse samplesheet
    Channel
        .fromPath(params.input)
        .splitCsv(
            header: false
        )
        .map { row ->
            tuple(
                file(row[0]).simpleName,
                file(row[0]),
                row[1]
            )
        }
        .set{ ch_fastqs }

    // Define lookahead parameter
    if (params.use_read_length) {
        /*
        // Get read length of dataset
        */
        GET_READ_LENGTH(
            ch_fastqs.first()
        )
        read_length = GET_READ_LENGTH.out.read_length.toInteger()
        lookahead = ((read_length - 2 * params.kmer_size) / 2).toInteger()

    } else {
        lookahead = params.lookahead
    }

    /*
    // Trim fastqs
    */
    TRIMGALORE(
        ch_fastqs
    )

    ch_fastqs = TRIMGALORE.out.fastq

    if (params.anchors_file) {

        SKIP_GET_ANCHORS(
            ch_fastqs,
            lookahead
        )

    } else {

        /*
        // Check if we are only using unmapped reads
        */
        if (params.unmapped) {
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
        // Process anchors and targets
        */
        ANCHORS(
            ch_fastqs,
            lookahead
        )

        abundant_anchors = ANCHORS.out.abundant_anchors

        if (params.run_type == "bulk") {
            // Perform analysis on anchors and targets
            ANALYZE_BULK(
                ch_fastqs,
                lookahead,
                abundant_anchors
            )

            anchor_target_counts = ANALYZE_BULK.out.anchor_target_counts
            anchor_scores = ANALYZE_BULK.out.anchor_scores

        } else if (params.run_type == "10X") {
            // Perform analysis on anchors and targets
            ANALYZE_10X(
                ch_fastqs,
                lookahead,
                abundant_anchors
            )

            anchor_target_counts = ANALYZE_10X.out.anchor_target_counts
            anchor_scores = ANALYZE_10X.out.anchor_scores
        }

        // Annotate anchors and targets
        ANNOTATE(
            anchor_target_counts,
            anchor_scores
        )

    }
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
