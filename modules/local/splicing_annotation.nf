
process SPLICING_ANNOTATION {

    tag "${index_name}"
    label 'process_low'
    container biocontainers/bedtools

    input:
    path anchor_fasta
    path target_fasta
    val index_bowtie // separate index for bowtie and bowtie2 should be made

    output:
    path "anchor_hits*tsv"  , emit: anchor_hits
    path "target_hits*tsv"  , emit: target_hits

    script:
    anchor_hits             = "anchor_hits_${index_name}.tsv"
    target_hits             = "target_hits_${index_name}.tsv"
    """
    
    # map the anchors and targets against the reference genome using bowtie

    echo -e "anchor\tanchor_hits_${index_name}\tanchor_hits_pos_${index_name}" >> ${anchor_hits}
    bowtie -f -n 0 -l 27 -x ${index_bowtie} -U ${anchor_fasta} --quiet 
        >> ${anchor_hits}

    echo -e "target\ttarget_hits_${index_name}\ttarget_hits_pos_${index_name}" >> ${target_hits}
    bowtie -f -n 0 -l 27 -x ${index_bowtie} -U ${target_fasta} --quiet 
        >> ${target_hits}

    
    # transform the samfile to bed
    sam2bed anchor_hits.sam > anchor_hits.bed
    sam2bed target_hits.sam > target_hits.bed
    cat anchor_hits.bed target_hits.bed > in.bed

    # use bedtools to get the closest exon
    bedtools sort -i in.bed > in.sorted.bed
    bedtools closest -d -t first -a in.sorted.bed -b /scratch/groups/horence/marieke/indexes/gtf_GRCh38/Homo_sapiens.GRCh38.103.exons.sorted.bed > out.bed

    # make a separated bed file for anchors/targets that are present in the exons
    python get_exon_seq.py

    # use bedtools to get the closest ss if anchor/target is in exon
    bedtools sort -i exon_in.bed > exon_in.sorted.bed
    bedtools closest -d -t first -a exon_in.sorted.bed -b /scratch/groups/horence/marieke/indexes/gtf_GRCh38/Homo_sapiens.GRCh38.103.ss.sorted.bed > exon_out.bed

    # report back
    python report_splicing.py
    """
}
