/* 
=================================================================
=================================================================

workflow
Date: 2025-02-28
Create: Jiwon

=================================================================
=================================================================
*/

workflow aaSuite_feature_similarity {
    take:
        dna_legacy_file_sheet
    main:

        println "workflow_name              = ${params.workflow_name}"
        println "dna_legacy_file_sheet      = ${params.dna_legacy_file_sheet}"
        println "dna_pairs_manifest_sheet   = ${params.dna_pairs_manifest_sheet}"
        println "aasuite_ver                = ${params.aasuite_ver}"
        println "aaclassifier_ver           = ${params.aaclassifier_ver}"
    
        println "genome                     = ${params.genome}"
        println "aa_gain                    = ${params.aa_gain}"
        println "aa_cnsize_min              = ${params.aa_cnsize_min}"
        println "aa_downsample              = ${params.aa_downsample}"

        // data_repo_12152022_AmpliconSuite is the same bundle as used in the Pan ecDNA
        bundle_ch = Channel.value(file("${params.nf_home}/references/${params.genome}/AmpliconArchitect/data_repo_12152022_AmpliconSuite"))

        sample_bams = make_samples( params.dna_legacy_file_sheet)
        // OUTPUT: aliquot, patient_barcode, project, tm_object_id, nm_object_id, AA_out_dir
       	//sample_bams.view()


        create_features_to_graph_out = create_features_to_graph(sample_bams)
        // OUTPUT: aliquot, patient_barcode, project, tm_object_id, nm_object_id, AA_out_dir, features_to_graph, env, cmd 
        //create_features_to_graph_out.view()

        merge_by_patient = create_features_to_graph_out
            .map { aliquot, patient_barcode, project, tm_object_id, nm_object_id, AA_out_dir, features_to_graph, env, cmd -> 
                println "Aliquot: ${aliquot}, Patient: ${patient_barcode}, Features File: ${features_to_graph}"  
                return [patient_barcode, features_to_graph]
            }
            .groupTuple(by: 0)

        //merge_by_patient.view()

        concat_features_to_graph_out = concat_features_to_graph(merge_by_patient)
        //concat_features_to_graph_out.view()
        similarity_input = concat_features_to_graph_out.map { patient_barcode, features_to_graph, cnv, cmd -> return [ patient_barcode, features_to_graph ] }
        similarity_input.view()

        run_similarity(similarity_input, bundle_ch).view()


}



/* 
=================================================================
=================================================================

Define processes

=================================================================
=================================================================
*/


// -----------------------------------------------------------------
//
// aaSuite-based AA
//
// -----------------------------------------------------------------
process create_features_to_graph {
    
    tag "${aliquot}"

    publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/aaSuite_feature_similarity/${patient_barcode}/create_features_to_graph/${aliquot}", pattern: "*.txt", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/aaSuite_feature_similarity/${patient_barcode}/create_features_to_graph/${aliquot}", pattern: "*{command,exitcode}*", mode: 'copy'

    //label "aasuite"
    

    input:
        tuple val(aliquot), val(patient_barcode), val(project), val(tm_object_id), val(nm_object_id), path(AA_out_dir)

    output:
        tuple val(aliquot), val(patient_barcode), val(project), val(tm_object_id), val(nm_object_id), path(AA_out_dir), path("${aliquot}_features_to_graph.txt"), path("env.txt"), path("*{command,exitcode}*",hidden:true)
    
    script:
    """
    /home/jiwon/HL-NF/nfs/aaSuite/extract_feature_to_graph.sh ${AA_out_dir} ${aliquot}
    
    ls -al -R . >> env.txt
    """


    
}   

process concat_features_to_graph {
    
    tag "${patient_barcode}"

    publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/aaSuite_feature_similarity/${patient_barcode}/concat_features_to_graph", pattern: "*.txt", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/aaSuite_feature_similarity/${patient_barcode}/concat_features_to_graph", pattern: "*{command,exitcode}*", mode: 'copy'

    //label "aasuite"
    

    input:
        tuple val(patient_barcode), path(features_to_graphs)

    output:
        tuple val(patient_barcode), path("${patient_barcode}_features_to_graph.txt"), path("env.txt"), path("*{command,exitcode}*",hidden:true)
    
    script:
    """
    echo "Patient barcode: ${patient_barcode}"
    output_file="${patient_barcode}_features_to_graph.txt"
    rm -f "\$output_file"

    for file in ${features_to_graphs}; do
        echo "Adding \$file to \$output_file"
        cat \$file >> \$output_file
    done \

    echo "Combined file created: \$output_file" \

    ls -al -R . >> env.txt
    """
    
}  


process run_similarity {
    
    tag "${patient_barcode}"

    label 'aasuite'

    publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/aaSuite_feature_similarity/${patient_barcode}/run_similarity/min_cn${min_cn}", pattern: "*.{tsv,txt}", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aasuite_ver}/${params.genome}/minCN${params.aa_gain}/cnsizeMin${params.aa_cnsize_min}/${params.aa_downsample}X/aaSuite_feature_similarity/${patient_barcode}/run_similarity/min_cn${min_cn}", pattern: "*{command,exitcode}*", mode: 'copy'

    input:
        tuple val(patient_barcode), path(features_to_graph)
        file(aa_data_repo)

    output:
        tuple val(patient_barcode), path("${patient_barcode}_feature_similarity_scores.tsv"), path("env.txt"), path("*{command,exitcode}*",hidden:true)

    script:
    """
    #!/bin/bash

    export AA_DATA_REPO=${aa_data_repo}

    echo "AA_DATA_REPO=\$AA_DATA_REPO" >> env.txt && \
    echo "AA_SRC=\$AA_SRC" >> env.txt && \

    python /home/programs/AmpliconClassifier-main/feature_similarity.py \
    --ref ${params.genome} \
    -f ${features_to_graph} \
    -o ${patient_barcode} \
    --min_cn ${params.min_cn} \
    --add_chr_tag \
    --required_classifications any \

    ls -al -R . >> env.txt
    """
} 

/*
Singularity> /home/programs/AmpliconClassifier-main/python amplicon_classifier.py -h
usage: amplicon_classifier.py [-h] [-i INPUT] [-c CYCLES] [-g GRAPH] --ref
                              {hg19,GRCh37,hg38,GRCh38,GRCh38_viral,mm10,GRCm38}
                              [--min_flow MIN_FLOW] [--min_size MIN_SIZE]
                              [-o O] [--plotstyle {grouped,individual,noplot}]
                              [--force] [--add_chr_tag] [--report_complexity]
                              [--verbose_classification] [--no_LC_filter]
                              [--exclude_bed EXCLUDE_BED]
                              [--decomposition_strictness DECOMPOSITION_STRICTNESS]
                              [--filter_similar] [-v]
*/

/* 
=================================================================
=================================================================

Define functions

=================================================================
=================================================================
*/

// Prepare pair map
def make_samples ( dna_legacy_file_sheet  ) {
    Channel
        .fromPath(dna_legacy_file_sheet)
        .splitCsv(header:true)
        .map{ row-> tuple(  row.aa_barcode,
                            row.patient_barcode,
                            row.project,
                            row.tm_object_id,
                            row.nm_object_id,
                            row.cp_output_dir,
                            row.action ) }
        .filter{    aliquot, patient_barcode, project, tm_object_id, nm_object_id, cp_output_dir, action ->
                    action == "run" }
        .map {      aliquot, patient_barcode, project, tm_object_id, nm_object_id, AA_ou_dir, action -> 
        return [    aliquot, patient_barcode, project, tm_object_id, nm_object_id, AA_ou_dir ] }

}







