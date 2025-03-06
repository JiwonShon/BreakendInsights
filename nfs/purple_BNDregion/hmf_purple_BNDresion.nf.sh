/*
Updated: 2024-12-12
Edit: 2025-01-16 by Jiwon 
Description: This version of aa_bp_seq_seek is for analyzing AA pair mode results from the cloud pipeline, modified to separate and analyze tm and nm in single mode.
*/

/*
=================================================================
=================================================================

workflow

=================================================================
=================================================================
*/

workflow hmf_purple_BNDresion {


    main:

        println "workflow_name = ${params.step}"
        println "aa_gain = ${params.aa_gain}"
        println "aa_cnsize_min = ${params.aa_cnsize_min}"
        println "aa_downsample = ${params.aa_downsample}"
        println "cohort_name = ${params.cohort_name}"

        R_ch = Channel.value(file("${params.nf_home}/R/nextflow"))
        py_ch = Channel.value(file("${params.nf_home}/python"))
        
        println "BND width = ${params.band_width}"
        println "profile_name = ${params.profile_name}"
        
        local_manifest = make_samples( params.dna_legacy_file_sheet )
        //local_manifest.view()

        extract_BNDregion_out = extract_BNDregion(local_manifest,py_ch)
		//OUTPUT: [output_barcode, tm_object_id, nm_object_id, repository, project, project, study, bed, env, cmds]
        //extract_BNDregion_out.view()

        purple_sv_bnd_input = find_vcf(extract_BNDregion_out)
        //purple_sv_bnd_input.view()

        purple_SV_BNDregion(purple_sv_bnd_input)
        //purple_SV_BNDregion.out.view()
}

///mnt/NAS7/home/jiwon/PanEcDNA2/hmf-dr-265-update6/hmf-dr-265-update6/LINX_PURPLE/DR-265-update6/${tm_object_id}/purple/${tm_object_id}.purple.sv.vcf.gz
//process purple_SV_BNDregion 

/*
=================================================================
=================================================================

Define processes

=================================================================
=================================================================
*/

// --------------------------------------------------------------
// 
//
// Purple - vcf target slicing
//
//
// --------------------------------------------------------------

process extract_BNDregion {
    
    tag "${output_barcode}"

    label 'bp_homology_python'

    publishDir "${params.scratch_dir}/results/purple_BND/${tm_object_id}/bed/${params.band_width}", pattern: "*.{bed,txt}", mode: 'copy'
    publishDir "${params.scratch_dir}/results/purple_BND/${tm_object_id}/bed/${params.band_width}", pattern: "*{command,exitcode}*", mode: 'copy'
    
    input:
        tuple val(output_barcode), val(tm_object_id), val(nm_object_id), val(repository), val(project), val(study), path(aa_files)
        file(py_ch)

    output:
        tuple val(output_barcode), val(tm_object_id), val(nm_object_id), val(repository), val(project), val(study), path("*.bed"), path("env.txt"), path("*{command,exitcode}*",hidden:true)

    script:
    """
    #!/bin/bash
    python3 ${py_ch}/BNDregion_AAinterval.py \
    --aa_out_path ${aa_files} \
    --aa_barcode ${output_barcode} \
    --out_dir . \
    --refbuild hg19 \
    --band_width ${params.band_width} 

    ls -al -R . >> env.txt
    """
}


process purple_SV_BNDregion {
    
    tag "${output_barcode}"

    label 'bcftools'

    publishDir "${params.scratch_dir}/results/purple_BND/${tm_object_id}/purple", pattern: "*.{vcf.gz,txt,tbi}", mode: 'copy'
    publishDir "${params.scratch_dir}/results/purple_BND/${tm_object_id}/purple", pattern: "*{command,exitcode}*", mode: 'copy'
    
    input:
        tuple val(output_barcode), val(tm_object_id), val(nm_object_id), val(repository), val(project), val(study), path(bed), path(vcf_file), path(vcf_idx_file)

    output:
        tuple val(output_barcode), val(tm_object_id), val(nm_object_id), val(repository), val(project), val(study), path("${tm_object_id}.purple.BND.sv.vcf.gz"), path("${tm_object_id}.purple.BND.sv.vcf.gz.tbi"), path("bcftools_version.txt"), path("env.txt"), path("*{command,exitcode}*",hidden:true)

    script:
    """
    #!/bin/bash
    bcftools -v >> bcftools_version.txt \

    bcftools view -R ${bed} ${vcf_file} -Oz -o "${tm_object_id}.purple.BND.sv.vcf.gz" \

    bcftools index -t "${tm_object_id}.purple.BND.sv.vcf.gz" \

    ls -al -R . >> env.txt
    """
}


/* 
=================================================================
=================================================================

Define functions

=================================================================
=================================================================
*/


def make_samples ( dna_legacy_file_sheet  ) {
    ch_samples = Channel
        .fromPath(dna_legacy_file_sheet)
        .splitCsv(header:true)
        .map{ row-> tuple(  row.aa_barcode,
                            row.tm_object_id,
                            row.nm_object_id,
                            row.repository,
                            row.project,
                            row.study,
                            row.cp_output_dir,
                            row.action) }
        .map{          output_barcode, tm_object_id, nm_object_id, repository, project, study, aa_files, action ->
            def aa_file = file(aa_files).exists() ? file(aa_files) : "File not found"
            return     [output_barcode, tm_object_id, nm_object_id, repository, project, study, aa_file, action ] }
        .filter{        output_barcode, tm_object_id, nm_object_id, repository, project, study, aa_file, action ->
                action == "run" }
        .map{           output_barcode, tm_object_id, nm_object_id, repository, project, study, aa_file, action ->
                return [output_barcode, tm_object_id, nm_object_id, repository, project, study, aa_file ]}
        .unique()
}


def find_vcf(ch_access_bed){
    ch_access_bed
        .map{                           output_barcode, tm_object_id, nm_object_id, repository, project, study, bed, env, cmds -> 
            def vcf_file = file("/mnt/NAS7/home/jiwon/PanEcDNA2/hmf-dr-265-update6/hmf-dr-265-update6/LINX_PURPLE/DR-265-update6/${tm_object_id}/purple/${tm_object_id}.purple.sv.vcf.gz")
            def vcf_idx_file = file("/mnt/NAS7/home/jiwon/PanEcDNA2/hmf-dr-265-update6/hmf-dr-265-update6/LINX_PURPLE/DR-265-update6/${tm_object_id}/purple/${tm_object_id}.purple.sv.vcf.gz.tbi")
            return vcf_file.exists() ? [output_barcode, tm_object_id, nm_object_id, repository, project, study, bed, vcf_file, vcf_idx_file] : null
        }
        .filter { it != null }
}
