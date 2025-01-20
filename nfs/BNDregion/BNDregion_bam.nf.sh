#!/usr/bin/env nextflow

workflow run_BNDregion_bed {

    take:
        workflow_name

    main:
        println "workflow_name = ${workflow_name}"
        println "BND width = ${params.band_width}"
        println "rbase_version = ${params.rbase_version}"
        println "profile_name = ${params.profile_name}"
        
        filtered_samples_ch = process_manifest_and_extract_samples( params.manifest_file )
        //filtered_samples_ch.view()

        local_manifest = ch_local_manifest( params.manifest_file, filtered_samples_ch ) 
        //local_manifest.view()
        extract_BNDregion(local_manifest)
        //extract_BNDregion.out.view()
/*
        if (params.repository_name == 'gcp') {
            ch_filtered_input = filter_empty_bed_files(extract_BNDregion.out.map{ output_barcode, tm_object_id, nm_object_id, repository, project, study, bed_file, env, cmds -> return[output_barcode, tm_object_id, nm_object_id, repository, project, study, bed_file]})
            //ch_filtered_input.view()
            upload_to_gcphmf(ch_filtered_input)
        } else {
            ch_filtered_input = filter_empty_bed_files(extract_BNDregion.out.map{ output_barcode, tm_object_id, nm_object_id, repository, project, study, bed_file, env, cmds -> return[output_barcode, tm_object_id, nm_object_id, repository, project, study, bed_file]})
            //ch_filtered_input.view()
            upload_to_aws(ch_filtered_input)      
        }
*/
}


workflow run_BNDregion_bam {

    take:
        workflow_name
        ch_access_bam_out

    main:
        println "workflow_name = ${workflow_name}"
        println "BNDregion_save_directory = ${params.nf_storage_workDir}"
        println "gatk_version = ${params.gatk_version}"        
        aws_igenome_base    = "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK"

        if (params.fasta_ref_verison == "Homo_sapiens_assembly19") {
                fasta_ch            = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/ref_genome/37/Homo_sapiens_assembly19.fasta",checkIfExists: true))
                fasta_index_ch      = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/ref_genome/37/Homo_sapiens_assembly19.fasta.fai",checkIfExists: true))
                fasta_dict_ch       = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/ref_genome/37/Homo_sapiens_assembly19.dict",checkIfExists: true))
        } else if (params.fasta_ref_verison == "Homo_sapiens.GRCh37.GATK.illumina") {
                fasta_ch            = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/Ref-Genome/Homo_sapiens.GRCh37.GATK.illumina.fasta",checkIfExists: true))
                fasta_index_ch      = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/Ref-Genome/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai",checkIfExists: true))
                fasta_dict_ch       = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/Ref-Genome/Homo_sapiens.GRCh37.GATK.illumina.dict",checkIfExists: true))
        } else {
                if (params.profile_name == "aws") {
                        fasta_ch            = Channel.value(file("${aws_igenome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta",checkIfExists: true))
                        fasta_index_ch      = Channel.value(file("${aws_igenome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.fai",checkIfExists: true))
                        fasta_dict_ch       = Channel.value(file("${aws_igenome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict",checkIfExists: true))
                } else {
                        gcp_genome_base    = "${params.cloud_bucket}/TCGA/ngi-igenomes/igenomes/Homo_sapiens/GATK"
                        fasta_ch            = Channel.value(file("${gcp_genome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta",checkIfExists: true))
                        fasta_index_ch      = Channel.value(file("${gcp_genome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.fai",checkIfExists: true))
                        fasta_dict_ch       = Channel.value(file("${gcp_genome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict",checkIfExists: true))
                }
        }

        local_manifest = ch_local_manifest_aws( params.manifest_file )
        pairs_bam = make_pair_bam(local_manifest, ch_access_bam_out)
        pairs_bam.view()
        ch_slice_pair_to_object_out = ch_slice_pair_to_object(pairs_bam)
        ch_slice_pair_to_object_out.view()

        if (params.repository_name == 'gcphmf') {
            slice_bam_pair_hmf(ch_slice_pair_to_object_out, fasta_ch, fasta_index_ch, fasta_dict_ch)
        } else {
            slice_bam_pair(ch_slice_pair_to_object_out, fasta_ch, fasta_index_ch, fasta_dict_ch)
            //slice_bam_pair.out.view()            
        }

}

workflow run_BNDregion_bam_hmf {

    take:
        workflow_name

    main:
        println "workflow_name = ${workflow_name}"
        println "BNDregion_save_directory = ${params.nf_storage_workDir}"
        println "gatk_version = ${params.gatk_version}"        
        println "gcphmf_app_credential_file = ${params.gcphmf_app_credential_file}"

        fasta_ch            = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/Ref-Genome/Homo_sapiens.GRCh37.GATK.illumina.fasta",checkIfExists: true))
        fasta_index_ch      = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/Ref-Genome/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai",checkIfExists: true))
        fasta_dict_ch       = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/Ref-Genome/Homo_sapiens.GRCh37.GATK.illumina.dict",checkIfExists: true))

        // this is for cloud mount version
        local_manifest = ch_local_manifest_hmf( params.manifest_file ) 
        
        //local_manifest = ch_local_manifest_hmf_ver_down( params.manifest_file )
        ch_slice_pair_to_object_out = ch_slice_pair_to_object_hmf(local_manifest)
        //ch_slice_pair_to_object_out.view()

        if (params.repository_name == 'gcphmf') {
            slice_bam_pair_hmf(ch_slice_pair_to_object_out, fasta_ch, fasta_index_ch, fasta_dict_ch)
        } else {
            slice_bam_pair_hmf(ch_slice_pair_to_object_out, fasta_ch, fasta_index_ch, fasta_dict_ch)
            //slice_bam_pair.out.view()            
        }

}

/*
workflow run_convert_to_cram {

    take:
        workflow_name
        ch_access_bam_out

    main:
        println "workflow_name = ${workflow_name}"
        println "BNDregion_save_directory = ${params.nf_storage_workDir}"
        println "gatk_version = ${params.gatk_version}"        
        aws_igenome_base    = "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK"

        if (params.fasta_ref_verison == "Homo_sapiens_assembly19") {
                fasta_ch            = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/ref_genome/37/Homo_sapiens_assembly19.fasta",checkIfExists: true))
                fasta_index_ch      = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/ref_genome/37/Homo_sapiens_assembly19.fasta.fai",checkIfExists: true))
                fasta_dict_ch       = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/ref_genome/37/Homo_sapiens_assembly19.dict",checkIfExists: true))
        } else {

                if (params.profile_name == "aws") {
                        fasta_ch            = Channel.value(file("${aws_igenome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta",checkIfExists: true))
                        fasta_index_ch      = Channel.value(file("${aws_igenome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.fai",checkIfExists: true))
                        fasta_dict_ch       = Channel.value(file("${aws_igenome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict",checkIfExists: true))
                } else {
                        gcp_genome_base    = "${params.cloud_bucket}/TCGA/ngi-igenomes/igenomes/Homo_sapiens/GATK"
                        fasta_ch            = Channel.value(file("${gcp_genome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta",checkIfExists: true))
                        fasta_index_ch      = Channel.value(file("${gcp_genome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.fai",checkIfExists: true))
                        fasta_dict_ch       = Channel.value(file("${gcp_genome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict",checkIfExists: true))
                }
        }


        local_manifest = ch_local_manifest_aws( params.manifest_file )
        pairs_bam = make_pair_bam(local_manifest, ch_access_bam_out)
        pairs_bam.view()

        bam_to_cram(pairs_bam, fasta_ch, fasta_index_ch, fasta_dict_ch)
        
        //bam_to_cram_out = bam_to_cram.out.map{output_barcode, study, repository, project, tm_object_id, nm_object_id,tm_cram_name,nm_cram_name, tm_bam, nm_bam, env, cmds -> return [output_barcode, study, repository, project, tm_object_id, nm_object_id,tm_cram_name,nm_cram_name, tm_bam, nm_bam]}

}
*/

workflow run_convert_to_cram {

    take:
        workflow_name

    main:
        println "workflow_name = ${workflow_name}"
        println "BNDregion_save_directory = ${params.nf_storage_workDir}"
        println "gatk_version = ${params.gatk_version}"        
        aws_igenome_base    = "s3://ngi-igenomes/igenomes/Homo_sapiens/GATK"

        if (params.fasta_ref_verison == "Homo_sapiens_assembly19") {
                fasta_ch            = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/ref_genome/37/Homo_sapiens_assembly19.fasta",checkIfExists: true))
                fasta_index_ch      = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/ref_genome/37/Homo_sapiens_assembly19.fasta.fai",checkIfExists: true))
                fasta_dict_ch       = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/ref_genome/37/Homo_sapiens_assembly19.dict",checkIfExists: true))
        } else {

                if (params.profile_name == "aws") { 
                        fasta_ch            = Channel.value(file("${aws_igenome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta",checkIfExists: true))
                        fasta_index_ch      = Channel.value(file("${aws_igenome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.fai",checkIfExists: true))
                        fasta_dict_ch       = Channel.value(file("${aws_igenome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict",checkIfExists: true))
                } else {
                        gcp_genome_base    = "${params.cloud_bucket}/TCGA/ngi-igenomes/igenomes/Homo_sapiens/GATK"
                        fasta_ch            = Channel.value(file("${gcp_genome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta",checkIfExists: true))
                        fasta_index_ch      = Channel.value(file("${gcp_genome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.fai",checkIfExists: true))
                        fasta_dict_ch       = Channel.value(file("${gcp_genome_base}/${params.genome}/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict",checkIfExists: true))
                }
        }


        cram_manifest = ch_convert_cram_manifest_gcp( params.manifest_file )
        cram_manifest.view()
            
        cram_manifest_out = cram_manifest  
            .map{ output_barcode, study, repository,  project, tm_object_id, nm_object_id, tm_file_path,  nm_file_path, coverage_cutoff, GB_size ->
                file_basename_tm = tm_file_path.split('/')[-1]
                file_basename_nm = nm_file_path.split('/')[-1]
                tm_bam = returnCloudBamPath_TCGA(file_basename_tm, tm_object_id)
                nm_bam = returnCloudBamPath_TCGA(file_basename_nm, nm_object_id)

                return [output_barcode, study, repository, project, tm_object_id, nm_object_id, tm_bam, nm_bam] 
            }
        //cram_manifest_out.view()
        cram_pair_to_object_out=cram_pair_to_object(cram_manifest_out)
        //cram_pair_to_object_out.view()
        bam_to_cram(cram_pair_to_object_out, fasta_ch, fasta_index_ch, fasta_dict_ch)
}


workflow run_no_need_rm_bam{
    
    take:
        workflow_name

    main:
        println "workflow_name = ${workflow_name}"

        rm_manifest = ch_rm_bam_manifest_gcp(params.manifest_file)
        //rm_manifest.view()

        rm_manifest_out = rm_manifest  
            .map{ output_barcode, study, repository,  project, tm_object_id, nm_object_id, tm_file_path,  nm_file_path, workflow, need_bam, action ->
                file_basename_tm = tm_file_path.split('/')[-1]
                file_basename_nm = nm_file_path.split('/')[-1]
                tm_bam = returnCloudBamPath_TCGA(file_basename_tm, tm_object_id)
                nm_bam = returnCloudBamPath_TCGA(file_basename_nm, nm_object_id)

                return [output_barcode, study, repository, project, tm_object_id, nm_object_id, tm_bam, nm_bam] 
            }
        rm_manifest_out.view()

        no_need_rm_bam(rm_manifest_out)
}



/* 
=================================================================
=================================================================

Define processes

=================================================================
=================================================================
*/

process extract_BNDregion {
    
    tag "${output_barcode}"

    label 'extract_BNDregion'

    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bed/${output_barcode}", pattern: "*.bed", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bed/${output_barcode}", pattern: "*{command,exitcode}*", mode: 'copy'
    
    input:
        tuple val(output_barcode), val(tm_object_id), val(nm_object_id), val(repository), val(project), val(study), path(aa_files)

    output:
        tuple val(output_barcode), val(tm_object_id), val(nm_object_id), val(repository), val(project), val(study), path("*.bed"), path("env.txt"), path("*{command,exitcode}*",hidden:true)

    script:
    """
    #!/bin/bash
    Rscript ${params.nf_home}/R/nextflow/BNDregion2Bed_v3.R \
    --aa_out_path ${aa_files} \
    --aa_barcode ${output_barcode} \
    --out_dir . \
    --refbuild hg19 \
    --band_width ${params.band_width} 

    ls -al -R . >> env.txt
    """
}


process upload_to_aws {
    
    tag "${output_barcode}"

    label 'upload_to_bucket'

    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bed/${output_barcode}", pattern: "*.bed", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bed/${output_barcode}", pattern: "*{command,exitcode}*", mode: 'copy'
 
    input:
        tuple val(output_barcode), val(tm_object_id), val(nm_object_id), val(repository), val(project), val(study), path(bed_file)

    output:
        tuple val(output_barcode), path("${output_barcode}_BNDregion.bed"), path("*{command,exitcode}*", hidden:true)

    script:
    """
    aws configure set aws_access_key_id {id}
    aws configure set aws_secret_access_key {id}

    aws s3 cp ${bed_file} ${params.cloud_bucket}/${params.profile_name}/scratch/results/${study}/${project}/BNDregion_bed/${output_barcode}/
    """
}

process upload_to_gcp {
    
    tag "${output_barcode}"

    label 'upload_to_bucket'

    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bed/${output_barcode}", pattern: "*.bed", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bed/${output_barcode}", pattern: "*{command,exitcode}*", mode: 'copy'
 
    input:
        tuple val(output_barcode), val(tm_object_id), val(nm_object_id), val(repository), val(project), val(study), path(bed_file)

    output:
        tuple val(output_barcode), path("${output_barcode}_BNDregion.bed"), path("*{command,exitcode}*", hidden:true)
    
    script:
    """
    gsutil -u {project_name} -m cp ${bed_file} ${params.cloud_bucket}/${params.profile_name}/scratch/results/${study}/${project}/BNDregion_bed/${output_barcode}/
    echo 0 >> finish_flag.txt

    """
}

process upload_to_gcphmf {
    
    tag "${output_barcode}"

    label 'upload_to_bucket'

    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bed/${output_barcode}", pattern: "*.bed", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bed/${output_barcode}", pattern: "*{command,exitcode}*", mode: 'copy'
 
    input:
        tuple val(output_barcode), val(tm_object_id), val(nm_object_id), val(repository), val(project), val(study), path(bed_file)

    output:
        tuple val(output_barcode), path("${output_barcode}_BNDregion.bed"), path("*{command,exitcode}*", hidden:true)
    
    script:
    """
    gsutil -u {project_name} -m cp ${bed_file} ${params.cloud_bucket}/${params.profile_name}/scratch/results/${study}/${project}/BNDregion_bed/${output_barcode}/
    echo 0 >> finish_flag.txt

    """
}

process slice_bam_pair {

    tag "${object_id}"

    label 'slice_bam_pair'
    
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bam/${output_barcode}", pattern: "*", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bam/${output_barcode}", pattern: "*{command,exitcode}*", mode: 'copy'
    //publishDir "${params.nf_storage_workDir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bam/${output_barcode}", pattern: "*", mode: 'copy'

    input:
        //tuple val(output_barcode), val(tm_object_id), val(nm_object_id), val(study), val(repository), val(project), path(tm_bam), path(tm_bai), path(tm_bam_header), path(nm_bam), path(nm_bai), path(nm_bam_header), path(bed_file)
        tuple val(output_barcode), val(study), val(repository), val(project), val(object_id), path(bam), path(bai), path(header), path(bed_file), val(tm_or_nm)
        file(genome_fasta)
        file(genome_index)
        file(genome_dict)

    output:
        //tuple val(output_barcode), val(study), val(repository), val(project), path(tm_bam_header), path(nm_bam_header), path("${tm_object_id}_BNDregion.bam"), path("${tm_object_id}_BNDregion.bai"), path("${nm_object_id}_BNDregion.bam"), path("${nm_object_id}_BNDregion.bai"), path("env.txt"), path("*{command,exitcode}*", hidden:true)
        tuple val(output_barcode), val(study), val(repository), val(project), path(header), path("${object_id}_BNDregion.bam"), path("${object_id}_BNDregion.bai"), path("env.txt"), path("*{command,exitcode}*", hidden:true)

    script:
    """
    #!/bin/bash
    gatk -version >> env.txt

    java -jar /opt/gatk/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar PrintReads \
        -R ${genome_fasta} \
        -I ${bam} \
        -L ${bed_file} \
        -O ${object_id}_BNDregion.bam \
        --disable-tool-default-read-filters


    echo 0 >> finish_flag.txt
    ls -al -R . >> env.txt
    """
}


process slice_bam_pair_hmf_downver {
    tag "${object_id}"

    label 'slice_bam_pair'
    
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion/${output_barcode}/${object_id}", pattern: "*{txt,bam}", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion/${output_barcode}/${object_id}", pattern: "*{command,exitcode}*", mode: 'copy'
    //publishDir "${params.nf_storage_workDir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion/${output_barcode}/${object_id}", pattern: "*", mode: 'copy'

    input:
        tuple val(output_barcode), val(study), val(repository), val(project), val(object_id), path(cram), path(crai), path(bed_file), val(tm_or_nm)
        file(genome_fasta)
        file(genome_index)
        file(genome_dict)

    output:
        tuple val(output_barcode), val(study), val(repository), val(project), path("${object_id}_cram_file_info.txt"), path("${object_id}_BNDregion.bam"), path("env.txt"), path("*{command,exitcode}*", hidden:true)

    script:
    """
    #!/bin/bash
    export GOOGLE_APPLICATION_CREDENTIALS=${params.gcphmf_app_credential_file}

    gsutil -u {project_name} ls -L ${params.cloud_bucket}/hmf_files/${params.project_name}/${params.profile_name}/${params.repository_name}/${params.genome}/${params.manifest_batch}/${object_id}/${cram} >> ${object_id}_cram_file_info.txt

    java -Xmx8G -cp /opt/hmftools/bamtools/bam-tools_v1.2.1.jar com.hartwig.hmftools.bamtools.slice.RegionSlicer \
    -output_dir ./ \
    -output_prefix ${object_id}_BNDregion \
    -threads 2 \
    -ref_genome ${genome_fasta} \
    -bam_file ${cram} \
    -regions_file ${bed_file} \
    -write_bam \
    -log_level INFO

    echo 0 >> finish_flag.txt
    ls -al -R . >> env.txt
    """
}

process slice_bam_pair_hmf {

    tag "${object_id}"

    label 'slice_bam_pair_bamtools_v13rc1'
    
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion/${output_barcode}/${object_id}", pattern: "*{txt,bam}", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion/${output_barcode}/${object_id}", pattern: "*{command,exitcode}*", mode: 'copy'

    input:
        tuple val(output_barcode), val(study), val(repository), val(project), val(object_id), path(cram), val(file_path), path(bed_file), val(tm_or_nm)
        file(genome_fasta)
        file(genome_index)
        file(genome_dict)

    output:
        tuple val(output_barcode), val(study), val(repository), val(project), path("*_cram_file_info.txt"), path("*.header.txt"), path("${object_id}_BNDregion.bam"), path("env.txt"), path("*{command,exitcode}*", hidden:true)

    script:
    """
    #!/bin/bash
    export GOOGLE_APPLICATION_CREDENTIALS=${params.gcphmf_app_credential_file}

    gsutil -u {project_name} ls -L ${file_path} >> ${cram}_cram_file_info.txt
    samtools view -H *.cram >> ${cram}.header.txt
    samtools index *.cram 

    java -Xmx50G -Xlog:gc -cp /opt/hmftools/bamtools/bam-tools_v1.3-rc.1.jar com.hartwig.hmftools.bamtools.slice.RegionSlicer \
        -output_dir ./ \
        -output_prefix ${object_id}_BNDregion \
        -threads 7 \
        -ref_genome ${genome_fasta} \
        -bam_file ${cram} \
        -regions_file ${bed_file} \
        -write_bam \
        -log_level INFO


    echo 0 >> finish_flag.txt
    ls -al -R . >> env.txt
    """
}

process slice_bam_pair_hmf_primtread {

    tag "${object_id}"

    label 'slice_bam_pair_gatk'
    
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion/${output_barcode}/${object_id}", pattern: "*{txt,cram,crai}", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion/${output_barcode}/${object_id}", pattern: "*{command,exitcode}*", mode: 'copy'
    //publishDir "${params.nf_storage_workDir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_cram/${output_barcode}", pattern: "*", mode: 'copy'

    input:
        tuple val(output_barcode), val(study), val(repository), val(project), val(object_id), path(cram), val(file_path), path(bed_file), val(tm_or_nm)
        file(genome_fasta)
        file(genome_index)
        file(genome_dict)

    output:
        tuple val(output_barcode), val(study), val(repository), val(project), path("*_cram_file_info.txt"), path("*.header.txt"), path("${object_id}_BNDregion.cram"), path("env.txt"), path("*{command,exitcode}*", hidden:true)

    script:
    """
    #!/bin/bash
    gatk -version >> env.txt

    export GOOGLE_APPLICATION_CREDENTIALS=${params.gcphmf_app_credential_file}

    gsutil -u {project_name} ls -L ${file_path} >> ${cram}_cram_file_info.txt
    samtools view -H *.cram >> ${cram}.header.txt
    samtools index *.cram 

    java -Xmx40G -jar /opt/gatk/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar PrintReads \
        -R ${genome_fasta} \
        -I ${cram} \
        -L ${bed_file} \
        -O ${object_id}_BNDregion.cram \
        --disable-tool-default-read-filters


    echo 0 >> finish_flag.txt
    ls -al -R . >> env.txt
    """
}

process copy_sliced_bam {

    tag "${output_barcode}"

    label 'access_bam_gcptcga'
    
    //publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bam_copy/${output_barcode}", pattern: "*.txt", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bam_copy/${output_barcode}", pattern: "*{command,exitcode}*", mode: 'copy'
    //publishDir "${params.nf_local_workDir}/BNDregion_bam/${output_barcode}", pattern: "*.{bam,bai}", mode: 'copy'

    input:
        tuple val(output_barcode), val(study), val(repository), val(project), path(T_bam), path(T_bai), path(N_bam), path(N_bai)

    output:
        tuple val(output_barcode), path(T_bam), path(N_bam), path("finish_flag.txt"), path("*{command,exitcode}*", hidden:true)

    script:
    """
    mkdir -p ${params.nf_storage_workDir}/BNDregion_bam/${output_barcode} && \

    echo "Copying from gs: ${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bam/${output_barcode}"
    echo "To local directory: ${params.nf_storage_workDir}/BNDregion_bam/${output_barcode}"
  
    gsutil -u "${params.gcp_project_id}" -m cp -r ${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/BNDregion_bam/${output_barcode}/* ${params.nf_storage_workDir}/BNDregion_bam/${output_barcode}/
    echo 0 >> finish_flag.txt

    """
}

process bam_to_cram {

    tag "${object_id}"

    label 'bam_to_cram'
    
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/bam_to_cram/${output_barcode}/${object_id}", pattern: "*.cram", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/bam_to_cram/${output_barcode}/${object_id}", pattern: "*{command,exitcode}*", mode: 'copy'
    //publishDir "${params.nf_storage_workDir}/results/${params.profile_name}/${study}/${repository}/${project}/bam_to_cram/${output_barcode}", pattern: "*", mode: 'copy'

    input:  
        tuple val(output_barcode), val(study), val(repository), val(project), val(object_id), path(bam), val(tm_or_nm)
        file(genome_fasta)
        file(genome_index)
        file(genome_dict)

    output:
        tuple val(output_barcode), val(study), val(repository), val(project), val(object_id), path("${bam.getBaseName()}.cram"), path("env.txt"), path("*{command,exitcode}*", hidden:true)

    script:
        def bamBaseName = bam.getBaseName()
        //println "Extracted BAM Base Name: ${bamBaseName}" 

        def cram_name = bamBaseName + '.cram'
        println "Generated CRAM Name: ${cram_name}" 

    """
    #!/bin/bash
    set -euo pipefail
    samtools --version >> env.txt
    echo "Converting $bam to $cram_name" >> env.txt

    samtools view -C --threads 6 -T ${genome_fasta} -o ${cram_name} ${bam}
    
    echo 0 >> finish_flag.txt
    ls -al -R . >> env.txt
    """
}


process transfer_and_cleanup {

    tag "${output_barcode}"

    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/transfer_and_cleanup/${output_barcode}", pattern: "*.txt", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/transfer_and_cleanup/${output_barcode}", pattern: "*{command,exitcode}*", mode: 'copy'

    input:
        tuple val(output_barcode), val(study), val(repository), val(project), val(tm_object_id), val(nm_object_id), val(tm_bam), val(nm_bam), val(tm_cram), val(nm_cram)

    output: 
        path("env.txt")
        path("*{command,exitcode}*", hidden:true)
    
    label 'transfer_file_to_bucket'

    script:
        tm_cram_name = tm_cram.split('/')[-1]
        nm_cram_name = nm_cram.split('/')[-1]
    """
    #!/bin/bash

    mkdir -p /mnt/storage/gcp_tcga/results/${study}/${repository}/${project}/cram/${tm_object_id} && \
    mkdir -p /mnt/storage/gcp_tcga/results/${study}/${repository}/${project}/cram/${nm_object_id} && \

    gsutil -u {project_name} cp ${tm_cram} /mnt/storage/gcp_tcga/results/${study}/${repository}/${project}/cram/${tm_object_id}/ && \
    gsutil -u {project_name} cp ${nm_cram} /mnt/storage/gcp_tcga/results/${study}/${repository}/${project}/cram/${nm_object_id}/ && \

    if [[ -f /mnt/storage/gcp_tcga/results/${study}/${repository}/${project}/cram/${tm_object_id}/${tm_cram_name} && -f /mnt/storage/gcp_tcga/results/${study}/${repository}/${project}/cram/${nm_object_id}/${nm_cram_name} ]]; then
        echo "CRAM files have been successfully saved to the local directory"

        gsutil -u {project_name} -m rm ${tm_bam} && \
        gsutil -u {project_name} -m rm ${nm_bam} && \

        gsutil -u {project_name} -m rm ${tm_cram} && \
        gsutil -u {project_name} -m rm ${nm_cram}

        echo "BAM and CRAM files have been successfully deleted from the remote bucket"
    else
        echo "Failed to copy CRAM files. Not deleting local files" >&2
        exit 1
    fi

    echo "Transfer and cleanup completed." >> env.txt
    ls -al -R . >> env.txt

    """
}



process no_need_rm_bam {

    tag "${output_barcode}"

    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/no_need_rm_bam/${output_barcode}", pattern: "*.txt", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/no_need_rm_bam/${output_barcode}", pattern: "*{command,exitcode}*", mode: 'copy'

    input:
        tuple val(output_barcode), val(study), val(repository), val(project), val(tm_object_id), val(nm_object_id), val(tm_bam), val(nm_bam)

    output: 
        path("env.txt")
        path("*{command,exitcode}*", hidden:true)
    
    label 'transfer_file_to_bucket'

    script:
    """
    #!/bin/bash

    gsutil -u {project_name} -m rm ${tm_bam} && \
    gsutil -u {project_name} -m rm ${nm_bam} && \


    echo "no_need_rm_bam completed." >> env.txt
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

def extract_key(manifest_path) {
    def parts = manifest_path.split('/')
    def filename = parts[-1]
    def key = filename.split('_manifest')[0]
    return key
}
def get_aaSuite_output_file(key) {
    if (key == 'tcga_hmftools_somatic_ss') {
        return params.tcga_aa_somatic_ss_file
    } else if (key == 'tcga_hmftools_somatic_ms') {
        return params.tcga_aa_somatic_ms_file
    } else if (key == 'pcawg_hmftools_somatic_ss') {
        return params.pcawg_aa_somatic_ss_file
    } else if (key=='hmf_hmftools_somatic_ss') {
        return params.hmf_aa_somatic_ss_file
    } else if (key=='hmf_hmftools_somatic_ms') {
        return params.hmf_aa_somatic_ms_file
    } else {
        throw new Exception("Unknown key: ${key}")
    }
}
def process_manifest_and_extract_samples(manifest_file) {
    def key = extract_key(manifest_file)
    def aaSuite_output_file = get_aaSuite_output_file(key)

    def samples_ch = Channel
        .fromPath(aaSuite_output_file)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.aa_barcode, row.cp_output_dir)
        }
        .map { output_barcode, aa_files ->
            return [output_barcode, aa_files]
        }

    return samples_ch
}


/*
def ch_local_manifest(manifest_file){
    samples_ch = Channel
            .fromPath(manifest_file)
            .splitCsv(header:true)
            .map{ row-> tuple(  row.aa_barcode,
                                row.tm_object_id,
                                row.tm_file_path,
                                row.tm_ega_file_id,
                                row.nm_object_id,
                                row.nm_file_path,
                                row.nm_ega_file_id,
                                row.repository,
                                row.patient_barcode,
                                row.project,
                                row.study,
                                row.action) }
            .filter{    aa_barcode, tm_object_id, tm_file_path, tm_ega_file_id, nm_object_id, nm_file_path, nm_ega_file_id, repository, patient_barcode, project, study, action->
                action == "run" }
            .map{   output_barcode, tm_object_id, tm_file_path, tm_ega_file_id, nm_object_id, nm_file_path, nm_ega_file_id, repository, patient_barcode, project, study, action-> 
                def aa_files = file("/mnt/NAS3/storage/EcDNA_Advanced/data/aaSuite_files/${study}/aaSuite_somatic_ss/v0.1344.2/GRCh37/minCN4.5/cnsizeMin50000/10X/calls/${output_barcode}/aaSuite_output/${output_barcode}/${output_barcode}_output")
                return[ output_barcode, tm_object_id,                                nm_object_id,                               repository,                  project, study,         aa_files]}
}
*/

def ch_local_manifest(manifest_file, filtered_samples_ch){
    samples_ch = Channel
            .fromPath(manifest_file)
            .splitCsv(header:true)
            .map{ row-> tuple(  row.aa_barcode,
                                row.tm_object_id,
                                row.tm_file_path,
                                row.tm_ega_file_id,
                                row.nm_object_id,
                                row.nm_file_path,
                                row.nm_ega_file_id,
                                row.repository,
                                row.patient_barcode,
                                row.project,
                                row.study,
                                row.action) }
            .filter{    aa_barcode, tm_object_id, tm_file_path, tm_ega_file_id, nm_object_id, nm_file_path, nm_ega_file_id, repository, patient_barcode, project, study, action->
                action == "run" }
            .map{        output_barcode, tm_object_id, tm_file_path, tm_ega_file_id, nm_object_id, nm_file_path, nm_ega_file_id, repository, patient_barcode, project, study, action-> 
                return [ output_barcode, tm_object_id,                               nm_object_id,                               repository,                  project, study        ]}
            .combine(
                filtered_samples_ch.map{ output_barcode, aa_files -> return [output_barcode, aa_files]},by:[0])
            .map{ output_barcode, tm_object_id, nm_object_id, repository, project, study, aa_files -> return [output_barcode, tm_object_id, nm_object_id, repository, project, study, aa_files]}
}


def ch_local_manifest_aws(manifest_file){
    samples_ch = Channel
            .fromPath(manifest_file)
            .splitCsv(header:true)
            .map{ row-> tuple(  row.aa_barcode,
                                row.tm_object_id,
                                row.tm_file_path,
                                row.tm_ega_file_id,
                                row.nm_object_id,
                                row.nm_file_path,
                                row.nm_ega_file_id,
                                row.repository,
                                row.patient_barcode,
                                row.project,
                                row.study,
                                row.action) }
            .filter{    aa_barcode, tm_object_id, tm_file_path, tm_ega_file_id, nm_object_id, nm_file_path, nm_ega_file_id, repository, patient_barcode, project, study, action->
                action == "run" }
            .map{    output_barcode, tm_object_id, tm_file_path, tm_ega_file_id, nm_object_id, nm_file_path, nm_ega_file_id, repository, patient_barcode, project, study, action-> 
            return [ output_barcode, tm_object_id,                               nm_object_id,                               repository,                  project, study]}
}

def make_pair_bam(local_manifest, ch_access_bam_out){
    local_manifest
        .map{        output_barcode, tm_object_id, nm_object_id, repository,  project, study -> 
            return [ output_barcode, tm_object_id, nm_object_id                             ]
        }
        .combine(
            ch_access_bam_out.map{  output_barcode,                         study, repository, project, tm_bam, tm_bai, tm_bam_header, nm_bam, nm_bai, nm_bam_header->
                            return [output_barcode,                         study, repository, project, tm_bam, tm_bai, tm_bam_header, nm_bam, nm_bai, nm_bam_header]
            }, by:[0] 
        )
        .map{        output_barcode, tm_object_id, nm_object_id, study, repository, project, tm_bam, tm_bai, tm_bam_header, nm_bam, nm_bai, nm_bam_header-> 
            //def bed_file = file("${params.cloud_bucket}/local/scratch/results/${study}/${project}/BNDregion_bed/${output_barcode}/${output_barcode}_BNDregion.bed", checkIfExists:true)
            def bed_file = file("${params.cloud_bucket}/local/scratch/results/${study}/${project}/BNDregion_bed/${output_barcode}/${output_barcode}_BNDregion.bed")
            return bed_file.exists() ? [output_barcode, tm_object_id, nm_object_id, study, repository, project, tm_bam, tm_bai, tm_bam_header, nm_bam, nm_bai, nm_bam_header, bed_file] : null
        }
        .filter { it != null }
}


def ch_slice_pair_to_object(ch_pair_bams) {
    def ch_tumor_bams = ch_pair_bams
        .map { output_barcode, tm_object_id, nm_object_id, study, repository, project, tm_bam, tm_bai, tm_bam_header, nm_bam, nm_bai, nm_bam_header, bed_file ->
            def tm_or_nm = 'tumor'
            return [ output_barcode, study, repository, project, tm_object_id, tm_bam, tm_bai, tm_bam_header, bed_file, tm_or_nm ]
        }
        .unique()

    def ch_normal_bams = ch_pair_bams
        .map { output_barcode, tm_object_id, nm_object_id, study, repository, project, tm_bam, tm_bai, tm_bam_header, nm_bam, nm_bai, nm_bam_header, bed_file ->
            def tm_or_nm = 'normal'
            return [ output_barcode, study, repository, project, nm_object_id, nm_bam, nm_bai, nm_bam_header, bed_file, tm_or_nm  ]
        }
        .unique()

    def merged_channel = ch_tumor_bams.concat(ch_normal_bams)

    return merged_channel
        .map { output_barcode, study, repository, project, object_id, bam, bai, header, bed_file, tm_or_nm ->
            return [
                output_barcode,
                study,
                repository,
                project,
                object_id,
                bam,
                bai,
                header,
                bed_file,
                tm_or_nm
            ]
        }
        .unique()
}


def filter_empty_bed_files(input_channel){
    return input_channel.filter { tuple ->
        def (output_barcode, tm_object_id, nm_object_id, repository, project, study, bed_file) = tuple
        return bed_file.size() > 0 }
}


def ch_local_manifest_gcp(manifest_file){
    samples_ch = Channel
            .fromPath(manifest_file)
            .splitCsv(header:true)
            .map{ row-> tuple(  row.aa_barcode,
                                row.tm_object_id,
                                row.tm_file_path,
                                row.tm_ega_file_id,
                                row.nm_object_id,
                                row.nm_file_path,
                                row.nm_ega_file_id,
                                row.repository,
                                row.patient_barcode,
                                row.project,
                                row.study,
                                row.action) }
            .filter{    aa_barcode, tm_object_id, tm_file_path, tm_ega_file_id, nm_object_id, nm_file_path, nm_ega_file_id, repository, patient_barcode, project, study, action->
                action == "run" }
            .map{    output_barcode, tm_object_id, tm_file_path, tm_ega_file_id, nm_object_id, nm_file_path, nm_ega_file_id, repository, patient_barcode, project, study, action-> 
            return [ output_barcode, study, repository,  project, tm_object_id, nm_object_id, tm_file_path,  nm_file_path ]}
}


def returnBamPath_TCGA(tg_bam_name, tg_object_id) {
    if (params.repository_name == 'gcp') {
        def cloud_file_path = "/tcga_files/${params.project_name}/gcptcga/${params.repository_name}/${params.genome}/${params.manifest_batch}/${tg_object_id}/${tg_bam_name}"
        return cloud_file_path
    } else {
        throw new IllegalArgumentException("this function does not work with profile_name ${params.profile_name}")
    }
}


def returnCloudBamPath_TCGA(tg_bam_name, tg_object_id) {
    if (params.repository_name == 'gcp') {
        def cloud_file_path = "${params.cloud_bucket}/tcga_files/${params.project_name}/gcptcga/${params.repository_name}/${params.genome}/${params.manifest_batch}/${tg_object_id}/${tg_bam_name}"
        return cloud_file_path
    } else {
        throw new IllegalArgumentException("this function does not work with profile_name ${params.profile_name}")
    }
}

def returnCloudCramPath_TCGA(tg_bam_name, output_barcode,study,repository,project) {
    if (params.repository_name == 'gcp') {
        def cloud_file_path = "${params.cloud_bucket}/${params.project_name}/scratch/results/gcptcga/${study}/${repository}/${project}/bam_to_cram/${output_barcode}/${tg_bam_name}"
        return cloud_file_path
    } else {
        throw new IllegalArgumentException("this function does not work with repository_name ${params.repository_name}")
    }
}


def ch_rm_bam_manifest_gcp(manifest_file){
    samples_ch = Channel
            .fromPath(manifest_file)
            .splitCsv(header:true)
            .map{ row-> tuple(  row.aa_barcode,
                                row.tm_object_id,
                                row.tm_file_path,
                                row.nm_object_id,
                                row.nm_file_path,
                                row.repository,
                                row.patient_barcode,
                                row.project,
                                row.study,
                                row.workflow,
                                row.need_bam,
                                row.action) }
            .filter{    aa_barcode, tm_object_id, tm_file_path, nm_object_id, nm_file_path, repository, patient_barcode, project, study, workflow,need_bam, action->
                action == "run" }
            .map{   output_barcode, tm_object_id, tm_file_path, nm_object_id, nm_file_path, repository, patient_barcode, project, study, workflow,need_bam, action-> 
            return [ output_barcode, study, repository,  project, tm_object_id, nm_object_id, tm_file_path,  nm_file_path, workflow, need_bam, action ]}
}

def ch_convert_cram_manifest_gcp(manifest_file){
    samples_ch = Channel
            .fromPath(manifest_file)
            .splitCsv(header:true)
            .map{ row-> tuple(  row.aa_barcode,
                                row.tm_object_id,
                                row.tm_file_path,
                                row.nm_object_id,
                                row.nm_file_path,
                                row.repository,
                                row.patient_barcode,
                                row.project,
                                row.study,
                                row.coverage_cutoff,
                                row.GB_size) }
            .filter{    aa_barcode, tm_object_id, tm_file_path, nm_object_id, nm_file_path, repository, patient_barcode, project, study,coverage_cutoff, GB_size->
                coverage_cutoff == ">=10X" }
            .map{   output_barcode, tm_object_id, tm_file_path, nm_object_id, nm_file_path, repository, patient_barcode, project, study,coverage_cutoff, GB_size-> 
            return [ output_barcode, study, repository,  project, tm_object_id, nm_object_id, tm_file_path,  nm_file_path, coverage_cutoff, GB_size  ]}
}



def ch_local_manifest_hmf(manifest_file){
    samples_ch = Channel
            .fromPath(manifest_file)
            .splitCsv(header:true)
            .map{ row-> tuple(  row.aa_barcode,
                                row.tm_object_id,
                                row.tm_file_path,
                                row.tm_crai,
                                row.nm_object_id,
                                row.nm_file_path,
                                row.nm_crai,
                                row.repository,
                                row.patient_barcode,
                                row.project,
                                row.study,
                                row.action) }
            .filter{     aa_barcode, tm_object_id, tm_file_path, tm_crai, nm_object_id, nm_file_path, nm_crai, repository, patient_barcode, project, study, action->
                action == "run" }
            .map{    output_barcode, tm_object_id, tm_file_path, tm_crai_path, nm_object_id, nm_file_path, nm_crai_path, repository, patient_barcode, project, study, action->
                def bed_file = file("${params.cloud_bucket}/local/scratch/results/${study}/${project}/BNDregion_bed/${output_barcode}/${output_barcode}_BNDregion.bed")
                def tm_cram = file("${tm_file_path}",checkIfExists: true)
                def nm_cram = file("${nm_file_path}",checkIfExists: true)
                def tm_crai = file("${tm_crai_path}",checkIfExists: true)
                def nm_crai = file("${nm_crai_path}",checkIfExists: true)
 
            return [ output_barcode, tm_object_id, nm_object_id, study, repository, project,  tm_cram, tm_file_path, tm_crai, nm_cram, nm_file_path, nm_crai, bed_file]}
}

def ch_local_manifest_hmf_ver_down(manifest_file){
    samples_ch = Channel
            .fromPath(manifest_file)
            .splitCsv(header:true)
            .map{ row-> tuple(  row.aa_barcode,
                                row.tm_object_id,
                                row.tm_cram,
                                row.tm_file_path,
                                row.tm_crai,
                                row.nm_object_id,
                                row.nm_cram,
                                row.nm_file_path,
                                row.nm_crai,
                                row.repository,
                                row.patient_barcode,
                                row.project,
                                row.study,
                                row.action) }
            .filter{     aa_barcode, tm_object_id, tm_cram, tm_file_path, tm_crai, nm_object_id, nm_cram, nm_file_path, nm_crai, repository, patient_barcode, project, study, action->
                action == "run" }
            .map{    output_barcode, tm_object_id, tm_cram, tm_file_path, tm_crai, nm_object_id, nm_cram, nm_file_path, nm_crai, repository, patient_barcode, project, study, action->
                def bed_file = file("${params.cloud_bucket}/local/scratch/results/${study}/${project}/BNDregion_bed/${output_barcode}/${output_barcode}_BNDregion.bed")
                def f_tm_cram =  file("${params.cloud_bucket}/hmf_files/${params.project_name}/${params.profile_name}/${params.repository_name}/${params.genome}/${params.manifest_batch}/${tm_object_id}/${tm_cram}",checkIfExists: true)
                def f_nm_cram =  file("${params.cloud_bucket}/hmf_files/${params.project_name}/${params.profile_name}/${params.repository_name}/${params.genome}/${params.manifest_batch}/${nm_object_id}/${nm_cram}",checkIfExists: true)
                def f_tm_crai =  file("${params.cloud_bucket}/hmf_files/${params.project_name}/${params.profile_name}/${params.repository_name}/${params.genome}/${params.manifest_batch}/${tm_object_id}/${tm_cram}.crai",checkIfExists: true)
                def f_nm_crai =  file("${params.cloud_bucket}/hmf_files/${params.project_name}/${params.profile_name}/${params.repository_name}/${params.genome}/${params.manifest_batch}/${nm_object_id}/${nm_cram}.crai",checkIfExists: true)
 
            return [ output_barcode, tm_object_id, nm_object_id, study, repository, project,  f_tm_cram, f_tm_crai, f_nm_cram, f_nm_crai, bed_file]}
}


def ch_slice_pair_to_object_hmf(ch_pair_bams) {
    def ch_tumor_bams = ch_pair_bams
        .map { output_barcode, tm_object_id, nm_object_id, study, repository, project,  tm_cram, tm_file_path, tm_crai, nm_cram, nm_file_path, nm_crai, bed_file ->
            def tm_or_nm = 'tumor'
            return [ output_barcode, study, repository, project, tm_object_id, tm_cram, tm_file_path, tm_crai, bed_file, tm_or_nm ]
        }
        .unique()

    def ch_normal_bams = ch_pair_bams
        .map { output_barcode, tm_object_id, nm_object_id, study, repository, project,  tm_cram, tm_file_path, tm_crai, nm_cram, nm_file_path, nm_crai, bed_file ->
            def tm_or_nm = 'normal'
            return [ output_barcode, study, repository, project, nm_object_id, nm_cram, nm_file_path, nm_crai, bed_file, tm_or_nm  ]
        }
        .unique()

    def merged_channel = ch_tumor_bams.concat(ch_normal_bams)

    return merged_channel
        .map { output_barcode, study, repository, project, object_id, cram, file_path, crai, bed_file, tm_or_nm ->
            return [
                output_barcode,
                study,
                repository,
                project,
                object_id,
                cram,
                file_path,
                bed_file,
                tm_or_nm
            ]
        }
        .unique()
}
