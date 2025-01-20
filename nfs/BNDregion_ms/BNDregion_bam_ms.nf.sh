#!/usr/bin/env nextflow

workflow run_BNDregion_bed_ms {

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
        extract_BNDregion_out = extract_BNDregion(local_manifest)
        
        emit_extract_BNDregion_out = extract_BNDregion_out.map{ output_barcode, tm_object_id, nm_object_id,patient_barcode, repository, project, study, bed_file, env, cmds -> return[output_barcode, tm_object_id, nm_object_id, patient_barcode, repository, project, study, bed_file]}
        //emit_extract_BNDregion_out.view()
    emit:
        ch_extract_BNDregion_out = emit_extract_BNDregion_out

}

workflow run_BNDregion_bed_ms_merge {
    take:
        workflow_name
        ch_extract_BNDregion_out

    main:
        println "workflow_name = ${workflow_name}"


        merge_patient_input = ch_merge_patient_input(ch_extract_BNDregion_out) 
        // [patient_barcode, repository, project, study]

        //merge_patient_input.view()

        merge_patient_input_output = merge_bed_files_by_patient(merge_patient_input)
        //merge_patient_input_output.view()
}




workflow run_merge_bed_to_cloud {
    take:
        workflow_name

    main:
        println "workflow_name = ${workflow_name}"
        
        ch_local_manifest_pair_bed_out = ch_local_manifest_pair_bed( params.manifest_file )
        ch_local_manifest_pair_bed_out.view()

        if (params.repository_name == 'gcp') {
            ch_filtered_input = filter_empty_bed_files(ch_local_manifest_pair_bed_out)
            //ch_filtered_input.view()

            upload_to_gcp(ch_filtered_input)

        } else {
            ch_filtered_input = filter_empty_bed_files(ch_local_manifest_pair_bed_out)
            //ch_filtered_input.view()
            upload_to_aws(ch_filtered_input)      
        }
}


workflow run_BNDregion_bam_ms {

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


        local_manifest = ch_local_manifest_cloud( params.manifest_file )
        pairs_bam = make_pair_bam(local_manifest, ch_access_bam_out)
        //pairs_bam.view() -> make_sliceBAM_input
        make_sliceBAM_input_out = ch_slice_pair_to_object_hmf(pairs_bam)
        slice_bam_pair(make_sliceBAM_input_out, fasta_ch, fasta_index_ch, fasta_dict_ch)
        //slice_bam_pair.out.view()

}

workflow run_BNDregion_bam_ms_hmf {

    take:
        workflow_name

    main:
        println "workflow_name = ${workflow_name}"
        //println "BNDregion_save_directory = ${params.nf_storage_workDir}"
        //println "gatk_version = ${params.gatk_version}"        
        println "gcphmf_app_credential_file = ${params.gcphmf_app_credential_file}"

        fasta_ch            = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/Ref-Genome/Homo_sapiens.GRCh37.GATK.illumina.fasta",checkIfExists: true))
        fasta_index_ch      = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/Ref-Genome/Homo_sapiens.GRCh37.GATK.illumina.fasta.fai",checkIfExists: true))
        fasta_dict_ch       = Channel.value(file("${params.cloud_bucket}/references/${params.genome}/HMFTools-Resources/Ref-Genome/Homo_sapiens.GRCh37.GATK.illumina.dict",checkIfExists: true))


        local_manifest = ch_local_manifest_hmf( params.manifest_file )

        make_sliceBAM_input_out = ch_slice_pair_to_object_hmf(local_manifest)

        slice_bam_pair_hmf(make_sliceBAM_input_out, fasta_ch, fasta_index_ch, fasta_dict_ch)
        //slice_bam_pair.out.view()

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

    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/ms_BNDregion_bed/${patient_barcode}", pattern: "*.bed", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/ms_BNDregion_bed/${patient_barcode}", pattern: "*{command,exitcode}*", mode: 'copy'
    
    input:
        tuple val(output_barcode), val(tm_object_id), val(nm_object_id), val(patient_barcode), val(repository), val(project), val(study), path(aa_files)

    output:
        tuple val(output_barcode), val(tm_object_id), val(nm_object_id), val(patient_barcode), val(repository), val(project), val(study), path("${output_barcode}_BNDregion.bed"), path("env.txt"), path("*{command,exitcode}*",hidden:true)

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


process merge_bed_files_by_patient {
    
    tag "${patient_barcode}"

    label 'extract_BNDregion'
    
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/ms_BNDregion_bed/${patient_barcode}", pattern: "*", mode: 'copy'

    input:
        tuple val(patient_barcode), val(repository), val(project), val(study)

    output:
        tuple val(patient_barcode), val(repository), val(project), val(study), path("${patient_barcode}_ms_BNDregion.bed")

    script:
    """
    #!/bin/bash

    # Find all .bed files in the specified directory
    find ${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/ms_BNDregion_bed/${patient_barcode}/ -type f -name '*.bed' > bed_files_list.txt

    if [ -s bed_files_list.txt ]; then
        echo "Found the following BED files:"
        cat bed_files_list.txt
        
        # Concatenate, sort, and deduplicate the BED files
        cat \$(cat bed_files_list.txt) | sort | uniq > ${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/ms_BNDregion_bed/${patient_barcode}/${patient_barcode}_ms_BNDregion.bed
    else
        echo "No .bed files found in the specified directory." >&2
        exit 1
    fi

    """
}



process upload_to_aws {
    
    tag "${patient_barcode}"

    label 'upload_to_bucket'

    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/ms_BNDregion_bed/${patient_barcode}", pattern: "*.bed", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/ms_BNDregion_bed/${patient_barcode}", pattern: "*{command,exitcode}*", mode: 'copy'
 
    input:

        tuple val(patient_barcode), val(repository), val(project), val(study), path(bed_file)

    output:
        tuple val(patient_barcode), path("${patient_barcode}_ms_BNDregion.bed"), path("*{command,exitcode}*", hidden:true)

    script:
    """
    aws configure set aws_access_key_id -
    aws configure set aws_secret_access_key -

    aws s3 cp ${bed_file} ${params.cloud_bucket}/${params.profile_name}/scratch/results/${study}/${project}/ms_BNDregion_bed/${patient_barcode}/
    """
}

process upload_to_gcp {
    
    tag "${patient_barcode}"

    label 'transfer_file_to_bucket'

    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/ms_BNDregion_bed/${patient_barcode}", pattern: "*.bed", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/ms_BNDregion_bed/${patient_barcode}", pattern: "*{command,exitcode}*", mode: 'copy'
 
    input:
        tuple val(patient_barcode), val(repository), val(project), val(study), path(bed_file)

    output:
        tuple val(patient_barcode), path("${patient_barcode}_ms_BNDregion.bed"), path("*{command,exitcode}*", hidden:true)
    
    script:
    """
    gsutil -u {project_name} -m cp ${bed_file} ${params.cloud_bucket}/${params.profile_name}/scratch/results/${study}/${project}/ms_BNDregion_bed/${patient_barcode}/
    echo 0 >> finish_flag.txt

    """
}



process slice_bam_pair {

    tag "${patient_barcode}"

    label 'slice_bam_pair'
    
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/ms_BNDregion_bam/${patient_barcode}/", pattern: "*", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/ms_BNDregion_bam/${patient_barcode}/", pattern: "*{command,exitcode}*", mode: 'copy'

    input:
        tuple val(patient_barcode), val(object_id), val(study), val(repository), val(project), path(bam), path(bai), path(bam_header), path(bed_file), val(tm_or_nm)
        file(genome_fasta)
        file(genome_index)
        file(genome_dict)

    output:
        tuple val(patient_barcode), val(object_id), val(study), val(repository), val(project), path(bam_header), path("${object_id}_ms_BNDregion.bam"), path("${object_id}_ms_BNDregion.bai"), path("env.txt"), path("*{command,exitcode}*", hidden:true)

    script:
    """
    #!/bin/bash
    gatk -version >> env.txt

    java -jar /opt/gatk/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar PrintReads \
        -R ${genome_fasta} \
        -I ${bam} \
        -L ${bed_file} \
        -O ${object_id}_ms_BNDregion.bam \
        --disable-tool-default-read-filters

    echo 0 >> finish_flag.txt
    ls -al -R . >> env.txt
    """
}

process slice_bam_pair_hmf {

    tag "${output_barcode}"

    label 'slice_bam_pair_bamtools_v13rc1'
    
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/ms_BNDregion_bam/${patient_barcode}/${output_barcode}/${object_id}", pattern: "*{txt,bam}", mode: 'copy'
    publishDir "${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/ms_BNDregion_bam/${patient_barcode}/${output_barcode}/${object_id}", pattern: "*{command,exitcode}*", mode: 'copy'

    input:
        tuple val(output_barcode), val(patient_barcode), val(study), val(repository), val(project), val(object_id), path(cram), val(file_path), path(bed_file), val(tm_or_nm)
        file(genome_fasta)
        file(genome_index)
        file(genome_dict)

    output:
        tuple val(output_barcode), val(patient_barcode), val(object_id), val(study), val(repository), val(project), path("*_cram_file_info.txt"), path("*.header.txt"), path("*.bam"), path("env.txt"), path("*{command,exitcode}*", hidden:true)

    script:
    """
    #!/bin/bash
    gsutil -u {project_name} ls -L ${file_path} >> ${cram}_cram_file_info.txt
    samtools view -H *.cram >> ${cram}.header.txt
    samtools index *.cram 

    java -Xmx50G -cp /opt/hmftools/bamtools/bam-tools_v1.3-rc.1.jar com.hartwig.hmftools.bamtools.slice.RegionSlicer \
        -output_dir ./ \
        -output_prefix ${object_id}_ms_BNDregion \
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
    } else if (key == 'hmf_hmftools_somatic_ms') {
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
                return [ output_barcode, tm_object_id,                               nm_object_id,                               repository, patient_barcode, project, study        ]}
            .combine(
                filtered_samples_ch.map{ output_barcode, aa_files -> return [output_barcode, aa_files]},by:[0])
            .map{ output_barcode, tm_object_id, nm_object_id, repository, patient_barcode, project, study, aa_files -> return [output_barcode, tm_object_id, nm_object_id, patient_barcode, repository, project, study, aa_files]}
}


def ch_local_manifest_pair_bed(manifest_file){
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
            return[ patient_barcode, repository, project, study]}
            .unique()
            .map{ patient_barcode, repository, project, study ->
                def ms_bed_file = file("${params.scratch_dir}/results/${params.profile_name}/${study}/${repository}/${project}/ms_BNDregion_bed/${patient_barcode}/${patient_barcode}_ms_BNDregion.bed", checkIfExists:true)
                return [patient_barcode, repository, project, study, ms_bed_file]
            }
}

def ch_local_manifest_cloud(manifest_file){
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
            return [ output_barcode, tm_object_id,                               nm_object_id, patient_barcode,              repository,                  project, study]}
}



def make_pair_bam(local_manifest, ch_access_bam_out){
    local_manifest
        .map{        output_barcode, tm_object_id, nm_object_id, patient_barcode, repository,  project, study -> 
            return [ output_barcode, tm_object_id, nm_object_id, patient_barcode                                ]
        }
        .combine(
            ch_access_bam_out.map{  output_barcode,                         study, repository, project, tm_bam, tm_bai, tm_bam_header, nm_bam, nm_bai, nm_bam_header->
                            return [output_barcode,                         study, repository, project, tm_bam, tm_bai, tm_bam_header, nm_bam, nm_bai, nm_bam_header]
            }, by:[0] 
        )
        .map{        output_barcode, tm_object_id, nm_object_id, patient_barcode, study, repository, project, tm_bam, tm_bai, tm_bam_header, nm_bam, nm_bai, nm_bam_header-> 
            def bed_file = file("${params.cloud_bucket}/local/scratch/results/${study}/${project}/ms_BNDregion_bed/${patient_barcode}/${patient_barcode}_ms_BNDregion.bed")
            return bed_file.exists() ? [output_barcode, tm_object_id, nm_object_id, patient_barcode, study, repository, project, tm_bam, tm_bai, tm_bam_header, nm_bam, nm_bai, nm_bam_header, bed_file] : null
        }
        .filter { it != null }
}




def ch_merge_patient_input(ch_extract_BNDregion_out){
    ch_extract_BNDregion_out
        .map{        output_barcode, tm_object_id, nm_object_id, patient_barcode, repository, project, study, bed_file -> 
            return [                                             patient_barcode, repository, project, study          ]
        }
        .unique()
}



def filter_empty_bed_files(input_channel){
    return input_channel.filter { tuple ->
        def (patient_barcode, repository, project, study, ms_bed_file) = tuple
        return ms_bed_file.size() > 0 }
}




def make_sliceBAM_input(ch_pair_bams) {
    def ch_tumor_bams = ch_pair_bams
        .map { output_barcode, tm_object_id, nm_object_id, patient_barcode, study, repository, project, tm_bam, tm_bai, tm_bam_header, nm_bam, nm_bai, nm_bam_header, bed_file ->
            def tm_or_nm = 'tumor'
            return [ patient_barcode, tm_object_id, study, repository, project, tm_bam, tm_bai, tm_bam_header, bed_file , tm_or_nm ]
        }
        .unique()

    def ch_normal_bams = ch_pair_bams
        .map { output_barcode, tm_object_id, nm_object_id, patient_barcode, study, repository, project, tm_bam, tm_bai, tm_bam_header, nm_bam, nm_bai, nm_bam_header, bed_file ->
            def tm_or_nm = 'normal'
            return [ patient_barcode, nm_object_id, study, repository, project, nm_bam, nm_bai, nm_bam_header, bed_file , tm_or_nm ]
        }
        .unique()

    def merged_channel = ch_tumor_bams.concat(ch_normal_bams)

    return merged_channel
        .map { patient_barcode, object_id, study, repository, project, bam, bai, bam_header, bed_file, tm_or_nm ->
            return [
                patient_barcode,
                object_id,
                study,
                repository,
                project,
                bam,
                bai,
                bam_header,
                bed_file,
                tm_or_nm
            ]
        }
        .unique()
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
                def bed_file = file("${params.cloud_bucket}/local/scratch/results/${study}/${project}/ms_BNDregion_bed/${patient_barcode}/${patient_barcode}_ms_BNDregion.bed")
                def tm_cram = file("${tm_file_path}",checkIfExists: true)
                def nm_cram = file("${nm_file_path}",checkIfExists: true)
                //def tm_crai = file("${tm_crai_path}",checkIfExists: true)
                //def nm_crai = file("${nm_crai_path}",checkIfExists: true)
 
            return [ output_barcode, patient_barcode, tm_object_id, nm_object_id, study, repository, project,  tm_cram, tm_file_path,  nm_cram, nm_file_path, bed_file]}
}

def ch_slice_pair_to_object_hmf(ch_pair_bams) {
    def ch_tumor_bams = ch_pair_bams
        .map { output_barcode, patient_barcode, tm_object_id, nm_object_id, study, repository, project,  tm_cram, tm_file_path, nm_cram, nm_file_path, bed_file ->
            def tm_or_nm = 'tumor'
            return [ output_barcode, study, repository, project, patient_barcode, tm_object_id, tm_cram, tm_file_path, bed_file, tm_or_nm ]
        }
        .unique()

    def ch_normal_bams = ch_pair_bams
        .map { output_barcode, patient_barcode, tm_object_id, nm_object_id, study, repository, project,  tm_cram, tm_file_path, nm_cram, nm_file_path, bed_file ->
            def tm_or_nm = 'normal'
            return [ output_barcode, study, repository, project, patient_barcode, nm_object_id, nm_cram, nm_file_path, bed_file, tm_or_nm  ]
        }
        .unique()

    def merged_channel = ch_tumor_bams.concat(ch_normal_bams)

    return merged_channel
        .map { output_barcode, study, repository, project, patient_barcode, object_id, cram, file_path, bed_file, tm_or_nm ->
            return [
                output_barcode,
                patient_barcode,
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
