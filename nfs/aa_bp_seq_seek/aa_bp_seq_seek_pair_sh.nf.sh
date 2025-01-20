/*
Updated: 2024-12-12
Edit: 2025-01-16 by Jiwon 
Description: This version of aa_bp_seq_seek is for analyzing AA pair mode results from the cloud pipeline, modified to separate and analyze tm and nm in single pair mode.
*/

/*
=================================================================
=================================================================

workflow

=================================================================
=================================================================
*/

workflow aa_bp_seq_seek_pair_sh {
	main:
        println "workflow_name = ${params.workflow_name}"
        println "aa_gain = ${params.aa_gain}"
        println "aa_cnsize_min = ${params.aa_cnsize_min}"
        println "aa_downsample = ${params.aa_downsample}"
        println "cohort_name = ${params.cohort_name}"

        R_ch = Channel.value(file("${params.nf_home}/R/nextflow"))
        py_ch = Channel.value(file("${params.nf_home}/python"))
        
        ch_svaba_dbsnp_vcf = Channel.value(file(params.svaba_dbsnp_vcf))
        if (params.cohort_name=="HMF") {
	        ch_genome_fasta = Channel.value(file(params.hmftools_ref_fasta))
	        ch_genome_index = Channel.value(file(params.hmftools_ref_fai))
	        ch_genome_dict 	= Channel.value(file(params.hmftools_ref_dict))
	        ch_genome_sa 	= Channel.value(file(params.hmftools_ref_sa))
	        ch_genome_bwt 	= Channel.value(file(params.hmftools_ref_bwt))
	        ch_genome_ann	= Channel.value(file(params.hmftools_ref_ann))        
	        ch_genome_amb 	= Channel.value(file(params.hmftools_ref_amb))
	        ch_genome_pac 	= Channel.value(file(params.hmftools_ref_pac))

        } else{
	        ch_genome_fasta = Channel.value(file(params.genome_fasta))
	        ch_genome_index = Channel.value(file(params.genome_index))
	        ch_genome_dict 	= Channel.value(file(params.genome_dict))
	        ch_genome_sa 	= Channel.value(file(params.genome_sa))
	        ch_genome_bwt 	= Channel.value(file(params.genome_bwt))
	        ch_genome_ann 	= Channel.value(file(params.genome_ann))
	        ch_genome_amb 	= Channel.value(file(params.genome_amb))
	        ch_genome_pac 	= Channel.value(file(params.genome_pac))
	        //ch_genome_alt 	= Channel.value(file(params.genome_alt))
        }

    // ==== indexing bam
    ch_pairs_bam = make_samples_only_bam ( params.dna_legacy_file_sheet )
   	// OUTPUT: [aliquot, study, project, repository, tm_object_id, tm_bam,         nm_object_id, nm_bam,         bed_path, aa_files ]
		ch_pairs = index_bam(ch_pairs_bam)
   	// OUTPUT: [aliquot, study, project, repository, tm_object_id, tm_bam, tm_bai, nm_object_id, nm_bam, nm_bai, bed_path, aa_files ]
		
		ch_input = ch_pairs
			.map{   aliquot, study, project, repository, tm_object_id, tm_bam, tm_bai, nm_object_id, nm_bam, nm_bai, bed_path, aa_files ->
			return[ aliquot, study, project, repository,                                                                       aa_files]}
		//OUTPUT: aliquot, study, project, repository, aa_files
		//ch_input.view()	

		// ==== aa_bp_ss: process 1
    count_aa_amp_num(ch_input)
		//OUTPUT: aliquot, study, project, repository, aa_files, aa_sum_file_path_txt, amp_num_txt, amp_num, cmds
		//count_aa_amp_num.out.view()

		// ==== aa_bp_ss: process 2-1 or 2-2 - check amplicon existance
		amp = count_aa_amp_num.out.map{ aliquot, study, project, repository, aa_files, aa_sum_file_path_txt, amp_num_txt, amp_num, cmds ->
                 				return[ aliquot, study, project, repository, aa_files,            					      amp_num ] }
		//amp.view()
    no_aa_amp(amp)
		// OUTPUT: aliquot, study, project, repository, aa_files, no_amplicon_flag, cmds
		
		make_input_for_breakpoints_to_bed(amp, R_ch)
		// OUTPUT: aliquot, study, project, repository, aa_files, amp_num, bpinput_txt, aainterval_txt, env, cmds
    	
		input_breakpoints_to_bed = make_input_for_breakpoints_to_bed.out
	      							    .map {  aliquot, study, project, repository, aa_files, amp_num, bpinput_txt, aainterval_txt, env, cmds -> 
	 							    return[     aliquot, study, project, repository, aa_files,          bpinput_txt, aainterval_txt]}
    get_aa_bp(input_breakpoints_to_bed)
		//OUTPUT: aliquot, study, project, repository, aa_files, bp_bed, env, cmds

	  // process 4-1 or 4-2
    output_breakpoints_to_bed = get_aa_bp.out
       				      .map{    aliquot, study, project, repository, aa_files, bp_bed_info, env, cmds -> 
	  	       	      return [ aliquot, study, project, repository,           bp_bed_info ]}
		//output_breakpoints_to_bed.view()


    // ==== aa_bp_ss: local assembly SvABA (remove:ch_genome_alt)
		run_svaba_ss_targeted_local_assembly_out=run_svaba_ss_targeted_local_assembly(ch_pairs, ch_genome_fasta, ch_genome_index, ch_genome_dict,
     			 		     				ch_genome_sa, ch_genome_bwt, ch_genome_ann, ch_genome_amb, ch_genome_pac,
         			 	     			ch_svaba_dbsnp_vcf)
		//OUTPUT: aliquot, study, project, repository, tm_object_id, nm_object_id, BND_bed, aa_files, svaba_somatic_sv_vcf, svaba_unfiltered_somatic_sv_vcf, svaba_germline_sv_vcf, svaba_unfiltered_germline_sv_vcf, svaba_output, bps_txt, env, cmds
		//run_svaba_ss_targeted_local_assembly_out.view()
        
	  // ==== aa_bp_ss: process 6 (preliminary process)

		input_align_aa_bp_to_svaba = output_breakpoints_to_bed
     						.combine(run_svaba_ss_targeted_local_assembly_out, by:[0,1,2,3])
								.map{    aliquot, study, project, repository, bp_bed_info, tm_object_id, nm_object_id, BND_bed, aa_files, svaba_somatic_sv_vcf, svaba_unfiltered_somatic_sv_vcf, svaba_germline_sv_vcf, svaba_unfiltered_germline_sv_vcf, svaba_output, bps_txt, env, cmds ->
    		        return [ aliquot, study, project, repository, bp_bed_info, tm_object_id, nm_object_id,          aa_files, svaba_somatic_sv_vcf, svaba_unfiltered_somatic_sv_vcf, svaba_germline_sv_vcf, svaba_unfiltered_germline_sv_vcf ]}
    //input_align_aa_bp_to_svaba.view()
    align_aa_bp_to_svaba(input_align_aa_bp_to_svaba, R_ch)


	  // process 7 ??
		//copy_script(ch_input, R_ch)


		// ==== POST_AA_BP_SS: vcf_to_bedpe
		input_vcf_to_bedpe = run_svaba_ss_targeted_local_assembly_out
    	       				      .map{    aliquot, study, project, repository, tm_object_id, nm_object_id, BND_bed, aa_files, svaba_somatic_sv_vcf, svaba_unfiltered_somatic_sv_vcf, svaba_germline_sv_vcf, svaba_unfiltered_germline_sv_vcf, svaba_output, bps_txt, env, cmds  -> 
	    			       	      return [ aliquot, study, project, repository, tm_object_id, nm_object_id,          aa_files, svaba_somatic_sv_vcf, svaba_unfiltered_somatic_sv_vcf, svaba_germline_sv_vcf, svaba_unfiltered_germline_sv_vcf,               bps_txt ]}
		vcf_to_bedpe(input_vcf_to_bedpe, py_ch)
		//OUTPUT:tuple val(aliquot_barcode), path("*.svaba.somatic.sv.vcf.bedpe.csv"), path("*{bedpe,log}"), path(bps_txt), path("env.txt"), path("*{command,exitcode}*",hidden:true)


		post_aabpss_input = align_aa_bp_to_svaba.out
			.map{ aliquot, study, project, repository, tm_object_id, nm_object_id, aa_files, svaba_somatic_sv_vcf, svaba_unfiltered_somatic_sv_vcf, somatic_matched_svaba_annotated_aa, germline_matched_svaba_annotated_aa, matched_svaba_aa, matched_unfiltered_svaba_annotated_aa, matched_unfiltered_svaba_aa, env, cmds -> 
			return [aliquot, study, project, repository, tm_object_id, nm_object_id, somatic_matched_svaba_annotated_aa ]}
			.combine(vcf_to_bedpe.output.map{aliquot, svaba_somatic_bedpe, files, bps_txt, env, cmds -> return[aliquot, svaba_somatic_bedpe, bps_txt]}, by:[0])
			.map{   aliquot, study, project, repository, tm_object_id, nm_object_id, somatic_matched_svaba_annotated_aa, svaba_somatic_bedpe, bps_txt -> 
			return [aliquot, study, project, repository, tm_object_id, nm_object_id, somatic_matched_svaba_annotated_aa, svaba_somatic_bedpe, bps_txt]}

		post_aabpss(post_aabpss_input, py_ch)
}


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
// Pair-mode analysis 
//
//
// --------------------------------------------------------------
process count_aa_amp_num {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/aa_amp_num", pattern: "*.{txt}", mode: 'copy'
	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/aa_amp_num", pattern: "*{command,exitcode}*", mode: 'copy'

	//label "get_aa_bp"

	input: 
		tuple val(aliquot_barcode), val(study), val(project), val(repository), path(aa_files)

	output:
		tuple val(aliquot_barcode), val(study), val(project), val(repository), path(aa_files), path("aa_summary_file_path.txt"), path("amplicon_num.txt"), env(amp_num), path("*{command,exitcode}*",hidden:true)

	shell:
		'''
		find  !{aa_files}/ -name "!{aliquot_barcode}_summary.txt" > ./aa_summary_file_path.txt
		cat ./aa_summary_file_path.txt | xargs -I{} sh -c 'head $0 -n 1 | cut -d " " -f3 > ./amplicon_num.txt' {} ;

		amp_num=`cat ./amplicon_num.txt`

		'''

}

process no_aa_amp {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/no_aa_amp", pattern: "*{flag}", mode: 'copy'
	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/no_aa_amp", pattern: "*{command,exitcode}*", mode: 'copy'
	
	label "get_aa_bp"
	
	input: 
		tuple val(aliquot_barcode), val(study), val(project), val(repository), path(aa_files), env(amp_num)

	output:
		tuple val(aliquot_barcode), val(study), val(project), val(repository), path(aa_files), path(no_amplicon_flag), path("*{command,exitcode}*",hidden:true)
	
	when:
		amp_num == "0"
	
	script:

		"""
		#!/bin/bash

		touch no_amplicon_flag

		"""

}

process make_input_for_breakpoints_to_bed {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/get_aa_bp", pattern: "*.{txt}", mode: 'copy'
	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/get_aa_bp", pattern: "*{command,exitcode}*", mode: 'copy'

	label "aa_bp_seq_seek"

	input: 
		tuple  val(aliquot_barcode), val(study), val(project), val(repository), path(aa_files), val(amp_num)
		file(R_dir)

	output:
		tuple val(aliquot_barcode), val(study), val(project), val(repository), path(aa_files), val(amp_num), path("${aliquot_barcode}_bp.input.txt"), path("${aliquot_barcode}_aa_summary_intervals.txt"), path("*env.txt"), path("*{command,exitcode}*",hidden:true)
	
	when:
		amp_num != "0"

	script:

		"""
		#!/bin/bash

			Rscript ${R_dir}/input_for_breakpoints_to_bed_script.R \
		    --aa_output_path ${aa_files} \
			--aliquot_barcode ${aliquot_barcode} \
		 	--output_path ./ && \
		 	ls -al -R ./ >> make_input_for_breakpoints_to_bed_env.txt
		

		"""

}


process get_aa_bp {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/get_aa_bp", pattern: "*.{txt,bed}", mode: 'copy'
	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/get_aa_bp", pattern: "*{command,exitcode}*", mode: 'copy'

	label "get_aa_bp"

	input:
		tuple val(aliquot_barcode), val(study), val(project), val(repository), path(aa_files), path(bpinput_txt), path(aainterval_txt)

	output:
		tuple val(aliquot_barcode), val(study), val(project), val(repository), path(aa_files), path("${aliquot_barcode}_bp.input_breakpoints.bed"), path("*env.txt"), path("*{command,exitcode}*",hidden:true)
	
	
	script:

		"""
		#!/bin/bash

		aainterval=`cat ${aainterval_txt}`

		python3 /home/breakpoints_to_bed.py \
		    -i ${bpinput_txt} \
        -r \$aainterval --add_chr_tag && \
        ls -al -R ./ >> get_aa_bp_env.txt
		"""

}


process no_aa_bp {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/no_aa_bp", pattern: "*{txt,flag}", mode: 'copy'
	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/no_aa_bp", pattern: "*{command,exitcode}*", mode: 'copy'

	label "get_aa_bp"

	input:
		tuple val(aliquot_barcode), val(study), val(project), val(repository), path(aa_files), path(bp_bed)

	output:
		tuple val(aliquot_barcode), val(study), val(project), val(repository), path(aa_files), path(no_aa_bp_flag), path("*{command,exitcode}*",hidden:true)
	
	when:
		bp_bed.size() == 0
	
	script:

		"""
		#!/bin/bash
		touch no_aa_bp_flag
		"""
}


// preliminary step
process align_aa_bp_to_svaba {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/align_aa_bp_to_svaba", pattern: "*.{txt,tsv}", mode: 'copy'
	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/align_aa_bp_to_svaba", pattern: "*{command,exitcode}*", mode: 'copy'

	label "aa_bp_seq_seek"

	input:
		tuple val(aliquot_barcode), val(study), val(project), val(repository), path(bp_bed), val(tm_object_id), val(nm_object_id), path(aa_files), path(svaba_somatic_sv_vcf), path(svaba_unfiltered_somatic_sv_vcf), path(svaba_germline_sv_vcf), path(svaba_unfiltered_germline_sv_vcf)
		file(R_dir)

	output:
		tuple val(aliquot_barcode), val(study), val(project), val(repository),               val(tm_object_id), val(nm_object_id), path(aa_files), path(svaba_somatic_sv_vcf), path(svaba_unfiltered_somatic_sv_vcf), path("*somatic.sv.vcf_matched_svaba_annotated_aa.tsv"), path("*germline.sv.vcf_matched_svaba_annotated_aa.tsv"), path("*_matched_svaba_aa.tsv"), path("*_matched_unfiltered_svaba_annotated_aa.tsv"), path("*_matched_unfiltered_svaba_aa.tsv"), path("env.txt"), path("*{command,exitcode}*",hidden:true)
	
	script:

		"""
    	#!/bin/bash
    	
    		Rscript ${R_dir}/matching_aa_svaba_ss.R \
		    --breakpoints_bed ${bp_bed} \
		    --aa_output_path ${aa_files} \
		    --svaba_sv_vcf ${svaba_somatic_sv_vcf} \
		    --svaba_unfiltered_sv_vcf ${svaba_unfiltered_somatic_sv_vcf} \
			  --aliquot_barcode ${aliquot_barcode} \
		   	--output_path ./ && \

    		Rscript ${R_dir}/matching_aa_svaba_ss.R \
		    --breakpoints_bed ${bp_bed} \
		    --aa_output_path ${aa_files} \
		    --svaba_sv_vcf ${svaba_germline_sv_vcf} \
		    --svaba_unfiltered_sv_vcf ${svaba_germline_sv_vcf} \
			  --aliquot_barcode ${aliquot_barcode} \
		 	  --output_path ./ && \

		 	ls -al -R ./ >> env.txt

    	"""
}


process copy_script {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/${object_id}/scripts", pattern: "*", mode: 'copy'

	input:
		val(aliquot_barcode)
		file(R_dir)

	output:
		tuple val(aliquot_barcode), path("*"), path("*{command,exitcode}*",hidden:true)
	
	when:
		amp_num != "0"
	
	script:

		"""
    	#!/bin/bash
    	
    		cp ${R_dir}/input_for_breakpoints_to_bed_script.R ./`date +%Y.%m.%d`_input_for_breakpoints_to_bed_script.R && \
    		cp ${R_dir}/input_for_svaba.R ./`date +%Y.%m.%d`_input_for_svaba.R && \
    		cp ${R_dir}/matching_aa_svaba.R ./`date +%Y.%m.%d`_matching_aa_svaba.R && \
    		cp ${params.nf_home}/nfs/aa_bp_seq_seek/${params.step}.nf.sh ./`date +%Y.%m.%d`_${params.step}.nf.sh
    	"""

}

process run_svaba_ss_targeted_local_assembly {
	
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/svaba_ss/2kb/local_assembly", pattern: "*{vcf,txt,bam,gz,log}", mode: 'copy'
	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/svaba_ss/2kb/local_assembly/bed", pattern: "*{bed}", mode: 'copy'
	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/svaba_ss/2kb/local_assembly", pattern: "*{command,exitcode}*", mode: 'copy'
	
	label "svaba"

	input:
		tuple val(aliquot_barcode), val(study), val(project), val(repository), val(tm_object_id), path(tm_bam), path(tm_bai), val(nm_object_id), path(nm_bam), path(nm_bai), path(BND_bed), path(aa_files)
		file(genome_fasta)
        file(genome_index)
        file(genome_fasta_dict)
        file(genome_sa) 
        file(genome_bwt) 
        file(genome_ann) 
        file(genome_amb) 
        file(genome_pac) 
        file(svaba_dbsnp_vcf)

	output:
		tuple val(aliquot_barcode), val(study), val(project), val(repository), val(tm_object_id), val(nm_object_id), path(BND_bed), path(aa_files),path("*.svaba.somatic.sv.vcf"), path("*.svaba.unfiltered.somatic.sv.vcf"),path("*.svaba.germline.sv.vcf"), path("*.svaba.unfiltered.germline.sv.vcf"), path("*{vcf,txt.gz,contigs.bam,log}"), path("*.bps.txt.gz"), path("env.txt"), path("*{command,exitcode}*",hidden:true)

	
	script:
		"""
    	#!/bin/bash

    	svaba run \
    	-t ${tm_bam} \
    	-n ${nm_bam} \
    	-p ${params.svaba_threads} \
    	-D ${svaba_dbsnp_vcf} \
    	-a ${aliquot_barcode} \
    	-G ${genome_fasta} \
	    -k ${BND_bed} && 

    	ls -al -R ./ >> env.txt

    	"""
}



process vcf_to_bedpe {
    
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/svaba_ss/2kb/local_assembly/bedpe", pattern: "*{bedpe,txt}", mode: 'copy'
	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/svaba_ss/2kb/local_assembly/bedpe", pattern: "*{command,exitcode}*", mode: 'copy'
	    			       	     
	label "bp_homology_python"

	input:
		tuple val(aliquot_barcode), val(study), val(project), val(repository), val(tm_object_id), val(nm_object_id), path(aa_files), path(svaba_somatic_sv_vcf), path(svaba_unfiltered_somatic_sv_vcf), path(svaba_germline_sv_vcf), path(svaba_unfiltered_germline_sv_vcf), path(bps_txt)
		file(py_dir)

	when:
		svaba_somatic_sv_vcf.exists() && svaba_somatic_sv_vcf.size() > 0
		
	output:
		tuple val(aliquot_barcode), path("*.svaba.somatic.sv.vcf.bedpe.csv"), path("*.bedpe"), path(bps_txt), path("env.txt"), path("*{command,exitcode}*",hidden:true)

    script:
    	
    	"""
    	python3 ${py_dir}/vcf_to_csv_and_bedpe_ss.py ${svaba_somatic_sv_vcf} "${aliquot_barcode}" &&\
    	python3 ${py_dir}/vcf_to_csv_and_bedpe_ss.py ${svaba_unfiltered_somatic_sv_vcf} "${aliquot_barcode}" &&\
    	python3 ${py_dir}/vcf_to_csv_and_bedpe_ss.py ${svaba_germline_sv_vcf} "${aliquot_barcode}" &&\
    	python3 ${py_dir}/vcf_to_csv_and_bedpe_ss.py ${svaba_unfiltered_germline_sv_vcf} "${aliquot_barcode}" &&\
	 	
	 	ls -al -R ./ >> env.txt
    	"""
}


process post_aabpss {
    
	tag "${aliquot_barcode}"

	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/svaba_ss/2kb/local_assembly/post_aabpss", pattern: "*.{csv,txt}", mode: 'copy'
	publishDir "${params.scratch_dir}/results/${params.workflow_name}/${params.aa_workflow}/${params.genome}/${study}/${repository}/${project}/${aliquot_barcode}/svaba_ss/2kb/local_assembly/post_aabpss", pattern: "*{command,exitcode}*", mode: 'copy'
	
	label "bp_homology_python"
	    			       	     
	input:
		tuple val(aliquot_barcode), val(study), val(project), val(repository), val(tm_object_id), val(nm_object_id), path(matched_svaba_annotated_aa), path(svaba_somatic_bedpe), path(bps_txt)
		file(py_dir)
	
	when:
		svaba_somatic_bedpe.exists() && svaba_somatic_bedpe.size() > 0
	
	output:
		tuple val(aliquot_barcode), path("*.csv"), path("env.txt"), path("*{command,exitcode}*",hidden:true)

    script:
    	
    	"""
    	gunzip -c ${bps_txt} > *.bps.txt
    	python3 ${py_dir}/post_aa_bp_ss.py ${aliquot_barcode} ${svaba_somatic_bedpe} *.bps.txt ${matched_svaba_annotated_aa}  ${params.gr4_paths[study]} ${params.amp_paths[study]} \
	 	
	 	  ls -al -R ./ >> env.txt
	
    	"""
}

process index_bam {
    
    tag "${aliquot_barcode}"
	label "bp_homology_python"

    input:
	    tuple val(aliquot_barcode), val(study), val(project), val(repository), val(tm_object_id), path(tm_bam), val(nm_object_id), path(nm_bam), path(bed_path), path(aa_file)

    output:
	    tuple val(aliquot_barcode), val(study), val(project), val(repository), val(tm_object_id), path(tm_bam), path("${tm_bam}.bai"), val(nm_object_id), path(nm_bam), path("${nm_bam}.bai"), path(bed_path), path(aa_file)

    script:
    """
    samtools index ${tm_bam}
    samtools index ${nm_bam}

 	  ls -al -R ./ >> env.txt

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
                            row.study,
                            row.project,
                            row.repository,
                            row.tm_object_id,
                            row.tm_file_path,
                            row.tm_bai_file_path,
                            row.nm_object_id,
                            row.nm_file_path,
                            row.nm_bai_file_path,
                            row.BNDregion_bed_file_path,
                            row.cp_output_dir,
                            row.action ) }
        .map{          aliquot, study, project, repository, tm_object_id, tm_bam_path, tm_bai_path, nm_object_id, nm_bam_path, nm_bai_path, bed_path, aa_files, action ->
            def tm_bam = file(tm_bam_path).exists() ? file(tm_bam_path): "File not found"
            def tm_bai = file(tm_bai_path).exists() ? file(tm_bai_path): "File not found"
            def nm_bam = file(nm_bam_path).exists() ? file(nm_bam_path): "File not found"
            def nm_bai = file(nm_bai_path).exists() ? file(nm_bai_path): "File not found"
            def aa_file = file(aa_files).exists() ? file(aa_files): "File not found"
	        return     [aliquot, study, project, repository, tm_object_id, tm_bam, tm_bai, nm_object_id, nm_bam, nm_bai, bed_path, aa_file, action ] }
        .filter{        aliquot, study, project, repository, tm_object_id, tm_bam, tm_bai, nm_object_id, nm_bam, nm_bai, bed_path, aa_file, action ->
			    action == "run" }
		.map{           aliquot, study, project, repository, tm_object_id, tm_bam, tm_bai, nm_object_id, nm_bam, nm_bai, bed_path, aa_file, action ->
				return [aliquot, study, project, repository, tm_object_id, tm_bam, tm_bai, nm_object_id, nm_bam, nm_bai, bed_path, aa_file ]}
        .unique()
}

def make_samples_only_bam ( dna_legacy_file_sheet  ) {
    ch_samples = Channel
        .fromPath(dna_legacy_file_sheet)
        .splitCsv(header:true)
        .map{ row-> tuple(  row.aa_barcode,
                            row.study,
                            row.project,
                            row.repository,
                            row.tm_object_id,
                            row.tm_file_path,
                            row.nm_object_id,
                            row.nm_file_path,
                            row.BNDregion_bed_file_path,
                            row.cp_output_dir,
                            row.action ) }
        .map{          aliquot, study, project, repository, tm_object_id, tm_bam_path, nm_object_id, nm_bam_path, bed_path, aa_files, action ->
            def tm_bam = file(tm_bam_path).exists() ? file(tm_bam_path): "File not found"
            def nm_bam = file(nm_bam_path).exists() ? file(nm_bam_path): "File not found"
            def aa_file = file(aa_files).exists() ? file(aa_files): "File not found"
	        return     [aliquot, study, project, repository, tm_object_id, tm_bam, nm_object_id, nm_bam, bed_path, aa_file, action ] }
        .filter{        aliquot, study, project, repository, tm_object_id, tm_bam, nm_object_id, nm_bam, bed_path, aa_file, action ->
			    action == "run" }
		.map{           aliquot, study, project, repository, tm_object_id, tm_bam, nm_object_id, nm_bam, bed_path, aa_file, action ->
				return [aliquot, study, project, repository, tm_object_id, tm_bam, nm_object_id, nm_bam, bed_path, aa_file ]}
        .unique()
}


def ch_pair_to_object(ch_pair_bams) {
    def ch_tumor_bams = ch_pair_bams
        .map { aliquot, study, project, repository, tm_object_id, tm_bam, tm_bai, nm_object_id, nm_bam, nm_bai, bed_path, aa_file ->
            def tm_or_nm = 'tumor'
            return [ aliquot, study, project, repository, tm_object_id, tm_bam, tm_bai, bed_path, tm_or_nm, aa_file ]
        }
        .unique()

    def ch_normal_bams = ch_pair_bams
        .map { aliquot, study, project, repository, tm_object_id, tm_bam, tm_bai, nm_object_id, nm_bam, nm_bai, bed_path, aa_file ->
            def tm_or_nm = 'normal'
            return [ aliquot, study, project, repository, nm_object_id, nm_bam, nm_bai, bed_path, tm_or_nm, aa_file  ]
        }
        .unique()

    def merged_channel = ch_tumor_bams.concat(ch_normal_bams)

    return merged_channel
        .map{        aliquot, study, project, repository, object_id, bam_path, bai_path, bed_path, tm_or_nm, aa_file ->
            return [ aliquot, study, project, repository, object_id, bam_path, bai_path, bed_path, tm_or_nm, aa_file ]}
        .unique()
}
