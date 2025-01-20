// ######################################################################
//
// workflow_name: run_BND_region_bam
//
// ######################################################################
include {run_BNDregion_bed} from './BNDregion/BNDregion_bam.nf.sh'
include {run_BNDregion_bam} from './BNDregion/BNDregion_bam.nf.sh'
include {run_convert_to_cram} from './BNDregion/BNDregion_bam.nf.sh'
include {run_no_need_rm_bam} from './BNDregion/BNDregion_bam.nf.sh'
include {run_BNDregion_bam_hmf} from './BNDregion/BNDregion_bam.nf.sh'


// ######################################################################
//
// workflow_name: run_BND_region_bam_ms
//
// ######################################################################
// ## local
include {run_BNDregion_bed_ms} from './BNDregion_ms/BNDregion_bam_ms.nf.sh'
include {run_BNDregion_bed_ms_merge} from './BNDregion_ms/BNDregion_bam_ms.nf.sh'
include {run_merge_bed_to_cloud} from './BNDregion_ms/BNDregion_bam_ms.nf.sh'
// ## cloud
include {run_BNDregion_bam_ms} from './BNDregion_ms/BNDregion_bam_ms.nf.sh'
include {run_BNDregion_bam_ms_hmf} from './BNDregion_ms/BNDregion_bam_ms.nf.sh'


// ######################################################################
//
// workflow_name: aa_bp_seq_seek
//
// ######################################################################
include { aa_bp_seq_seek } from './aa_bp_seq_seek/aa_bp_seq_seek.nf.sh' addParams(workflow_name: "aa_bp_seq_seek", aa_workflow:"aaSuite_somatic_ss", aa_gain: 4.5, aa_cnsize_min: 50000, aa_downsample: 10)
include { aa_bp_seq_seek_pair_sh } from './aa_bp_seq_seek/aa_bp_seq_seek_pair_sh.nf.sh' addParams(workflow_name: "aa_bp_seq_seek_pair_sh", aa_workflow:"aaSuite_somatic_ss", aa_gain: 4.5, aa_cnsize_min: 50000, aa_downsample: 10)

// ======================================================================
// ======================================================================
//
// workflow
//
// ======================================================================
// ======================================================================
workflow {
    
    step2List = [
    'aa_bp_seq_seek']


    if (['aa_bp_seq_seek'].contains(params.step)) aa_bp_seq_seek(dna_bams)
    if (['aa_bp_seq_seek_pair_sh'].contains(params.step)) aa_bp_seq_seek_pair_sh()

//this part for sliced bam & collect bed file
    if (['transfer_controlled_data_to_gcp_bucket'].contains(params.step)) {
    
        transfer_controlled_data_to_gcp_bucket("transfer_controlled_data_to_gcp_bucket")
        //transfer_controlled_data_to_gcp_bucket.out.view()
    
    }
    
    if (['run_get_tmnm_pair_bams'].contains(params.step)) {
    
        run_get_tmnm_pair_bams("run_get_tmnm_pair_bams")
        //run_get_tmnm_pair_bams.out.view()
    
    }
    
    if (['hmftools_pipeline'].contains(params.step)) {
    
        run_get_tmnm_pair_bams("run_get_tmnm_pair_bams")
        //run_get_tmnm_pair_bams.out.paired_bam_out.view()
    
        hmftools_pipeline("hmftools_pipeline", run_get_tmnm_pair_bams.out.paired_bam_out)
    }
    
    if (['run_BNDregion_bed'].contains(params.step)) {
        run_BNDregion_bed("run_BNDregion_bed")
    }
    
    if (['run_BNDregion_bam'].contains(params.step)) {
        run_get_tmnm_pair_bams("run_get_tmnm_pair_bams")
        run_BNDregion_bam("run_BNDregion_bam", run_get_tmnm_pair_bams.out.paired_bam_out)
    }
    
    
    if (['run_convert_to_cram'].contains(params.step)) {
        //run_get_tmnm_pair_bams("run_get_tmnm_pair_bams")
        
        run_convert_to_cram("run_convert_to_cram")
    }
    
    if (['run_no_need_rm_bam'].contains(params.step)) {
        run_no_need_rm_bam("run_no_need_rm_bam")
    }
    
    
    if (['run_BNDregion_bed_ms'].contains(params.step)) {
        run_BNDregion_bed_ms("run_BNDregion_bed_ms")
    }
    
    if (['run_BNDregion_bed_ms_merge'].contains(params.step)) {
        run_BNDregion_bed_ms("run_BNDregion_bed_ms")
    
        run_BNDregion_bed_ms_merge("run_BNDregion_bed_ms_merge", run_BNDregion_bed_ms.out.ch_extract_BNDregion_out)
    }
    
    if (['run_merge_bed_to_cloud'].contains(params.step)) {
        run_merge_bed_to_cloud("run_merge_bed_to_cloud")
    }
    
    if (['run_BNDregion_bam_ms'].contains(params.step)) {
        run_get_tmnm_pair_bams("run_get_tmnm_pair_bams")
    
        run_BNDregion_bam_ms("run_BNDregion_bam_ms", run_get_tmnm_pair_bams.out.paired_bam_out)
    }
    
    //    if (['run_BNDregion_bam_hmf'].contains(params.step)) {
    //        run_BNDregion_bam_hmf("run_BNDregion_bam_hmf")
    //    }
    
    if (['run_BNDregion_bam_hmf'].contains(params.step)) {
        run_BNDregion_bam_hmf("run_BNDregion_bam_hmf")
    }
    
    if (['run_BNDregion_bam_ms_hmf'].contains(params.step)) {
        run_BNDregion_bam_ms_hmf("run_BNDregion_bam_ms_hmf")
    }
}
