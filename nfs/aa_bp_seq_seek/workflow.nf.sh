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

}
