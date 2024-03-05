version 1.0

task run_susie {

	input {
    String script_dir
    File script = "~{script_dir}/sum_stats_susie.r"
  }

  command <<<
    Rscript ~{script}
  >>>

	runtime {
		docker: "quay.io/thedevilinthedetails/work/susie:v0.12.35"
		dx_timeout: '4h'
		memory: mem
		continueOnReturnCode: [0, 79] # allow for timeouts, this will be handled in the retry workflow
	}

}

workflow finemap_region {

  input {
    File plink_gwas # must have columns #CHROM POS REF ALT BETA SE P MACH_R2 A1_FREQ ERRCODE, other columns are optional
    File? precalculated_ld
    Int n_samples
  }

  # get new sample list

  # take existing fine-mapping regions, expand by 500kb in each direction

  # rerun GWAS in those regions

  # reconstruct fine-mapping regions, check to see if they are out of bounds, if so, truncate

  # choose how to filter variants 

  # create LD

  # write finemap input file

  # call finemap

  # call susie
}
