version 1.0

import "tasks.wdl"

workflow test_ld {

  bgen chr1_bgen = {
    "bgen": "../ukbiobank/array_imputed/ukb_imp_chr1_v3.bgen",
    "index": "../ukbiobank/array_imputed/ukb_imp_chr1_v3.bgen.bgi"
  }
  File mfi = "../ukbiobank/array_imputed/ukb_mfi_chr1_v3.txt"
  pfiles chr1_pfiles = {
    "pgen": "../ukbiobank/array_imputed/pfile_converted/chr1.pgen",
    "psam": "../ukbiobank/array_imputed/pfile_converted/chr1.psam",
    "pvar": "../ukbiobank/array_imputed/pfile_converted/chr1.pvar",
  }
  File samples = "data/platelet_count/combined_unrelated.sample"
  Int n_samples = 50000
  region bounds = {
    "chrom": 1,
    "start": 0,
    "end": 136062 #1000 variants in
  }

  call tasks.get_training_samples { input : 
    samples = samples,
    num_samples = n_samples
  }

  call tasks.ldstore { input : 
    input_bgen = chr1_bgen,
    mfi = mfi,
    bounds = bounds,
    n_samples = n_samples,

    sample_file = get_training_samples.training_samples
  }

  call tasks.plink_ld { input :
    input_pfiles = chr1_pfiles,
    bounds = bounds,
    sample_file = get_training_samples.training_samples
  }
}
