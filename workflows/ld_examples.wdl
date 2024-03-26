version 1.0

import "tasks.wdl"

workflow ld_examples {
  call tasks.get_training_samples { input :
    samples = "data/platelet_count/combined_unrelated.sample",
    num_samples = 50000,
  }

  call tasks.get_validation_subsample_replicates { input :
    samples = "data/platelet_count/combined_unrelated.sample",
    num_samples = 50000,
    replicate = 1
  }

  Array[region] many_bounds = [
    {
      "chrom": 1,
      "start": 15250110,
      "end": 15506160
    },
    {
      "chrom": 1,
      "start": 31116886,
      "end": 31366886
    },
    {
      "chrom": 1,
      "start": 32207350,
      "end": 32487347
    },
    {
      "chrom": 1,
      "start": 44514722,
      "end": 44791297 
    },
    {
      "chrom": 1,
      "start": 47545525,
      "end": 47832027
    }
  ]

  scatter (chrom in range(22)) {
    pfiles pfiles_ = {
      "pgen": "../ukbiobank/array_imputed/pfile_converted/chr~{chrom+1}.pgen",
      "pvar": "../ukbiobank/array_imputed/pfile_converted/chr~{chrom+1}.pvar",
      "psam": "../ukbiobank/array_imputed/pfile_converted/chr~{chrom+1}.psam",
    }
  }

  scatter (bounds in many_bounds) {
    Int chrom_p_one = bounds.chrom - 1
    call tasks.plink_ld as training_ld { input :
      input_pfiles = pfiles_[chrom_p_one],
      bounds = bounds,

      sample_file = get_training_samples.training_samples,
      prefix = "~{bounds.chrom}_~{bounds.start}_~{bounds.end}_training"
    }

    call tasks.plink_ld as validation_ld { input :
      input_pfiles = pfiles_[chrom_p_one],
      bounds = bounds,

      sample_file = get_validation_subsample_replicates.replicate_samples,
      prefix = "~{bounds.chrom}_~{bounds.start}_~{bounds.end}_validation"
    }
  }

  output {
    Array[File] training_lds = training_ld.ld
    Array[File] training_logs = training_ld.log
    Array[File] training_freqs = training_ld.freq
    Array[File] validation_lds = validation_ld.ld
    Array[File] validation_logs = validation_ld.log
    Array[File] validation_freqs = validation_ld.freq
  }
}
