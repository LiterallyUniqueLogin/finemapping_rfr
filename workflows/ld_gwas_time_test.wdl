version 1.0

import "tasks.wdl"
import "prep_ukb_inputs.wdl"
import "../../ukbiobank/workflow/gwas_wdl/gwas_tasks.wdl"

workflow ld_gwas_time_test {

  call prep_ukb_inputs.prep_ukb_inputs { input :
    project_ukb_script_dir = 'scripts/ukb_prep',
    old_ukb_script_dir = '../ukbiobank',
    phenotype_name = 'platelet_count',
    is_binary = false,
    phenotype_id = 30080,
    categorical_covariate_names = ['platelet_count_device_id'],
    categorical_covariate_ids = ['30083'],
  }

  call tasks.get_training_samples { input :
    samples = prep_ukb_inputs.samples_for_phenotype_no_header,
    num_samples = 50000,
  }

  scatter (gwas_chrom_zero_indexed in range(22)) {
    Int gwas_chrom_one_indexed = gwas_chrom_zero_indexed + 1

    pfiles pfiles_ = {
      'pgen': '../ukbiobank/array_imputed/pfile_converted/chr~{gwas_chrom_one_indexed}.pgen',
      'pvar': '../ukbiobank/array_imputed/pfile_converted/chr~{gwas_chrom_one_indexed}.pvar',
      'psam': '../ukbiobank/array_imputed/pfile_converted/chr~{gwas_chrom_one_indexed}.psam',
    }

    call tasks.plink_gwas as chromosomal_plink_gwas { input :
      input_pfiles = pfiles_,
      phenotype_name = 'platelet_count',
      pheno_covar_file = prep_ukb_inputs.unrelated_pheno_df,
      quantile_normalize = true,
      sample_file = get_training_samples.training_samples,
      prefix = 'platelet_count_training_gwas_chr~{gwas_chrom_one_indexed}'
    }
  }

  call tasks.concatenate_tsvs as plink_gwas { input :
    tsvs = chromosomal_plink_gwas.gwas,
    prefix = 'platelet_count_training_gwas'
  }

  call gwas_tasks.generate_finemapping_regions { input :
    script_dir = '../ukbiobank',
    chr_lens = '../ukbiobank/misc_data/genome/chr_lens.txt',
    phenotype = 'platelet_count',
    snp_assoc_results = plink_gwas.tsv,
    prefix = 'platelet_count_training_'
  }
  Array[Array[String]] finemapping_regions_tsv_with_header = read_tsv(generate_finemapping_regions.data)

  scatter (region_idx in range(length(finemapping_regions_tsv_with_header) - 1)) {
    Int region_idx_plus_one = region_idx + 1
    region bounds = {
      'chrom': finemapping_regions_tsv_with_header[region_idx_plus_one][0],
      'start': finemapping_regions_tsv_with_header[region_idx_plus_one][1],
      'end': finemapping_regions_tsv_with_header[region_idx_plus_one][2],
    }

    Int region_chrom_zero_indexed = bounds.chrom - 1
    call tasks.plink_ld { input :
      input_pfiles = pfiles_[region_chrom_zero_indexed],
      bounds = bounds,

      sample_file = get_training_samples.training_samples,
      prefix = '~{bounds.chrom}_~{bounds.start}_~{bounds.end}'
    }
  }

  output {
    File pheno_and_covars = prep_ukb_inputs.unrelated_pheno_df
    Array[File] gwas_stdout_logs = chromosomal_plink_gwas.stdout_log
    Array[File] gwas_stderr_logs = chromosomal_plink_gwas.stderr_log
    File gwas = plink_gwas.tsv
    File data = generate_finemapping_regions.data
    Array[File] lds = plink_ld.ld
    Array[File] ld_logs = plink_ld.log
    Array[File] ld_freqs = plink_ld.freq
  }
}
