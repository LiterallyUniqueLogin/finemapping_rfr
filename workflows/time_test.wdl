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

  scatter (chrom_zero_indexed in [0]) { # range(22)) {
    Int chrom_one_indexed = chrom_zero_indexed + 1

    pfiles pfiles_ = {
      'pgen': '../ukbiobank/array_imputed/pfile_converted/chr~{chrom_one_indexed}.pgen',
      'pvar': '../ukbiobank/array_imputed/pfile_converted/chr~{chrom_one_indexed}.pvar',
      'psam': '../ukbiobank/array_imputed/pfile_converted/chr~{chrom_one_indexed}.psam',
    }

    call tasks.plink_gwas as chromosomal_plink_gwas { input :
      input_pfiles = pfiles_,
      phenotype_name = 'platelet_count',
      pheno_covar_file = prep_ukb_inputs.unrelated_pheno_df,
      quantile_normalize = true,
      sample_file = get_training_samples.training_samples,
      prefix = 'platelet_count_training_gwas_chr~{chrom_one_indexed}'
    }

    call gwas_tasks.generate_finemapping_regions { input :
      script_dir = '../ukbiobank',
      chr_lens = '../ukbiobank/misc_data/genome/chr_lens.txt',
      phenotype = 'platelet_count',
      snp_assoc_results = chromosomal_plink_gwas.gwas,
      prefix = 'platelet_count_training_chr~{chrom_one_indexed}_'
    }
    Array[Array[String]] finemapping_regions_tsv_with_header = read_tsv(generate_finemapping_regions.data)

    scatter (region_idx in range(length(finemapping_regions_tsv_with_header) - 1)) {
      Int region_idx_plus_one = region_idx + 1
      Int finemapping_region_starts = finemapping_regions_tsv_with_header[region_idx_plus_one][1]
      Int finemapping_region_ends = finemapping_regions_tsv_with_header[region_idx_plus_one][2]
      String out_ld_files_in_order = 'platelet_count_training_ld_~{chrom_one_indexed}_~{finemapping_region_starts}_~{finemapping_region_ends}.unphased.vcor1'
      String out_vars_files_in_order = 'platelet_count_training_ld_~{chrom_one_indexed}_~{finemapping_region_starts}_~{finemapping_region_ends}.unphased.vcor1.vars'
      String out_log_files_in_order = 'platelet_count_training_ld_~{chrom_one_indexed}_~{finemapping_region_starts}_~{finemapping_region_ends}.log'
    }

    call tasks.plink_chromosomal_ld_many_regions { input :
      input_pfiles = pfiles_,
      chrom = chrom_one_indexed,
      starts = finemapping_region_starts,
      ends = finemapping_region_ends,

      sample_file = get_training_samples.training_samples,
      prefix = 'platelet_count_training_ld',

      out_ld_files_in_order = out_ld_files_in_order,
      out_vars_files_in_order = out_vars_files_in_order,
      out_log_files_in_order = out_log_files_in_order
    }

    scatter (region_idx_again in range(length(finemapping_regions_tsv_with_header) - 1)) {
      Int region_idx_again_plus_one = region_idx_again + 1
      region bounds = {
        "chrom": chrom_one_indexed,
        "start": finemapping_regions_tsv_with_header[region_idx_again_plus_one][1],
        "end": finemapping_regions_tsv_with_header[region_idx_again_plus_one][2],
      }

      call tasks.convert_gwas_results_to_finemap_input { input :
        gwas_results = chromosomal_plink_gwas.gwas,
        bounds = bounds,
        is_binary = false 
      }

      call tasks.finemap { input :
        n_samples = 50000,
        input_z = convert_gwas_results_to_finemap_input.finemap_input,
        ld = plink_chromosomal_ld_many_regions.ld[region_idx_again],
        max_causal_snps = 30,
        prefix = "platelet_count_training_~{bounds.chrom}_~{bounds.start}_~{bounds.end}_"
      }
      # TODO check no max causal snps, check converged

      call tasks.convert_gwas_results_to_susie_input { input :
        bounds = bounds,
        gwas_results = chromosomal_plink_gwas.gwas,
        is_binary = false
      }

      call tasks.susie { input :
        script_dir = "scripts",
        vars = convert_gwas_results_to_susie_input.vars,
        effect_sizes = convert_gwas_results_to_susie_input.effect_sizes,
        effect_standard_errors = convert_gwas_results_to_susie_input.effect_standard_errors,
        correlation_matrix = plink_chromosomal_ld_many_regions.ld[region_idx_again],
        n_samples = 50000,
        L = 30,
        prefix = "platelet_count_training_~{bounds.chrom}_~{bounds.start}_~{bounds.end}_"
      }
      # TODO check no max causal snps, check converged
    }
  }

  call tasks.concatenate_tsvs as plink_gwas { input :
    tsvs = chromosomal_plink_gwas.gwas,
    prefix = 'platelet_count_training_gwas'
  }

  call tasks.concatenate_tsvs as finemapping_regions { input :
    tsvs = generate_finemapping_regions.data,
    prefix = "platelet_count_training_finemapping_regions"
  }

  output {
    File pheno_and_covars = prep_ukb_inputs.unrelated_pheno_df

    Array[File] gwas_stdout_logs = chromosomal_plink_gwas.stdout_log
    Array[File] gwas_stderr_logs = chromosomal_plink_gwas.stderr_log
    File gwas = plink_gwas.tsv

    File data = finemapping_regions.tsv

    Array[File] ld_logs = flatten(plink_chromosomal_ld_many_regions.log)

    Array[File] finemap_snp_file = flatten(finemap.snp_file)
    Array[File] finemap_log_sss = flatten(finemap.log_sss)
    Array[File] finemap_config = flatten(finemap.config)
    Array[Array[File]] finemap_creds = flatten(finemap.creds)
    Array[File] finemap_input_z = flatten(finemap.finemap_input_z)

    Array[File] susie_lbf = flatten(susie.lbf)
    Array[File] susie_lbf_variable = flatten(susie.lbf_variable)
    Array[File] susie_sigma2 = flatten(susie.sigma2)
    Array[File] susie_V = flatten(susie.V)
    Array[File] susie_converged = flatten(susie.converged)
    Array[File] susie_lfsr = flatten(susie.lfsr)
    Array[File] susie_requested_coverage = flatten(susie.requested_coverage)
    Array[File] susie_alpha = flatten(susie.alpha)
    Array[File] susie_vars = flatten(susie.vars)
    Array[Array[File]] susie_CSs = flatten(susie.CSs)
  }
}
