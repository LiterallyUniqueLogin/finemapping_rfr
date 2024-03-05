version 1.0

import "../../ukbiobank/workflow/gwas_wdl/gwas_tasks.wdl"
import "../../ukbiobank/workflow/expanse_wdl/expanse_tasks.wdl"
import "../../ukbiobank/workflow/expanse_wdl/expanse_files.wdl"

task load_continuous_phenotype {
  input {
    String script_dir
    File script = "~{script_dir}/load_init_assessment_phenotype.py"

    String outname = "pheno"
    File sc
    File samples 
    String phenotype_name
    File fam_file
    File pcs_fname
    Int n_pcs
    File assessment_ages

    Array[String] categorical_covariate_names
    Array[File] categorical_covariate_scs
  }

  output {
    File df = "~{outname}.tab"
  }

  command <<<
    ~{script} \
      ~{outname} \
      ~{samples} \
      ~{sc} \
      ~{phenotype_name} \
      ~{fam_file} \
      ~{pcs_fname} \
      ~{n_pcs} \
      ~{assessment_ages} \
      --categorical-covar-names ~{sep=" " categorical_covariate_names} \
      --categorical-covar-files ~{sep=" " categorical_covariate_scs}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/python_data:v1.0"
    shortTask: true
    dx_timeout: "1h"
    memory: "2GB"
  }
}

task subset_df_to_sample_list {
  input {
    String outname = "pheno"
    File in_df # first column must be id
    File samples 
  }

  output {
    File df = "~{outname}.tab"
  }

  command <<<
    python -c '
    import polars as pl
    pl.read_csv("~{in_df}", separator="\t", new_columns=["id"]).join(
      pl.read_csv("~{samples}", separator="\t", new_columns=["id"]).select("id"),
      how="inner",
      on="id"
    ).write_csv(
      "~{outname}.tab",
      separator="\t"
    )
    '
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/python_data:v1.0"
    shortTask: true
    dx_timeout: "1h"
    memory: "2GB"
  }
}

task create_numpy_from_df {
  input {
    File df # first column must be id
  }

  output {
    File npy = "~{basename(df, '.tab')}.npy"
  }

  command <<<
    python -c '
    import numpy as np
    np.save(
      "~{basename(df, '.tab')}.npy",
      np.genfromtxt(
          "~{df}",
          delimiter="\t",
          skip_header = 1,
      )
    )
    '
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/python_data:v1.0"
    shortTask: true
    dx_timeout: "1h"
    memory: "2GB"
  }
}

task cut_samples_from_first_column {
  input {
    File df
  }

  output {
    File samples = "~{basename(df, '.tab')}.samples"
  }

  command <<<
    cut -f 1 -d $'\t' ~{df} > ~{basename(df, '.tab')}.samples
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/python_data:v1.0"
    shortTask: true
    dx_timeout: "1h"
    memory: "2GB"
  }
}

# --------------- Begin workflow ---------------------

workflow prep_samples_and_phenotype {

  input {
    String project_ukb_script_dir
    String old_ukb_script_dir

    String phenotype_name
    Boolean is_binary

    # should define phenotype_id or premade_pheno_df
    Int? phenotype_id
    Array[String] categorical_covariate_names = []
    Array[Int] categorical_covariate_ids = []

    File? premade_pheno_df # first col ID, second pheno, remaining are covars

    Boolean transform = true

    # If specified, must contain all samples of all ethnicities that you want included
    # (so any samples not included will be omitted)
    # samples that fail QC will still be removed
    # analyses will still be split per ethnicity
    # each ethnicity's sample list will still be shrunk to remove related participants
    File? subpop_sample_list

    Int n_pcs = 8
  }

  call expanse_files.files

  call expanse_tasks.extract_field as white_brits { input:
    script_dir = old_ukb_script_dir,
    id = 22006
  }

  call expanse_tasks.extract_field as ethnicity_self_report { input :
    script_dir = old_ukb_script_dir,
    id = 21000
  }

  call expanse_tasks.extract_field as sex_aneuploidy { input:
    script_dir = old_ukb_script_dir,
    id = 22019
  }

  call expanse_tasks.extract_field as genetic_sex { input:
    script_dir = old_ukb_script_dir,
    id = 22001
  }

  call expanse_tasks.extract_field as reported_sex { input:
    script_dir = old_ukb_script_dir,
    id = 31
  }

  call expanse_tasks.extract_field as kinship_count { input:
    script_dir = old_ukb_script_dir,
    id = 22021
  }

  call expanse_tasks.extract_field as assessment_ages { input :
    script_dir = old_ukb_script_dir,
    id = 21003
  }

  call expanse_tasks.extract_field as pcs { input :
    script_dir = old_ukb_script_dir,
    id = 22009
  }

  if (defined(phenotype_id)) {
    call expanse_tasks.extract_field as phenotype { input :
      script_dir = old_ukb_script_dir,
      id = select_first([phenotype_id])
    }

    scatter (categorical_covariate_id in categorical_covariate_ids) {
      call expanse_tasks.extract_field as categorical_covariates { input :
        script_dir = old_ukb_script_dir,
        id = categorical_covariate_id
      }
    }
  }

  call gwas_tasks.write_sample_list as white_brits_sample_list { input:
    script_dir = old_ukb_script_dir,
    sc = white_brits.data
  }

  call gwas_tasks.write_sample_list as sex_aneuploidy_sample_list { input:
    script_dir = old_ukb_script_dir,
    sc = sex_aneuploidy.data
  }

  call gwas_tasks.sex_mismatch_sample_list { input:
    script_dir = old_ukb_script_dir,
    sc_genetic_sex = genetic_sex.data,
    sc_reported_sex = reported_sex.data
  }
  
  call gwas_tasks.write_sample_list as low_genotyping_quality_sample_list { input:
    script_dir = old_ukb_script_dir,
    sc = kinship_count.data,
    value = -1
  }

  call gwas_tasks.qced_sample_list as qced_sample_list { input:
    script_dir = old_ukb_script_dir,
    unqced_sample_list = white_brits_sample_list.data,
    withdrawn_sample_list = files.withdrawn_sample_list,
    sex_aneuploidy_sample_list = sex_aneuploidy_sample_list.data,
    sex_mismatch_sample_list = sex_mismatch_sample_list.data,
    low_genotyping_quality_sample_list = low_genotyping_quality_sample_list.data,
    subpop_sample_list = subpop_sample_list
  }

  if (defined(phenotype_id)) {
    call load_continuous_phenotype { input:
      script_dir = project_ukb_script_dir,
      outname = "qced_white_brits_~{phenotype_name}",
      sc = select_first([phenotype.data]),
      samples = qced_sample_list.data,
      phenotype_name = phenotype_name,
      fam_file = files.fam_file,
      pcs_fname = pcs.data,
      n_pcs = n_pcs,
      assessment_ages = assessment_ages.data,
      categorical_covariate_names = categorical_covariate_names,
      categorical_covariate_scs = select_first([categorical_covariates.data])
    }
  }
  if (!defined(phenotype_id)) {
    call subset_df_to_sample_list as subsetted_pheno_related { input :
      outname = "qced_white_brits_~{phenotype_name}",
      in_df = select_first([premade_pheno_df]),
      samples = qced_sample_list.data
    }
  }
  
  # regardless of continuous or binary, get the outputs and move on
  File pheno_df_ = select_first([subsetted_pheno_related.df, load_continuous_phenotype.df])

  call cut_samples_from_first_column { input :
    df = pheno_df_
  }

  if (!is_binary) {
    call gwas_tasks.unrelated_samples as not_binary_phenotype_unrelated_samples { input:
      script_dir = old_ukb_script_dir,
      PRIMUS_command = 'run_PRIMUS.pl',
      kinship = files.kinship,
      sample_list = cut_samples_from_first_column.samples,
      prefix = "unrelated_qced_white_brits_~{phenotype_name}"
    }
  }
  if (is_binary) {
    call create_numpy_from_df { input :
      df = pheno_df_
    }

    call gwas_tasks.unrelated_samples as binary_phenotype_unrelated_samples { input:
      script_dir = old_ukb_script_dir,
      PRIMUS_command = 'run_PRIMUS.pl',
      kinship = files.kinship,
      sample_list = cut_samples_from_first_column.samples,
      binary_pheno_data = create_numpy_from_df.npy,
      prefix = "unrelated_qced_white_brits_~{phenotype_name}"
    }
  }

  File samples_for_phenotype_ = select_first([
    not_binary_phenotype_unrelated_samples.data,
    binary_phenotype_unrelated_samples.data
  ])

  call subset_df_to_sample_list as subsetted_pheno_unrelated { input :
    in_df = pheno_df_,
    samples = samples_for_phenotype_,
    outname = "unrelated_qced_white_brits_~{phenotype_name}"
  }

  if (transform && !is_binary) {
    call create_numpy_from_df as unrelated_npy { input :
      df = subsetted_pheno_unrelated.df
    }

    call gwas_tasks.transform_trait_values { input :
      script_dir = old_ukb_script_dir,
      pheno_data = unrelated_npy.npy,
      prefix = "rin_unrelated_qced_white_brits_~{phenotype_name}"
    }
  }

  output {
    File related_pheno_df = pheno_df_ # including related samples
    File unrelated_pheno_df = subsetted_pheno_unrelated.df # unrelated

    # sample lists subset to those with the specified phenotype
    File samples_for_phenotype = samples_for_phenotype_ # unrelated, qced and takes into account subpop if specified

    File? transformed_trait_values_npy = transform_trait_values.data
  }
}
