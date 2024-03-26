version 1.0

struct region {
  Int chrom
  Int start
  Int end
}

struct bgen {
  File bgen
  File index
}

struct pfiles {
  File pgen
  File psam
  File pvar
}

## helper task

task concatenate_tsvs {
  input {
    # this can be larger than the max number of arguments for a command on the command line,
		# but cannot be too long for a string variable or an array in bash
    Array[File]+ tsvs
    String prefix = "out"
  }

  output {
    File tsv = "~{prefix}.tab"
  }

  command <<<
    bash -c '
      tsvs=(~{sep=" " tsvs})
      if (( $(for (( i=0; i<${#tsvs[@]}; i++ )); do
        head -1 "${tsvs[$i]}"
      done | uniq -c | wc -l) != 1 )) ; then
        echo "Different headers" 2>&1
        exit 1
      fi

      head -1 ~{tsvs[0]} > ~{prefix}.tab
      for (( i=0; i<${#tsvs[@]}; i++ )); do
        tail -n +2 "${tsvs[$i]}" >> ~{prefix}.tab
      done
    '
  >>>

  runtime {
    docker: "ubuntu:jammy-20240212"
    dx_timeout: "4h"
    memory: "2GB"
  }
}


# sample files described in tasks in this file are
# one sample per line, no header

#task finemap {
#
#  input {
#    
#  }
#
#  output {
#
#  }
#
#  command <<<
#
#  >>>
#
#  runtime {
#    docker: "quay.io/thedevilinthedetails/work/finemap:v1.4.2"
#    dx_timeout:
#    memory:
#  }
#}

task susie {

  input {
    String script_dir
    File script = "~{script_dir}/sum_stats_susie.r"

    File effect_sizes
    File effect_standard_errors
    File correlation_matrix
    Int n_samples
    Float phenotype_variance
    Int L

    String prefix = ''
  }

  command <<<
    Rscript ~{script} \
      ~{effect_sizes} \
      ~{effect_standard_errors} \
      ~{correlation_matrix} \
      ~{n_samples} \
      ~{phenotype_variance} \
      ~{L}
  >>>

  output {
    File lbf = "~{prefix}lbf.tab"
    File lbf_variable = "~{prefix}lbf_variable.tab"
    File sigma2 = "~{prefix}sigma2.txt"
    File V = "~{prefix}V.tab"
    File converged = "~{prefix}converged.txt"
    File lfsr = "~{prefix}lfsr.tab"
    File requested_coverage = "~{prefix}requested_coverage.txt"
    File alpha = "~{prefix}alpha.tab"
    File colnames = "~{prefix}colnames.txt"
    Array[File] CSs = glob("~{prefix}cs*.txt")
  }

  runtime {
    docker: "quay.io/thedevilinthedetails/work/susie:v0.12.35"
    dx_timeout: '4h'
    memory: '10GB'
  }
}

#task plink_ld_many_regions {
#  input {
#    pfiles input_pfiles
#    Array[region] bounds
#
#    File? sample_file
#  }
#
#  output {
#    File ld = "plink2.unphased.vcor1"
#    File log = "plink2.log"
#  }
#
#  command <<<
#    ~{if defined(bounds) then "printf '~{select_first([bounds]).chrom}\\t~{select_first([bounds]).start}\\t~{select_first([bounds]).end + 1}\\n' > region.bed" else ""}
#    ~{if defined(sample_file) then "{ printf 'FID\\tIID\\n' ; awk '{ print $1 \"\\t\" $1 }' ~{sample_file} ; } > plink.sample" else "" }
#
#    plink2 \
#      --pfile $(echo '~{input_pfiles.pgen}' | sed -e 's/\.pgen$//') \
#      ~{if defined(bounds) then "--extract bed1 region.bed" else ""} \
#      ~{if defined(sample_file) then "--keep plink.sample" else ""} \
#      --r-unphased square ref-based \
#      --memory 55000 \
#      --threads 28
#  >>>
#
#  runtime {
#    docker: "quay.io/thedevilinthedetails/work/plink:v2.00a6LM_AVX2_Intel_5_Feb_2024"
#    # with 50k samples takes ~12sec for 1k variants, ~15m for 10k variants
#    # with 5k samples takes ~2sec for 1k variants, ~2m8s for 10k variants
#    dx_timeout: "48h"
#    memory: "56GB"
#    cpus: 28
#  }
#
#}

task plink_ld {
  # confirmed that this returns correct values: 2024/03/01
  # using bgen_reader to read instead and then np.corrcoef produces 0.00356942
  # for chr1 first and second vars instead of 0.00358023 from plink
  # so there's some precision difference one of these somewhere,
  # but no such difference for first and third vars with corr 0.227007

  input {
    pfiles input_pfiles
    region? bounds

    File? sample_file
    String prefix = 'plink2'
  }

  output {
    File ld = "~{prefix}.unphased.vcor1"
    File freq = "~{prefix}.afreq"
    File vars = "~{prefix}.unphased.vcor1.vars"
    File log = "~{prefix}.log"
  }

  command <<<
    ~{if defined(bounds) then "printf '~{select_first([bounds]).chrom}\\t~{select_first([bounds]).start}\\t~{select_first([bounds]).end + 1}\\n' > region.bed" else ""}
    ~{if defined(sample_file) then "{ printf 'FID\\tIID\\n' ; awk '{ print $1 \"\\t\" $1 }' ~{sample_file} ; } > plink.sample" else "" }

    plink2 \
      --pfile $(echo '~{input_pfiles.pgen}' | sed -e 's/\.pgen$//') \
      ~{if defined(bounds) then "--extract bed1 region.bed" else ""} \
      --maf 0.0005 \
      ~{if defined(sample_file) then "--keep plink.sample" else ""} \
      --r-unphased square ref-based \
      --memory 55000 \
      --threads 28 \
      --freq \
      --out ~{prefix}
  >>>

  ## --maf 0.0002 culls to about 1/3 of the variants (MAC 20 for 50k individuals)
  ## --maf 0.005 culls to between 1/4-1/5 of the variants (MAC 50 for 50k individuals)
  ## --maf 0.001 culls to about 1/6 of the variants (MAC 100 for 50k individuals)


  runtime {
    docker: "quay.io/thedevilinthedetails/work/plink:v2.00a6LM_AVX2_Intel_5_Feb_2024"
    # with 50k samples takes ~12sec for 1k variants, ~15m for 10k variants - got to 2 min in a later run, unclear why
    # with 5k samples takes ~2sec for 1k variants, ~2m8s for 10k variants
    # costs $0.15 for a genome-wide run for platelet count at an MAF threshold corresponding to MAC 50
    dx_timeout: "8h"
    memory: "56GB"
    cpus: 28
  }
}

#task ld_with_covars {
#
#  input {
#
#  }
#
#  output {
#
#  }
#
#  command <<<
#
#  >>>
#
#  runtime {
#    docker: "quay.io/thedevilinthedetails/work/"
#    dx_timeout:
#    memory:
#  }
#}

task plink_gwas {

  input {
    pfiles input_pfiles
    region? bounds
    String phenotype_name

    File pheno_covar_file
    Boolean quantile_normalize

    File? sample_file

    String prefix = "plink2"
  }

  output {
    File gwas = "~{prefix}.~{phenotype_name}.glm.linear"
    File stdout_log = "~{prefix}.stdout"
    File stderr_log = "~{prefix}.stderr"
  }

  command <<<
    ~{if defined(bounds) then "printf '~{select_first([bounds]).chrom}\\t~{select_first([bounds]).start}\\t~{select_first([bounds]).end + 1}\\n' > region.bed" else ""}
    ~{if defined(sample_file) then "{ printf 'FID\\tIID\\n' ; awk '{ print $1 \"\\t\" $1 }' ~{sample_file} ; } > plink.sample" else "" }

    # replace id column with two identical columns FID and FID
    {
      printf "FID\tIID\t" ; head -1 ~{pheno_covar_file} | cut -f 2-
      paste <(cut -f 1 ~{pheno_covar_file} | tail -n +2) <(tail -n +2 ~{pheno_covar_file})
    } > plink_formatted_df.tab

    plink2 \
      ~{if defined(sample_file) then "--keep plink.sample" else ""} \
      --no-psam-pheno \
      --pheno plink_formatted_df.tab \
      --pheno-name ~{phenotype_name} \
      --covar-name $(head -1 plink_formatted_df.tab | cut -f 4- ) \
      --covar-variance-standardize \
      ~{if quantile_normalize then "--pheno-quantile-normalize" else ""} \
      --pfile $(echo '~{input_pfiles.pgen}' | sed -e 's/\.pgen$//') \
      ~{if defined(bounds) then "--extract bed1 region.bed" else ""} \
      --maf 0.0005 \
      --glm omit-ref pheno-ids hide-covar cols=+machr2,-test,-nobs \
      --memory 55000 \
      --threads 28 \
      --out ~{prefix} \
      > ~{prefix}.stdout 2> ~{prefix}.stderr
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/plink:v2.00a5LM_AVX2_Intel_23_Sep_2023"
    dx_timeout: "24h"
    memory: "56GB"
    cpus: 28
    # costs $0.79 for each run genome-wide with MAC < 20 filter and 50k samples
    # so would expect the cost to be 3/4.5 that amount with MAF equivalent to MAC < 50 filter, i.e. $0.53
  }
}

task sample_betas {

  input {
    String script_dir
    File script = "~{script_dir}/scripts/sample_betas.py"

    File original_betas # one value per line, corresponding to each variant
    File ld # whitespace delimited file of square float array written as text
    Float sigma
    Int num_samples
    Int n_replicates
  }

  output {
    File new_betas = "replicate_betas.tab" # len(original_betas) x n_replicates
  }

  command <<<
    python ~{script} ~{original_betas} ~{ld}  ~{sigma} ~{num_samples} ~{n_replicates}
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/numpy:v1.26.4"
    dx_timeout: "2h"
    memory: "10GB"
  }
}

task get_training_samples {

  input {
    File samples
    Int num_samples
  }

  output {
    File training_samples = "training.sample" 
  }

  command <<<
    head -n ~{num_samples} ~{samples} > training.sample
  >>>

  runtime {
    docker: "ubuntu:jammy-20240212"
    dx_timeout: "30m"
    memory: "2GB"
  }
}

task get_validation_subsample_replicates {
  
  input {
    File samples
    Int num_samples
    Int replicate 
  }

  output {
    File replicate_samples = "replicate_~{replicate}.sample"
  }

  command <<<
    # from https://www.gnu.org/software/coreutils/manual/html_node/Random-sources.html
    get_seeded_random()
    {
      seed="$1"
      openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
        </dev/zero 2>/dev/null
    }

    # cut off the first num_samples samples from the file which will have been used for training, then select randomly from the rest
    shuf -n ~{num_samples} \
      --random-source=<(get_seeded_random ~{replicate}) \
      <(tail -n +~{num_samples + 1} ~{samples}) \
      > replicate_~{replicate}.sample
  >>>

  runtime {
    #docker: "ubuntu:jammy-20240212" would need a docker with openssl
    dx_timeout: "30m"
    memory: "2GB"
  }
}
