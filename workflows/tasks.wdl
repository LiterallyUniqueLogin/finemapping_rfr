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
#
#task susie {
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

task ldstore {

  # can't get this to run - it won't match the variants I'm writing the z file to those
  # its loading from the bgen but doesn't say why.

  input {
    bgen input_bgen
    File mfi
    region bounds # need to make this optional
    Int n_samples

    File? sample_file # one sample ID per line, no header?
  }

  output {
    File ld = "ld"
  }

  command <<<
    # write z file
    {
      printf "rsid chromosome position allele1 allele2\n" ;
      awk '{ if (($3 >= ~{bounds.start}) && ($3 <= ~{bounds.end})) { print $2 " " ~{bounds.chrom} " " $3 " " $4 " " $5 } }' ~{mfi}
    } > file.z

    ~{if defined(sample_file) then "ln ~{sample_file} sample" else "" }
    ~{if defined(sample_file) then "ln ~{sample_file} incl" else "" }
    echo "z;bgen;bgi;ld;n_samples~{if defined(sample_file) then ";sample;incl" else ""}
    file.z;~{input_bgen.bgen};~{input_bgen.index};ld;~{n_samples}~{if defined(sample_file) then ";sample;incl" else ""}" > master

    ldstore \
      --write-text \
      --read-only-bgen \
      --in-files master \
      --n-threads 28
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ldstore:v2.0"
    dx_timeout: "48h"
    memory: "56GB"
    cpus: 28
  }
}

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
  }

  output {
    File ld = "plink2.unphased.vcor1"
  }

  command <<<
    ~{if defined(bounds) then "printf 'chr~{select_first([bounds]).chrom}\\t~{select_first([bounds]).start}\\t~{select_first([bounds]).end + 1}\\n' > region.bed" else ""}
    ~{if defined(sample_file) then "{ printf 'FID\\tIID\\n' ; awk '{ print $1 \"\\t\" $1 }' ~{sample_file} ; } > plink.sample" else "" }

    plink2 \
      --pfile $(echo '~{input_pfiles.pgen}' | sed -e 's/\.pgen$//') \
      ~{if defined(bounds) then "--extract bed1 region.bed" else ""} \
      ~{if defined(sample_file) then "--keep plink.sample" else ""} \
      --r-unphased square ref-based \
      --memory 55000 \
      --threads 28
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/plink:v2.00a6LM_AVX2_Intel_5_Feb_2024"
    dx_timeout: "48h"
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

task regular_plink_gwas {

  input {
    pfiles input_pfiles
    region? bounds
    String phenotype_name

    File pheno_covar_file

    File? sample_file
  }

  output {
    File gwas = "plink2.~{phenotype_name}.glm.linear"
    File stdout_log = "plink.stdout"
    File stderr_log = "plink.stderr"
  }

  command <<<
    ~{if defined(bounds) then "printf 'chr~{select_first([bounds]).chrom}\\t~{select_first([bounds]).start}\\t~{select_first([bounds]).end + 1}\\n' > region.bed" else ""}
    ~{if defined(sample_file) then "{ printf 'FID\\tIID\\n' ; awk '{ print $1 \"\\t\" $1 }' ~{sample_file} ; } > plink.sample" else "" }

    plink2 \
      ~{if defined(sample_file) then "--keep plink.sample" else ""} \
      --no-psam-pheno \
      --pheno ~{pheno_covar_file} \
      --pheno-name ~{phenotype_name} \
      --pfile $(echo '~{input_pfiles.pgen}' | sed -e 's/\.pgen$//') \
      ~{if defined(bounds) then "--extract bed1 region.bed" else ""} \
      --mac 20 \
      --glm omit-ref pheno-ids hide-covar cols=+machr2,-test,-nobs \
      --memory 55000 \
      --threads 28 \
      > plink.stdout 2> plink.stderr
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/plink:v2.00a5LM_AVX2_Intel_23_Sep_2023"
    dx_timeout: "48h"
    memory: "56GB"
    cpus: 28
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
    docker: "ubuntu:jammy-20240212"
    dx_timeout: "30m"
    memory: "2GB"
  }
}
