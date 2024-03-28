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

task convert_gwas_results_to_finemap_input {
  input {
    region bounds
    File gwas_results

    Boolean is_binary
  }

  output {
    File finemap_input = "finemap_input.z"
  }

  command <<<
    python -c '
    import ast
    import polars as pl
    pl.scan_csv(
        "~{gwas_results}", separator="\t", null_values="NA"
    ).filter(
        (pl.col("#CHROM") == ~{bounds.chrom}) &
        (pl.col("POS") >= ~{bounds.start}) &
        (pl.col("POS") <= ~{bounds.end}) &
        (pl.col("ERRCODE") == ".")
    ).select([
        ("SNP_" + pl.col("#CHROM").cast(str) + "_" + pl.col("POS").cast(str) + "_" + pl.col("REF") + "_" + pl.col("ALT")).alias("rsid"),
        ("0" + pl.col("#CHROM").cast(str)).str.slice(-2).alias("chromosome"),
        pl.col("POS").alias("position"),
        pl.col("REF").alias("allele1"),
        pl.col("ALT").alias("allele2"),
        pl.lit("nan").alias("maf"),
        (pl.col("BETA") if not ast.literal_eval("~{is_binary}".capitalize()) else pl.col("OR").log()).alias("beta"),
        (pl.col("SE") if not ast.literal_eval("~{is_binary}".capitalize()) else pl.col("LOG(OR)_SE")).alias("se"),
    ]).collect().write_csv("finemap_input.z", separator=" ")
    '
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/python_data:v1.0"
    dx_timeout: "30m"
    memory: "6GB"
  }
}

task finemap {

  input {
    Int n_samples
    File input_z 
    File ld

    Int max_causal_snps
    String prefix = ""
  }

  output {
    File snp_file = "finemapping/~{prefix}finemap_output.snp"
    File log_sss = "finemapping/~{prefix}finemap_output.log_sss"
    File config = "finemapping/~{prefix}finemap_output.config"
    Array[File] creds = glob("finemapping/~{prefix}finemap_output.cred*")
    File finemap_input_z = "finemapping/~{prefix}finemap_input.z"
  }


  command <<<
    ######## write/link inputs
    {
      echo 'z;ld;snp;config;cred;log;n_samples'
      echo 'finemap_input.z;all_variants.ld;finemap_output.snp;finemap_output.config;finemap_output.cred;finemap_output.log;~{n_samples}'
    } > finemap_input.master

    ln ~{input_z} finemap_input.z
    ln ~{ld} all_variants.ld

    ###### run
    finemap \
      --sss \
      --in-files finemap_input.master \
      --log \
      --n-threads 2 \
      --n-causal-snps ~{max_causal_snps}

    mkdir finemapping
    for file in finemap_input.z finemap_output.cred* finemap_output.config finemap_output.log_sss finemap_output.snp ; do
      ln $file finemapping/~{prefix}$file
    done
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/finemap:v1.4.2"
    dx_timeout: "8h"
    memory: "4GB"
    cpus: 2
  }
}

task convert_gwas_results_to_susie_input {
  input {
    region bounds
    File gwas_results

    Boolean is_binary
  }

  output {
    File vars = 'vars.txt'
    File effect_sizes = "effect_sizes.txt"
    File effect_standard_errors = "effect_standard_errors.txt"
  }

  command <<<
    python -c '   
    import ast
    import polars as pl
    snps = pl.scan_csv(
        "~{gwas_results}", separator="\t", null_values="NA"
    ).filter(
        (pl.col("#CHROM") == ~{bounds.chrom}) &
        (pl.col("POS") >= ~{bounds.start}) &
        (pl.col("POS") <= ~{bounds.end}) &
        (pl.col("ERRCODE") == ".")
    ).select([
        ("SNP_" + pl.col("#CHROM").cast(str) + "_" + pl.col("POS").cast(str) + "_" + pl.col("REF") + "_" + pl.col("ALT")).alias("id"),
        (pl.col("BETA") if not ast.literal_eval("~{is_binary}".capitalize()) else pl.col("OR").log()).alias("beta"),
        (pl.col("SE") if not ast.literal_eval("~{is_binary}".capitalize()) else pl.col("LOG(OR)_SE")).alias("se"),
    ]).collect()
    snps[["id"]].write_csv("vars.txt", include_header=False)
    snps[["beta"]].write_csv("effect_sizes.txt", include_header=False)
    snps[["se"]].write_csv("effect_standard_errors.txt", include_header=False)
    '
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/python_data:v1.0"
    dx_timeout: "30m"
    memory: "6GB"
  }
}

task susie {

  input {
    String script_dir
    File script = "~{script_dir}/sum_stats_susie.r"

    File vars
    File effect_sizes
    File effect_standard_errors
    File correlation_matrix

    Int n_samples
    Int L

    String prefix = ''
  }

  output {
   File lbf = "finemapping/~{prefix}lbf.tab"
   File lbf_variable = "finemapping/~{prefix}lbf_variable.tab"
   File sigma2 = "finemapping/~{prefix}sigma2.txt"
   File V = "finemapping/~{prefix}V.tab"
   File converged = "finemapping/~{prefix}converged.txt"
   File lfsr = "finemapping/~{prefix}lfsr.tab"
   File requested_coverage = "finemapping/~{prefix}requested_coverage.txt"
   File alpha = "finemapping/~{prefix}alpha.tab"
   File vars =  "finemapping/~{prefix}vars.txt"
   Array[File] CSs = glob("finemapping/~{prefix}cs*.txt")
  }

  command <<<
    Rscript ~{script} \
      "~{prefix}" \
      ~{effect_sizes} \
      ~{effect_standard_errors} \
      ~{correlation_matrix} \
      ~{n_samples} \
      ~{L}

    mkdir finemapping
    for file in "~{prefix}"*{lbf.tab,lbf_variable.tab,sigma2.txt,V.tab,converged.txt,lfsr.tab,requested_coverage.txt,alpha.tab,cs*.txt} ; do
      ln $file finemapping/$file
    done
    ln ~{vars} finemapping/~{prefix}vars.txt
  >>>


  runtime {
    docker: "quay.io/thedevilinthedetails/work/susie:v0.12.35"
    dx_timeout: '4h'
    memory: '10GB'
  }
}

task plink_chromosomal_ld_many_regions {

  input {
    pfiles input_pfiles
    Int chrom
    Array[Int] starts
    Array[Int] ends

    File? sample_file
    String prefix = 'plink2'

    Array[String] out_ld_files_in_order
    Array[String] out_vars_files_in_order
    Array[String] out_log_files_in_order
  }

  output {
    Array[File] ld = out_ld_files_in_order #glob("*.unphased.vcor1")
    Array[File] vars = out_vars_files_in_order #glob("*.unphased.vcor1.vars")
    Array[File] log = out_log_files_in_order #glob("*.log")
  }

  command <<<
    ~{if defined(sample_file) then "{ printf 'FID\\tIID\\n' ; awk '{ print $1 \"\\t\" $1 }' ~{sample_file} ; } > plink.sample" else "" }

    starts=(~{sep=" " starts})
    ends=(~{sep=" " ends})
    for (( region_idx=0; region_idx<${#starts[@]}; region_idx++ )); do
      full_prefix="~{prefix}_~{chrom}_${starts[region_idx]}_${ends[region_idx]}"
      printf "~{chrom}\t${starts[region_idx]}\t$(( ${ends[region_idx]} + 1 ))\n" > "$full_prefix"_region.bed

      plink2 \
        --pfile $(echo '~{input_pfiles.pgen}' | sed -e 's/\.pgen$//') \
        --extract bed1 "$full_prefix"_region.bed \
        --maf 0.0005 \
        ~{if defined(sample_file) then "--keep plink.sample" else ""} \
        --r-unphased square ref-based \
        --memory 55000 \
        --threads 28 \
        --out "$full_prefix"
      
      mv "$full_prefix".unphased.vcor1 "$full_prefix".unphased.vcor1.tab
      sed -e 's/\t/ /g' "$full_prefix".unphased.vcor1.tab > "$full_prefix".unphased.vcor1
    done
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/plink:v2.00a6LM_AVX2_Intel_5_Feb_2024"
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

    mv "$prefix".unphased.vcor1 "$prefix".unphased.vcor1.tab
    sed -e 's/\t/ /g' "$prefix".unphased.vcor1.tab > "$prefix".unphased.vcor1
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
