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

task finemap {

  input {
    
  }

  output {

  }

  command <<<

  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/finemap:v1.4.2"
    dx_timeout:
    memory:
  }
}

task susie {

  input {

  }

  output {

  }

  command <<<

  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/"
    dx_timeout:
    memory:
  }
}

task ldstore {

  input {
    bgen input_bgen
    File mfi
    Int chrom
    Int start
    Int end
    Int n_samples

    File? sample_file # one sample ID per line, no header?
  }

  output {
    File ld = "ld"
  }

  command <<<
    # write z file
    {
      printf "rsid\tchromosome\tposition\tallele1\tallele2\n" ;
      awk 'if (($3 >= ~{start}) && ($3 <= ~{end})) { print $2 "\t" ~{chrom} "\t" $3 "\t" $4 "\t" $5 }' ~{mfi}
    } > file.z

    echo "z;bgen;bgi;ld;n_samples~{if defined(sample_file) then ";incl" else ""}
    file.z;~{input_bgen.bgen};{input_bgen}.index;ld;~{n_samples}~{";" + sample_file}" > master

    ldstore \
      --write-text \
      --read-only-bgen \
      --in-files master \
      --memory 55 \
      --n-threads 28 \
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/ldstore:v2.0"
    dx_timeout: "48h"
    memory: "56GB"
    cpus: 28
  }
}

task plink_ld {

  input {
    pfiles input_pfiles
    region? range

    File? sample_file
  }

  output {
    File ld = "plink2.unphased.vcor1"
  }

  command <<<
    ~{if is_defined(range) then "printf 'chr~{range.chrom}\t~{range.start}\t~{range.end + 1}\n' > region.bed" else ""}

    plink2 \
      --pfile $(echo '~{input_pfiles.pgen}' | sed -e 's/\.pgen$//') \
      ~{if is_defined(range) then "--extract bed1 region.bed \
      ~{"--keep " + sample_file} \
      --r-unphased square ref-based \
      --memory 55000 \
      --threads 28
  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/plink:v2.00a5LM_AVX2_Intel_23_Sep_2023"
    dx_timeout: "48h"
    memory: "56GB"
    cpus: 28
  }
}

task ld_with_covars {

  input {

  }

  output {

  }

  command <<<

  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/"
    dx_timeout:
    memory:
  }
}

task regular_plink_gwas {

  input {
    pfiles input_pfiles
    region? range
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
    ~{if is_defined(range) then "printf 'chr~{range.chrom}\t~{range.start}\t~{range.end + 1}\n' > region.bed" else ""}

    plink2 \
      ~{"--keep " + sample_file} \
      --no-psam-pheno \
      --pheno ~{pheno_covar_file} \
      --pheno-name ~{phenotype_name} \
      --pfile $(echo '~{input_pfiles.pgen}' | sed -e 's/\.pgen$//') \
      ~{if is_defined(range) then "--extract bed1 region.bed" else ""} \
      --mac 20 \
      --glm omit-ref pheno-ids hide-covar cols=-test,-nobs \ # TODO add columns
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

  }

  output {

  }

  command <<<

  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/"
    dx_timeout:
    memory:
  }
}

task ld_with_covars {

  input {

  }

  output {

  }

  command <<<

  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/"
    dx_timeout:
    memory:
  }
}

