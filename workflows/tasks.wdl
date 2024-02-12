version 1.0

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

task z_file_for_ldstore {
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
    File z_file
    Int n_samples

    File? sample_file # one sample ID per line, no header?
  }

  output {
    File ld = "ld"
  }

  command <<<
    echo "z;bgen;bgi;ld;n_samples~{if defined(sample_file) then ";incl" else ""}
    ~{z_file};~{input_bgen.bgen};{input_bgen}.index;ld;~{n_samples}~{";" + sample_file}" > master

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
    File range

    File? sample_file
  }

  output {
    
  }

  command <<<
    plink2 \
      --pfile $(echo '~{input_pfiles.pgen}' | sed -e 's/\.pgen$//') \
      ~{"--keep " + sample_file} \
      --extract bed1 
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

  }

  output {

  }

  command <<<

  >>>

  runtime {
    docker: "quay.io/thedevilinthedetails/work/plink:v2.00a5LM_AVX2_Intel_23_Sep_2023"
    dx_timeout:
    memory:
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

