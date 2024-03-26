version 1.0

import "tasks.wdl"

workflow validation_replicates {

  scatter (replicate in range(1)) {
   
    call tasks.get_validation_subsample_replicates 

    scatter (chromosome in range(22)) {
      call gwas_regions

      scatter (gwas_region in gwas_regions.region) {
        call gwas
      }

      call collect_gwas

      call calculate_ld_matrices

      scatter (ld_matrix in calculate_ld_matrices.ld_matrices) {

      }
    }
  }
}
