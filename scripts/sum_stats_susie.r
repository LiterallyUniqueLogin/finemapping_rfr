library(susieR)
library(argparse)

parser = ArgumentParser()
parser$add_argument('effect_sizes')
parser$add_argument('effect_standard_errors')
parser$add_argument('correlation_matrix')
parser$add_argument('n_samples')
parser$add_argument('phenotype_variance')
parser$add_argument('L')
args = parser$parse_args()

results = susie_rss(
	R=read.table(args$correlation_matrix, sep="\t"),
	n=args$n_samples,
	bhat=read.table(args$effect_sizes),
	shat=read.table(args$effect_standard_errors),
	var_y=args$phenotype_variance,
	estimate_residual_variance = TRUE,
	check_prior = TRUE,
	min_abs_corr = 0,
	L=args$L,
	coverage = 0.9,
	max_iter = 500
)

write.table(results$alpha, paste(dir, '/alpha.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # get pips per variable by 1-apply(1-alpha, 1, prod)
write.table(results$lbf, paste(dir, '/lbf.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) #log bayes factor of each cs
write.table(results$lbf_variable, paste(dir, '/lbf_variable.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) #lbf of each variable per cs
write.table(results$sigma2, paste(dir, '/sigma2.txt', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # estimated residual posterior variance of y
write.table(results$V, paste(dir, '/V.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # estimated prior variance of each CS
write.table(results$converged, paste(dir, '/converged.txt', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # TRUE or FALSE, did SuSiE converge within the specified tolerance
write.table(susie_get_lfsr(results), paste(dir, '/lfsr.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # local false sign rate for each cs - probability that we've got the sign of the effect wrong

sink(paste(dir, '/requested_coverage.txt', sep=''))
cat(results$sets$requested_coverage, '\n')
sink()
count = 0
for (id in results$sets$cs_index)  {
  count = count+1
  cs_id = paste('L', id, sep='')
  sink(paste(dir, '/cs', id, '.txt', sep=''))
  cat(results$sets$cs[[cs_id]], sep=' ', '\n')
  cat(results$sets$coverage[count], '\n')
  cat(results$sets$purity[cs_id, 'min.abs.corr'],
      results$sets$purity[cs_id, 'mean.abs.corr'],
      results$sets$purity[cs_id, 'median.abs.corr'], '\n')
  sink()
}

