#' Simulate genotypes influencing a set of phenotypes
#'
#' One genotype has an influence on one phenotype. A counfounder has a random effect on all phenotypes. 
#'
#' @param n_phen How many phenotypes
#' @param n_id How many individuals
#' @param phen_mean The means for each phenotype
#' @param phen_var The variances for each phenotype
#' @param phen_h2 The heritabilities for each phenotype
#' @param genotype_af The allele frequencies of each genotype
#'
#' @export
#' @return List, where geno is the matrix of genotypes and phen is a data frame of the phenotypes
simulate_starting_conditions <- function(n_phen, n_id, phen_mean, phen_var, phen_h2, genotype_af) 
{
	stopifnot(length(phen_mean) == n_phen)
	stopifnot(length(phen_var) == n_phen)
	stopifnot(length(phen_h2) == n_phen)
	stopifnot(length(genotype_af) == n_phen)

	confounder <- rnorm(n_id)
	confounder_var <- runif(n_phen, 0, (1 - phen_h2) * phen_var)
	residual_var <- (phen_var - phen_h2 * phen_var - confounder_var)

	colSums(rbind(residual_var, confounder_var, phen_h2 * phen_var))

	geno <- sapply(genotype_af, function(x) rbinom(n_id, 2, x))
	b_gx <- sqrt(phen_h2 * phen_var / (2 * genotype_af * (1- genotype_af)))

	l <- list()
	for(i in 1:n_phen)
	{
		l[[i]] <- phen_mean[i] + (geno[,i] - mean(geno[,i])) * b_gx[i] + confounder * sqrt(confounder_var[i]) + rnorm(n_id, 0, sqrt(residual_var[i]))
	}

	sapply(l, mean)

	return(list(
		geno=geno,
		phen=do.call(cbind, l) %>% as_tibble()
	))
}


# example
set.seed(1234)

n_phen = 3
n_id = 1000
phen_mean = c(1,2,3)
phen_var = c(10,20,30)
phen_h2 = c(0.2, 0.3, 0.4)
genotype_af = c(0.2, 0.3, 0.4)
dat <- simulate_starting_conditions(n_phen=n_phen, n_id=n_id, phen_mean=phen_mean, phen_var=phen_var, phen_h2=phen_h2, genotype_af=genotype_af)


# check
colMeans(dat$geno)/2
genotype_af

colMeans(dat$phen)
phen_mean

apply(dat$phen, 2, var)
phen_var

