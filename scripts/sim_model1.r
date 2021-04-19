source("dynamics_model1.r")

resfile <- commandArgs(T)

parameters <- expand.grid(
    n_id = c(1000),
    scenario = 1,
    time_start = 0,
    time_end = 200,
    time_steps = 1,
    analysis_timepoint = 200,
    per_timepoint = c(5, 100),
    per_effect = c(-0.1, 0.2),
    disease_threshold = 0.45,
    a = 0.18,
    b = 0.27,
    c = 0.0747,
    e = 0.18,
    f = 0.27,
    g = 0.0585,
    pval_threshold = 1e-4,
    sim=1:1000
)

starting_condition_parameters <- tibble(
    scenario = 1,
    n_phen = 3,
    variables = c("X", "K", "P"),
    phen_mean = c(0.74,1.2,1),
    phen_var = c(0.08,0.08,0.08),
    phen_h2 = c(0.2, 0.3, 0.4),
    genotype_af = c(0.2, 0.3, 0.4)
)

l <- list()
for(i in 1:nrow(parameters))
{
    l[[i]] <- simulation(parameters[i,], starting_condition_parameters)$res
}

result <- bind_rows(l)

save(result, file=resfile)

