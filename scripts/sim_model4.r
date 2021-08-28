library(parallel)
<<<<<<< HEAD
source("dynamics_model4.r")

resfile <- commandArgs(T)

consts <- expand.grid(
    n_id = c(1000),
    scenario = 1,
    time_start = 0,
    time_end = 100,
    time_steps = 1,
    parm_n = 2,
    per_timepoint = c(5, 100),
    per_effect = c(-0.1, 0.2),
    disease_threshold = 0.5,
=======
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
>>>>>>> 4eb3a3ae05592178e7000bfb41e57d78ff150528
    pval_threshold = 1e-4,
    sim=1:100
)

<<<<<<< HEAD
cond <- tibble(
    scenario = 1,
    n_phen = 4,
    variables = c("p", "a", "m", "f"),
    phen_mean = c(1,7.5,3,3),
    phen_var = c(1,7.5,3,3),
    phen_h2 = c(0.8, 0.8, 0.8, 0.8),
    genotype_af = c(0.4, 0.4, 0.4, 0.4)
)

consts %>% unlist() %>% sum() %>% set.seed()


cond_p <- tibble(
    scenario = 1,
    n_phen = 14,
    variables = c("Pbp", "Ppp", "Pfp", "App", "Aph", "Afp", "Afh", "Mbp", "Mpp", "Mph", "Fdam", "gp", "gm", "gf"),
    par_mean = c(0.25,200,25,75,1.5,80,2.5,0.05,25,5,0.05,5,5,2.5),
    par_var = c(0.25,200,25,75,1.5,80,2.5,0.05,25,5,0.05,5,5,2.5),
    par_h2 = c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8),
    genotype_af = c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4)
)


l <- mclapply(1:nrow(consts), function(i)
{
    simulation(cond_p, cond, consts[i,])$res
=======
starting_condition_parameters <- tibble(
    scenario = 1,
    n_phen = 3,
    variables = c("X", "K", "P"),
    phen_mean = c(0.74,1.2,1),
    phen_var = c(0.08,0.08,0.08),
    phen_h2 = c(0.4, 0.4, 0.4),
    genotype_af = c(0.2, 0.3, 0.4)
)

l <- mclapply(1:nrow(parameters), function(i)
{
    simulation(parameters[i,], starting_condition_parameters)$res
>>>>>>> 4eb3a3ae05592178e7000bfb41e57d78ff150528
}, mc.cores=16)

result <- bind_rows(l)

save(result, file=resfile)

