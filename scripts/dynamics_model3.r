library(data.table)
library(tidyverse)
library(pbapply)
library(simulateGP)
library(TwoSampleMR)
require(deSolve)

#solving an ode

phosph <- function(Time, State, parms)
{
    with(as.list(c(State,parms)),{
        dp = (Pbp + Ppp*p^n/(1+p^n) + Pfp*f^n/(1+f^n))*1/(1+a^n) - gp*p
        da = App*p^n/(Aph^n+p^n) + Afp*f^n/(Afh^n+f^n) - a
        dm = Mbp + Mpp*p^n/(Mph^n+p^n) - gm*m
        df = m + Fdam - gf*f

        res<-list(c(dp,da,dm,df))

        return(res)
    })
}

simulate_dynamics <- function(starting_conditions, params, n = 2, disease_threshold, time_start, time_end, time_steps)
{

    out_full <- pblapply(1:nrow(starting_conditions), function(i)
    {
        tmp <- starting_conditions[i,]
        tmp1 <- params[i,]

        # Starting concentrations
        p0 = as.numeric(tmp[1])
        a0 = as.numeric(tmp[2])
        m0 = as.numeric(tmp[3])
        f0 = as.numeric(tmp[4])

        # Parameters
        Pbp = as.numeric(tmp1[1])
        Ppp = as.numeric(tmp1[2])
        Pfp = as.numeric(tmp1[3])
        App = as.numeric(tmp1[4])
        Aph = as.numeric(tmp1[5])
        Afp = as.numeric(tmp1[6])
        Afh = as.numeric(tmp1[7])
        Mbp = as.numeric(tmp1[8])
        Mpp = as.numeric(tmp1[9])
        Mph = as.numeric(tmp1[10])
        Fdam = as.numeric(tmp1[11])
        gp = as.numeric(tmp1[12])
        gm = as.numeric(tmp1[13])
        gf = as.numeric(tmp1[14])

        odesol <- ode(func=phosph,y=c(p=p0,a=a0,m=m0,f=f0), parms=c(Pbp=Pbp, Ppp=Ppp, Pfp=Pfp, App=App, Aph=Aph, Afp=Afp, Afh=Afh, Mbp=Mbp, Mpp=Mpp, Mph=Mph, Fdam=Fdam, gp=gp, gm=gm, gf=gf, n=n), times=seq(0, 100, by = 1))
        
        sst <- odesol[nrow(odesol),"m"]
        if (sst>=disease_threshold) {
            ds <- 1
        } else {
            ds <- 0
        }
        tb <- cbind(odesol,p0=p0,a0=a0,m0=m0,f0=f0,D=ds) %>%
            as_tibble() %>%
            mutate(id=i)
        return(tb)
    }) %>% bind_rows()
}



simulation <- function(params, starting_condition_parameters)
{
    params %>% unlist() %>% sum() %>% set.seed()
    cond <- subset(starting_condition_parameters, scenario==params$scenario)

    message("Generating starting conditions")
    starting_conditions <- simulate_starting_conditions(
        phen_names=cond$variables, 
        n_id=params$n_id, 
        phen_mean=cond$phen_mean, 
        phen_var=cond$phen_var, 
        phen_h2=cond$phen_h2, 
        genotype_af=cond$genotype_af
    )
    
    starting_conditions$phen <- starting_conditions$phen + abs(min(starting_conditions$phen)) + 0.05

    message("General simulation")
    dyn <- simulate_dynamics(
        starting_conditions = starting_conditions$phen,
        Pbp = params$parm_Pbp,
        Ppp = params$parm_Ppp,
        Pfp = params$parm_Pfp,
        App = params$parm_App,
        Aph = params$parm_Aph,
        Afp = params$parm_Afp,
        Afh = params$parm_Afh,
        Mbp = params$parm_Mbp,
        Mpp = params$parm_Mpp,
        Mph = params$parm_Mph,
        Fdam = params$parm_Fdam,
        gp = params$parm_gp,
        gm = params$parm_gm,
        gf = params$parm_gf,
        n = params$parm_n,
        disease_threshold = params$disease_threshold,
        time_start = params$time_start,
        time_end = params$time_end,
        time_steps = params$time_steps
    )

    message("Analysis")
    # Analyse

    # 1. What instruments are available 

    dyn1 <- subset(dyn, time == params$analysis_timepoint)
    mrres <- mr_analysis(geno=starting_conditions$geno, out = dyn1, exposures=c("p", "a", "f"), outcomes="m")
    obsres <- obs_analysis(out = dyn1, exposures=starting_condition_parameters$variables, outcomes="D")
    res <- bind_rows(mrres, obsres) %>%
        rename(beta=b)
    res <- bind_cols(params, res)
    return(list(starting_conditions=starting_conditions, dyn=dyn, res=res))
}




# =================== Starting values and parameters ======================


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
simulate_starting_conditions <- function(phen_names, n_id, phen_mean, phen_var, phen_h2, genotype_af) 
{
    set.seed(1234)
    n_phen <- length(phen_names)
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
    phen <- do.call(cbind, l) %>% as_tibble()
    names(phen) <- phen_names
    return(list(
        geno=geno,
        phen=phen
    ))
}


simulate_params <- function(par_names, n_id, par_mean, par_var, par_h2, genotype_af) 
{
    set.seed(1234)
    n_par <- length(par_names)
    stopifnot(length(par_mean) == n_par)
    stopifnot(length(par_var) == n_par)
    stopifnot(length(par_h2) == n_par)
    stopifnot(length(genotype_af) == n_par)

    confounder <- rnorm(n_id)
    confounder_var <- runif(n_par, 0, (1 - par_h2) * par_var)
    residual_var <- (par_var - par_h2 * par_var - confounder_var)
    
    colSums(rbind(residual_var, confounder_var, par_h2 * par_var))
    
    geno <- sapply(genotype_af, function(x) rbinom(n_id, 2, x))
    b_gx <- sqrt(par_h2 * par_var / (2 * genotype_af * (1- genotype_af)))
    
    l <- list()
    for(i in 1:n_par)
    {
        l[[i]] <- par_mean[i] + (geno[,i] - mean(geno[,i])) * b_gx[i] + confounder * sqrt(confounder_var[i]) + rnorm(n_id, 0, sqrt(residual_var[i]))
    }
    par <- do.call(cbind, l) %>% as_tibble()
    names(par) <- par_names
    return(list(
        geno=geno,
        par=par
    ))
}





# MR analysis

mr_analysis <- function(geno, out, exposures, outcomes, instrument_threshold=1e-4)
{
    analyses <- expand.grid(exposure=exposures, outcome=outcomes, stringsAsFactors=FALSE)
    d <- lapply(1:nrow(analyses), function(i){
        simulateGP::merge_exp_out(
            simulateGP::gwas(out[[analyses$exposure[i]]], geno),
            simulateGP::gwas(out[[analyses$outcome[i]]], geno),
            xname=analyses$exposure[i],
            yname=analyses$outcome[i]
        ) %>% filter(pval.exposure < instrument_threshold) %>%
        TwoSampleMR::mr(., method_list=c("mr_wald_ratio", "mr_ivw")) %>%
        mutate(method="MR") %>%
        dplyr::select(-id.exposure, -id.outcome)
    }) %>% bind_rows() %>%
        rename(mr_method=method)
    return(d)
}

# Observational analysis

obs_analysis <- function(out, exposures, outcomes)
{
    analyses <- expand.grid(exposure=exposures, outcome=outcomes, stringsAsFactors=FALSE)
    d <- lapply(1:nrow(analyses), function(i){
        temp <- tibble(
            exposure = out[[analyses$exposure[i]]],
            outcome = out[[analyses$outcome[i]]]
        )
        summary(glm(outcome ~ exposure, temp, family="binomial")) %>%
        coefficients() %>%
        {
            tibble(
                exposure=analyses$exposure[i],
                outcome=analyses$outcome[i],
                b=.[2,1],
                se=.[2,2],
                pval=.[2,4],
                method = "Obs"
            )
        } %>% return()
    }) %>% bind_rows()
    return(d)
}

per_analysis <- function(out_null, out_per, tp)
{
    lapply(names(out_per), function(i)
    {
        dper <- subset(out_per[[i]], time==tp)
        dnull <- subset(out_null, time==tp)
        phendiff <- mean(dper[[i]]) - mean(dnull[[i]])
        dper$T <- 1
        dnull$T <- 0
        d <- bind_rows(dper, dnull)
        mod <- summary(glm(D ~ T, data=d, family="binomial"))
        res <- tibble(
            exposure=i, 
            outcome="D",
            b = coefficients(mod)[2,1] / phendiff,
            se = coefficients(mod)[2,2] / phendiff,
            pval = coefficients(mod)[2,4],
            method="RCT"
        )
        return(res)
    }) %>% bind_rows()
}


## plotting


recreate_data <- function(dyn, dyn_per)
{
    tp <- min(dyn_per$time)
    dyn <- subset(dyn, time < tp)
    return(rbind(dyn, dyn_per))
}


plot_dynamics <- function(dyn, phen, w)
{
    dyn<-as.data.frame(dyn)
    phen<-as.data.frame(phen)
    #sb <- subset(dyn, subset=(Xt==phen[w,1]&Kt==phen[w,2]&Pt==phen[w,3]))
    sb <- subset(dyn, subset=(id==w))
    #print(sb)
    matplot(sb[,1],sb[,2:5], type="l",lty=1,col = 1:4)
    abline(h=sb$Ys, col="black",lwd=1, lty=2)
    legend("right",  colnames(sb[,2:7]),col=seq_len(6),cex=0.8,fill=seq_len(6),bty = "n", xpd=TRUE)
}

