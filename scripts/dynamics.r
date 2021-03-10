library(data.table)
library(dplyr)
library(pbapply)
library(simulateGP)
library(TwoSampleMR)
require(deSolve)


# =================== Starting values ======================


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


# ===================== Dynamics simulation ==========================


###########################


#solving an ode


phosph <- function(Time,State,Params)
{
    with(as.list(c(State,Params)),{
    dX = -a*X*K + b*XK + g*YP
    dY = c*XK - e*Y*P + f*YP
    dK = -a*X*K + (b+c)*XK
    dXK = a*X*K - (b+c)*XK
    dP = -e*Y*P + (f+g)*YP
    dYP = e*Y*P - (f+g)*YP
    
    res<-list(c(dX,dY,dK,dXK,dP,dYP))
    
    return(res)
    })
}

fs <- function(x) {
    if (Im(x)>1e-13) {
        x <- NaN
    } else {
        x <- Re(x)
    }
    return(x)
}

steadystate <- function(a,b,c,e,f,g,Xt,Kt,Pt) 
{
    Kp <- (b+c)/a
    Kf <- (g+f)/e
    h1 <- Kp/Xt
    h2 <- Kf/Xt
    h3 <- c/g
    h4 <- Kt/Pt
    h5 <- Pt/Xt
    
    u <- (1+h1+h2)/(1+h3)+h4*h5+h5/h3
    v <- ((1+h2)*h4*h5+(1+h1)*h5/h3)/(1+h3)+h4*h5*h5/h3
    w <- h4*h5*h5/h3/(1+h3)
    
    D <- sqrt(as.complex(3*((4*v-u^2)*v^2-w*(18*v*u-27*w-4*u^3))))
    DD <- as.complex(4*(3*D-27*w-u*(2*u^2-9*v)))^(1/3)
    
    x1 <- 2*(v-1/3*u^2)/DD + 1/3*(u-1/2*DD)
    
    x1 <-fs(x1)
    
    q1 <- 1-h3*h4
    q2 <- h1-x1*(h3+1)+1+h3*h4*(h2+x1*(h3+1)-1)
    q3 <- h2*h3*h4*(1-x1*(h3+1))
    y <- 1/(2*q1)*(q2-sqrt(as.complex(q2^2-4*q1*q3)))
    
    y <- fs(y)
    
    y <- y*Xt
    
    return(y)
}


simulate_dynamics <- function(starting_conditions, parm_a, parm_b, parm_c, parm_e, parm_f, parm_g, disease_threshold, time_start, time_end, time_steps)
{
    starting_conditions <- lapply(starting_conditions, function(x){
        x + abs(min(x)) + 0.05
    }) %>% bind_cols()
    out_full <- pblapply(1:nrow(starting_conditions), function(i)
    {
        tmp <- starting_conditions[i,]
        Xt = tmp$V1
        Kt = tmp$V2
        Pt = tmp$V3
        odesol <- ode(func=phosph,y=c(X=Xt,Y=0,K=Kt,XK=0,P=Pt,YP=0),parms=c(a = parm_a,b = parm_b,c = parm_c,e = parm_e,f = parm_f,g = parm_g),times=seq(time_start, time_end, by = time_steps))
        sst <- steadystate(parm_a,parm_b,parm_c,parm_e,parm_f,parm_g,Xt,Kt,Pt)
        if (sst>=disease_threshold) {
            ds <- 1
        } else {
            ds <- 0
        }
        o <- cbind(odesol,Xt=Xt,Kt=Kt,Pt=Pt,Ys=sst,D=ds) %>%
            as_tibble() %>%
            mutate(id=i)
        return(o)

    }) %>% bind_rows()
    return(out_full)
}


mr_analysis <- function(geno, out, exposures, outcomes, pval_threshold=1e-4)
{
    analyses <- expand.grid(exposure=exposures, outcome=outcomes)
    d <- lapply(1:nrow(analyses), function(i){
        simulateGP::merge_exp_out(
            simulateGP::gwas(out[[analyses$exposure[i]]], geno),
            simulateGP::gwas(out[[analyses$outcome[i]]], geno),
            xname=analyses$exposure[i],
            yname=analyses$outcome[i]
        ) %>% filter(pval.exposure < pval_threshold) %>%
        TwoSampleMR::mr(., method_list=c("mr_wald_ratio", "mr_ivw"))
    }) %>% bind_rows()
}



simulation <- function(params, starting_condition_parameters)
{
    cond <- subset(starting_condition_parameters, scenario==params$scenario)
    starting_conditions <- simulate_starting_conditions(
        n_id=params$n_id, 
        n_phen=cond$n_phen[1], 
        phen_mean=cond$phen_mean, 
        phen_var=cond$phen_var, 
        phen_h2=cond$phen_h2, 
        genotype_af=cond$genotype_af
    )
    dyn <- simulate_dynamics(
        starting_conditions = starting_conditions$phen,
        parm_a = params$a,
        parm_b = params$b,
        parm_c = params$c,
        parm_e = params$e,
        parm_f = params$f,
        parm_g = params$g,
        disease_threshold = params$disease_threshold,
        time_start = params$time_start,
        time_end = params$time_end,
        time_steps = params$time_steps
    )
    # simulate_perturbation_reaction()
    dyn1 <- subset(dyn, time == params$analysis_timepoint)
    res <- mr_analysis(geno=starting_conditions$geno, out = dyn1, exposures=c("Xt", "Kt", "Pt"), outcomes="Ys")
    return(list(starting_conditions=starting_conditions, dyn=dyn, mr_res=res))
}

