library(data.table)
library(tidyverse)
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


#solving ode system

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

#getting rid of small complex values

fs <- function(x) {
    if (Im(x)>1e-13) {
        x <- NaN
    } else {
        x <- Re(x)
    }
    return(x)
}

#find steady state concentration in analytical form

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


# simulate dynamics without perturbation

simulate_dynamics <- function(starting_conditions, parm_a, parm_b, parm_c, parm_e, parm_f, parm_g, disease_threshold, time_start, time_end, time_steps)
{
    #starting_conditions <- lapply(starting_conditions, function(x){
    #    x + abs(min(x)) + 0.05
    #}) %>% bind_cols()
    out_full <- pblapply(1:nrow(starting_conditions), function(i)
    {
        tmp <- starting_conditions[i,]

        # Starting concentrations
        X0 = tmp$V1
        Y0 = 0
        K0 = tmp$V2
        XK0 = 0
        P0 = tmp$V3
        YP0 = 0

        # Total concentrations
        Xt = X0
        Kt = K0
        Pt = P0

        odesol <- ode(func=phosph,y=c(X=X0,Y=Y0,K=K0,XK=XK0,P=P0,YP=YP0),parms=c(a = parm_a,b = parm_b,c = parm_c,e = parm_e,f = parm_f,g = parm_g),times=seq(time_start, time_end, by = time_steps))
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

# simulate dynamics with perturbation

simulate_dynamics_per <- function(starting_conditions, parm_a, parm_b, parm_c, parm_e, parm_f, parm_g, disease_threshold, time_start, time_end, time_steps)
{
    #starting_conditions <- lapply(starting_conditions, function(x){
    #    x + abs(min(x)) + 0.05
    #}) %>% bind_cols()
    out_full <- pblapply(1:nrow(starting_conditions), function(i)
    {
        tmp <- starting_conditions[i,]

        # Starting concentrations
        X0 = tmp$X
        Y0 = tmp$Y
        K0 = tmp$K
        XK0 = tmp$XK
        P0 = tmp$P
        YP0 = tmp$YP

        # Total concentrations
        Xt = X0+Y0+XK0+YP0
        Kt = K0+XK0
        Pt = P0+YP0

        odesol <- ode(func=phosph,y=c(X=X0,Y=Y0,K=K0,XK=XK0,P=P0,YP=YP0),parms=c(a = parm_a,b = parm_b,c = parm_c,e = parm_e,f = parm_f,g = parm_g),times=seq(time_start, time_end, by = time_steps))
        sst <- steadystate(parm_a,parm_b,parm_c,parm_e,parm_f,parm_g,Xt,Kt,Pt)
        ds <- 
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

# MR analysis

mr_analysis <- function(geno, out, exposures, outcomes, instrument_threshold=1e-4)
{
    analyses <- expand.grid(exposure=exposures, outcome=outcomes, stringsAsFactors=FALSE)
    d <- lapply(1:nrow(analyses), function(i){
        simulateGP::merge_exp_out(
            simulateGP::gwas(out[[analyses$exposure[i]]], geno),
            simulateGP::gwas(out[[analyses$outcome[i]]], geno, logistic=TRUE),
            xname=analyses$exposure[i],
            yname=analyses$outcome[i]
        ) %>% filter(pval.exposure < instrument_threshold) %>%
        TwoSampleMR::mr(., method_list=c("mr_wald_ratio", "mr_ivw")) %>%
        rename(mr_method=method) %>%
        mutate(method="MR") %>%
        dplyr::select(-id.exposure, -id.outcome)
    }) %>% bind_rows()
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


# ===================== Dynamics simulation with perturbation ==========================


simulation_per <- function(params, starting_condition_parameters, dyn)
{
    # starting conditions

    cond <- subset(starting_condition_parameters, scenario==params$scenario)

    # Dynamics after the perturbation
    dyn_per <- lapply(1:nrow(starting_condition_parameters), function(i)
    {
        message("perturbing ", starting_condition_parameters$variables[i])
        dyn1 <- subset(dyn, time == params$per_timepoint)
        v <- starting_condition_parameters$variables[i]
        dyn1[[v]] <- dyn1[[v]] * (1 + params$per_effect)
        dyn_per <- simulate_dynamics_per(
            starting_conditions = dyn1,
            parm_a = params$a,
            parm_b = params$b,
            parm_c = params$c,
            parm_e = params$e,
            parm_f = params$f,
            parm_g = params$g,
            disease_threshold = params$disease_threshold,
            time_start = params$per_timepoint,
            time_end = params$time_end,
            time_steps = params$time_steps
        )
        return(dyn_per)
    })
    names(dyn_per) <- starting_condition_parameters$variables
    return(dyn_per)
}

## overall simulation

simulation <- function(params, starting_condition_parameters)
{
    cond <- subset(starting_condition_parameters, scenario==params$scenario)

    message("Generating starting conditions")
    starting_conditions <- simulate_starting_conditions(
        n_id=params$n_id, 
        n_phen=cond$n_phen[1], 
        phen_mean=cond$phen_mean, 
        phen_var=cond$phen_var, 
        phen_h2=cond$phen_h2, 
        genotype_af=cond$genotype_af
    )
    
    starting_conditions$phen <- starting_conditions$phen + abs(min(starting_conditions$phen)) + 0.05

    message("General simulation")
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

    message("Perturbations")
    dyn_per <- simulation_per(
        params,
        starting_condition_parameters,
        dyn
    )

    message("Analysis")
    # Analyse
    dyn1 <- subset(dyn, time == params$analysis_timepoint)
    mrres <- mr_analysis(geno=starting_conditions$geno, out = dyn1, exposures=starting_condition_parameters$variables, outcomes="D")
    obsres <- obs_analysis(out = dyn1, exposures=starting_condition_parameters$variables, outcomes="D")
    perres <- per_analysis(dyn, dyn_per, params$per_timepoint)
    res <- bind_rows(mrres, obsres, perres) %>%
        rename(beta=b)
    res <- bind_cols(params, res)
    return(list(starting_conditions=starting_conditions, dyn=dyn, dyn_per=dyn_per, res=res))
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
    matplot(sb[,1],sb[,2:7], type="l",lty=1,col = 1:6)
    abline(h=sb$Ys, col="black",lwd=1, lty=2)
    legend("right",  colnames(sb[,2:7]),col=seq_len(6),cex=0.8,fill=seq_len(6),bty = "n", xpd=TRUE)
}

