library(data.table)
library(tidyverse)
library(pbapply)
library(simulateGP)
library(TwoSampleMR)
require(deSolve)
library("scatterplot3d")
require(deSolve)
library("rgl")


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




# ===================== Dynamics simulation ==========================


#solving ode system

ode_sys <- function(Time, State, parms)
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

## steady-state: root is the condition where sum of |rates of change|
## is very small



rootfun <- function (Time, State, parms) {
        dstate <- unlist(ode_sys(Time, State, parms)) # rate of change vector return(sum(abs(dstate)) - 1e-10)
        return(sum(abs(dstate)) - 1e-10)
}


# simulate dynamics without perturbation

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

        odesol <- ode(func=ode_sys,y=c(p=p0,a=a0,m=m0,f=f0), parms=c(Pbp=Pbp, Ppp=Ppp, Pfp=Pfp, App=App, Aph=Aph, Afp=Afp, Afh=Afh, Mbp=Mbp, Mpp=Mpp, Mph=Mph, Fdam=Fdam, gp=gp, gm=gm, gf=gf, n=n), times=seq(0, 100, by = 1))
        
        ic <- list(p0=p0, a0=a0, m0=m0, f0=f0)

        pars <- list(Pbp=Pbp, Ppp=Ppp, Pfp=Pfp, App=App, Aph=Aph, Afp=Afp, Afh=Afh, Mbp=Mbp, Mpp=Mpp, Mph=Mph, Fdam=Fdam, gp=gp, gm=gm, gf=gf, n=n)

        tout <- c(0, 1e10)

        State <- c(p = ic$p0, a = ic$a0, m = ic$m0, f = ic$f0)

        out <- suppressMessages(suppressWarnings(lsodar(State, tout, ode_sys, pars, rootfun = rootfun, verbose=FALSE)))

        sst=as.data.frame(t(out[2,]))

        if (sst$p>=disease_threshold) {
            ds <- 1
        } else {
            ds <- 0
        }
        tb <- cbind(odesol,p0=p0,a0=a0,m0=m0,f0=f0,ps=sst$p,as=sst$a,ms=sst$m,fs=sst$f,D=ds) %>%
            as_tibble() %>%
            mutate(id=i)
        return(tb)
    }) %>% bind_rows()
}

<<<<<<< HEAD
=======
# simulate dynamics with perturbation!!!!!!!!!!!!!!!!!!!!

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

        odesol <- ode(func=ode_sys,y=c(X=X0,Y=Y0,K=K0,XK=XK0,P=P0,YP=YP0),parms=c(a = parm_a,b = parm_b,c = parm_c,e = parm_e,f = parm_f,g = parm_g),times=seq(time_start, time_end, by = time_steps))
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

>>>>>>> 4eb3a3ae05592178e7000bfb41e57d78ff150528
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
<<<<<<< HEAD
    if (nrow(d)==0) {
        d<-data.frame(outcome=NA,exposure=NA,mr_method=NA,nsnp=NA,b=NA,se=NA,pval=NA,method="MR")
    }
=======
>>>>>>> 4eb3a3ae05592178e7000bfb41e57d78ff150528
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


<<<<<<< HEAD
simulation_per <- function(params, starting_condition_parameters, consts, dyn)
{
    # starting conditions

    cond <- subset(starting_condition_parameters, scenario==consts$scenario)

    pert_v <- c(starting_condition_parameters$variables,colnames(params$par))

    # Dynamics after the perturbation
    dyn_per <- lapply(1:length(pert_v), function(i)
    {
        message("perturbing ", pert_v[i])
        dyn1 <- subset(dyn, time == consts$per_timepoint)
        v <- pert_v[i]

        if (v %in% c("p","a","m","f")) {
            dyn1[[v]] <- dyn1[[v]] * (1 + consts$per_effect)
            dyn_per <- simulate_dynamics(
                starting_conditions = as.data.frame(dyn1[,c(2:5)]),
                params = params$par,
                n = consts$parm_n,
                disease_threshold = consts$disease_threshold,
                time_start = consts$time_start,
                time_end = consts$time_end,
                time_steps = consts$time_steps
            )
        } else {
            pars = as.data.frame(params$par)
            pars[v] <- pars[v] * (1 + consts$per_effect)
            dyn_per <- simulate_dynamics(
                starting_conditions = as.data.frame(dyn1[,c(2:5)]),
                params = pars,
                n = consts$parm_n,
                disease_threshold = consts$disease_threshold,
                time_start = consts$time_start,
                time_end = consts$time_end,
                time_steps = consts$time_steps
            )
        }

        return(dyn_per)
    })
    names(dyn_per) <- pert_v
=======
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
>>>>>>> 4eb3a3ae05592178e7000bfb41e57d78ff150528
    return(dyn_per)
}

## overall simulation

simulation <- function(cond_p, starting_condition_parameters, consts)
{
    cond <- subset(starting_condition_parameters, scenario==consts$scenario)

    message("Generating starting conditions")
    starting_conditions <- simulate_starting_conditions(
        phen_names=cond$variables, 
        n_id=consts$n_id, 
        phen_mean=cond$phen_mean, 
        phen_var=cond$phen_var, 
        phen_h2=cond$phen_h2, 
        genotype_af=cond$genotype_af
    )
    
    starting_conditions$phen <- t(t(starting_conditions$phen) + abs(apply(starting_conditions$phen,2,min))+apply(starting_conditions$phen,2,max)*0.001)

    message("Generating parameters")
    params <- simulate_params(
        n_id=consts$n_id, 
        par_names=cond_p$variables, 
        par_mean=cond_p$par_mean, 
        par_var=cond_p$par_var, 
        par_h2=cond_p$par_h2, 
        genotype_af=cond_p$genotype_af
    )

    params$par <- t(t(params$par) + abs(apply(params$par,2,min))+apply(params$par,2,max)*0.001)


    message("General simulation")
    dyn <- simulate_dynamics(
        starting_conditions = starting_conditions$phen,
        params = params$par,
        n = consts$parm_n,
        disease_threshold = consts$disease_threshold,
        time_start = consts$time_start,
        time_end = consts$time_end,
        time_steps = consts$time_steps
    )

<<<<<<< HEAD
    message("Perturbations")
    dyn_per <- simulation_per(
        params,
        starting_condition_parameters,
        consts,
        dyn
    )

    message("Analysis")
    # Analyse
    dyn1 <- subset(dyn, time == consts$time_end)
    # remove columns for dynamical p,a,m,f and rename ps,as,ms,fs -> p,a,m,f for the follwoing analysis
    dyn1 <- dyn1[,-c(2:5)]
    names(dyn1)[6:9] <- c("p","a","m","f")
    mrres <- mr_analysis(geno=starting_conditions$geno, out = dyn1, exposures=starting_condition_parameters$variables, outcomes="D")
    obsres <- obs_analysis(out = dyn1, exposures=starting_condition_parameters$variables, outcomes="D")
    perres <- per_analysis(dyn, dyn_per, consts$per_timepoint)
    res <- bind_rows(mrres, obsres, perres) %>%
        rename(beta=b)
    res <- bind_cols(consts, res)
    return(list(starting_conditions=starting_conditions, dyn=dyn, dyn_per=dyn_per, res=res))
    #return(list(starting_conditions=starting_conditions, dyn=dyn))
=======

    #message("Analysis")
    # Analyse
    #dyn1 <- subset(dyn, time == params$analysis_timepoint)
    #mrres <- mr_analysis(geno=starting_conditions$geno, out = dyn1, exposures=starting_condition_parameters$variables, outcomes="D")
    #obsres <- obs_analysis(out = dyn1, exposures=starting_condition_parameters$variables, outcomes="D")
    #perres <- per_analysis(dyn, dyn_per, params$per_timepoint)
    #res <- bind_rows(mrres, obsres, perres) %>%
    #    rename(beta=b)
    #res <- bind_cols(params, res)
    #return(list(starting_conditions=starting_conditions, dyn=dyn, dyn_per=dyn_per, res=res))
    return(list(starting_conditions=starting_conditions, dyn=dyn))
>>>>>>> 4eb3a3ae05592178e7000bfb41e57d78ff150528
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

plot_3d_1 <- function(dyn)
{
    dyn<-as.data.frame(dyn)
    t3 <- c()
    for(w in 1:1000)
    {
        t1 <- subset(dyn, subset=(id==w))
        t2 <- t1[1,]
        t3 <- rbind(t3,t2)
    }
    t3 <- t3[complete.cases(t3),]
    colors<-rgb(0,t3$ps/max(t3$ps),0)
    scatterplot3d(log10(t3[,c(11:13)]), pch = 16, color=colors)
}

rgl_init <- function(new.device = FALSE, bg = "white", width = 640) { 
    if( new.device | rgl.cur() == 0 ) {
        rgl.open()
        par3d(windowRect = 50 + c( 0, 0, width, width ) )
        rgl.bg(color = bg )
    }
    rgl.clear(type = c("shapes", "bboxdeco"))
    rgl.viewpoint(theta = 15, phi = 20, zoom = 0.7)
}

rgl_add_axes <- function(x, y, z, axis.col = "grey",
                         xlab = "", ylab="", zlab="", show.plane = TRUE, 
                         show.bbox = FALSE, bbox.col = c("#333377","black"))
{ 
    
    lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
    # Add axes
    xlim <- lim(x); ylim <- lim(y); zlim <- lim(z)
    rgl.lines(xlim, c(0, 0), c(0, 0), color = axis.col)
    rgl.lines(c(0, 0), ylim, c(0, 0), color = axis.col)
    rgl.lines(c(0, 0), c(0, 0), zlim, color = axis.col)
    
    # Add a point at the end of each axes to specify the direction
    axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0), 
                  c(0, 0, zlim[2]))
    rgl.points(axes, color = axis.col, size = 3)
    
    # Add axis labels
    rgl.texts(axes, text = c(xlab, ylab, zlab), color = axis.col,
              adj = c(0.5, -0.8), size = 2)
    
    # Add plane
    if(show.plane) 
        xlim <- xlim/1.1; zlim <- zlim /1.1
    rgl.quads( x = rep(xlim, each = 2), y = c(0, 0, 0, 0),
               z = c(zlim[1], zlim[2], zlim[2], zlim[1]))
    
    # Add bounding box decoration
    if(show.bbox){
        rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, 
                 emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
                 xlen = 3, ylen = 3, zlen = 3) 
    }
}

plot_3d_2 <- function(dyn)
{
    dyn<-as.data.frame(dyn)
    t3 <- c()
    for(w in 1:1000)
    {
        t1 <- subset(dyn, subset=(id==w))
        t2 <- t1[1,]
        t3 <- rbind(t3,t2)
    }
    t3 <- t3[complete.cases(t3),]
    colors<-rgb(0,t3$ps/max(t3$ps),0)

    x <- as <- log10(t3$as)
    y <- ms <- log10(t3$ms)
    z <- fs <- log10(t3$fs)

    rgl_init()
    rgl.spheres(x, y, z, r = 0.2, color = colors) 
    rgl_add_axes(x, y, z, show.bbox = TRUE, xlab ="as", ylab = "ms", zlab = "fs")
    aspect3d(1,1,1)
    
}

