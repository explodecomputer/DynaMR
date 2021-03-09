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
library(dplyr)

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
phen_mean = c(0.74,1.2,1)
phen_var = c(0.15,0.12,0.17)
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

dat$V1 <- dat$V1 + min(dat$V1)

#outfile <- file.path(output, "phen.txt")
#write.table(dat, file="phen.txt", row=F, col=T, qu=F)
#save(dat, file="phen.rdata")







# ===================== Dynamics simulation ==========================


###########################

library(data.table)

#gp <- fread("phen.txt", header=TRUE)
#load("phen.rdata")
gp <- dat

phen <- as.data.frame(gp$phen)
phen <- phen + abs(min(phen))+0.05

#solving an ode

require(deSolve)

phosph=function(Time,State,Params){
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




out <- c()

for (i in 1:nrow(phen)){
        tmp <- phen[i,]
        Xt = tmp[,1]
        Kt = tmp[,2]
        Pt = tmp[,3]
        a <- 0.18
        b <- 0.27
        c <- 0.0747
        e <- 0.18
        f <- 0.27
        g <- 0.0585
        odesol <- ode(func=phosph,y=c(X=Xt,Y=0,K=Kt,XK=0,P=Pt,YP=0),parms=c(a = a,b = b,c = c,e = e,f = f,g = g),times=seq(0, 200, by = 1))
        sst <- steadystate(a,b,c,e,f,g,Xt,Kt,Pt)
        if (sst>=0.45) {
                ds <- 1
        } else {
                ds <- 0
        }
        tb <- cbind(odesol,Xt=Xt,Kt=Kt,Pt=Pt,Ys=sst,D=ds)
        out <- rbind(out,tb)
}

out <- as.data.frame(out)

# Plot
w <- 4
sb <- subset(out, subset=(Xt==phen[w,1]&Kt==phen[w,2]&Pt==phen[w,3]))
matplot(sb[,1],sb[,2:7], type="l",lty=1,col = 1:6)
abline(h=sb[,11], col="black",lwd=1, lty=2)
legend("bottomright",  colnames(sb[,2:7]),col=seq_len(6),cex=0.8,fill=seq_len(6),bty = "n", xpd=TRUE)

###