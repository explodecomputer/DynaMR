library(dplyr)
library(rootSolve)
library(data.table)


# ===================== Dynamics simulation ==========================


###########################

RJ2 <- function(i2, pars){
        x = M2-r2t+p2t*(M1/(n1*i2)+1/n1+1)
        d = r2t*M2
        eps = 2.2204e-16
        er_ap = d/(x^2)
        er_ex = eps*(x^2)/(2*d)
        
        if (er_ap < er_ex) {
                w2 = d/x - 1/x*(d/x)^2
        } else {
                w2 = -x/2 + sqrt(x^2+4*d)/2
        }
        
        A = w2
        return(A)
}



#solving an ode

require(deSolve)

ode_sys=function(Time,State,Params){
with(as.list(c(State,Params)),{
        
        s33 = s3p^2/M13
        s35 = s3p*s5p/M14
        s55 = s5p^2/M15
        p3s3p = p3t*s3p/(M9+s3p)
        p5s5p = p5t*s5p/(M12+s5p)
        
        w2 = RJ2(i2, pars)
        
        s3 = (1 - s3p - 2*s33 - s35 - p3s3p)/(1 + w2/M7 + Q6)
        s5 = (s5t - s5p - 2*s55 - s35 - p5s5p)/(1 + w2/M10 + Q21)
        
        w2s3 = w2*s3/M7
        w6s3 = Q6*s3
        w2s5 = w2*s5/M10
        w21s5 = Q21*s5
        
        ds3p = m7*(n4*w2s3 + n5*w6s3 - p3s3p)
        ds5p = m16*(n6*w2s5 + n7*w21s5 - p5s5p)
        
        res<-list(c(ds3p,ds5p))
        
        return(res)
        })
}


# =================== Starting values ======================
s3p0 <- 0
s5p0 <- 0

# =================== Parameters ======================
# Receptor
r2t = 0.0027 # R2t/S3t
p2t = 0.0027 # P2t/S3t
M1 = 0.1333 # (f2+f3)/(f1*S3t) : phosphorylation
M2 = 4.2474e-05 # (f5+f6)/(f4*S3t) : dephosphorylation
n1 = 118.4200 # f2/f5 : phosphorylation/dephosphorylation

# Total amounts of STAT5 and STAT phosphatases
s5t = 0.0247 # S5t/S3t
p3t = 2.5924 # P3t/S3t
p5t = 0.0012 # P5t/S3t

# STAT3
M7 = 3.6348e-04 # (a2+a3)/(a1*S3t) : phosphorylation by IL-2
M9 = 47.7140 # (a8+a9)/(a7*S3t) : dephosphorylation
M13 = 19.1540 # a20/(a19*S3t) : homodimer formation
n4 = 0.1987 # a2/a8 : phosphorylation by IL-2/dephosphorylation
n5 = 1.5589 # a5/a8 : phosphorylation by IL-6/dephosphorylation
m7 = 1.005 # a8/(a2+a3) : rate of dephosphorylation

# Heterodimer
M14 = 0.1002 # a22/(a21*S3t) : S35 formation

# STAT5
M10 = 0.0047 # (a11+a12)/(a10*S3t) : phosphorylation by IL-2
M12 = 1.9586e+03 # (a17+a18)/(a16*S3t) : dephosphorylation
M15 = 0.3787 # a24/(a23*S3t) : homodimer formation
n6 = 5.5056 # a11/a17 : phosphorylation by IL-2/dephosphorylation
n7 = 0.0322 # a14/a17 : phosphorylation by IL-6/dephosphorylation
m16 = 1.01 # a17/(a2+a3) : rate of dephosphorylation

# IL-6 signal
Q6 = 0.0014 # RpJ6*a4/(a5+a6)

# IL-21 signal
Q21 = 3.4207e-11 # RpJ21*a13/(a14+a15)

# Cytokine production
ggt = 0.8949 # Ggt/S3t
mp1t = 0.0034 # Mp1t/S3t
M18 = 9.9130 # (k4+k5)/(k3*S3t) : inactivation by MP
M19 = 4.3540 # k2/(k1*S3t) : STAT55 complex formation with Gg
n8 = 0.0101 # l5/k4 : rate of IFN-g production/rate of inactivation by MP, <1

g10t = 7.2184 # G10t/S3t
mp2t = 0.5933 # Mp2t/S3t
M20 = 0.0149 # (k11+k12)/(k10*S3t) : inactivation by MP
M21 = 0.0138 # k7/(k6*S3t) : STAT33 complex formation with G10
M22 = 0.1885 # k9/(k8*S3t) : SP1 complex formation with G10
n9 = 0.0191 # l6/k11 : rate of IL-10 production/rate of inactivation by MP, <1

# CD46
cd46 = 0.6826 # CD46/S3t
sp1t = 33.1420 # Sp1t/S3t
M16 = 8.9570e-06 # l2/(l1*S3t) : sensitivity to IL-2
M17 = 0.1071 # l4/(l3*S3t)

# Initial IL-2
i2 = 1e-10

# ======== Dynamics ================

odesol_0 <- ode(func=ode_sys,y=c(s3p=s3p0,s5p=s5p0),parms=c(r2t, p2t, M1, M2, n1, s5t, p3t, p5t, M7, M9, M13, n4, n5, m7, M14, M10, M12, M15, n6, n7, m16, Q6, Q21, ggt, mp1t, M18, M19, n8, g10t, mp2t, M20, M21, M22, n9, cd46, sp1t, M16, M17, i2), times=seq(0, 200, by = 1))
odesol_1 <- as.data.frame(odesol_0)
odesol_1$s5p <- odesol_1$s5p/s5t

odesol_1$s33 <- odesol_1$s3p^2/M13
odesol_1$s35 <- odesol_1$s3p*odesol_1$s5p/M14
odesol_1$s55 <- (odesol_1$s5p^2/M15)*s5t

sp1 = sp1t*i2/(M16+i2)*cd46/(M17+cd46)

odesol_1$ifn <- M18/(mp1t/(n8*ggt*odesol_1$s55*s5t/(M19+odesol_1$s55*s5t))-1)
odesol_1$il10 <- M20/(mp2t/(n9*g10t*(odesol_1$s33/(M21+odesol_1$s33)+sp1/(M22+sp1)-odesol_1$s33/(M21+odesol_1$s33)*sp1/(M22+sp1)))-1)

# Plot dynamics
matplot(odesol_1[,1],odesol_1[,7:8], type="l",lty=1,col = 1:6)
#abline(h=sb[,12], col="black",lwd=1, lty=2)
legend("bottomright",  colnames(odesol_1[,7:8]),col=seq_len(6),cex=0.8,fill=seq_len(6),bty = "n", xpd=TRUE)

# ======== Steady-states ================

## root is the condition where sum of |rates of change|
## is very small
rootfun <- function (Time, State, Params) {
        dstate <- unlist(ode_sys(Time, State, Params)) # rate of change vector return(sum(abs(dstate)) - 1e-10)
        return(sum(abs(dstate)) - 1e-10)
}

ic <- list(s3p0=s3p0, s5p0=s5p0)

tout <- c(0, 1e10)

State <- c(s3p = ic$s3p0, s5p = ic$s5p0)

# IL-2 is varied between 1e-10 to 1
r = c()
for (q in seq(-10, 0, length.out=200)) {
        i2 = 10^q
        pars <- list(r2t=r2t, p2t=p2t, M1=M1, M2=M2, n1=n1, s5t=s5t, p3t=p3t, p5t=p5t, M7=M7, M9=M9, M13=M13, n4=n4, n5=n5, m7=m7, M14=M14, M10=M10, M12=M12, M15=M15, n6=n6, n7=n7, m16=m16, Q6=Q6, Q21=Q21, ggt=ggt, mp1t=mp1t, M18=M18, M19=M19, n8=n8, g10t=g10t, mp2t=mp2t, M20=M20, M21=M21, M22=M22, n9=n9, cd46=cd46, sp1t=sp1t, M16=M16, M17=M17, i2=i2)
        out = lsodar(State, tout, ode_sys, pars, rootfun = rootfun)
        sst = as.data.frame(t(out[2,]))
        
        sst$s5p <- sst$s5p/s5t
        
        sst$s33 <- sst$s3p^2/M13
        sst$s35 <- sst$s3p*sst$s5p/M14
        sst$s55 <- (sst$s5p^2/M15)*s5t
        
        sp1 = sp1t*i2/(M16+i2)*cd46/(M17+cd46)
        
        sst$ifn <- M18/(mp1t/(n8*ggt*sst$s55*s5t/(M19+sst$s55*s5t))-1)
        sst$il10 <- M20/(mp2t/(n9*g10t*(sst$s33/(M21+sst$s33)+sp1/(M22+sp1)-sst$s33/(M21+sst$s33)*sp1/(M22+sp1)))-1)
        
        r = rbind(r,data.frame(i2=i2,s3ps=sst$s3p,s5ps=sst$s5p,s33s=sst$s33,s35s=sst$s35,s55s=sst$s55,ifns=sst$ifn,il10s=sst$il10))
}


# Plot steady-states
# STATs
matplot(log10(r[,1]),r[,2:3], type="l",lty=1,col = 1:6)
legend("bottomright",  colnames(odesol_1[,2:3]),col=seq_len(6),cex=0.8,fill=seq_len(6),bty = "n", xpd=TRUE)

# STAT dimers
matplot(log10(r[,1]),r[,c(4,6)], type="l",lty=1,col = 1:6)
legend("bottomright",  colnames(odesol_1[,c(4,6)]),col=seq_len(6),cex=0.8,fill=seq_len(6),bty = "n", xpd=TRUE)

# Cytokines
matplot(log10(r[,1]),r[,7:8], type="l",lty=1,col = 1:6)
legend("bottomright",  colnames(odesol_1[,7:8]),col=seq_len(6),cex=0.8,fill=seq_len(6),bty = "n", xpd=TRUE)



