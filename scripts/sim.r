parameters <- expand.grid(
    n_id = c(1000),
    scenario = 1,
    time_start = 0,
    time_end = 200,
    time_steps = 1,
    analysis_timepoint = 200,
    per_timepoint = 30,
    disease_threshold = 0.45,
    a = 0.18,
    b = 0.27,
    c = 0.0747,
    e = 0.18,
    f = 0.27,
    g = 0.0585,
    pval_threshold = 1e-4
)

starting_condition_parameters <- tibble(
    scenario = 1,
    n_phen = 3,
    phen_mean = c(0.74,1.2,1),
    phen_var = c(0.08,0.08,0.08),
    phen_h2 = c(0.2, 0.3, 0.4),
    genotype_af = c(0.2, 0.3, 0.4)
)



params <- parameters[1,]
source("dynamics.r")




# ---------- MR analysis without perturbation ------------

plot_dynamics <- function(dyn, phen, w)
{
    dyn<-as.data.frame(dyn)
    phen<-as.data.frame(phen)
    #sb <- subset(dyn, subset=(Xt==phen[w,1]&Kt==phen[w,2]&Pt==phen[w,3]))
    sb <- subset(dyn, subset=(id==w))
    #print(sb)
    png("rplot.png") 
    matplot(sb[,1],sb[,2:7], type="l",lty=1,col = 1:6)
    abline(h=sb$Ys, col="black",lwd=1, lty=2)
    legend("right",  colnames(sb[,2:7]),col=seq_len(6),cex=0.8,fill=seq_len(6),bty = "n", xpd=TRUE)
    dev.off() 
}

out <- simulation(params, starting_condition_parameters)

plot_dynamics(out$dyn, out$starting_conditions$phen, 4)



# ---------- Dyn analysis without perturbation ------------

per_X <- 0
per_K <- 0
per_P <- 0

out_per <- simulation_per(params, starting_condition_parameters, per_X, per_K, per_P)

plot_dynamics(out_per$dyn, out_per$starting_conditions$phen, 4)

dis <- out_per$sts
dis_n <- cbind(dis,T=0)



# ---------- Dyn analysis with perturbation ------------

per_X <- 0
per_K <- 0
per_P <- 0.2

out_per <- simulation_per(params, starting_condition_parameters, per_X, per_K, per_P)

plot_dynamics(out_per$dyn, out_per$starting_conditions$phen, 4)

dis <- out_per$sts
dis_p <- cbind(dis,T=1)


dis_all <- rbind(dis_n,dis_p)


summary(glm(D ~ T, data = dis_all, family="binomial"))









### old stuff



out <- subset(out_full, time==max(time))

# example
set.seed(1234)

n_phen = 3
n_id = 1000
dat <- simulate_starting_conditions(n_phen=n_phen, n_id=n_id, phen_mean=phen_mean, phen_var=phen_var, phen_h2=phen_h2, genotype_af=genotype_af)


# check
colMeans(dat$geno)/2
genotype_af

colMeans(dat$phen)
phen_mean

apply(dat$phen, 2, var)
phen_var

dat$V1 <- dat$V1 + min(dat$V1)

table(dat$phen < 0)

hist(dat$phen$V3)
hist(dat$phen$phen)
#outfile <- file.path(output, "phen.txt")
#write.table(dat, file="phen.txt", row=F, col=T, qu=F)
#save(dat, file="phen.rdata")


#gp <- fread("phen.txt", header=TRUE)
#load("phen.rdata")
gp <- dat

phen <- as.data.frame(gp$phen)
phen <- phen + abs(min(phen))+0.05

