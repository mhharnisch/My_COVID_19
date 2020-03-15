### COVID 19 Simulation ####
setwd("C:/Users/Stephan/Documents/GitHub/My_COVID_19/")

# Functions: --------------------------------------------------------------

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}





###### Empirics #######
# Read Data ---------------------------------------------------------------
e_dat_con <- read.csv("../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv")
e_dat_dea <- read.csv("../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Deaths.csv")
e_dat_rec <- read.csv("../COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Recovered.csv")



# First analysis ----------------------------------------------------------
country <- 'Germany'

plot(as.numeric(e_dat_con[e_dat_con$Country.Region==country,5:ncol(e_dat_con)]), ylab = 'Confirmed', xlab = "Date")
lines(as.numeric(e_dat_dea[e_dat_dea$Country.Region==country,5:ncol(e_dat_con)]), ylab = 'Confirmed', xlab = "Date", type = 'b', col = 2)
lines(as.numeric(e_dat_rec[e_dat_rec$Country.Region==country,5:ncol(e_dat_con)]), ylab = 'Confirmed', xlab = "Date", type = 'b', col = 3)

plot(as.numeric(e_dat_con[e_dat_con$Country.Region==country,5:ncol(e_dat_con)]), ylab = 'Confirmed', xlab = "Date", log = 'y')
lines(as.numeric(e_dat_dea[e_dat_dea$Country.Region==country,5:ncol(e_dat_con)]), ylab = 'Confirmed', xlab = "Date", type = 'b', col = 2)
lines(as.numeric(e_dat_rec[e_dat_rec$Country.Region==country,5:ncol(e_dat_con)]), ylab = 'Confirmed', xlab = "Date", type = 'b', col = 3)


# Estimate daily rate:

growth_rates <- function(x){
  return(x[2:length(x)]/x[1:(length(x)-1)]-1)
}


con_GER <- as.numeric(e_dat_con[e_dat_con$Country.Region==country,5:ncol(e_dat_con)])
gr_GER <- growth_rates(as.numeric(e_dat_con[e_dat_con$Country.Region==country,5:ncol(e_dat_con)]))

plot(gr_GER, type = 'h')

hist(gr_GER, breaks = seq(0,3,0.1))

plot(density(gr_GER[34:length(gr_GER)]))
mean(gr_GER[34:length(gr_GER)])

mean(gr_GER[10:length(gr_GER)], na.rm = T)


# Breakpoint estimation ---------------------------------------------------

log_gr_GER <- log(gr_GER)

gr_dat <- data.frame(t = 1:length(gr_GER), log_gr_GER = log_gr_GER, gr_GER = gr_GER)
gr_dat$log_gr_GER[is.finite(gr_dat$log_gr_GER)==F] <- NA

dat  <- data.frame(t = 1:length(con_GER), con_GER = as.numeric(con_GER))
dat$con_GER[dat$con_GER == 0] <- NA


(lm1 <- lm(log(con_GER)~t, data = dat))
exp(lm1$coefficients[2])-1

# With breakpoint at T

lm_res <- list()

i<- 1
for(T.break in 1:20){
  for(T.break2 in 21:nrow(dat)){
    lm_res[[i]] <- lm(log(con_GER)~ t + I(t>T.break)*t + I(t>T.break2)*t, data = dat)
    i<-i+1
  }
}

model_selct <- which.max(unlist(lapply(lm_res, function(x) summary(x)$r.squared)))

par(mar=c(4.1,4.1,1.1,1.1))
plot(con_GER)
lines(as.numeric(names(exp(lm_res[[model_selct]]$fitted.values))),exp(lm_res[[model_selct]]$fitted.values), col = 2, type = 'l')
legend('topleft', legend = c('confirmed cases', 'tri sectional log linear fit'), pch = c(1,NA), lty = c(NA, 1), col = c(1,2))

plot(con_GER, log = 'y')
lines(as.numeric(names(exp(lm_res[[model_selct]]$fitted.values))), exp(lm_res[[model_selct]]$fitted.values), col = 2, type = 'l')
legend('topleft', legend = c('confirmed cases', 'tri sectional log linear fit'), pch = c(1,NA), lty = c(NA, 1), col = c(1,2))

# Estimated growth rate:

sum_best_model <- summary(lm_res[[model_selct]])

a <- c(0,1,0,0,1,1)

daily_growth <- exp(a%*%sum_best_model$coefficients[,1])-1

log(2)/log(1+daily_growth)

# y = (1+gr)^t
# log(y) = t*log(1+gr)

###### Simulation #####
# Parameters: -------------------------------------------------------------

# population size:
n <- 2000

# initially infected
i0 <- 10

# number of days
Tmax <- 200

# initial average number of "close" encounters per person per day:
ane <- 80

### "Lockdown" parameters: 
# after measures are introduced reduction by factor
red_fac <- 0.50
ane_2 <- ane*(1-red_fac)
# day of introducing measures:
t_lockdown <- 200

# incubation period:
Tinc <- 10

# symptomatic period:
Tsym <-  10

# infectious period:
Tinf <- Tinc + Tsym

# town clusters:
n_town <- 1
within_town <- 1

# expected family size
fam_size <- 2
within_fam <- 1

# probability of overall infecting "R0" (expected value)
R0 <- 3.0

# probability of infecting in a single encounter
p <- R0/(ane*Tinf)


# Overview of infection distributions:
par(mar=c(4.1,4.1,2.1,1.1))
plot(seq(0,ane*Tinf,1), dbinom(seq(0,ane*Tinf,1), ane*Tinf, prob = p), type = 'b', xlab = 'Number of infected', ylab = 'Probability', main = 'Probability of infecting X persons in infectious period', cex.main = 0.8)
grid()

ane*Tinf*p

plot(seq(0,ane,1), dbinom(seq(0,ane,1), ane, prob = p), type = 'b', xlab = 'Number of infected', ylab = 'Probability', main = 'Probability of infecting X persons on one day', cex.main = 0.8)
grid()

lines(seq(0,ane_2,1), dbinom(seq(0,ane_2,1), ane_2, prob = p), type = 'b', xlab = 'Number of infected', ylab = 'Probability', main = 'Probability of infecting X persons on one day', cex.main = 0.8, col = 2)
grid()


# simulation parameters: --------------------------------------------------

simu <- 100
simu_res <- list()

library(progress)
pb <- progress_bar$new(total = simu)

# Set up simulation dataframes: -------------------------------------------
for(j in 1:simu){

#Data on who is infected:
inf_dat <- as.data.frame(matrix(0, nrow = n, ncol = Tmax))

#initially infected:
inf_dat[1:i0, 1] <- 1

#Individual number of meetings (constant troughout)
ine <- rpois(n, lambda = ane)

ine_2 <- round(ine*(1-red_fac))

### Probabilities of meeting with clusters for family & town

# matrix for relationship:
cluster_matrix <- matrix(1/n, ncol = n, nrow = n)


# Generate town clusters:
for(i in 1:(n_town)){
  ind_1 <- round((i-1)*n/n_town)+1
  ind_2 <- round((i)*n/n_town)
  cluster_matrix[ind_1:ind_2,ind_1:ind_2] <- cluster_matrix[ind_1:ind_2,ind_1:ind_2]*within_town
}

# Generate family clusters:
i<-1
while(i<=n){
  temp_fam_size<-rpois(n = 1, fam_size)
  
  ind_1 <- i
  
  if(i+temp_fam_size-1<=n){
    ind_2 <- i+temp_fam_size-1
  }else{
    ind_2 <- n
  }
  
  cluster_matrix[ind_1:ind_2,ind_1:ind_2] <- cluster_matrix[ind_1:ind_2,ind_1:ind_2]*within_fam
  
  i <- i + temp_fam_size
}

#Set diagonal to 0:
diag(cluster_matrix) <- 0

#renormalize matrix:
cluster_matrix_norm <- cluster_matrix/rowSums(cluster_matrix)

###407 481 4820

#list of who infects who (depending on round):
#(first dimension is time, second dimension is individual, thrid dimension who gets infected)
inf_trans <- rep(list(rep(list(NA),n)),Tmax)

# function determing who gets infected:
trans_func <- function(ine, p, n, m_probs){
  ntrans <- rbinom(1, ine, p)
  if(ntrans == 0){return(NA)} else{return(sample(1:n, ntrans, prob = m_probs))}
}


# Start infection rounds:
for(t in 2:Tmax){

# previous round infected:
it_1_ind <- which(inf_dat[,t-1]>0)

if(length(it_1_ind)>0){

# determine all transmission:
if(t<= t_lockdown){
  inf_trans[[t]][it_1_ind] <- lapply(it_1_ind, function(x) trans_func(ine[[x]],p,n, m_probs = cluster_matrix_norm[x,]))
}else{
  inf_trans[[t]][it_1_ind] <- lapply(it_1_ind, function(x) trans_func(ine_2[[x]],p,n, m_probs = cluster_matrix_norm[x,]))
}



#continue previous round infected:
inf_dat[,t] <- inf_dat[,t-1] 
inf_dat[it_1_ind,t] <- inf_dat[it_1_ind,t]+1

#set newly infected persons in current round (only if not NA or positive yet):
indis <- unlist(inf_trans[[t]][it_1_ind])
if(any(is.finite(indis))){
  #Check if already NA or infected:
  indi_sub <- which(inf_dat[indis,t]==0)
  if(length(indi_sub)>0){
    inf_dat[indis[indi_sub],t] <- 1
  }
}

#set people who are healthy again to NA
indis2 <- which(inf_dat[,t]>Tinf)
if(any(is.finite(indis2))){inf_dat[indis2,t:Tmax] <- NA}
}
# else{
#   inf_dat[,t] <- inf_dat[,t]
# }

}

simu_res[[j]] <- inf_dat
pb$tick()

}

###### Evaluate ######

## Across simulations:

# Total infected:
tot_inf <- unlist(lapply(simu_res, function(x) sum(rowSums(x, na.rm = T)>0)))

mean(tot_inf)
sd(tot_inf)

hist(tot_inf)

# Simultaneaous infected:
sim_inf <- simplify2array(lapply(simu_res, function(x) colSums(x>0, na.rm = T)))
dim(sim_inf)

sim_inf_mean <- rowMeans(sim_inf)
sim_inf_quant <- apply(sim_inf, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))


plot(sim_inf_quant[1,], type = 'l', ylim = c(0,n), lty = 2, xlim = c(0,200))
lines(sim_inf_quant[2,], type = 'l', ylim = c(0,n))
lines(sim_inf_quant[3,], type = 'l', ylim = c(0,n), lty = 2)
abline(v=t_lockdown, col = 2, lty = 2)

matplot(sim_inf, type = 'l', col = add.alpha('blue', 0.3), lty = 1, ylim = c(0,n))
lines(sim_inf_mean, lwd = 2, col =2)

abline(v=t_lockdown, col = 2, lty = 2)

# Cumulative infected:
cum_inf <- simplify2array(lapply(simu_res, function(x) colSums((is.na(x)|x>0), na.rm = T)))

cum_inf_mean <- rowMeans(cum_inf)
cum_inf_quant <- apply(cum_inf, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))

plot(cum_inf_quant[2,], type = 'l', ylim = c(0,n), lty = 1, xlab = 'Days', ylab = 'Cumulative infected')
grid()
polygon(c(1:Tmax, Tmax:1), c(cum_inf_quant[1,],rev(cum_inf_quant[3,])), border = NA, col = add.alpha('gray', 0.5))
lines(cum_inf_quant[2,], type = 'l', ylim = c(0,n), lty = 1, lwd = 2)
abline(v=t_lockdown, col = 2, lty = 2)

matplot(cum_inf,type = 'l', col = add.alpha('blue', 0.3), lty = 1, ylim = c(0, n))
grid()
lines(cum_inf_mean, lwd = 2, col = 2)
abline(v=t_lockdown, col = 2, lty = 2)

plot(cum_inf_mean/n, lwd = 2, col = 2, type = 'l')
lines(for_daddy_1000/1000, lwd = 2, col = 1, type = 'l', lty = 2)

cum_inf_mean/n
for_daddy_1000/1000

#Mean infection rate:
plot(2:Tmax, cum_inf_mean[2:Tmax]/cum_inf_mean[1:(Tmax-1)], xlab = 'Day', ylab = 'Infection rate', type = 'h')


# Newly infected:
new_inf <- simplify2array(lapply(simu_res, function(x) colSums(x ==1, na.rm = T)))

new_inf_mean <- rowMeans(new_inf)
new_inf_quant <- apply(new_inf, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))

plot(new_inf_quant[2,], type = 'l', ylim = c(0,max(new_inf_quant[3,])), lty = 1, xlab = 'Days', ylab = 'newly infected')
grid()
polygon(c(1:Tmax, Tmax:1), c(new_inf_quant[1,],rev(new_inf_quant[3,])), border = NA, col = add.alpha('gray', 0.5))
lines(new_inf_quant[2,], type = 'l', ylim = c(0,n), lty = 1, lwd = 2)

matplot(new_inf, type = 'l', col = add.alpha('blue', 0.3), lty = 1)
lines(new_inf_mean, lwd = 2, col =2)




plot(colSums(inf_dat==1, na.rm = T), type = 'h')
#Immune:
plot(colSums(is.na(inf_dat), na.rm = T), type = 's')
# Total infected:
sum(rowSums(inf_dat, na.rm = T)>0)



