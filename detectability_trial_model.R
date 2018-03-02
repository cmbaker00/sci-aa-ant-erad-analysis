library(gdata)
library(greta)
library(bayesplot)


this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
getwd()

detect_data = as.matrix(read.csv("detect_trial.csv"))
detect_data = detect_data[rowSums(detect_data,na.rm = TRUE)>0,] # only include sites where there was a detection

persist_data = read.csv("Persistance_data.csv")
pre_treat = persist_data$pre.treatment #sites with ants before treatment
post_treat = persist_data$post.2.treatments #sites with ants after treatment

treat_data_sites = (pre_treat == 0 | pre_treat == 1) & (post_treat == 0 | post_treat == 1) # finding sites with data pre and post treatment

# only keep data when there is pre and post data
pre_treat = pre_treat[treat_data_sites]
post_treat = post_treat[treat_data_sites]

# only keep data where ants were detected pre treatment
post_treat = post_treat[pre_treat==1]
pre_treat = pre_treat[pre_treat==1]

# Data as great arrays
Y = as_data(detect_data[is.nan(detect_data)==FALSE])
Z = as_data(as.numeric(as.matrix(post_treat)))

# Define variables
p = uniform(0,1)
psi = uniform(0,1)

# Define distributions for the data
distribution(Y) = binomial(1,p)
distribution(Z) = binomial(1,p*psi^2)

# Create the greta model
m <- model(p,psi)

# Plot the model
plot(m)

# Run 3 MCMC chains
draws_list <- replicate(3, mcmc(m))
draws <- coda::mcmc.list(draws_list)

# Check convergence
sirt::mcmc.list.descriptives(draws)[, "Rhat"]

# Plot trace and intervals
mcmc_trace(draws)
mcmc_intervals(draws)
s <- summary(draws)
s
