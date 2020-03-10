# Extended Oscillations Function Source
# By Hannah De los Santos
# Code description: Contains all the funcitons for extended harmonic oscillator work, in order to have less confusion between scripts.

# combined models ----

echo_lin_const <- function(a1,gam1,omega1,phi1,slo1,y_shift1,y_shift2,t, s1, s2){
  return(((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+slo1*t+y_shift1)*s1)+(y_shift2*s2))
}

echo_lin_lin <- function(a1,gam1,omega1,phi1,slo1,y_shift1,slo2,y_shift2,t, s1, s2){
  return(((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+slo1*t+y_shift1)*s1)+((slo2*t+y_shift2)*s2))
}

echo_lin_exp <- function(a1,gam1,omega1,phi1,slo1,y_shift1,a2,gamma2,y_shift2,t, s1, s2){
  return(((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+slo1*t+y_shift1)*s1)+((a2*exp(gamma2*t)+y_shift2)*s2))
}

echo_lin_exp_lin <- function(a1,gam1,omega1,phi1,slo1,y_shift1,a2,gamma2,slo2,y_shift2,t, s1, s2){
  return(((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+slo1*t+y_shift1)*s1)+((a2*exp(gamma2*t)+slo2*t+y_shift2)*s2))
}

echo_lin_echo <- function(a1,gam1,omega1,phi1,slo1,y_shift1,a2,gam2,omega2,phi2,y_shift2,t, s1, s2){
  return(((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+slo1*t+y_shift1)*s1)+((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+y_shift1)*s1))
}

echo_const <- function(a1,gam1,omega1,phi1,y_shift1,y_shift2,t, s1, s2){
  return(((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+y_shift1)*s1)+(y_shift2*s2))
}

echo_lin <- function(a1,gam1,omega1,phi1,y_shift1,slo2,y_shift2,t, s1, s2){
  return(((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+y_shift1)*s1)+((slo2*t+y_shift2)*s2))
}

echo_exp <- function(a1,gam1,omega1,phi1,y_shift1,a2,gamma2,y_shift2,t, s1, s2){
  return(((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+y_shift1)*s1)+((a2*exp(gamma2*t)+y_shift2)*s2))
}

echo_exp_lin <- function(a1,gam1,omega1,phi1,y_shift1,a2,gamma2,slo2,y_shift2,t, s1, s2){
  return(((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+y_shift1)*s1)+((a2*exp(gamma2*t)+slo2*t+y_shift2)*s2))
}

echo_echo <- function(a1,gam1,omega1,phi1,y_shift1,a2,gam2,omega2,phi2,y_shift2,t, s1, s2){
  return(((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+y_shift1)*s1)+((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+y_shift1)*s1))
}

exp_lin_const <- function(a1,gamma1,slo1,y_shift1,y_shift2,t, s1, s2){
  return(((a1*exp(gamma1*t)+slo1*t+y_shift1)*s1)+(y_shift2*s2))
}

exp_lin_lin <- function(a1,gamma1,slo1,y_shift1,slo2,y_shift2,t, s1, s2){
  return(((a1*exp(gamma1*t)+slo1*t+y_shift1)*s1)+((slo2*t+y_shift2)*s2))
}

exp_lin_exp <- function(a1,gamma1,slo1,y_shift1,a2,gamma2,y_shift2,t, s1, s2){
  return(((a1*exp(gamma1*t)+slo1*t+y_shift1)*s1)+((a2*exp(gamma2*t)+y_shift2)*s2))
}

exp_lin_exp_lin <- function(a1,gamma1,slo1,y_shift1,a2,gamma2,slo2,y_shift2,t, s1, s2){
  return(((a1*exp(gamma1*t)+slo1*t+y_shift1)*s1)+((a2*exp(gamma2*t)+slo2*t+y_shift2)*s2))
}

exp_const <- function(a1,gamma1,y_shift1,y_shift2,t, s1, s2){
  return(((a1*exp(gamma1*t)+y_shift1)*s1)+(y_shift2*s2))
}

exp_lin <- function(a1,gamma1,y_shift1,slo2,y_shift2,t, s1, s2){
  return(((a1*exp(gamma1*t)+y_shift1)*s1)+((slo2*t+y_shift2)*s2))
}

exp_exp <- function(a1,gamma1,y_shift1,a2,gamma2,y_shift2,t, s1, s2){
  return(((a1*exp(gamma1*t)+y_shift1)*s1)+((a2*exp(gamma2*t)+y_shift2)*s2))
}

lin_const <- function(slo1,y_shift1,y_shift2,t, s1, s2){
  return(((slo1*t+y_shift1)*s1)+(y_shift2*s2))
}

lin_lin <- function(slo1,y_shift1,slo2,y_shift2,t, s1, s2){
  return(((slo1*t+y_shift1)*s1)+((slo2*t+y_shift2)*s2))
}

const_const <- function(y_shift1,y_shift2,t, s1, s2){
  return((y_shift1*s1)+(y_shift2*s2))
}

# other functions ----

# function to represent damped oscillator with phase and equilibrium shift formula
# inputs:
#  a: Amplitude
#  gam: Amplitude.Change.Coefficient (amount of damping/driving)
#  omega: Radial frequency
#  phi: Phase Shift (radians)
#  y_shift: Equilibrium shift
#  t: time
# outputs:
#  result of inputs into formula
alt_form <- function(a,gam,omega,phi,y_shift,t){
  return (a*exp(-1*gam*(t)/2)*cos((omega*t)+phi)+y_shift)
}

# echo with linear trend
echo_lin_mod <- function(a,gam,omega,phi,slo,y_shift,t){
  return (a*exp(-1*gam*(t)/2)*cos((omega*t)+phi)+slo*t+y_shift)
}

# constant model
const_mod <- function(y_shift,t){
  return(rep(y_shift, length(t)))
}

# linear model
lin_mod <- function(slo,y_shift,t){
  return(slo*t+y_shift)
}

# exponential model with linear trend
exp_lin_mod <- function(a,gamma,slo,y_shift,t){
  return(a*exp(gamma*t)+slo*t+y_shift)
}

# exponential model
exp_mod <- function(a,gamma,y_shift,t){
  return(a*exp(gamma*t)+y_shift)
}

echo_lin_joint <- function(a1,a2, gam1,gam2, omega1, phi1,phi2, slo1,slo2, y_shift1,y_shift2, t, s1, s2){
  return ((a1*s1+a2*s2)*exp(-1*(gam1*s1+gam2*s1)*(t)/2)*cos((omega1*t)+(phi1*s1+phi2*s2))+((slo1*s1 + slo2*s2)*t)+(y_shift1*s1+y_shift2*s2))
}

echo_joint <- function(a1,a2, gam1,gam2, omega1, phi1,phi2, y_shift1,y_shift2, t, s1, s2){
  return ((a1*s1+a2*s2)*exp(-1*(gam1*s1+gam2*s2)*(t)/2)*cos((omega1*t)+(phi1*s1+phi2*s2))+(y_shift1*s1+y_shift2*s2))
}

# omega and chang are in hours
echo_joint_res_constr <- function(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, y_shift1,y_shift2, t, s1, s2, resol){
  # time_constr <- 0
  # if (abs((2*pi/omega1)-(2*pi/omega2))>resol[1]){
  #   time_constr <- Inf
  # }
  
  # return ((a1*s1+a2*s2)*exp(-1*(gam1*s1+gam2*s2)*(t)/2)*cos(((omega1*s1+((omega1+chang)*s2))*t)+(phi1*s1+phi2*s2))+(y_shift1*s1+y_shift2*s2))
  return ((a1*s1+a2*s2)*exp(-1*(gam1*s1+gam2*s2)*(t)/2)*cos(((2*pi/(omega1*s1+((omega1+chang)*s2)))*t)+(phi1*s1+phi2*s2))+(y_shift1*s1+y_shift2*s2))
}

# omega and chang are in hours
echo_lin_joint_res_constr <- function(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, slo1,slo2, y_shift1,y_shift2, t, s1, s2, resol){
  # time_constr <- 0
  # if (abs((2*pi/omega1)-(2*pi/omega2))>resol[1]){
  #   time_constr <- Inf
  # }
  
  # return ((a1*s1+a2*s2)*exp(-1*(gam1*s1+gam2*s2)*(t)/2)*cos(((omega1*s1+((omega1+chang)*s2))*t)+(phi1*s1+phi2*s2))+(y_shift1*s1+y_shift2*s2))
  return ((a1*s1+a2*s2)*exp(-1*(gam1*s1+gam2*s2)*(t)/2)*cos(((2*pi/(omega1*s1+((omega1+chang)*s2)))*t)+(phi1*s1+phi2*s2))+((slo1*s1 + slo2*s2)*t)+(y_shift1*s1+y_shift2*s2))
}

jac_echo_joint <- function(a1,a2, gam1,gam2, omega1, phi1,phi2, y_shift1,y_shift2, t, s1, s2){
  echo_expr <- expression((a1*s1+a2*s2)*exp(-1*(gam1*s1+gam2*s2)*(t)/2)*cos((omega1*t)+(phi1*s1+phi2*s2))+(y_shift1*s1+y_shift2*s2))
  
  # needs to be in the same order as the initial guess
  jac_val <- c(
    eval(D(echo_expr, "gam1")),
    eval(D(echo_expr, "a1")),
    eval(D(echo_expr, "omega1")),
    eval(D(echo_expr, "phi1")),
    eval(D(echo_expr, "y_shift1")),
    eval(D(echo_expr, "gam2")),
    eval(D(echo_expr, "a2")),
    eval(D(echo_expr, "phi2")),
    eval(D(echo_expr, "y_shift2"))
  )
  
  return(jac_val)
}

jac_exp_lin_joint <- function(a,gam,slo,y_shift,t){
  exp_lin_expr <- expression(a*exp(gam*t)+slo*t+y_shift)
  
  # needs to be in the same order as the initial guess
  jac_val <- c(
    eval(D(exp_lin_expr, "a")),
    eval(D(exp_lin_expr, "gam")),
    eval(D(exp_lin_expr, "slo")),
    # eval(D(exp_lin_expr, "y_shift"))
    rep(1, length(t))
  )
  
  return(jac_val)
}


fcn <- function(p, t, s1, s2, w, y, fcall, jcall){
  sqrt(w)*(y - do.call("fcall", c(list(t = t, s1 =s1, s2 = s2), as.list(p))))
}

fcn_one_rep <- function(p, t, s1, s2, y, fcall, jcall){
  (y - do.call("fcall", c(list(t = t, s1 =s1, s2 = s2), as.list(p))))
}

# fcn_constr <- function(p, t, s1, s2, w, y, fcall, jcall){
#   sqrt(w)*(y - do.call("fcall", c(list(t = t, s1 =s1, s2 = s2, resol=resol), as.list(p))))
# }

# fcn.jac <- function(p, t, s1, s2, w, y, fcall, jcall){
#   -do.call("jcall", c(list(t = t, s1 =s1, s2 = s2), as.list(p)))
# } 

fcn_single <- function(p, t, w, y, fcall, jcall){
  sqrt(w)*(y - do.call("fcall", c(list(t = t), as.list(p))))
}

fcn_single_one_rep <- function(p, t, y, fcall, jcall){
  (y - do.call("fcall", c(list(t = t), as.list(p))))
}


# fcn.jac_else <- function(p, t, w, y, fcall, jcall){
#   -do.call("jcall", c(list(t = t), as.list(p)))
# }

# function to calculate the average of replicates
# inputs:
#  current_gene: row number of gene being examined
#  num_reps: number of replicates
# outputs:
#  a vector of means of gene expressions for replicates
avg_rep <- function(current_gene, num_reps){
  return(rbind(sapply(seq(2,ncol(genes), by = num_reps), function(x) mean(unlist(genes[current_gene,c(x:(num_reps-1+x))]), na.rm = TRUE))))
}

# function to calculate the average of all replicates. Used primarily for cases with multiple replicates.
# inputs:
#  num_reps: number of replicates
# outputs:
#  a matrix of means of gene expressions for replicates
avg_all_rep <- function(num_reps, genes, timen){
  # originally based on heat map code, but it will work fine here
  
  #get matrix of just the relative expression over time
  hm_mat <- as.matrix(genes[,2:ncol(genes)])
  
  #if there are replicates, average the relative expression for each replicate
  mtx_reps <- list() # to store actual replicate matrix
  mtx_count <- list() # to store how many are NA
  for (i in 1:num_reps){
    mtx_reps[[i]] <- hm_mat[, seq(i,ncol(hm_mat), by=num_reps)]
    mtx_count[[i]] <- is.na(mtx_reps[[i]])
    mtx_reps[[i]][is.na(mtx_reps[[i]])] <- 0
  }
  repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat))+num_reps # to store how many we should divide by
  hm_mat <- matrix(0L,ncol = length(timen),nrow = nrow(hm_mat)) # to store the final result
  for (i in 1:num_reps){
    hm_mat <- hm_mat + mtx_reps[[i]] # sum the replicates
    repmtx <- repmtx - mtx_count[[i]] # how many replicates are available for each time point
  }
  repmtx[repmtx==0] <- NA # to avoid division by 0 and induce NAs if there are no time points available
  hm_mat <- hm_mat/repmtx
  
  return(hm_mat)
}

# function to calculate the variance of replicates at a certain time point
# inputs:
#  x: time point
#  current_gene: row number of gene being examined
#  num_reps: number of replicates
# outputs:
#  the variance of replicates at a certain time point
calc_var <- function(x, current_gene, num_reps, genes){
  eps <- 1e-7 # slight adjustment for 0 variance case
  std2 <- var(unlist(genes[current_gene,c(x:(num_reps-1+x))]), na.rm = TRUE) # calc variance
  if (is.na(std2)){ # possible bug for NA
    std2 <- 1
  }
  return(1/(std2+eps))
}

# function to calculate the weights for replicate fitting (these weights are the inverse of the variance at each time point)
# inputs:
#  current_gene: row number of gene being examined
#  num_reps: number of replicates
# outputs:
#  the variance of replicates at a certain time point
calc_weights <- function(current_gene, num_reps, genes){
  return(sapply(seq(2,ncol(genes), by = num_reps), function(x) calc_var(x,current_gene,num_reps, genes)))
}

# function to find confidence intervals by the standard jackknifing method (leave one out)
# inputs:
#  temp: data frame with time points, y values, and weights (if used)
#  parameters: parameters found by full nlsLM fit
#  num_reps: number of replicates
#  start_param: list of starting parameters for fits
# outputs:
#  confidence intervals for each parameter, low and then high
jackknife <- function(temp, parameters, num_reps, start_param){
  # preallocate spaces
  all_fits <- list()
  all_pars <- vector(mode = "list", 5)
  n <- nrow(temp)
  
  # go through and compute fit for each leave one out
  for (i in 1:nrow(temp)){
    temp_edit <- temp[-i,]
    
    if (num_reps==1){
      all_fits[[i]] <- suppressWarnings(nlsLM(y ~ alt_form(a,gam,omega,phi,y_shift,t),
                                              data=temp_edit,
                                              start=, start_param,
                                              lower=c(-Inf, -Inf, high, -Inf, min(temp_edit$y)),
                                              upper=c(Inf, Inf, low, Inf, max(temp_edit$y)),
                                              control = nls.lm.control(maxiter = 1000, maxfev = 2000,
                                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)))
    } else {
      #fitting
      all_fits[[i]] <- suppressWarnings(nlsLM(y ~ alt_form(a,gam,omega,phi,y_shift,t),
                                              data=temp_edit,
                                              start=, start_param,
                                              lower=c(-Inf, -Inf, high, -Inf, min(temp_edit$y)),
                                              upper=c(Inf, Inf, low, Inf, max(temp_edit$y)),
                                              control = nls.lm.control(maxiter = 1000, maxfev = 2000,
                                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6),
                                              weights = w))
    }
    
    # store all parameters
    for (p in 1:5){
      all_pars[[p]] <- c(all_pars[[p]],all_fits[[i]]$m$getAllPars()[p])
    }
  }
  
  # compute jackknifed confidence intervals
  par_names <- c("gam","a","omega","phi","y_shift")
  ci_names <- c(paste0("ci_low_",par_names), paste0("ci_high_",par_names))
  ci_int <- 1:(2*5)
  for (p in 1:5){
    ps <- (n*parameters[p])-((n-1)*all_pars[[p]])
    v <- 1/n/(n-1)*sum((ps-(sum(ps)/n))^2)
    ci_int[p] <- mean(ps) - qt(0.975, n-1)*sqrt(var(ps)/n)
    ci_int[p+5] <- mean(ps) + qt(0.975, n-1)*sqrt(var(ps)/n)
  }
  names(ci_int) <- ci_names
  
  return(ci_int)
}


# function to calculate the variance of replicates at a certain time point for bootstrapping
# inputs:
#  x: time point
#  boot_gene: gene time course being examined
#  num_reps: number of replicates
# outputs:
#  the variance of replicates at a certain time point
calc_var_boot <- function(x, boot_gene, num_reps){
  eps <- 1e-7 # slight adjustment for 0 variance case
  std2 <- var(unlist(boot_gene[c(x:(num_reps+x))]), na.rm = TRUE) # calc variance
  if (is.na(std2)){ # possible bug for NA
    std2 <- 1
  }
  return(1/(std2+eps))
}

# function to calculate the weights for replicate fitting for bootstrapping (these weights are the inverse of the variance at each time point)
# inputs:
#  boot_gene: gene time course being examined
#  num_reps: number of replicates
# outputs:
#  the variance of replicates at a certain time point
calc_weights_boot <- function(boot_gene, num_reps){
  return(sapply(seq(1,length(boot_gene), by = num_reps), function(x) calc_var_boot(x,boot_gene,num_reps)))
}

# function to find confidence intervals by the standard jackknifing method (leave one out)
# inputs:
#  temp: data frame with time points, y values, and weights (if used)
#  fit: original fit data
#  parameters: parameters found by full nlsLM fit
#  num_reps: number of replicates
#  start_param: list of starting parameters for fits
#  current_gene: row number of gene being examined
#  seed: random seed to set for reproducibility
# outputs:
#  confidence intervals for each parameter, low and then high
bootstrap <- function(temp, fit, start_param, num_reps, current_gene, seed){
  #this copy should be avoided if you have big data, but I don't have time right now:
  df <- temp
  df$fitted <- fitted(fit)
  df$resid <- residuals(fit)
  fun <- function(df, inds) {
    # get the new bootstrapped values
    df$bootGPP <- df$fitted + df$resid[inds]
    
    
    tryCatch({
      if (num_reps > 1){
        # get new weights based on these
        w <- rep(calc_weights_boot(df$bootGPP, num_reps), each = num_reps)
        df$w <- w[!is.na(genes[current_gene,-1])]
        
        suppressWarnings(coef(nlsLM(bootGPP ~ alt_form(a,gam,omega,phi,y_shift,t),
                                    data=df,
                                    start=start_param,
                                    lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
                                    upper=c(Inf, Inf, low, Inf, max(temp$y)),
                                    control = nls.lm.control(maxiter = 1000, maxfev = 2000,
                                                             ftol=1e-6, ptol=1e-6, gtol=1e-6),
                                    weights = w)))
      } else {
        suppressWarnings(coef(nlsLM(bootGPP ~ alt_form(a,gam,omega,phi,y_shift,t),
                                    data=df,
                                    start=start_param,
                                    lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
                                    upper=c(Inf, Inf, low, Inf, max(temp$y)),
                                    control = nls.lm.control(maxiter = 1000, maxfev = 2000,
                                                             ftol=1e-6, ptol=1e-6, gtol=1e-6))))
      }
    }, error = function(e){
      return(unlist(start_param))
    })
  }
  
  set.seed(seed)
  b <- boot(df, fun, R = 999)
  # get the 95% confidence intervals
  par_names <- c("gam","a","omega","phi","y_shift")
  ci_names <- c(paste0("ci_low_",par_names), paste0("ci_high_",par_names))
  ci_int <- 1:(2*5)
  for (p in 1:5){
    bci <- boot.ci(b, index=p, type = "perc")
    
    if (!is.null(bci)){
      ci_int[p] <- bci$percent[4]
      ci_int[p+5] <- bci$percent[5]
    } else { # if there's no variation in the parameter - fixed interval
      ci_int[p] <- ci_int[p+5] <- b$t0[3]
    }
  }
  names(ci_int) <- ci_names
  
  
  return(ci_int)
}

# Function to calculate the parameters for the extended harmonic oscillator equation for a specific gene.
#
#  current_gene: row number of current gene we want to calculate parameters for
#  timen: time points for dataset
#  resol: resolution of time points
#  num_reps: number of replicates
#  tied: if replicate data, whether the replicates are related (paired) or not (unpaired)
#  is_smooth: boolean that indicates whether data should be smoothed or not
#  is_weighted: if there is smoothing, is it weighted (1,2,1) smoothing, or unweighted smoothing (1,1,1)
#  low: the highest frequency we are looking for, in radians (lowest period)
#  high: the lowest frequency we are looking for, in radians (highest period)
#  rem_unexpr: logical indicating whether genes with less than rem_unexpr_amt % expression should not be considered
#  rem_unexpr_amt: percentage of expression for which genes should not be considered
#  run_conf: boolean of whether or not to run confidence intervals
#  which_conf: string of which type of confidence interval to compute
#  harm_cut: postive number indicating the cutoff for a gene to be considered harmonic
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
#  seed: number for random seed to fix for bootstrapping for confidence intervals
#  results, a data frame which contains:
#   gene: gene name
#   conv: did the fit converge, or descriptor of type of data (constant, unexpressed, etc.)
#   iter: number of iterations
#   gamma: forcing coefficient value for fit
#   type_gam: Type of oscillation (damped, forced, etc.)
#   amplitude: Amplitude value for fit
#   omega: Radial frequency for fit
#   period: Period for fit (in time units)
#   phase.shift: Phase shift for fit (radians)
#   hours.shift: Phase shift for fit (hours)
#   tau: Kendall's tau between original and fitted values
#   y_shift: Equilibrium shift for fit
#   pval: P-value calculated based on Kendall's tau
#   (ci_int: confidence for each parameter)
#   original.values: original values for gene
#   fitted.values: fitted values for gene
calculate_param <- function(current_gene,timen,resol,num_reps,tied,is_smooth=FALSE,is_weighted=FALSE,low,high,rem_unexpr=FALSE,rem_unexpr_amt=70,run_conf = F, which_conf = "Bootstrap", harm_cut = .03, over_cut = .15, seed = 30){
  
  if (run_conf){
    ci_int <- rep(NA, 10)
    par_names <- c("gam","a","omega","phi","y_shift")
    ci_names <- c(paste0("ci_low_",par_names), paste0("ci_high_",par_names))
    names(ci_int) <- ci_names
  }
  
  gene_n <- as.character(genes[current_gene,1]) # gene name
  # first we need to check whether or not the gene is just a straight line
  if (!is_deviating(current_gene)){ # one replicate
    results <- data.frame(gene = gene_n, conv = "No Deviation", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA,  tau = NA, pval = NA, stringsAsFactors = FALSE)
    
    if (run_conf){
      results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
    } else {
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
    }
    return(results)
  }
  
  # then we need to check if < threshold % are expressed (if desired)
  if (rem_unexpr){
    if (rem_unexpr_vect[current_gene]){
      results <- data.frame(gene = gene_n, conv = "Unexpressed", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
      if (run_conf){
        results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
      } else {
        results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
      }
      return(results)
    }
  }
  
  # calculating averages for initial values
  if (num_reps == 1){
    y_val <- rbind(as.numeric(as.character(t(genes[current_gene,c(2:ncol(genes))])))) # all the gene values
  } else{
    y_val <- avg_genes[current_gene,] # starting values determined by average of replicates
  }
  
  tryCatch({ # throw exception upon error
    # calculate the amount of peaks
    peaks <- c(); # vector of peak values
    peaks_time <- c(); # vector of peak timen
    counting <- 1; # counter
    {
      if (resol <= ((1/12)+10^-8)){ # 17 hour surround
        mod <- 102
        for(i in (mod+1):(length(y_val)-mod)){
          # deal with complete missingness
          if (suppressWarnings(max(y_val[(i-mod):(i+mod)], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
            next
          }
          # otherwise continue as normal
          if (y_val[i] == max(y_val[(i-mod):(i+mod)], na.rm = TRUE)){
            peaks[counting] <- y_val[i]
            peaks_time[counting] <- timen[i]
            counting <- counting+1
          }
        }
      } else if (resol <= ((1/6)+10^-8)){ # 15 hour surround
        mod <- 45
        for(i in (mod+1):(length(y_val)-mod)){
          # deal with complete missingness
          if (suppressWarnings(max(y_val[(i-mod):(i+mod)], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
            next
          }
          # otherwise continue as normal
          if (y_val[i] == max(y_val[(i-mod):(i+mod)], na.rm = TRUE)){
            peaks[counting] <- y_val[i]
            peaks_time[counting] <- timen[i]
            counting <- counting+1
          }
        }
      } else if (resol <= ((1/4)+10^-8)){ # 13 hour surround
        mod <- 26
        for(i in (mod+1):(length(y_val)-mod)){
          # deal with complete missingness
          if (suppressWarnings(max(y_val[(i-mod):(i+mod)], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
            next
          }
          # otherwise continue as normal
          if (y_val[i] == max(y_val[(i-mod):(i+mod)], na.rm = TRUE)){
            peaks[counting] <- y_val[i]
            peaks_time[counting] <- timen[i]
            counting <- counting+1
          }
        }
      } else if (resol <= ((1/2)+10^-8)){ # 11 hour surround
        mod <- 11
        for(i in (mod+1):(length(y_val)-mod)){
          # deal with complete missingness
          if (suppressWarnings(max(y_val[(i-mod):(i+mod)], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
            next
          }
          # otherwise continue as normal
          if (y_val[i] == max(y_val[(i-mod):(i+mod)], na.rm = TRUE)){
            peaks[counting] <- y_val[i]
            peaks_time[counting] <- timen[i]
            counting <- counting+1
          }
        }
      } else if (resol <= 1){
        # go through gene values and find maximum as compared to 8 surrounding values
        # finding peaks for first 4 points
        # deal with complete missingness
        # if (suppressWarnings(max(y_val[1:9], na.rm = TRUE)) != -Inf){
        #   for (i in 1:4){
        #     # otherwise continue as normal
        #     if (y_val[i] == max(y_val[1:9], na.rm = TRUE)){
        #       peaks[counting] <- y_val[i]
        #       peaks_time[counting] <- timen[i]
        #       counting <- counting+1
        #     }
        #   }
        # }
        
        for(i in 5:(length(y_val)-4)){
          # deal with complete missingness
          if (suppressWarnings(max(y_val[i-4],y_val[i-3],y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],y_val[i+3],y_val[i+4], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
            next
          }
          # otherwise continue as normal
          if (y_val[i] == max(y_val[i-4],y_val[i-3],y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],y_val[i+3],y_val[i+4], na.rm = TRUE)){
            peaks[counting] <- y_val[i]
            peaks_time[counting] <- timen[i]
            counting <- counting+1
          }
        }
        
        # finding peaks for last 4 points
        # deal with complete missingness
        # if (suppressWarnings(max(y_val[(length(y_val)-8):length(y_val)], na.rm = TRUE)) != -Inf){
        #   for (i in (length(y_val)-3):length(y_val)){
        #     # otherwise continue as normal
        #     if (y_val[i] == max(y_val[(length(y_val)-8):length(y_val)], na.rm = TRUE)){
        #       peaks[counting] <- y_val[i]
        #       peaks_time[counting] <- timen[i]
        #       counting <- counting+1
        #     }
        #   }
        # }
      } else if (resol <=2){
        # go through gene values and find maximum as compared to six surrounding values
        # finding peaks for first 3 points
        # deal with complete missingness
        # if (suppressWarnings(max(y_val[1:7], na.rm = TRUE)) != -Inf){
        #   for (i in 1:3){
        #     # otherwise continue as normal
        #     if (y_val[i] == max(y_val[1:7], na.rm = TRUE)){
        #       peaks[counting] <- y_val[i]
        #       peaks_time[counting] <- timen[i]
        #       counting <- counting+1
        #     }
        #   }
        # }
        for(i in 4:(length(y_val)-3)){
          # deal with complete missingness
          if (suppressWarnings(max(y_val[i-3],y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],y_val[i+3], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
            next
          }
          # otherwise continue as normal
          if (y_val[i] == max(y_val[i-3],y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],y_val[i+3], na.rm = TRUE)){
            peaks[counting] <- y_val[i]
            peaks_time[counting] <- timen[i]
            counting <- counting+1
          }
        }
        # finding peaks for last 3 points
        # deal with complete missingness
        # if (suppressWarnings(max(y_val[(length(y_val)-6):length(y_val)], na.rm = TRUE)) != -Inf){
        #   for (i in (length(y_val)-2):length(y_val)){
        #     # otherwise continue as normal
        #     if (y_val[i] == max(y_val[(length(y_val)-6):length(y_val)], na.rm = TRUE)){
        #       peaks[counting] <- y_val[i]
        #       peaks_time[counting] <- timen[i]
        #       counting <- counting+1
        #     }
        #   }
        # }
      } else if (resol <= 4){
        # finding peaks for first 2 points
        # deal with complete missingness
        # if (suppressWarnings(max(y_val[1:5], na.rm = TRUE)) != -Inf){
        #   for (i in 1:2){
        #     # otherwise continue as normal
        #     if (y_val[i] == max(y_val[1:5], na.rm = TRUE)){
        #       peaks[counting] <- y_val[i]
        #       peaks_time[counting] <- timen[i]
        #       counting <- counting+1
        #     }
        #   }
        # }
        # go through gene values and find maximum as compared to four surrounding values
        for(i in 3:(length(y_val)-2)){
          # to deal with complete missingness
          if(suppressWarnings(max(y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],na.rm = TRUE))== -Inf  | is.na(y_val[i])){
            next
          }
          if (y_val[i] == max(y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],na.rm = TRUE)){
            peaks[counting] <- y_val[i]
            peaks_time[counting] <- timen[i]
            counting <- counting+1
          }
        }
        
        # finding peaks for last 3 points
        # deal with complete missingness
        # if (suppressWarnings(max(y_val[(length(y_val)-4):length(y_val)], na.rm = TRUE)) != -Inf){
        #   for (i in (length(y_val)-1):length(y_val)){
        #     # otherwise continue as normal
        #     if (y_val[i] == max(y_val[(length(y_val)-4):length(y_val)], na.rm = TRUE)){
        #       peaks[counting] <- y_val[i]
        #       peaks_time[counting] <- timen[i]
        #       counting <- counting+1
        #     }
        #   }
        # }
        
      } else{
        # finding peaks for first point
        # deal with complete missingness
        # if (suppressWarnings(max(y_val[1:3], na.rm = TRUE)) != -Inf){
        #   for (i in 1){
        #     # otherwise continue as normal
        #     if (y_val[i] == max(y_val[1:3], na.rm = TRUE)){
        #       peaks[counting] <- y_val[i]
        #       peaks_time[counting] <- timen[i]
        #       counting <- counting+1
        #     }
        #   }
        # }
        
        # go through gene values and find maximum as compared to two surrounding values
        for(i in 2:(length(y_val)-1)){
          # to deal with complete missingness
          if(suppressWarnings(max(y_val[i-2],y_val[i-1],y_val[i],y_val[i+1],y_val[i+2],na.rm = TRUE))== -Inf  | is.na(y_val[i])){
            next
          }
          if (y_val[i] == max(y_val[i-1],y_val[i],y_val[i+1], na.rm = TRUE)){
            peaks[counting] <- y_val[i]
            peaks_time[counting] <- timen[i]
            counting <- counting+1
          }
        }
        
        # finding peaks for last 3 points
        # deal with complete missingness
        # if (suppressWarnings(max(y_val[(length(y_val)-2):length(y_val)], na.rm = TRUE)) != -Inf){
        #   for (i in length(y_val)){
        #     # otherwise continue as normal
        #     if (y_val[i] == max(y_val[(length(y_val)-2):length(y_val)], na.rm = TRUE)){
        #       peaks[counting] <- y_val[i]
        #       peaks_time[counting] <- timen[i]
        #       counting <- counting+1
        #     }
        #   }
        # }
      }
    }
    # calculate starting amplitude, y_shift
    y0 <- mean(y_val,na.rm = TRUE) # intial value for the equilibrium shift
    if (y0 < 10^-10 && y0 > -10^-10){
      y0 <- 10^-8 # avoiding 0 mean, which makes gradient singular
    }
    x0 <- min(timen) # the x start parameter
    a0 <- max(y_val,na.rm = TRUE) - y0 # mean(y_val) # initial guess for amplitude
    
    # intial value for gamma
    if (length(peaks)==0){ # if there are no peaks, we account for that
      gam0 <- 0
    } else if (which.max(peaks)==1){ # if the highest peak is first, then damping is likely
      if (length(peaks)>1){
        # trying to figure out gamma based on logarithmic decrement
        n <- peaks_time[2]-peaks_time[1]
        log_dec <- (1/n)*log(abs(peaks[1]/peaks[2]))
        gam0 <- 1/(sqrt(1+((2*pi/log_dec)^2)))
      } else{ # use a guess
        gam0 <- .01
      }
    } else{ # otherwise driving is likely
      if (length(peaks)>1){
        # trying to figure out gamma based on logarithmic decrement
        n <- peaks_time[2]-peaks_time[1]
        log_dec <- (1/n)*log(abs(peaks[2]/peaks[1]))
        gam0 <- -1*1/(sqrt(1+((2*pi/log_dec)^2)))
      } else{ # use a guess
        gam0 <- -.01
      }
    }
    
    # let frequency depend on amount of peaks = (length(timen)*resol/(no of peaks+1 [accounts for phase shift])
    if (length(peaks) == 0){
      if (high == -Inf || low == Inf){
        w0 <- 2*pi/(length(timen)*resol/2)
      } else{
        # want to get their actual integer period values
        highfix <- (high/2/pi)^-1
        lowfix <- (low/2/pi)^-1
        w0 <- 2*pi/(length(timen)*resol/((highfix+lowfix)/2))
      }
    } else if (length(peaks) == 1){ # phase shift causes only one peak to appear
      w0 <- 2*pi/(length(timen)*resol/(length(peaks)+1))
    } else{
      w0 <- 2*pi/(length(timen)*resol/(length(peaks)))
    }
    
    # can't be outside the specified parameters
    if (w0 > low){
      w0 <- low
    } else if (w0 < high){
      w0 <- high
    }
    
    
    # we estimate our phase shift on the second and third nonmissing points for accuracy
    # if you have less than 3 points nonmissing, I have no hope for you
    second <- which(!is.na(y_val))[2]
    third <- which(!is.na(y_val))[3]
    min_i <- 0;
    min_vect <- rep(10000, length(0:11))
    for (i in 0:11){
      # if the phase shift at both the second and third gene values are the smallest value available for the fitted value
      min_vect[i+1] <- sum(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[second])-y_val[second],
                           alt_form(a0,gam0,w0,(i*pi/6),y0,timen[third])-y_val[third])
    }
    # phi0 <- min_i*pi/6 # intial value for phase shift
    phi0 <- (which.min(abs(min_vect))-1)*pi/6 # intial value for phase shift
    
    start_param <- list(gam=gam0,a=a0,omega=w0,phi=phi0,y_shift=y0)
    temp <- data.frame()
    if (num_reps == 1){ # one replicate
      # put the timen into a data frame
      temp <- data.frame(y=t(y_val),t=timen)
      temp <- temp[!is.na(temp$y),] # remove any missing data points
      
      # fitting
      oscillator.fit <- suppressWarnings(nlsLM(y ~ alt_form(a,gam,omega,phi,y_shift,t),
                                               data=temp,
                                               start=start_param,
                                               lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
                                               upper=c(Inf, Inf, low, Inf, max(temp$y)),
                                               control = nls.lm.control(maxiter = 1000, maxfev = 2000,
                                                                        ftol=1e-6, ptol=1e-6, gtol=1e-6)))
      
    } else{ # multiple replicates
      #put the timen and data point into a data frame
      weights <- calc_weights(current_gene,num_reps, genes)
      temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),t=cbind(rep(timen,each = num_reps)),w=cbind(rep(weights,each = num_reps)))
      temp <- temp[!is.na(temp$y),] # remove any missing data points
      
      #fitting
      oscillator.fit <- suppressWarnings(nlsLM(y ~ alt_form(a,gam,omega,phi,y_shift,t),
                                               data=temp,
                                               start=start_param,
                                               lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
                                               upper=c(Inf, Inf, low, Inf, max(temp$y)),
                                               control = nls.lm.control(maxiter = 1000, maxfev = 2000,
                                                                        ftol=1e-6, ptol=1e-6, gtol=1e-6),
                                               weights = w))
    }
    
    did_conv <- oscillator.fit$convInfo$isConv # did the fit converge
    num_iter <- oscillator.fit$convInfo$finIter # amount of iterations
    
    parameters <- coef(oscillator.fit) #extract parameter estimates
    
    # alt_form parameters:
    # the parameters go in the order of: gam,a,omega,phi,y_shift
    gam <- parameters[1]
    a <- parameters[2]
    omega <- parameters[3]
    phi <- parameters[4]
    y_shift <- parameters[5]
    
    if (run_conf){
      if (which_conf == "Bootstrap"){
        ci_int <- bootstrap(temp, oscillator.fit, start_param, num_reps, current_gene, seed)
      } else {
        ci_int <- jackknife(temp, parameters, num_reps, start_param)
      }
    }
    
    # calculating whether (over)damped, (over)forced, harmonic
    if (gam < -over_cut){
      type_gam <- "Overexpressed"
    } else if (gam <= -harm_cut){
      type_gam <- "Forced"
    } else if (gam <= harm_cut){
      type_gam <- "Harmonic"
    } else if (gam <= over_cut){
      type_gam <- "Damped"
    } else{
      type_gam <- "Repressed"
    }
    
    # calculating the phase shift in terms of period (omega inverse of period)
    if (!is.na(a)){ # all param will either be na or not
      if (a >= 0){ # positive amplitudes
        if (phi > 0){ # shift to the left
          frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
          dist_peak <- frac_part*(2*pi/omega) # distance from first peak
          phase_hours <- (2*pi/omega)-dist_peak
        } else { # shift to the right
          frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
          dist_peak <- frac_part*(2*pi/omega) # distance from first peak
          phase_hours <- abs(dist_peak)
        }
      } else { # negative ampltitudes
        if (phi > 0){ # shift to the left
          frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
          dist_peak <- frac_part*(2*pi/omega) # distance from first peak
          if (abs(frac_part) < .5){
            phase_hours <- (2*pi/omega)-dist_peak - (2*pi/omega/2)
          } else {
            phase_hours <- (2*pi/omega)-dist_peak + (2*pi/omega/2)
          }
        } else { # shift to the right
          frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
          dist_peak <- frac_part*(2*pi/omega) # distance from first peak
          if (abs(frac_part) < .5){
            phase_hours <- abs(dist_peak) + (2*pi/omega/2)
          } else {
            phase_hours <- abs(dist_peak) - (2*pi/omega/2)
          }
        }
      }
    } else {
      phase_hours <- NA
    }
    
    # calculate p-value
    ref_wave <- (alt_form(a,gam,omega,phi,y_shift,timen)) # fitted values
    all_pred <- rep(ref_wave, each=num_reps)[!is.na(unlist(genes[current_gene,-1]))]
    testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
    pval <- testing$p.value
    tau <- testing$estimate
    
    # list of parameters and other resulting values
    results <- data.frame(gene = gene_n, conv = did_conv, iter = num_iter, gamma = gam, type_gam = type_gam,amplitude = a, omega = omega, period = (2*pi/omega), phase.shift = phi,hours.shifted = phase_hours, y_shift=y_shift, tau = tau, pval = pval, stringsAsFactors = FALSE)
    if (!run_conf){
      if (num_reps == 1){
        results <- cbind(results, y_val, rbind(ref_wave))
      } else {
        results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(ref_wave))
      }
    } else {
      if (num_reps == 1){
        results <- cbind(results, rbind(ci_int), y_val, rbind(ref_wave))
      } else {
        results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(ref_wave))
      }
    }
    return (results)
  }, error = function(e){ # if there's failure in convergence
    results <- data.frame(gene = gene_n, conv = NA, iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
    if (!run_conf){
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
    } else {
      results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
    }
    
    return (results)
  })
}

calculate_param_update <- function(current_gene,timen,resol,num_reps,tied,is_smooth=FALSE,is_weighted=FALSE,low,high,rem_unexpr=FALSE,rem_unexpr_amt=70,run_conf = F, which_conf = "Bootstrap", harm_cut = .03, over_cut = .15, seed = 30){
  
  if (run_conf){
    ci_int <- rep(NA, 10)
    par_names <- c("gam","a","omega","phi","y_shift")
    ci_names <- c(paste0("ci_low_",par_names), paste0("ci_high_",par_names))
    names(ci_int) <- ci_names
  }
  
  gene_n <- as.character(genes[current_gene,1]) # gene name
  # first we need to check whether or not the gene is just a straight line
  if (!is_deviating(current_gene)){ # one replicate
    results <- data.frame(gene = gene_n, conv = "No Deviation", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA,  tau = NA, pval = NA, stringsAsFactors = FALSE)
    
    if (run_conf){
      results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
    } else {
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
    }
    return(results)
  }
  
  # then we need to check if < threshold % are expressed (if desired)
  if (rem_unexpr){
    if (rem_unexpr_vect[current_gene]){
      results <- data.frame(gene = gene_n, conv = "Unexpressed", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
      if (run_conf){
        results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
      } else {
        results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
      }
      return(results)
    }
  }
  
  # calculating averages for initial values
  if (num_reps == 1){
    y_val <- rbind(as.numeric(as.character(t(genes[current_gene,c(2:ncol(genes))])))) # all the gene values
  } else{
    y_val <- avg_genes[current_gene,] # starting values determined by average of replicates
  }
  
  tryCatch({ # throw exception upon error
    start_param <- calc_start_echo(y_val)
    temp <- data.frame()
    if (num_reps == 1){ # one replicate
      # put the timen into a data frame
      temp <- data.frame(y=t(y_val),t=timen)
      temp <- temp[!is.na(temp$y),] # remove any missing data points
      
      # fitting
      oscillator.fit <- suppressWarnings(nlsLM(y ~ alt_form(a,gam,omega,phi,y_shift,t),
                                               data=temp,
                                               start=start_param,
                                               lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
                                               upper=c(Inf, Inf, low, Inf, max(temp$y)),
                                               control = nls.lm.control(maxiter = 1000, maxfev = 2000,
                                                                        ftol=1e-6, ptol=1e-6, gtol=1e-6)))
      
    } else{ # multiple replicates
      #put the timen and data point into a data frame
      weights <- calc_weights(current_gene,num_reps, genes)
      temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),t=cbind(rep(timen,each = num_reps)),w=cbind(rep(weights,each = num_reps)))
      temp <- temp[!is.na(temp$y),] # remove any missing data points
      
      #fitting
      # oscillator.fit <- suppressWarnings(nlsLM(y ~ alt_form(a,gam,omega,phi,y_shift,t),
      #                                          data=temp,
      #                                          start=start_param,
      #                                          lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
      #                                          upper=c(Inf, Inf, low, Inf, max(temp$y)),
      #                                          control = nls.lm.control(maxiter = 1000, maxfev = 2000,
      #                                                                   ftol=1e-6, ptol=1e-6, gtol=1e-6),
      #                                          weights = w))
      
      #fitting
      oscillator_init <- nls.lm(par = start_param,
                            fn = fcn_single,
                            fcall = alt_form,
                            lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
                            upper=c(Inf, Inf, low, Inf, max(temp$y)),
                            t = temp$t,
                            y = temp$y,
                            w = temp$w,
                            control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                     ftol=1e-6, ptol=1e-6, gtol=1e-6)
      )
      
      coeff <- coef(oscillator_init)
      
      oscillator.fit <- nls_edit(
        formula = y ~ exp_mod(a, gam, y_shift, t),
        data=temp,
        start = list(a=coeff[1],gam=coeff[2], y_shift=coeff[3]),
        control = nls.control(maxiter = 0)
      )
    }
    
    did_conv <- oscillator.fit$convInfo$isConv # did the fit converge
    num_iter <- oscillator.fit$convInfo$finIter # amount of iterations
    
    parameters <- coef(oscillator.fit) #extract parameter estimates
    
    # alt_form parameters:
    # the parameters go in the order of: gam,a,omega,phi,y_shift
    gam <- parameters[1]
    a <- parameters[2]
    omega <- parameters[3]
    phi <- parameters[4]
    y_shift <- parameters[5]
    
    if (run_conf){
      if (which_conf == "Bootstrap"){
        ci_int <- bootstrap(temp, oscillator.fit, start_param, num_reps, current_gene, seed)
      } else {
        ci_int <- jackknife(temp, parameters, num_reps, start_param)
      }
    }
    
    # calculating whether (over)damped, (over)forced, harmonic
    if (gam < -over_cut){
      type_gam <- "Overexpressed"
    } else if (gam <= -harm_cut){
      type_gam <- "Forced"
    } else if (gam <= harm_cut){
      type_gam <- "Harmonic"
    } else if (gam <= over_cut){
      type_gam <- "Damped"
    } else{
      type_gam <- "Repressed"
    }
    
    # calculating the phase shift in terms of period (omega inverse of period)
    if (!is.na(a)){ # all param will either be na or not
      if (a >= 0){ # positive amplitudes
        if (phi > 0){ # shift to the left
          frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
          dist_peak <- frac_part*(2*pi/omega) # distance from first peak
          phase_hours <- (2*pi/omega)-dist_peak
        } else { # shift to the right
          frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
          dist_peak <- frac_part*(2*pi/omega) # distance from first peak
          phase_hours <- abs(dist_peak)
        }
      } else { # negative ampltitudes
        if (phi > 0){ # shift to the left
          frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
          dist_peak <- frac_part*(2*pi/omega) # distance from first peak
          if (abs(frac_part) < .5){
            phase_hours <- (2*pi/omega)-dist_peak - (2*pi/omega/2)
          } else {
            phase_hours <- (2*pi/omega)-dist_peak + (2*pi/omega/2)
          }
        } else { # shift to the right
          frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
          dist_peak <- frac_part*(2*pi/omega) # distance from first peak
          if (abs(frac_part) < .5){
            phase_hours <- abs(dist_peak) + (2*pi/omega/2)
          } else {
            phase_hours <- abs(dist_peak) - (2*pi/omega/2)
          }
        }
      }
    } else {
      phase_hours <- NA
    }
    
    # calculate p-value
    ref_wave <- (alt_form(a,gam,omega,phi,y_shift,timen)) # fitted values
    all_pred <- rep(ref_wave, each=num_reps)[!is.na(unlist(genes[current_gene,-1]))]
    testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
    pval <- testing$p.value
    tau <- testing$estimate
    
    # list of parameters and other resulting values
    results <- data.frame(gene = gene_n, conv = did_conv, iter = num_iter, gamma = gam, type_gam = type_gam,amplitude = a, omega = omega, period = (2*pi/omega), phase.shift = phi,hours.shifted = phase_hours, y_shift=y_shift, tau = tau, pval = pval, stringsAsFactors = FALSE)
    if (!run_conf){
      if (num_reps == 1){
        results <- cbind(results, y_val, rbind(ref_wave))
      } else {
        results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(ref_wave))
      }
    } else {
      if (num_reps == 1){
        results <- cbind(results, rbind(ci_int), y_val, rbind(ref_wave))
      } else {
        results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(ref_wave))
      }
    }
    return (results)
  }, error = function(e){ # if there's failure in convergence
    results <- data.frame(gene = gene_n, conv = NA, iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
    if (!run_conf){
      results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
    } else {
      results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
    }
    
    return (results)
  })
}

calculate_param_both <- function(current_gene,timen,resol,num_reps,tied,is_smooth=FALSE,is_weighted=FALSE,low,high,rem_unexpr=FALSE,rem_unexpr_amt=70,run_conf = F, which_conf = "Bootstrap", harm_cut = .03, over_cut = .15, seed = 30, start_type = "orig"){
  
  if (run_conf){
    ci_int <- rep(NA, 10)
    par_names <- c("gam","a","omega","phi","y_shift")
    ci_names <- c(paste0("ci_low_",par_names), paste0("ci_high_",par_names))
    names(ci_int) <- ci_names
  }
  
  gene_n <- as.character(genes[current_gene,1]) # gene name
  
  # FOR NOW, NO CHECKING
  # first we need to check whether or not the gene is just a straight line
  # if (!is_deviating(current_gene)){ # one replicate
  #   results <- data.frame(gene = gene_n, conv = "No Deviation", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA,  tau = NA, pval = NA, stringsAsFactors = FALSE)
  #   
  #   if (run_conf){
  #     results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
  #   } else {
  #     results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
  #   }
  #   return(results)
  # }
  
  # then we need to check if < threshold % are expressed (if desired)
  if (rem_unexpr){
    if (rem_unexpr_vect[current_gene]){
      results <- data.frame(gene = gene_n, conv = "Unexpressed", iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
      if (run_conf){
        results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
      } else {
        results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
      }
      return(results)
    }
  }
  
  # calculating averages for initial values
  if (num_reps == 1){
    y_val_rna <- rbind(as.numeric(as.character(t(genes_rna[current_gene,c(2:ncol(genes))])))) # all the gene values
    y_val_pro <- rbind(as.numeric(as.character(t(genes_pro[current_gene,c(2:ncol(genes))])))) # all the gene values
  } else{
    y_val_pro <- avg_genes_pro[current_gene,] # starting values determined by average of replicates
    y_val_rna <- avg_genes_rna[current_gene,] # starting values determined by average of replicates
  }
  
  # tryCatch({ # throw exception upon error
    if (start_type == "orig"){
      rna_start <- calc_start_echo(y_val_rna)
      pro_start <- calc_start_echo(y_val_pro)
    } else {
      if (!is.na(rna_tot$Amplitude.Change.Coefficient[current_gene]) 
          & rna_tot$`BH Adj P-Value`[current_gene] < .0005
          & rna_tot$`Oscillation Type` != "Overexpressed"
          & rna_tot$`Oscillation Type` != "Repressed"){
        rna_start <- list("gam" = rna_tot$Amplitude.Change.Coefficient[current_gene],
                          "a" = rna_tot$Initial.Amplitude[current_gene],
                          "omega" = rna_tot$Radian.Frequency[current_gene],
                          "phi" = rna_tot$`Phase Shift`[current_gene],
                          "y_shift" = rna_tot$`Equilibrium Value`[current_gene])
      } else {
        rna_start <- calc_start_echo(y_val_rna)
      }
      
      
      if (!is.na(pro_tot$Amplitude.Change.Coefficient[current_gene]) &
          pro_tot$`BH Adj P-Value`[current_gene] < .0005
          & rna_tot$`Oscillation Type` != "Overexpressed"
          & rna_tot$`Oscillation Type` != "Repressed"){
        pro_start <- list("gam" = pro_tot$Amplitude.Change.Coefficient[current_gene],
                          "a" = pro_tot$Initial.Amplitude[current_gene],
                          "omega" = pro_tot$Radian.Frequency[current_gene],
                          "phi" = pro_tot$`Phase Shift`[current_gene],
                          "y_shift" = pro_tot$`Equilibrium Value`[current_gene])
      } else {
        pro_start <- calc_start_echo(y_val_pro)
      }
    }
    # we're going to use the starting value of rna for the joint period, respective for everything else
    # now we need to format it into a way that is nice
    names(rna_start) <- paste0(names(rna_start),"1")
    names(pro_start) <- paste0(names(pro_start),"2")
    start_param <- c(rna_start,pro_start)
    # start_param[["y_shift1"]] <- unlist(alph_rna[current_gene])
    # start_param[["y_shift2"]] <- unlist(alph_pro[current_gene])
    # start_param[["slo1"]] <- unlist(beta_rna[current_gene])
    # start_param[["slo2"]] <- unlist(beta_pro[current_gene])
    start_param$omega2 <- NULL
    
    y_stack <- c(unlist(genes_rna[current_gene,-1]), unlist(genes_pro[current_gene,-1]))
    t_stack <- c(rep(timen, each=num_reps), (rep(timen, each = num_reps)))
    weights <- c(rep(calc_weights_boot(unlist(genes_rna[current_gene,-1]), num_reps), each = num_reps),
                 rep(calc_weights_boot(unlist(genes_pro[current_gene,-1]), num_reps), each = num_reps))
    
    t_rep <- length(rep(timen, each=num_reps))
    # ASSUME NO MISSING DATA ATM
    s1 <- rep(c(1,0), c(t_rep,t_rep))
    s2 <- rep(c(0,1), c(t_rep,t_rep))
    
    # CONSTRAINTS ARE IN THE ORDER OF THE START PARAM
    lower_constr <- upper_constr <- c(rep(1, length(start_param)))
    names(lower_constr) <- names(upper_constr) <- names(start_param)
    # specific constr
    # lower_constr <- c( -.3,-(3*min(abs(y_stack))), high, -4*pi, min(y_stack[1:(t_rep)]),
    #                    -.3,-(3*min(abs(y_stack))), -4*pi, min(y_stack[((t_rep)+1):length(y_stack)]))
    # upper_constr <- c( .3, (3*max(abs(y_stack))), low, 4*pi, max(y_stack[1:t_rep]),
    #                    .3, (3*min(abs(y_stack))), 4*pi, max(y_stack[((t_rep)+1):length(y_stack)]))
    # back to normal constraints
    lower_constr <- c( -Inf,-Inf, high, -Inf, min(y_stack[1:(t_rep)]),
                       -Inf,-Inf, -Inf, min(y_stack[((t_rep)+1):length(y_stack)]))
    upper_constr <- c( Inf, Inf, low, Inf, max(y_stack[1:t_rep]),
                       Inf, Inf, Inf, max(y_stack[((t_rep)+1):length(y_stack)]))
    
    temp <- data.frame()
    if (num_reps == 1){ # one replicate
      # put the timen into a data frame
      temp <- data.frame(y=t(y_val),t=timen)
      temp <- temp[!is.na(temp$y),] # remove any missing data points
      
      # fitting
      oscillator.fit <- suppressWarnings(nlsLM(y ~ alt_form(a,gam,omega,phi,y_shift,t),
                                               data=temp,
                                               start=start_param,
                                               lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
                                               upper=c(Inf, Inf, low, Inf, max(temp$y)),
                                               control = nls.lm.control(maxiter = 1000, maxfev = 2000,
                                                                        ftol=1e-6, ptol=1e-6, gtol=1e-6)))
      
    } else{ # multiple replicates
      #put the timen and data point into a data frame
      # weights <- calc_weights(current_gene,num_reps, genes)
      # temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),t=cbind(rep(timen,each = num_reps)),w=cbind(rep(weights,each = num_reps)))
      temp <- data.frame(y = y_stack, t =t_stack, w= weights, s1 = s1, s2 = s2)
      temp <- temp[!is.na(temp$y),] # remove any missing data points-NONE IN EXAMPLE
      
      # oscillator.fit <- suppressWarnings(nlsLM(y ~ alt_form(a2,gam2,omega2,phi2,y_shift2,t),
      #                                          data=temp,
      #                                          start=pro_start,
      #                                          lower=c(-Inf, -Inf, high, -Inf, min(temp$y)),
      #                                          upper=c(Inf, Inf, low, Inf, max(temp$y)),
      #                                          control = nls.lm.control(maxiter = 1000, maxfev = 2000,
      #                                                                   ftol=1e-6, ptol=1e-6, gtol=1e-6),
      #                                          weights = w))
      oscillator.fit <- nls.lm(par = start_param,
                          fn = fcn,
                          jac = fcn.jac,
                          fcall = echo_joint,
                          jcall = jac_echo_joint,
                          lower = lower_constr,
                          upper = upper_constr,
                          t = temp$t,
                          y = temp$y,
                          s1 = temp$s1,
                          s2 = temp$s2,
                          w = temp$w,
                          control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                   ftol=1e-6, ptol=1e-6, gtol=1e-6)
      )
      
      par <- coef(oscillator.fit) #extract parameter estimates
      
      if (is.na(par[["gam1"]])){
        print(paste("here we go again", current_gene))
        
        rna_start <- calc_start_echo(y_val_rna)
        pro_start <- calc_start_echo(y_val_pro)
        names(rna_start) <- paste0(names(rna_start),"1")
        names(pro_start) <- paste0(names(pro_start),"2")
        start_param <- c(rna_start,pro_start)
        # start_param[["y_shift1"]] <- unlist(alph_rna[current_gene])
        # start_param[["y_shift2"]] <- unlist(alph_pro[current_gene])
        # start_param[["slo1"]] <- unlist(beta_rna[current_gene])
        # start_param[["slo2"]] <- unlist(beta_pro[current_gene])
        start_param$omega2 <- NULL
        
        oscillator.fit <- nls.lm(par = start_param,
                                 fn = fcn,
                                 jac = fcn.jac,
                                 fcall = echo_joint,
                                 jcall = jac_echo_joint,
                                 lower = lower_constr,
                                 upper = upper_constr,
                                 t = temp$t,
                                 y = temp$y,
                                 s1 = temp$s1,
                                 s2 = temp$s2,
                                 w = temp$w,
                                 control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                          ftol=1e-6, ptol=1e-6, gtol=1e-6)
        )
      }
      
      # oscillator.fit <- suppressWarnings(nlsLM(y~echo_joint(a1,a2, gam1,gam2, omega1, phi1,phi2, y_shift1,y_shift2, t, s1, s2),
      #                                          jac = jac_echo_joint(a1,a2, gam1,gam2, omega1, phi1,phi2, y_shift1,y_shift2, t, s1, s2),
      #                                          data=temp,
      #                                          start=start_param,
      #                                          lower=lower_constr,
      #                                          upper=upper_constr,
      #                                          control = nls.lm.control(maxiter = 1000, maxfev = 2000,
      #                                                                   ftol=1e-6, ptol=1e-6, gtol=1e-6),
      #                                          weights = w))
      
      
    }
    
    # did_conv <- oscillator.fit$convInfo$isConv # did the fit converge
    num_iter <- oscillator.fit$niter # amount of iterations
    
    par <- coef(oscillator.fit) #extract parameter estimates
    res <- as.data.frame(t(par))
    
    # calculate pvalues
    fit <- echo_joint(res$a1,res$a2, res$gam1,res$gam2, res$omega1, res$phi1,res$phi2, res$y_shift1,res$y_shift2, temp$t, temp$s1, temp$s2)
    res$both_pval <- cor.test(fit, temp$y, method = "kendall")$p.value
    res$rna_pval <- cor.test(fit[1:(t_rep)], temp$y[1:(t_rep)], method = "kendall")$p.value
    res$pro_pval <- cor.test(fit[((t_rep)+1):length(y_stack)], temp$y[((t_rep)+1):length(y_stack)], method = "kendall")$p.value
    
    res <- cbind(data.frame("gene_n"=gene_n, stringsAsFactors = F), res)
    return(res)
    
    # alt_form parameters:
    # the parameters go in the order of: gam,a,omega,phi,y_shift
    gam <- parameters[1]
    a <- parameters[2]
    omega <- parameters[3]
    phi <- parameters[4]
    y_shift <- parameters[5]
    
    # NEED TO REWRITE BOOTSTRAPPING
    if (run_conf){
      if (which_conf == "Bootstrap"){
        ci_int <- bootstrap(temp, oscillator.fit, start_param, num_reps, current_gene, seed)
      } else {
        ci_int <- jackknife(temp, parameters, num_reps, start_param)
      }
    }
    
    # calculating whether (over)damped, (over)forced, harmonic
    if (gam < -over_cut){
      type_gam <- "Overexpressed"
    } else if (gam <= -harm_cut){
      type_gam <- "Forced"
    } else if (gam <= harm_cut){
      type_gam <- "Harmonic"
    } else if (gam <= over_cut){
      type_gam <- "Damped"
    } else{
      type_gam <- "Repressed"
    }
    
    # calculating the phase shift in terms of period (omega inverse of period)
    if (!is.na(a)){ # all param will either be na or not
      if (a >= 0){ # positive amplitudes
        if (phi > 0){ # shift to the left
          frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
          dist_peak <- frac_part*(2*pi/omega) # distance from first peak
          phase_hours <- (2*pi/omega)-dist_peak
        } else { # shift to the right
          frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
          dist_peak <- frac_part*(2*pi/omega) # distance from first peak
          phase_hours <- abs(dist_peak)
        }
      } else { # negative ampltitudes
        if (phi > 0){ # shift to the left
          frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
          dist_peak <- frac_part*(2*pi/omega) # distance from first peak
          if (abs(frac_part) < .5){
            phase_hours <- (2*pi/omega)-dist_peak - (2*pi/omega/2)
          } else {
            phase_hours <- (2*pi/omega)-dist_peak + (2*pi/omega/2)
          }
        } else { # shift to the right
          frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
          dist_peak <- frac_part*(2*pi/omega) # distance from first peak
          if (abs(frac_part) < .5){
            phase_hours <- abs(dist_peak) + (2*pi/omega/2)
          } else {
            phase_hours <- abs(dist_peak) - (2*pi/omega/2)
          }
        }
      }
    } else {
      phase_hours <- NA
    }
    
    # calculate p-value
    ref_wave <- (alt_form(a,gam,omega,phi,y_shift,timen)) # fitted values
    all_pred <- rep(ref_wave, each=num_reps)[!is.na(unlist(genes[current_gene,-1]))]
    testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
    pval <- testing$p.value
    tau <- testing$estimate
    
    # list of parameters and other resulting values
    results <- data.frame(gene = gene_n, conv = did_conv, iter = num_iter, gamma = gam, type_gam = type_gam,amplitude = a, omega = omega, period = (2*pi/omega), phase.shift = phi,hours.shifted = phase_hours, y_shift=y_shift, tau = tau, pval = pval, stringsAsFactors = FALSE)
    if (!run_conf){
      if (num_reps == 1){
        results <- cbind(results, y_val, rbind(ref_wave))
      } else {
        results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(ref_wave))
      }
    } else {
      if (num_reps == 1){
        results <- cbind(results, rbind(ci_int), y_val, rbind(ref_wave))
      } else {
        results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(ref_wave))
      }
    }
    return (results)
  # }, error = function(e){ # if there's failure in convergence
  #   results <- data.frame(gene = gene_n, conv = NA, iter = NA, gamma = NA, type_gam = NA, amplitude = NA, omega = NA, period = NA, phase.shift = NA, hours.shifted = NA, y_shift=NA, tau = NA, pval = NA, stringsAsFactors = FALSE)
  #   if (!run_conf){
  #     results <- cbind(results, rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
  #   } else {
  #     results <- cbind(results, rbind(ci_int), rbind(as.numeric(as.character(t(genes[current_gene,-1])))), rbind(rep(NA,length(timen))))
  #   }
  #   
  #   return (results)
  # })
}

run_echo <- function(genes, timen, resol,num_reps,tied,is_smooth=FALSE,is_weighted=FALSE,low,high,rem_unexpr=FALSE,rem_unexpr_amt=70,run_conf = F, which_conf = "Bootstrap", harm_cut = .03, over_cut = .15, seed = 30){
  # prepare for parallelism
  cores <- detectCores() # dectect how many processors
  cl <- makeCluster(cores[1]-1) # not to overload your computer, need one for OS
  registerDoSNOW(cl)
  
  # making a progress bar
  if (nrow(genes)==1){
    print(paste("Percentage finished out of",nrow(genes)-1,"expression:"))
  } else {
    print(paste("Percentage finished out of",nrow(genes),"expressions:"))
  }
  pb <- txtProgressBar(max = nrow(genes), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # where we put the result
  total_results <- foreach (i=1:nrow(genes), .combine = rbind, .packages=c('minpack.lm',"boot"),.options.snow = opts,.export=ls(envir=globalenv())) %dopar% {
    calculate_param_update(i, timen, resol, num_reps, tied = tied, is_smooth = is_smooth, is_weighted = is_weighted,low = low,high = high,rem_unexpr = rem_unexpr, rem_unexpr_amt = rem_unexpr_amt, run_conf = run_conf, which_conf = which_conf, harm_cut = harm_cut, over_cut = over_cut, seed = seed)
  }
  
  close(pb)
  
  stopCluster(cl) # stop using the clusters
  
  # renaming columns of the final results
  if (!run_conf){
    colnames(total_results) <- c("Gene Name","Convergence","Iterations","Amplitude.Change.Coefficient","Oscillation Type","Initial.Amplitude","Radian.Frequency","Period","Phase Shift","Hours Shifted","Equilibrium Value", "Tau", "P-Value", paste(rep("Original TP",length(rep(timen, each = num_reps))),rep(timen, each = num_reps),rep(".",length(rep(timen, each = num_reps))),rep(c(1:num_reps), length(timen)),sep=""), paste(rep("Fitted TP",length(timen)),timen,sep=""))
  } else {
    conf_int_names <- c("CI.AC.Coeff","CI.Init.Amp","CI.Rad.Freq","CI.Phase.Shift","CI.Eq.Val")
    conf_int_names <- c(paste0(conf_int_names,".Low"), paste0(conf_int_names,".High"))
    colnames(total_results) <- c("Gene Name","Convergence","Iterations","Amplitude.Change.Coefficient","Oscillation Type","Initial.Amplitude","Radian.Frequency","Period","Phase Shift","Hours Shifted","Equilibrium Value", "Tau", "P-Value", conf_int_names, paste(rep("Original TP",length(rep(timen, each = num_reps))),rep(timen, each = num_reps),rep(".",length(rep(timen, each = num_reps))),rep(c(1:num_reps), length(timen)),sep=""), paste(rep("Fitted TP",length(timen)),timen,sep=""))
  }
  # remove the fake row I added if there is only one gene
  if (nrow(genes) == 1){
    total_results <- total_results[-nrow(total_results),]
  }
  
  # add slope
  beta <- rep(NA, nrow(genes))
  total_results <- cbind(total_results[,c(1:11)],`Slope` = beta, total_results[,c(12:ncol(total_results))])
  
  adjusted_p_val_us <- p.adjust(unlist(total_results$`P-Value`), method = "BH") # benjamini-hochberg adjust p-values
  # add BH adjusted pvalue
  total_results <- cbind(total_results[,c(1:14)],`BH Adj P-Value` = adjusted_p_val_us, total_results[,c(15:ncol(total_results))])
  
  # adding the benjamini-hochberg-yekutieli p-value adjustment
  total_results <- cbind(total_results[,c(1:15)],`BY Adj P-Value` = p.adjust(unlist(total_results$`P-Value`), method = "BY"), total_results[,c(16:ncol(total_results))])
  
  
  
  return(total_results)
}

# CLEAN THIS UP
smooth_y_val <- function(y_val){
  # smoothing the starting averages - weighted or unweighted?
  #get matrix of just the relative expression over time
  all_reps <- matrix(c(y_val), nrow = 1, ncol = length(y_val))
  
  is_weighted <- T
  
  # weighted averaging
  
  #if there are replicates, average the relative expression for each replicate
  center_reps <- list() # to store actual replicate matrix
  mtx_count <- list() # to store how many are NA
  for (i in 1:1){
    center_reps[[i]] <- all_reps[, seq(i,ncol(all_reps), by=1)]*2
    mtx_count[[i]] <- is.na(center_reps[[i]])*2
    center_reps[[i]][is.na(center_reps[[i]])] <- 0
  }
  
  repmtx_l <- list() # store amount to divide by for each rep
  repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(all_reps))+4 # to store how many we should divide by
  repmtx[,c(1,ncol(repmtx))] <- repmtx[,c(1,ncol(repmtx))] - 1 # we only divide by 3 at edges
  
  
  # sum the replicates
  left <- c(rep(0,1),center_reps[[i]][-length(center_reps[[i]])]/2) # left shifted matrix
  right <- c(center_reps[[i]][-1]/2,rep(0,1)) # right shifted matrix
  center_reps[[i]] <- left + center_reps[[i]] + right
  
  # figure out how many replicates are actually available for each time point
  left_na <- c(rep(0,1),mtx_count[[i]][-length(mtx_count[[i]])]/2) # left shifted matrix
  right_na <- c(mtx_count[[i]][-1]/2,rep(0,1)) # right shifted matrix
  repmtx_l[[i]] <- repmtx - left_na - mtx_count[[i]] - right_na
  # to avoid division by 0 and induce NAs if there are no time points available
  repmtx_l[[i]][repmtx_l[[i]]==0] <- NA
  
  dat <- all_reps
  x <- 0
  dat[,seq(1+x,ncol(all_reps),by=1)] <- center_reps[[x+1]]/repmtx_l[[x+1]]
  dat[is.na(all_reps)] <- NA # do not impute missing values
  
  return(dat[1,])
}

calc_peaks <- function(y_val, resol, y0, timen){
  # figure out the resolution modifier
  # note to self: as long as we stay under 24 hours, we should be good
  mod <- 0
  if (resol <= ((1/12)+10^-8)){ # 17 hour surround
    mod <- 102
  } else if (resol <= ((1/6)+10^-8)){ # 15 hour surround
    mod <- 45
  } else if (resol <= ((1/4)+10^-8)){ # 13 hour surround
    mod <- 26
  } else if (resol <= ((1/2)+10^-8)){ # 11 hour surround
    mod <- 11
  } else if (resol <= 1){ # 8 hour surround
    mod <- 4
  } else if (resol <=2){ # 6 hour surround
    mod <- 3
  } else if (resol <= 4){ # 4 hour surround
    mod <- 2
  } else{
    mod <- 1
  }
  
  # calculate the amount of peaks
  peaks <- c(); # vector of peak values
  peaks_time <- c(); # vector of peak timen
  counting <- 1; # counter
  
  for(i in (mod+1):(length(y_val)-mod)){
    # deal with complete missingness
    if (suppressWarnings(max(y_val[(i-mod):(i+mod)], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
      next
    }
    # otherwise continue as normal
    if (y_val[i] == max(y_val[(i-mod):(i+mod)], na.rm = TRUE)){
      peaks[counting] <- y_val[i]
      peaks_time[counting] <- timen[i]
      counting <- counting+1
    }
  }
  
  # if (length(peaks) <= 1){ # phase might be out of sync
  troughs <- c()
  troughs_time <- c()
  counting <- 1
  
  for(i in (mod+1):(length(y_val)-mod)){
    # deal with complete missingness
    if (suppressWarnings(min(y_val[(i-mod):(i+mod)], na.rm = TRUE)) == -Inf | is.na(y_val[i])){
      next
    }
    # otherwise continue as normal
    if (y_val[i] == min(y_val[(i-mod):(i+mod)], na.rm = TRUE)){
      troughs[counting] <- y_val[i]
      troughs_time[counting] <- timen[i]
      counting <- counting+1
    }
  }
  
  # old criterion:
  # if we have now found more troughts than peaks, we want to make those our peaks
  peaks_per_diff <- 0
  troughs_per_diff <- 0
  if (length(peaks)>1 & length(troughs)>1){
    # new criterion: which peaks are more evenly distributed?
    peaks_per <- length(timen)*resol/length(peaks)
    troughs_per <- length(timen)*resol/length(troughs)
    
    # get avg per from the timen
    peaks_avg_per <- mean(sapply(seq(length(peaks_time),2,-1), 
                                 function(x){peaks_time[x]-peaks_time[x-1]}))
    troughs_avg_per <- mean(sapply(seq(length(troughs_time),2,-1), 
                                   function(x){troughs_time[x]-troughs_time[x-1]}))
    
    # see which is probably more consistent
    peaks_per_diff <- abs(peaks_per-peaks_avg_per)
    troughs_per_diff <- abs(troughs_per-troughs_avg_per)
  }
  
  if ((length(peaks)<=1 & length(troughs)>length(peaks)) | 
      (troughs_per_diff < peaks_per_diff)){  
    # flip the troughs for accurate calculations later
    # absolute value the difference from the midline
    peaks <- abs(y0 - troughs)
    peaks_time <- troughs_time
  } else {
    # absolute value the difference from the midline
    peaks <- abs(peaks - y0)
  }
  # it's possible for this not to work if you don't have more than one
  # oscillation over your time course
  # }
  
  return(list("peaks"=peaks, "peaks_time"=peaks_time))
  
}

calc_start_echo <- function(y_val){
  timen <- timen
  y_val <- smooth_y_val(y_val)
  
  # calculate starting y_shift
  y0 <- mean(y_val,na.rm = TRUE) # intial value for the equilibrium shift
  if (y0 < 10^-10 && y0 > -10^-10){
    y0 <- 10^-8 # avoiding 0 mean, which makes gradient singular
  }
  
  # calculate peaks
  peaks_list <- calc_peaks(y_val, resol, y0, timen)
  peaks <- peaks_list$peaks
  peaks_time <- peaks_list$peaks_time
  
  # calc starting amplitude
  if (length(peaks) >0){
    a0 <- abs(peaks[1]) # y0 already removed
  } else {
    a0 <- max(y_val,na.rm = TRUE) - y0
  }
  
  # intial value for gamma
  if (length(peaks)==0){ # if there are no peaks, we account for that
    gam0 <- 0
  } else if (which.max(peaks)==1){ # if the highest peak is first, then damping is likely
    if (length(peaks)>1){
      # trying to figure out gamma based on logarithmic decrement
      n <- 1
      # n <- peaks_time[2]-peaks_time[1]
      log_dec <- (1/n)*log(abs(peaks[1]/peaks[2]))
      gam0 <- 1/(sqrt(1+((2*pi/log_dec)^2)))/2
    } else{ # use a guess
      gam0 <- .01
    }
  } else{ # otherwise driving is likely
    if (length(peaks)>1){
      # trying to figure out gamma based on logarithmic decrement
      n <- 1 
      # n <- peaks_time[2]-peaks_time[1]
      log_dec <- (1/n)*log(abs(peaks[2]/peaks[1]))
      gam0 <- -1*1/(sqrt(1+((2*pi/log_dec)^2)))/2
    } else{ # use a guess
      gam0 <- -.01
    }
  }
  
  if (gam0 < -.15){
    gam0 <- -.15
  } else if (gam0 > .15){
    gam0 <- .15
  }
  
  # let frequency depend on amount of peaks = (length(timen)*resol/(no of peaks+1 [accounts for phase shift])
  if (length(peaks) == 0){
    if (high == -Inf || low == Inf){
      w0 <- 2*pi/(length(timen)*resol/2)
    } else{
      # want to get their actual integer period values
      highfix <- (high/2/pi)^-1
      lowfix <- (low/2/pi)^-1
      w0 <- 2*pi/(length(timen)*resol/((highfix+lowfix)/2))
    }
  } else if (length(peaks) == 1){ # phase shift causes only one peak to appear
    w0 <- 2*pi/(length(timen)*resol/(length(peaks)+1))
  } else{
    w0 <- 2*pi/(length(timen)*resol/(length(peaks)))
  }
  
  # can't be outside the specified parameters
  if (w0 > low){
    w0 <- low
  } else if (w0 < high){
    w0 <- high
  }
  
  
  # we estimate our phase shift on the second and third nonmissing points for accuracy
  # if you have less than 3 points nonmissing, I have no hope for you
  # beginning
  # second <- which(!is.na(y_val))[2]
  # third <- which(!is.na(y_val))[3]
  beg <- which(!is.na(y_val))[2]
  # middle
  mid <- which(!is.na(y_val))[floor(length(timen)/2)]
  # end
  en <- which(!is.na(y_val))[length(timen)-1]
  
  # higher fidelity guess for higher amounts of points
  # this still works for 3 points, 
  beg1 <- which(!is.na(y_val))[2]
  beg2 <- which(!is.na(y_val))[3]
  # middle
  mid1 <- which(!is.na(y_val))[floor(length(timen)/2)]
  mid2 <- which(!is.na(y_val))[floor(length(timen)/2)+1]
  # end
  en1 <- which(!is.na(y_val))[length(timen)-2]
  en2 <- which(!is.na(y_val))[length(timen)-1]
  
  min_i <- 0;
  min_vect <- rep(Inf, length(0:11))
  for (i in 0:11){
    # if the phase shift at both the second and third gene values are the smallest value available for the fitted value
    # min_vect[i+1] <- sum(abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[second])-y_val[second]),
    #                      abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[third])-y_val[third]))
    
    # min_vect[i+1] <- sum(abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[beg])-y_val[beg]),
    #                      abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[mid])-y_val[mid]),
    #                      abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[en])-y_val[en]))
    min_vect[i+1] <- sum(abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[beg1])-y_val[beg1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[mid1])-y_val[mid1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[en1])-y_val[en1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[beg2])-y_val[beg2]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[mid2])-y_val[mid2]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[en2])-y_val[en2]))
  }
  # phi0 <- min_i*pi/6 # intial value for phase shift
  phi0 <- (which.min(abs(min_vect))-1)*pi/6 # intial value for phase shift
  
  start_param <- list(gam=gam0,a=a0,omega=w0,phi=phi0,y_shift=y0)
  
  return(start_param)
}

calc_start_echo_lin <- function(y_val, slo0, y_shift0){
  # subtract baseline
  y_val <- y_val - lin_mod(slo0,y_shift0,timen)
  
  # smooth y_val
  y_val <- smooth_y_val(y_val)
  
  # calculate starting y_shift
  y0 <- mean(y_val,na.rm = TRUE) # intial value for the equilibrium shift
  if (y0 < 10^-10 && y0 > -10^-10){
    y0 <- 10^-8 # avoiding 0 mean, which makes gradient singular
  }
  
  # calculate peaks
  peaks_list <- calc_peaks(y_val, resol, y0, timen)
  peaks <- peaks_list$peaks
  peaks_time <- peaks_list$peaks_time
  
  # calc starting amplitude
  if (length(peaks) >0){
    a0 <- abs(peaks[1]) # y0 already removed from peaks
  } else {
    a0 <- max(y_val,na.rm = TRUE) - y0
  }
  
  # intial value for gamma
  if (length(peaks)==0){ # if there are no peaks, we account for that
    gam0 <- 0
  } else if (which.max(peaks)==1){ # if the highest peak is first, then damping is likely
    if (length(peaks)>1){
      # trying to figure out gamma based on logarithmic decrement
      n <- 1
      # n <- peaks_time[2]-peaks_time[1]
      log_dec <- (1/n)*log(abs(peaks[1]/peaks[2]))
      gam0 <- 1/(sqrt(1+((2*pi/log_dec)^2)))/2
    } else{ # use a guess
      gam0 <- .01
    }
  } else{ # otherwise driving is likely
    if (length(peaks)>1){
      # trying to figure out gamma based on logarithmic decrement
      n <- 1 
      # n <- peaks_time[2]-peaks_time[1]
      log_dec <- (1/n)*log(abs(peaks[2]/peaks[1]))
      gam0 <- -1*1/(sqrt(1+((2*pi/log_dec)^2)))/2
    } else{ # use a guess
      gam0 <- -.01
    }
  }
  
  if (gam0 < -.15){
    gam0 <- -.15
  } else if (gam0 > .15){
    gam0 <- .15
  }
  
  # let frequency depend on amount of peaks = (length(timen)*resol/(no of peaks+1 [accounts for phase shift])
  if (length(peaks) == 0){
    if (high == -Inf || low == Inf){
      w0 <- 2*pi/(length(timen)*resol/2)
    } else{
      # want to get their actual integer period values
      highfix <- (high/2/pi)^-1
      lowfix <- (low/2/pi)^-1
      w0 <- 2*pi/(length(timen)*resol/((highfix+lowfix)/2))
    }
  } else if (length(peaks) == 1){ # phase shift causes only one peak to appear
    w0 <- 2*pi/(length(timen)*resol/(length(peaks)+1))
  } else{
    w0 <- 2*pi/(length(timen)*resol/(length(peaks)))
  }
  
  # can't be outside the specified parameters
  if (w0 > low){
    w0 <- low
  } else if (w0 < high){
    w0 <- high
  }
  
  
  # we estimate our phase shift on the second and third nonmissing points for accuracy
  # if you have less than 3 points nonmissing, I have no hope for you
  # beginning
  # second <- which(!is.na(y_val))[2]
  # third <- which(!is.na(y_val))[3]
  beg <- which(!is.na(y_val))[2]
  # middle
  mid <- which(!is.na(y_val))[floor(length(timen)/2)]
  # end
  en <- which(!is.na(y_val))[length(timen)-1]
  
  # higher fidelity guess for higher amounts of points
  # this still works for 3 points, 
  beg1 <- which(!is.na(y_val))[2]
  beg2 <- which(!is.na(y_val))[3]
  # middle
  mid1 <- which(!is.na(y_val))[floor(length(timen)/2)]
  mid2 <- which(!is.na(y_val))[floor(length(timen)/2)+1]
  # end
  en1 <- which(!is.na(y_val))[length(timen)-2]
  en2 <- which(!is.na(y_val))[length(timen)-1]
  
  min_i <- 0;
  min_vect <- rep(10000, length(0:11))
  for (i in 0:11){
    # if the phase shift at both the second and third gene values are the smallest value available for the fitted value
    # min_vect[i+1] <- sum(abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[second])-y_val[second]),
    #                      abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[third])-y_val[third]))
    
    # min_vect[i+1] <- sum(abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[beg])-y_val[beg]),
    #                      abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[mid])-y_val[mid]),
    #                      abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[en])-y_val[en]))
    
    min_vect[i+1] <- sum(abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[beg1])-y_val[beg1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[mid1])-y_val[mid1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[en1])-y_val[en1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[beg2])-y_val[beg2]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[mid2])-y_val[mid2]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[en2])-y_val[en2]))
  }
  # phi0 <- min_i*pi/6 # intial value for phase shift
  phi0 <- (which.min(abs(min_vect))-1)*pi/6 # intial value for phase shift
  
  start_param <- list(gam=gam0,a=a0,omega=w0,phi=phi0,slo=slo0,y_shift=y_shift0)
  
  return(start_param)
}

get_type_gam <- function(gam){
  if (gam < -.15){
    type_gam <- "Overexpressed"
  } else if (gam <= -.03){
    type_gam <- "Forced"
  } else if (gam <= .03){
    type_gam <- "Harmonic"
  } else if (gam <= .15){
    type_gam <- "Damped"
  } else{
    type_gam <- "Repressed"
  }
  return(type_gam)
}

# get hours shifted
get_hs <- function(a,phi,omega){
  # calculating the phase shift in terms of period (omega inverse of period)
  if (!is.na(a)){ # all param will either be na or not
    if (a >= 0){ # positive amplitudes
      if (phi > 0){ # shift to the left
        frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
        dist_peak <- frac_part*(2*pi/omega) # distance from first peak
        phase_hours <- (2*pi/omega)-dist_peak
      } else { # shift to the right
        frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
        dist_peak <- frac_part*(2*pi/omega) # distance from first peak
        phase_hours <- abs(dist_peak)
      }
    } else { # negative ampltitudes
      if (phi > 0){ # shift to the left
        frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
        dist_peak <- frac_part*(2*pi/omega) # distance from first peak
        if (abs(frac_part) < .5){
          phase_hours <- (2*pi/omega)-dist_peak - (2*pi/omega/2)
        } else {
          phase_hours <- (2*pi/omega)-dist_peak + (2*pi/omega/2)
        }
      } else { # shift to the right
        frac_part <- (phi/(2*pi)) - trunc(phi/(2*pi))
        dist_peak <- frac_part*(2*pi/omega) # distance from first peak
        if (abs(frac_part) < .5){
          phase_hours <- abs(dist_peak) + (2*pi/omega/2)
        } else {
          phase_hours <- abs(dist_peak) - (2*pi/omega/2)
        }
      }
    }
  } else {
    phase_hours <- NA
  }
  
  return(phase_hours)
}

# calculate echo joint with inequality constr
calc_param_echo_joint_constr <- function(current_gene, timen, resol, num_reps, start_type, res){
  gene_n <- as.character(genes_rna[current_gene,1]) # gene name
  
  # calculating averages for initial values
  if (num_reps == 1){
    y_val_rna <- as.numeric(genes_rna[current_gene,-1])
    y_val_pro <- as.numeric(genes_pro[current_gene,-1]) # all the gene values
  } else{
    y_val_pro <- avg_genes_pro[current_gene,] # starting values determined by average of replicates
    y_val_rna <- avg_genes_rna[current_gene,] # starting values determined by average of replicates
  }
  
  if (start_type == "orig"){
    rna_start <- calc_start_echo(y_val_rna)
    pro_start <- calc_start_echo(y_val_pro)
  } else {
    if (!is.na(rna_tot$Amplitude.Change.Coefficient[current_gene]) 
        & rna_tot$`BH Adj P-Value`[current_gene] < .0005
        & rna_tot$`Oscillation Type` != "Overexpressed"
        & rna_tot$`Oscillation Type` != "Repressed"){
      rna_start <- list("gam" = rna_tot$Amplitude.Change.Coefficient[current_gene],
                        "a" = rna_tot$Initial.Amplitude[current_gene],
                        "omega" = rna_tot$Radian.Frequency[current_gene],
                        "phi" = rna_tot$`Phase Shift`[current_gene],
                        "y_shift" = rna_tot$`Equilibrium Value`[current_gene])
    } else {
      rna_start <- calc_start_echo(y_val_rna)
    }
    
    
    if (!is.na(pro_tot$Amplitude.Change.Coefficient[current_gene]) &
        pro_tot$`BH Adj P-Value`[current_gene] < .0005
        & rna_tot$`Oscillation Type` != "Overexpressed"
        & rna_tot$`Oscillation Type` != "Repressed"){
      pro_start <- list("gam" = pro_tot$Amplitude.Change.Coefficient[current_gene],
                        "a" = pro_tot$Initial.Amplitude[current_gene],
                        "omega" = pro_tot$Radian.Frequency[current_gene],
                        "phi" = pro_tot$`Phase Shift`[current_gene],
                        "y_shift" = pro_tot$`Equilibrium Value`[current_gene])
    } else {
      pro_start <- calc_start_echo(y_val_pro)
    }
  }
  # we're going to use the starting value of rna for the joint period, respective for everything else
  # now we need to format it into a way that is nice
  names(rna_start) <- paste0(names(rna_start),"1")
  names(pro_start) <- paste0(names(pro_start),"2")
  start_param <- c(rna_start,pro_start)
  # start_param[["y_shift1"]] <- unlist(alph_rna[current_gene])
  # start_param[["y_shift2"]] <- unlist(alph_pro[current_gene])
  # start_param[["slo1"]] <- unlist(beta_rna[current_gene])
  # start_param[["slo2"]] <- unlist(beta_pro[current_gene])
  names(start_param)[names(start_param)=="omega2"] <- "chang"
  start_param$omega1 <- 2*pi/start_param$omega1
  start_param$chang <- 0
  
  # i am trying something
  start_joint <- start_param
  start_joint[names(start_joint)!="chang"] <- as.list(res[,names(start_joint)[names(start_joint)!="chang"]])
  start_joint$omega1 <- 2*pi/start_joint$omega1
  
  y_stack <- c(unlist(genes_rna[current_gene,-1]), unlist(genes_pro[current_gene,-1]))
  t_stack <- c(rep(timen, each=num_reps), (rep(timen, each = num_reps)))
  if (num_reps > 1){
    weights <- c(rep(calc_weights_boot(unlist(genes_rna[current_gene,-1]), num_reps), each = num_reps),
                 rep(calc_weights_boot(unlist(genes_pro[current_gene,-1]), num_reps), each = num_reps))
  }
  
  t_rep <- length(rep(timen, each=num_reps))
  # ASSUME NO MISSING DATA ATM
  s1 <- rep(c(1,0), c(t_rep,t_rep))
  s2 <- rep(c(0,1), c(t_rep,t_rep))
  
  # CONSTRAINTS ARE IN THE ORDER OF THE START PARAM
  lower_constr <- upper_constr <- c(rep(1, length(start_param)))
  names(lower_constr) <- names(upper_constr) <- names(start_param)
  # specific constr
  # lower_constr <- c( -.3,-(3*min(abs(y_stack))), high, -4*pi, min(y_stack[1:(t_rep)]),
  #                    -.3,-(3*min(abs(y_stack))), -4*pi, min(y_stack[((t_rep)+1):length(y_stack)]))
  # upper_constr <- c( .3, (3*max(abs(y_stack))), low, 4*pi, max(y_stack[1:t_rep]),
  #                    .3, (3*min(abs(y_stack))), 4*pi, max(y_stack[((t_rep)+1):length(y_stack)]))
  # back to normal constraints
  lower_constr <- c( -.15,-Inf, 2*pi/low, -Inf, min(y_stack[1:(t_rep)]),
                     -.15,-Inf, (-resol+1e-12), -Inf, min(y_stack[((t_rep)+1):length(y_stack)]))
  upper_constr <- c( .15, Inf, 2*pi/high, Inf, max(y_stack[1:t_rep]),
                     .15, Inf, resol-1e-12, Inf, max(y_stack[((t_rep)+1):length(y_stack)]))
  
  temp <- data.frame()
  coeff <- c()
  best <- ""
  num_iter <- 1000
  if (num_reps == 1){ # one replicate
    #put the timen and data point into a data frame
    temp <- data.frame(y = y_stack, t =t_stack, s1 = s1, s2 = s2, resol = rep(resol,length(y_stack)))
    temp <- temp[!is.na(temp$y),] # remove any missing data points-NONE IN EXAMPLE
    
    # USING NLS.LM
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_one_rep,
                              # jac = fcn.jac,
                              fcall = echo_joint_res_constr,
                              # jcall = jac_echo_joint,
                              lower = lower_constr,
                              upper = upper_constr,
                              t = temp$t,
                              y = temp$y,
                              s1 = temp$s1,
                              s2 = temp$s2,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff.orig <- coeff <- coef(oscillator_init) #extract parameter estimates
    iter.orig <- oscillator_init$niter
    
    oscillator.fit.orig <- nls_edit(
      formula = y ~ echo_joint_res_constr(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, y_shift1,y_shift2, t, s1, s2, resol),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], y_shift1=coeff[5],
                   gam2=coeff[6], a2=coeff[7], chang = coeff[8], phi2 = coeff[9], y_shift2=coeff[10]
      ),
      control = nls.control(maxiter = 0)
    )
    
    # USING NLS.LM
    oscillator_init <- nls.lm(par = start_joint,
                              fn = fcn_one_rep,
                              # jac = fcn.jac,
                              fcall = echo_joint_res_constr,
                              # jcall = jac_echo_joint,
                              lower = lower_constr,
                              upper = upper_constr,
                              t = temp$t,
                              y = temp$y,
                              s1 = temp$s1,
                              s2 = temp$s2,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff.joint <- coeff <- coef(oscillator_init) #extract parameter estimates
    iter.joint <- oscillator_init$niter
    
    oscillator.fit.joint <- nls_edit(
      formula = y ~ echo_joint_res_constr(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, y_shift1,y_shift2, t, s1, s2, resol),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], y_shift1=coeff[5],
                   gam2=coeff[6], a2=coeff[7], chang = coeff[8], phi2 = coeff[9], y_shift2=coeff[10]
      ),
      control = nls.control(maxiter = 0)
    )
    
  } else{ # multiple replicates
    #put the timen and data point into a data frame
    temp <- data.frame(y = y_stack, t =t_stack, w= weights, s1 = s1, s2 = s2, resol = rep(resol,length(y_stack)))
    temp <- temp[!is.na(temp$y),] # remove any missing data points-NONE IN EXAMPLE
    
    # USING NLS.LM
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn,
                              # jac = fcn.jac,
                              fcall = echo_joint_res_constr,
                              # jcall = jac_echo_joint,
                              lower = lower_constr,
                              upper = upper_constr,
                              t = temp$t,
                              y = temp$y,
                              s1 = temp$s1,
                              s2 = temp$s2,
                              w = temp$w,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff.orig <- coeff <- coef(oscillator_init) #extract parameter estimates
    iter.orig <- oscillator_init$niter
    
    # USING COBYLA
    
    # oscillator_init <-suppressMessages(cobyla(x0=as.numeric(start_param),
    #               fn=obj_constr_fun,
    #               # eval_grad_f=eval_grad_f0,
    #               lower =as.numeric(lower_constr),
    #               upper =as.numeric(upper_constr),
    #               hin =echo_joint_constr,
    #               # eval_jac_g_ineq =eval_jac_g0,
    #               s1 = temp$s1,
    #               s2 = temp$s2,
    #               t = temp$t,
    #               y = temp$y,
    #               w = temp$w,
    #               resol = resol,
    #               control = list(maxeval = 1000)))
    
    # if (oscillator_init$convergence == 5){
    #   print(oscillator_init$message)
    # }
    
    # coeff <- oscillator_init$par
    
    oscillator.fit.orig <- nls_edit(
      formula = y ~ echo_joint_res_constr(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, y_shift1,y_shift2, t, s1, s2, resol),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], y_shift1=coeff[5],
                   gam2=coeff[6], a2=coeff[7], chang = coeff[8], phi2 = coeff[9], y_shift2=coeff[10]
      ),
      control = nls.control(maxiter = 0)
    )
    
    # USING NLS.LM
    oscillator_init <- nls.lm(par = start_joint,
                              fn = fcn,
                              # jac = fcn.jac,
                              fcall = echo_joint_res_constr,
                              # jcall = jac_echo_joint,
                              lower = lower_constr,
                              upper = upper_constr,
                              t = temp$t,
                              y = temp$y,
                              s1 = temp$s1,
                              s2 = temp$s2,
                              w = temp$w,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff.joint <- coeff <- coef(oscillator_init) #extract parameter estimates
    iter.joint <- oscillator_init$niter
    
    oscillator.fit.joint <- nls_edit(
      formula = y ~ echo_joint_res_constr(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, y_shift1,y_shift2, t, s1, s2, resol),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], y_shift1=coeff[5],
                   gam2=coeff[6], a2=coeff[7], chang = coeff[8], phi2 = coeff[9], y_shift2=coeff[10]
      ),
      control = nls.control(maxiter = 0)
    )
    
  }
  
  aic_res <- AIC(oscillator.fit.orig, oscillator.fit.joint)
  ord_start <- c("orig","joint")
  oscillator.fit <- list(oscillator.fit.orig, oscillator.fit.joint)[[which.min(aic_res$AIC)]]
  best <- ord_start[which.min(aic_res$AIC)]
  coeff <- list(coeff.orig, coeff.joint)[[which.min(aic_res$AIC)]]
  num_iter <- c(iter.orig, iter.joint)[which.min(aic_res$AIC)]
  
  # did_conv <- oscillator.fit$convInfo$isConv # did the fit converge
  # num_iter <- oscillator_init$niter # amount of iterations
  # num_iter <- oscillator_init$iter # amount of iterations
  
  res <- data.frame("gene_n"=gene_n,
                    "gam1" =  coeff[1],
                    "a1" =  coeff[2],
                    "omega1" =  coeff[3],
                    "per1" =  coeff[3],
                    "phi1" =  coeff[4],
                    "y_shift1" = coeff[5],
                    "gam2" =  coeff[6],
                    "a2" =  coeff[7],
                    "chang" =  coeff[8],
                    "per2" =  (coeff[3]+coeff[8]),
                    "phi2" =  coeff[9],
                    "y_shift2" = coeff[10],
                    "iter" = num_iter,
                    "best" = best,
                    stringsAsFactors = F)
  
  # save data for the final dataframe
  # get oscillation type
  type_gam_rna <- get_type_gam(coeff[1])
  type_gam_pro <- get_type_gam(coeff[6])
  
  # get hours shifted
  phase_hours_rna <- get_hs(coeff[2],coeff[4],2*pi/coeff[3])
  phase_hours_pro <- get_hs(coeff[7],coeff[9],2*pi/(coeff[3]+coeff[8]))
  
  fin <- c("a1" =  coeff[2],
           "gam1" =  coeff[1],
           "osc1" = type_gam_rna,
           "omega1" =  2*pi/coeff[3],
           "per1" =  coeff[3],
           "phi1" =  coeff[4],
           "hs1" = phase_hours_rna,
           "y_shift1" = coeff[5],
           "a2" =  coeff[7],
           "gam2" =  coeff[6],
           "osc2" = type_gam_pro,
           "omega2" =  2*pi/(coeff[3]+coeff[8]),
           "per2" =  (coeff[3]+coeff[8]),
           "phi2" =  coeff[9],
           "hs2" = phase_hours_pro,
           "y_shift2" = coeff[10])
  
  # adjust names correctly
  echo_param_n <- c("Initial_Amplitude", "AC_Coefficient", "Oscillation_Type", "Radian_Frequency", "Period", "Phase_Shift", "Hours_Shifted", "Equilibrium_Value")
  echo_joint_param_n <- c(paste0(echo_param_n, "_RNA"), paste0(echo_param_n,"_Protein"))
  
  names(fin) <- echo_joint_param_n
  
  # generating p-values
  # joint pvalue
  all_pred <- oscillator.fit$m$predict() # fitted values
  testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
  joint_pval <- testing$p.value
  
  # rna pvalue
  rna_val <- genes_rna[current_gene,-1]
  rna_val <- rna_val[!is.na(rna_val)]
  testing <- suppressWarnings(cor.test(all_pred[1:length(rna_val)],rna_val, method = "kendall"))
  rna_pval <- testing$p.value
  
  # pro pvalue
  pro_val <- genes_pro[current_gene,-1]
  pro_val <- pro_val[!is.na(pro_val)]
  testing <- suppressWarnings(cor.test(all_pred[(length(rna_val)+1):length(all_pred)],pro_val, method = "kendall"))
  pro_pval <- testing$p.value
  
  
  res <- list("mod_name" = "echo_joint",
              "model" = oscillator.fit,
              "param" = res,
              "joint_pval" = joint_pval,
              "rna_pval" = rna_pval,
              "pro_pval" = pro_pval,
              "final" = fin,
              "type_gam_rna" = type_gam_rna,
              "type_gam_pro" = type_gam_pro)
  
  return(res)
}

# calculate echo lin joint with inequality constr - CHECK
calc_param_echo_lin_joint_constr <- function(current_gene, timen, resol, num_reps, start_type, res){
  gene_n <- as.character(genes_rna[current_gene,1]) # gene name
  
  # calculating averages for initial values
  if (num_reps == 1){
    y_val_rna <- as.numeric(genes_rna[current_gene,-1])
    y_val_pro <- as.numeric(genes_pro[current_gene,-1]) # all the gene values
  } else{
    y_val_pro <- avg_genes_pro[current_gene,] # starting values determined by average of replicates
    y_val_rna <- avg_genes_rna[current_gene,] # starting values determined by average of replicates
  }
  
  
  # calculating linear model to remove for starting values
  # calculate linear model parameters
  lin.fit <- lm(as.numeric(genes_rna[current_gene,-1])~(rep(timen, each = num_reps)))
  coeff <- as.numeric(lin.fit$coefficients)
  y_shift0_rna <- coeff[1]
  slo0_rna <- coeff[2]
  
  lin.fit <- lm(as.numeric(genes_pro[current_gene,-1])~(rep(timen, each = num_reps)))
  coeff <- as.numeric(lin.fit$coefficients)
  y_shift0_pro <- coeff[1]
  slo0_pro <- coeff[2]
  
  rna_start <- calc_start_echo_lin(y_val_rna, slo0_rna, y_shift0_rna)
  pro_start <- calc_start_echo_lin(y_val_pro, slo0_pro, y_shift0_pro)
  # we're going to use the starting value of rna for the joint period, respective for everything else
  # now we need to format it into a way that is nice
  names(rna_start) <- paste0(names(rna_start),"1")
  names(pro_start) <- paste0(names(pro_start),"2")
  start_param <- c(rna_start,pro_start)
  
  names(start_param)[names(start_param)=="omega2"] <- "chang"
  start_param$omega1 <- 2*pi/start_param$omega1
  start_param$chang <- 0
  
  # i am trying something
  start_joint <- start_param
  start_joint[names(start_joint)!="chang"] <- as.list(res[,names(start_joint)[names(start_joint)!="chang"]])
  start_joint$omega1 <- 2*pi/start_joint$omega1
  
  
  
  y_stack <- c(unlist(genes_rna[current_gene,-1]), unlist(genes_pro[current_gene,-1]))
  t_stack <- c(rep(timen, each=num_reps), (rep(timen, each = num_reps)))
  if (num_reps > 1){
    weights <- c(rep(calc_weights_boot(unlist(genes_rna[current_gene,-1]), num_reps), each = num_reps),
                 rep(calc_weights_boot(unlist(genes_pro[current_gene,-1]), num_reps), each = num_reps))
  }
  t_rep <- length(rep(timen, each=num_reps))
  # ASSUME NO MISSING DATA ATM
  s1 <- rep(c(1,0), c(t_rep,t_rep))
  s2 <- rep(c(0,1), c(t_rep,t_rep))
  
  # CONSTRAINTS ARE IN THE ORDER OF THE START PARAM
  lower_constr <- upper_constr <- c(rep(1, length(start_param)))
  names(lower_constr) <- names(upper_constr) <- names(start_param)
  # specific constr
  # lower_constr <- c( -.3,-(3*min(abs(y_stack))), high, -4*pi, min(y_stack[1:(t_rep)]),
  #                    -.3,-(3*min(abs(y_stack))), -4*pi, min(y_stack[((t_rep)+1):length(y_stack)]))
  # upper_constr <- c( .3, (3*max(abs(y_stack))), low, 4*pi, max(y_stack[1:t_rep]),
  #                    .3, (3*min(abs(y_stack))), 4*pi, max(y_stack[((t_rep)+1):length(y_stack)]))
  # back to normal constraints
  lower_constr <- c( -.15,-Inf, 2*pi/low, -Inf, -Inf, min(y_stack[1:(t_rep)]),
                     -.15,-Inf, (-resol+1e-12), -Inf, -Inf, min(y_stack[((t_rep)+1):length(y_stack)]))
  upper_constr <- c( .15, Inf, 2*pi/high, Inf, Inf, max(y_stack[1:t_rep]),
                     .15, Inf, resol-1e-12, Inf, Inf, max(y_stack[((t_rep)+1):length(y_stack)]))
  
  temp <- data.frame()
  coeff <- c()
  best <- ""
  num_iter <- 1000
  if (num_reps == 1){ # one replicate
    #put the timen and data point into a data frame
    temp <- data.frame(y = y_stack, t =t_stack, s1 = s1, s2 = s2)
    temp <- temp[!is.na(temp$y),] # remove any missing data points-NONE IN EXAMPLE
    
    # USING NLS.LM
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_one_rep,
                              # jac = fcn.jac,
                              fcall = echo_lin_joint_res_constr,
                              # jcall = jac_echo_joint,
                              lower = lower_constr,
                              upper = upper_constr,
                              t = temp$t,
                              y = temp$y,
                              s1 = temp$s1,
                              s2 = temp$s2,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff.orig <- coeff <- coef(oscillator_init) #extract parameter estimates
    iter.orig <- oscillator_init$niter
    
    oscillator.fit.orig <- nls_edit(
      formula = y ~ echo_lin_joint_res_constr(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, slo1,slo2, y_shift1,y_shift2, t, s1, s2, resol),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], slo1 = coeff[5], y_shift1=coeff[6],
                   gam2=coeff[7], a2=coeff[8], chang = coeff[9], phi2 = coeff[10],  slo2 = coeff[11], y_shift2=coeff[12]
      ),
      control = nls.control(maxiter = 0)
    )
    
    # USING NLS.LM
    oscillator_init <- nls.lm(par = start_joint,
                              fn = fcn_one_rep,
                              # jac = fcn.jac,
                              fcall = echo_lin_joint_res_constr,
                              # jcall = jac_echo_joint,
                              lower = lower_constr,
                              upper = upper_constr,
                              t = temp$t,
                              y = temp$y,
                              s1 = temp$s1,
                              s2 = temp$s2,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff.joint <- coeff <- coef(oscillator_init) #extract parameter estimates
    iter.joint <- oscillator_init$niter
    
    oscillator.fit.joint <- nls_edit(
      formula = y ~ echo_lin_joint_res_constr(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, slo1,slo2, y_shift1,y_shift2, t, s1, s2, resol),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], slo1 = coeff[5], y_shift1=coeff[6],
                   gam2=coeff[7], a2=coeff[8], chang = coeff[9], phi2 = coeff[10],  slo2 = coeff[11], y_shift2=coeff[12]
      ),
      control = nls.control(maxiter = 0)
    )
    
  } else{ # multiple replicates
    #put the timen and data point into a data frame
    temp <- data.frame(y = y_stack, t =t_stack, w= weights, s1 = s1, s2 = s2)
    temp <- temp[!is.na(temp$y),] # remove any missing data points-NONE IN EXAMPLE
    
    # USING NLS.LM
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn,
                              # jac = fcn.jac,
                              fcall = echo_lin_joint_res_constr,
                              # jcall = jac_echo_joint,
                              lower = lower_constr,
                              upper = upper_constr,
                              t = temp$t,
                              y = temp$y,
                              s1 = temp$s1,
                              s2 = temp$s2,
                              w = temp$w,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff.orig <- coeff <- coef(oscillator_init) #extract parameter estimates
    iter.orig <- oscillator_init$niter
    
    oscillator.fit.orig <- nls_edit(
      formula = y ~ echo_lin_joint_res_constr(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, slo1,slo2, y_shift1,y_shift2, t, s1, s2, resol),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], slo1 = coeff[5], y_shift1=coeff[6],
                   gam2=coeff[7], a2=coeff[8], chang = coeff[9], phi2 = coeff[10],  slo2 = coeff[11], y_shift2=coeff[12]
      ),
      control = nls.control(maxiter = 0)
    )
    
    # USING NLS.LM
    oscillator_init <- nls.lm(par = start_joint,
                              fn = fcn,
                              # jac = fcn.jac,
                              fcall = echo_lin_joint_res_constr,
                              # jcall = jac_echo_joint,
                              lower = lower_constr,
                              upper = upper_constr,
                              t = temp$t,
                              y = temp$y,
                              s1 = temp$s1,
                              s2 = temp$s2,
                              w = temp$w,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff.joint <- coeff <- coef(oscillator_init) #extract parameter estimates
    iter.joint <- oscillator_init$niter
    
    oscillator.fit.joint <- nls_edit(
      formula = y ~ echo_lin_joint_res_constr(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, slo1,slo2, y_shift1,y_shift2, t, s1, s2, resol),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], slo1 = coeff[5], y_shift1=coeff[6],
                   gam2=coeff[7], a2=coeff[8], chang = coeff[9], phi2 = coeff[10],  slo2 = coeff[11], y_shift2=coeff[12]
      ),
      control = nls.control(maxiter = 0)
    )
  }
  
  aic_res <- AIC(oscillator.fit.orig, oscillator.fit.joint)
  ord_start <- c("orig","joint")
  oscillator.fit <- list(oscillator.fit.orig, oscillator.fit.joint)[[which.min(aic_res$AIC)]]
  best <- ord_start[which.min(aic_res$AIC)]
  coeff <- list(coeff.orig, coeff.joint)[[which.min(aic_res$AIC)]]
  num_iter <- c(iter.orig, iter.joint)[which.min(aic_res$AIC)]
  
  # did_conv <- oscillator.fit$convInfo$isConv # did the fit converge
  # num_iter <- oscillator_init$niter # amount of iterations
  # num_iter <- oscillator_init$iter # amount of iterations
  
  res <- data.frame("gene_n"=gene_n,
                    "gam1" =  coeff[1],
                    "a1" =  coeff[2],
                    "omega1" =  coeff[3],
                    "per1" =  coeff[3],
                    "phi1" =  coeff[4],
                    "slo1" = coeff[5],
                    "y_shift1" = coeff[6],
                    "gam2" =  coeff[7],
                    "a2" =  coeff[8],
                    "chang" =  coeff[9],
                    "per2" =  (coeff[3]+coeff[9]),
                    "phi2" =  coeff[10],
                    "slo2" = coeff[11],
                    "y_shift2" = coeff[12],
                    "iter" = num_iter,
                    "best" = best,
                    stringsAsFactors = F)
  
  # save data for the final dataframe
  # get oscillation type
  type_gam_rna <- get_type_gam(coeff[1])
  type_gam_pro <- get_type_gam(coeff[7])
  
  # get hours shifted
  phase_hours_rna <- get_hs(coeff[2],coeff[4],2*pi/coeff[3])
  phase_hours_pro <- get_hs(coeff[8],coeff[10],2*pi/(coeff[3]+coeff[9]))
  
  fin <- c("a1" =  coeff[2],
           "gam1" =  coeff[1],
           "osc1" = type_gam_rna,
           "omega1" =  2*pi/coeff[3],
           "per1" =  coeff[3],
           "phi1" =  coeff[4],
           "hs1" = phase_hours_rna,
           "slo1" = coeff[5],
           "y_shift1" = coeff[6],
           "a2" =  coeff[8],
           "gam2" =  coeff[7],
           "osc2" = type_gam_pro,
           "omega2" =  2*pi/(coeff[3]+coeff[9]),
           "per2" =  (coeff[3]+coeff[9]),
           "phi2" =  coeff[10],
           "hs2" = phase_hours_pro,
           "slo2" = coeff[11],
           "y_shift2" = coeff[12])
  
  echo_lin_param_n <- c("Initial_Amplitude", "AC_Coefficient", "Oscillation_Type", "Radian_Frequency", "Period", "Phase_Shift", "Hours_Shifted", "Slope", "Equilibrium_Value")
  echo_lin_joint_param_n <- c(paste0(echo_lin_param_n, "_RNA"), paste0(echo_lin_param_n,"_Protein"))
  names(fin) <- echo_lin_joint_param_n
  
  # generating p-values
  # joint pvalue
  all_pred <- oscillator.fit$m$predict() # fitted values
  testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
  joint_pval <- testing$p.value
  
  # rna pvalue
  rna_val <- genes_rna[current_gene,-1]
  rna_val <- rna_val[!is.na(rna_val)]
  testing <- suppressWarnings(cor.test(all_pred[1:length(rna_val)],rna_val, method = "kendall"))
  rna_pval <- testing$p.value
  
  # pro pvalue
  pro_val <- genes_pro[current_gene,-1]
  pro_val <- pro_val[!is.na(pro_val)]
  testing <- suppressWarnings(cor.test(all_pred[(length(rna_val)+1):length(all_pred)],pro_val, method = "kendall"))
  pro_pval <- testing$p.value
  
  
  res <- list("mod_name" = "echo_lin_joint",
              "model" = oscillator.fit,
              "param" = res,
              "joint_pval" = joint_pval,
              "rna_pval" = rna_pval,
              "pro_pval" = pro_pval,
              "final" = fin,
              "type_gam_rna" = type_gam_rna,
              "type_gam_pro" = type_gam_pro)
  
  return(res)
}


# calculate echo joint parameters - getting initial starting points
calc_param_echo_joint <- function(current_gene, timen, resol, num_reps, start_type){
  gene_n <- as.character(genes_rna[current_gene,1]) # gene name
  
  # calculating averages for initial values
  if (num_reps == 1){
    y_val_rna <- as.numeric(genes_rna[current_gene,-1])
    y_val_pro <- as.numeric(genes_pro[current_gene,-1]) # all the gene values
  } else{
    y_val_pro <- avg_genes_pro[current_gene,] # starting values determined by average of replicates
    y_val_rna <- avg_genes_rna[current_gene,] # starting values determined by average of replicates
  }
  
  if (start_type == "orig"){
    rna_start <- calc_start_echo(y_val_rna)
    pro_start <- calc_start_echo(y_val_pro)
  } else {
    if (!is.na(rna_tot$Amplitude.Change.Coefficient[current_gene]) 
        & rna_tot$`BH Adj P-Value`[current_gene] < .0005
        & rna_tot$`Oscillation Type` != "Overexpressed"
        & rna_tot$`Oscillation Type` != "Repressed"){
      rna_start <- list("gam" = rna_tot$Amplitude.Change.Coefficient[current_gene],
                        "a" = rna_tot$Initial.Amplitude[current_gene],
                        "omega" = rna_tot$Radian.Frequency[current_gene],
                        "phi" = rna_tot$`Phase Shift`[current_gene],
                        "y_shift" = rna_tot$`Equilibrium Value`[current_gene])
    } else {
      rna_start <- calc_start_echo(y_val_rna)
    }
    
    
    if (!is.na(pro_tot$Amplitude.Change.Coefficient[current_gene]) &
        pro_tot$`BH Adj P-Value`[current_gene] < .0005
        & rna_tot$`Oscillation Type` != "Overexpressed"
        & rna_tot$`Oscillation Type` != "Repressed"){
      pro_start <- list("gam" = pro_tot$Amplitude.Change.Coefficient[current_gene],
                        "a" = pro_tot$Initial.Amplitude[current_gene],
                        "omega" = pro_tot$Radian.Frequency[current_gene],
                        "phi" = pro_tot$`Phase Shift`[current_gene],
                        "y_shift" = pro_tot$`Equilibrium Value`[current_gene])
    } else {
      pro_start <- calc_start_echo(y_val_pro)
    }
  }
  # we're going to use the starting value of rna for the joint period, respective for everything else
  # now we need to format it into a way that is nice
  names(rna_start) <- paste0(names(rna_start),"1")
  names(pro_start) <- paste0(names(pro_start),"2")
  start_param <- c(rna_start,pro_start)
  # start_param[["y_shift1"]] <- unlist(alph_rna[current_gene])
  # start_param[["y_shift2"]] <- unlist(alph_pro[current_gene])
  # start_param[["slo1"]] <- unlist(beta_rna[current_gene])
  # start_param[["slo2"]] <- unlist(beta_pro[current_gene])
  start_param$omega2 <- NULL
  
  y_stack <- c(unlist(genes_rna[current_gene,-1]), unlist(genes_pro[current_gene,-1]))
  t_stack <- c(rep(timen, each=num_reps), (rep(timen, each = num_reps)))
  if (num_reps > 1){
    weights <- c(rep(calc_weights_boot(unlist(genes_rna[current_gene,-1]), num_reps), each = num_reps),
                 rep(calc_weights_boot(unlist(genes_pro[current_gene,-1]), num_reps), each = num_reps))
  }
  
  t_rep <- length(rep(timen, each=num_reps))
  # ASSUME NO MISSING DATA ATM
  s1 <- rep(c(1,0), c(t_rep,t_rep))
  s2 <- rep(c(0,1), c(t_rep,t_rep))
  
  # CONSTRAINTS ARE IN THE ORDER OF THE START PARAM
  lower_constr <- upper_constr <- c(rep(1, length(start_param)))
  names(lower_constr) <- names(upper_constr) <- names(start_param)
  # specific constr
  # lower_constr <- c( -.3,-(3*min(abs(y_stack))), high, -4*pi, min(y_stack[1:(t_rep)]),
  #                    -.3,-(3*min(abs(y_stack))), -4*pi, min(y_stack[((t_rep)+1):length(y_stack)]))
  # upper_constr <- c( .3, (3*max(abs(y_stack))), low, 4*pi, max(y_stack[1:t_rep]),
  #                    .3, (3*min(abs(y_stack))), 4*pi, max(y_stack[((t_rep)+1):length(y_stack)]))
  # back to normal constraints
  lower_constr <- c( -.15,-Inf, high, -Inf, min(y_stack[1:(t_rep)]),
                     -.15,-Inf, -Inf, min(y_stack[((t_rep)+1):length(y_stack)]))
  upper_constr <- c( .15, Inf, low, Inf, max(y_stack[1:t_rep]),
                     .15, Inf, Inf, max(y_stack[((t_rep)+1):length(y_stack)]))
  
  temp <- data.frame()
  coeff <- c()
  if (num_reps == 1){ # one replicate
    #put the timen and data point into a data frame
    temp <- data.frame(y = y_stack, t =t_stack, s1 = s1, s2 = s2)
    temp <- temp[!is.na(temp$y),] # remove any missing data points-NONE IN EXAMPLE
    
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_one_rep,
                              # jac = fcn.jac,
                              fcall = echo_joint,
                              # jcall = jac_echo_joint,
                              lower = lower_constr,
                              upper = upper_constr,
                              t = temp$t,
                              y = temp$y,
                              s1 = temp$s1,
                              s2 = temp$s2,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff <- coef(oscillator_init) #extract parameter estimates
    
    oscillator.fit <- nls_edit(
      formula = y ~ echo_joint(a1,a2, gam1,gam2, omega1, phi1,phi2, y_shift1,y_shift2, t, s1, s2),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], y_shift1=coeff[5],
                   gam2=coeff[6], a2=coeff[7], phi2 = coeff[8], y_shift2=coeff[9]
      ),
      control = nls.control(maxiter = 0)
    )
    
  } else{ # multiple replicates
    #put the timen and data point into a data frame
    temp <- data.frame(y = y_stack, t =t_stack, w= weights, s1 = s1, s2 = s2)
    temp <- temp[!is.na(temp$y),] # remove any missing data points-NONE IN EXAMPLE
    
    oscillator_init <- nls.lm(par = start_param,
                             fn = fcn,
                             # jac = fcn.jac,
                             fcall = echo_joint,
                             # jcall = jac_echo_joint,
                             lower = lower_constr,
                             upper = upper_constr,
                             t = temp$t,
                             y = temp$y,
                             s1 = temp$s1,
                             s2 = temp$s2,
                             w = temp$w,
                             control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                      ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff <- coef(oscillator_init) #extract parameter estimates
    
    oscillator.fit <- nls_edit(
      formula = y ~ echo_joint(a1,a2, gam1,gam2, omega1, phi1,phi2, y_shift1,y_shift2, t, s1, s2),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], y_shift1=coeff[5],
                   gam2=coeff[6], a2=coeff[7], phi2 = coeff[8], y_shift2=coeff[9]
      ),
      control = nls.control(maxiter = 0)
    )
  }
  
  # did_conv <- oscillator.fit$convInfo$isConv # did the fit converge
  num_iter <- oscillator_init$niter # amount of iterations
  
  res <- data.frame("gene_n"=gene_n,
                    "gam1" =  coeff[1],
                    "a1" =  coeff[2],
                    "omega1" =  coeff[3],
                    "phi1" =  coeff[4],
                    "y_shift1" = coeff[5],
                    "gam2" =  coeff[6],
                    "a2" =  coeff[7],
                    "phi2" =  coeff[8],
                    "y_shift2" = coeff[9],
                    stringsAsFactors = F)
  
  res <- list("mod_name" = "echo_joint",
              "model" = oscillator.fit,
              "param" = res)
  return(res)
}

# calculate echo lin joint parameters - getting initial starting points
calc_param_echo_lin_joint <- function(current_gene, timen, resol, num_reps, start_type){
  gene_n <- as.character(genes_rna[current_gene,1]) # gene name
  
  # calculating averages for initial values
  if (num_reps == 1){
    y_val_rna <- as.numeric(genes_rna[current_gene,-1])
    y_val_pro <- as.numeric(genes_pro[current_gene,-1]) # all the gene values
  } else{
    y_val_pro <- avg_genes_pro[current_gene,] # starting values determined by average of replicates
    y_val_rna <- avg_genes_rna[current_gene,] # starting values determined by average of replicates
  }
  
  # calculating linear model to remove for starting values
  # calculate linear model parameters
  lin.fit <- lm(as.numeric(genes_rna[current_gene,-1])~(rep(timen, each = num_reps)))
  coeff <- as.numeric(lin.fit$coefficients)
  y_shift0_rna <- coeff[1]
  slo0_rna <- coeff[2]
  
  lin.fit <- lm(as.numeric(genes_pro[current_gene,-1])~(rep(timen, each = num_reps)))
  coeff <- as.numeric(lin.fit$coefficients)
  y_shift0_pro <- coeff[1]
  slo0_pro <- coeff[2]
  
  rna_start <- calc_start_echo_lin(y_val_rna, slo0_rna, y_shift0_rna)
  pro_start <- calc_start_echo_lin(y_val_pro, slo0_pro, y_shift0_pro)
  
  # we're going to use the starting value of rna for the joint period, respective for everything else
  # now we need to format it into a way that is nice
  names(rna_start) <- paste0(names(rna_start),"1")
  names(pro_start) <- paste0(names(pro_start),"2")
  start_param <- c(rna_start,pro_start)
  start_param$omega2 <- NULL
  
  y_stack <- c(unlist(genes_rna[current_gene,-1]), unlist(genes_pro[current_gene,-1]))
  t_stack <- c(rep(timen, each=num_reps), (rep(timen, each = num_reps)))
  if (num_reps > 1){
    weights <- c(rep(calc_weights_boot(unlist(genes_rna[current_gene,-1]), num_reps), each = num_reps),
                 rep(calc_weights_boot(unlist(genes_pro[current_gene,-1]), num_reps), each = num_reps))
  }
  t_rep <- length(rep(timen, each=num_reps))
  # ASSUME NO MISSING DATA ATM
  s1 <- rep(c(1,0), c(t_rep,t_rep))
  s2 <- rep(c(0,1), c(t_rep,t_rep))
  
  # CONSTRAINTS ARE IN THE ORDER OF THE START PARAM
  lower_constr <- upper_constr <- c(rep(1, length(start_param)))
  names(lower_constr) <- names(upper_constr) <- names(start_param)
  # specific constr
  # lower_constr <- c( -.3,-(3*min(abs(y_stack))), high, -4*pi, min(y_stack[1:(t_rep)]),
  #                    -.3,-(3*min(abs(y_stack))), -4*pi, min(y_stack[((t_rep)+1):length(y_stack)]))
  # upper_constr <- c( .3, (3*max(abs(y_stack))), low, 4*pi, max(y_stack[1:t_rep]),
  #                    .3, (3*min(abs(y_stack))), 4*pi, max(y_stack[((t_rep)+1):length(y_stack)]))
  # back to normal constraints
  lower_constr <- c( -.15,-Inf, high, -Inf, -Inf,  min(y_stack[1:(t_rep)]),
                     -.15,-Inf, -Inf, -Inf, min(y_stack[((t_rep)+1):length(y_stack)]))
  upper_constr <- c( .15, Inf, low, Inf, Inf, max(y_stack[1:t_rep]),
                     .15, Inf, Inf, Inf, max(y_stack[((t_rep)+1):length(y_stack)]))
  
  temp <- data.frame()
  coeff <- c()
  if (num_reps == 1){ # one replicate
    #put the timen and data point into a data frame
    temp <- data.frame(y = y_stack, t =t_stack, s1 = s1, s2 = s2)
    temp <- temp[!is.na(temp$y),] # remove any missing data points-NONE IN EXAMPLE
    
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_one_rep,
                              # jac = fcn.jac,
                              fcall = echo_lin_joint,
                              # jcall = jac_echo_joint,
                              lower = lower_constr,
                              upper = upper_constr,
                              t = temp$t,
                              y = temp$y,
                              s1 = temp$s1,
                              s2 = temp$s2,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff <- coef(oscillator_init) #extract parameter estimates
    
    oscillator.fit <- nls_edit(
      formula = y ~ echo_lin_joint(a1,a2, gam1,gam2, omega1, phi1,phi2, slo1,slo2, y_shift1,y_shift2, t, s1, s2),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], slo1 = coeff[5], y_shift1=coeff[6],
                   gam2=coeff[7], a2=coeff[8], phi2 = coeff[9], slo2 = coeff[10], y_shift2=coeff[11]
      ),
      control = nls.control(maxiter = 0)
    )
  } else{ # multiple replicates
    #put the timen and data point into a data frame
    temp <- data.frame(y = y_stack, t =t_stack, w= weights, s1 = s1, s2 = s2)
    temp <- temp[!is.na(temp$y),] # remove any missing data points-NONE IN EXAMPLE
    
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn,
                              # jac = fcn.jac,
                              fcall = echo_lin_joint,
                              # jcall = jac_echo_joint,
                              lower = lower_constr,
                              upper = upper_constr,
                              t = temp$t,
                              y = temp$y,
                              s1 = temp$s1,
                              s2 = temp$s2,
                              w = temp$w,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff <- coef(oscillator_init) #extract parameter estimates
    
    oscillator.fit <- nls_edit(
      formula = y ~ echo_lin_joint(a1,a2, gam1,gam2, omega1, phi1,phi2, slo1,slo2, y_shift1,y_shift2, t, s1, s2),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], slo1 = coeff[5], y_shift1=coeff[6],
                   gam2=coeff[7], a2=coeff[8], phi2 = coeff[9], slo2 = coeff[10], y_shift2=coeff[11]
      ),
      control = nls.control(maxiter = 0)
    )
  }
  
  # did_conv <- oscillator.fit$convInfo$isConv # did the fit converge
  num_iter <- oscillator_init$niter # amount of iterations
  
  res <- data.frame("gene_n"=gene_n,
                    "gam1" =  coeff[1],
                    "a1" =  coeff[2],
                    "omega1" =  coeff[3],
                    "phi1" =  coeff[4],
                    "slo1" = coeff[5],
                    "y_shift1" = coeff[6],
                    "gam2" =  coeff[7],
                    "a2" =  coeff[8],
                    "slo2" = coeff[9],
                    "phi2" =  coeff[10],
                    "y_shift2" = coeff[11],
                    stringsAsFactors = F)
  
  res <- list("mod_name" = "echo_lin_joint",
              "model" = oscillator.fit,
              "param" = res)
  return(res)
}

# calculate echo parameters
calc_param_echo <- function(current_gene, timen, resol, num_reps, genes, avg_genes){

  gene_n <- as.character(genes[current_gene,1]) # gene name
  
  # calculating averages for initial values
  if (num_reps == 1){
    y_val <- as.numeric(genes[current_gene,-1]) # all the gene values
  } else{
    y_val <- avg_genes[current_gene,] # starting values determined by average of replicates
  }
  
  start_param <- calc_start_echo(y_val)
  temp <- data.frame()
  if (num_reps == 1){ # one replicate
    # put the timen into a data frame
    temp <- data.frame(y=y_val,t=timen)
    temp <- temp[!is.na(temp$y),] # remove any missing data points
    
    # fitting
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_single_one_rep,
                              fcall = alt_form,
                              lower=c(-.15, -Inf, high, -Inf, min(temp$y)),
                              upper=c(.15, Inf, low, Inf, max(temp$y)),
                              t = temp$t,
                              y = temp$y,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff <- coef(oscillator_init)
    
    oscillator.fit <- nls_edit(
      formula = y ~ alt_form(a, gam, omega, phi, y_shift, t),
      data=temp,
      start = list(gam=coeff[1], a=coeff[2], omega = coeff[3], phi = coeff[4], y_shift=coeff[5]),
      control = nls.control(maxiter = 0)
    )
    
  } else{ # multiple replicates
    #put the timen and data point into a data frame
    weights <- calc_weights(current_gene,num_reps, genes)
    temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),t=cbind(rep(timen,each = num_reps)),w=cbind(rep(weights,each = num_reps)))
    temp <- temp[!is.na(temp$y),] # remove any missing data points
    
    
    #fitting
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_single,
                              fcall = alt_form,
                              lower=c(-.15, -Inf, high, -Inf, min(temp$y)),
                              upper=c(.15, Inf, low, Inf, max(temp$y)),
                              t = temp$t,
                              y = temp$y,
                              w = temp$w,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff <- coef(oscillator_init)
    
    oscillator.fit <- nls_edit(
      formula = y ~ alt_form(a, gam, omega, phi, y_shift, t),
      data=temp,
      start = list(gam=coeff[1], a=coeff[2], omega = coeff[3], phi = coeff[4], y_shift=coeff[5]),
      control = nls.control(maxiter = 0)
    )
  }
  
  did_conv <- oscillator.fit$convInfo$isConv # did the fit converge
  num_iter <- oscillator.fit$convInfo$finIter # amount of iterations
  
  coeff <- coef(oscillator.fit) #extract parameter estimates
  
  # alt_form parameters:
  # the parameters go in the order of: gam,a,omega,phi,y_shift
  
  res <- data.frame("gene_n"=gene_n,
                    "gam" =  coeff[1],
                    "a" =  coeff[2],
                    "omega" =  coeff[3],
                    "phi" =  coeff[4],
                    "y_shift" = coeff[5],
                    stringsAsFactors = F)
  
  # save data for the final dataframe
  # get oscillation type
  type_gam <- get_type_gam(coeff[1])
  
  phase_hours <- get_hs(coeff[2],coeff[4],coeff[3])
  
  fin <- c(coeff[2], coeff[1], type_gam, coeff[3], 2*pi/coeff[3], coeff[4], phase_hours, coeff[5])
  echo_param_n <- c("Initial_Amplitude", "AC_Coefficient", "Oscillation_Type", "Radian_Frequency", "Period", "Phase_Shift", "Hours_Shifted", "Equilibrium_Value")
  names(fin) <- echo_param_n
  
  # generating p-values
  all_pred <- oscillator.fit$m$predict() # fitted values
  testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
  pval <- testing$p.value
  
  res <- list("mod_name" = "echo",
              "model" = oscillator.fit,
              "param" = res,
              "pval" = pval,
              "final" = fin,
              "type_gam" = type_gam)
  
  return(res)
}

# calculate echo_lin parameters
calc_param_echo_lin <- function(current_gene, timen, resol, num_reps, genes, avg_genes){
  
  gene_n <- as.character(genes[current_gene,1]) # gene name
  
  # calculating linear model to remove for starting values
  # calculate linear model parameters
  lin.fit <- lm(as.numeric(genes[current_gene,-1])~(rep(timen, each = num_reps)))
  coeff <- lin.fit$coefficients
  y_shift0 <- coeff[1]
  slo0 <- coeff[2]
  
  # calculating averages for initial values
  if (num_reps == 1){
    y_val <- as.numeric(genes[current_gene,-1]) # all the gene values
  } else{
    y_val <- avg_genes[current_gene,] # starting values determined by average of replicates
  }
  
  # get starting parameters
  start_param <- calc_start_echo_lin(y_val, slo0, y_shift0)
  temp <- data.frame()
  if (num_reps == 1){ # one replicate
    # put the timen into a data frame
    temp <- data.frame(y=y_val,t=timen)
    temp <- temp[!is.na(temp$y),] # remove any missing data points
    
    # fitting
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_single_one_rep,
                              fcall = echo_lin_mod,
                              lower=c(-.15, -Inf, high, -Inf, -Inf, min(temp$y)),
                              upper=c(.15, Inf, low, Inf, Inf, max(temp$y)),
                              t = temp$t,
                              y = temp$y,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff <- coef(oscillator_init)
    
    oscillator.fit <- nls_edit(
      formula = y ~ echo_lin_mod(a, gam, omega, phi, slo, y_shift, t),
      data=temp,
      start = list(gam=coeff[1], a=coeff[2], omega = coeff[3], phi = coeff[4], slo = coeff[5], y_shift=coeff[6]),
      control = nls.control(maxiter = 0)
    )
    
  } else{ # multiple replicates
    #put the timen and data point into a data frame
    weights <- calc_weights(current_gene,num_reps, genes)
    temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),t=cbind(rep(timen,each = num_reps)),w=cbind(rep(weights,each = num_reps)))
    temp <- temp[!is.na(temp$y),] # remove any missing data points
    
    
    #fitting
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_single,
                              fcall = echo_lin_mod,
                              lower=c(-.15, -Inf, high, -Inf, -Inf, min(temp$y)),
                              upper=c(.15, Inf, low, Inf, Inf, max(temp$y)),
                              t = temp$t,
                              y = temp$y,
                              w = temp$w,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    coeff <- coef(oscillator_init)
    
    oscillator.fit <- nls_edit(
      formula = y ~ echo_lin_mod(a, gam, omega, phi, slo, y_shift, t),
      data=temp,
      start = list(gam=coeff[1], a=coeff[2], omega = coeff[3], phi = coeff[4], slo = coeff[5], y_shift=coeff[6]),
      control = nls.control(maxiter = 0)
    )
  }
  
  did_conv <- oscillator.fit$convInfo$isConv # did the fit converge
  num_iter <- oscillator.fit$convInfo$finIter # amount of iterations
  
  coeff <- coef(oscillator.fit) #extract parameter estimates
  
  # alt_form parameters:
  # the parameters go in the order of: gam,a,omega,phi,y_shift
  
  res <- data.frame("gene_n"=gene_n,
                    "gam" =  coeff[1],
                    "a" =  coeff[2],
                    "omega" =  coeff[3],
                    "phi" =  coeff[4],
                    "slo" = coeff[5],
                    "y_shift" = coeff[6],
                    stringsAsFactors = F)
  
  # save data for the final dataframe
  # get oscillation type
  gam <- coeff[1]
  type_gam <- get_type_gam(gam)
  
  a <- coeff[2]
  phi <- coeff[4]
  omega <- coeff[3]
  phase_hours <- get_hs(a,phi,omega)
  
  fin <- c(coeff[2], coeff[1], type_gam, coeff[3], 2*pi/coeff[3], coeff[4], phase_hours, coeff[5], coeff[6])
  names(fin) <- echo_lin_param_n <- c("Initial_Amplitude", "AC_Coefficient", "Oscillation_Type", "Radian_Frequency", "Period", "Phase_Shift", "Hours_Shifted", "Slope", "Equilibrium_Value")
  
  # generating p-values
  all_pred <- oscillator.fit$m$predict() # fitted values
  testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
  pval <- testing$p.value
  
  res <- list("mod_name" = "echo_lin",
              "model" = oscillator.fit,
              "param" = res,
              "pval" = pval,
              "final" = fin,
              "type_gam" = type_gam)
  
  return(res)
}

# calculate linear model
calc_param_lin <- function(current_gene, timen, resol, num_reps, genes, avg_genes){
  # calculate linear model parameters
  temp <- data.frame(y=as.numeric(genes[current_gene,-1]),t=(rep(timen, each = num_reps)))
  temp <- temp[!is.na(temp$y),] # remove any missing data points
  
  lin.fit <- lm(temp$y~temp$t)
  coeff <- lin.fit$coefficients
  y_shift <- coeff[1]
  slo <- coeff[2]
  
  res <- as.data.frame(t(c(slo,y_shift)))
  
  gene_n <- genes[current_gene,1]
  res <- cbind(data.frame("gene_n"=gene_n, stringsAsFactors = F), res)
  
  fin <- c(coeff[2], coeff[1])
  names(fin) <- lin_param_n <- c("Slope", "Equilibrium_Value")
  
  # generating p-values
  all_pred <- lin.fit$fitted.values # fitted values
  testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
  pval <- testing$p.value
  
  # generate slope p-value
  sumfit <- summary(lin.fit)
  slo_pval <- sumfit$coefficients[2,4]
  
  res <- list("mod_name" = "lin",
              "model" = lin.fit,
              "param" = res,
              "pval" = pval,
              "slo_pval" = slo_pval,
              "final" = fin)
  
  return(res)
}

# claculate constant model - DO NOT USE
calc_param_const <- function(current_gene, timen, resol, num_reps, genes, avg_genes){
  gene_n <- genes[current_gene,1]
  # calculate constant model
  # result is mean, but we need the lm object
  const.fit <- lm(as.numeric(genes[current_gene,-1]) ~ 1)
  coeff <- const.fit$coefficients
  res <- data.frame("gene_n"=gene_n, 
                    "y_shift" = coeff[1],
                    stringsAsFactors = F)
  
  fin <- coeff[1]
  names(fin) <- const_param_n <- c("Equilibrium_Value")
  
  # # generating p-values
  # # fitted values
  # all_pred <- const.fit$fitted.values
  # testing <- suppressWarnings(cor.test(all_pred,as.numeric(genes[current_gene,-1]), method = "kendall"))
  # pval <- testing$p.value
  
  pval <- summary(const.fit)$coefficients[[4]]
  
  res <- list("mod_name" = "const",
              "model" = const.fit,
              "param" = res,
              "pval" = pval,
              "final" = fin)
  
  return(res)
}

# calculate exponential model
calc_param_exp <- function(current_gene, timen, resol, num_reps, genes, avg_genes){
  # calculating averages for initial values
  if (num_reps == 1){
    y_val <- as.numeric(genes[current_gene,-1]) # all the gene values
  } else{
    y_val <- avg_genes[current_gene,] # starting values determined by average of replicates
  }
  
  # fitting based on actual values (assume positive amplitude)
  # shift values to get rid of y shift
  y_shift_e0 <- min(y_val,na.rm = T)
  y_edit <- y_val
  y_edit <- y_val-y_shift_e0
  # make it the next worst value, after the edit, so it doesn't skew the points (because log(0)=-Inf)
  y_edit[which.min(y_val)] <- abs(min(y_edit[-which.min(y_val)], na.rm = T))
  # do the fitting and assignment
  exp_mod_coeff <- lm(log(y_edit)~timen)$coefficients
  gam_e0 <- exp_mod_coeff[2]
  a_e0 <- exp(exp_mod_coeff[1])
  
  # fitting based on negative values
  y_shift_e0_neg <- min(-y_val,na.rm = T)
  y_edit <- -y_val
  y_edit <- -y_val-y_shift_e0_neg
  # make it the next worst value, after the edit, so it doesn't skew the points (because log(0)=-Inf)
  y_edit[which.min(-y_val)] <- abs(min(y_edit[-which.min(-y_val)], na.rm = T))
  # do the fitting and assignment
  exp_mod_coeff <- lm(log(y_edit)~timen)$coefficients
  gam_e0_neg <- exp_mod_coeff[2]
  a_e0_neg <- -exp(exp_mod_coeff[1])
  y_shift_e0_neg <- -y_shift_e0_neg
  
  
  
  if (num_reps == 1){
    temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),t=cbind(rep(timen,each = num_reps)))
    temp <- temp[!is.na(temp$y),] # remove any missing data points
    
    
    #fitting
    exp_mod_pos <- nls.lm(par = list(a=a_e0,gam=gam_e0, y_shift=y_shift_e0),
                           fn = fcn_single_one_rep,
                           fcall = exp_mod,
                           t = temp$t,
                           y = temp$y,
                           control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                    ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    exp_mod_neg <- nls.lm(par = list(a=a_e0_neg,gam=gam_e0_neg, y_shift=y_shift_e0_neg),
                          fn = fcn_single_one_rep,
                          fcall = exp_mod,
                          t = temp$t,
                          y = temp$y,
                          control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                   ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
  } else {
    weights <- calc_weights(current_gene,num_reps, genes)
    temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),t=cbind(rep(timen,each = num_reps)),w=cbind(rep(weights,each = num_reps)))
    temp <- temp[!is.na(temp$y),] # remove any missing data points
    
    #fitting
    exp_mod_pos <- nls.lm(par = list(a=a_e0,gam=gam_e0, y_shift=y_shift_e0),
                          fn = fcn_single,
                          fcall = exp_mod,
                          t = temp$t,
                          y = temp$y,
                          w = temp$w,
                          control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                   ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
    
    exp_mod_neg <- nls.lm(par = list(a=a_e0_neg,gam=gam_e0_neg, y_shift=y_shift_e0_neg),
                           fn = fcn_single,
                           fcall = exp_mod,
                           t = temp$t,
                           y = temp$y,
                           w = temp$w,
                           control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                    ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
  }
  coeff_pos <- coef(exp_mod_pos)
  coeff_neg <- coef(exp_mod_neg)
  
  exp_mod.fit_pos <- nls_edit(
    formula = y ~ exp_mod(a, gam, y_shift, t),
    data=temp,
    start = list(a=coeff_pos[1],gam=coeff_pos[2], y_shift=coeff_pos[3]),
    control = nls.control(maxiter = 0)
  )
  
  exp_mod.fit_neg <- nls_edit(
    formula = y ~ exp_mod(a, gam, y_shift, t),
    data=temp,
    start = list(a=coeff_neg[1],gam=coeff_neg[2], y_shift=coeff_neg[3]),
    control = nls.control(maxiter = 0)
  )
  
  
  aic_res <- AIC(exp_mod.fit_pos, exp_mod.fit_neg)
  exp_mod.fit<- list(exp_mod.fit_pos, exp_mod.fit_neg)[[which.min(aic_res$AIC)]]
  coeff <- list(coeff_pos, coeff_neg)[[which.min(aic_res$AIC)]]
  
  gene_n <- genes[current_gene,1]
  res <- data.frame("gene_n"=gene_n,
                    "a" = coeff[1],
                    "gam" = coeff[2],
                    "y_shift" = coeff[3],
                    stringsAsFactors = F)
  
  fin <- coeff
  names(fin) <- exp_param_n <- c("Initial_Amplitude", "Growth_Rate", "Equilibrium_Value")
  
  # generating p-values
  all_pred <- exp_mod.fit$m$predict() # fitted values
  testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
  pval <- testing$p.value
  
  
  res <- list("mod_name" = "exp",
              "model" = exp_mod.fit,
              "param" = res,
              "pval" = pval,
              "final" = fin)
  
  return(res)
}

# calculate exponential linear model - WRONG, DO NOT USE
calc_param_exp_lin <- function(current_gene, timen, resol, num_reps, genes, avg_genes){
  # calculating averages for initial values
  if (num_reps == 1){
    y_val <- as.numeric(genes[current_gene,-1]) # all the gene values
  } else{
    y_val <- avg_genes[current_gene,] # starting values determined by average of replicates
  }
  
  coeff <- lm(as.numeric(genes[current_gene,-1])~(rep(timen, each = num_reps)))$coefficients
  y_shift_e0 <- coeff[1]
  slo_e0 <- coeff[2]
  
  
  y_edit <- y_lin_rem  <- y_val-lin_mod(slo_e0, y_shift_e0, timen)
  poss_shift <- (min(y_edit, na.rm = T))
  # need to shift everything for the exponential, as per ushe, if negative
  y_edit <- y_edit-poss_shift
  # make it the next worst value, after the edit, so it doesn't skew the points (because log(0)=-Inf)
  y_edit[which.min(y_lin_rem)] <- abs(min(y_edit[-which.min(y_lin_rem)], na.rm = T))
  y_shift_e0 <- y_shift_e0 + poss_shift
  
  exp_mod_coeff <- lm(log(y_edit)~timen)$coefficients
  gam_e0 <- exp_mod_coeff[2]
  a_e0 <- exp(exp_mod_coeff[1])
  
  
  # fitting based on actual values (assume positive amplitude)
  # shift values to get rid of y shift
  y_shift_e0 <- min(y_val,na.rm = T)
  y_edit <- y_val
  y_edit <- y_val-y_shift_e0
  # make it the next worst value, after the edit, so it doesn't skew the points (because log(0)=-Inf)
  y_edit[which.min(y_val)] <- abs(min(y_edit[-which.min(y_val)], na.rm = T))
  # do the fitting and assignment
  exp_mod_coeff <- lm(log(y_edit)~timen)$coefficients
  gam_e0 <- exp_mod_coeff[2]
  a_e0 <- exp(exp_mod_coeff[1])
  
  # fitting based on negative values
  y_shift_e0_neg <- min(-y_val,na.rm = T)
  y_edit <- -y_val
  y_edit <- -y_val-y_shift_e0_neg
  # make it the next worst value, after the edit, so it doesn't skew the points (because log(0)=-Inf)
  y_edit[which.min(-y_val)] <- abs(min(y_edit[-which.min(-y_val)], na.rm = T))
  # do the fitting and assignment
  exp_mod_coeff <- lm(log(y_edit)~timen)$coefficients
  gam_e0_neg <- exp_mod_coeff[2]
  a_e0_neg <- -exp(exp_mod_coeff[1])
  y_shift_e0_neg <- -y_shift_e0_neg
  
  
  weights <- calc_weights(current_gene,num_reps, genes)
  temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),t=cbind(rep(timen,each = num_reps)),w=cbind(rep(weights,each = num_reps)))
  temp <- temp[!is.na(temp$y),] # remove any missing data points
  
  #fitting
  
  exp_lin_init <- nls.lm(par = list(a=a_e0,gam=gam_e0,slo=slo_e0, y_shift=y_shift_e0),
                           fn = fcn_single,
                           # jac = fcn.jac_else,
                           fcall = exp_lin_mod,
                           # jcall = jac_exp_lin_joint,
                           # lower = lower_constr,
                           # upper = upper_constr,
                           t = temp$t,
                           y = temp$y,
                           w = temp$w,
                           control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                    ftol=1e-6, ptol=1e-6, gtol=1e-6)
  )
  
  coeff <- coef(exp_lin_init)
  
  exp_lin.fit <- nls_edit(
    formula = y ~ exp_lin_mod(a, gam, slo, y_shift, t),
    data=temp,
    start = list(a=coeff[1],gam=coeff[2],slo=coeff[3], y_shift=coeff[4]),
    control = nls.control(maxiter = 0, warnOnly = T)
  )
  
  gene_n <- genes[current_gene,1]
  res <- data.frame("gene_n"=gene_n,
                    "a" = coeff[1],
                    "gam" = coeff[2],
                    "slo" = coeff[3],
                    "y_shift" = coeff[4],
                    stringsAsFactors = F)
  
  fin <- coeff
  names(fin) <- exp_lin_param_n <- c("Initial_Amplitude", "Growth_Rate", "Slope", "Equilibrium_Value")
  
  # generating p-values
  all_pred <- exp_lin.fit$m$predict() # fitted values
  testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
  pval <- testing$p.value
  
  res <- list("mod_name" = "exp_lin",
              "model" = exp_lin.fit,
              "param" = res,
              "pval" = pval,
              "final" = fin)
  
  return(res)
}

# ORIGINAL FUNCTION --- DO NOT USE
choose_best_mod_ORIG <- function(timen, resol, num_reps){
  model_n <- c("const", "lin", "exp", "exp_lin", "echo", "echo_lin")
  aic_res <- data.frame(matrix(0, nrow = 0, ncol = 13))
  colnames(aic_res) <- c("best_mod_orig", "const_val", "lin_val", "exp_val", "exp_lin_val", "echo_val", "echo_lin_val", "echo_joint_val", "echo_lin_joint_val", "free_joint_val", "best_mod_after", "best_mod_pval_orig", "best_mod_pval_after")
  # all_const_val <- all_lin_val <- all_exp_val <- all_exp_lin_val <- all_echo_val <- all_echo_lin_val <- list()  
  best_mods <- list()
  for (current_gene in 1:nrow(genes)){
    if (current_gene%%100 == 0 ){
      print(current_gene)
    }
    
    # calculate models
    # all_const_val[[current_gene]] <- const_val <- calc_param_const(current_gene, timen, resol, num_reps)
    # all_lin_val[[current_gene]] <- lin_val <- calc_param_lin(current_gene, timen, resol, num_reps)
    # all_exp_val[[current_gene]] <- exp_val <- calc_param_exp(current_gene, timen, resol, num_reps)
    # all_exp_lin_val[[current_gene]] <- exp_lin_val <- calc_param_exp_lin(current_gene, timen, resol, num_reps)
    # all_echo_val[[current_gene]] <- echo_val <- calc_param_echo(current_gene, timen, resol, num_reps)
    # all_echo_lin_val[[current_gene]] <- echo_lin_val <- calc_param_echo_lin(current_gene, timen, resol, num_reps)
    
    const_val <- calc_param_const(current_gene, timen, resol, num_reps)
    lin_val <- calc_param_lin(current_gene, timen, resol, num_reps)
    exp_val <- calc_param_exp(current_gene, timen, resol, num_reps)
    exp_lin_val <- calc_param_exp_lin(current_gene, timen, resol, num_reps)
    echo_val <- calc_param_echo(current_gene, timen, resol, num_reps)
    echo_lin_val <- calc_param_echo_lin(current_gene, timen, resol, num_reps)
    
    # calculate AICs
    aic_vals <- AIC(const_val$model, lin_val$model, exp_val$model, exp_lin_val$model, echo_val$model, echo_lin_val$model)
    pvals <- c(const_val$pval, lin_val$pval, exp_val$pval, exp_lin_val$pval, echo_val$pval, echo_lin_val$pval)
    
    # save results
    best_mods[[current_gene]] <- list(const_val, lin_val, exp_val, exp_lin_val, echo_val, echo_lin_val)[[which.min(aic_vals$AIC)]]
    
    aic_df <- data.frame(rbind(c(0, aic_vals$AIC, NA, NA, NA, NA, 0, 0)))
    aic_df[1,1] <- aic_df[1,11] <- model_n[which.min(aic_vals$AIC)]
    aic_df[1,12] <- aic_df[1,13] <- pvals[which.min(aic_vals$AIC)]
    colnames(aic_df) <- colnames(aic_res)
    aic_res <- rbind(aic_res, aic_df)
  }
  
  # i think we're going to use the aic, but i will doublecheck
  table(aic_res$best_mod)
  
  all_mod <- list("best_mods" = best_mods,
                  # "const" = all_const_val,
                  # "lin" = all_lin_val,
                  # "exp" = all_exp_val,
                  # "exp_lin" = all_exp_lin_val,
                  # "echo" = all_echo_val,
                  # "echo_lin" = all_echo_lin_val,
                  "aic_res" = aic_res
                  )
  
  return(all_mod)
}

# HAS EXP LIN FUNCTION --- DO NOT USE
choose_best_mod_W_EXP_LIN <- function(current_gene, timen, resol, num_reps, genes, avg_genes){
  model_n <- c("const", "lin", "exp", "exp_lin", "echo", "echo_lin")
  
  const_val <- calc_param_const(current_gene, timen, resol, num_reps, genes, avg_genes)
  lin_val <- calc_param_lin(current_gene, timen, resol, num_reps, genes, avg_genes)
  exp_val <- calc_param_exp(current_gene, timen, resol, num_reps, genes, avg_genes)
  exp_lin_val <- calc_param_exp_lin(current_gene, timen, resol, num_reps, genes, avg_genes)
  echo_val <- calc_param_echo(current_gene, timen, resol, num_reps, genes, avg_genes)
  echo_lin_val <- calc_param_echo_lin(current_gene, timen, resol, num_reps, genes, avg_genes)
  
  # calculate AICs
  aic_vals <- AIC(const_val$model, lin_val$model, exp_val$model, exp_lin_val$model, echo_val$model, echo_lin_val$model)
  pvals <- c(const_val$pval, lin_val$pval, exp_val$pval, exp_lin_val$pval, echo_val$pval, echo_lin_val$pval)
  
  # save results
  best_mod <- list(const_val, lin_val, exp_val, exp_lin_val, echo_val, echo_lin_val)[[which.min(aic_vals$AIC)]]
  
  aic_df <- data.frame(rbind(c(0, aic_vals$AIC, NA, NA, NA, NA, 0, 0, NA, NA)))
  aic_df[1,1] <- aic_df[1,11] <- model_n[which.min(aic_vals$AIC)]
  aic_df[1,12] <- aic_df[1,13] <- pvals[which.min(aic_vals$AIC)]
  aic_df[1,14] <- aic_df[1,15] <- if (!is.null(best_mod$type_gam)) best_mod$type_gam else ""
  colnames(aic_df) <- c("best_mod_orig", "const_val", "lin_val", "exp_val", "exp_lin_val", "echo_val", "echo_lin_val", "echo_joint_val", "echo_lin_joint_val", "free_joint_val", "best_mod_after", "best_mod_pval_orig", "best_mod_pval_after", "best_mod_type_gam_orig", "best_mod_type_gam_after")
  
  all_mod <- list("best_mod" = best_mod, "aic_df" = aic_df)
  
  return(all_mod)
}

# HAS CONST FUNCTION --- DO NOT USE
choose_best_mod_W_CONST <- function(current_gene, timen, resol, num_reps, genes, avg_genes){
  model_n <- c("const", "lin", "exp", "echo", "echo_lin")
  
  const_val <- calc_param_const(current_gene, timen, resol, num_reps, genes, avg_genes)
  lin_val <- calc_param_lin(current_gene, timen, resol, num_reps, genes, avg_genes)
  exp_val <- calc_param_exp(current_gene, timen, resol, num_reps, genes, avg_genes)
  echo_val <- calc_param_echo(current_gene, timen, resol, num_reps, genes, avg_genes)
  echo_lin_val <- calc_param_echo_lin(current_gene, timen, resol, num_reps, genes, avg_genes)
  
  # calculate AICs
  aic_vals <- AIC(const_val$model, lin_val$model, exp_val$model, echo_val$model, echo_lin_val$model)
  pvals <- c(const_val$pval, lin_val$pval, exp_val$pval, echo_val$pval, echo_lin_val$pval)
  
  # save results
  best_mod <- list(const_val, lin_val, exp_val, echo_val, echo_lin_val)[[which.min(aic_vals$AIC)]]
  
  aic_df <- data.frame(rbind(c(0, aic_vals$AIC, NA, NA, NA, NA, 0, 0, NA, NA)))
  aic_df[1,1] <- aic_df[1,10] <- model_n[which.min(aic_vals$AIC)]
  aic_df[1,11] <- aic_df[1,12] <- pvals[which.min(aic_vals$AIC)]
  aic_df[1,13] <- aic_df[1,14] <- if (!is.null(best_mod$type_gam)) best_mod$type_gam else ""
  colnames(aic_df) <- c("best_mod_orig", "const_val", "lin_val", "exp_val", "echo_val", "echo_lin_val", "echo_joint_val", "echo_lin_joint_val", "free_joint_val", "best_mod_after", "best_mod_pval_orig", "best_mod_pval_after", "best_mod_type_gam_orig", "best_mod_type_gam_after")
  
  all_mod <- list("best_mod" = best_mod, "aic_df" = aic_df)
  
  return(all_mod)
}

# function to calculate corrected AIC
AICc <- function(object, ..., n){
  aics <- AIC(object,...)
  
  aicc <- sapply(1:nrow(aics), function(x){aics[x,"AIC"] + ((2*(aics$df[x]-1)*(aics$df[x]-1+1))/(n-(aics$df[x]-1)-1))})
  aics$BIC <- aicc
  
  return(aics)
}

get_mod_fits <- function(best_mod, timen, resol){
  prm <- suppressWarnings(as.numeric(best_mod$param))
  
  fit_val <-
    if (best_mod$mod_name == "echo"){
      alt_form(prm[3],prm[2],prm[4],prm[5],prm[6], timen)
    } else if (best_mod$mod_name == "echo_lin"){
      echo_lin_mod(prm[3],prm[2],prm[4],prm[5],prm[6],prm[7], timen)
    } else if (best_mod$mod_name == "exp"){
      exp_mod(prm[2],prm[3],prm[4],timen)
    } else if (best_mod$mod_name == "lin") {
      lin_mod(prm[2], prm[3], timen)
    } else {
      s1 <- rep(c(1,0), c(length(timen), length(timen)))
      s2 <- rep(c(0,1), c(length(timen), length(timen)))
      if (best_mod$mod_name == "echo_joint"){
        
        echo_joint_res_constr(
          prm[3],prm[9], # a
          prm[2],prm[8], # gam
          prm[4],prm[10], # omega
          prm[6],prm[12], # phi
          prm[7],prm[13], # yshift
          timen,s1,s2,resol)
      } else if (best_mod$mod_name == "echo_lin_joint"){
        echo_lin_joint_res_constr(
          prm[3],prm[10], # a
          prm[2],prm[9], # gam
          prm[4],prm[11], # omega
          prm[6],prm[13], # phi
          prm[7],prm[14], # slo
          prm[8],prm[15], # yshift
          timen,s1,s2,resol)
      }
    }
  
  return(fit_val)
}

choose_best_mod <- function(current_gene, timen, resol, num_reps, genes, avg_genes){
  model_n <- c("lin", "exp", "echo", "echo_lin")
  
  lin_val <- calc_param_lin(current_gene, timen, resol, num_reps, genes, avg_genes)
  exp_val <- calc_param_exp(current_gene, timen, resol, num_reps, genes, avg_genes)
  echo_val <- calc_param_echo(current_gene, timen, resol, num_reps, genes, avg_genes)
  echo_lin_val <- calc_param_echo_lin(current_gene, timen, resol, num_reps, genes, avg_genes)
  
  # calculate AICs
  # n <- length(as.numeric(genes[current_gene,-1]))
  
  aic_vals <- BIC(lin_val$model, exp_val$model, echo_val$model, echo_lin_val$model)
  pvals <- c(lin_val$pval, exp_val$pval, echo_val$pval, echo_lin_val$pval)
  
  # save results
  best_mod <- list(lin_val, exp_val, echo_val, echo_lin_val)[[which.min(aic_vals$BIC)]]
  
  aic_df <- data.frame(rbind(c(0, aic_vals$BIC, NA, NA, NA, NA, 0, 0, NA, NA)))
  aic_df[1,1] <- aic_df[1,9] <- model_n[which.min(aic_vals$BIC)]
  aic_df[1,10] <- aic_df[1,11] <- pvals[which.min(aic_vals$BIC)]
  aic_df[1,12] <- aic_df[1,13] <- if (!is.null(best_mod$type_gam)) best_mod$type_gam else ""
  colnames(aic_df) <- c("best_mod_orig", "lin_val", "exp_val", "echo_val", "echo_lin_val", "echo_joint_val", "echo_lin_joint_val", "free_joint_val", "best_mod_after", "best_mod_pval_orig", "best_mod_pval_after", "best_mod_type_gam_orig", "best_mod_type_gam_after")
  
  all_mod <- list("best_mod" = best_mod, "aic_df" = aic_df)
  
  return(all_mod)
}

# created in order to put everything in parallel -- carrying around dataframes
multi_pipeline <- function(current_gene, timen, resol, num_reps, final_df_row, genes_rna, genes_pro, rem_unexpr_combine){
  # preallocating
  final_df <- final_df_row
  
  # remove unexpressed genes
  if (rem_unexpr_combine[current_gene]){
    final_df$Best_Model_RNA <- final_df$Best_Model_Protein <- "Unexpressed"
    return(list(
      "final_df" = final_df,
      "aic_rna" = data.frame(),
      "aic_pro" = data.frame()
    ))
  }
  
  # remove genes with less than 3 datapoints
  if (sum(!is.na(genes_rna[,-1])) < 3 | sum(!is.na(genes_pro[,-1])) < 3){
    final_df$Best_Model_RNA <- final_df$Best_Model_Protein <- "Unexpressed"
    return(list(
      "final_df" = final_df,
      "aic_rna" = data.frame(),
      "aic_pro" = data.frame()
    ))
  }
  
  # get names of parameters in order
  # model name to pretty mapping
  mod_name_map <- c("echo_lin" = "ECHO Linear",
                    "echo" = "ECHO",
                    "echo_lin_joint" = "ECHO Linear Joint",
                    "echo_joint" = "ECHO Joint",
                    "exp_lin" = "Exponential Linear",
                    "exp" = "Exponential",
                    "lin" = "Linear",
                    "const" = "Constant")
  
  # parameter names
  echo_lin_param_n <- c("Initial_Amplitude", "AC_Coefficient", "Oscillation_Type", "Radian_Frequency", "Period", "Phase_Shift", "Hours_Shifted", "Slope", "Equilibrium_Value")
  echo_lin_joint_param_n <- c(paste0(echo_lin_param_n, "_RNA"), paste0(echo_lin_param_n,"_Protein"))
  echo_param_n <- c("Initial_Amplitude", "AC_Coefficient", "Oscillation_Type", "Radian_Frequency", "Period", "Phase_Shift", "Hours_Shifted", "Equilibrium_Value")
  echo_joint_param_n <- c(paste0(echo_param_n, "_RNA"), paste0(echo_param_n,"_Protein"))
  # exp_lin_param_n <- c("Initial_Amplitude", "Growth_Rate", "Slope", "Equilibrium_Value")
  exp_param_n <- c("Initial_Amplitude", "Growth_Rate", "Equilibrium_Value")
  lin_param_n <- c("Slope", "Equilibrium_Value")
  # const_param_n <- c("Equilibrium_Value")
  
  # put all the specific parameter names in a list
  param_map <- list(
    "echo_lin" = echo_lin_param_n,
    "echo_lin_joint" = echo_lin_joint_param_n,
    "echo" = echo_param_n,
    "echo_joint" = echo_joint_param_n,
    # "exp_lin" = exp_lin_param_n,
    "exp" = exp_param_n,
    "lin" = lin_param_n#,
    # "const" = const_param_n
  )
  
  # so, we start by getting the best model for RNA and protein
  # rna:
  genes <- genes_rna
  avg_genes <- avg_genes_rna
  rna_mods <- choose_best_mod(current_gene, timen, resol, num_reps, genes_rna, avg_genes_rna)
  aic_rna <- rna_mods$aic_df # save the aic
  
  # protein
  genes <- genes_pro
  avg_genes <- avg_genes_pro
  pro_mods <- choose_best_mod(current_gene, timen, resol, num_reps, genes_pro, avg_genes_pro)
  aic_pro <- pro_mods$aic_df # save the aic
  
  # now, are both models oscillatory, or is one oscillatory and one is not?
  send_joint <- which((aic_rna$best_mod_orig == "echo" | aic_pro$best_mod_orig == "echo" | 
                   aic_rna$best_mod_orig == "echo_lin" | aic_pro$best_mod_orig == "echo_lin") &
    !(aic_rna$best_mod_orig == "echo" & aic_pro$best_mod_orig == "echo") &
    !(aic_rna$best_mod_orig == "echo_lin" & aic_pro$best_mod_orig == "echo_lin") &
    !(aic_rna$best_mod_orig == "echo_lin" & aic_pro$best_mod_orig == "echo") &
    !(aic_rna$best_mod_orig == "echo" & aic_pro$best_mod_orig == "echo_lin"))
  
  if (length(send_joint) > 0){
    # first, we need to get the joint model
    
    # get what the best model is supposed to be
    
    if ("echo" %in% c(aic_rna$best_mod_orig[1], aic_pro$best_mod_orig[1])){
      # echo joint model
      ej_mod <- calc_param_echo_joint(current_gene, timen, resol, num_reps, "orig")
      res <- ej_mod$param #extract parameter estimates
      
      # echo joint model with constraints
      ejc_mod <- calc_param_echo_joint_constr(current_gene, timen, resol, num_reps, "orig", res)
    } else {
      # echo lin joint model
      ej_mod <- calc_param_echo_lin_joint(current_gene, timen, resol, num_reps, "orig")
      res <- ej_mod$param #extract parameter estimates
      
      # echo lin joint model with constraints
      ejc_mod <- calc_param_echo_lin_joint_constr(current_gene, timen, resol, num_reps, "orig", res)
    }
    
    # compare with best free joint model
    rna_mod <- rna_mods$best_mod
    pro_mod <- pro_mods$best_mod
    
    best_mod <- reeval_best_mod(ejc_mod, rna_mod, pro_mod, current_gene, timen, num_reps, resol)
    
    if (best_mod$mod_name != "free"){ # joint model is best model
      joint_mod <- best_mod$joint_mod
      # we want to substitute them in
      aic_rna$best_mod_after[1] <- aic_pro$best_mod_after[1] <- best_mod$mod_name
      aic_rna[1,paste0(best_mod$mod_name,"_val")] <- 
        aic_pro[1,paste0(best_mod$mod_name,"_val")] <- best_mod$aic_joint
      aic_rna$best_mod_pval_after[1] <- joint_mod$rna_pval
      aic_pro$best_mod_pval_after[1] <- joint_mod$pro_pval
      
      aic_rna$free_joint_val[1] <- aic_pro$free_joint_val[1] <- best_mod$aic_free
      aic_rna$best_mod_type_gam_after[1] <- joint_mod$type_gam_rna
      aic_pro$best_mod_type_gam_after[1] <- joint_mod$type_gam_pro
      
      # update the final dataframes
      
      # add the designation and the model names
      final_df$Best_Model_RNA[1] <- final_df$Best_Model_Protein[1] <- mod_name_map[[best_mod$mod_name]]
      final_df$P_Value_RNA[1] <- aic_rna$best_mod_pval_after[1]
      final_df$P_Value_Protein[1] <- aic_pro$best_mod_pval_after[1]
      final_df$P_Value_Joint[1] <- joint_mod$joint_pval
      
      # add the parameters
      final_df[1, param_map[[joint_mod$mod_name]]] <- joint_mod$final
      
      final_df[1,(33+(length(timen)*num_reps)+1):(33+(length(timen)*num_reps)+length(timen))] <- 
        get_mod_fits(joint_mod, timen, resol)[1:length(timen)]
      final_df[1,(33+(length(timen)*num_reps*2)+length(timen)+1):ncol(final_df)] <- 
        get_mod_fits(joint_mod, timen, resol)[(length(timen)+1):(length(timen)*2)]
    } else {
      aic_rna$best_mod_after[1] <- aic_pro$best_mod_after[1] <- "free"
      aic_rna[1,paste0(best_mod$mod_name,"_val")] <- 
        aic_pro[1,paste0(best_mod$mod_name,"_val")] <- best_mod$aic_joint
      aic_rna$free_joint_val[1] <- aic_pro$free_joint_val[1] <- best_mod$aic_free
      
      # add the designation and the model names
      final_df$Best_Model_RNA[1] <- mod_name_map[aic_rna$best_mod_orig[1]]
      final_df$Best_Model_Protein[1] <- mod_name_map[aic_pro$best_mod_orig[1]]
      final_df$P_Value_RNA[1] <- aic_rna$best_mod_pval_orig[1]
      final_df$P_Value_Protein[1] <- aic_pro$best_mod_pval_orig[1]
      
      # add the linear slope p-value, if it exists
      if (final_df$Best_Model_RNA[1] == "Linear"){
        final_df$P_Value_Linear_Slope_RNA[1] <- rna_mods$best_mod$slo_pval
      }
      if (final_df$Best_Model_Protein[1] == "Linear"){
        final_df$P_Value_Linear_Slope_Protein[1] <- pro_mods$best_mod$slo_pval
      }
      
      # add the parameters
      final_df[1, paste0(param_map[[aic_rna$best_mod_orig[1]]],"_RNA")] <-
        rna_mods$best_mod$final
      final_df[1, paste0(param_map[[aic_pro$best_mod_orig[1]]],"_Protein")] <-
        pro_mods$best_mod$final
      
      # add the fitted values
      final_df[1,(33+(length(timen)*num_reps)+1):(33+(length(timen)*num_reps)+length(timen))] <- 
        get_mod_fits(rna_mods$best_mod, timen, resol)
      final_df[1,(33+(length(timen)*num_reps*2)+length(timen)+1):ncol(final_df)] <- 
        get_mod_fits(pro_mods$best_mod, timen, resol)
    }
  } else {
    # add the designation and the model names
    final_df$Best_Model_RNA[1] <- mod_name_map[aic_rna$best_mod_orig[1]]
    final_df$Best_Model_Protein[1] <- mod_name_map[aic_pro$best_mod_orig[1]]
    final_df$P_Value_RNA[1] <- aic_rna$best_mod_pval_orig[1]
    final_df$P_Value_Protein[1] <- aic_pro$best_mod_pval_orig[1]
    
    # add the linear slope p-value, if it exists
    if (final_df$Best_Model_RNA[1] == "Linear"){
      final_df$P_Value_Linear_Slope_RNA[1] <- rna_mods$best_mod$slo_pval
    }
    if (final_df$Best_Model_Protein[1] == "Linear"){
      final_df$P_Value_Linear_Slope_Protein[1] <- pro_mods$best_mod$slo_pval
    }
    
    # add the parameters
    final_df[1, paste0(param_map[[aic_rna$best_mod_orig[1]]],"_RNA")] <-
      rna_mods$best_mod$final
    final_df[1, paste0(param_map[[aic_pro$best_mod_orig[1]]],"_Protein")] <-
      pro_mods$best_mod$final
    
    # add the fitted values
    final_df[1,(33+(length(timen)*num_reps)+1):(33+(length(timen)*num_reps)+length(timen))] <- 
      get_mod_fits(rna_mods$best_mod, timen, resol)
    final_df[1,(33+(length(timen)*num_reps*2)+length(timen)+1):ncol(final_df)] <- 
      get_mod_fits(pro_mods$best_mod, timen, resol)
  }
  
  fin <- list(
    "final_df" = final_df,
    "aic_rna" = aic_rna,
    "aic_pro" = aic_pro
  )
  
  return(fin)
  
}


get_res <- function(gene_n, timen, joint_tot, rna_tot, pro_tot, start_joint, start_echo){
  
  current_gene <- which(genes_rna$Gene.Name == gene_n)
  timen <- timen
  y_stack <- c(unlist(genes_rna[current_gene,-1]), unlist(genes_pro[current_gene,-1]))
  t_stack <- c(rep(timen, each=num_reps), (rep(timen, each = num_reps)))
  weights <- c(rep(calc_weights_boot(unlist(genes_rna[current_gene,-1]), num_reps), each = num_reps),
               rep(calc_weights_boot(unlist(genes_pro[current_gene,-1]), num_reps), each = num_reps))
  
  t_rep <- length(rep(timen, each=num_reps))
  # ASSUME NO MISSING DATA ATM
  s1 <- rep(c(1,0), c(t_rep,t_rep))
  s2 <- rep(c(0,1), c(t_rep,t_rep))
  
  lower_constr <- c( -Inf,-Inf, high, -Inf, min(y_stack[1:(t_rep)]),
                     -Inf,-Inf, -Inf, min(y_stack[((t_rep)+1):length(y_stack)]))
  upper_constr <- c( Inf, Inf, low, Inf, max(y_stack[1:t_rep]),
                     Inf, Inf, Inf, max(y_stack[((t_rep)+1):length(y_stack)]))
  temp <- data.frame(y = y_stack, t =t_stack, w= weights, s1 = s1, s2 = s2)
  temp <- temp[!is.na(temp$y),] # remove any missing data points-NONE IN EXAMPLE
  
  tryCatch({
    
    osc.joint <- nlsLM(formula = y ~ echo_joint(a1,a2, gam1,gam2, omega1, phi1,phi2, y_shift1,y_shift2, t, s1, s2),
                       data = temp,
                       start = start_joint,
                       control = nls.lm.control(maxiter = 0),
                       weights = temp$w,
                       lower = lower_constr,
                       upper = upper_constr)
    
    lower_constr <- c( -Inf,-Inf, high, -Inf, min(y_stack[1:(t_rep)]),
                       -Inf,-Inf, high, -Inf, min(y_stack[((t_rep)+1):length(y_stack)]))
    upper_constr <- c( Inf, Inf, low, Inf, max(y_stack[1:t_rep]),
                       Inf, Inf, low, Inf, max(y_stack[((t_rep)+1):length(y_stack)]))
    osc.full <- nlsLM(formula = y ~ echo_full(a1,a2, gam1,gam2, omega1,omega2, phi1,phi2, y_shift1,y_shift2, t, s1, s2),
                      data = temp,
                      start = start_echo,
                      control = nls.lm.control(maxiter = 0),
                      weights = temp$w,
                      lower = lower_constr,
                      upper = upper_constr)
    
    # get log likelihood
    res <- lrtest(osc.full, osc.joint)
    
    return(res)
  }, error = function(err){
    res <- list("Pr(>Chisq)"=c(NA,NA), "LogLik"=c(NA,NA))
    
    return(res)
  })
}

compare_echo <- function(joint_mod, full_mod_rna, full_mod_pro,
                         current_gene, timen, resol, num_reps){
  # set up the vector for data
  y_stack <- c(unlist(genes_rna[current_gene,-1]), unlist(genes_pro[current_gene,-1]))
  t_stack <- c(rep(timen, each=num_reps), (rep(timen, each = num_reps)))
  weights <- c(rep(calc_weights_boot(unlist(genes_rna[current_gene,-1]), num_reps), each = num_reps),
               rep(calc_weights_boot(unlist(genes_pro[current_gene,-1]), num_reps), each = num_reps))
  
  t_rep <- length(rep(timen, each=num_reps))
  s1 <- rep(c(1,0), c(t_rep,t_rep))
  s2 <- rep(c(0,1), c(t_rep,t_rep))
  
  #put the timen and data point into a data frame
  temp <- data.frame(y = y_stack, t =t_stack, w= weights, s1 = s1, s2 = s2)
  temp <- temp[!is.na(temp$y),] # remove any missing data points-NONE IN EXAMPLE
  
  
  start_joint <- joint_mod$param
  start_joint$gene_n <- NULL
  
  start_echo <- c(full_mod_rna$param, full_mod_pro$param)
  names(start_echo) <- c("rna_name","gam1", "a1","omega1","phi1","y_shift1",
                          "pro_name","gam2", "a2","omega2","phi2","y_shift2")
  start_echo$rna_name <- start_echo$pro_name <- NULL
  
  joint_mod <- calc_param_echo_joint(current_gene, timen, resol, num_reps, start_type)
  
  # getting the log likelihood
  osc.full <-nls_edit(formula = y ~ echo_full(a1,a2, gam1,gam2, omega1,omega2, phi1,phi2, y_shift1,y_shift2, t, s1, s2),
                    data = temp,
                    start = start_echo,
                    control = nls.control(maxiter = 0)
                    )
  osc.joint <- joint_mod$model
  
  osc_test <- lrtest(osc.full,osc.joint)
  
  winner <- ""
  if (osc_test$`Pr(>Chisq)`[2] < .05){
    # the two fits are significantly different -- take the one with the higher logLik
    if (logLik(osc.joint) > logLik(osc.full)){
      winner <- "joint"
    } else {
      winner <- "full"
    }
  } else {
    # the fits are basically the same, so take the parsimonious one
    winner <- "joint"
  }
  
  # get kendall's tau
  # full_fit <- echo_full(start_echo$a1, start_echo$a2, start_echo$gam1, start_echo$gam2,
  #                       start_echo$omega1, start_echo$omega2, start_echo$phi1, start_echo$phi2,
  #                       start_echo$y_shift1, start_echo$y_shift2, temp$t, temp$s1, temp$s2)
  # kenp_full <- cor.test(full_fit, temp$y, method = "kendall")$p.value
  # joint_fit <- echo_joint(start_echo$a1, start_echo$a2, start_echo$gam1, start_echo$gam2,
  #                         start_echo$omega1, start_echo$phi1, start_echo$phi2,
  #                         start_echo$y_shift1, start_echo$y_shift2, temp$t, temp$s1, temp$s2)
  # kenp_joint <- cor.test(joint_fit, temp$y, method = "kendall")$p.value
  # 
  # return(list("gene_n" = gene_n,
  #             "pr"=res$`Pr(>Chisq)`[2],
  #             "logLik_full" = res$LogLik[1],
  #             "logLik_joint" = res$LogLik[2],
  #             "kenp_full" = kenp_full,
  #             "kenp_joint" = kenp_joint))
  
  return(winner)
}

get_free_comb_mod <- function(rna_mod, pro_mod, current_gene, timen, num_reps, resol){
  model_n <- c("echo_lin","echo","exp_lin","exp","lin","const")
  
  t_stack <- c(rep(timen, each=num_reps), (rep(timen, each = num_reps)))
  weights <- c(rep(calc_weights_boot(unlist(genes_rna[current_gene,-1]), num_reps), each = num_reps),
               rep(calc_weights_boot(unlist(genes_pro[current_gene,-1]), num_reps), each = num_reps))
  
  t_rep <- length(rep(timen, each=num_reps))
  # ASSUME NO MISSING DATA ATM
  s1 <- rep(c(1,0), c(t_rep,t_rep))
  s2 <- rep(c(0,1), c(t_rep,t_rep))
  
  if (which(model_n == rna_mod$mod_name) <= which(model_n == pro_mod$mod_name)){
    rna_ord <- 1
    pro_ord <- 2
    
    first_mod <- rna_mod
    second_mod <- pro_mod
    
    y_stack <- c(unlist(genes_rna[current_gene,-1]), unlist(genes_pro[current_gene,-1]))
    
    mod_name <- paste0(rna_mod$mod_name,"_",pro_mod$mod_name)
  } else {
    rna_ord <- 2
    pro_ord <- 1
    
    first_mod <- pro_mod
    second_mod <- rna_mod
    
    y_stack <- c(unlist(genes_pro[current_gene,-1]), unlist(genes_rna[current_gene,-1]))
    
    mod_name <- paste0(pro_mod$mod_name,"_",rna_mod$mod_name)
  }
  temp <- data.frame(y = y_stack, t =t_stack, s1 = s1, s2 = s2, resol = rep(resol,length(y_stack)))
  temp <- temp[!is.na(temp$y),] # remove any missing data points-NONE IN EXAMPLE
  
  
  coeff1 <- as.numeric(first_mod$param[-1])
  coeff2 <- as.numeric(second_mod$param[-1])
  
  # model choice - very long switch statement!
  free.fit <- {
    switch(mod_name,
           "echo_lin_const" = {
             nls_edit(
               formula = y ~ echo_lin_const(a1,gam1,omega1,phi1,slo1,y_shift1,y_shift2,t, s1, s2),
               data=temp,
               start = list(gam1=coeff1[1], a1=coeff1[2], omega1 = coeff1[3], phi1 = coeff1[4], slo1 = coeff1[5], y_shift1=coeff1[6],
                            y_shift2 = coeff2[1]
               ),
               control = nls.control(maxiter = 0)
             )
           },
           "echo_lin_lin" = {
             nls_edit(
               formula = y ~ echo_lin_lin(a1,gam1,omega1,phi1,slo1,y_shift1,slo2,y_shift2,t, s1, s2),
               data=temp,
               start = list(gam1=coeff1[1], a1=coeff1[2], omega1 = coeff1[3], phi1 = coeff1[4], slo1 = coeff1[5], y_shift1=coeff1[6],
                            slo2 = coeff2[1], y_shift2 = coeff2[2]
               ),
               control = nls.control(maxiter = 0)
             )
           },
           "echo_lin_exp" = {
             nls_edit(
               formula = y ~ echo_lin_exp(a1,gam1,omega1,phi1,slo1,y_shift1,a2,gamma2,y_shift2,t, s1, s2),
               data=temp,
               start = list(gam1=coeff1[1], a1=coeff1[2], omega1 = coeff1[3], phi1 = coeff1[4], slo1 = coeff1[5], y_shift1=coeff1[6],
                            a2= coeff2[1],gamma2 = coeff2[2], y_shift2 = coeff2[3]
               ),
               control = nls.control(maxiter = 0)
             )
           },
           "echo_lin_exp_lin" = {
             nls_edit(
               formula = y ~ echo_lin_exp_lin(a1,gam1,omega1,phi1,slo1,y_shift1,a2,gamma2,slo2,y_shift2,t, s1, s2),
               data=temp,
               start = list(gam1=coeff1[1], a1=coeff1[2], omega1 = coeff1[3], phi1 = coeff1[4], slo1 = coeff1[5], y_shift1=coeff1[6],
                            a2= coeff2[1],gamma2 = coeff2[2],slo2 = coeff2[3], y_shift2 = coeff2[4]
               ),
               control = nls.control(maxiter = 0)
             )
           },
           "echo_lin_echo" = {
             nls_edit(
               formula = y ~ echo_lin_echo(a1,gam1,omega1,phi1,slo1,y_shift1,a2,gam2,omega2,phi2,y_shift2,t, s1, s2),
               data=temp,
               start = list(gam1=coeff1[1], a1=coeff1[2], omega1 = coeff1[3], phi1 = coeff1[4], slo1 = coeff1[5], y_shift1=coeff1[6],
                            gam2=coeff2[1], a2=coeff2[2], omega2 = coeff2[3], phi2 = coeff2[4], y_shift2=coeff2[5]
               ),
               control = nls.control(maxiter = 0)
             )
           },
         "echo_const" = {
           nls_edit(
             formula = y ~ echo_const(a1,gam1,omega1,phi1,y_shift1,y_shift2,t, s1, s2),
             data=temp,
             start = list(gam1=coeff1[1], a1=coeff1[2], omega1 = coeff1[3], phi1 = coeff1[4], y_shift1=coeff1[5],
                          y_shift2 = coeff2[1]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "echo_lin" = {
           nls_edit(
             formula = y ~ echo_lin(a1,gam1,omega1,phi1,y_shift1,slo2,y_shift2,t, s1, s2),
             data=temp,
             start = list(gam1=coeff1[1], a1=coeff1[2], omega1 = coeff1[3], phi1 = coeff1[4], y_shift1=coeff1[5],
                          slo2 = coeff2[1], y_shift2 = coeff2[2]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "echo_exp" = {
           nls_edit(
             formula = y ~ echo_exp(a1,gam1,omega1,phi1,y_shift1,a2,gamma2,y_shift2,t, s1, s2),
             data=temp,
             start = list(gam1=coeff1[1], a1=coeff1[2], omega1 = coeff1[3], phi1 = coeff1[4], y_shift1=coeff1[5],
                          a2= coeff2[1],gamma2 = coeff2[2], y_shift2 = coeff2[3]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "echo_exp_lin" = {
           nls_edit(
             formula = y ~ echo_exp_lin(a1,gam1,omega1,phi1,y_shift1,a2,gamma2,slo2,y_shift2,t, s1, s2),
             data=temp,
             start = list(gam1=coeff1[1], a1=coeff1[2], omega1 = coeff1[3], phi1 = coeff1[4], y_shift1=coeff1[5],
                          a2= coeff2[1],gamma2 = coeff2[2],slo2 = coeff2[3], y_shift2 = coeff2[4]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "echo_echo" = {
           nls_edit(
             formula = y ~ echo_echo(a1,gam1,omega1,phi1,y_shift1,a2,gam2,omega2,phi2,y_shift2,t, s1, s2),
             data=temp,
             start = list(gam1=coeff1[1], a1=coeff1[2], omega1 = coeff1[3], phi1 = coeff1[4], y_shift1=coeff1[5],
                          gam2=coeff2[1], a2=coeff2[2], omega2 = coeff2[3], phi2 = coeff2[4], y_shift2=coeff2[5]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "exp_lin_const" = {
           nls_edit(
             formula = y ~ exp_lin_const(a1,gamma1,slo1,y_shift1,y_shift2,t, s1, s2),
             data=temp,
             start = list(a1= coeff1[1],gamma1 = coeff1[2],slo1 = coeff1[3], y_shift1 = coeff2[4],
                          y_shift2 = coeff2[1]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "exp_lin_lin" = {
           nls_edit(
             formula = y ~ exp_lin_lin(a1,gamma1,slo1,y_shift1,slo2,y_shift2,t, s1, s2),
             data=temp,
             start = list(a1= coeff1[1],gamma1 = coeff1[2],slo1 = coeff1[3], y_shift1 = coeff2[4],
                          slo2 = coeff2[1], y_shift2 = coeff2[2]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "exp_lin_exp" = {
           nls_edit(
             formula = y ~ exp_lin_exp(a1,gamma1,slo1,y_shift1,a2,gamma2,y_shift2,t, s1, s2),
             data=temp,
             start = list(a1= coeff1[1],gamma1 = coeff1[2],slo1 = coeff1[3], y_shift1 = coeff2[4],
                          a2= coeff2[1],gamma2 = coeff2[2], y_shift2 = coeff2[3]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "exp_lin_exp_lin" = {
           nls_edit(
             formula = y ~ exp_lin_exp_lin(a1,gamma1,slo1,y_shift1,a2,gamma2,slo2,y_shift2,t, s1, s2),
             data=temp,
             start = list(a1= coeff1[1],gamma1 = coeff1[2],slo1 = coeff1[3], y_shift1 = coeff2[4],
                          a2= coeff2[1],gamma2 = coeff2[2],slo2 = coeff2[3], y_shift2 = coeff2[4]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "exp_const" = {
           nls_edit(
             formula = y ~ exp_const(a1,gamma1,y_shift1,y_shift2,t, s1, s2),
             data=temp,
             start = list(a1= coeff1[1],gamma1 = coeff1[2], y_shift1 = coeff1[3],
                          y_shift2 = coeff2[1]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "exp_lin" = {
           nls_edit(
             formula = y ~ exp_lin(a1,gamma1,y_shift1,slo2,y_shift2,t, s1, s2),
             data=temp,
             start = list(a1= coeff1[1],gamma1 = coeff1[2], y_shift1 = coeff1[3],
                          slo2 = coeff2[1], y_shift2 = coeff2[2]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "exp_exp" = {
           nls_edit(
             formula = y ~ exp_exp(a1,gamma1,y_shift1,a2,gamma2,y_shift2,t, s1, s2),
             data=temp,
             start = list(a1= coeff1[1],gamma1 = coeff1[2], y_shift1 = coeff1[3],
                          a2= coeff2[1],gamma2 = coeff2[2], y_shift2 = coeff2[3]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "lin_const" = {
           nls_edit(
             formula = y ~ lin_const(slo1,y_shift1,y_shift2,t, s1, s2),
             data=temp,
             start = list(slo1 = coeff1[1], y_shift1 = coeff1[2],
                          y_shift2 = coeff2[1]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "lin_lin" = {
           nls_edit(
             formula = y ~ lin_lin(slo1,y_shift1,slo2,y_shift2,t, s1, s2),
             data=temp,
             start = list(slo1 = coeff1[1], y_shift1 = coeff1[2],
                          slo2 = coeff2[1], y_shift2 = coeff2[2]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "const_const" = {
           nls_edit(
             formula = y ~ const_const(y_shift1,y_shift2,t, s1, s2),
             data=temp,
             start = list(y_shift1 = coeff1[1],
                          y_shift2 = coeff2[1]
             ),
             control = nls.control(maxiter = 0)
           )
         }
         )
  }
  
  return(free.fit)
}

reeval_best_mod <- function(joint_mod, rna_mod, pro_mod, current_gene, timen, num_reps, resol){
  free.fit <- get_free_comb_mod(rna_mod, pro_mod, current_gene, timen, num_reps, resol)
  
  # calculate AICs
  # n <- length(as.numeric(genes[current_gene,-1]))
  
  
  mod_ord <- c("free","joint")
  aic <- BIC(free.fit, joint_mod$model)
  
  if (mod_ord[which.min(aic$BIC)] == "free"){
    return(
      list(
        "mod_name" = "free",
        "aic_free" = aic$BIC[1],
        "aic_joint" = aic$BIC[2]
      )
    )
  } else {
    return(
      list(
        "mod_name" = joint_mod$mod_name,
        "joint_mod" = joint_mod,
        "aic_free" = aic$BIC[1],
        "aic_joint" = aic$BIC[2]
      )
    )
  }
}

# does not work with replicates at the moment
interpolate_pts <- function(current_gene, timen, resol, num_reps, genes){
  intr_pts <- as.numeric(genes[current_gene,-1])
  
  if (resol > 2){
    intr_pts <- approx(timen, intr_pts, xout = seq(timen[1],timen[length(timen)],2))
  }
  
  return(intr_pts)
}

#UNNEEDED
# function for determining whether gene values are deviating or constant for full matrix
# inputs:
#  current_gene: row number of current gene we want to calculate parameters for
# outputs:
#  boolean vector if there is deviation within the gene values
is_deviating_all <- function(){
  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])
  
  # vector of standard deviations
  all_row_stdev <- sqrt(rowSums((all_reps - rowMeans(all_reps,na.rm = TRUE))^2, na.rm = TRUE)/(dim(all_reps)[2] - 1))
  
  # return vector if there is deviation within the gene values
  return (all_row_stdev != 0)
}

# function for determining whether gene values are unexpressed (less than 70% expressed (i.e., below specified cutoff)) for full matrix
# inputs:
#  rem_unexpr_amt: what percentage of the time course must be expressed
#  rem_unexpr_amt_below: cutoff for expression
#  boolean if there is 70% expression
genes_unexpressed_all <- function(genes, rem_unexpr_amt, rem_unexpr_amt_below){
  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])
  # get how many genes are expressed for each gene
  tot_expressed <- rowSums(abs(all_reps) > abs(rem_unexpr_amt_below),na.rm = TRUE)
  
  # return true if amount is less than threshold
  return(tot_expressed <= (ncol(all_reps)*rem_unexpr_amt))
}


# function for determining whether gene values are deviating or constant (one replicate)
# inputs:
#  current_gene: row number of current gene we want to calculate parameters for
# outputs:
#  boolean if there is deviation within the gene values
is_deviating <- function(current_gene){
  if (all(is.na(genes[current_gene,-1])) | sum(!is.na(genes[current_gene,-1]))==1){
    stdev <- 0
  } else {
    stdev <- sd(genes[current_gene,-1], na.rm = TRUE)
  }
  return (stdev != 0)
}

# DEPRECIATED
# function for determining whether gene values are deviating or constant (multiple replicates)
# inputs:
#  current_gene: row number of current gene we want to calculate parameters for
#  num_reps: number of replicates
# outputs:
#  boolean if there is deviation within the gene values
is_deviating_rep <- function(current_gene,num_reps){
  y_val <- avg_rep(current_gene,num_reps)
  stdev <- sd(y_val,na.rm = TRUE)
  return (stdev != 0)
}

# function for determining whether gene values are unexpressed (less than 70% expressed (i.e., not 0))
# inputs:
#  current_gene: row number of current gene we want to calculate parameters for
# outputs:
#  boolean if there is 70% expression
genes_unexpressed <- function(current_gene,rem_unexpr_amt){
  y_val <- genes[current_gene,!as.logical(is.na(genes[current_gene,]))]
  y_val <- y_val[,-1]
  tot_expressed <- sum(y_val != 0,na.rm = TRUE)
  return(tot_expressed <= (length(y_val)*rem_unexpr_amt))
}

# function to adjust pvals according to the benjamini-hochberg criterion
# inputs:
#  pvals: vector of pvalues
# outputs:
#  BH adjusted pvalues
adjust_p_values <- function(pvals){
  return (p.adjust(unlist(pvals), method = "BH"))
}

# function to calculate a rolling average (3-wide window) for a set of data (smoothing) for paired
# (tied) replicates
# inputs:
#  is_weighted: logical if there is smoothing, is it weighted (1,2,1) smoothing,
#    or unweighted smoothing (1,1,1)
#  num_reps: number of replicates
# outputs:
#  smoothed data
smoothing_all_tied <- function(is_weighted, num_reps, genes, avg_genes, timen){
  # originally based on heat map code, but it will work fine here
  
  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])
  
  if (!is_weighted){
    
    #if there are replicates, average the relative expression for each replicate
    center_reps <- list() # to store actual replicate matrix
    mtx_count <- list() # to store how many are NA
    for (i in 1:num_reps){
      center_reps[[i]] <- all_reps[, seq(i,ncol(all_reps), by=num_reps)]
      mtx_count[[i]] <- is.na(center_reps[[i]])
      center_reps[[i]][is.na(center_reps[[i]])] <- 0
    }
    
    repmtx_l <- list() # store amount to divide by for each rep
    repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(all_reps))+3 # to store how many we should divide by
    repmtx[,c(1,ncol(repmtx))] <- repmtx[,c(1,ncol(repmtx))] - 1 # we only divide by 2 at edges
  } else{ # weighted averaging
    
    #if there are replicates, average the relative expression for each replicate
    center_reps <- list() # to store actual replicate matrix
    mtx_count <- list() # to store how many are NA
    for (i in 1:num_reps){
      center_reps[[i]] <- all_reps[, seq(i,ncol(all_reps), by=num_reps)]*2
      mtx_count[[i]] <- is.na(center_reps[[i]])*2
      center_reps[[i]][is.na(center_reps[[i]])] <- 0
    }
    
    repmtx_l <- list() # store amount to divide by for each rep
    repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(all_reps))+4 # to store how many we should divide by
    repmtx[,c(1,ncol(repmtx))] <- repmtx[,c(1,ncol(repmtx))] - 1 # we only divide by 3 at edges
    
  }
  for (i in 1:num_reps){
    # sum the replicates
    left <- cbind(matrix(0,nrow(all_reps),1),center_reps[[i]][,-ncol(center_reps[[i]])]/2) # left shifted matrix
    right <- cbind(center_reps[[i]][,-1]/2,matrix(0,nrow(genes),1)) # right shifted matrix
    center_reps[[i]] <- left + center_reps[[i]] + right
    
    # figure out how many replicates are actually available for each time point
    left_na <- cbind(matrix(0,nrow(all_reps),1),mtx_count[[i]][,-ncol(mtx_count[[i]])]/2) # left shifted matrix
    right_na <- cbind(mtx_count[[i]][,-1]/2,matrix(0,nrow(genes),1)) # right shifted matrix
    repmtx_l[[i]] <- repmtx - left_na - mtx_count[[i]] - right_na
    # to avoid division by 0 and induce NAs if there are no time points available
    repmtx_l[[i]][repmtx_l[[i]]==0] <- NA
  }
  
  dat <- genes
  for (x in 0:(num_reps-1)){
    dat[,seq(2+x,ncol(genes),by=num_reps)] <- center_reps[[x+1]]/repmtx_l[[x+1]]
  }
  dat[is.na(genes)] <- NA # do not impute missing values
  return(dat)
}

# function to calculate a rolling average (3-wide window) for a set of data (smoothing) for paired
# (tied) replicates
# inputs:
#  is_weighted: logical if there is smoothing, is it weighted (1,2,1) smoothing,
#    or unweighted smoothing (1,1,1)
#  num_reps: number of replicates
# outputs:
#  smoothed data
smoothing_all_untied <- function(is_weighted, num_reps, genes, avg_genes, timen){
  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])
  
  if (!is_weighted){
    
    #if there are replicates, average the relative expression for each replicate
    center_reps <- list() # to store actual replicate matrix
    mtx_count <- list() # to store how many are
    
    side_reps <- avg_genes
    mtx_side_count <- is.na(side_reps)
    side_reps[is.na(side_reps)] <- 0
    for (i in 1:num_reps){
      center_reps[[i]] <- all_reps[, seq(i,ncol(all_reps), by=num_reps)]
      mtx_count[[i]] <- is.na(center_reps[[i]])
      center_reps[[i]][is.na(center_reps[[i]])] <- 0
    }
    
    repmtx_l <- list() # store amount to divide by for each rep
    repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(all_reps))+3 # to store how many we should divide by
    repmtx[,c(1,ncol(repmtx))] <- repmtx[,c(1,ncol(repmtx))] - 1 # we only divide by 3 at edges
  } else{ # weighted averaging
    
    #if there are replicates, average the relative expression for each replicate
    center_reps <- list() # to store actual replicate matrix
    mtx_count <- list() # to store how many are
    
    side_reps <- avg_genes
    mtx_side_count <- is.na(side_reps)
    side_reps[is.na(side_reps)] <- 0
    
    for (i in 1:num_reps){
      center_reps[[i]] <- all_reps[, seq(i,ncol(all_reps), by=num_reps)]*2
      mtx_count[[i]] <- is.na(center_reps[[i]])*2
      center_reps[[i]][is.na(center_reps[[i]])] <- 0
    }
    
    repmtx_l <- list() # store amount to divide by for each rep
    repmtx <- matrix(0L,ncol = length(timen),nrow = nrow(all_reps))+4 # to store how many we should divide by
    repmtx[,c(1,ncol(repmtx))] <- repmtx[,c(1,ncol(repmtx))] - 1 # we only divide by 3 at edges
    
  }
  for (i in 1:num_reps){
    # sum the replicates
    left <- cbind(matrix(0,nrow(genes),1),side_reps[,-ncol(side_reps)]) # left shifted matrix
    right <- cbind(side_reps[,-1],matrix(0,nrow(genes),1)) # right shifted matrix
    center_reps[[i]] <- left + center_reps[[i]] + right
    
    # figure out how many replicates are actually available for each time point
    left_na <- cbind(matrix(0,nrow(all_reps),1),mtx_side_count[,-ncol(mtx_side_count)]) # left shifted matrix
    right_na <- cbind(mtx_side_count[,-1],matrix(0,nrow(genes),1)) # right shifted matrix
    repmtx_l[[i]] <- repmtx - left_na - mtx_count[[i]] - right_na
    # to avoid division by 0 and induce NAs if there are no time points available
    repmtx_l[[i]][repmtx_l[[i]]==0] <- NA
  }
  
  dat <- genes # assigning to dataframe to return
  for (x in 0:(num_reps-1)){
    dat[,seq(2+x,ncol(genes),by=num_reps)] <- center_reps[[x+1]]/repmtx_l[[x+1]]
  }
  dat[is.na(genes)] <- NA # do not impute missing values
  return(dat)
}

# function to normalize expressions in a matricized manner, by row
# thanks to: https://stackoverflow.com/questions/25099825/row-wise-variance-of-a-matrix-in-r
# outputs:
#   normalized expression data
normalize_all <- function(genes){
  
  #get matrix of just the relative expression over time
  all_reps <- as.matrix(genes[,2:ncol(genes)])
  
  # vector of row means
  all_row_mean <- rowMeans(all_reps, na.rm = TRUE)
  # vector of standard deviations
  all_row_stdev <- sqrt(rowSums((all_reps - rowMeans(all_reps,na.rm = TRUE))^2, na.rm = TRUE)/(dim(all_reps)[2] - 1))
  
  # get the matrix of normalized expressions
  all_reps_normal <- (all_reps - all_row_mean)/all_row_stdev
  # if standard deviation is 0, imposes NA, so this expression shouldn't be considered anyway
  # and is now constant
  all_reps_normal[is.na(all_reps_normal)] <- 0
  
  # create dataframe with normalized expressions
  dat <- genes
  dat[,-1] <- all_reps_normal
  dat[is.na(genes)] <- NA # do not impute missing values
  
  return(list("dat"=dat, "means"=all_row_mean, "stdevs"=all_row_stdev))
}

# DEPRECIATED
# function to calculate a rolling average (3-wide window) for a set of data (smoothing) for paired
# (tied) replicates
# smooths each replicate separately
# inputs:
#  current_gene: row number of current gene we want to calculate parameters for
#  is_weighted: logical if there is smoothing, is it weighted (1,2,1) smoothing,
#    or unweighted smoothing (1,1,1)
#  num_reps: number of replicates
# outputs:
#  smoothed data
smoothing_tied <- function(current_gene, is_weighted, num_reps){
  ldat <- lapply(c(0:(num_reps - 1)), function (x) { # return each replicate in a list of smoothed data
    if (is_weighted){ # (1,2,1) weighted average
      center.dat <- sapply(seq(num_reps+2+x,ncol(genes)-num_reps,by= num_reps), function(z)
        weighted.mean(genes[current_gene,seq(z-num_reps,z+num_reps,by = num_reps)],c(1,2,1),na.rm=TRUE))
      first <- weighted.mean(genes[current_gene,c(2+x,2+num_reps+x)],c(2,1),na.rm=TRUE)
      last <- weighted.mean(genes[current_gene,c(ncol(genes)-(2*num_reps)+1+x,ncol(genes)-num_reps+x+1)],c(1,2),na.rm=TRUE)
    }
    else{ # (1,1,1) average
      center.dat <- sapply(seq(num_reps+2+x,ncol(genes)-num_reps,by= num_reps), function(z)
        weighted.mean(genes[current_gene,seq(z-num_reps,z+num_reps,by = num_reps)],c(1,1,1),na.rm=TRUE))
      first <- weighted.mean(genes[current_gene,c(2+x,2+num_reps+x)],c(1,1),na.rm=TRUE)
      last <- weighted.mean(genes[current_gene,c(ncol(genes)-(2*num_reps)+1+x,ncol(genes)-num_reps+x+1)],c(1,1),na.rm=TRUE)
    }
    return (c(first,center.dat,last))}
  )
  dat <- genes[current_gene,] # original data
  for (x in 0:(num_reps-1)){ # replace original data by smoothed data, keeping na values
    ldat[[x+1]][is.na(genes[current_gene,seq(2+x,ncol(genes),by=num_reps)])] <- NA
    dat[,seq(2+x,ncol(genes),by=num_reps)] <- ldat[[x+1]]
  }
  
  return(dat)
}



# DEPRECIATED
# function to calculate a rolling average (3-wide window) for a set of data (smoothing) for unpaired
# (untied) replicates
# smooths each replicate with the following scheme: left average, data of specific replicate, right average
# inputs:
#  current_gene: row number of current gene we want to calculate parameters for
#  is_weighted: logical if there is smoothing, is it weighted (1,2,1) smoothing,
#    or unweighted smoothing (1,1,1)
#  num_reps: number of replicates
# outputs:
#  smoothed data
smoothing_untied <- function(current_gene, is_weighted, num_reps){
  y_val <- avg_rep(current_gene,num_reps) # calculating the average of the technical reps
  y_val <- c(genes[current_gene,1],y_val) # adding the gene name for similarity of index
  
  ldat <- lapply(c(0:(num_reps - 1)), function (x) { # return each replicate in a list of smoothed data
    if (is_weighted){ # (1,2,1) weighted average
      center.dat <- sapply(seq(num_reps+2+x,ncol(genes)-num_reps,by= num_reps), function(z)
        weighted.mean(c(as.numeric(y_val[z-num_reps]),genes[current_gene,z],as.numeric(y_val[z+num_reps])),c(1,2,1),na.rm=TRUE))
      first <- weighted.mean(c(genes[current_gene,2+x],as.numeric(y_val[2+num_reps+x])),c(2,1),na.rm=TRUE)
      last <- weighted.mean(c(as.numeric(y_val[ncol(genes)-(2*num_reps)+1+x]),genes[current_gene,ncol(genes)-num_reps+x+1]),c(1,2),na.rm=TRUE)
    }
    else{ # (1,1,1) average
      center.dat <- sapply(seq(num_reps+2+x,ncol(genes)-num_reps,by= num_reps), function(z)
        weighted.mean(c(as.numeric(y_val[z-num_reps]),genes[current_gene,z],as.numeric(y_val[z+num_reps])),c(1,1,1),na.rm=TRUE))
      first <- weighted.mean(c(genes[current_gene,2+x],as.numeric(y_val[2+num_reps+x])),c(1,1),na.rm=TRUE)
      last <- weighted.mean(c(as.numeric(y_val[ncol(genes)-(2*num_reps)+1+x]),genes[current_gene,ncol(genes)-num_reps+x+1]),c(1,1),na.rm=TRUE)
    }
    return (c(first,center.dat,last))}
  )
  dat <- genes[current_gene,] # original data
  for (x in 0:(num_reps-1)){ # replace original data by smoothed data, keeping na values
    ldat[[x+1]][is.na(genes[current_gene,seq(2+x,ncol(genes),by=num_reps)])] <- NA
    dat[,seq(2+x,ncol(genes),by=num_reps)] <- ldat[[x+1]]
  }
  
  return(dat)
}

# DEPRECIATED
de_linear_trend <- function(current_gene, time_begin, time_end, resol,num_reps,timen,
                            rem_unexpr_vect){
  if (!rem_unexpr_vect[current_gene]){
    gene_n <- as.character(genes[current_gene,1]) # gene name
    rep_timen <- rep(timen,each=num_reps)
    y_val <- as.numeric(as.character(t(genes[current_gene,-1]))) # all the y values
    #do linear regression
    trend_test <- lm((y_val) ~ rep_timen)
    coeff <- trend_test$coefficients # resulting coefficients
    
    # detrend the data
    adjusted_y_val <- y_val - (coeff[1] + rep_timen*coeff[2])
    df2 <- cbind(data.frame("Gene.Name" = gene_n),rbind(adjusted_y_val))
    colnames(df2) <- colnames(genes)
    return (df2)
  } else {
    return (genes[current_gene,])
  }
  
}

# function to remove linear trend from all data
de_linear_trend_all <- function(timen,num_reps,tied){
  all_rep <- as.matrix(genes[,-1]) # y values for linear fit
  if (!tied){ # if they're not paired, we just fit an aggregate data model
    
    # x values for linear fit
    xrow <- rep(timen,each=num_reps)
    xmtx <- matrix(rep(xrow,each=nrow(all_rep)),nrow = nrow(all_rep))
    
    # covariance
    cov <- rowSums((all_rep-rowMeans(all_rep,na.rm = TRUE))*(xmtx-rowMeans(xmtx)),na.rm = TRUE)
    # variance
    var <- rowSums((xmtx - rowMeans(xmtx))^2,na.rm = TRUE)
    
    beta <- cov/var
    alph <- rowMeans(all_rep,na.rm = TRUE)-(beta*rowMeans(xmtx))
    
    df <- all_rep-(alph+(beta*xmtx)) # linear fit
  } else { # we have to do the models separately for each replicate
    # preallocate matrix where we put results
    df <- matrix(NA, nrow = dim(all_rep)[1], ncol = dim(all_rep)[2])
    
    # x values for linear fit
    xmtx <- matrix(rep(timen,each=nrow(all_rep)),nrow = nrow(all_rep))
    # preallocating to store slopes
    beta.df <- data.frame(matrix(NA, nrow(genes), num_reps))
    for (i in 1:num_reps){
      each_rep <- all_rep[,seq(i,ncol(all_rep),by=num_reps)]
      
      # covariance
      cov <- rowSums((each_rep-rowMeans(each_rep,na.rm = TRUE))*(xmtx-rowMeans(xmtx)),na.rm = TRUE)
      # variance
      var <- rowSums((xmtx - rowMeans(xmtx))^2,na.rm = TRUE)
      
      beta.df[,i] <- beta <- cov/var
      alph <- rowMeans(each_rep,na.rm = TRUE)-(beta*rowMeans(xmtx))
      
      df[,seq(i,ncol(all_rep),by=num_reps)] <- each_rep -(alph+(beta*xmtx)) # linear fit
    }
    # get average slope for each gene
    beta <- rowMeans(beta.df, na.rm = T)
  }
  # get the data frame correctly set up for returning
  res_df <- genes
  res_df[,-1] <- df
  
  # now we return the slope and the altered expressions
  res_list <- list("res_df" = res_df, "beta" = beta, "alph" = alph)
  
  return (res_list)
}

# THE FOLLOWING ARE VERY SPECIFICALLY EDITED FOR OUR PROBLEM

# for getting an object of class nls
nlsModel_edit <- function (form, data, start, wts) 
{
  thisEnv <- environment()
  env <- new.env(parent = environment(form))
  for (i in names(data)) {
    assign(i, data[[i]], envir = env)
  }
  ind <- as.list(start)
  parLength <- 0
  for (i in names(ind)) {
    temp <- start[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp, envir = env)
    ind[[i]] <- parLength + seq(along = start[[i]])
    parLength <- parLength + length(start[[i]])
  }
  useParams <- rep(TRUE, parLength)
  lhs <- eval(form[[2]], envir = env)
  rhs <- eval(form[[3]], envir = env)
  resid <- lhs - rhs
  dev <- sum(resid^2)
  if (is.null(attr(rhs, "gradient"))) {
    # the gradient doesn't matter, we're not evaluating anything
    # getRHS.noVarying <- function() numericDeriv(form[[3]], 
                                                # names(ind), env)
    # getRHS <- getRHS.noVarying
    # rhs <- getRHS()
    attr(rhs, "gradient") <- matrix(1, nrow = length(data$y), ncol = length(names(ind)))
  }
  else {
    getRHS.noVarying <- function() eval(form[[3]], envir = env)
    getRHS <- getRHS.noVarying
  }
  dimGrad <- dim(attr(rhs, "gradient"))
  marg <- length(dimGrad)
  if (marg > 0) {
    gradSetArgs <- vector("list", marg + 1)
    for (i in 2:marg) gradSetArgs[[i]] <- rep(TRUE, dimGrad[i - 
                                                              1])
    useParams <- rep(TRUE, dimGrad[marg])
  }
  else {
    gradSetArgs <- vector("list", 2)
    useParams <- rep(TRUE, length(attr(rhs, "gradient")))
  }
  npar <- length(useParams)
  gradSetArgs[[1]] <- (~attr(ans, "gradient"))[[2]]
  gradCall <- switch(length(gradSetArgs) - 1, call("[", gradSetArgs[[1]], 
                                                   gradSetArgs[[2]]), call("[", gradSetArgs[[1]], gradSetArgs[[2]], 
                                                                           gradSetArgs[[2]]), call("[", gradSetArgs[[1]], gradSetArgs[[2]], 
                                                                                                   gradSetArgs[[2]], gradSetArgs[[3]]), call("[", gradSetArgs[[1]], 
                                                                                                                                             gradSetArgs[[2]], gradSetArgs[[2]], gradSetArgs[[3]], 
                                                                                                                                             gradSetArgs[[4]]))
  getRHS.varying <- function() {
    ans <- getRHS.noVarying()
    attr(ans, "gradient") <- eval(gradCall)
    ans
  }
  QR <- qr(attr(rhs, "gradient"))
  qrDim <- min(dim(QR$qr))
  # if (QR$rank < qrDim) 
  #   stop("singular gradient matrix at initial parameter estimates")
  getPars.noVarying <- function() unlist(setNames(lapply(names(ind), 
                                                         get, envir = env), names(ind)))
  getPars.varying <- function() unlist(setNames(lapply(names(ind), 
                                                       get, envir = env), names(ind)))[useParams]
  getPars <- getPars.noVarying
  internalPars <- getPars()
  setPars.noVarying <- function(newPars) {
    assign("internalPars", newPars, envir = thisEnv)
    for (i in names(ind)) {
      assign(i, unname(newPars[ind[[i]]]), envir = env)
    }
  }
  setPars.varying <- function(newPars) {
    internalPars[useParams] <- newPars
    for (i in names(ind)) {
      assign(i, unname(internalPars[ind[[i]]]), envir = env)
    }
  }
  setPars <- setPars.noVarying
  on.exit(remove(i, data, parLength, start, temp, m))
  m <- list(resid = function() resid, fitted = function() rhs, 
            formula = function() form, deviance = function() dev, 
            gradient = function() attr(rhs, "gradient"), conv = function() {
              rr <- qr.qty(QR, resid)
              sqrt(sum(rr[1:npar]^2)/sum(rr[-(1:npar)]^2))
            }, incr = function() qr.coef(QR, resid), setVarying = function(vary = rep(TRUE, 
                                                                                      length(useParams))) {
              assign("useParams", if (is.character(vary)) {
                temp <- logical(length(useParams))
                temp[unlist(ind[vary])] <- TRUE
                temp
              } else if (is.logical(vary) && length(vary) != length(useParams)) stop("setVarying : vary length must match length of parameters") else {
                vary
              }, envir = thisEnv)
              gradCall[[length(gradCall)]] <<- useParams
              if (all(useParams)) {
                assign("setPars", setPars.noVarying, envir = thisEnv)
                assign("getPars", getPars.noVarying, envir = thisEnv)
                assign("getRHS", getRHS.noVarying, envir = thisEnv)
                assign("npar", length(useParams), envir = thisEnv)
              } else {
                assign("setPars", setPars.varying, envir = thisEnv)
                assign("getPars", getPars.varying, envir = thisEnv)
                assign("getRHS", getRHS.varying, envir = thisEnv)
                assign("npar", length((1:length(useParams))[useParams]), 
                       envir = thisEnv)
              }
            }, setPars = function(newPars) {
              setPars(newPars)
              assign("resid", lhs - assign("rhs", getRHS(), envir = thisEnv), 
                     envir = thisEnv)
              assign("dev", sum(resid^2), envir = thisEnv)
              assign("QR", qr(attr(rhs, "gradient")), envir = thisEnv)
              return(QR$rank < min(dim(QR$qr)))
            }, getPars = function() getPars(), getAllPars = function() getPars(), 
            getEnv = function() env, trace = function() cat(format(dev), 
                                                            ": ", format(getPars()), "\n"), Rmat = function() qr.R(QR), 
            predict = function(newdata = list(), qr = FALSE) {
              eval(form[[3]], as.list(newdata), env)
            })
  class(m) <- "nlsModel"
  m
}

# for getting an object of class nls
nls_edit <- function (formula, data = parent.frame(), start, control = nls.control(), 
                      algorithm = c("default", "plinear", "port"), trace = FALSE, 
                      subset, weights, na.action, model = FALSE, lower = -Inf, 
                      upper = Inf, ...) 
{
  formula <- as.formula(formula)
  algorithm <- match.arg(algorithm)
  if (!is.list(data) && !is.environment(data)) 
    stop("'data' must be a list or an environment")
  mf <- cl <- match.call()
  varNames <- all.vars(formula)
  if (length(formula) == 2L) {
    formula[[3L]] <- formula[[2L]]
    formula[[2L]] <- 0
  }
  form2 <- formula
  form2[[2L]] <- 0
  varNamesRHS <- all.vars(form2)
  mWeights <- missing(weights)
  pnames <- if (missing(start)) {
    if (!is.null(attr(data, "parameters"))) {
      names(attr(data, "parameters"))
    }
    else {
      cll <- formula[[length(formula)]]
      if (is.symbol(cll)) {
        cll <- substitute(S + 0, list(S = cll))
      }
      fn <- as.character(cll[[1L]])
      if (is.null(func <- tryCatch(get(fn), error = function(e) NULL))) 
        func <- get(fn, envir = parent.frame())
      if (!is.null(pn <- attr(func, "pnames"))) 
        as.character(as.list(match.call(func, call = cll))[-1L][pn])
    }
  }
  else names(start)
  env <- environment(formula)
  if (is.null(env)) 
    env <- parent.frame()
  if (length(pnames)) 
    varNames <- varNames[is.na(match(varNames, pnames))]
  lenVar <- function(var) tryCatch(length(eval(as.name(var), 
                                               data, env)), error = function(e) -1L)
  if (length(varNames)) {
    n <- vapply(varNames, lenVar, 0)
    if (any(not.there <- n == -1L)) {
      nnn <- names(n[not.there])
      if (missing(start)) {
        if (algorithm == "plinear") 
          stop("no starting values specified")
        warning("No starting values specified for some parameters.\n", 
                "Initializing ", paste(sQuote(nnn), collapse = ", "), 
                " to '1.'.\n", "Consider specifying 'start' or using a selfStart model", 
                domain = NA)
        start <- setNames(as.list(rep_len(1, length(nnn))), 
                          nnn)
        varNames <- varNames[i <- is.na(match(varNames, 
                                              nnn))]
        n <- n[i]
      }
      else stop(gettextf("parameters without starting value in 'data': %s", 
                         paste(nnn, collapse = ", ")), domain = NA)
    }
  }
  else {
    if (length(pnames) && any((np <- sapply(pnames, lenVar)) == 
                              -1)) {
      message(sprintf(ngettext(sum(np == -1), "fitting parameter %s without any variables", 
                               "fitting parameters %s without any variables"), 
                      paste(sQuote(pnames[np == -1]), collapse = ", ")), 
              domain = NA)
      n <- integer()
    }
    else stop("no parameters to fit")
  }
  respLength <- length(eval(formula[[2L]], data, env))
  if (length(n) > 0L) {
    varIndex <- n%%respLength == 0
    if (is.list(data) && diff(range(n[names(n) %in% names(data)])) > 
        0) {
      mf <- data
      if (!missing(subset)) 
        warning("argument 'subset' will be ignored")
      if (!missing(na.action)) 
        warning("argument 'na.action' will be ignored")
      if (missing(start)) 
        start <- getInitial(formula, mf)
      startEnv <- new.env(hash = FALSE, parent = environment(formula))
      for (i in names(start)) assign(i, start[[i]], envir = startEnv)
      rhs <- eval(formula[[3L]], data, startEnv)
      n <- NROW(rhs)
      wts <- if (mWeights) 
        rep_len(1, n)
      else eval(substitute(weights), data, environment(formula))
    }
    else {
      vNms <- varNames[varIndex]
      if (any(nEQ <- vNms != make.names(vNms))) 
        vNms[nEQ] <- paste0("`", vNms[nEQ], "`")
      mf$formula <- as.formula(paste("~", paste(vNms, 
                                                collapse = "+")), env = environment(formula))
      mf$start <- mf$control <- mf$algorithm <- mf$trace <- mf$model <- NULL
      mf$lower <- mf$upper <- NULL
      mf[[1L]] <- quote(stats::model.frame)
      mf <- eval.parent(mf)
      n <- nrow(mf)
      mf <- as.list(mf)
      wts <- if (!mWeights) 
        model.weights(mf)
      else rep_len(1, n)
    }
    if (any(wts < 0 | is.na(wts))) 
      stop("missing or negative weights not allowed")
  }
  else {
    varIndex <- logical()
    mf <- list(0)
    wts <- numeric()
  }
  if (missing(start)) 
    start <- getInitial(formula, mf)
  for (var in varNames[!varIndex]) mf[[var]] <- eval(as.name(var), 
                                                     data, env)
  varNamesRHS <- varNamesRHS[varNamesRHS %in% varNames[varIndex]]
  m <- switch(algorithm, plinear = nlsModel_edit.plinear(formula, 
                                                         mf, start, wts), port = nlsModel_edit(formula, mf, start, 
                                                                                               wts, upper), nlsModel_edit(formula, mf, start, wts))
  ctrl <- nls.control()
  if (!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  if (algorithm != "port") {
    if (!identical(lower, -Inf) || !identical(upper, +Inf)) {
      warning("upper and lower bounds ignored unless algorithm = \"port\"")
      cl$lower <- NULL
      cl$upper <- NULL
    }
    if (ctrl$maxiter > 0){
      convInfo <- .Call(C_nls_iter, m, ctrl, trace)
    } else {
      convInfo <- list("isConv" = F, "finIter" = 0, "finTol" = NA, "stopCode" = 9,
                       "stopMessage"="The number of iterations has reached maxiter.")
    }
    nls.out <- list(m = m, convInfo = convInfo, data = substitute(data), 
                    call = cl)
  }
  else {
    pfit <- nls_port_fit(m, start, lower, upper, control, 
                         trace, give.v = TRUE)
    iv <- pfit[["iv"]]
    msg.nls <- port_msg(iv[1L])
    conv <- (iv[1L] %in% 3:6)
    if (!conv) {
      msg <- paste("Convergence failure:", msg.nls)
      if (ctrl$warnOnly) 
        warning(msg)
      else stop(msg)
    }
    v. <- port_get_named_v(pfit[["v"]])
    cInfo <- list(isConv = conv, finIter = iv[31L], finTol = v.[["NREDUC"]], 
                  nEval = c(`function` = iv[6L], gradient = iv[30L]), 
                  stopCode = iv[1L], stopMessage = msg.nls)
    cl$lower <- lower
    cl$upper <- upper
    nls.out <- list(m = m, data = substitute(data), call = cl, 
                    convInfo = cInfo, convergence = as.integer(!conv), 
                    message = msg.nls)
  }
  nls.out$call$algorithm <- algorithm
  nls.out$call$control <- ctrl
  nls.out$call$trace <- trace
  nls.out$na.action <- attr(mf, "na.action")
  nls.out$dataClasses <- attr(attr(mf, "terms"), "dataClasses")[varNamesRHS]
  if (model) 
    nls.out$model <- mf
  if (!mWeights) 
    nls.out$weights <- wts
  nls.out$control <- control
  class(nls.out) <- "nls"
  nls.out
}
