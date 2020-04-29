# MOSAIC Function Source
# By Hannah De los Santos
# Code description: Contains all the funcitons for extended harmonic oscillator work, in order to have less confusion between scripts.

# combined models ----

# function to represent the free model of ECHO linear and linear models
# inputs:
#  a1: Amplitude, ECHO linear
#  gam1: Amplitude.Change.Coefficient (amount of damping/forcing), ECHO linear
#  omega1: Radial frequency, ECHO linear
#  phi1: Phase Shift (radians), ECHO linear
#  slo1: Slope, ECHO linear
#  y_shift1: Equilibrium shift, ECHO linear
#  slo2: Slope, linear
#  y_shift2: Equilibrium shift, linear
#  t: time
#  s1: binary vector, with 1's corresponding to rna data
#  s2: binary vector, with 1's corresponding to protein data
# outputs:
#  numerical vector of model outputs
echo_lin_lin <- function(a1,gam1,omega1,phi1,slo1,y_shift1,slo2,y_shift2,t, s1, s2){
  return(((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+slo1*t+y_shift1)*s1)+((slo2*t+y_shift2)*s2))
}

# function to represent the free model of ECHO linear and exponential models
# inputs:
#  a1: Amplitude, ECHO linear
#  gam1: Amplitude.Change.Coefficient (amount of damping/forcing), ECHO linear
#  omega1: Radial frequency, ECHO linear
#  phi1: Phase Shift (radians), ECHO linear
#  slo1: Slope, ECHO linear
#  y_shift1: Equilibrium shift, ECHO linear
#  a2: Amplitude, exponential
#  gamma2: Growth Rate, exponential
#  y_shift2: Equilibrium shift, exponential
#  t: time
#  s1: binary vector, with 1's corresponding to rna data
#  s2: binary vector, with 1's corresponding to protein data
# outputs:
#  numerical vector of model outputs
echo_lin_exp <- function(a1,gam1,omega1,phi1,slo1,y_shift1,a2,gamma2,y_shift2,t, s1, s2){
  return(((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+slo1*t+y_shift1)*s1)+((a2*exp(gamma2*t)+y_shift2)*s2))
}

# function to represent the free model of ECHO and linear models
# inputs:
#  a1: Amplitude, ECHO
#  gam1: Amplitude.Change.Coefficient (amount of damping/forcing), ECHO
#  omega1: Radial frequency, ECHO
#  phi1: Phase Shift (radians), ECHO
#  y_shift1: Equilibrium shift, ECHO
#  slo2: Slope, linear
#  y_shift2: Equilibrium shift, linear
#  t: time
#  s1: binary vector, with 1's corresponding to rna data
#  s2: binary vector, with 1's corresponding to protein data
# outputs:
#  numerical vector of model outputs
echo_lin <- function(a1,gam1,omega1,phi1,y_shift1,slo2,y_shift2,t, s1, s2){
  return(((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+y_shift1)*s1)+((slo2*t+y_shift2)*s2))
}

# function to represent the free model of ECHO and exponential models
# inputs:
#  a1: Amplitude, ECHO
#  gam1: Amplitude.Change.Coefficient (amount of damping/forcing), ECHO
#  omega1: Radial frequency, ECHO
#  phi1: Phase Shift (radians), ECHO
#  y_shift1: Equilibrium shift, ECHO
#  a2: Amplitude, exponential
#  gamma2: Growth Rate, exponential
#  y_shift2: Equilibrium shift, exponential
#  t: time
#  s1: binary vector, with 1's corresponding to rna data
#  s2: binary vector, with 1's corresponding to protein data
# outputs:
#  numerical vector of model outputs
echo_exp <- function(a1,gam1,omega1,phi1,y_shift1,a2,gamma2,y_shift2,t, s1, s2){
  return(((a1*exp(-1*gam1*(t)/2)*cos((omega1*t)+phi1)+y_shift1)*s1)+((a2*exp(gamma2*t)+y_shift2)*s2))
}

# models ----

# function to represent ECHO model
# inputs:
#  a: Amplitude
#  gam: Amplitude.Change.Coefficient (amount of damping/forcing)
#  omega: Radial frequency
#  phi: Phase Shift (radians)
#  y_shift: Equilibrium shift
#  t: time
# outputs:
#  numerical vector of model outputs
alt_form <- function(a,gam,omega,phi,y_shift,t){
  return (a*exp(-1*gam*(t)/2)*cos((omega*t)+phi)+y_shift)
}

# function to represent ECHO Linear model
# inputs:
#  a: Amplitude
#  gam: Amplitude.Change.Coefficient (amount of damping/forcing)
#  omega: Radial frequency
#  phi: Phase Shift (radians)
#  slo: Slope
#  y_shift: Equilibrium shift
#  t: time
# outputs:
#  numerical vector of model outputs
echo_lin_mod <- function(a,gam,omega,phi,slo,y_shift,t){
  return (a*exp(-1*gam*(t)/2)*cos((omega*t)+phi)+slo*t+y_shift)
}

# function to represent linear model
# inputs:
#  slo: Slope
#  y_shift: Equilibrium shift
#  t: time
# outputs:
#  numerical vector of model outputs
lin_mod <- function(slo,y_shift,t){
  return(slo*t+y_shift)
}

# function to represent exponential model
# inputs:
#  a: Amplitude
#  gamma: Growth Rate
#  y_shift: Equilibrium shift
#  t: time
# outputs:
#  numerical vector of model outputs
exp_mod <- function(a,gamma,y_shift,t){
  return(a*exp(gamma*t)+y_shift)
}

# function to represent the ECHO linear completely joint model (no slack for protein)
# inputs:
#  a1: Amplitude, rna
#  a2: Amplitude, protein
#  gam1: Amplitude.Change.Coefficient (amount of damping/forcing), rna
#  gam2: Amplitude.Change.Coefficient (amount of damping/forcing), protein
#  omega1: Radial frequency, both
#  phi1: Phase Shift (radians), rna
#  phi2: Phase Shift (radians), protein
#  slo1: Slope, rna
#  slo2: Slope, protein
#  y_shift1: Equilibrium shift, rna
#  y_shift2: Equilibrium shift, protein
#  t: time
#  s1: binary vector, with 1's corresponding to rna data
#  s2: binary vector, with 1's corresponding to protein data
# outputs:
#  numerical vector of model outputs
echo_lin_joint <- function(a1,a2, gam1,gam2, omega1, phi1,phi2, slo1,slo2, y_shift1,y_shift2, t, s1, s2){
  return ((a1*s1+a2*s2)*exp(-1*(gam1*s1+gam2*s1)*(t)/2)*cos((omega1*t)+(phi1*s1+phi2*s2))+((slo1*s1 + slo2*s2)*t)+(y_shift1*s1+y_shift2*s2))
}

# function to represent the ECHO completely joint model (no slack for protein)
# inputs:
#  a1: Amplitude, rna
#  a2: Amplitude, protein
#  gam1: Amplitude.Change.Coefficient (amount of damping/forcing), rna
#  gam2: Amplitude.Change.Coefficient (amount of damping/forcing), protein
#  omega1: Radial frequency, both
#  phi1: Phase Shift (radians), rna
#  phi2: Phase Shift (radians), protein
#  y_shift1: Equilibrium shift, rna
#  y_shift2: Equilibrium shift, protein
#  t: time
#  s1: binary vector, with 1's corresponding to rna data
#  s2: binary vector, with 1's corresponding to protein data
# outputs:
#  numerical vector of model outputs
echo_joint <- function(a1,a2, gam1,gam2, omega1, phi1,phi2, y_shift1,y_shift2, t, s1, s2){
  return ((a1*s1+a2*s2)*exp(-1*(gam1*s1+gam2*s2)*(t)/2)*cos((omega1*t)+(phi1*s1+phi2*s2))+(y_shift1*s1+y_shift2*s2))
}

# function to represent the ECHO joint model (with slack for protein)
# inputs:
#  a1: Amplitude, rna
#  a2: Amplitude, protein
#  gam1: Amplitude.Change.Coefficient (amount of damping/forcing), rna
#  gam2: Amplitude.Change.Coefficient (amount of damping/forcing), protein
#  omega1: Period (hours)
#  change: Change in period for protein (hours)
#  phi1: Phase Shift (radians), rna
#  phi2: Phase Shift (radians), protein
#  y_shift1: Equilibrium shift, rna
#  y_shift2: Equilibrium shift, protein
#  t: time
#  s1: binary vector, with 1's corresponding to rna data
#  s2: binary vector, with 1's corresponding to protein data
#  resol: resolution
# outputs:
#  numerical vector of model outputs
echo_joint_res_constr <- function(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, y_shift1,y_shift2, t, s1, s2, resol){
  return ((a1*s1+a2*s2)*exp(-1*(gam1*s1+gam2*s2)*(t)/2)*cos(((2*pi/(omega1*s1+((omega1+chang)*s2)))*t)+(phi1*s1+phi2*s2))+(y_shift1*s1+y_shift2*s2))
}

# function to represent the ECHO linear joint model (with slack for protein)
# inputs:
#  a1: Amplitude, rna
#  a2: Amplitude, protein
#  gam1: Amplitude.Change.Coefficient (amount of damping/forcing), rna
#  gam2: Amplitude.Change.Coefficient (amount of damping/forcing), protein
#  omega1: Period (hours)
#  change: Change in period for protein (hours)
#  phi1: Phase Shift (radians), rna
#  phi2: Phase Shift (radians), protein
#  slo1: Slope, rna
#  slo2: Slope, protein
#  y_shift1: Equilibrium shift, rna
#  y_shift2: Equilibrium shift, protein
#  t: time
#  s1: binary vector, with 1's corresponding to rna data
#  s2: binary vector, with 1's corresponding to protein data
#  resol: resolution
# outputs:
#  numerical vector of model outputs
echo_lin_joint_res_constr <- function(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, slo1,slo2, y_shift1,y_shift2, t, s1, s2, resol){
  return ((a1*s1+a2*s2)*exp(-1*(gam1*s1+gam2*s2)*(t)/2)*cos(((2*pi/(omega1*s1+((omega1+chang)*s2)))*t)+(phi1*s1+phi2*s2))+((slo1*s1 + slo2*s2)*t)+(y_shift1*s1+y_shift2*s2))
}

# optimization support functions ----

# function to calculate optimization objective for weighted nls (with replicates)
# inputs:
#  p: paramters
#  t: time points
#  s1: binary vector, with 1's corresponding to rna data
#  s2: binary vector, with 1's corresponding to protein data
#  w: weights
#  y: experimental data
#  fcall: function call (model)
#  jcall: jacobian call
# output:
#  numeric vector of objective function
fcn <- function(p, t, s1, s2, w, y, fcall, jcall){
  sqrt(w)*(y - do.call("fcall", c(list(t = t, s1 =s1, s2 = s2), as.list(p))))
}

# function to calculate optimization objective for non-weighted nls (with no replicates)
# inputs:
#  p: paramters
#  t: time points
#  s1: binary vector, with 1's corresponding to rna data
#  s2: binary vector, with 1's corresponding to protein data
#  y: experimental data
#  fcall: function call (model)
#  jcall: jacobian call
# output:
#  numeric vector of objective function
fcn_one_rep <- function(p, t, s1, s2, y, fcall, jcall){
  (y - do.call("fcall", c(list(t = t, s1 =s1, s2 = s2), as.list(p))))
}

# function to calculate optimization objective for weighted nls (with replicates) for one omics type
# inputs:
#  p: paramters
#  t: time points
#  w: weights
#  y: experimental data
#  fcall: function call (model)
#  jcall: jacobian call
# output:
#  numeric vector of objective function
fcn_single <- function(p, t, w, y, fcall, jcall){
  sqrt(w)*(y - do.call("fcall", c(list(t = t), as.list(p))))
}

# function to calculate optimization objective for non-weighted nls (with no replicates) for one omics type
# inputs:
#  p: paramters
#  t: time points
#  y: experimental data
#  fcall: function call (model)
#  jcall: jacobian call
# output:
#  numeric vector of objective function
fcn_single_one_rep <- function(p, t, y, fcall, jcall){
  (y - do.call("fcall", c(list(t = t), as.list(p))))
}

# supports ----

# function to calculate the average of all replicates. Used primarily for cases with multiple replicates.
# inputs:
#  num_reps: number of replicates
#  genes: dataframe of gene expressions (first column is names, other columns are numeric)
#  timen: numeric vector of time points
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
#  genes: dataframe of gene expressions (first column is names, other columns are numeric)
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
#  genes: dataframe of gene expressions (first column is names, other columns are numeric)
# outputs:
#  the variance of replicates at a certain time point
calc_weights <- function(current_gene, num_reps, genes){
  return(sapply(seq(2,ncol(genes), by = num_reps), function(x) calc_var(x,current_gene,num_reps, genes)))
}

# function to calculate the variance of replicates at a certain time point for one gene
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

# function to calculate the weights for replicate fitting for one gene (these weights are the inverse of the variance at each time point)
# inputs:
#  boot_gene: gene time course being examined
#  num_reps: number of replicates
# outputs:
#  the variance of replicates at a certain time point
calc_weights_boot <- function(boot_gene, num_reps){
  return(sapply(seq(1,length(boot_gene), by = num_reps), function(x) calc_var_boot(x,boot_gene,num_reps)))
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
# calculate_param <- function(current_gene,timen,resol,num_reps,tied,is_smooth=FALSE,is_weighted=FALSE,low,high,rem_unexpr=FALSE,rem_unexpr_amt=70,run_conf = F, which_conf = "Bootstrap", harm_cut = .03, over_cut = .15, seed = 30)

# function to smooth, with weighting, a given averaged expression
# inputs:
#  y_val: vector of averaged expressions
# outputs:
#  smoothed y_val
smooth_y_val <- function(y_val){
  # smoothing the starting averages - weighted or unweighted?
  #get matrix of just the relative expression over time
  all_reps <- matrix(c(y_val), nrow = 1, ncol = length(y_val))

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
  # averaging to get smoothed replicates
  dat[,seq(1+x,ncol(all_reps),by=1)] <- center_reps[[x+1]]/repmtx_l[[x+1]]
  dat[is.na(all_reps)] <- NA # do not impute missing values

  return(dat[1,])
}

# function to calculate the "peaks" vector for starting points for echo-type models, "peaks" can be peaks or troughs
# inputs:
#  y_val: vector of averaged expressions
#  resol: resolution
#  y0: initial equilibrium shift
#  timen: numeric vector of time points
# outputs:
#  list with numeric vector of peaks, and the time points for those peaks
calc_peaks <- function(y_val, resol, y0, timen){
  # figure out the resolution modifier (must have at least 24 hour data)
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
  # go through all time points and calculate peaks (max values)
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

  # calculate amount of troughs
  troughs <- c() # vector of trough values
  troughs_time <- c() # vector of trough times
  counting <- 1 # counter
  # go through all time points and calculate peaks (max values)
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

  peaks_per_diff <- 0
  troughs_per_diff <- 0
  if (length(peaks)>1 & length(troughs)>1){
    # new criterion: which peaks are more evenly distributed?
    # what should the period be, based on the amount of peaks/troughs seen in the time course
    peaks_per <- length(timen)*resol/length(peaks)
    troughs_per <- length(timen)*resol/length(troughs)

    # get avg period, calculated as the difference in time between successive peaks
    peaks_avg_per <- mean(sapply(seq(length(peaks_time),2,-1),
                                 function(x){peaks_time[x]-peaks_time[x-1]}))
    troughs_avg_per <- mean(sapply(seq(length(troughs_time),2,-1),
                                   function(x){troughs_time[x]-troughs_time[x-1]}))

    # calculate the difference between the theoretical and estimated period
    peaks_per_diff <- abs(peaks_per-peaks_avg_per)
    troughs_per_diff <- abs(troughs_per-troughs_avg_per)
  }

  # if there is only one peak, or if our troughs are more consistent, we choose troughs
  if ((length(peaks)<=1 & length(troughs)>length(peaks)) |
      (troughs_per_diff < peaks_per_diff)){
    # flip the troughs for accurate calculations later
    # absolute value the difference from the midline
    peaks <- abs(y0 - troughs)
    peaks_time <- troughs_time
  } else {# else we choose peaks
    # absolute value the difference from the midline
    peaks <- abs(peaks - y0)
  }
  # it's possible for this not to work if you don't have more than one
  # oscillation over your time course
  # }

  return(list("peaks"=peaks, "peaks_time"=peaks_time))

}

# function to calculate the starting points for echo models
# inputs:
#  y_val: vector of averaged expressions
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
# outputs:
#  list with all starting points for parameters (amplitude, ac coeff, frequency, phase shift, equilibrium shift)
calc_start_echo <- function(y_val, over_cut){
  timen <- timen
  # smooth the expression
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
  } else{ # otherwise forcing is likely
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

  # restricting ac coeff to not be overexpressed/repressed
  if (gam0 < -over_cut){
    gam0 <- -over_cut
  } else if (gam0 > over_cut){
    gam0 <- over_cut
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


  # we estimate our phase shift on the second and third nonmissing sets of points for accuracy
  # if you have less than 3 points nonmissing, I have no hope for you

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

  min_vect <- rep(Inf, length(0:11))
  for (i in 0:11){
    # if the phase shift at the beginning, middle, and end are the smallest value available for the fitted value
    min_vect[i+1] <- sum(abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[beg1])-y_val[beg1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[mid1])-y_val[mid1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[en1])-y_val[en1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[beg2])-y_val[beg2]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[mid2])-y_val[mid2]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[en2])-y_val[en2]))
  }
  phi0 <- (which.min(abs(min_vect))-1)*pi/6 # intial value for phase shift

  # form all starting parameters into a nice list
  start_param <- list(gam=gam0,a=a0,omega=w0,phi=phi0,y_shift=y0)

  return(start_param)
}

# function to calculate the starting points for echo linear models
# inputs:
#  y_val: vector of averaged expressions
#  slo0: initial slope
#  y_shift0: initial equilibrium value
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
# outputs:
#  list with all starting points for parameters (amplitude, ac coeff, frequency, phase shift, slope, equilibrium shift)
calc_start_echo_lin <- function(y_val, slo0, y_shift0, over_cut){
  # subtract sloped baseline
  y_val <- y_val - lin_mod(slo0,y_shift0,timen)

  # smooth y_val
  y_val <- smooth_y_val(y_val)

  # calculate starting y_shift -- only used here for calculations
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
  } else{ # otherwise forcing is likely
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

  # restricting ac coeff to not be overexpressed/repressed
  if (gam0 < over_cut){
    gam0 <- -over_cut
  } else if (gam0 > over_cut){
    gam0 <- over_cut
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
  # higher fidelity guess for higher amounts of points
  # this still works for 3 points
  # beginning
  beg1 <- which(!is.na(y_val))[2]
  beg2 <- which(!is.na(y_val))[3]
  # middle
  mid1 <- which(!is.na(y_val))[floor(length(timen)/2)]
  mid2 <- which(!is.na(y_val))[floor(length(timen)/2)+1]
  # end
  en1 <- which(!is.na(y_val))[length(timen)-2]
  en2 <- which(!is.na(y_val))[length(timen)-1]

  min_vect <- rep(10000, length(0:11))
  for (i in 0:11){
    # if the phase shift at both the second and third gene values are the smallest value available for the fitted value
    min_vect[i+1] <- sum(abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[beg1])-y_val[beg1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[mid1])-y_val[mid1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[en1])-y_val[en1]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[beg2])-y_val[beg2]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[mid2])-y_val[mid2]),
                         abs(alt_form(a0,gam0,w0,(i*pi/6),y0,timen[en2])-y_val[en2]))
  }
  phi0 <- (which.min(abs(min_vect))-1)*pi/6 # intial value for phase shift

  # form all starting parameters into a nice list
  start_param <- list(gam=gam0,a=a0,omega=w0,phi=phi0,slo=slo0,y_shift=y_shift0)

  return(start_param)
}

# function to categorize ac coefficient into categories
# inputs:
#  gam: ac coefficient
#  harm_cut: postive number indicating the cutoff for a gene to be considered harmonic
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
# outputs:
#  string, ac coefficient categories
get_type_gam <- function(gam, harm_cut, over_cut){
  # categorize gamma in ac coeff categories
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
  return(type_gam)
}

# function to calculate hours shifted (phase shift, in hours) for ECHO-type rhythms
# inputs:
#  a: initial amplitude
#  phi: phase shift, radians
#  omega: radian frequency
# outputs:
#  numeric, phase shift, in hours, relative to the 0 specified in the time course
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

# function for determining whether gene values are unexpressed (less than 70% expressed (i.e., below specified cutoff)) for full matrix
# inputs:
#  rem_unexpr_amt: what percentage of the time course must be expressed
#  rem_unexpr_amt_below: cutoff for expression
# outputs:
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

# function for getting an object of class nlsModel, EDITED VERY SPECIFICALLY FOR MOSAIC PROBLEM
# inputs:
#  form: formula specifying model
#  data: dataframe containing the data to be used
#  start: named initial values for parameters
#  wts: weights for data
# outputs:
#  object of class nlsModel
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

# function for getting an object of class nls, EDITED VERY SPECIFICIALLY FOR MOSAIC PROBLEM
# inputs:
#  formula: a nonlinear model formula including variables and parameters. Will be coerced to a formula if necessary.
#  data: an optional data frame in which to evaluate the variables in formula and weights. Can also be a list or an environment, but not a matrix.
# start: a named list or named numeric vector of starting estimates. When start is missing (and formula is not a self-starting model, see selfStart), a very cheap guess for start is tried (if algorithm != "plinear").
# control: an optional list of control settings. See nls.control for the names of the settable control values and their effect.
# algorithm: character string specifying the algorithm to use. The default algorithm is a Gauss-Newton algorithm. Other possible values are "plinear" for the Golub-Pereyra algorithm for partially linear least-squares models and "port" for the 'nl2sol' algorithm from the Port library - see the references. Can be abbreviated.
# trace: logical value indicating if a trace of the iteration progress should be printed. Default is FALSE. If TRUE the residual (weighted) sum-of-squares and the parameter values are printed at the conclusion of each iteration. When the "plinear" algorithm is used, the conditional estimates of the linear parameters are printed after the nonlinear parameters. When the "port" algorithm is used the objective function value printed is half the residual (weighted) sum-of-squares.
# subset: an optional vector specifying a subset of observations to be used in the fitting process.
# weights: an optional numeric vector of (fixed) weights. When present, the objective function is weighted least squares.
# na.action: a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset. The 'factory-fresh' default is na.omit. Value na.exclude can be useful.
# model: logical. If true, the model frame is returned as part of the object. Default is FALSE.
# lower, upper: vectors of lower and upper bounds, replicated to be as long as start. If unspecified, all parameters are assumed to be unconstrained. Bounds can only be used with the "port" algorithm. They are ignored, with a warning, if given for other algorithms.
# ...: Additional optional arguments. None are used at present.
# outputs:
#  a list containing:
# m: an nlsModel object incorporating the model.
# data: the expression that was passed to nls as the data argument. The actual data values are present in the environment of the m component.
# call: the matched call with several components, notably algorithm.
# na.action: the "na.action" attribute (if any) of the model frame.
# dataClasses: the "dataClasses" attribute (if any) of the "terms" attribute of the model frame.
# model: if model = TRUE, the model frame.
# weights: if weights is supplied, the weights.
# convInfo: a list with convergence information.
# control: the control list used, see the control argument.
# convergence, message: for an algorithm = "port" fit only, a convergence code (0 for convergence) and message. To use these is deprecated, as they are available from convInfo now.
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


# calculate model parameters ----

# function to calculate echo joint model parameters with inequality constraints (slack for protein)
# inputs:
#  current_gene: index for gene that we're calculating parameters for
#  timen: time points vector
#  resol: resolution
#  num_reps: number of replicates
#  start_type: "orig", only calculate starting points from scratch
#  res: vector of parameter estimates from a completely joint echo model
#  harm_cut: postive number indicating the cutoff for a gene to be considered harmonic
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
#  genes_rna: dataframe of rna gene expressions (first column is names, other columns are numeric)
#  genes_protein: dataframe of protein gene expressions (first column is names, other columns are numeric)
#  avg_genes_rna: matrix of averaged gene expressions, over replicates, for rna
#  avg_genes_pro: matrix of averaged gene expressions, over replicates, for protein
# outputs:
#  list which contains: model name, fit model object, parameters, joint p value, rna p value, protein p value, data frame of parameters for output, ac coefficient category for rna and protein
calc_param_echo_joint_constr <- function(current_gene, timen, resol, num_reps, start_type = "orig", res, harm_cut, over_cut, genes_rna, genes_pro, avg_genes_rna, avg_genes_pro){
  gene_n <- as.character(genes_rna[current_gene,1]) # gene name

  # calculating averages for initial values
  if (num_reps == 1){
    y_val_rna <- as.numeric(genes_rna[current_gene,-1])
    y_val_pro <- as.numeric(genes_pro[current_gene,-1]) # all the gene values
  } else{
    y_val_pro <- avg_genes_pro[current_gene,] # starting values determined by average of replicates
    y_val_rna <- avg_genes_rna[current_gene,] # starting values determined by average of replicates
  }

  # calculate first starting points for rna and protein, using the original method ("orig")
  rna_start <- calc_start_echo(y_val_rna, over_cut)
  pro_start <- calc_start_echo(y_val_pro, over_cut)

  # we're going to use the starting value of rna for the joint period, respective for everything else
  # now we need to format it into a way that is nice (renaming)
  names(rna_start) <- paste0(names(rna_start),"1")
  names(pro_start) <- paste0(names(pro_start),"2")
  start_param <- c(rna_start,pro_start)
  names(start_param)[names(start_param)=="omega2"] <- "chang"
  # changing starting period to be in hours
  start_param$omega1 <- 2*pi/start_param$omega1
  start_param$chang <- 0

  # create second starting points -- based on the results of a completely joint model
  start_joint <- start_param
  # change the values
  start_joint[names(start_joint)!="chang"] <- as.list(res[,names(start_joint)[names(start_joint)!="chang"]])
  # changing starting joint period to be in hours
  start_joint$omega1 <- 2*pi/start_joint$omega1

  # putting together expression values and time points
  y_stack <- c(unlist(genes_rna[current_gene,-1]), unlist(genes_pro[current_gene,-1]))
  t_stack <- c(rep(timen, each=num_reps), (rep(timen, each = num_reps)))
  # calculate weights for time points, if there are replicates
  if (num_reps > 1){
    weights <- c(rep(calc_weights_boot(unlist(genes_rna[current_gene,-1]), num_reps), each = num_reps),
                 rep(calc_weights_boot(unlist(genes_pro[current_gene,-1]), num_reps), each = num_reps))
  }

  # time points, repeated for each replicate
  t_rep <- length(rep(timen, each=num_reps))
  # binary vectors to indicate rna and protein in concatenation
  s1 <- rep(c(1,0), c(t_rep,t_rep)) # rna
  s2 <- rep(c(0,1), c(t_rep,t_rep)) # protein

  # constraints are in the order of the start parameters
  # preallocate constraints
  lower_constr <- upper_constr <- c(rep(1, length(start_param)))
  # name the constraints
  names(lower_constr) <- names(upper_constr) <- names(start_param)
  # assign upper and lower constraints
  lower_constr <- c( -over_cut,-Inf, 2*pi/low, -Inf, min(y_stack[1:(t_rep)]),
                     -over_cut,-Inf, (-resol+1e-12), -Inf, min(y_stack[((t_rep)+1):length(y_stack)]))
  upper_constr <- c( over_cut, Inf, 2*pi/high, Inf, max(y_stack[1:t_rep]),
                     over_cut, Inf, resol-1e-12, Inf, max(y_stack[((t_rep)+1):length(y_stack)]))

  # preallocate
  temp <- data.frame() # df for data
  coeff <- c() # coefficient
  best <- "" # best model
  num_iter <- 1000 # max number of iterations
  if (num_reps == 1){ # one replicate
    #put the timen and data point into a data frame
    temp <- data.frame(y = y_stack, t =t_stack, s1 = s1, s2 = s2, resol = rep(resol,length(y_stack)))
    temp <- temp[!is.na(temp$y),] # remove any missing data points

    # first: fitting with original starting points
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_one_rep,
                              fcall = echo_joint_res_constr,
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
    iter.orig <- oscillator_init$niter # extract amount of iterations

    # put coefficient estimates into an nls object
    oscillator.fit.orig <- nls_edit(
      formula = y ~ echo_joint_res_constr(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, y_shift1,y_shift2, t, s1, s2, resol),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], y_shift1=coeff[5],
                   gam2=coeff[6], a2=coeff[7], chang = coeff[8], phi2 = coeff[9], y_shift2=coeff[10]
      ),
      control = nls.control(maxiter = 0)
    )

    # second: fitting with starting points from a completely joint model
    oscillator_init <- nls.lm(par = start_joint,
                              fn = fcn_one_rep,
                              fcall = echo_joint_res_constr,
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
    iter.joint <- oscillator_init$niter # extract amount of iterations

    # put coefficient estimates into an nls object
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
    temp <- temp[!is.na(temp$y),] # remove any missing data points

    # first: fitting with original starting points
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn,
                              fcall = echo_joint_res_constr,
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
    iter.orig <- oscillator_init$niter # extract amount of iterations

    # put coefficient estimates into an nls object
    oscillator.fit.orig <- nls_edit(
      formula = y ~ echo_joint_res_constr(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, y_shift1,y_shift2, t, s1, s2, resol),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], y_shift1=coeff[5],
                   gam2=coeff[6], a2=coeff[7], chang = coeff[8], phi2 = coeff[9], y_shift2=coeff[10]
      ),
      control = nls.control(maxiter = 0)
    )

    # second: fitting with starting points from a completely joint model
    oscillator_init <- nls.lm(par = start_joint,
                              fn = fcn,
                              fcall = echo_joint_res_constr,
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
    iter.joint <- oscillator_init$niter # extract amount of iterations

    # put coefficient estimates into an nls object
    oscillator.fit.joint <- nls_edit(
      formula = y ~ echo_joint_res_constr(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, y_shift1,y_shift2, t, s1, s2, resol),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], y_shift1=coeff[5],
                   gam2=coeff[6], a2=coeff[7], chang = coeff[8], phi2 = coeff[9], y_shift2=coeff[10]
      ),
      control = nls.control(maxiter = 0)
    )

  }

  # compare the results of the starting points using aic
  aic_res <- AIC(oscillator.fit.orig, oscillator.fit.joint)
  ord_start <- c("orig","joint") # order of starting point fits
  # choose fit based on which has the lowest aic
  oscillator.fit <- list(oscillator.fit.orig, oscillator.fit.joint)[[which.min(aic_res$AIC)]] # nls object
  best <- ord_start[which.min(aic_res$AIC)] # best starting points
  coeff <- list(coeff.orig, coeff.joint)[[which.min(aic_res$AIC)]] # coefficient
  num_iter <- c(iter.orig, iter.joint)[which.min(aic_res$AIC)] # number of iterations

  # coefficients and name, with some other options
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
  type_gam_rna <- get_type_gam(coeff[1], harm_cut, over_cut)
  type_gam_pro <- get_type_gam(coeff[6], harm_cut, over_cut)

  # get hours shifted
  phase_hours_rna <- get_hs(coeff[2],coeff[4],2*pi/coeff[3])
  phase_hours_pro <- get_hs(coeff[7],coeff[9],2*pi/(coeff[3]+coeff[8]))

  # final df row vector
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
  # add those name
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

  # results list
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

# function to calculate echo linear joint model parameters with inequality constraints (slack for protein)
# inputs:
#  current_gene: index for gene that we're calculating parameters for
#  timen: time points vector
#  resol: resolution
#  num_reps: number of replicates
#  start_type: "orig", only calculate starting points from scratch
#  res: vector of parameter estimates from a completely joint echo linear model
#  harm_cut: postive number indicating the cutoff for a gene to be considered harmonic
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
#  genes_rna: dataframe of rna gene expressions (first column is names, other columns are numeric)
#  genes_protein: dataframe of protein gene expressions (first column is names, other columns are numeric)
#  avg_genes_rna: matrix of averaged gene expressions, over replicates, for rna
#  avg_genes_pro: matrix of averaged gene expressions, over replicates, for protein
# outputs:
#  list which contains: model name, fit model object, parameters, joint p value, rna p value, protein p value, data frame of parameters for output, ac coefficient category for rna and protein
calc_param_echo_lin_joint_constr <- function(current_gene, timen, resol, num_reps, start_type = "orig", res, harm_cut, over_cut, genes_rna, genes_pro, avg_genes_rna, avg_genes_pro){
  gene_n <- as.character(genes_rna[current_gene,1]) # gene name

  # calculating averages for initial values
  if (num_reps == 1){
    y_val_rna <- as.numeric(genes_rna[current_gene,-1])
    y_val_pro <- as.numeric(genes_pro[current_gene,-1]) # all the gene values
  } else{
    y_val_pro <- avg_genes_pro[current_gene,] # starting values determined by average of replicates
    y_val_rna <- avg_genes_rna[current_gene,] # starting values determined by average of replicates
  }

  # calculating linear model to remove for starting values for rna and protein
  # rna
  # calculate linear model parameters
  lin.fit <- lm(as.numeric(genes_rna[current_gene,-1])~(rep(timen, each = num_reps)))
  coeff <- as.numeric(lin.fit$coefficients)
  y_shift0_rna <- coeff[1]
  slo0_rna <- coeff[2]
  # protein
  # calculate linear model parameters
  lin.fit <- lm(as.numeric(genes_pro[current_gene,-1])~(rep(timen, each = num_reps)))
  coeff <- as.numeric(lin.fit$coefficients)
  y_shift0_pro <- coeff[1]
  slo0_pro <- coeff[2]

  # calculate starting points using the original method
  rna_start <- calc_start_echo_lin(y_val_rna, slo0_rna, y_shift0_rna, over_cut)
  pro_start <- calc_start_echo_lin(y_val_pro, slo0_pro, y_shift0_pro, over_cut)

  # we're going to use the starting value of rna for the joint period, respective for everything else
  # now we need to format it into a way that is nice
  names(rna_start) <- paste0(names(rna_start),"1")
  names(pro_start) <- paste0(names(pro_start),"2")
  start_param <- c(rna_start,pro_start)
  names(start_param)[names(start_param)=="omega2"] <- "chang"
  # changing starting period to be in hours
  start_param$omega1 <- 2*pi/start_param$omega1
  start_param$chang <- 0

  # create second starting points -- based on the results of a completely joint model
  start_joint <- start_param
  # change the values
  start_joint[names(start_joint)!="chang"] <- as.list(res[,names(start_joint)[names(start_joint)!="chang"]])
  # changing starting joint period to be in hours
  start_joint$omega1 <- 2*pi/start_joint$omega1

  # putting together expression values and time points
  y_stack <- c(unlist(genes_rna[current_gene,-1]), unlist(genes_pro[current_gene,-1]))
  t_stack <- c(rep(timen, each=num_reps), (rep(timen, each = num_reps)))
  # calculate weights for time points, if there are replicates
  if (num_reps > 1){
    weights <- c(rep(calc_weights_boot(unlist(genes_rna[current_gene,-1]), num_reps), each = num_reps),
                 rep(calc_weights_boot(unlist(genes_pro[current_gene,-1]), num_reps), each = num_reps))
  }

  # time points, repeated for each replicate
  t_rep <- length(rep(timen, each=num_reps))
  # binary vectors to indicate rna and protein in concatenation
  s1 <- rep(c(1,0), c(t_rep,t_rep)) # rna
  s2 <- rep(c(0,1), c(t_rep,t_rep)) # protein

  # constraints are in the order of the start parameters
  # preallocate constraints
  lower_constr <- upper_constr <- c(rep(1, length(start_param)))
  # name the constraints
  names(lower_constr) <- names(upper_constr) <- names(start_param)
  # assign upper and lower constraints
  lower_constr <- c( -over_cut,-Inf, 2*pi/low, -Inf, -Inf, min(y_stack[1:(t_rep)]),
                     -over_cut,-Inf, (-resol+1e-12), -Inf, -Inf, min(y_stack[((t_rep)+1):length(y_stack)]))
  upper_constr <- c( over_cut, Inf, 2*pi/high, Inf, Inf, max(y_stack[1:t_rep]),
                     over_cut, Inf, resol-1e-12, Inf, Inf, max(y_stack[((t_rep)+1):length(y_stack)]))

  # preallocate
  temp <- data.frame() # df for data
  coeff <- c() # coefficient
  best <- "" # best model
  num_iter <- 1000 # max number of iterations
  if (num_reps == 1){ # one replicate
    #put the timen and data point into a data frame
    temp <- data.frame(y = y_stack, t =t_stack, s1 = s1, s2 = s2)
    temp <- temp[!is.na(temp$y),] # remove any missing data points

    # first: fitting with original starting points
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_one_rep,
                              fcall = echo_lin_joint_res_constr,
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
    iter.orig <- oscillator_init$niter # extract amount of iterations

    # put coefficient estimates into an nls object
    oscillator.fit.orig <- nls_edit(
      formula = y ~ echo_lin_joint_res_constr(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, slo1,slo2, y_shift1,y_shift2, t, s1, s2, resol),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], slo1 = coeff[5], y_shift1=coeff[6],
                   gam2=coeff[7], a2=coeff[8], chang = coeff[9], phi2 = coeff[10],  slo2 = coeff[11], y_shift2=coeff[12]
      ),
      control = nls.control(maxiter = 0)
    )

    # second: fitting with starting points from a completely joint model
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
    iter.joint <- oscillator_init$niter # extract amount of iterations

    # put coefficient estimates into an nls object
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
    temp <- temp[!is.na(temp$y),] # remove any missing data points

    # first: fitting with original starting points
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
    iter.orig <- oscillator_init$niter # extract amount of iterations

    # put coefficient estimates into an nls object
    oscillator.fit.orig <- nls_edit(
      formula = y ~ echo_lin_joint_res_constr(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, slo1,slo2, y_shift1,y_shift2, t, s1, s2, resol),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], slo1 = coeff[5], y_shift1=coeff[6],
                   gam2=coeff[7], a2=coeff[8], chang = coeff[9], phi2 = coeff[10],  slo2 = coeff[11], y_shift2=coeff[12]
      ),
      control = nls.control(maxiter = 0)
    )

    # second: fitting with starting points from a completely joint model
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
    iter.joint <- oscillator_init$niter # extract amount of iterations

    # put coefficient estimates into an nls object
    oscillator.fit.joint <- nls_edit(
      formula = y ~ echo_lin_joint_res_constr(a1,a2, gam1,gam2, omega1,chang, phi1,phi2, slo1,slo2, y_shift1,y_shift2, t, s1, s2, resol),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], slo1 = coeff[5], y_shift1=coeff[6],
                   gam2=coeff[7], a2=coeff[8], chang = coeff[9], phi2 = coeff[10],  slo2 = coeff[11], y_shift2=coeff[12]
      ),
      control = nls.control(maxiter = 0)
    )
  }

  # compare the results of the starting points using aic
  aic_res <- AIC(oscillator.fit.orig, oscillator.fit.joint)
  ord_start <- c("orig","joint") # order of starting point fits
  # choose fit based on which has the lowest aic
  oscillator.fit <- list(oscillator.fit.orig, oscillator.fit.joint)[[which.min(aic_res$AIC)]] # nls object
  best <- ord_start[which.min(aic_res$AIC)] # best starting points
  coeff <- list(coeff.orig, coeff.joint)[[which.min(aic_res$AIC)]] # coefficient
  num_iter <- c(iter.orig, iter.joint)[which.min(aic_res$AIC)] # number of iterations

  # coefficients and name, with some other options
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
  type_gam_rna <- get_type_gam(coeff[1], harm_cut, over_cut)
  type_gam_pro <- get_type_gam(coeff[7], harm_cut, over_cut)

  # get hours shifted
  phase_hours_rna <- get_hs(coeff[2],coeff[4],2*pi/coeff[3])
  phase_hours_pro <- get_hs(coeff[8],coeff[10],2*pi/(coeff[3]+coeff[9]))

  # final df row vector
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

  # adjust names correctly
  echo_lin_param_n <- c("Initial_Amplitude", "AC_Coefficient", "Oscillation_Type", "Radian_Frequency", "Period", "Phase_Shift", "Hours_Shifted", "Slope", "Equilibrium_Value")
  echo_lin_joint_param_n <- c(paste0(echo_lin_param_n, "_RNA"), paste0(echo_lin_param_n,"_Protein"))
  # add those names
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

  # results list
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

# function to calculate echo joint model parameters for a completely joint model (no slack for protein), for the purpose of getting starting points for the slacked model
# inputs:
#  current_gene: index for gene that we're calculating parameters for
#  timen: time points vector
#  resol: resolution
#  num_reps: number of replicates
#  start_type: "orig", only calculate starting points from scratch
#  harm_cut: postive number indicating the cutoff for a gene to be considered harmonic
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
#  genes_rna: dataframe of rna gene expressions (first column is names, other columns are numeric)
#  genes_protein: dataframe of protein gene expressions (first column is names, other columns are numeric)
#  avg_genes_rna: matrix of averaged gene expressions, over replicates, for rna
#  avg_genes_pro: matrix of averaged gene expressions, over replicates, for protein
# outputs:
#  list which contains: model name, fit model object, and parameters
calc_param_echo_joint <- function(current_gene, timen, resol, num_reps, start_type = "orig", harm_cut, over_cut, genes_rna, genes_pro, avg_genes_rna, avg_genes_pro){
  gene_n <- as.character(genes_rna[current_gene,1]) # gene name

  # calculating averages for initial values
  if (num_reps == 1){
    y_val_rna <- as.numeric(genes_rna[current_gene,-1])
    y_val_pro <- as.numeric(genes_pro[current_gene,-1]) # all the gene values
  } else{
    y_val_pro <- avg_genes_pro[current_gene,] # starting values determined by average of replicates
    y_val_rna <- avg_genes_rna[current_gene,] # starting values determined by average of replicates
  }

  # calculate first starting points for rna and protein, using the original method ("orig")
  rna_start <- calc_start_echo(y_val_rna, over_cut)
  pro_start <- calc_start_echo(y_val_pro, over_cut)

  # we're going to use the starting value of rna for the joint period, respective for everything else
  # now we need to format it into a way that is nice (renaming)
  names(rna_start) <- paste0(names(rna_start),"1")
  names(pro_start) <- paste0(names(pro_start),"2")
  start_param <- c(rna_start,pro_start)
  start_param$omega2 <- NULL

  # putting together expression values and time points
  y_stack <- c(unlist(genes_rna[current_gene,-1]), unlist(genes_pro[current_gene,-1]))
  t_stack <- c(rep(timen, each=num_reps), (rep(timen, each = num_reps)))
  # calculate weights for time points, if there are replicates
  if (num_reps > 1){
    weights <- c(rep(calc_weights_boot(unlist(genes_rna[current_gene,-1]), num_reps), each = num_reps),
                 rep(calc_weights_boot(unlist(genes_pro[current_gene,-1]), num_reps), each = num_reps))
  }

  # time points, repeated for each replicate
  t_rep <- length(rep(timen, each=num_reps))
  # binary vectors to indicate rna and protein in concatenation
  s1 <- rep(c(1,0), c(t_rep,t_rep)) # rna
  s2 <- rep(c(0,1), c(t_rep,t_rep)) # protein

  # constraints are in the order of the start parameters
  # preallocate constraints
  lower_constr <- upper_constr <- c(rep(1, length(start_param)))
  # assign upper and lower constraints
  names(lower_constr) <- names(upper_constr) <- names(start_param)
  lower_constr <- c( -over_cut,-Inf, high, -Inf, min(y_stack[1:(t_rep)]),
                     -over_cut,-Inf, -Inf, min(y_stack[((t_rep)+1):length(y_stack)]))
  upper_constr <- c( over_cut, Inf, low, Inf, max(y_stack[1:t_rep]),
                     over_cut, Inf, Inf, max(y_stack[((t_rep)+1):length(y_stack)]))

  # preallocate
  temp <- data.frame() # df for data
  coeff <- c() # coefficient
  if (num_reps == 1){ # one replicate
    #put the timen and data point into a data frame
    temp <- data.frame(y = y_stack, t =t_stack, s1 = s1, s2 = s2)
    temp <- temp[!is.na(temp$y),] # remove any missing data points

    # fitting with original starting points
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_one_rep,
                              fcall = echo_joint,
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

    # put coefficient estimates into an nls object
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
    temp <- temp[!is.na(temp$y),] # remove any missing data points

    # fitting with original starting points
    oscillator_init <- nls.lm(par = start_param,
                             fn = fcn,
                             fcall = echo_joint,
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

    # put coefficient estimates into an nls object
    oscillator.fit <- nls_edit(
      formula = y ~ echo_joint(a1,a2, gam1,gam2, omega1, phi1,phi2, y_shift1,y_shift2, t, s1, s2),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], y_shift1=coeff[5],
                   gam2=coeff[6], a2=coeff[7], phi2 = coeff[8], y_shift2=coeff[9]
      ),
      control = nls.control(maxiter = 0)
    )
  }

  # parameters
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

  # results list
  res <- list("mod_name" = "echo_joint",
              "model" = oscillator.fit,
              "param" = res)
  return(res)
}

# function to calculate echo linear joint model parameters for a completely joint model (no slack for protein), for the purpose of getting starting points for the slacked model
# inputs:
#  current_gene: index for gene that we're calculating parameters for
#  timen: time points vector
#  resol: resolution
#  num_reps: number of replicates
#  start_type: "orig", only calculate starting points from scratch
#  harm_cut: postive number indicating the cutoff for a gene to be considered harmonic
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
#  genes_rna: dataframe of rna gene expressions (first column is names, other columns are numeric)
#  genes_protein: dataframe of protein gene expressions (first column is names, other columns are numeric)
#  avg_genes_rna: matrix of averaged gene expressions, over replicates, for rna
#  avg_genes_pro: matrix of averaged gene expressions, over replicates, for protein
# outputs:
#  list which contains: model name, fit model object, and parameters
calc_param_echo_lin_joint <- function(current_gene, timen, resol, num_reps, start_type = "orig", harm_cut, over_cut, genes_rna, genes_pro, avg_genes_rna, avg_genes_pro){
  gene_n <- as.character(genes_rna[current_gene,1]) # gene name

  # calculating averages for initial values
  if (num_reps == 1){
    y_val_rna <- as.numeric(genes_rna[current_gene,-1])
    y_val_pro <- as.numeric(genes_pro[current_gene,-1]) # all the gene values
  } else{
    y_val_pro <- avg_genes_pro[current_gene,] # starting values determined by average of replicates
    y_val_rna <- avg_genes_rna[current_gene,] # starting values determined by average of replicates
  }

  # calculating linear model to remove for starting values for rna and protein
  # rna
  # calculate linear model parameters
  lin.fit <- lm(as.numeric(genes_rna[current_gene,-1])~(rep(timen, each = num_reps)))
  coeff <- as.numeric(lin.fit$coefficients)
  y_shift0_rna <- coeff[1]
  slo0_rna <- coeff[2]
  # protein
  # calculate linear model parameters
  lin.fit <- lm(as.numeric(genes_pro[current_gene,-1])~(rep(timen, each = num_reps)))
  coeff <- as.numeric(lin.fit$coefficients)
  y_shift0_pro <- coeff[1]
  slo0_pro <- coeff[2]

  # calculate starting points using the original method
  rna_start <- calc_start_echo_lin(y_val_rna, slo0_rna, y_shift0_rna, over_cut)
  pro_start <- calc_start_echo_lin(y_val_pro, slo0_pro, y_shift0_pro, over_cut)

  # we're going to use the starting value of rna for the joint period, respective for everything else
  # now we need to format it into a way that is nice
  names(rna_start) <- paste0(names(rna_start),"1")
  names(pro_start) <- paste0(names(pro_start),"2")
  start_param <- c(rna_start,pro_start)
  start_param$omega2 <- NULL

  # putting together expression values and time points
  y_stack <- c(unlist(genes_rna[current_gene,-1]), unlist(genes_pro[current_gene,-1]))
  t_stack <- c(rep(timen, each=num_reps), (rep(timen, each = num_reps)))
  # calculate weights for time points, if there are replicates
  if (num_reps > 1){
    weights <- c(rep(calc_weights_boot(unlist(genes_rna[current_gene,-1]), num_reps), each = num_reps),
                 rep(calc_weights_boot(unlist(genes_pro[current_gene,-1]), num_reps), each = num_reps))
  }

  # time points, repeated for each replicate
  t_rep <- length(rep(timen, each=num_reps))
  # binary vectors to indicate rna and protein in concatenation
  s1 <- rep(c(1,0), c(t_rep,t_rep)) # rna
  s2 <- rep(c(0,1), c(t_rep,t_rep)) # protein

  # constraints are in the order of the start parameters
  # preallocate constraints
  lower_constr <- upper_constr <- c(rep(1, length(start_param)))
  # name the constraints
  names(lower_constr) <- names(upper_constr) <- names(start_param)
  # assign upper and lower constraints
  lower_constr <- c( -over_cut,-Inf, high, -Inf, -Inf,  min(y_stack[1:(t_rep)]),
                     -over_cut,-Inf, -Inf, -Inf, min(y_stack[((t_rep)+1):length(y_stack)]))
  upper_constr <- c( over_cut, Inf, low, Inf, Inf, max(y_stack[1:t_rep]),
                     over_cut, Inf, Inf, Inf, max(y_stack[((t_rep)+1):length(y_stack)]))

  # preallocate
  temp <- data.frame() # df for data
  coeff <- c() # coefficient
  if (num_reps == 1){ # one replicate
    #put the timen and data point into a data frame
    temp <- data.frame(y = y_stack, t =t_stack, s1 = s1, s2 = s2)
    temp <- temp[!is.na(temp$y),] # remove any missing data points

    # fitting with original starting points
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_one_rep,
                              fcall = echo_lin_joint,
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

    # put coefficient estimates into an nls object
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
    temp <- temp[!is.na(temp$y),] # remove any missing data points

    # fitting with original starting points
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn,
                              fcall = echo_lin_joint,
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

    # put coefficient estimates into an nls object
    oscillator.fit <- nls_edit(
      formula = y ~ echo_lin_joint(a1,a2, gam1,gam2, omega1, phi1,phi2, slo1,slo2, y_shift1,y_shift2, t, s1, s2),
      data=temp,
      start = list(gam1=coeff[1], a1=coeff[2], omega1 = coeff[3], phi1 = coeff[4], slo1 = coeff[5], y_shift1=coeff[6],
                   gam2=coeff[7], a2=coeff[8], phi2 = coeff[9], slo2 = coeff[10], y_shift2=coeff[11]
      ),
      control = nls.control(maxiter = 0)
    )
  }

  # parameters
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

  # results list
  res <- list("mod_name" = "echo_lin_joint",
              "model" = oscillator.fit,
              "param" = res)
  return(res)
}

# function to calculate echo model parameters
# inputs:
#  current_gene: index for gene that we're calculating parameters for
#  timen: time points vector
#  resol: resolution
#  num_reps: number of replicates
#  genes: dataframe of gene expressions (first column is names, other columns are numeric)
#  avg_genes: matrix of averaged gene expressions, over replicates
#  harm_cut: postive number indicating the cutoff for a gene to be considered harmonic
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
# outputs:
#  list which contains: model name, fit model object, parameters, p value, data frame of parameters for output, ac coefficient category
calc_param_echo <- function(current_gene, timen, resol, num_reps, genes, avg_genes, harm_cut, over_cut, genes_rna, genes_pro, avg_genes_rna, avg_genes_pro){
  gene_n <- as.character(genes[current_gene,1]) # gene name

  # calculating averages for initial values
  if (num_reps == 1){
    y_val <- as.numeric(genes[current_gene,-1]) # all the gene values
  } else{
    y_val <- avg_genes[current_gene,] # starting values determined by average of replicates
  }

  # calculate starting values
  start_param <- calc_start_echo(y_val, over_cut)

  # preallocate df for data
  temp <- data.frame()
  if (num_reps == 1){ # one replicate
    # put the timen and data points into a data frame
    temp <- data.frame(y=y_val,t=timen)
    temp <- temp[!is.na(temp$y),] # remove any missing data points

    # fitting
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_single_one_rep,
                              fcall = alt_form,
                              lower=c(-over_cut, -Inf, high, -Inf, min(temp$y)),
                              upper=c(over_cut, Inf, low, Inf, max(temp$y)),
                              t = temp$t,
                              y = temp$y,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )

    # get coefficients
    coeff <- coef(oscillator_init)

    # put coefficient estimates into an nls object
    oscillator.fit <- nls_edit(
      formula = y ~ alt_form(a, gam, omega, phi, y_shift, t),
      data=temp,
      start = list(gam=coeff[1], a=coeff[2], omega = coeff[3], phi = coeff[4], y_shift=coeff[5]),
      control = nls.control(maxiter = 0)
    )

  } else{ # multiple replicates
    #put the timen and data points into a data frame, and weights
    weights <- calc_weights(current_gene,num_reps, genes)
    temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),t=cbind(rep(timen,each = num_reps)),w=cbind(rep(weights,each = num_reps)))
    temp <- temp[!is.na(temp$y),] # remove any missing data points

    # fitting
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_single,
                              fcall = alt_form,
                              lower=c(-over_cut, -Inf, high, -Inf, min(temp$y)),
                              upper=c(over_cut, Inf, low, Inf, max(temp$y)),
                              t = temp$t,
                              y = temp$y,
                              w = temp$w,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )

    # get coefficients
    coeff <- coef(oscillator_init)

    # put coefficient estimates into an nls object
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
  type_gam <- get_type_gam(coeff[1], harm_cut, over_cut)

  # get hours shifted
  phase_hours <- get_hs(coeff[2],coeff[4],coeff[3])

  # make final dataframe row
  fin <- c(coeff[2], coeff[1], type_gam, coeff[3], 2*pi/coeff[3], coeff[4], phase_hours, coeff[5])
  echo_param_n <- c("Initial_Amplitude", "AC_Coefficient", "Oscillation_Type", "Radian_Frequency", "Period", "Phase_Shift", "Hours_Shifted", "Equilibrium_Value")
  names(fin) <- echo_param_n

  # generating p-values
  all_pred <- oscillator.fit$m$predict() # fitted values
  testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
  pval <- testing$p.value

  # results list
  res <- list("mod_name" = "echo",
              "model" = oscillator.fit,
              "param" = res,
              "pval" = pval,
              "final" = fin,
              "type_gam" = type_gam)

  return(res)
}

# function to calculate echo linear model parameters
# inputs:
#  current_gene: index for gene that we're calculating parameters for
#  timen: time points vector
#  resol: resolution
#  num_reps: number of replicates
#  genes: dataframe of gene expressions (first column is names, other columns are numeric)
#  avg_genes: matrix of averaged gene expressions, over replicates
#  harm_cut: postive number indicating the cutoff for a gene to be considered harmonic
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
# outputs:
#  list which contains: model name, fit model object, parameters, p value, data frame of parameters for output, ac coefficient category
calc_param_echo_lin <- function(current_gene, timen, resol, num_reps, genes, avg_genes, harm_cut, over_cut){
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
  start_param <- calc_start_echo_lin(y_val, slo0, y_shift0, over_cut)

  # preallocate df for data
  temp <- data.frame()
  if (num_reps == 1){ # one replicate
    # put the timen into a data frame
    temp <- data.frame(y=y_val,t=timen)
    temp <- temp[!is.na(temp$y),] # remove any missing data points

    # fitting
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_single_one_rep,
                              fcall = echo_lin_mod,
                              lower=c(-over_cut, -Inf, high, -Inf, -Inf, min(temp$y)),
                              upper=c(over_cut, Inf, low, Inf, Inf, max(temp$y)),
                              t = temp$t,
                              y = temp$y,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )

    # get coefficients
    coeff <- coef(oscillator_init)

    # put coefficient estimates into an nls object
    oscillator.fit <- nls_edit(
      formula = y ~ echo_lin_mod(a, gam, omega, phi, slo, y_shift, t),
      data=temp,
      start = list(gam=coeff[1], a=coeff[2], omega = coeff[3], phi = coeff[4], slo = coeff[5], y_shift=coeff[6]),
      control = nls.control(maxiter = 0)
    )

  } else{ # multiple replicates (add weights)
    #put the timen and data point into a data frame
    weights <- calc_weights(current_gene,num_reps, genes)
    temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),t=cbind(rep(timen,each = num_reps)),w=cbind(rep(weights,each = num_reps)))
    temp <- temp[!is.na(temp$y),] # remove any missing data points


    #fitting
    oscillator_init <- nls.lm(par = start_param,
                              fn = fcn_single,
                              fcall = echo_lin_mod,
                              lower=c(-over_cut, -Inf, high, -Inf, -Inf, min(temp$y)),
                              upper=c(over_cut, Inf, low, Inf, Inf, max(temp$y)),
                              t = temp$t,
                              y = temp$y,
                              w = temp$w,
                              control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                       ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )

    # get coefficients
    coeff <- coef(oscillator_init)

    # put coefficient estimates into an nls object
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
  type_gam <- get_type_gam(gam, harm_cut, over_cut)

  # get hours shifted
  a <- coeff[2]
  phi <- coeff[4]
  omega <- coeff[3]
  phase_hours <- get_hs(a,phi,omega)

  # row for final data frame
  fin <- c(coeff[2], coeff[1], type_gam, coeff[3], 2*pi/coeff[3], coeff[4], phase_hours, coeff[5], coeff[6])
  names(fin) <- echo_lin_param_n <- c("Initial_Amplitude", "AC_Coefficient", "Oscillation_Type", "Radian_Frequency", "Period", "Phase_Shift", "Hours_Shifted", "Slope", "Equilibrium_Value")

  # generating p-values
  all_pred <- oscillator.fit$m$predict() # fitted values
  testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
  pval <- testing$p.value

  # results list
  res <- list("mod_name" = "echo_lin",
              "model" = oscillator.fit,
              "param" = res,
              "pval" = pval,
              "final" = fin,
              "type_gam" = type_gam)

  return(res)
}

# function to calculate linear model parameters
# inputs:
#  current_gene: index for gene that we're calculating parameters for
#  timen: time points vector
#  resol: resolution
#  num_reps: number of replicates
#  genes: dataframe of gene expressions (first column is names, other columns are numeric)
#  avg_genes: matrix of averaged gene expressions, over replicates
# outputs:
#  list which contains: model name, fit model object, parameters, fit p value, slope p value, data frame of parameters for output
calc_param_lin <- function(current_gene, timen, resol, num_reps, genes, avg_genes){
  # calculate linear model parameters
  # put data in a df
  temp <- data.frame(y=as.numeric(genes[current_gene,-1]),t=(rep(timen, each = num_reps)))
  temp <- temp[!is.na(temp$y),] # remove any missing data points

  # get linear fit
  lin.fit <- lm(temp$y~temp$t)
  # get coefficients
  coeff <- lin.fit$coefficients
  y_shift <- coeff[1]
  slo <- coeff[2]

  # get results dataframe
  res <- as.data.frame(t(c(slo,y_shift)))
  # add gene name
  gene_n <- genes[current_gene,1]
  res <- cbind(data.frame("gene_n"=gene_n, stringsAsFactors = F), res)

  # make final dataframe row
  fin <- c(coeff[2], coeff[1])
  names(fin) <- lin_param_n <- c("Slope", "Equilibrium_Value")

  # generating p-values
  all_pred <- lin.fit$fitted.values # fitted values
  testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
  pval <- testing$p.value

  # generate slope p-value
  sumfit <- summary(lin.fit)
  slo_pval <- sumfit$coefficients[2,4]

  # results list
  res <- list("mod_name" = "lin",
              "model" = lin.fit,
              "param" = res,
              "pval" = pval,
              "slo_pval" = slo_pval,
              "final" = fin)

  return(res)
}

# function to calculate exponential model parameters
# inputs:
#  current_gene: index for gene that we're calculating parameters for
#  timen: time points vector
#  resol: resolution
#  num_reps: number of replicates
#  genes: dataframe of gene expressions (first column is names, other columns are numeric)
#  avg_genes: matrix of averaged gene expressions, over replicates
# outputs:
#  list which contains: model name, fit model object, parameters, p value, data frame of parameters for output
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

  # fitting based on negative values (assuming negative ampltide)
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
    # set up data in df
    temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),t=cbind(rep(timen,each = num_reps)))
    temp <- temp[!is.na(temp$y),] # remove any missing data points


    #first: fitting with starting points based on positive amplitude assumption
    exp_mod_pos <- nls.lm(par = list(a=a_e0,gam=gam_e0, y_shift=y_shift_e0),
                           fn = fcn_single_one_rep,
                           fcall = exp_mod,
                           t = temp$t,
                           y = temp$y,
                           control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                    ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )

    #second: fitting with starting points based on negative amplitude assumption
    exp_mod_neg <- nls.lm(par = list(a=a_e0_neg,gam=gam_e0_neg, y_shift=y_shift_e0_neg),
                          fn = fcn_single_one_rep,
                          fcall = exp_mod,
                          t = temp$t,
                          y = temp$y,
                          control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                   ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )
  } else { #replicates
    # calculate weights
    weights <- calc_weights(current_gene,num_reps, genes)
    # set up data in df
    temp <- data.frame(y=cbind(unlist(genes[current_gene,-1])),t=cbind(rep(timen,each = num_reps)),w=cbind(rep(weights,each = num_reps)))
    temp <- temp[!is.na(temp$y),] # remove any missing data points

    #first: fitting with starting points based on positive amplitude assumption
    exp_mod_pos <- nls.lm(par = list(a=a_e0,gam=gam_e0, y_shift=y_shift_e0),
                          fn = fcn_single,
                          fcall = exp_mod,
                          t = temp$t,
                          y = temp$y,
                          w = temp$w,
                          control = nls.lm.control(nprint=0, maxiter = 1000, maxfev = 2000,
                                                   ftol=1e-6, ptol=1e-6, gtol=1e-6)
    )

    #second: fitting with starting points based on negative amplitude assumption
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

  # get coefficients
  coeff_pos <- coef(exp_mod_pos)
  coeff_neg <- coef(exp_mod_neg)

  # put coefficient estimates into an nls object
  # positive amplitude assumption
  exp_mod.fit_pos <- nls_edit(
    formula = y ~ exp_mod(a, gam, y_shift, t),
    data=temp,
    start = list(a=coeff_pos[1],gam=coeff_pos[2], y_shift=coeff_pos[3]),
    control = nls.control(maxiter = 0)
  )

  # negative amplitude assumption
  exp_mod.fit_neg <- nls_edit(
    formula = y ~ exp_mod(a, gam, y_shift, t),
    data=temp,
    start = list(a=coeff_neg[1],gam=coeff_neg[2], y_shift=coeff_neg[3]),
    control = nls.control(maxiter = 0)
  )

  # choose final fit based on lowest aic
  aic_res <- AIC(exp_mod.fit_pos, exp_mod.fit_neg)
  exp_mod.fit<- list(exp_mod.fit_pos, exp_mod.fit_neg)[[which.min(aic_res$AIC)]] # fitted object
  coeff <- list(coeff_pos, coeff_neg)[[which.min(aic_res$AIC)]] # coefficient

  gene_n <- genes[current_gene,1] # gene name
  # parameters
  res <- data.frame("gene_n"=gene_n,
                    "a" = coeff[1],
                    "gam" = coeff[2],
                    "y_shift" = coeff[3],
                    stringsAsFactors = F)

  # final df row
  fin <- coeff
  names(fin) <- exp_param_n <- c("Initial_Amplitude", "Growth_Rate", "Equilibrium_Value")

  # generating p-values
  all_pred <- exp_mod.fit$m$predict() # fitted values
  testing <- suppressWarnings(cor.test(all_pred,temp$y, method = "kendall"))
  pval <- testing$p.value

  # results list
  res <- list("mod_name" = "exp",
              "model" = exp_mod.fit,
              "param" = res,
              "pval" = pval,
              "final" = fin)

  return(res)
}

# mosaic run functions ----

# function to get the fit of the best model chosen
# inputs:
#  best_mod: a list containing the model parameters and model name
#  timen: time points vector
#  resol: resolution
# outputs:
#  vector of fitted values
get_mod_fits <- function(best_mod, timen, resol){
  # turn parameters into a vector
  prm <- suppressWarnings(as.numeric(best_mod$param))

  # get the fitted values for each model, based on the parameters
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

# function to choose the best model (echo, echo linear, linear, exponential) for each dataset (not including joint model)
# inputs:
#  current_gene: index for gene that we're calculating parameters for
#  timen: time points vector
#  resol: resolution
#  num_reps: number of replicates
#  genes: dataframe of gene expressions (first column is names, other columns are numeric)
#  avg_genes: matrix of averaged gene expressions, over replicates
#  harm_cut: postive number indicating the cutoff for a gene to be considered harmonic
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
# outputs:
#  list containing the best model object (see the calc_param functions), and the dataframe with bic values, model names, pvalues, ac category
choose_best_mod <- function(current_gene, timen, resol, num_reps, genes, avg_genes, harm_cut, over_cut){
  # model names
  model_n <- c("lin", "exp", "echo", "echo_lin")

  # values from fitting each of the models with each of the gene data
  lin_val <- calc_param_lin(current_gene, timen, resol, num_reps, genes, avg_genes)
  exp_val <- calc_param_exp(current_gene, timen, resol, num_reps, genes, avg_genes)
  echo_val <- calc_param_echo(current_gene, timen, resol, num_reps, genes, avg_genes, harm_cut, over_cut)
  echo_lin_val <- calc_param_echo_lin(current_gene, timen, resol, num_reps, genes, avg_genes, harm_cut, over_cut)

  # get the bic values for each of the fits (aic is a misnomer based on original tests)
  aic_vals <- BIC(lin_val$model, exp_val$model, echo_val$model, echo_lin_val$model)
  # get the pvalues for each of the fits
  pvals <- c(lin_val$pval, exp_val$pval, echo_val$pval, echo_lin_val$pval)

  # save results -- choose the best model based on the lowest bic
  best_mod <- list(lin_val, exp_val, echo_val, echo_lin_val)[[which.min(aic_vals$BIC)]]

  # save all the data based on the lowest bic
  aic_df <- data.frame(rbind(c(0, aic_vals$BIC, NA, NA, NA, NA, 0, 0, NA, NA)))
  aic_df[1,1] <- aic_df[1,9] <- model_n[which.min(aic_vals$BIC)] # best model name
  aic_df[1,10] <- aic_df[1,11] <- pvals[which.min(aic_vals$BIC)] # best model pvalue
  aic_df[1,12] <- aic_df[1,13] <- if (!is.null(best_mod$type_gam)) best_mod$type_gam else "" # best ac coeff category
  colnames(aic_df) <- c("best_mod_orig", "lin_val", "exp_val", "echo_val", "echo_lin_val", "echo_joint_val", "echo_lin_joint_val", "free_joint_val", "best_mod_after", "best_mod_pval_orig", "best_mod_pval_after", "best_mod_type_gam_orig", "best_mod_type_gam_after") # names

  # results list
  all_mod <- list("best_mod" = best_mod, "aic_df" = aic_df)

  return(all_mod)
}

# function to run the full mosaic pipeline
# inputs:
#  current_gene: index for gene that we're calculating parameters for
#  timen: time points vector
#  resol: resolution
#  num_reps: number of replicates
#  final_df_row: row for final mosaic data frame, unfilled
#  genes_rna: dataframe of rna gene expressions (first column is names, other columns are numeric)
#  genes_protein: dataframe of protein gene expressions (first column is names, other columns are numeric)
#  avg_genes_rna: matrix of averaged gene expressions, over replicates, for rna
#  avg_genes_pro: matrix of averaged gene expressions, over replicates, for protein
#  rem_unexpr_combine: a logical vector indicating whether the gene is unexpressed in rna OR protein
#  harm_cut: postive number indicating the cutoff for a gene to be considered harmonic
#  over_cut: postive number indicating the cutoff for a gene to be considered repressed/overexpressed
# outputs:
#  list containing the final filled mosaic data frame row, and the bic dataframes for rna and protein (before and after joint modeling results)
multi_pipeline <- function(current_gene, timen, resol, num_reps, final_df_row, genes_rna, genes_pro, avg_genes_rna, avg_genes_pro, rem_unexpr_combine, harm_cut = .03, over_cut = .15){
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

  # remove genes with less than 3 datapoints -- should be unexpressed
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
  exp_param_n <- c("Initial_Amplitude", "Growth_Rate", "Equilibrium_Value")
  lin_param_n <- c("Slope", "Equilibrium_Value")

  # put all the specific parameter names in a list
  param_map <- list(
    "echo_lin" = echo_lin_param_n,
    "echo_lin_joint" = echo_lin_joint_param_n,
    "echo" = echo_param_n,
    "echo_joint" = echo_joint_param_n,
    "exp" = exp_param_n,
    "lin" = lin_param_n
  )

  # so, we start by getting the best model for RNA and protein
  # rna:
  genes <- genes_rna
  avg_genes <- avg_genes_rna
  rna_mods <- choose_best_mod(current_gene, timen, resol, num_reps, genes_rna, avg_genes_rna, harm_cut, over_cut)
  aic_rna <- rna_mods$aic_df # save the bic (aic is a misnomer based on earlier tests)

  # protein
  genes <- genes_pro
  avg_genes <- avg_genes_pro
  pro_mods <- choose_best_mod(current_gene, timen, resol, num_reps, genes_pro, avg_genes_pro, harm_cut, over_cut)
  aic_pro <- pro_mods$aic_df # save the bic (aic is a misnomer based on earlier tests)

  # now, are both models oscillatory, or is one oscillatory and one is not?
  send_joint <- which((aic_rna$best_mod_orig == "echo" | aic_pro$best_mod_orig == "echo" |
                   aic_rna$best_mod_orig == "echo_lin" | aic_pro$best_mod_orig == "echo_lin") &
    !(aic_rna$best_mod_orig == "echo" & aic_pro$best_mod_orig == "echo") &
    !(aic_rna$best_mod_orig == "echo_lin" & aic_pro$best_mod_orig == "echo_lin") &
    !(aic_rna$best_mod_orig == "echo_lin" & aic_pro$best_mod_orig == "echo") &
    !(aic_rna$best_mod_orig == "echo" & aic_pro$best_mod_orig == "echo_lin"))

  if (length(send_joint) > 0){ # if we have a joint model
    # first, we need to get the joint model

    # get what the oscillatory model to jointly model on is supposed to be (echo or echo linear)
    if ("echo" %in% c(aic_rna$best_mod_orig[1], aic_pro$best_mod_orig[1])){
      # get one choice for starting points, using a completely joint model
      # echo completely joint model
      ej_mod <- calc_param_echo_joint(current_gene, timen, resol, num_reps, "orig", harm_cut, over_cut, genes_rna, genes_pro, avg_genes_rna, avg_genes_pro)
      res <- ej_mod$param #extract parameter estimates

      # final fit: echo joint model with slack for protein
      ejc_mod <- calc_param_echo_joint_constr(current_gene, timen, resol, num_reps, "orig", res, harm_cut, over_cut, genes_rna, genes_pro, avg_genes_rna, avg_genes_pro)
    } else {
      # get one choice for starting points, using a completely joint model
      # echo lin completely joint model
      ej_mod <- calc_param_echo_lin_joint(current_gene, timen, resol, num_reps, "orig", harm_cut, over_cut, genes_rna, genes_pro, avg_genes_rna, avg_genes_pro)
      res <- ej_mod$param #extract parameter estimates

      # final fit: echo linear joint model with slack for protein
      ejc_mod <- calc_param_echo_lin_joint_constr(current_gene, timen, resol, num_reps, "orig", res, harm_cut, over_cut, genes_rna, genes_pro, avg_genes_rna, avg_genes_pro)
    }

    # then, compare with best free joint model
    # get best free joint model for rna and protein
    rna_mod <- rna_mods$best_mod
    pro_mod <- pro_mods$best_mod

    # compare joint to free model
    best_mod <- reeval_best_mod(ejc_mod, rna_mod, pro_mod, current_gene, timen, num_reps, resol, genes_rna, genes_pro)

    if (best_mod$mod_name != "free"){ # joint model is best model
      joint_mod <- best_mod$joint_mod
      # we want to substitute the joint model into the bic dataframes, updating the best model
      aic_rna$best_mod_after[1] <- aic_pro$best_mod_after[1] <- best_mod$mod_name # name
      aic_rna[1,paste0(best_mod$mod_name,"_val")] <-
        aic_pro[1,paste0(best_mod$mod_name,"_val")] <- best_mod$aic_joint # bic
      aic_rna$best_mod_pval_after[1] <- joint_mod$rna_pval # rna pvalue
      aic_pro$best_mod_pval_after[1] <- joint_mod$pro_pval # protein pvalue

      aic_rna$free_joint_val[1] <- aic_pro$free_joint_val[1] <- best_mod$aic_free # free model name
      aic_rna$best_mod_type_gam_after[1] <- joint_mod$type_gam_rna # ac coefficient category, rna
      aic_pro$best_mod_type_gam_after[1] <- joint_mod$type_gam_pro # ac coefficient category, protein

      # update the final mosaic dataframes

      # add the designation and the model names
      final_df$Best_Model_RNA[1] <- final_df$Best_Model_Protein[1] <- mod_name_map[[best_mod$mod_name]]
      final_df$P_Value_RNA[1] <- aic_rna$best_mod_pval_after[1]
      final_df$P_Value_Protein[1] <- aic_pro$best_mod_pval_after[1]
      final_df$P_Value_Joint[1] <- joint_mod$joint_pval

      # add the parameters
      final_df[1, param_map[[joint_mod$mod_name]]] <- joint_mod$final

      # add the fits
      final_df[1,(33+(length(timen)*num_reps)+1):(33+(length(timen)*num_reps)+length(timen))] <-
        get_mod_fits(joint_mod, timen, resol)[1:length(timen)]
      final_df[1,(33+(length(timen)*num_reps*2)+length(timen)+1):ncol(final_df)] <-
        get_mod_fits(joint_mod, timen, resol)[(length(timen)+1):(length(timen)*2)]
    } else { # the free model was best
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
  } else { # there is no reason to joint model
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

  # results list
  fin <- list(
    "final_df" = final_df,
    "aic_rna" = aic_rna,
    "aic_pro" = aic_pro
  )

  return(fin)

}

# function to combine respectively fitted models into one "free" model
# inputs:
#  rna_mod: best model object for rna (see the calc_param functions)
#  pro_mod: best model object for protein (see the calc_param functions)
#  current_gene: index for gene that we're calculating parameters for
#  timen: time points vector
#  num_reps: number of replicates
#  resol: resolution
#  genes_rna: dataframe of rna gene expressions (first column is names, other columns are numeric)
#  genes_protein: dataframe of protein gene expressions (first column is names, other columns are numeric)
# outputs:
#  nls model object of combined free models
get_free_comb_mod <- function(rna_mod, pro_mod, current_gene, timen, num_reps, resol, genes_rna, genes_pro){
  # model names
  model_n <- c("echo_lin","echo","exp_lin","exp","lin","const")

  # get the concatenated time points and weights
  t_stack <- c(rep(timen, each=num_reps), (rep(timen, each = num_reps)))
  weights <- c(rep(calc_weights_boot(unlist(genes_rna[current_gene,-1]), num_reps), each = num_reps),
               rep(calc_weights_boot(unlist(genes_pro[current_gene,-1]), num_reps), each = num_reps))

  # time points, repeated for each replicate
  t_rep <- length(rep(timen, each=num_reps))
  # binary vectors to indicate rna and protein in concatenation
  s1 <- rep(c(1,0), c(t_rep,t_rep)) # rna
  s2 <- rep(c(0,1), c(t_rep,t_rep)) # protein

  # figure out which omics type is the assumed oscillatory one
  # rna oscillatory
  if (which(model_n == rna_mod$mod_name) <= which(model_n == pro_mod$mod_name)){
    # order: rna is concatenated first
    rna_ord <- 1
    pro_ord <- 2

    # order models
    first_mod <- rna_mod
    second_mod <- pro_mod

    # concatenate data accordingly
    y_stack <- c(unlist(genes_rna[current_gene,-1]), unlist(genes_pro[current_gene,-1]))

    # model names put together
    mod_name <- paste0(rna_mod$mod_name,"_",pro_mod$mod_name)
  } else { # protein oscillatory
    # order: protein is concatenated first
    rna_ord <- 2
    pro_ord <- 1

    # order models
    first_mod <- pro_mod
    second_mod <- rna_mod

    # concatenate data accordingly
    y_stack <- c(unlist(genes_pro[current_gene,-1]), unlist(genes_rna[current_gene,-1]))

    # model names put together
    mod_name <- paste0(pro_mod$mod_name,"_",rna_mod$mod_name)
  }
  # put together the data
  temp <- data.frame(y = y_stack, t =t_stack, s1 = s1, s2 = s2, resol = rep(resol,length(y_stack)))
  temp <- temp[!is.na(temp$y),] # remove any missing data points

  # get the coefficients
  coeff1 <- as.numeric(first_mod$param[-1]) # oscillatory
  coeff2 <- as.numeric(second_mod$param[-1]) # non-oscillatory

  # model choice - very long switch statement!
  # nls model object of the free model
  free.fit <- {
    switch(mod_name,
           "echo_lin_lin" = { # echo linear, linear concatenation
             nls_edit(
               formula = y ~ echo_lin_lin(a1,gam1,omega1,phi1,slo1,y_shift1,slo2,y_shift2,t, s1, s2),
               data=temp,
               start = list(gam1=coeff1[1], a1=coeff1[2], omega1 = coeff1[3], phi1 = coeff1[4], slo1 = coeff1[5], y_shift1=coeff1[6],
                            slo2 = coeff2[1], y_shift2 = coeff2[2]
               ),
               control = nls.control(maxiter = 0)
             )
           },
           "echo_lin_exp" = { # echo linear, exponential concatenation
             nls_edit(
               formula = y ~ echo_lin_exp(a1,gam1,omega1,phi1,slo1,y_shift1,a2,gamma2,y_shift2,t, s1, s2),
               data=temp,
               start = list(gam1=coeff1[1], a1=coeff1[2], omega1 = coeff1[3], phi1 = coeff1[4], slo1 = coeff1[5], y_shift1=coeff1[6],
                            a2= coeff2[1],gamma2 = coeff2[2], y_shift2 = coeff2[3]
               ),
               control = nls.control(maxiter = 0)
             )
           },
         "echo_lin" = { # echo, linear concatenation
           nls_edit(
             formula = y ~ echo_lin(a1,gam1,omega1,phi1,y_shift1,slo2,y_shift2,t, s1, s2),
             data=temp,
             start = list(gam1=coeff1[1], a1=coeff1[2], omega1 = coeff1[3], phi1 = coeff1[4], y_shift1=coeff1[5],
                          slo2 = coeff2[1], y_shift2 = coeff2[2]
             ),
             control = nls.control(maxiter = 0)
           )
         },
         "echo_exp" = { # echo, exponential concatenation
           nls_edit(
             formula = y ~ echo_exp(a1,gam1,omega1,phi1,y_shift1,a2,gamma2,y_shift2,t, s1, s2),
             data=temp,
             start = list(gam1=coeff1[1], a1=coeff1[2], omega1 = coeff1[3], phi1 = coeff1[4], y_shift1=coeff1[5],
                          a2= coeff2[1],gamma2 = coeff2[2], y_shift2 = coeff2[3]
             ),
             control = nls.control(maxiter = 0)
           )
         }
         )
  }

  return(free.fit)
}

# function to compare the fitted joint model to the fitted "free" model
# inputs:
#  joint model:  best joint model object (see the calc_param functions)
#  rna_mod: best model object for rna (see the calc_param functions)
#  pro_mod: best model object for protein (see the calc_param functions)
#  current_gene: index for gene that we're calculating parameters for
#  timen: time points vector
#  num_reps: number of replicates
#  resol: resolution
#  genes_rna: dataframe of rna gene expressions (first column is names, other columns are numeric)
#  genes_protein: dataframe of protein gene expressions (first column is names, other columns are numeric)
# outputs:
#  list containing which fitted model is better, the BICs of both, and the model (if joint is best)
reeval_best_mod <- function(joint_mod, rna_mod, pro_mod, current_gene, timen, num_reps, resol, genes_rna, genes_pro){
  # get the concatenated "free" model
  free.fit <- get_free_comb_mod(rna_mod, pro_mod, current_gene, timen, num_reps, resol, genes_rna, genes_pro)

  # compare free to joint models with BIC
  mod_ord <- c("free","joint")
  aic <- BIC(free.fit, joint_mod$model)

  # choose the model with the lowest BIC
  if (mod_ord[which.min(aic$BIC)] == "free"){
    # free is best
    return(
      list(
        "mod_name" = "free",
        "aic_free" = aic$BIC[1],
        "aic_joint" = aic$BIC[2]
      )
    )
  } else {
    # joint is best
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

