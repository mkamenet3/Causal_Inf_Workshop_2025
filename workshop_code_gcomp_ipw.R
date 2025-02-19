######################################################################################################################
# Author: Alex Keil
# Program: workshop_code.sas
# Date: 20180826
# Project: ISEE 2018 workshop: Causal inference foundations and applications in environmental health sciences
# Tasks:
# Data in: miners, coalplant
# Description: Carry out parametric g-formula analysis and strucgtural nested model analysis of data simulated to
#   emulate an occupational study of uranium miners looking at the effects of occupational exposures to radon
#   and all cause mortality, as well as a population study of air pollutants that arise from coal-fired power
#   plants
# Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
######################################################################################################################

# important packages
# install with: 
# install.packages(c('ggplot2', 'readr', 'dplyr', 'survival', 'tibble', 'boot', 'Hmisc', 'zoo'), repos = 'https://cran.mtu.edu')
library(ggplot2)
library(readr)
library(dplyr)
library(survival)
library(tibble)
library(boot)
library(Hmisc)
library(zoo)
# reading in data
# set working directory to the 2025_NCI_causal folder: e.g. setwd("~/Documents/2025_NCI_causal")
coalplant = read_csv("analyses/data/coalplant.csv", show_col_types = FALSE)
miners = read_csv("analyses/data/miners.csv", show_col_types = FALSE)
outpath = "analyses/output/"

################################################################################
#  Coal plant example - parametric g-formula with point exposure
################################################################################
# MDI (mental development index) appears to be roughly normally distributed
#  so we will model using a linear model

mdimod = glm(mdi ~ as*urbanicity + be*urbanicity + cd*urbanicity + as*be + as*cd + be*cd, data=coalplant)

# Say we know based on local emissions data and discussions with atmospheric
#  chemists that 91% of As, 96% of Be, and 45% of Cd ambient levels come
#  from a local coal-fired power plant


gformula_means <- function(ymodel, data, ...){
   # function to calculate mean MDI under a given intervention
   require('dplyr')
   # if only 'ymodel' and 'data' are supplied, then natural course is fit
  postinterventionX <- mutate(data,...)
  mean(predict(ymodel, newdata=postinterventionX))
}

#point estimates for mean MDI under each intervention
nc = gformula_means(mdimod, data=coalplant)
no_coal = gformula_means(mdimod, data=coalplant, 
                         cd=cd*(1-0.5), as=as*(1-0.96), be=be*(1-0.91))
no_cd = gformula_means(mdimod, data=coalplant, cd=0)
no_as = gformula_means(mdimod, data=coalplant, as=0)
no_be = gformula_means(mdimod, data=coalplant, be=0)
# expected population mean MDI
print(c(nc=nc, no_coal=no_coal, no_cd=no_cd, no_as=no_as, no_be=no_be), 4)

# bootstrap to get confidence intervals on effects (boot function also gives point estimates from original data)
gformula_meandiffs <- function(data, index){
  mcsample = data[index,]
  bootmod = glm(mdi ~ as*urbanicity + be*urbanicity + cd*urbanicity + as*be + as*cd + be*cd, data=mcsample)
  nc = gformula_means(bootmod, data=mcsample)
  no_coal = gformula_means(bootmod, data=mcsample, cd=cd*(1-0.5), as=as*(1-0.96), be=be*(1-0.91))
  no_cd = gformula_means(bootmod, data=mcsample, cd=0)
  no_as = gformula_means(bootmod, data=mcsample, as=0)
  no_be = gformula_means(bootmod, data=mcsample, be=0)
  c(nc_mdi=nc, shutdown=no_coal-nc, attrmdi_cd=no_cd-nc, attrmdi_as=no_as-nc, attrmdi_be=no_be-nc)
}

set.seed(12321)
bootsamples = boot(data=coalplant, statistic=gformula_meandiffs, R=500)

# effect estimates with confidence intervals
se = apply(bootsamples$t, 2, sd)
print(cbind(estimate=bootsamples$t0, lowerCI=bootsamples$t0-1.96*se, upperCI=bootsamples$t0+1.96*se), 3)

################################################################################
#  Uranium miners examples (longitudinal data)
################################################################################



################################################################################
#
#
# g-computation/g-formula
#
#
################################################################################

get_coefs <- function(data){
  # pooled linear model for log(exposure), while at work
  data$lograd = with(data, ifelse(atwork==1, log(rad), NA))
  mod_e = glm(lograd ~ outtime + cum_rad_lag1 + smoker, data=filter(data, atwork==1))
  # pooled logistic model for leaving work (at the end of the year), while at work
  mod_w = glm(leavework ~ outtime + I(outtime^2) + rad_lag1 + cum_rad_lag1 + smoker, data=filter(data, atwork==1 | leavework==1), family=binomial(link=logit))
  # pooled logistic model for death
  mod_d = glm(dead ~ outtime + atwork + cum_rad + I(cum_rad*cum_rad) + smoker + smoker*cum_rad, data=data, family=binomial(link=logit))
  # when troubleshooting code, it is helpful to have a time-only model (which should fit the natural course well)
  #mod_e = glm(lograd ~ outtime + I(outtime^2), data=filter(data, atwork==1))
  #mod_w = glm(leavework ~ outtime + I(outtime^2), data=filter(data, atwork==1 | leavework==1), family=binomial(link=logit))
  #mod_d = glm(dead ~ outtime + I(outtime^2) + I(outtime^3), data=data, family=binomial(link=logit))
  list(mod_e=mod_e, mod_w=mod_w, mod_d=mod_d)
}

gformula_risks <- function(intervention=NULL, pseudo_cohort_bl, 
                           endtime, mod_e, mod_w, mod_d, seed=NULL){
  # esimate intervention specific risks using a Monte Carlo algorithm
  # this function assumes ordering by ID, TIME
  require(dplyr)
  N <- dim(pseudo_cohort_bl)[1]
  ncols <- dim(pseudo_cohort_bl)[2]
  # create container to hold simulated results
  pseudo_cohort <- slice(pseudo_cohort_bl,rep(1:N, each=endtime))
  pseudo_cohort$gfid <- rep(seq(1, N, 1), each=endtime) # new id for each pseudo individual
  
  # initialize values
  pseudo_cohort$intime <- rep(seq(0, (endtime-1), 1), times=N)
  pseudo_cohort$outtime <- pseudo_cohort$intime+1
  pseudo_cohort$atwork <- 1
  pseudo_cohort$durwork <- 1
  pseudo_cohort$leavework <- 0
  pseudo_cohort$rad <- 0
  pseudo_cohort$rad_lag1 <- 0
  pseudo_cohort$cum_rad_lag1 <- 0
  pseudo_cohort$durwork_lag1 <- 0
  pseudo_cohort$dead <- 0
  pseudo_cohort$cum_hazard <- 0 # new cumulative hazard
  pseudo_cohort$cum_incidence <- 0 # new cumulative incidence
  # error term variance for exposure model
  #var_e = sum(mod_e$residuals^2)/mod_e$df.residual # default R approach (unbiased SD - better approach)
  var_e = mean(mod_e$residuals^2) # to match with SAS version (unbiased variance)
  set.seed(seed)
  # limits on log exposure (keep simulated values in range of observed)
  max_e = 7.321276
  
  # simulate time-varying values
  # R is slow with looping, so this simulation is vectorized as much as possible
  for(t in seq(1, endtime, 1)){
    # index that keeps track of time
    idx = which(pseudo_cohort$outtime == t)
    #update all values of lagged variables (easy to forget!)
    if(t > 1){
      pseudo_cohort[idx,'cum_rad_lag1'] = pseudo_cohort[idx-1,'cum_rad']
      pseudo_cohort[idx,'rad_lag1'] = pseudo_cohort[idx-1,'rad']
    }    
    ############ 
    # leaving work (time varying confounder)
    ############
    if(t==1) widx = 1:length(idx) # everyone is at work at first time point
    if(t > 1){
       widx = which(pseudo_cohort[idx-1,'atwork']==1)
       # index keeping track of which workers still work
       pseudo_cohort[idx[widx],'leavework'] <- rbinom(n=length(widx), size=1, 
                                                      prob=predict(mod_w, newdata = pseudo_cohort[idx[widx],], 
                                                                   type = 'response'))
       # if worker didn't leave, then assume stay at work
       pseudo_cohort[idx,'atwork'] <- (pseudo_cohort[idx,'leavework']==0 & pseudo_cohort[idx-1,'atwork']==1)
       # update at work index to account for workers who left
       widx = which(pseudo_cohort[idx,'atwork']==1)
       pseudo_cohort[idx,'durwork'] <- pseudo_cohort[idx-1,'durwork'] + pseudo_cohort[idx,'atwork']
       pseudo_cohort[idx,'durwork_lag1'] <- pseudo_cohort[idx-1,'durwork']
    }
    ############
    # exposure/interventions
    ############
    # exposure is assumed log-normally distributed (=predicted value plus draw from the residuals)
    meanlogr = predict(mod_e, newdata = pseudo_cohort[idx[widx],])
    logr = meanlogr + rnorm(n=length(widx), mean=0, sd=sqrt(var_e))
    pseudo_cohort[idx[widx],'rad'] <- exp(pmin(log(max_e), logr))
    # exposure under interventions
    if(typeof(intervention) != "NULL"){
      if(is.numeric(intervention)){
        # static, deterministic intervention, e.g. set exposure = 0
        pseudo_cohort[idx,'rad'] <- intervention
      } else{
      if(is.character(intervention)){
        # dynamic interventions are defined by character strings: e.g. 'rad>0.3' means we intervene to prevent radon exposure from being > 0.3
        #  this line creates an index of rows at the current time which are in violation of the intervention
        #  exposure will be redrawn for these observations until they are in compliance.
        #  This is equivalent to drawing from a truncated log-normal distribution
        viol_idx <- with(pseudo_cohort[idx,], which(eval(parse(text=intervention))))
        while(length(viol_idx)>0){
          # accept/rejection sampling to get draws from truncated log-normal
          # this is how we (I) assume an intervention would be implemented, but
          # we could imagine other ways (e.g. a hard cap on exposure)
          meanlogr = predict(mod_e, newdata = pseudo_cohort[idx[viol_idx],])
          lograd = meanlogr + rnorm(n=length(viol_idx), mean=0, sd=sqrt(var_e))
          pseudo_cohort[idx[viol_idx],'rad'] <- exp(lograd)
          # check whether any are still in violation of the distribution, if so, repeat loop
          viol_idx <- with(pseudo_cohort[idx,], which(eval(parse(text=intervention))))
        }}} # end dynamic intervention
    } # end all interventions on exposure
    if(t > 1){
      pseudo_cohort[idx,'cum_rad'] = pseudo_cohort[idx-1,'cum_rad'] + pseudo_cohort[idx,'rad']
    } else pseudo_cohort[idx,'cum_rad'] = pseudo_cohort[idx,'rad']
    ############
    # death (simulate discrete hazard as in Taubman et al 2009, rather than 1/0 as in Keil et al 2014)
    #  see note at bottom of program about late entry - in the case of late entry
    #  the typical approach is to simulate actual death (1/0 variable) and then
    #  estimate survival in the pseudo-population using a survival curve estimator
    #  that can account for late entry
    ############
    pseudo_cohort[idx,'dead'] = predict(mod_d, newdata = pseudo_cohort[idx,], type = 'response')
    if(t > 1){
      # product limit estimator of cumulative incidence (note this is applied on the individual basis)
      pseudo_cohort[idx,'cum_incidence'] = pseudo_cohort[idx-1,'cum_incidence'] + 
        (1-pseudo_cohort[idx-1,'cum_incidence'])*pseudo_cohort[idx,'dead']
    } else{
      pseudo_cohort[idx,'cum_incidence'] = pseudo_cohort[idx,'dead']
    }
    #alternative estimator of the cumulative incidence: Breslow estimator
    #if(t > 1){
    #  # kaplain-meier estimator of cumulative incidence (note this is applied on the individual basis)
    #  pseudo_cohort[idx,'cum_hazard'] = pseudo_cohort[idx-1,'cum_hazard'] + pseudo_cohort[idx,'dead']
    #} else{
    #  pseudo_cohort[idx,'cum_hazard'] = pseudo_cohort[idx,'dead']
    #}
    #  pseudo_cohort[idx,'cum_incidence'] = 1-exp(-pseudo_cohort[idx,'cum_hazard'])
    # could also use Breslow estimator, but need aalen-johansen estimator if there are competing risks
    # the 'cum_incidence' variable is an estimate of the individual risk. We then
    # average that over the population below to get the marginal risk
  } # end time loop
  pseudo_cohort
}# end function

# observed data
baseline_data <- filter(miners, intime==0)
N = dim(baseline_data)[1]


# large sample from observed baseline data
mc_iter = 20000
set.seed(12325)
mcdata <- slice(baseline_data, sample(1:N, size = mc_iter, replace=TRUE))

# example of getting the cumulative incidence curves
mods = get_coefs(data=miners)
nc = gformula_risks(intervention=NULL, pseudo_cohort_bl=mcdata, endtime=20, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=12325)
unex = gformula_risks(intervention=0, pseudo_cohort_bl=mcdata, endtime=20, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=12325)
exp3 = gformula_risks(intervention='rad>0.3', pseudo_cohort_bl=mcdata, endtime=20, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=12325)
exp1 = gformula_risks(intervention='rad>0.1', pseudo_cohort_bl=mcdata, endtime=20, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=12325)

# comparing natural course and observed, cumulative incidence
# we get population risk/cumulative incidence by the mean of the individual risks
obscifit <- survfit(Surv(intime, outtime, dead)~1, data=miners)
gfcinc <- with(nc, tibble(time=seq(0,20, 1), ci=c(0, tapply(cum_incidence, outtime, mean))))
gfciunex <- with(unex, tibble(time=seq(0,20, 1), ci=c(0, tapply(cum_incidence, outtime, mean))))
gfciexp3 <- with(exp3, tibble(time=seq(0,20, 1), ci=c(0, tapply(cum_incidence, outtime, mean))))
gfciexp1 <- with(exp1, tibble(time=seq(0,20, 1), ci=c(0, tapply(cum_incidence, outtime, mean))))
obsnc <- tibble(time=obscifit$time, ci=1-obscifit$surv)

# time ratio for comparison with SNM
# total expected years of life in natural course versus intervention to limit exposure to be below 0.1
# with no censoring this calculation is easy! The total expected number of years lived is just
# the survival function summed over all the person-time.
elifenc = with(nc, sum(1-cum_incidence))
elifeexp1 = with(exp1, sum(1-cum_incidence))
elifeexp1/elifenc # [1] 1.476379
# recall g-estimate was 1.36 (1.09, 1.67)


ggplot(aes(x=time, y=ci), data=obsnc) + 
  geom_line(aes(color="Observed")) +
  geom_line(aes(x=time, y=ci, color="Natural course"), data=gfcinc) +
  geom_line(aes(x=time, y=ci, color="Unexposed"), data=gfciunex) +
  geom_line(aes(x=time, y=ci, color="Limit: 0.3"), data=gfciexp3) +
  geom_line(aes(x=time, y=ci, color="Limit: 0.1"), data=gfciexp1) +
  scale_y_continuous(name="Cumulative incidence") +
  scale_x_continuous(name="Time", limits=c(0,20)) +
  scale_color_discrete(name="") +
  theme_classic() +
  theme(legend.position = c(0.99,0.01),legend.justification = c(1,0))+
  ggtitle(expression(paste("G-formula estimates of cumulative incidence")))
ggsave(paste0(outpath, 'r_gformulaCI.png'), width=5, height=4)
  


# now we turn the g-formula algorithm to get risk difference into a function in order
# to get bootstrap variance
gformula_effectestimates <- function(data=baseline_data, index=NULL, fdata=miners, endtime=10, seed=NULL){
  # for bootstrap samples, when needed: first take a random sample of the 
  #  study IDs, with replacement
  # 0) preprocessing: take bootstrap sample if bootstrapping, otherwise use the observed data
    if(is.null(index)){
      # if index is NULL, then just use original data and do nothing here
    } else{
      # bootstrap sample of the baseline data using the vector 'index'
      data = slice(data, index)
      # resample from the full data according to the id variable in the baseline data
      idindex = tibble(reps = table(data$id), id = as.numeric(names(table(data$id))))
      fdata = merge(idindex, fdata, by='id')
      fdata = slice(fdata, rep(1:dim(fdata)[1], times=fdata$reps))
    }
  # 1) fit models to full data/bootstrap sample of the full data
    mods = get_coefs(data=fdata)
  # 2) take large MC sample from baseline data
    mc_iter = 20000
    set.seed(seed)
    mcdata <- slice(data, sample(1:N, size = mc_iter, replace=TRUE))
  # 3) simulate probability distributions in large MC sample using model fits from step 1
    nc = gformula_risks(intervention=NULL,        pseudo_cohort_bl=mcdata, endtime=endtime, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=seed)
    unex = gformula_risks(intervention=0,         pseudo_cohort_bl=mcdata, endtime=endtime, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=seed)
    exp3 = gformula_risks(intervention='rad>0.3', pseudo_cohort_bl=mcdata, endtime=endtime, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=seed)
    exp1 = gformula_risks(intervention='rad>0.1', pseudo_cohort_bl=mcdata, endtime=endtime, 
                    mod_e=mods$mod_e, mod_w=mods$mod_w, mod_d=mods$mod_d, seed=seed)
    # note: I change the seed for each of these to emulate SAS approach - this is not required
  # 4) summarize data
    ci_nc = as.numeric(with(nc, tapply(cum_incidence, outtime, mean)))
    ci_unex = as.numeric(with(unex, tapply(cum_incidence, outtime, mean)))
    ci_exp3 = as.numeric(with(exp3, tapply(cum_incidence, outtime, mean)))
    ci_exp1 = as.numeric(with(exp1, tapply(cum_incidence, outtime, mean)))
   # return a vector of risk differences and log-risk ratios at time t=endtime
   # could also calculate for all times
     c(
       rd_exp3_10=c(ci_exp3-ci_nc)[endtime],
       rd_exp1_10=c(ci_exp1-ci_nc)[endtime],
       rd_unex_10=c(ci_unex-ci_nc)[endtime],
       lrr_exp3_10=log(c(ci_exp3/ci_nc)[endtime]),
       lrr_exp1_10=log(c(ci_exp1/ci_nc)[endtime]),
       lrr_unex_10=log(c(ci_unex/ci_nc)[endtime])
     )
}

# point estimate for risk differences at t=10 (also given by the 'boot' function below)
system.time(gf_est <- gformula_effectestimates(data=baseline_data, fdata=miners, endtime=10, seed=NULL))
# took about 3 seconds on my machine
# risk difference, log RR estimates (will vary slightly by seed value due to simulation approach)
gf_est

# getting confidence intervals: non-parametric bootstrap
# note that to get valid SE, we should use 200+ bootstrap samples
# for percentile based confidence intervals we should use 1000+
# for an approach that takes 3 seconds to run, 200 bootstrap samples = 10 minutes
#nbootsamples = 200
#system.time(boot_samples <- boot(data = baseline_data, fdata=miners, statistic = gformula_effectestimates, R = nbootsamples, endtime=10))
#   user  system elapsed 
#674.059 227.130 901.453 


nbootsamples = 10
system.time(boot_samples <- boot(data = baseline_data, fdata=miners, statistic = gformula_effectestimates, R = nbootsamples, endtime=10))
boot_samples$t0

# point estimate and confidence intervals
est = boot_samples$t0
lci =  boot_samples$t0 - 1.96*apply(boot_samples$t, 2, sd)
uci =  boot_samples$t0 + 1.96*apply(boot_samples$t, 2, sd)
cat(paste0(names(gf_est), ":", round(est, 3), " (", round(lci, 3), ", ", round(uci, 3), ")\n"))



# Further considerations:
# Administrative censoring: generally, we would stop follow-up of the pseudo-cohort
#   at the time of administrative censoring, but if we assume we have correct models
#   then follow-up could be, in principle, extended.
# Censoring by loss-to-follow-up: under assumption of correct model specification
#   the approach to censoring in g-computation is to simulate individuals covariate
#   and event history after the time of censorin
# Competing risks: The discrete, cause specific hazard from g-computation will be
#   consistent under correct model specification, even with competing risks. 
#   We can also "account" for competing risks by using an estimate of the cumulative
#   incidence that allows for competing events, like the Aalen-Johansen estimator.
#   If we use the Kaplan-Meier/product limit estimator with competing risks, we estimate the 
#   "conditional" risk rather than the unconditional risk. The unconditional risk
#   is typically what we want.
# Late entry: notice that the pseudo-cohort has all individuals enter at the same
#   time. This time scale is "time since hire" or "time on study" if this study
#   included individuals hired prior to the start of the study. We could late
#   enter people for a different time scale (e.g. age), but would have to calculate
#   risk differently - we would have to simulate values for death and then use an
#   estimator of the cumulative incidence that allows for late entry, such as
#   the Kaplan meier estimator
# Last: Functional form of all of the models matter. We could fit, for example,
#   machine learning or data adaptive approaches





################################################################################
#
#
# Inverse probability weighting
#
#
################################################################################



ipw_risks = function(data=miners){
  
  limits=c(Inf, 0.3, 0.1)
  #1b) create copy of the dataset for each limit
  clones <- do.call(rbind,
                    lapply(1:length(limits), function(x) {  
                      df = data                  
                      df$limit = as.numeric(limits[x])
                      df$cloneid = paste0(x, "_", df$id)
                      df                    
                    }
                    ))
  #2) artificially censor data (keeping first observation that is censored)
  cens_data <- clones %>%
    group_by(cloneid) %>%
    mutate(
      cens = pmin(1,cumsum(rad > limit)),
      drop = pmin(2, cumsum(cens))
    ) %>%
    group_by() %>%
    filter(
      drop < 2
    )
  
  # 3) create 1/0 weights to use for confounding/censoring during follow-up (alternative is to create new datasets)
  cens_data <- group_by(cens_data, cloneid) %>%
    mutate(
      one = 1,
      nobs = cumsum(one),
      fobs____ = as.numeric(nobs == 1)
    ) %>%
    mutate(
      conf_weight = as.numeric(nobs  == 1), # this duplicates fobs____ but is kept for clarity
      fu_weight = as.numeric((nobs  > 1) & (atwork == 1)),    
      dconf = 0,
      dcens = 0,
      intervened = 0,
    ) %>%
    select(-c(one, nobs)) %>%
    ungroup()
  # check: which(tapply(data$x, data$id, max) < limit & tapply(data$atwork, data$id, sum)>3)[2] # index 427
  # check: names(tapply(data$x, data$id, max))[427]
  # check: print(select(cens_data, c(id, x, limit, conf_weight, fu_weight, cens, fobs____)) %>% filter(id==427), n=50)
  # check: print(select(filter(data, id==427), c(id, x, atwork)))
  
  
  #4) fit censoring models (can include pre-baseline exposure in practice)
  
  # fit models if there is more than a small amount of censoring (seems to work either way, but this avoids convergence problems)
  for (l in limits){
    limidx = which(cens_data$limit == l)
    tempdat = cens_data[limidx,]
    tempdat$mxl = tempdat$cum_rad_lag1 / (tempdat$durwork_lag1 + as.numeric(tempdat$outtime==1))
    outkn0 = attr(rcspline.eval(filter(tempdat, conf_weight==1)$outtime, nk = 4), "knots")
    outkn = attr(rcspline.eval(filter(tempdat, fu_weight==1)$outtime, nk = 4), "knots")
    
    #yearkn0 = attr(rcspline.eval(filter(tempdat, conf_weight==1)$year, nk = 4), "knots")
    #agekn = attr(rcspline.eval(filter(tempdat, fu_weight==1)$age, nk = 4), "knots")
    #yearkn = attr(rcspline.eval(filter(tempdat, fu_weight==1)$year, nk = 4), "knots")
    #
    if(sum(tempdat[tempdat$conf_weight==1,"cens"]) > 10){
      summary(confnmod <- glm(cens ~ 1 , data = tempdat, weight=conf_weight, family=binomial()))
      summary(confdmod <- glm(cens ~ smoker, data = tempdat, weight=conf_weight, family=binomial()))
      # need to use base R operations here
      # only get non-zero predictions if "conf_weight" == 1 (eligible to be "censored" at baseline)
      cens_data[limidx,"nconf"] = tempdat$conf_weight*as.numeric(predict(confnmod, type="response"))
      cens_data[limidx,"dconf"] = tempdat$conf_weight*as.numeric(predict(confdmod, type="response"))
    } else{
      cens_data[limidx,"nconf"] = 0
      cens_data[limidx,"dconf"] = 0
    }
    if(sum(tempdat[tempdat$fu_weight==1,"cens"]) > 1){
      summary(censnmod <- glm(cens ~ outtime + rcspline.eval(outtime, knots=outkn), data = tempdat, weight=fu_weight, family=binomial()))
      summary(censdmod <- glm(cens ~  outtime + rcspline.eval(outtime, knots=outkn) + rad_lag1 + cum_rad_lag1 + smoker, data = tempdat, weight=fu_weight, family=binomial()))
      # only get non-zero predictions if "fu_weight" == 1 (eligible to be censored during follow-up)
      cens_data[limidx,"ncens"] = tempdat$fu_weight*as.numeric(predict(censnmod, type="response")) 
      cens_data[limidx,"dcens"] = tempdat$fu_weight*as.numeric(predict(censdmod, type="response")) 
    } else{
      cens_data[limidx,"ncens"] = 0
      cens_data[limidx,"dcens"] = 0
    }
  }
  print(str(cens_data))
  browser()
  # create final weights
  combined_wtd_data <- cens_data %>% 
    mutate(
      wtcontr = case_when(
        # weight: (stabilized) inverse probability of NOT being censored
        # note regarding unstabilized weights: Cain et al 2010 show that weight stabilization is not guaranteed to reduce variance
        ((fobs____ == 1) & (atwork==1)) ~ (1-cens)*(1-nconf)/(1-dconf),
        ((fobs____ == 0) & (atwork==1)) ~ (1-cens)*(1-ncens)/(1-dcens),
        .default=1
      )
    ) %>%
    #select(id,cloneid,intime,outtime,dead,rad,wtcontr,rad,cum_rad,atwork,cens,limit, smoker) %>% # optional step here to reduce comp. burden
    group_by(cloneid) %>%
    mutate(
      ipw = cumprod(wtcontr)
    ) %>%
    group_by() %>%
    mutate(
      # weight truncation at 10.0 following Cain et al 2010 (should be a user option)
      ipw_trunc = pmin(ipw, 10.0)
    )
  
  # check: summary(select(combined_wtd_data, c(ipw, ipw_trunc)))
  # check: print(select(combined_wtd_data, c(id, age, atwork, x, cens, wtcontr, ipw, ipw_trunc)) %>% filter(id==427), n=10)
  
  # check mean weights by intervention
  # note mean is taken across all possible observations, even those with weights = 0
 # N = nrow(data)
 # Nid = length(unique(data$id))
 # wtdx = data.frame(
 #   limit = limits,
 #   mean_cumx = do.call(c,lapply(limits, function(x) mean(combined_wtd_data[combined_wtd_data$limit == x  & combined_wtd_data$cens == 0,]$cumx))),
 #   n_conf = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$fobs____==1,]$cens))),
 #   n_cens = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$fobs____==0,]$cens))),
 #   n_d1w = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$ipw*combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$d1))),
 #   n_d2w = do.call(c,lapply(limits, function(x) sum(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$ipw*combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$d2))),
 #   mean_ipw = do.call(c,lapply(limits, function(x) mean(combined_wtd_data[combined_wtd_data$limit == x & combined_wtd_data$cens == 0,]$ipw)))
 # )
  

  ##########################################################################################
  # marginal structural (IP weighted) Cox model: effects of limit vs. ref
  # robust variance is done by ID level, but observations are done by CLONEID level
  # the robust variance gives valid standard errors for the hazard ratio
  ##########################################################################################
  combined_wtd_data$limitf = factor(combined_wtd_data$limit, levels=c("Inf", "0.3", "0.1"), labels=c("Natural course", "Limit: 0.3", "Limit: 0.1"))
  msm = coxph(Surv(intime, outtime, dead)~limitf, data=filter(combined_wtd_data, ipw>0), 
               id=cloneid, weight=ipw, cluster=id, x=FALSE, y=FALSE)
  ##########################################################################################
  # IP weighted survival curve: risk under each limit
  ##########################################################################################
  # Survival curves for each limit
 browser()
   tt0 = survfit(Surv(intime, outtime, dead)~limitf, data=combined_wtd_data, weight=ipw, id=cloneid)
  print(tt0)
  # return a vector of risk differences and log-risk ratios at time t=endtime
  
  results = data.frame(
    time=tt0$time, 
    risk = 1-tt0$surv,
    int = do.call(c,sapply(1:length(tt0$strata), function(x) rep(names(tt0$strata)[x], tt0$strata[x])))
  ) %>%
    mutate(
      ci0 = case_when(
        int == "limitf=Natural course" ~ risk,
        .default = as.numeric(NA)
      ),
      cip1 = case_when(
        int == "limitf=Limit: 0.1" ~ risk,
        .default = as.numeric(NA)
      ),
      cip3 = case_when(
        int == "limitf=Limit: 0.3" ~ risk,
        .default = as.numeric(NA)
      )
    ) %>%
    select(time, ci0, cip1, cip3) %>%
    rbind(data.frame(time=0, ci0=0, cip1=0, cip3=0)) %>%
    arrange(time) %>%
    na.locf() %>%
    mutate(
      rd_exp1_10 = cip1-ci0,
      rd_exp3_10 = cip3-ci0,
      lrr_exp1_10 = log(cip1/ci0),
      lrr_exp3_10 = log(cip3/ci0),
    ) %>%
    select(
      time, ci0, rd_exp1_10, rd_exp3_10, lrr_exp1_10, lrr_exp3_10
    )
    
  results
}


ipw_effectestimates <- function(data=baseline_data, index=NULL, fdata=miners, endtime=10, seed=NULL){
  # for bootstrap samples, when needed: first take a random sample of the 
  #  study IDs, with replacement
  # 0) preprocessing: take bootstrap sample if bootstrapping, otherwise use the observed data
  if(is.null(index)){
    # if index is NULL, then just use original data and do nothing here
  } else{
    # bootstrap sample of the baseline data using the vector 'index'
    data = slice(data, index)
    # resample from the full data according to the id variable in the baseline data
    idindex = tibble(reps = table(data$id), id = as.numeric(names(table(data$id))))
    #fdata = merge(idindex, fdata, by='id', sort=FALSE)
    #fdata = right_join(fdata, idindex, by=c('id'))
    fdata = left_join(idindex, fdata, by=c('id'))
    #sortindex = rep(1:dim(fdata)[1], times=fdata$reps)
    
    idxlist = sapply(1:dim(fdata)[1], function(x) rep(x, fdata$reps[x]))
    maxcopies = max(sapply(idxlist, length))
    sortindex_full = do.call(c,lapply(1:maxcopies, function(z) sapply(idxlist, function(x) x[z])))
    sortindex = sortindex_full[!is.na(sortindex_full)]
    fdata = slice(fdata, sortindex)
    fdata$oldid = fdata$id
    fdata$id = cumsum(fdata$intime==0) # need to keep track of new ID or bootstrap observations will get combined
  }
  ipw_boot = ipw_risks(data=fdata)
  unlist(tail(ipw_boot[ipw_boot$time<=endtime,], n=1))
}
  
# the entire risk curve
system.time(ipw_est <- ipw_risks(data=miners))



ggplot(aes(x=time), data=ipw_est) + 
  geom_line(aes(x=time, y=ci0, color="Natural course"), data=ipw_est) +
  geom_line(aes(x=time, y=ci0+rd_exp3_10, color="Limit: 0.3"), data=ipw_est) +
  geom_line(aes(x=time, y=ci0+rd_exp1_10, color="Limit: 0.1"), data=ipw_est) +
  scale_y_continuous(name="Cumulative incidence") +
  scale_x_continuous(name="Time", limits=c(0,20)) +
  scale_color_discrete(name="") +
  theme_classic() +
  theme(legend.position = c(0.99,0.01),legend.justification = c(1,0))+
  ggtitle(expression(paste("IPW estimates of cumulative incidence")))
ggsave(paste0(outpath, 'r_ipwCI.png'), width=5, height=4)



# getting confidence intervals: non-parametric bootstrap
# note that to get valid SE, we should use 200+ bootstrap samples
# for percentile based confidence intervals we should use 1000+
# for an approach that takes 3 seconds to run, 200 bootstrap samples = 10 minutes
#nbootsamples = 200
#system.time(boot_samples_ipw <- boot(miners, statistic = ipw_effectestimates, R = nbootsamples, endtime=10))
#    user  system elapsed 
# 120.876   5.029 126.876 


nbootsamples = 10
system.time(boot_samples_ipw <- boot(miners, statistic = ipw_effectestimates, R = nbootsamples, endtime=10))
boot_samples_ipw$t0

# point estimate and confidence intervals
est_ipw = boot_samples_ipw$t0
lci_ipw =  boot_samples_ipw$t0 - 1.96*apply(boot_samples_ipw$t, 2, sd)
uci_ipw =  boot_samples_ipw$t0 + 1.96*apply(boot_samples_ipw$t, 2, sd)
cat(paste0(names(ipw_est), ":", round(est_ipw, 3), " (", round(lci_ipw, 3), ", ", round(uci_ipw, 3), ")\n"))
