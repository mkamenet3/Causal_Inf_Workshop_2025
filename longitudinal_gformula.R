knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE,cache=TRUE)

#Load in libraries
library(ggplot2) #for plotting
library(dplyr) #for data cleaning
library(survival) #for survival analysis
library(boot) #for bootstrapping



miners <- read.csv("data/miners.csv") %>%
  mutate(lograd = ifelse(atwork==1, log(rad), NA))
str(miners)


summary(miners[, 2:13])



baseline_data <- filter(miners, intime==0) 
N = nrow(baseline_data)


# large sample from observed baseline data
mc_iter = 20000
set.seed(12325)
mcdata <- slice(baseline_data, sample(1:N, size = mc_iter, replace=TRUE))


#Model 1: pooled linear model for log(exposure), while at work
mod_e = lm(lograd ~ outtime + cum_rad_lag1 + smoker, data=filter(miners, atwork==1))


# pooled logistic model for leaving work (at the end of the year), while at work
mod_w = glm(leavework ~ outtime + I(outtime^2) + rad_lag1 + cum_rad_lag1 + smoker, data=filter(miners, atwork==1 | leavework==1), family=binomial())
 


mod_d = glm(dead ~ outtime + atwork + cum_rad + I(cum_rad*cum_rad) + smoker + smoker*cum_rad, data=miners, family=binomial(link=logit))
  
################################################################################################################################
N <- nrow(mcdata)
ncols <- ncol(mcdata)
pseudo_cohort <- slice(mcdata,rep(1:N, each=20))
#create a new identifier
pseudo_cohort$gfid <- rep(seq(1, N, 1), each=20) # new id for each pseudo individual
# initialize values
pseudo_cohort$intime <- rep(seq(0, (20-1), 1), times=N)
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

# set error variance for exposure model
var_e = sum(mod_e$residuals^2)/mod_e$df.residual 


# limits on log exposure (keep simulated values in range of observed)
max_e = 7.321276
################################################################################################################################




#for the natural course, there is no intervention so we set it to NULL
intervention <- NULL

head(pseudo_cohort)#check out the first 6 observations of the pseudo-cohort and confirm that initial values were initialized correctly

################################################################################################################################
set.seed(12325)
for(t in 1:20){
  # index that keeps track of time
  idx = which(pseudo_cohort$outtime == t)
  #update all values of lagged variables (easy to forget!)
  if(t > 1){
    #lagging the data because of disease latency
    pseudo_cohort[idx,'cum_rad_lag1'] = pseudo_cohort[idx-1,'cum_rad']
    pseudo_cohort[idx,'rad_lag1'] = pseudo_cohort[idx-1,'rad']
  }    
  ############ 
  # leaving work (time varying confounder)
  ############
  if(t==1){
    widx = 1:length(idx) # everyone is at work at first time point
  } 
  if(t > 1){
    #simulation below is forward in time, did someone leave work during this time point?
    widx = which(pseudo_cohort[idx-1,'atwork']==1) #work index
    # index keeping track of which workers still work
    ## simulate here from known probability distribution
    pseudo_cohort[idx[widx],'leavework'] <- rbinom(n=length(widx), size=1, 
                                                   prob=predict(mod_w, newdata = pseudo_cohort[idx[widx],], 
                                                                type = 'response'))
    #for everyone else still alive at time t, do they leave work at this time?
    # if worker didn't leave, then assume stay at work
    pseudo_cohort[idx,'atwork'] <- (pseudo_cohort[idx,'leavework']==0 & pseudo_cohort[idx-1,'atwork']==1)
    # update at work index to account for workers who left
    widx = which(pseudo_cohort[idx,'atwork']==1)
    #durwork: keeping track of work duration
    pseudo_cohort[idx,'durwork'] <- pseudo_cohort[idx-1,'durwork'] + pseudo_cohort[idx,'atwork']
    pseudo_cohort[idx,'durwork_lag1'] <- pseudo_cohort[idx-1,'durwork']
  }
  ############
  # exposure/interventions
  ############
  # exposure is assumed log-normally distributed (=predicted value plus draw from the residuals)
  # these are predicted values from the upstream process
  meanlogr = predict(mod_e, newdata = pseudo_cohort[idx[widx],]) 
  logr = meanlogr + rnorm(n=length(widx), mean=0, sd=sqrt(var_e))
  pseudo_cohort[idx[widx],'rad'] <- exp(pmin(log(max_e), logr))
  # exposure under different interventions/regimes
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
        #  
        #  Below is rejection sampling; if value is above the regime threshold,
        #  then reject it and sample again. This is more or less equivalent to sampling
        #  from a truncated distribution at regime limit
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
  #below is final prediction of death or not (if not, then keep calculating
  #cumulative incidence)
  #here, we are making mortality predictions from the updated data. All people
  #alive at time t and making predictions based on the simulated data
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
  #  # Kaplan-meier estimator of cumulative incidence (note this is applied on the individual basis)
  #  pseudo_cohort[idx,'cum_hazard'] = pseudo_cohort[idx-1,'cum_hazard'] + pseudo_cohort[idx,'dead']
  #} else{
  #  pseudo_cohort[idx,'cum_hazard'] = pseudo_cohort[idx,'dead']
  #}
  #  pseudo_cohort[idx,'cum_incidence'] = 1-exp(-pseudo_cohort[idx,'cum_hazard'])
  # could also use Breslow estimator, but need aalen-johansen estimator if there are competing risks
  # the 'cum_incidence' variable is an estimate of the individual risk. We then
  # average that over the population below to get the marginal risk
} # end time loop
head(pseudo_cohort)#check out first 6 observations with new predicted values
################################################################################################################################

cuminc_nc0 <- with(pseudo_cohort, tibble(time=seq(0,20, 1), ci=c(0, tapply(cum_incidence, outtime, mean)))) 
knitr::kable(cuminc_nc0, caption = "Cumulative incidence: unexposed regime (hard code)")


################################################################################################################################
gformula_risks <- function(intervention=NULL, pseudo_cohort_bl, 
                           endtime, mod_e, mod_w, mod_d, seed=NULL){
  N <- nrow(pseudo_cohort_bl)
  ncols <- ncol(pseudo_cohort_bl)
  pseudo_cohort <- slice(pseudo_cohort_bl,rep(1:N, each=endtime))
  pseudo_cohort$gfid <- rep(seq(1, N, 1), each=endtime) 
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
  pseudo_cohort$cum_hazard <- 0 
  pseudo_cohort$cum_incidence <- 0 
  var_e = mean(mod_e$residuals^2) 
  set.seed(seed)
  max_e = 7.321276
  for(t in seq(1, endtime, 1)){
    idx = which(pseudo_cohort$outtime == t)
    if(t > 1){
      pseudo_cohort[idx,'cum_rad_lag1'] = pseudo_cohort[idx-1,'cum_rad']
      pseudo_cohort[idx,'rad_lag1'] = pseudo_cohort[idx-1,'rad']
    }    
    if(t==1) widx = 1:length(idx) 
    if(t > 1){
      widx = which(pseudo_cohort[idx-1,'atwork']==1)
      pseudo_cohort[idx[widx],'leavework'] <- rbinom(n=length(widx), size=1, 
                                                     prob=predict(mod_w, newdata = pseudo_cohort[idx[widx],], 
                                                                  type = 'response'))
      pseudo_cohort[idx,'atwork'] <- (pseudo_cohort[idx,'leavework']==0 & pseudo_cohort[idx-1,'atwork']==1)
      widx = which(pseudo_cohort[idx,'atwork']==1)
      pseudo_cohort[idx,'durwork'] <- pseudo_cohort[idx-1,'durwork'] + pseudo_cohort[idx,'atwork']
      pseudo_cohort[idx,'durwork_lag1'] <- pseudo_cohort[idx-1,'durwork']
    }
    meanlogr = predict(mod_e, newdata = pseudo_cohort[idx[widx],])
    logr = meanlogr + rnorm(n=length(widx), mean=0, sd=sqrt(var_e))
    pseudo_cohort[idx[widx],'rad'] <- exp(pmin(log(max_e), logr))
    if(typeof(intervention) != "NULL"){
      if(is.numeric(intervention)){
        pseudo_cohort[idx,'rad'] <- intervention
      } else{
        if(is.character(intervention)){
          viol_idx <- with(pseudo_cohort[idx,], which(eval(parse(text=intervention))))
          while(length(viol_idx)>0){
            meanlogr = predict(mod_e, newdata = pseudo_cohort[idx[viol_idx],])
            lograd = meanlogr + rnorm(n=length(viol_idx), mean=0, sd=sqrt(var_e))
            pseudo_cohort[idx[viol_idx],'rad'] <- exp(lograd)
            viol_idx <- with(pseudo_cohort[idx,], which(eval(parse(text=intervention))))
          }}}
    } 
    if(t > 1){
      pseudo_cohort[idx,'cum_rad'] = pseudo_cohort[idx-1,'cum_rad'] + pseudo_cohort[idx,'rad']
    } else pseudo_cohort[idx,'cum_rad'] = pseudo_cohort[idx,'rad']
    pseudo_cohort[idx,'dead'] = predict(mod_d, newdata = pseudo_cohort[idx,], type = 'response')
    if(t > 1){
      pseudo_cohort[idx,'cum_incidence'] = pseudo_cohort[idx-1,'cum_incidence'] + 
        (1-pseudo_cohort[idx-1,'cum_incidence'])*pseudo_cohort[idx,'dead']
      #individual-level risks calculated using Kaplan-Meier formula and then summed and averaged
    } else{
      pseudo_cohort[idx,'cum_incidence'] = pseudo_cohort[idx,'dead']
    }
  } 
  pseudo_cohort
}

################################################################################################################################


nc = gformula_risks(intervention=NULL, pseudo_cohort_bl=mcdata, endtime=20, 
                    mod_e=mod_e, mod_w=mod_w, mod_d=mod_d, seed=12325)
gf_cuminc_nc <- with(nc, tibble(time=seq(0,20, 1), ci=c(0, tapply(cum_incidence, outtime, mean)))) 
knitr::kable(gf_cuminc_nc, caption = "Cumulative incidence: natural course regime")





################################################################################################################################
unex = gformula_risks(intervention=0, 
                      pseudo_cohort_bl=mcdata, 
                      endtime=20, 
                      mod_e=mod_e, mod_w=mod_w, mod_d=mod_d, 
                      seed=12325)
gf_cuminc_unex <- with(unex, tibble(time=seq(0,20, 1), ci=c(0, tapply(cum_incidence, outtime, mean)))) 
knitr::kable(gf_cuminc_unex, caption = "Cumulative incidence: unexposed regime")

exp3 = gformula_risks(intervention='rad>0.3', 
                      pseudo_cohort_bl=mcdata,
                      endtime=20, 
                      mod_e=mod_e, mod_w=mod_w, mod_d=mod_d, 
                      seed=12325)
gf_cuminc_exp3 <- with(exp3, tibble(time=seq(0,20, 1), ci=c(0, tapply(cum_incidence, outtime, mean)))) 
knitr::kable(gf_cuminc_exp3, caption = "Cumulative incidence: radon limit = 0.3 exposure regime")

exp1 = gformula_risks(intervention='rad>0.1', 
                      pseudo_cohort_bl=mcdata, 
                      endtime=20, 
                      mod_e=mod_e, mod_w=mod_w, mod_d=mod_d, 
                      seed=12325)
gf_cuminc_exp1 <- with(exp1, tibble(time=seq(0,20, 1), ci=c(0, tapply(cum_incidence, outtime, mean)))) 
knitr::kable(gf_cuminc_exp1, caption = "Cumulative incidence: radon limit = 0.1 exposure regime")


ggplot() +
  geom_line(aes(x=time, y=ci, color="Natural course"), data=gf_cuminc_nc) +
  geom_line(aes(x=time, y=ci, color="Unexposed"), data=gf_cuminc_unex ) +
  geom_line(aes(x=time, y=ci, color="Limit: 0.3"), data=gf_cuminc_exp3 ) +
  geom_line(aes(x=time, y=ci, color="Limit: 0.1"), data=gf_cuminc_exp1 ) +
  scale_y_continuous(name="Cumulative incidence") +
  scale_x_continuous(name="Time", limits=c(0,20)) +
  scale_color_discrete(name="") +
  theme_classic() +
  theme(legend.position = c(0.99,0.01),legend.justification = c(1,0))+
  ggtitle(expression(paste("G-formula estimates of cumulative incidence")))



obscifit <- survfit(Surv(intime, outtime, dead)~1, data=miners)
obsnc <- tibble(time=obscifit$time, ci=1-obscifit$surv)



ggplot() +
  geom_line(aes(x=time, y=ci, color="Natural course"), data=gf_cuminc_nc) +
  geom_line(aes(x=time, y=ci, color="Unexposed"), data=gf_cuminc_unex ) +
  geom_line(aes(x=time, y=ci, color="Limit: 0.3"), data=gf_cuminc_exp3 ) +
  geom_line(aes(x=time, y=ci, color="Limit: 0.1"), data=gf_cuminc_exp1 ) + 
  geom_line(aes(x=time, y=ci, color="Observed"), data=obsnc) +
  scale_y_continuous(name="Cumulative incidence") +
  scale_x_continuous(name="Time", limits=c(0,20)) +
  scale_color_discrete(name="") +
  theme_classic() +
  theme(legend.position = c(0.99,0.01),legend.justification = c(1,0))+
  ggtitle(expression(paste("G-formula estimates of cumulative incidence")))

# now we turn the g-formula algorithm to get risk difference into a function in order to get bootstrap variance
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
  mod_e = glm(lograd ~ outtime + cum_rad_lag1 + smoker, data=filter(data, atwork==1))
  # pooled logistic model for leaving work (at the end of the year), while at work
  mod_w = glm(leavework ~ outtime + I(outtime^2) + rad_lag1 + cum_rad_lag1 + smoker, 
              data=filter(fdata, atwork==1 | leavework==1), family=binomial(link=logit))
  # pooled logistic model for death
  mod_d = glm(dead ~ outtime + atwork + cum_rad + I(cum_rad*cum_rad) + smoker + smoker*cum_rad, 
              data=fdata, family=binomial(link=logit))
  
  # 2) take large MC sample from baseline data
  mc_iter = 20000
  set.seed(seed)
  mcdata <- slice(data, sample(1:N, size = mc_iter, replace=TRUE))
  # 3) simulate probability distributions in large MC sample using model fits from step 1
  nc = gformula_risks(intervention=NULL,        
                      pseudo_cohort_bl=mcdata, 
                      endtime=endtime, 
                      mod_e=mod_e, mod_w=mod_w, mod_d=mod_d, 
                      seed=seed)
  unex = gformula_risks(intervention=0,         
                        pseudo_cohort_bl=mcdata, 
                        endtime=endtime, 
                        mod_e=mod_e, mod_w=mod_w, mod_d=mod_d, 
                        seed=seed)
  exp3 = gformula_risks(intervention='rad>0.3', 
                        pseudo_cohort_bl=mcdata, 
                        endtime=endtime, 
                        mod_e=mod_e, mod_w=mod_w, mod_d=mod_d, 
                        seed=seed)
  exp1 = gformula_risks(intervention='rad>0.1', 
                        pseudo_cohort_bl=mcdata, endtime=endtime, 
                        mod_e=mod_e, mod_w=mod_w, mod_d=mod_d, 
                        seed=seed)
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
################################################################################################################################
gf_est <- gformula_effectestimates(data=baseline_data, fdata=miners, endtime=10, seed=NULL)
knitr::kable(gf_est, caption="1 Bootstrapped Resample: Risk Differences")


nbootsamples = 10
boot_samples <- boot(data = baseline_data, fdata=miners, statistic = gformula_effectestimates, R = nbootsamples, endtime=10)
knitr::kable(boot_samples$t0, caption="10 Bootstrapped Resamples: Risk Differences")

################################################################################################################################
est = boot_samples$t0 #the risk differences
lci =  boot_samples$t0 - 1.96*apply(boot_samples$t, 2, sd)
uci =  boot_samples$t0 + 1.96*apply(boot_samples$t, 2, sd)

dat1 <- cbind.data.frame(regime = c("Exp3","Exp1","Unex","Exp3","Exp1","Unex"),
                        type = c(rep("Risk Difference",3), rep("Log RR",3)),
                        RD = est, LB=lci, UB = uci) %>%
  mutate(regime = factor(regime, levels=c("Unex", "Exp1","Exp3"))) 
  
dat1 %>%
  knitr::kable(caption = "Risk Differences/Log Relative Risks and 95% Bootstrapped Confidence Intervals (R=10) Under Different Treatment Regimes at t=10 Using G-Formula")


#plot risk difference
ggplot(data=subset(dat1, type=="Risk Difference")) + 
  geom_point(aes(x=regime, y=RD)) +
  geom_errorbar(aes(ymin = LB, ymax=UB, x=regime))+
  facet_wrap(.~type) +
  ylab("Log RR/Risk Difference") +
  theme_bw() +
  ggtitle("G-Formula-Based Risk Differences Predictions Under Different Treatment Regimes") +
  guides(color='none')
  




#make table
#exponentiate the log RR to get RRs
dat1 %>%
  filter(type=="Log RR") %>%
  dplyr::rename(RR = RD) %>%
  mutate(across(c(RR,LB,UB), exp)) %>%
  print() %>%
ggplot() +
  geom_point(aes(x=regime, y=RR)) +
  geom_errorbar(aes(ymin = LB, ymax=UB, x=regime))+
  facet_wrap(.~type) +
  ylab("Risk Ratio") +
  theme_bw() +
  ggtitle("G-Formula-Based Risk Ratio Predictions Under Different Treatment Regimes") +
  guides(color='none') +
  geom_hline(yintercept = 1, linetype="dashed")




boot_samples200 <- boot(data = baseline_data, fdata=miners, statistic = gformula_effectestimates, R = 200, endtime=10)

est200 = boot_samples200$t0 #the risk differences
lci200 =  boot_samples200$t0 - 1.96*apply(boot_samples200$t, 2, sd)
uci200 =  boot_samples200$t0 + 1.96*apply(boot_samples200$t, 2, sd)

dat <- cbind.data.frame(regime = c("Exp3","Exp1","Unex","Exp3","Exp1","Unex"),
                        type = c(rep("Risk Difference",3), rep("Log RR",3)),
                        RD = est200, LB=lci200, UB = uci200) 
knitr::kable(dat, caption = "Risk Differences/Log Relative Risks and 95% Bootstrapped Confidence Intervals (R=200) Under Different Treatment Regimes at t=10 Using G-Formula")

ggplot(data=dat) + 
  geom_point(aes(x=regime, y=RD, color=type)) +
  geom_errorbar(aes(ymin = LB, ymax=UB, x=regime, color=type))+
  facet_wrap(.~type) +
  ylab("Log RR/Risk Difference") +
  theme_bw() +
  ggtitle("G-Formula-Based Risk Differences/Log Relative Risk Predictions Under Different Treatment Regimes") +
  theme(legend.position = "bottom")


# NA

knitr::purl(input = "longitudinal_gformula.Rmd", output = "longitudinal_gformula.R",documentation = 0)
