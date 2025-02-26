knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)


#Load in libraries
library(ggplot2) #for plotting
library(dplyr) #for data cleaning
library(survival) #for survival analysis
library(boot) #for bootstrapping
library(ggcorrplot)


coalplant <- read.csv("data/coalplant.csv")
str(coalplant)



exposures <- c("as", "be","cd")
corr <- cor(coalplant[, exposures])
ggcorrplot(corr)


mdimod = lm(mdi ~ as*urbanicity + be*urbanicity + cd*urbanicity + as*be + as*cd + be*cd, data=coalplant)
plot(mdimod)
summary(mdimod)



mean(predict(mdimod))

#create new data set for prediction
dat_nc <- coalplant %>%
    mutate(as = as*(1-0.00), #no reduction in as
         be = be*(1-0.00), #no reduction in be
         cd = cd*(1-0.00)) #no reduction in cd
  
#predict mean MDI 
(nc <- mean(predict(mdimod, newdata = dat_nc)))
  

#create new data set based on new intervention
dat_no_coal <- coalplant %>%
  mutate(as = as*(1-0.96),
         be = be*(1-0.91),
         cd = cd*(1-0.45))
#predict mean MDI 
(no_coal <- mean(predict(mdimod, newdata = dat_no_coal)))
  
  

dat_no_as <- coalplant %>%
  mutate(as = as*0)
#predict mean MDI 
(no_as <- mean(predict(mdimod, newdata = dat_no_as)))

dat_no_be <- coalplant %>%
  mutate(be = be*0)
#predict mean MDI 
(no_be <- mean(predict(mdimod, newdata = dat_no_be)))

dat_no_cd <- coalplant %>%
  mutate(cd = cd*0)
#predict mean MDI 
(no_cd <- mean(predict(mdimod, newdata = dat_no_cd)))

c(naturalcourse = nc, 
  no_coal = no_coal, 
  no_arsenic = no_as, 
  no_beryllium = no_be, 
  no_cadmium = no_cd) %>%
  knitr::kable(col.names = c("Model","Predicted Mean MDI"))


no_coal-nc

################################################################################################################################
bootfunc <- function(data, index){
  mcsample = data[index,] #resampled data set
  bootmod = lm(mdi ~ as*urbanicity + be*urbanicity + cd*urbanicity + as*be + as*cd + be*cd, data=mcsample) #model
  
  #scenarios
  ###natural course
  dat_nc <- mcsample %>%
    mutate(as = as*(1-0.00), #no reduction in as
           be = be*(1-0.00), #no reduction in be
           cd = cd*(1-0.00)) #no reduction in cd
  nc <- mean(predict(bootmod, newdata = dat_nc))
  
  ###no coal-fire power plant
  dat_no_coal <- mcsample  %>%
  mutate(as = as*(1-0.96),
         be = be*(1-0.91),
         cd = cd*(1-0.45))
  no_coal <- mean(predict(bootmod, newdata = dat_no_coal))
  
  ###no arsenic
  dat_no_as <- mcsample  %>%
  mutate(as = as*0)
  no_as <- mean(predict(bootmod, newdata = dat_no_as))
  
  ###no beryllium
  dat_no_be <- mcsample %>%
  mutate(be = be*0)
  no_be <- mean(predict(bootmod, newdata = dat_no_be))
  
  ###no cadmium
  dat_no_cd <- mcsample  %>%
  mutate(cd = cd*0)
  no_cd <- mean(predict(bootmod, newdata = dat_no_cd))
  
  #combine estiamtes together to export
  c(nc_mdi=nc, shutdown=no_coal-nc, attrmdi_cd=no_cd-nc, attrmdi_as=no_as-nc, attrmdi_be=no_be-nc)
  
}
################################################################################################################################


set.seed(2)
bootsamples = boot(data=coalplant, statistic=bootfunc, R=500)
bootsamples


bootsamples$t[1:20,]

se = apply(bootsamples$t, 2, sd)
se

print(cbind(estimate=bootsamples$t0, lowerCI=bootsamples$t0-1.96*se, upperCI=bootsamples$t0+1.96*se), 3)


