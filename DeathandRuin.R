cat("\014")
chunkNumber <- 1000
chunkSize <- 1000
samplesForGraphs <- 1000
yearsReturns <- 65
seed1 <- 0114
cores <- 3


#Libraries
library(gdata)
library(survival)
library(cmprsk)
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(cmprsk)
library(reshape2)


#Files
men_risk_death <- tbl_df(read.xls("Table02.xlsx", skip = 2)) %>% #CDC US 2011 Male life table
  select(Age..years., qx) %>%
  mutate(age = row_number() - 1) %>%
  filter(age >= 65)

women_risk_death <- tbl_df(read.xls("Table03.xlsx", skip = 2)) %>% #CDC US 2011 Female life table
  select(Age..years., qx) %>%
  mutate(age = row_number() - 1) %>%
  filter(age >= 65)


#Make random male, female, and joint times of death packaged as a "data_frame"
make.frame <- function (nc, sseed, type) {
  set.seed(sseed)
  cohort <- data_frame(age = rep(65, nc), dead = rep(0, nc), id = 1:nc)
  mcohort <- transmute(cohort, male_age = age, male_dead = dead, id = id)
  fcohort <- transmute(cohort, female_age = age, female_dead = dead, id = id)
  
  for(ageyr in 65:99) {
    dead <- filter(mcohort, male_dead == 1)
    alive <- filter(mcohort, male_dead == 0)
    
    aliveNum <- length(alive$male_dead)
    riskDeath <- filter(men_risk_death, age == ageyr)$qx
    
    alive <- alive %>%
      mutate(male_age = male_age + 1,
             male_dead = rbinom(aliveNum, 1, riskDeath))
    
    mcohort <- bind_rows(dead, alive) %>%
      mutate(passseed = sseed)
  }
  
  for(ageyr in 65:99) {
    dead <- filter(fcohort, female_dead == 1)
    alive <- filter(fcohort, female_dead == 0)
    
    aliveNum <- length(alive$female_dead)
    riskDeath <- filter(women_risk_death, age == ageyr)$qx
    
    alive <- alive %>%
      mutate(female_age = female_age + 1,
             female_dead = rbinom(aliveNum, 1, riskDeath))
    
    fcohort <- bind_rows(dead, alive)
  }
  
  cohort <- arrange(inner_join(mcohort, fcohort, by = "id"), id) %>%
    mutate(female_tail = floor(rgamma(nc, shape = 1, scale = 2.3)),
           male_tail = floor(rgamma(nc, shape = 1, scale = 2.1)),
           male_dead = ifelse(male_age == 100, 1, male_dead),
           female_dead = ifelse(female_age == 100, 1, female_dead),
           male_age = ifelse(male_age == 100, male_age + male_tail, male_age),
           female_age = ifelse(female_age == 100, female_age + female_tail, female_age))
  
  if (sum(cohort$male_dead) != nc) warning("men alive")
  if (sum(cohort$female_dead) != nc) warning("women alive")
  
  if (type == "male") {
    cohort <- transmute(cohort, age_death = male_age, id = id, seed = sseed)
  }
  
  if (type == "female") {
    cohort <- transmute(cohort, age_death = female_age, id = id, seed = sseed)
  }
  
  if (type == "joint") {
    cohort <- transmute(cohort, age_death = pmax(male_age, female_age), id = id, seed = sseed)
  }
  
  return(cohort)
}


#Add Portfolio failure
add.pfailure <- function(cohort, meanreturns, sdreturns, initialEgg, annualWithdrawal, returns = F) {
  set.seed(cohort$passseed[1])
  nc <- length(cohort$age_death)
  matrixOfReturns <- matrix(rlnorm(nc * yearsReturns, log(meanreturns), sdreturns), ncol = yearsReturns)
  
  cohort <- mutate(cohort, balance = initialEgg, failure_year = 999)
  for(i in 1:yearsReturns) {
    cohort <- mutate(cohort, balance = balance - annualWithdrawal,
                     failure_year = ifelse(failure_year == 999 & balance < 0, i + 65, failure_year),
                     balance = balance*matrixOfReturns[id, i],
                     yearevent = pmin(age_death, failure_year),
                     event = ifelse(failure_year < age_death, 1, 0))
  }
  cohort <- select(cohort, -i, -balance, -seed)
  
  if(returns == T) {
    cohort <- bind_cols(cohort, mutate_each(data.frame(matrixOfReturns), funs(mone = . - 1)))
  }
  return(cohort)
}


#Death Sims
male <- data_frame(seed = seed1:(seed1 + chunkNumber - 1)) %>%
  group_by(seed) %>%
  do(cohort = make.frame(chunkSize, .$seed, "male"),
     type = rep("Single Man", chunkSize))
female <- data_frame(seed = seed1:(seed1 + chunkNumber - 1)) %>%
  group_by(seed) %>%
  do(cohort = make.frame(chunkSize, .$seed, "female"),
     type = rep("Single Woman", chunkSize))
both <- data_frame(seed = seed1:(seed1 + chunkNumber - 1)) %>%
  group_by(seed) %>%
  do(cohort = make.frame(chunkSize, .$seed, "joint"),
     type = rep("Joint Man/Woman", chunkSize))

three.categories <- bind_rows(male, female, both) %>% unnest()

three.fit <- survfit(Surv(age_death, rep(1, length(id))) ~ type, data = three.categories)

strat <- vector("integer", 0)
for(i in 1:3) strat <- append(strat, rep(names(three.fit$strata[i]), three.fit$strata[i]))
real.death.plot <- data_frame(time = three.fit$time, surv = three.fit$surv, stratum = strat) %>%
  mutate(stratum = substr(stratum, 6, 99))
rm(strat)

filter(real.death.plot, time == 95)

pmort <- ggplot(aes(x = time, y = surv), data = real.death.plot) +
  geom_step(aes(linetype = stratum)) +
  theme_classic(base_size = 12) +
  theme(text=element_text(family="Times")) +
  theme(legend.position = c(0.8, 0.8),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_x_continuous("Years of Age", limits = c(65, 105)) +
  scale_y_continuous("Proportion of Retirees Living") +
  labs(linetype = "Household Type:")
print(pmort)
ggsave("mortality.png", dpi = 400, width = 4, height = 3)

#Deterministic cases
#Run the simulation
singlesim <- list()
listCounter <- 0
for(meanReturns in seq(1.02, 1.06, 0.01)) {
  for(propWithdrawal in seq(0.03, 0.05, 0.01)) {
    listCounter <- listCounter + 1
    #Container of make.frames
    container <- data_frame(seed = seed1:(seed1 + chunkNumber - 1)) %>%
      group_by(seed) %>%
      do(cohort = make.frame(chunkSize, .$seed, "male"))
    
    #Add portfolio failures
    container$cohort <- mclapply(X = as.list(container$cohort), FUN = add.pfailure,
                                 meanreturns = meanReturns,
                                 sdreturns = 0.11,
                                 initialEgg = 1000000,
                                 annualWithdrawal = 1000000*propWithdrawal,
                                 mc.cores = cores)
    
    singlesim[[listCounter]] <- mutate(container, meanReturns = meanReturns, propWithdrawal = propWithdrawal)
  }
}

fits.groups <- bind_rows(singlesim) %>%
  unnest() %>%
  group_by(meanReturns, propWithdrawal) %>%
  do(fit = survfit(Surv(.$yearevent, .$event) ~ 1, data = .))

unnested.pile <- fits.groups %>%
  group_by(meanReturns, propWithdrawal) %>%
  do(time = unlist(.$fit, recursive = F)$time,
     surv = unlist(.$fit, recursive = F)$surv,
     lower = unlist(.$fit, recursive = F)$lower,
     upper = unlist(.$fit, recursive = F)$upper) %>%
  unnest() %>%
  mutate(Return = meanReturns, Withdrawal = propWithdrawal)

pmatrix <- ggplot(data = filter(unnested.pile, Return != 1.06), aes(x = time, y = surv)) +
  facet_grid(Return ~ Withdrawal, labeller = label_both) +
  theme_classic(base_size = 8) +
  theme(axis.title = element_text(size=12)) +
  theme(text=element_text(family="Times")) +
  scale_y_continuous("Proportion of Portfolios Surviving \nAmong Living Retirees",
                     limits = c(0, 1)) +
  scale_x_continuous(limits = c(65, 110), "Years of Age") +
  geom_step() +
  geom_ribbon(aes(ymax = upper, ymin = lower))
print(pmatrix)
ggsave("matrix.png", dpi = 400, width = 4, height = 3)

#Timing of risk
container <- data_frame(seed = seed1:(seed1 + 100 - 1)) %>%
  group_by(seed) %>%
  do(cohort = make.frame(chunkSize, .$seed, "male"))

container$cohort <- mclapply(X = as.list(container$cohort), FUN = add.pfailure,
                             meanreturns = 1.04,
                             sdreturns = 0.11,
                             initialEgg = 1000000,
                             annualWithdrawal = 1000000*0.04,
                             returns = T,
                             mc.cores = cores)

ds <- unnest(container) %>% mutate(trimyearevent = yearevent - 65)

xnam <- paste0("X", 1:65)
fmla <- as.formula(paste("Surv(trimyearevent, event) ~ ", paste(xnam, collapse= " + ")))


nonparam.model <- coxph(fmla, data = ds, robust = T)

estimates <- data_frame(year = 1:65,
                        estimate = exp(coef(nonparam.model)),
                        lower = exp(confint(nonparam.model))[, 1],
                        upper = exp(confint(nonparam.model))[, 2])

ggplot(aes(x = year, y = estimate), data = estimates) +
  geom_ribbon(aes(ymax = upper, ymin = lower)) +
  geom_line()

# cmprsk
ci <- cuminc(ds$yearevent, ds$event, cencode = 2)

ci.plot <- bind_rows(data_frame(time = ci$`1 0`$time, est = ci$`1 0`$est, var = ci$`1 0`$var,
                                label = "Death"),
                     data_frame(time = ci$`1 1`$time, est = ci$`1 1`$est, var = ci$`1 1`$var,
                                label = "Ruin")) %>%
  mutate(lower = est - (1.96 * var**0.5),
         upper = est + (1.96 * var**0.5))

pcmp <- ggplot(aes(x = time, y = est), data = ci.plot) +
  geom_step(aes(linetype = label)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, alpha = label)) +
  theme_classic(base_size = 12) +
  theme(text=element_text(family="Times")) +
  theme(legend.position = c(0.2, 0.8),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_x_continuous("Years of Age") +
  scale_y_continuous("Proportion of Retirees with Event") +
  coord_cartesian(xlim = c(65, 110), ylim = c(0, 1)) +
  labs(linetype = "Event:") +
  guides(alpha = F)
print(pcmp)
ggsave("cmprisk.png", dpi = 400, width = 4, height = 3)


#Comparison of measures
su <- survfit(Surv(yearevent, event) ~ 1, data = ds)
us <- survfit(Surv(failure_year, rep(1, length(id))) ~ 1, data = ds)
dp <- filter(ci.plot, label == "Ruin")

three.measures <- bind_rows(data_frame(time = su$time, risk = 1 - su$surv,
                                       type = "Conditional"),
                            data_frame(time = us$time, risk = 1 - us$surv,
                                       type = "Without death"),
                            data_frame(time = dp$time, risk = dp$est,
                                       type = "Absolute"))

pmeasure <- ggplot(aes(x = time, y = risk), data = three.measures) +
  geom_step(aes(linetype = type)) +
  theme_classic(base_size = 12) +
  theme(text=element_text(family="Times")) +
  theme(legend.position = c(0.2, 0.8),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  scale_x_continuous("Years of Age") +
  scale_y_continuous("Proportion of Ruined Retirees") +
  coord_cartesian(xlim = c(65, 110), ylim = c(0, 0.4)) +
  labs(linetype = "Event:")
print(pmeasure)
ggsave("measures.png", dpi = 400, width = 4, height = 3)

#Three cmprisk scenarios
three.cmps <- bind_rows(singlesim) %>%
  unnest() %>%
  filter(meanReturns == 1.04)

ci.three <- cuminc(three.cmps$yearevent, three.cmps$event, cencode = 2, group = three.cmps$propWithdrawal)

cit.plot <- bind_rows(data_frame(time = ci.three$`0.03 0`$time, risk = ci.three$`0.03 0`$est,
                                 withdrawal = "3% Annual Withdrawal", event = "Death"),
                      data_frame(time = ci.three$`0.03 1`$time, risk = ci.three$`0.03 1`$est,
                                 withdrawal = "3% Annual Withdrawal", event = "Ruin"),
                      data_frame(time = ci.three$`0.04 0`$time, risk = ci.three$`0.04 0`$est,
                                 withdrawal = "4% Annual Withdrawal", event = "Death"),
                      data_frame(time = ci.three$`0.04 1`$time, risk = ci.three$`0.04 1`$est,
                                 withdrawal = "4% Annual Withdrawal", event = "Ruin"),
                      data_frame(time = ci.three$`0.05 0`$time, risk = ci.three$`0.05 0`$est,
                                 withdrawal = "5% Annual Withdrawal", event = "Death"),
                      data_frame(time = ci.three$`0.05 1`$time, risk = ci.three$`0.05 1`$est,
                                 withdrawal = "5% Annual Withdrawal", event = "Ruin"))

p3cmp <- ggplot(aes(x = time, y = risk), data = cit.plot) +
  facet_grid(. ~ withdrawal) +
  geom_step(aes(linetype = event)) +
  theme_classic(base_size = 9) +
  theme(text=element_text(family="Times")) +
  theme(legend.position = c(0.92, 0.3),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  theme(axis.title = element_text(size=12)) +
  scale_x_continuous("Years of Age") +
  scale_y_continuous("Proportion of Single Men with Event") +
  coord_cartesian(xlim = c(65, 110), ylim = c(0, 1)) +
  labs(linetype = "Event:")
print(p3cmp)
ggsave("threecmprsk.png", dpi = 400, width = 4, height = 3)

#
# Graph Rates of Ruin for Our model, Milevsky Formula, Bengen-style model
# 

withdraw <- c(.025,.03,.035,.04,.045,.05,.055)
ourModel <- c(0,.0097,.026,.056,.1,.161,.233)
milevsky <- c(.0229,.0411,.0657,.0966,.133,.174,.2195)
bengen <- c(.0012,.0083,.028,.063,.134,.23,.3353)
data <- data.frame(withdraw,ourModel,milevsky,bengen)
names(data) <- c("withdraw","Kaplan-Meier","Milevsky","Fixed Lifespan")

data_long <- melt(data, id="withdraw")  # convert to long format
colnames(data_long) <- c("Annual Withdrawal","Method","pRuin")

pother <- (ggplot(data=data_long,
                  aes(x=data_long$`Annual Withdrawal`, y=pRuin, linetype=Method)) +
             scale_fill_discrete(labels=c("Our Model","Milevsky","Bengen Model")) +
             geom_line() +
             # ggtitle("Probability of Ruin by\n Three Different Estimating Methods") +
             labs(x="Annual Withdrawal Rate", y="Probability of Ruin") +
             theme(text=element_text(family="Times")) +
             theme(axis.title = element_text(size=12)) + 
             theme(legend.position = c(0.2, 0.8), legend.text = element_text(size = 10),
                   legend.title = element_text(size = 10)) +
             scale_y_continuous(breaks = seq(0,.5,by = .05)) )
print(pother)

ggsave (paste("Figure 2.png",sep=""),dpi=600, width = 4, height = 3)


