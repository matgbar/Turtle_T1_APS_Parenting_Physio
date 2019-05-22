###############################################################################
library(tidyverse)
library(yarrr)
library(nFactors)
library(brms)
library(rstan)
library(bayesplot)
library(mice)
library(mclust)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
###############################################################################


#Organizing data for parenting cluster analysis 
#Prior analyses revealed unidimensionality of FAS 
#Evidence that there were likely 2-factors for CRPR 
#CRPR items in current sample did not all load as expected according to scale instructions
#Used PCA to derive more "pure" scores of nurt & rest in this sample (performed by LD in SPSS)
###############################################################################
DF_kmeans<-dat_miss[,c("ID",
                       "FAS_avg_MR85",
                       "REST_PCA", 
                       "NURT_PCA")]

pairs.panels(DF_kmeans[,2:4])

#Cluster analysis to identify parenting profiles
###############################################################################
BIC_vals<-mclustBIC(na.omit(DF_kmeans[,2:4]))
plot(BIC_vals)
summary(BIC_vals)

mod<-Mclust(na.omit(DF_kmeans[,2:4]), x = BIC_vals)
summary(mod, parameters = TRUE)

plot(mod, what = "classification")

dat.tmp<-na.omit(DF_kmeans)
dat.tmp$class.new<-paste0("Class", mod$classification)

POS_T1_new<-merge(dat.tmp[,c(1,5)], dat_comb, by = "ID")

POS_T1_new$RSA_vid<-rowMeans(POS_T1_new[,c("RSA15_Video1",
                                           "RSA15_Video2",
                                           "RSA15_Video3")])
POS_T1_new$c.RSA_vid<-POS_T1_new$RSA_vid - mean(POS_T1_new$RSA_vid, na.rm=TRUE)
POS_T1_new$class.new<-as.factor(POS_T1_new$class.new) #apparently the marginal_effects function is too dumb to know what to do with a character value

psych::describe(POS_T1_new)

#Modeling Group Play as a function of Physio
###############################################################################

#Binomial - Intercept only model
#------------------------------------------------------------------------------
Grp_Form_int<-bf(
  Grp.tot | trials(Obs.tot) ~ 1 + (1|ID)
) + binomial()

system.time(
  Grp_int<-brm(formula = Grp_Form_int, 
              data = dat_mod,
              warmup = 9000, 
              iter = 10000, 
              chains = 3, 
              control = list(adapt_delta = .99, 
                             max_treedepth = 15))
)

#Initial Covariate Model: 
#------------------------------------------------------------------------------
#Missingness models need to be created for predictors
#There are only three cases missing data on school hours - will use intercept to predict
Grp_Form_CV<-bf(
  Grp.tot | trials(Obs.tot) ~ 1 + (1|ID)+
    CHILD_age_yrs + CHILD_sex + mi(schl_hrs_wkly)
) + binomial()

Schl_hrs_miss<-bf(
  schl_hrs_wkly | mi() ~ 1
) + gaussian()

system.time(
  Grp_CV<-brm(formula = Grp_Form_CV + 
                Schl_hrs_miss + 
                set_rescor(rescor = FALSE), 
              data = POS_T1_new,
              warmup = 9000, 
              iter = 10000, 
              chains = 3, 
              control = list(adapt_delta = .99, 
                             max_treedepth = 15))
)

summary(Grp_CV)

#-----------------------------------------------------------------------------
Grp_Form_CV_par<-bf(
  Grp.tot | trials(Obs.tot) ~ 1 + class.new + (1|ID)+
    CHILD_age_yrs + CHILD_sex + mi(schl_hrs_wkly)
) + binomial()

Schl_hrs_miss<-bf(
  schl_hrs_wkly | mi() ~ 1
) + gaussian()

system.time(
  Grp_CV_par<-brm(formula = Grp_Form_CV_par + 
                Schl_hrs_miss + 
                set_rescor(rescor = FALSE), 
              data = POS_T1_new,
              warmup = 9000, 
              iter = 10000, 
              chains = 3, 
              control = list(adapt_delta = .99, 
                             max_treedepth = 15))
)

summary(Grp_CV_par)

###############################################################################
###############################################################################
#RSA Models - Predicting Group Play: 
###############################################################################
###############################################################################

POS_T1_new$RSA_vid<-rowMeans(POS_T1_new[,c("RSA15_Video1",
                                           "RSA15_Video2",
                                           "RSA15_Video3")])

POS_T1_new$c.RSA_vid<-POS_T1_new$RSA_vid - mean(POS_T1_new$RSA_vid, na.rm=TRUE)
POS_T1_new$D_RSA_clown<-POS_T1_new$RSA15_Clown-POS_T1_new$RSA_vid
POS_T1_new$D_RSA_intro<-POS_T1_new$RSA15_Intro-POS_T1_new$RSA_vid
POS_T1_new$D_RSA_kids<-POS_T1_new$RSA15_Kids-POS_T1_new$RSA_vid


#Group play predicted by aggreagate video RSA
#Need to be able to include missingness model for RSA video predictor
#---------------------------------------------------------------------
POS_T1_new$miss_RSA_vid<-as.numeric(is.na(POS_T1_new$c.RSA_vid))
fit_RSA_Miss<-glm(miss_RSA_vid ~ 1 + 
                    BIQ_avg_MR85 + 
                    CBCL.anxd + 
                    CBCL.withd +
                    CBCL.soma +
                    CHILD_age_yrs + 
                    CHILD_sex, 
                  data = POS_T1_new, 
                  family = "binomial")

summary(fit_RSA_Miss)
step<-MASS::stepAIC(fit_RSA_Miss)
step$anova

#Using a stepwise procedure - looks as though CBCL.anxd is predictive of missingness

#Complete model with missingness component for missing RSA videos
#------------------------------------------------------------------------------
Grp_Form_RSA_vid<-bf(
  Grp.tot | trials(Obs.tot) ~ 1 + mi(c.RSA_vid)+
    (1|ID)+
    CHILD_age_yrs + CHILD_sex + mi(schl_hrs_wkly)
) + binomial()

bf_miss_RSA<-bf(
  c.RSA_vid | mi() ~ 1 + CBCL.anxd
) + gaussian()

system.time(
  Grp_RSA_vid<-brm(formula = Grp_Form_RSA_vid + 
                     bf_miss_RSA + 
                     Schl_hrs_miss + 
                     set_rescor(rescor = FALSE),
                   data = POS_T1_new,
                   warmup = 9000, 
                   iter = 10000, 
                   chains = 3, 
                   control = list(adapt_delta = .99, 
                                  max_treedepth = 15))
)

summary(Grp_RSA_vid)

tmp<-data.frame(c.RSA_vid = sample(seq(min(POS_T1_new$c.RSA_vid, na.rm = TRUE),
                                         max(POS_T1_new$c.RSA_vid, na.rm = TRUE), 
                                         by = .01), 
                                     size = 5000, replace = TRUE), 
                Obs.tot = rep(1, 5000),
                ID = sample(POS_T1_new$ID, size = 5000, replace = TRUE), 
                CHILD_age_yrs = sample(POS_T1_new$CHILD_age_yrs, size = 5000, replace = TRUE), 
                CHILD_sex = sample(POS_T1_new$CHILD_sex, size = 5000, replace = TRUE), 
                schl_hrs_wkly = sample(na.omit(POS_T1_new$schl_hrs_wkly), size = 5000, replace = TRUE), 
                CBCL.anxd = sample(POS_T1_new$CBCL.anxd, size = 5000, replace = TRUE)) 

pp_check(Grp_RSA_vid, resp = "Grptot", nsamples = 100)

fit2fitnew<-fitted(Grp_RSA_vid, 
                   newdata = tmp, 
                   re_formula = NA, 
                   resp = "Grptot")

tmp$Fit_vals<-fit2fitnew[,1]

g.fit.vid<-ggplot(data= tmp,
                    aes(x = c.RSA_vid, 
                        y = Fit_vals))+
  stat_smooth(lwd = 2,
              se = FALSE)+
  geom_point(alpha = .25)+
  labs(title="Predicted Group Play Proportions as a Function of RSA during Video Task", 
       x="RSA during Video Task", 
       y="Proportion of Group Play")+
  theme_bw()

g.fit.vid

g.raw.vid<-ggplot(data = POS_T1_new)+
  geom_point(size = 2.5, 
             shape = 3, 
             color = "red",
             aes(x = c.RSA_vid,
                 y = Grp.prop))+
  geom_smooth(data = POS_T1_new, 
              lwd = 2,
              color = "red",
              aes(x = c.RSA_vid,
                  y = Grp.prop),
              method = "lm", 
              se = FALSE)+
  labs(title="Raw Group Play Proportions as a Function of RSA during Video Task", 
       x="RSA during Video Task", 
       y="Proportion of Group Play")+
  theme_bw()

g.raw.vid

#Prediction missingness on change in RSA during Intro 
#------------------------------------------------------------------------------
POS_T1_new$miss_RSA_intro<-as.numeric(is.na(POS_T1_new$D_RSA_intro))
fit_RSA_Miss<-glm(miss_RSA_intro ~ 1 + 
                    BIQ_avg_MR85 + 
                    CBCL.anxd + 
                    CBCL.withd +
                    CBCL.soma +
                    CHILD_age_yrs + 
                    CHILD_sex, 
                  data = POS_T1_new, 
                  family = "binomial")

summary(fit_RSA_Miss)
step<-MASS::stepAIC(fit_RSA_Miss)
step$anova

#Complete model for change in RSA during Intro task (relative to aggregate vids)
#Includes a missingness model in which missigness depends ?????
#------------------------------------------------------------------------------
Grp_Form_RSA_intro<-bf(
  Grp.tot | trials(Obs.tot) ~ 1 + mi(D_RSA_intro)+
    (1|ID)+
    CHILD_age_yrs + CHILD_sex + mi(schl_hrs_wkly)
) + binomial()

bf_miss_RSA<-bf(
  D_RSA_intro | mi() ~ 1 + CBCL.anxd
) + gaussian()

system.time(
  Grp_RSA_intro<-brm(formula = Grp_Form_RSA_intro + 
                       bf_miss_RSA + 
                       Schl_hrs_miss + 
                       set_rescor(rescor = FALSE),
                     data = POS_T1_new,
                     prior = c(set_prior("normal(0,3)", class = "b"),
                               set_prior("normal(1.36, 1)", class = "sd", resp = "Grptot"), 
                               set_prior("normal(-.74, 3)", class = "Intercept", resp = "DRSAintro")),
                     warmup = 19000, 
                     iter = 20000, 
                     chains = 3, 
                     control = list(adapt_delta = .99, 
                                    max_treedepth = 15))
)

summary(Grp_RSA_intro)

tmp<-data.frame(D_RSA_intro = sample(seq(min(POS_T1_new$D_RSA_intro, na.rm = TRUE),
                                        max(POS_T1_new$D_RSA_intro, na.rm = TRUE), 
                                        by = .01), 
                                    size = 5000, replace = TRUE), 
                Obs.tot = rep(1, 5000),
                ID = sample(POS_T1_new$ID, size = 5000, replace = TRUE), 
                CHILD_age_yrs = sample(POS_T1_new$CHILD_age_yrs, size = 5000, replace = TRUE), 
                CHILD_sex = sample(POS_T1_new$CHILD_sex, size = 5000, replace = TRUE), 
                schl_hrs_wkly = sample(na.omit(POS_T1_new$schl_hrs_wkly), size = 5000, replace = TRUE), 
                CBCL.anxd = sample(POS_T1_new$CBCL.anxd, size = 5000, replace = TRUE)) 

pp_check(Grp_RSA_intro, resp = "Grptot", nsamples = 100)

fit2fitnew<-fitted(Grp_RSA_intro, 
                   newdata = tmp, 
                   re_formula = NA, 
                   resp = "Grptot")

tmp$Fit_vals<-fit2fitnew[,1]


g.fit.intro<-ggplot(data= tmp,
                    aes(x = D_RSA_intro, 
                        y = Fit_vals))+
  stat_smooth(lwd = 2,
              se = FALSE)+
  geom_point(alpha = .25)+
  labs(title=expression("Predicted Group Play Proportions as a Function of"~Delta~"RSA during Intro Task"), 
       x=expression(Delta~"RSA during Intro Task"), 
       y="Proportion of Group Play")+
  theme_bw()

g.fit.intro

g.raw.intro<-ggplot(data = POS_T1_new)+
  geom_point(size = 2.5, 
             shape = 3, 
             color = "red",
             aes(x = D_RSA_intro,
                 y = Grp.prop))+
  geom_smooth(data = POS_T1_new, 
              lwd = 2,
              color = "red",
              aes(x = D_RSA_intro,
                  y = Grp.prop),
              method = "lm", 
              se = FALSE)+
  labs(title=expression("Raw Group Play Proportions as a Function of"~Delta~"RSA during Intro Task"), 
       x=expression(Delta~"RSA during Intro Task"), 
       y="Proportion of Group Play")+
  theme_bw()

g.raw.intro

#Prediction missingness on change in RSA during kids 
#---------------------------------------------------------------------
POS_T1_new$miss_RSA_kids<-as.numeric(is.na(POS_T1_new$D_RSA_kids))
fit_RSA_Miss<-glm(miss_RSA_kids ~ 1 + 
                    BIQ_avg_MR85 + 
                    CBCL.anxd + 
                    CBCL.withd +
                    CBCL.soma +
                    CHILD_age_yrs + 
                    CHILD_sex, 
                  data = POS_T1_new, 
                  family = "binomial")

summary(fit_RSA_Miss)
step<-MASS::stepAIC(fit_RSA_Miss)
step$anova

#Complete model for change in RSA during Clown task (relative to aggregate vids)
#Includes a missingness mode in which missigness depends on intercept only (see stepwise above)
#------------------------------------------------------------------------------
Grp_Form_RSA_kids<-bf(
  Grp.tot | trials(Obs.tot) ~ 1 + mi(D_RSA_kids) +
    (1|ID) +
    CHILD_age_yrs + CHILD_sex + mi(schl_hrs_wkly)
) + binomial()

bf_miss_RSA<-bf(
  D_RSA_kids | mi() ~ 1
) + gaussian()

system.time(
  Grp_RSA_kids<-brm(formula = Grp_Form_RSA_kids + 
                       bf_miss_RSA + 
                       Schl_hrs_miss + 
                       set_rescor(rescor = FALSE),
                     data = POS_T1_new,
                    prior = c(set_prior("normal(0,3)", class = "b"),
                              set_prior("normal(1.36, 1)", class = "sd", resp = "Grptot")), 
                     warmup = 9000, 
                     iter = 10000, 
                     chains = 3, 
                     control = list(adapt_delta = .99, 
                                    max_treedepth = 15))
)

summary(Grp_RSA_kids)

tmp<-data.frame(D_RSA_kids = sample(seq(min(POS_T1_new$D_RSA_kids, na.rm = TRUE),
                                       max(POS_T1_new$D_RSA_kids, na.rm = TRUE), 
                                       by = .01), 
                                   size = 5000, replace = TRUE), 
                Obs.tot = rep(1, 5000),
                ID = sample(POS_T1_new$ID, size = 5000, replace = TRUE), 
                CHILD_age_yrs = sample(POS_T1_new$CHILD_age_yrs, size = 5000, replace = TRUE), 
                CHILD_sex = sample(POS_T1_new$CHILD_sex, size = 5000, replace = TRUE), 
                schl_hrs_wkly = sample(na.omit(POS_T1_new$schl_hrs_wkly), size = 5000, replace = TRUE)) 

pp_check(Grp_RSA_kids, resp = "Grptot", nsamples = 100)

fit2fitnew<-fitted(Grp_RSA_kids, 
                   newdata = tmp, 
                   re_formula = NA, 
                   resp = "Grptot")

tmp$Fit_vals<-fit2fitnew[,1]

g.fit.kids<-ggplot(data= tmp,
                    aes(x = D_RSA_kids, 
                        y = Fit_vals))+
  stat_smooth(lwd = 2,
              se = FALSE)+
  geom_point(alpha = .25)+
  labs(title=expression("Predicted Group Play Proportions as a Function of"~Delta~"RSA during Kids Task"), 
       x=expression(Delta~"RSA during Kids Task"), 
       y="Proportion of Group Play")+
  theme_bw()

g.fit.kids

g.raw.kids<-ggplot(data = POS_T1_new)+
  geom_point(size = 2.5, 
             shape = 3, 
             color = "red",
             aes(x = D_RSA_kids,
                 y = Grp.prop))+
  geom_smooth(data = POS_T1_new, 
              lwd = 2,
              color = "red",
              aes(x = D_RSA_kids,
                  y = Grp.prop),
              method = "lm", 
              se = FALSE)+
  labs(title=expression("Raw Group Play Proportions as a Function of"~Delta~"RSA during Kids Task"), 
       x=expression(Delta~"RSA during Kids Task"), 
       y="Proportion of Group Play")+
  theme_bw()

g.raw.kids



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Adding pareting profiles as predictors
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------



#Complete model with missingness component for missing RSA videos
#------------------------------------------------------------------------------
Grp_Form_RSA_vid_par<-bf(
  Grp.tot | trials(Obs.tot) ~ 1 + mi(c.RSA_vid) * class.new +
    (1|ID)+
    CHILD_age_yrs + CHILD_sex + mi(schl_hrs_wkly)
) + binomial()

bf_miss_RSA<-bf(
  c.RSA_vid | mi() ~ 1 + CBCL.anxd
) + gaussian()

system.time(
  Grp_RSA_vid_par<-brm(formula = Grp_Form_RSA_vid_par + 
                     bf_miss_RSA + 
                     Schl_hrs_miss + 
                     set_rescor(rescor = FALSE),
                   data = POS_T1_new,
                   warmup = 9000, 
                   iter = 10000, 
                   chains = 3, 
                   control = list(adapt_delta = .99, 
                                  max_treedepth = 15))
)

summary(Grp_RSA_vid_par)

POS_T1_new$class.new2<-ifelse(POS_T1_new$class.new=="Class1", 
                              "Nurturing", 
                              "High Nurturing")

tmp<-data.frame(c.RSA_vid = sample(seq(min(POS_T1_new$c.RSA_vid, na.rm = TRUE),
                                       max(POS_T1_new$c.RSA_vid, na.rm = TRUE), 
                                       by = .01), 
                                   size = 10000, replace = TRUE), 
                class.new = sample(POS_T1_new$class.new, size = 10000, replace = TRUE),
                Obs.tot = rep(1, 10000),
                ID = sample(POS_T1_new$ID, size = 10000, replace = TRUE), 
                CHILD_age_yrs = sample(POS_T1_new$CHILD_age_yrs, size = 10000, replace = TRUE), 
                CHILD_sex = sample(POS_T1_new$CHILD_sex, size = 10000, replace = TRUE), 
                schl_hrs_wkly = sample(na.omit(POS_T1_new$schl_hrs_wkly), size = 10000, replace = TRUE), 
                CBCL.anxd = sample(POS_T1_new$CBCL.anxd, size = 10000, replace = TRUE)) 

pp_check(Grp_RSA_vid_par, resp = "Grptot", nsamples = 100)

fit2fitnew<-fitted(Grp_RSA_vid_par, 
                   newdata = tmp, 
                   re_formula = NA, 
                   resp = "Grptot")

tmp$Fit_vals<-fit2fitnew[,1]
tmp$class.new2<-ifelse(tmp$class.new=="Class1", "Nurturing", "High Nurturing")

g.fit.vid_par<-ggplot(data= tmp,
                      aes(x = c.RSA_vid, 
                          y = Fit_vals, 
                          color = class.new2))+
  stat_smooth(lwd = 2,
              se = FALSE, 
              aes(group = class.new2))+
  geom_point(alpha = .25, 
             aes(color = class.new2))+
  labs(title=expression("Predicted Group Play Proportions as a Function of RSA during Video Task"), 
       x=expression("RSA during Video Task"), 
       y="Proportion of Group Play", 
       color = "Parent Profile")+theme_bw()

g.fit.vid_par

#------------------------------------------------------------------------------
Grp_Form_RSA_intro_par<-bf(
  Grp.tot | trials(Obs.tot) ~ 1 + mi(D_RSA_intro) * class.new +
    (1|ID)+
    CHILD_age_yrs + CHILD_sex + mi(schl_hrs_wkly)
) + binomial()

bf_miss_RSA<-bf(
  D_RSA_intro | mi() ~ 1 + CBCL.anxd
) + gaussian()

system.time(
  Grp_RSA_intro_par<-brm(formula = Grp_Form_RSA_intro_par + 
                           bf_miss_RSA + 
                           Schl_hrs_miss + 
                           set_rescor(rescor = FALSE),
                         data = POS_T1_new,
                         prior = c(set_prior("normal(0,3)", class = "b"),
                                   set_prior("normal(1.36, 1)", class = "sd", resp = "Grptot"), 
                                   set_prior("normal(-.74, 3)", class = "Intercept", resp = "DRSAintro")),
                         warmup = 14000, 
                         iter = 15000, 
                         chains = 3, 
                         control = list(adapt_delta = .99, 
                                        max_treedepth = 15))
)

summary(Grp_RSA_intro_par)  #Also needs more iterations

tmp<-data.frame(D_RSA_intro = sample(seq(min(POS_T1_new$D_RSA_intro, na.rm = TRUE),
                                         max(POS_T1_new$D_RSA_intro, na.rm = TRUE), 
                                         by = .01), 
                                     size = 10000, replace = TRUE), 
                class.new = sample(POS_T1_new$class.new, size = 10000, replace = TRUE),
                Obs.tot = rep(1, 10000),
                ID = sample(POS_T1_new$ID, size = 10000, replace = TRUE), 
                CHILD_age_yrs = sample(POS_T1_new$CHILD_age_yrs, size = 10000, replace = TRUE), 
                CHILD_sex = sample(POS_T1_new$CHILD_sex, size = 10000, replace = TRUE), 
                schl_hrs_wkly = sample(na.omit(POS_T1_new$schl_hrs_wkly), size = 10000, replace = TRUE), 
                CBCL.anxd = sample(POS_T1_new$CBCL.anxd, size = 10000, replace = TRUE)) 

pp_check(Grp_RSA_intro_par, resp = "Grptot", nsamples = 100)

fit2fitnew<-fitted(Grp_RSA_intro_par, 
                   newdata = tmp, 
                   re_formula = NA, 
                   resp = "Grptot")

tmp$Fit_vals<-fit2fitnew[,1]
tmp$class.new2<-ifelse(tmp$class.new=="Class1", "Nurturing", "High Nurturing")

g.fit.intro_par<-ggplot(data= tmp,
                       aes(x = D_RSA_intro, 
                           y = Fit_vals, 
                           color = class.new2))+
  stat_smooth(lwd = 2,
              se = FALSE, 
              aes(group = class.new2))+
  geom_point(alpha = .25, 
             aes(color = class.new2))+
  labs(title=expression("Predicted Group Play Proportions as a Function of"~Delta~"RSA during Intro Task"), 
       x=expression(Delta~"RSA during Intro Task"), 
       y="Proportion of Group Play", 
       color = "Parent Profile")

g.fit.intro_par

#------------------------------------------------------------------------------
Grp_Form_RSA_kids_par<-bf(
  Grp.tot | trials(Obs.tot) ~ 1 + mi(D_RSA_kids) * class.new +
    (1|ID)+
    CHILD_age_yrs + CHILD_sex + mi(schl_hrs_wkly)
) + binomial()

bf_miss_RSA<-bf(
  D_RSA_kids | mi() ~ 1
) + gaussian()

system.time(
  Grp_RSA_kids_par<-brm(formula = Grp_Form_RSA_kids_par + 
                           bf_miss_RSA + 
                           Schl_hrs_miss + 
                           set_rescor(rescor = FALSE),
                         data = POS_T1_new,
                        prior = c(set_prior("normal(0,3)", class = "b"),
                                  set_prior("normal(1.36, 3)", class = "sd", resp = "Grptot"), 
                                  set_prior("normal(.06, 3)", class = "Intercept", resp = "DRSAkids"), 
                                  set_prior("normal(28.96, 3)", class = "Intercept", resp = "schlhrswkly")),
                         warmup = 39000, 
                         iter = 40000, 
                         chains = 3, 
                         control = list(adapt_delta = .99, 
                                        max_treedepth = 15))
)

summary(Grp_RSA_kids_par) #Required lots of iterations!!!

tmp<-data.frame(D_RSA_kids = sample(seq(min(POS_T1_new$D_RSA_kids, na.rm = TRUE),
                                         max(POS_T1_new$D_RSA_kids, na.rm = TRUE), 
                                         by = .01), 
                                     size = 10000, replace = TRUE), 
                class.new = sample(POS_T1_new$class.new, size = 10000, replace = TRUE),
                Obs.tot = rep(1, 10000),
                ID = sample(POS_T1_new$ID, size = 10000, replace = TRUE), 
                CHILD_age_yrs = sample(POS_T1_new$CHILD_age_yrs, size = 10000, replace = TRUE), 
                CHILD_sex = sample(POS_T1_new$CHILD_sex, size = 10000, replace = TRUE), 
                schl_hrs_wkly = sample(na.omit(POS_T1_new$schl_hrs_wkly), size = 10000, replace = TRUE), 
                CBCL.anxd = sample(POS_T1_new$CBCL.anxd, size = 10000, replace = TRUE)) 

pp_check(Grp_RSA_kids_par, resp = "Grptot", nsamples = 100)

fit2fitnew<-fitted(Grp_RSA_kids_par, 
                   newdata = tmp, 
                   re_formula = NA, 
                   resp = "Grptot")

tmp$Fit_vals<-fit2fitnew[,1]
tmp$class.new2<-ifelse(tmp$class.new=="Class1", "Nurturing", "High Nurturing")

g.fit.kids_par<-ggplot(data= tmp,
              aes(x = D_RSA_kids, 
                  y = Fit_vals, 
                  color = class.new2))+
  stat_smooth(lwd = 2,
              se = FALSE, 
              aes(group = class.new2))+
  geom_point(alpha = .25, 
             aes(color = class.new2))+
  labs(title=expression("Predicted Group Play Proportions as a Function of"~Delta~"RSA during Kids Task"), 
       x=expression(Delta~"RSA during Kids Task"), 
       y="Proportion of Group Play", 
       color = "Parent Profile")

g.fit.kids_par

###############################################################################
#plots for poster
dat.tmp$class.new2<-ifelse(dat.tmp$class.new=="Class1", "Nurturing", "High Nurturing")
#Scatterplot with distributions for parenting class: 
g21<-ggplot(data = dat.tmp, 
           aes(x = FAS_avg_MR85, 
               y = REST_PCA, 
               group = class.new2, 
               color = class.new2, 
               shape = class.new2))+
  geom_point()+
  labs(y = "", 
       x = "", 
       color = "Parenting Profile", 
       shape = "Parenting Profile")+
  stat_ellipse()+
  theme_bw()
g21

g31<-ggplot(data = dat.tmp, 
            aes(x = FAS_avg_MR85, 
                y = NURT_PCA, 
                group = class.new2, 
                color = class.new2, 
                shape = class.new2))+
  geom_point()+
  labs(y = "", 
       x = "", 
       color = "Parenting Profile", 
       shape = "Parenting Profile")+
  stat_ellipse()+
  theme_bw()
g31

g32<-ggplot(data = dat.tmp, 
            aes(x = REST_PCA, 
                y = NURT_PCA, 
                group = class.new2, 
                color = class.new2, 
                shape = class.new2))+
  geom_point()+
  labs(y = "", 
       x = "", 
       color = "Parenting Profile", 
       shape = "Parenting Profile")+
  stat_ellipse()+
  theme_bw()

g32

legend <- cowplot::get_legend(g21)

g11<-ggplot(data = dat.tmp, 
            aes(x = FAS_avg_MR85, 
                y = class.new2, 
                color = class.new2, 
                fill = class.new2))+
  geom_density_ridges()+
  theme_bw()+
  theme(axis.text.y = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = "none")+
  labs(y = "", 
       x = "", 
       title = "Accommodation")
g11

g22<-ggplot(data = dat.tmp, 
            aes(x = REST_PCA, 
                y = class.new2, 
                color = class.new2, 
                fill = class.new2))+
  geom_density_ridges()+
  theme_bw()+
  theme(axis.text.y = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = "none")+
  labs(y = "", 
       x = "", 
       title = "Restrictiveness")
g22

g33<-ggplot(data = dat.tmp, 
            aes(x = NURT_PCA, 
                y = class.new2, 
                color = class.new2, 
                fill = class.new2))+
  geom_density_ridges()+
  theme_bw()+
  theme(axis.text.y = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = "none")+
  labs(y = "", 
       x = "", 
       title = "Nurturance")
g33

g12<-ggplot(data = dat.tmp, 
            aes(y = FAS_avg_MR85, 
                x = REST_PCA, 
                group = class.new2, 
                fill = class.new2, 
                shape = class.new2))+
  stat_density2d(geom="tile", aes(alpha=..density..), contour=FALSE)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(y="", 
       x="")

g12

g13<-ggplot(data = dat.tmp, 
            aes(y = FAS_avg_MR85, 
                x = NURT_PCA, 
                group = class.new2, 
                fill = class.new2, 
                shape = class.new2))+
  stat_density2d(geom="tile", aes(alpha=..density..), contour=FALSE)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(y="", 
       x="")

g13

g23<-ggplot(data = dat.tmp, 
            aes(y = REST_PCA, 
                x = NURT_PCA, 
                group = class.new2, 
                fill = class.new2, 
                shape = class.new2))+
  stat_density2d(geom="tile", aes(alpha=..density..), contour=FALSE)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(y="", 
       x="")

g23


png(paste0(graphics.folder, "/Grid_Plot_pareting.png"), 
    res = 1200, 
    units = "in", 
    width = 10, 
    height = 8)
cowplot::plot_grid(g11, g12, g13, ggplot()+theme_void(),
                   g21 + theme(legend.position="none"), g22, g23, legend,
                   g31 + theme(legend.position="none"), g32 + theme(legend.position="none"), g33, ggplot()+theme_void(),
                   align = 'v',
                   hjust = -1,
                   nrow = 3, 
                   axis = 'l'
)
dev.off()

#------------------------------------------------------------------------------
png(paste0(graphics.folder, "/Main_RSA_models.png"), 
    res = 1200, 
    units = "in", 
    width = 14, 
    height = 10)
cowplot::plot_grid(g.fit.vid, g.raw.vid, 
                   g.fit.intro, g.raw.intro,
                   g.fit.kids, g.raw.kids,
                   align = 'h',
                   hjust = -1,
                   nrow = 3, 
                   axis = 'b')
dev.off()
#------------------------------------------------------------------------------
png(paste0(graphics.folder, "/Parenting_x_RSA_models.png"), 
    res = 1200, 
    units = "in", 
    width = 14, 
    height = 10)
cowplot::plot_grid(g.fit.vid_par+theme_bw(), 
                   g.fit.intro_par+theme_bw(), 
                   g.fit.kids_par+theme_bw(),
                   align = 'v',
                   hjust = -1,
                   nrow = 3, 
                   axis = 'l')
dev.off()


