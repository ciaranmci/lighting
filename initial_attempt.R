# The purpose of this R project is to recreate the lighting-crime models that
# were written in MLWIN, in R.
#
# ********************
# ** OVERDISPERSION **
# ********************
# Paul raised concerns that the variates were overdispersed. If we are more 
# concerned with accounting for overdispersion rather than being interested in
# the describing the overdispersion, then we could just add a random intercept
# covariate that indicates the observation ID for observations. This is called
# the "observation-level random effect" (ORLE) approach (see https://peerj.com/articles/1114.pdf).
#
# Harrison 2014 and 2015 (https://peerj.com/articles/616.pdf;
# https://peerj.com/articles/1114.pdf) concluded:
# 1. For Poisson count data, ORLE works well when the overdispersion is caused
# by random noise in the counts, or by non-independence (a.k.a.aggregation) of
# the count in the count data.
#     - It is plausible that overdispersion of our `PropDarkCrime` data is due
#       to random noise in the counts.
#     - I assume non-independence is highly likely for crimes, but I would like
#       to cite research on, for example, crime sprees to justify my assumption.
# 2. For Poisson count data, ORLE also doesn't work well when the overdispersion
# is caused by excess zeros.
#     - Zero-inflation does not seem to be an issue for the DarknessCrime_sum`
#       and `DaylightCrime_sum` variables.
# 3. For binomial proportion data, ORLE works well when the overdispersion is
# caused by random noise in the proportions 
# 4. For binomial proportion data, ORLE doesn't work well when the overdispersion
# is caused by each observation arising from binomial distributions with
# possibly-different probabilities of success (i.e. when the proportion data is
# beta-binomial distributed).
#     - To tease out whether this might be the case for our proportions, we have
#       to clarify what distribution each observation is of. Are we saying that
#       each observation is a crime from a criminal's distribution, or is it a
#       crime from an MSOA's distribution? Fundamentally, we are asking: taking
#       a given crime, was the probability that the crime occurred at night
#       entirely independent of other crime's probability, or are all crimes
#       similarly probable of occurring at night? One side of the argument is
#       that a crime is a crime, so darkness should affect them all the same. The
#       other side of the argument is that our dataset is made up of many types
#       of crime, each of which have different affordances. I feel more
#       comfortable with the second view, which says that different types of
#       crimes are permitted to have different probabilities of occurring in the
#       hours of darkness. The overdispersion term in a beta-binomial regression
#       model works with the grand mean probability of a crime occurring in 
#       darkness to inform the shape and location of a beta distribution of
#       permissible probabilities; the larger the value of the overdispersion
#       term, the greater the degree of overdispersion in the data.
#     - Harrison's one caveat with using beta-binomial models is that they 
#       perform poorly if the data are not beta-binomial distributed.
#
# Given Harrison's findings, I suggest that ORLE is the way to go for the models
# of count data, and beta-binomial models should be used for models of proportions.
#
# R functions available for fitting beta-binomial models are:
# 1. spaMM::HLfit(), but I don't know what the `HLmethod` argument is.
# 2. glmmADMB::glmmadmb()
# 3. PROreg::BBmm(), but I don't understand most of the arguments.
# 4. runJAGS::run.jags(), but that is a mad Bayesian world through an emulator!
#
# Finally, resources for discussion about quasi-binomial and quasi-Poisson models:
#   1. https://online.stat.psu.edu/stat504/lesson/7/7.3
#   2. https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#fitting-models-with-overdispersion
#
#

###############
# Requisites. #
###############
# ----
if (!("pacman" %in% installed.packages()[,"Package"])) install.packages("pacman")
library(pacman)
pacman::p_load(tidyverse, haven, R2MLwiN, MASS, glmmADMB)
install.packages("glmmADMB",
                 repos=c("http://glmmadmb.r-forge.r-project.org/repos",
                         getOption("repos")),
                 type="source")
library(glmmADMB)
#
# file.show(system.file("demo", "UserGuide09.R", package = "R2MLwiN"))
#
# ----

##############
# Load data. #
##############
# ----
spssData <-
  haven::read_sav(file.choose()) %>%
  dplyr::mutate(across(everything(), as.vector)) %>%
  dplyr::mutate_at("MSOAN112", as.factor)
# ----

#######################
# Model specification # 
#######################
# ----
# Several models need to be fitted, varying in outcome variable and model type
#
#     OUTCOMES / VARIATES
#     -------------------
# 1. PropDarkCrime (the ratio of the number of crimes in darkness to all crimes)
#     numerator = ...
#
#
#     DISTRIBUTION
#     ------------
# For count outcome, Paul uses an extra Poisson distribution so that the 
# variance of the count is free to be estimated.
# To model a count outcome using the `RMLwiN::runMLwiN()` call, set the
# distribution argument to D = 'Binomial', and variate in formula to 
# logit(count), and estoptions = list(extra = TRUE).
#
# For proportion outcome, Paul uses an extra binomial distribution so that the
# variance of the proportion is free to be estimated.
# To model a proportion outcome using the `RMLwiN::runMLwiN()` call, set the
# distribution argument to D = 'Binomial', and the variate in formula to 
# logit(numerator, denominator), and estoptions = list(extra = TRUE).
#
#
#     ESTIMATION PROCEDURES
#     ---------------------
# The various estimating procedures for fitting the models are:
# 1. MQL1 = first-order marginal quasi-likelihood
# 2. PQL2 = second-order penalised quasi-likelihood
# 3. PQL6 = sixth-order penalised quasi-likelihood
# 
#                     **************************
#                    **** RMLwiN::runMLwiN() ****
#                     **************************
# The `RMLwiN::runMLwiN()` call can fit models using these estimating procedures.
# In the RMLwIN package, the nonlinear estimating procedure is specified by the 
# `nonlinear` option within the `estoptions` argument to the `runMLwiN()` call.
# For example, to specify a first-order marginal quasi-likelihood estimation,
# the argument must be `runMLwiN(..., estoptions = list(nonlinear = c(N = 0, M = 1), ...)`.
# The `N` sub-argument specifies either marginal (N = 0) or penalised (N = 1)
# quasi-likelihood linearization. The `M` sub-argument specifies either first-order
# (M = 1) or second-order (M = 2) approximation.
#
# When running the PQL estimating procedures, MLWinN recommend first running the
# MQL and then submitting that first pass to the PQL. In the example below, note
# that attributes of the MQL-estimated model are used as starting values for the
# PQL-estimated model via the `startval` sub-argument in the `estoptions`
# argument:
#
# mod_MQL <-
#   runMLwiN(
#     logit(variate) ~ 1 + covariate1_lvl1 + covariate2_lvl1 + (1 | covariate3_lvl2),
#     D = "Binomial"
#           )
#
# mod_PQL2 <-
#   runMLwiN(
#     logit(variate) ~ 1 + <fixed_covariate_i + (<random_slope_i> | <random_intercept_i>),
#     D = "Binomial",
#     estoptions = list(
#                     nonlinear = c(N = 1, M = 2),
#                     startval = list(
#                                   FP.b = mod_MQL@FP,
#                                   FP.v = mod_MQL@FP.cov,
#                                   RP.b = mod_MQL@RP,
#                                   RP.v = mod_MQL@RP.cov
#                                   )
#                     )
#             )
#
#
#                     ***********************
#                    **** MASS::glmmPQL() ****
#                     ***********************
# Unlike `RMLwiN::runMLwiN()`, the `MASS::glmmPQL()` fits a PQL model straight
# away with much cleaner syntax:
#     MASS::glmmPQL(
#       fixed = <variate> ~ <fixed_covariate_i>,
#       random = ~ <random_slope_i> | <random_intercept_i>,
#       family = binomial
#       )
# 
#
#                     *********************
#                    **** lme4::glmer() ****
#                     *********************
# The `lme4::glmer()` function can be used to fit a observation-level random
# effect model:
#     lme4::glmer(
#         cbind(<indicator of success>,
#               <indicator of trials> - <indicator of success>
#               ) ~ 
#          <fixed_covariate_i> + (<random_slope_i> | <random_intercept_i>) + (1 | <obs_id>),
#         family = binomial(logit)
#         )
#
#
#     FIXED EFFECTS
#     -------------
# 1.  `DiffLfromMnBy100` (Number of lamps minus its MSOA mean number of lamps)
# 2.  `DiffMSOA_MnLfromGrandMnBy100` (The MSOA's arithmetic mean number of lamps
#                                   minus the grand mean number of lamps)
# 3.  `MonthMidWk_2` (February indicator variable)
# 4.  `MonthMidWk_3` (March indicator variable)
# 5.  `MonthMidWk_4` (April indicator variable)
# 6.  `MonthMidWk_5` (May indicator variable)
# 7.  `MonthMidWk_6` (June indicator variable)
# 8.  `MonthMidWk_7` (July indicator variable)
# 9.  `MonthMidWk_8` (August indicator variable)
# 10. `MonthMidWk_9` (Septemeber indicator variable)
# 11. `MonthMidWk_10` (October indicator variable)
# 12. `MonthMidWk_11` (November indicator variable)
# 13. `MonthMidWk_12` (December indicator variable)
# 14. `NewYrWk` (Indicator for the week of New Year)
# 15. `GoodFriWk` (Indicator for the week of Good Friday)
# 16. `EasterWk` (Indicator for the week of Easter)
# 17. `MayDayWk` (Indicator for the week of May Day bank holiday)
# 18. `SpringBankWk` (Indicator for the week of SPring bank holiday)
# 19. `SummerBankWk` (Indicator for the week of Summer bank holiday)
# 20. `XmasWk` (Indicator for the week of Christmas)
# 21. `DifferenceInOffsets` (Difference between the natural logarithm of a given
#                            weeks' duration of darkness as a fraction of 24 hrs, 
#                            and the natural logarithm of a given weeks' duration
#                            of daylight as a fraction of 24 hrs)
# 22. `T_Dec4` (`T_Dec` raised to the fourth power; see RANDOM EFFECTS for
#               details of `T_Dec`)
# 23. `T_Dec5` (`T_Dec` raised to the fifth power)
#
#
#     RANDOM EFFECTS - intercepts - level 1
#     -------------------------------------
# 1. `cons` (constant, a vector of 1s)
#
#     RANDOM EFFECTS - intercepts - level 2
#     -------------------------------------
# 1. `MSOAN112` (Middle Layer Super Output Areas up to #112)
#
#     RANDOM EFFECTS - slopes - level 1
#     ---------------------------------
# NONE
#
#     RANDOM EFFECTS - slopes - level 2
#     ---------------------------------
# 1. `T_Dec`  (The decadal time at the middle of the given week, relative to the
#              10-year window under study; used as a scaling variable to scale
#              the coefficients)
# 2. `T_Dec2` (`T_Dec` raised to the second power)
# 3. `T_Dec3` (`T_Dec` raised to the third power)





# glmmADMB::glmmadmb()
# ## Just respect the nested data structure; ignore overdispersion.
mod_ADMB_1 <- 
  glmmADMB::glmmadmb(
    cbind(DarknessCrime_sum, SumDarkAndDaylight - DarknessCrime_sum) ~ 1 + (1 | MSOAN112),
    family = "betabinomial",
    data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111)),
    debug = TRUE
  )

# Poisson regression with:
# 1. Random intercept for `MSOAN112`, to respect nested data structure.
# 2. Random intercept for `CaseID`, to account for overdispersion via OLRE.
mod_GLMR_1 <-
  lme4::glmer(
    DarknessCrime_sum ~
      1 +
      (1 | CaseID) +
      (1 | MSOAN112),
    family = "poisson",
    data = spssData
  )
# RESULTS:
# - AIC = 223,159.5
# - BIC = 223,185.9
# - Variance of `CaseID` = 0.2787
#
# Next step is to add fixed effects to see if the fit improves.
mod_GLMR_2 <-
  lme4::glmer(
    DarknessCrime_sum ~
      1 +
      # FIXED EFFECTS.
      # ## Level 1
      DiffLfromMnBy100 +
      MonthMidWk_2 +
      MonthMidWk_3 +
      MonthMidWk_4 +
      MonthMidWk_5 +
      MonthMidWk_6 +
      MonthMidWk_7 +
      MonthMidWk_8 +
      MonthMidWk_9 +
      MonthMidWk_10 +
      MonthMidWk_11 +
      MonthMidWk_12 +
      NewYrWk +
      GoodFriWk +
      EasterWk +
      MayDayWk +
      SpringBankWk +
      SummerBankWk +
      XmasWk +
      DifferenceInOffsets +
      T_Dec4 +
      T_Dec5 +
      # ## Level 2.
      DiffMSOA_MnLfromGrandMnBy100 +
      # RANDOM EFFECTS
      (1 | CaseID) +  # for overdispersion.
      (1 | MSOAN112), # to respect nested data structure.
    family = "poisson",
    data = spssData
  )





R2MLwiN::runMLwiN(
  # OUTCOME.
  logit(DarknessCrime_sum, SumDarkAndDaylight) ~ 
    # FIXED EFFECTS.
    DiffLfromMnBy100 +
    DiffMSOA_MnLfromGrandMnBy100 +
    MonthMidWk_2 +
    MonthMidWk_3 +
    MonthMidWk_4 +
    MonthMidWk_5 +
    MonthMidWk_6 +
    MonthMidWk_7 +
    MonthMidWk_8 +
    MonthMidWk_9 +
    MonthMidWk_10 +
    MonthMidWk_11 +
    MonthMidWk_12 +
    NewYrWk +
    GoodFriWk +
    EasterWk +
    MayDayWk +
    SpringBankWk +
    SummerBankWk +
    XmasWk +
    DifferenceInOffsets +
    T_Dec4 +
    T_Dec5 +
    # LEVEL 2 RANDOM EFFECTS intercept and slopes, correlated
    (T_Dec + T_Dec2 + T_Dec3| MSOAN112) +
    # LEVEL 1 RANDOM EFFECTS intercept.
    (1 | cons),
  # DISTRIBUTION.
  D = "Binomial",
  estoptions = list(
                    # DISTRIBUTION.
                    extra = TRUE,
                    # ESTIMATION.
                    nonlinear = c(N = 1, M = 2),
                    startval = list(
                                  FP.b = mymodel4@FP,
                                  FP.v = mymodel4@FP.cov,
                                  RP.b = mymodel4@RP,
                                  RP.v = mymodel4@RP.cov
                                  )
                    ),
  data = spssData)

# ----


#############
# Model fit #
#############
# ----
# https://www.rdocumentation.org/packages/MuMIn/versions/1.47.5/topics/QAIC
# ----