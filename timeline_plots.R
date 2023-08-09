# The script is for producing timeline plots of both day and darkness crimes
# over time, and overlaying the expectations of increasingly-complicated
# Poisson regression GLMs.


###############
# Requisites. #
###############
# ----
if (!("pacman" %in% installed.packages()[,"Package"])) install.packages("pacman")
library(pacman)
pacman::p_load(tidyverse, haven)
# ----


##############
# Load data. #
##############
# ----
if (!exists(spssData))
{
  spssData <-
    haven::read_sav(file.choose()) %>%
    mutate(across(everything(), as.vector))
}
# It seems that there is one MSOA that has distinctly different crime counts.
# It is MSOA = 111. It has 103 crimes in the hours of darkness within the first
# week, compared to the next highest of 19.
spssData %>%
  dplyr::select(MSOAN112, DarknessCrime_sum, WkNoStartFrom1) %>%
  dplyr::filter(WkNoStartFrom1 == 1) %>%
  dplyr::group_by(MSOAN112) %>%
  mutate(max_value = max(DarknessCrime_sum)) %>%
  dplyr::arrange(-max_value)

# MSOA = 5 also shows some unusually high values later on in the period.
spssData %>%
  dplyr::select(MSOAN112, DarknessCrime_sum, WkNoStartFrom1) %>%
  dplyr::filter(!MSOAN112 %in% c(111), WkNoStartFrom1 > 280, WkNoStartFrom1 < 300) %>%
  dplyr::group_by(MSOAN112) %>%
  mutate(max_value = max(DarknessCrime_sum)) %>%
  dplyr::arrange(-max_value)


# I will remove it from the data set in plots.
# ----


################################################################################
# Plot base:
# In the following plot, I randomly sample ten MSOAs and show their sum of crimes
# in darkness hours, over the data collection period. This plot will form the
# base of all future plots, which will overlay the expectations from increasingly
# complicated statistical models of the sum of crimes in darkness hours.
################################################################################
# ----
# Make additional dataset.
random_selection_of_MSOAs <-
  spssData %>%
   dplyr::select(MSOAN112) %>%
   unique() %>%
   dplyr::sample_n(10) %>%
   unlist() %>%
   as.numeric()
# Plot.
(p_base <-
  spssData %>%
    dplyr::select(DarknessCrime_sum, MSOAN112, WkNoStartFrom1) %>%
    dplyr::filter(!MSOAN112 %in% c(5, 111)) %>%
    dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs) %>%
    ggplot() +
    geom_line(aes(x = WkNoStartFrom1, y = DarknessCrime_sum), color = "darkgrey") +
    facet_wrap(vars(MSOAN112), ncol = 5) +
    labs(title = 'Number of crimes in the hours of darkness for a random sample of 10 MSOAs.',
         x = 'Week number', y = 'Weekly count of crimes in darkness hours')
)
# ----


################################################################################
# Plot 1:
# Data structure:
# - Grand mean
#
# The most-basic model is the arithmetic average sum of crimes in darkness hours, 
# In the following figure, I will overlay this value over each plot. I will also
# note the fit statistic, AIC.
################################################################################
# ----
# Make additional dataset.
grand_mean <-
  spssData %>% dplyr::select(DarknessCrime_sum, MSOAN112) %>%
  dplyr::group_by(MSOAN112) %>%
  dplyr::summarise(msoa_mean = mean(DarknessCrime_sum)) %>%
  dplyr::select(msoa_mean) %>% as.matrix() %>% mean()
# Fit model.
if (!exists("mod_GLMR_1"))
{
  mod_GLMR_1 <-
    glm(DarknessCrime_sum ~ 1,
        family = "poisson",
        data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
pred_GLMR_1 <-
  mod_GLMR_1$coefficients %>% exp()

# Plot.
(p1 <-
    p_base +
    geom_hline(yintercept = pred_GLMR_1, linewidth = 2) +
    labs(subtitle = paste0('Thick black line is the grand mean.\n',
                           'AIC: ', round(mod_GLMR_1$aic),
                           '\nBIC: ', round(BIC(mod_GLMR_1)) )
         )
)
# ----


################################################################################
# Plot 2:
# Data structure:
# - Grand mean
# - MSOA means, i.e. random intercept for MSOA.
#
# The next most-basic model respects the nested structure of the crimes within
# MSOAs, by fitting a random intercept for each MSOA.
################################################################################
# ----
# Make additional dataset.
# ## Poisson regression with:
# ## 1. Random intercept for `MSOAN112`, to respect nested data structure.
if (!exists(mod_GLMR_2))
{
  mod_GLMR_2 <-
    lme4::glmer(
      DarknessCrime_sum ~
        1 +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
pred_GLMR_2 <-
  predict(mod_GLMR_2) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Plot
(p2 <-
    p_base +
    geom_line(data = pred_GLMR_2, aes(x = WkNoStartFrom1, y = .), linewidth = 2) +
    labs(subtitle =
           paste0('Thick black line is the MSOA-specific mean.\n',
                  'Random effect: MSOA. \n',
                  'AIC: ', round(AIC(mod_GLMR_2)), ' (previous AIC: ', round(AIC(mod_GLMR_1)), ')',
                  '\nBIC: ', round(BIC(mod_GLMR_2)), ' (previous BIC: ', round(BIC(mod_GLMR_1)), ')'
                  )
         )
    )
# ----


################################################################################
# Plot 3:
# Data structure:
# - Grand mean
# - MSOA means, i.e. random intercept for MSOA.
# - overdispersion.
#
# We note that there is considerable overdispersion, as evidenced by the variance
# of the count variate being multiple-times larger than the arithmetic mean. We,
# therefore, attempt to account for this overdispersion by taking an
# "observation-level random effect" (ORLE) approach.
# (see https://peerj.com/articles/1114.pdf).
################################################################################
# ----
# Make additional dataset.
# ## Poisson regression with:
# ## 1. Random intercept for `MSOAN112`, to respect nested data structure.
# ## 2. Random intrecept for 'CaseID', to account for overdispersion.
if (!exists(mod_GLMR_3))
{
  mod_GLMR_3 <-
    lme4::glmer(
      DarknessCrime_sum ~
        1 +
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
pred_GLMR_3 <-
  predict(mod_GLMR_3) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Plot
(p3 <-
    p_base +
    geom_line(data = pred_GLMR_3, aes(x = WkNoStartFrom1, y = .)) +
    labs(subtitle =
           paste0('Thick black line is the MSOA-specific mean.\n',
                  'Random effect: MSOA, CaseID. \n',
                  'AIC: ', round(AIC(mod_GLMR_3)), ' (previous AIC: ', round(AIC(mod_GLMR_2)), ')',
                  '\nBIC: ', round(BIC(mod_GLMR_3)), ' (previous BIC: ', round(BIC(mod_GLMR_2)), ')'
           )
    )
)
# ----


################################################################################
# Plot 4:
# Data structure:
# - Grand mean
# - MSOA means, i.e. random intercept for MSOA.
# - overdispersion.
# Covariates:
# - Month
#
# The next model include covariates for each month, in an attempt to represent
# seasonality.
################################################################################
# ----
# Make additional dataset.
# ## Poisson regression with:
# ## 1. Random intercept for `MSOAN112`, to respect nested data structure.
# ## 2. Random intrecept for 'CaseID', to account for overdispersion.
if (!exists(mod_GLMR_4))
{
  mod_GLMR_4 <-
    lme4::glmer(
      DarknessCrime_sum ~
        1 +
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
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
pred_GLMR_4 <-
  predict(mod_GLMR_4) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)

# Plot.
(p4 <-
    p_base + 
    geom_line(data = pred_GLMR_4, aes(x = WkNoStartFrom1, y = ., group = MSOAN112)) +
    labs(subtitle =
           paste0('Thick black line is the MSOA-specific expected sum.\n',
                  'Random effect: MSOA, CaseID. | Fixed effect: Month \n',
                  'AIC: ', round(AIC(mod_GLMR_4)), ' (previous AIC: ', round(AIC(mod_GLMR_3)), ')',
                  '\nBIC: ', round(BIC(mod_GLMR_4)), ' (previous BIC: ', round(BIC(mod_GLMR_3)), ')'
           )
    )
)
# ----


################################################################################
# Plot 5:
# Data structure:
# - Grand mean
# - MSOA means, i.e. random intercept for MSOA.
# - overdispersion
# Covariates:
# - Month
# - Holidays
#
# The next model includes covariates for holiday periods, which might exacerbate 
# domestic crimes while owners are away.
################################################################################
# ----
# Make additional dataset.
# ## Poisson regression with:
# ## 1. Random intercept for `MSOAN112`, to respect nested data structure.
# ## 2. Random intrecept for 'CaseID', to account for overdispersion.
if (!exists(mod_GLMR_5))
{
  mod_GLMR_5 <-
    lme4::glmer(
      DarknessCrime_sum ~
        1 +
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
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
pred_GLMR_5 <-
  predict(mod_GLMR_5) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)

# Plot.
(p5 <-
    p_base + 
    geom_line(data = pred_GLMR_5, aes(x = WkNoStartFrom1, y = ., group = MSOAN112)) +
    labs(subtitle =
           paste0('Thick black line is the MSOA-specific expected sum.\n',
                  'Random effect: MSOA, CaseID. | Fixed effect: Month, Holidays. \n',
                  'AIC: ', round(AIC(mod_GLMR_5)), ' (previous AIC: ', round(AIC(mod_GLMR_4)), ')',
                  '\nBIC: ', round(BIC(mod_GLMR_5)), ' (previous BIC: ', round(BIC(mod_GLMR_4)), ')'
           )
    )
)
# ----


################################################################################
# Plot 6:
# Data structure:
# - Grand mean
# - MSOA means, i.e. random intercept for MSOA.
# - overdispersion.
# Covariates:
# - Month
# - Holidays
# - Dark/Light difference
#
# The next model includes a covariate to represent the changing proportion of
# daylight and darkness, over the year. Specifically, it is the difference
# between the natural logarithm of a given weeks' duration of darkness as a
# fraction of 24 hrs, and the natural logarithm of a given weeks' duration of
# daylight as a fraction of 24 hrs).
################################################################################
# ----
# Make additional dataset.
# ## Poisson regression with:
# ## 1. Random intercept for `MSOAN112`, to respect nested data structure.
# ## 2. Random intrecept for 'CaseID', to account for overdispersion.
if (!exists(mod_GLMR_6))
{
  mod_GLMR_6 <-
    lme4::glmer(
      DarknessCrime_sum ~
        1 +
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
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
pred_GLMR_6 <-
  predict(mod_GLMR_6) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)

# Plot.
(p6 <-
    p_base + 
    geom_line(data = pred_GLMR_6, aes(x = WkNoStartFrom1, y = ., group = MSOAN112)) +
    labs(subtitle =
           paste0('Thick black line is the MSOA-specific expected sum.\n',
                  'Random effect: MSOA, CaseID. | Fixed effect: Month, Holidays, Dark/Light difference. \n',
                  'AIC: ', round(AIC(mod_GLMR_6)), ' (previous AIC: ', round(AIC(mod_GLMR_5)), ')',
                  '\nBIC: ', round(BIC(mod_GLMR_6)), ' (previous BIC: ', round(BIC(mod_GLMR_5)), ')'
           )
    )
)
# ----


# These models are just so off; such bad fits. It will not be insightful to 
# evaluate the contribution of a covariate.


