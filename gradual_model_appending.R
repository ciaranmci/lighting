# The script is for producing timeline plots of both day and darkness crimes
# over time, and overlaying the expectations of increasingly-complicated
# Poisson regression GLMs.




### A NEW STORYLINE ###
#
# My new storyline is inspired by ASBATES's blogs on the use of GAMs for time
# series.
# I will start by using a lme4::glmer() model to fit `DarknessCrime_sum` with
# the main lamp covariate, DiffLfromMnBy100.
# I will interpret the coefficient but note the poor fit in the plot. I will
# then explain that the subsequent models represent progressive attempts to
# improve fit. We will watch the coefficient for DiffLfromMnBy100 as the models
# progress.
#
# 2. The second model will address the concern of overdispersion that is evident
#    from the variance of the count variate is many times the mean.
# 3. The third model will respect the data structure by including a random
#    intercept for `MSOAN112`.
# 4. The fourth model will check to see if a random slope for MSOAN112 also
#    improves fit.
# 5. The fifth model will use gamm4::gamm4() to include a joint spline for year
#    and month to account for a creeping annual trend and seasonality. I'll need
#    to create new variables for this. Both will have to be monotonically-increasing
#    integer variables, over the entire study period. Note that when predicting
#    using gamm4::gamm4(), you have to predict from the .$gam object in the list.
# 6. The sixth model will add a binary fixed effect variable representing
#    holidays. Unlike Paul's setup, I'll combine all binary holiday variables
#    into one sparse variable with 1s when the crime occurred during a holiday
#    weekend.
#
# At this point, the model will have accounted for overdispersion, respected the
# nested structure of the data, smoothly accounted any possible annual creep,
# smoothly accounted for seasonality, and will have accounted for holiday spikes.
# The last variable to include is `DifferenceInOffsets`, which represents the 
# changing proportion of light and dark.
#
# 7. The seventh model will include `DifferenceInOffsets` as a spline.
# 
# If model fit is still bad, I will include lag variables up to order four, which
# represents a month worth of history. Thus:
# 
# 8. The eighth model will include four lags of `DarknessCrime_sum`, lagging by
#    one, two, three, and four weeks.
#



###############
# Requisites. #
###############
# ----
if (!("pacman" %in% installed.packages()[,"Package"])) install.packages("pacman")
library(pacman)
pacman::p_load(tidyverse, haven, R2MLwiN, splines)
# ----


##############
# Load data. #
##############
# ----
if (!exists("spssData"))
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


# Add other variables.
if(!"month_idx" %in% colnames(spssData))
{
  spssData$month_idx <-
    spssData %>%
    dplyr::select(starts_with("MonthMidWk_")) %>%
    tidyr::unite(col = unitedMonths, sep = "") %>%
    dplyr::mutate(
      month_idx = dplyr::case_when(
        . == '00000000000' ~ '1',
        . == '10000000000' ~ '2',
        . == '01000000000' ~ '3',
        . == '00100000000' ~ '4',
        . == '00010000000' ~ '5',
        . == '00001000000' ~ '6',
        . == '00000100000' ~ '7',
        . == '00000010000' ~ '8',
        . == '00000001000' ~ '9',
        . == '00000000100' ~ '10',
        . == '00000000010' ~ '11',
        . == '00000000001' ~ '12',
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(month_idx) %>%
    unlist() %>%
    as.numeric()
} 
if(!"year_idx" %in% colnames(spssData))
{
  spssData <-
    spssData %>%
    group_by(MSOAN112) %>% 
    dplyr::select(CaseID, month_idx) %>%
    mutate(increment = 1 + cumsum( (month_idx == 12) & (lead(month_idx) == 1) )) %>% 
    dplyr::mutate(
      year_idx = dplyr::lag(increment, n = 1L)
    ) %>% 
    replace_na(list(year_idx=1)) %>%
    dplyr::ungroup() %>%
    dplyr::select(CaseID, year_idx) %>%  
    dplyr::right_join(spssData, by = 'CaseID')
}
if(!"T_Dec_gm" %in% colnames(spssData))
{
  spssData <-
    spssData %>%
    dplyr::mutate(
      T_Dec_gm = scale(T_Dec, scale = F),
      T_Dec_gm_pw2 = scale(T_Dec, scale = F)^2,
      T_Dec_gm_pw3 = scale(T_Dec, scale = F)^3,
      T_Dec_gm_pw4 = scale(T_Dec, scale = F)^4,
      T_Dec_gm_pw5 = scale(T_Dec, scale = F)^5
    )
}
if(!"PropDarkCrime_lag1" %in% colnames(spssData))
{
  spssData <-
    spssData %>%
    dplyr::group_by(MSOAN112) %>%
    dplyr::mutate(
      PropDarkCrime_lag1 = dplyr::lag(PropDarkCrime, n = 1L),
      PropDarkCrime_lag2 = dplyr::lag(PropDarkCrime, n = 2L)
    )
}
# ----


############################################
# Randomly select some MSOAs for plotting. #
############################################
# ----
random_selection_of_MSOAs <-
  spssData %>%
  dplyr::select(MSOAN112) %>%
  unique() %>%
  dplyr::sample_n(10) %>%
  unlist() %>%
  as.numeric()
# ----

#########################
# Define some constants #
#########################
# ----
text_str_prop <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the proportion of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
text_str_count <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
text_lrtest <- 
  paste0("A log-likelihood ratio test suggests that the inclusion of the ",
         "covariate is statistically significant at any reasonable value:")
text_plot_bad_fit <-
  "...plotting the expected values from the model show a very poor fit:"
text_plot_good_fit <-
  "...plotting the expected values from the model show a good fit:"
# ----


################################################################################
#####                                                                      #####
#####                        VARIATE: PropDarkCrime                        #####
#####                                                                      #####
################################################################################

################################################################################
# Plot base:
# In the following plot, I randomly sample ten MSOAs and show their proportion
# crimes in darkness, over the data collection period. This plot will form the
# base of all future plots, which will overlay the expectations from increasingly
# complicated statistical models of the sum of crimes in darkness hours.
#
# I also fit the base Binomial model for the proportion of crimes in darkness,
# fit using a basic GLM, stas::glm().
################################################################################
# ----
# Fit null model.
if (!exists("prop_mod_GLMR_0"))
{
  prop_mod_GLMR_0 <-
    glm(
      cbind(DarknessCrime_sum, SumDarkAndDaylight) ~ 1,
      family = binomial("logit"),
      data = spssData)
}
prop_pred_GLMR_0 <-
  predict(prop_mod_GLMR_0) %>% exp() %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>% dplyr::select(WkNoStartFrom1, MSOAN112) ) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Plot.
(p_base_prop <-
    spssData %>%
    dplyr::select(PropDarkCrime, MSOAN112, WkNoStartFrom1) %>%
    dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs) %>%
    ggplot() +
    geom_line(aes(x = WkNoStartFrom1, y = PropDarkCrime), color = "darkgrey") +
    facet_wrap(vars(MSOAN112), ncol = 5) +
    labs(title = 'Proportion of crimes in the hours of darkness for a random sample of 10 MSOAs.',
         x = 'Week number', y = 'Weekly proportion of crimes in darkness hours')
)
# ----
p_base_prop + geom_line(data = prop_pred_GLMR_0, aes(x = WkNoStartFrom1, y = .), color = "black")
prop_mod_GLMR_0 %>% summary()

################################################################################
# Model 1:
# Data structure:
# - Grand mean
#
# The most-basic model is the arithmetic average proportion of crimes in darkness
# In the following figure, I will overlay this value over each plot. I will also
# note the fit statistic, AIC.
################################################################################
# ----
# Fit model.
if (!exists("prop_mod_GLMR_1"))
{
  prop_mod_GLMR_1 <-
    glm(
      cbind(DarknessCrime_sum, SumDarkAndDaylight) ~ 
        WithinMSOALightingTerm +
        BetweenAreaLightingTerm,
      family = binomial("logit"),
      data = spssData)
}
prop_pred_GLMR_1 <-
  predict(prop_mod_GLMR_1) %>% exp() %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112)) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)

# Make plot.
p1_prop <-
  p_base_prop +
  geom_line(data = prop_pred_GLMR_1, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle = paste0('Black line is the model including number of lights.\n',
                         'AIC: ', round(prop_mod_GLMR_1$aic),
                         '\nBIC: ', round(BIC(prop_mod_GLMR_1)) )
  )

# Interpret coefficient and fit.
print_vals <-
  round(prop_mod_GLMR_1$coefficients %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_str_prop, print_vals))
print(text_lrtest)
lmtest::lrtest(prop_mod_GLMR_1)
print(text_plot_bad_fit)
p1_prop # Show plot.
print(
  paste0("The subsequent models represent progressive attempts to improve fit.",
         " I will note the coefficient for `WithinMSOALightingTerm` as the models progress.")
)
# ----
prop_mod_GLMR_1 %>% summary()



################################################################################
# Model 2:
# Data structure:
# - Grand mean
# - Random intercept for MSOA
#
# Respect the data structure by including a random intercept for `MSOAN112`.
################################################################################
# ----
# Fit model.
if (!exists("prop_mod_GLMR_2"))
{
  prop_mod_GLMR_2 <-
    lme4::glmer(
      cbind(DarknessCrime_sum, SumDarkAndDaylight) ~ 
        WithinMSOALightingTerm +
        BetweenAreaLightingTerm +
        (1 | MSOAN112),
      family = binomial("logit"),
      data = spssData)
}
prop_pred_GLMR_2 <-
  predict(prop_mod_GLMR_2) %>% exp() %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112)) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)

# Make plot.
p2_prop <-
  p_base_prop +
  geom_line(data = prop_pred_GLMR_2, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle = paste0('Black line is the model including number of lights.\n',
                         'AIC: ', round(AIC(prop_mod_GLMR_2)),
                         '\nBIC: ', round(BIC(prop_mod_GLMR_2)) )
  )

# Interpret coefficient and fit.
print_vals <-
  round(prop_mod_GLMR_2$coefficients %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_str_prop, print_vals))
print(text_lrtest)
lmtest::lrtest(prop_mod_GLMR_1, prop_mod_GLMR_2)
print(text_plot_bad_fit)
p2_prop # Show plot.
# ----
prop_mod_GLMR_2 %>% summary()


################################################################################
# Model 3:
# Data structure:
# - Grand mean
# - Random intercept for MSOA
#
# Respect the data structure by including a random intercept for `MSOAN112`.
################################################################################
# ----
# Fit model.
if (!exists("prop_mod_GLMR_3"))
{
  prop_mod_GLMR_3 <-
    lme4::glmer(
      cbind(DarknessCrime_sum, SumDarkAndDaylight) ~ 
        WithinMSOALightingTerm +
        BetweenAreaLightingTerm +
        (1 | MSOAN112),
      family = binomial("logit"),
      offset = DifferenceInOffsets,
      data = spssData)
}
prop_pred_GLMR_3 <-
  predict(prop_mod_GLMR_3) %>% exp() %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112)) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)

# Make plot.
p3_prop <-
  p_base_prop +
  geom_line(data = prop_pred_GLMR_3, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle = paste0('Black line is the model including number of lights.\n',
                         'AIC: ', round(AIC(prop_mod_GLMR_3)),
                         '\nBIC: ', round(BIC(prop_mod_GLMR_3)) )
  )

# Interpret coefficient and fit.
print_vals <-
  round(prop_mod_GLMR_3@beta %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = text_str_prop, print_vals))
print(text_lrtest)
lmtest::lrtest(prop_mod_GLMR_2, prop_mod_GLMR_3)
print(text_plot_good_fit)
p3_prop # Show plot.
# ----
prop_mod_GLMR_3 %>% summary()

################################################################################
# Model 3:
# Data structure:
# - Grand mean
# - Random intercept for MSOA
# - Offset for darkness hours
#
# Add offset indicating the difference in dark and daylight proportions of a day.
################################################################################
# ----
# Fit model.
if (!exists("prop_mod_GLMR_3"))
{
  prop_mod_GLMR_3 <-
    lme4::glmer(
      cbind(DarknessCrime_sum, SumDarkAndDaylight) ~ 
        WithinMSOALightingTerm +
        BetweenAreaLightingTerm +
        (1 | MSOAN112),
      family = binomial("logit"),
      offset = DifferenceInOffsets,
      data = spssData)
}
prop_pred_GLMR_3 <-
  predict(prop_mod_GLMR_3) %>% exp() %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112)) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)

# Make plot.
p3_prop <-
  p_base_prop +
  geom_line(data = prop_pred_GLMR_3, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle = paste0('Black line is the model including number of lights.\n',
                         'AIC: ', round(AIC(prop_mod_GLMR_3)),
                         '\nBIC: ', round(BIC(prop_mod_GLMR_3)) )
  )

# Interpret coefficient and fit.
print_vals <-
  round(prop_mod_GLMR_3@beta %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = text_str_prop, print_vals))
print(text_lrtest)
lmtest::lrtest(prop_mod_GLMR_2, prop_mod_GLMR_3)
print(text_plot_good_fit)
p3_prop # Show plot.
# ----
prop_mod_GLMR_3 %>% summary()


################################################################################
# Model 4:
# Data structure:
# - Grand mean
# - Random intercept for MSOA
# - Offset for darkness hours
#
# Fixed effects:
# - Indicators of holiday weeks.
################################################################################
# ----
# Fit model.
if (!exists("prop_mod_GLMR_4"))
{
  prop_mod_GLMR_4 <-
    lme4::glmer(
      cbind(DarknessCrime_sum, SumDarkAndDaylight) ~ 
        WithinMSOALightingTerm +
        BetweenAreaLightingTerm +
        NewYrWk +
        GoodFriWk +
        EasterWk +
        MayDayWk +
        SpringBankWk +
        SummerBankWk +
        XmasWk +
        (1 | MSOAN112),
      family = binomial("logit"),
      offset = DifferenceInOffsets,
      data = spssData)
}
prop_pred_GLMR_4 <-
  predict(prop_mod_GLMR_4) %>% exp() %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112)) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)

# Make plot.
p4_prop <-
  p_base_prop +
  geom_line(data = prop_pred_GLMR_4, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle = paste0('Black line is the model including number of lights.\n',
                         'AIC: ', round(AIC(prop_mod_GLMR_4)),
                         '\nBIC: ', round(BIC(prop_mod_GLMR_4)) )
  )

# Interpret coefficient and fit.
print_vals <-
  round(prop_mod_GLMR_4@beta %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = text_str_prop, print_vals))
print(text_lrtest)
lmtest::lrtest(prop_mod_GLMR_3, prop_mod_GLMR_4)
print(text_plot_good_fit)
p4_prop # Show plot.
# ----
prop_mod_GLMR_4 %>% summary()

################################################################################
# Model 5:
# Data structure:
# - Grand mean
# - Random intercept for MSOA
# - Offset for darkness hours
#
# Fixed effects:
# - Indicators of holiday weeks.
# - Indicators for months.
################################################################################
# ----
# Fit model.
if (!exists("prop_mod_GLMR_5"))
{
  prop_mod_GLMR_5 <-
    lme4::glmer(
      cbind(DarknessCrime_sum, SumDarkAndDaylight) ~ 
        WithinMSOALightingTerm +
        BetweenAreaLightingTerm +
        NewYrWk +
        GoodFriWk +
        EasterWk +
        MayDayWk +
        SpringBankWk +
        SummerBankWk +
        XmasWk +
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
        (1 | MSOAN112),
      family = binomial("logit"),
      offset = DifferenceInOffsets,
      data = spssData)
}
prop_pred_GLMR_5 <-
  predict(prop_mod_GLMR_5) %>% exp() %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112)) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)

# Make plot.
p5_prop <-
  p_base_prop +
  geom_line(data = prop_pred_GLMR_5, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle = paste0('Black line is the model including number of lights.\n',
                         'AIC: ', round(AIC(prop_mod_GLMR_5)),
                         '\nBIC: ', round(BIC(prop_mod_GLMR_5)) )
  )

# Interpret coefficient and fit.
print_vals <-
  round(prop_mod_GLMR_5@beta %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = text_str_prop, print_vals))
print(text_lrtest)
lmtest::lrtest(prop_mod_GLMR_4, prop_mod_GLMR_5)
print(text_plot_good_fit)
p5_prop # Show plot.
# ----
prop_mod_GLMR_5 %>% summary()

################################################################################
# Model 6:
# Data structure:
# - Grand mean
# - Random intercept for MSOA
# - Offset for darkness hours
#
# Fixed effects:
# - Indicators of holiday weeks.
# - Indicators for months.
# - Lag terms for the outcome.
#
################################################################################
# ----
# Fit model.
if (!exists("prop_mod_GLMR_6"))
{
  prop_mod_GLMR_6 <-
    lme4::glmer(
      cbind(DarknessCrime_sum, SumDarkAndDaylight) ~ 
        WithinMSOALightingTerm +
        BetweenAreaLightingTerm +
        NewYrWk +
        GoodFriWk +
        EasterWk +
        MayDayWk +
        SpringBankWk +
        SummerBankWk +
        XmasWk +
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
        PropDarkCrime_lag1 +
        PropDarkCrime_lag2 +
        (1 | MSOAN112),
      family = binomial("logit"),
      offset = DifferenceInOffsets,
      data = spssData)
}
prop_pred_GLMR_6 <-
  predict(prop_mod_GLMR_6, newdata = spssData) %>%
  exp() %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>% dplyr::select(WkNoStartFrom1, MSOAN112)) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)

# Make plot.
p6_prop <-
  p_base_prop +
  geom_line(data = prop_pred_GLMR_6 %>% dplyr::filter(!is.na(.)),
            aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle = paste0('Black line is the model including number of lights.\n',
                         'AIC: ', round(AIC(prop_mod_GLMR_6)),
                         '\nBIC: ', round(BIC(prop_mod_GLMR_6)) )
  )

# Interpret coefficient and fit.
print_vals <-
  round(prop_mod_GLMR_6@beta %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = text_str_prop, print_vals))
print(text_lrtest)
lmtest::lrtest(prop_mod_GLMR_5, prop_mod_GLMR_6)
print("I can't figure out why the plot is shifted")
p6_prop # Show plot.
# ----
prop_mod_GLMR_6 %>% summary()

################################################################################
# Model 7:
# Data structure:
# - Grand mean
# - Random intercept for MSOA
# - Offset for darkness hours
#
# Fixed effects:
# - Indicators of holiday weeks.
# - Indicators for months.
# - Lag terms for the outcome.
# - Polynomials for decimal time throughout the study period.
#
# Random effects:
# - Polynomials for decimal time throughout the study period.
#
################################################################################
# ----
# Fit model.
if (!exists("prop_mod_GLMR_7"))
{
  prop_mod_GLMR_7 <-
    lme4::glmer(
      cbind(DarknessCrime_sum, SumDarkAndDaylight) ~ 
        WithinMSOALightingTerm +
        BetweenAreaLightingTerm +
        NewYrWk +
        GoodFriWk +
        EasterWk +
        MayDayWk +
        SpringBankWk +
        SummerBankWk +
        XmasWk +
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
        PropDarkCrime_lag1 +
        PropDarkCrime_lag2 +
        T_Dec_gm_pw4 +
        T_Dec_gm_pw5 +
        (T_Dec + T_Dec_gm_pw2 + T_Dec_gm_pw3 | MSOAN112),
      family = binomial("logit"),
      offset = DifferenceInOffsets,
      data = spssData)
}
prop_pred_GLMR_7 <-
  predict(prop_mod_GLMR_7,newdata = spssData) %>%
  exp() %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>% dplyr::select(WkNoStartFrom1, MSOAN112)) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)

# Make plot.
p7_prop <-
  p_base_prop +
  geom_line(data = prop_pred_GLMR_7 %>% dplyr::filter(!is.na(.)),
            aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle = paste0('Black line is the model including number of lights.\n',
                         'AIC: ', round(AIC(prop_mod_GLMR_7)),
                         '\nBIC: ', round(BIC(prop_mod_GLMR_7)) )
  )

# Interpret coefficient and fit.
print_vals <-
  round(prop_mod_GLMR_7@beta %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = text_str_prop, print_vals))
print(text_lrtest)
lmtest::lrtest(prop_mod_GLMR_6, prop_mod_GLMR_7)
print("I can't figure out why the plot is shifted")
p7_prop # Show plot.
# ----
prop_mod_GLMR_7 %>% summary()


p_base_prop +
  geom_point(data =
               spssData %>%
               dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs),
             aes(x = WkNoStartFrom1, y = N_LampsChanged/1000),
             color = 'black') +
  scale_y_continuous(
    name =  "Cumulative # lamps installed",
    sec.axis = sec_axis(~.*1000))


## Beta-binomial model 
# 
# Some notes on an alternative approach using gamm4:
# https://stackoverflow.com/questions/34495473/model-selection-with-beta-and-quassi-families-using-gamm4
#
prop_mod_GLMR_7 <-
  spaMM::HLfit(
    formula = cbind(DarknessCrime_sum, SumDarkAndDaylight) ~ 
      WithinMSOALightingTerm +
      BetweenAreaLightingTerm +
      NewYrWk +
      GoodFriWk +
      EasterWk +
      MayDayWk +
      SpringBankWk +
      SummerBankWk +
      XmasWk +
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
      PropDarkCrime_lag1 +
      PropDarkCrime_lag2 +
      T_Dec_gm_pw4 +
      T_Dec_gm_pw5 +
      (T_Dec + T_Dec_gm_pw2 + T_Dec_gm_pw3 | MSOAN112) +
      offset(DifferenceInOffsets)
    ,data = spssData %>% dplyr::filter(SumDarkAndDaylight!=0)
    ,family=binomial()
    ,rand.family=Beta()
    ,HLmethod="HL(0,0)"
    #,resid.model = ~1
    #,REMLformula = NULL
    #,verbose = c(trace = FALSE)
    #,control.HLfit = list()
    #,control.glm = list()
    #,init.HLfit = list()
    #,ranFix = list()
    #,etaFix = list()
    #,prior.weights = NULL
    #,processed = NULL
    )



################################################################################
# Model 8:
# Data structure:
# - Overdispersion not modelled.
# - Random intercept for MSOA.
# Splines:
# - Spline for `T_Dec`
# Autoregression:
# - Variate lagged by one.
# Offset:
# - `DifferenceInOffsets`
# Covariates:
# - Month indicators
# - Holiday indicators
# - WithinMSOALightingTerm
# - BetweenAreaLightingTerm
################################################################################
# ----
# Fit model.
if (!exists("prop_mod_GAMM_8"))
{
  prop_mod_GAMM_8fac <-
    gamm4::gamm4(
      cbind(DarknessCrime_sum, SumDarkAndDaylight) ~
        WithinMSOALightingTerm +
        BetweenAreaLightingTerm +
        as.factor(NewYrWk) +
        as.factor(GoodFriWk) +
        as.factor(EasterWk) +
        as.factor(MayDayWk) +
        as.factor(SpringBankWk) +
        as.factor(SummerBankWk) +
        as.factor(XmasWk) +
        as.factor(MonthMidWk_2) +
        as.factor(MonthMidWk_3) +
        as.factor(MonthMidWk_4) +
        as.factor(MonthMidWk_5) +
        as.factor(MonthMidWk_6) +
        as.factor(MonthMidWk_7) +
        as.factor(MonthMidWk_8) +
        as.factor(MonthMidWk_9) +
        as.factor(MonthMidWk_10) +
        as.factor(MonthMidWk_11) +
        as.factor(MonthMidWk_12) +
        PropDarkCrime_lag1 +
        s(T_Dec) +
        offset(DifferenceInOffsets),
      random = ~ (1 | MSOAN112),
      family = "binomial",
      data = spssData %>% dplyr::filter(!is.na(PropDarkCrime_lag1))
    )
}
prop_pred_GAMM_8fac <-
  predict(prop_mod_GAMM_8fac$gam, newdata = spssData) %>%
  exp() %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>% dplyr::select(WkNoStartFrom1, MSOAN112)) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p8fac_prop <-
  p_base_prop +
  geom_line(data = prop_pred_GAMM_8fac %>% dplyr::filter(!is.na(.)),
            aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle = paste0('Black line is the model including number of lights.\n',
                         'AIC: ', round(AIC(prop_mod_GAMM_8fac$mer)),
                         '\nBIC: ', round(BIC(prop_mod_GAMM_8fac$mer)) )
  )
# Interpret coefficient and fit.
print_vals <-
  round(prop_mod_GAMM_8$mer@beta[1:3] %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = text_str_prop, print_vals))
p8fac_prop # Show plot.
# ----
prop_mod_GAMM_8$gam %>% summary()
prop_mod_GAMM_8$mer %>% summary()
termplot(prop_mod_GAMM_8fac$gam, all.terms = T)


################################################################################
# Model 9:
# Data structure:
# - Overdispersion not modelled.
# - Random intercept for MSOA.
# Splines:
# - None
# Autoregression:
# - Variate lagged by one.
# Offset:
# - `DifferenceInOffsets`
# Covariates:
# - Month indicators
# - Holiday indicators
# - WithinMSOALightingTerm
# - BetweenAreaLightingTerm
################################################################################
# ----
# Fit model.
if (!exists("prop_mod_GAMM_9"))
{
  prop_mod_GAMM_9 <-
    gamm4::gamm4(
      cbind(DarknessCrime_sum, SumDarkAndDaylight) ~
        WithinMSOALightingTerm +
        BetweenAreaLightingTerm +
        NewYrWk +
        GoodFriWk +
        EasterWk +
        MayDayWk +
        SpringBankWk +
        SummerBankWk +
        XmasWk +
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
        PropDarkCrime_lag1 +
        offset(DifferenceInOffsets),
      random = ~ (1 | MSOAN112),
      family = "binomial",
      data = spssData %>% dplyr::filter(!is.na(PropDarkCrime_lag1))
    )
}
prop_pred_GAMM_9 <-
  predict(prop_mod_GAMM_8$gam, newdata = spssData) %>%
  exp() %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>% dplyr::select(WkNoStartFrom1, MSOAN112)) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p9_prop <-
  p_base_prop +
  geom_line(data = prop_pred_GAMM_9 %>% dplyr::filter(!is.na(.)),
            aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle = paste0('Black line is the model including number of lights.\n',
                         'AIC: ', round(AIC(prop_mod_GAMM_9$mer)),
                         '\nBIC: ', round(BIC(prop_mod_GAMM_9$mer)) )
  )
# Interpret coefficient and fit.
print_vals <-
  round(prop_mod_GAMM_9$mer@beta[1:3] %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = text_str_prop, print_vals))
p9_prop # Show plot.
# ----
prop_mod_GAMM_9$gam %>% summary()
prop_mod_GAMM_9$mer %>% summary()




################################################################################
# Model 10:
# Data structure:
# - Overdispersion modelled using a beta-binomial link function, in glmmADMB.
# - Random intercept for MSOA.
# Splines:
# - Spline for `T_Dec`, setting the degrees of freedom to 5 because the freely-
#   estimated degrees of freedom using gamm4::gamm4() was 5.242.
# Autoregression:
# - Variate lagged by one.
# Offset:
# - `DifferenceInOffsets`
# Covariates:
# - Month indicators
# - Holiday indicators
# - WithinMSOALightingTerm
# - BetweenAreaLightingTerm
################################################################################
# ----
# Fit model.
if (!exists("prop_mod_GAMM_10") )
{
  prop_mod_GAMM_10 <-
    glmmADMB::glmmadmb(
      cbind(DarknessCrime_sum, SumDarkAndDaylight) ~
        WithinMSOALightingTerm +
        BetweenAreaLightingTerm +
        as.factor(NewYrWk) +
        as.factor(GoodFriWk) +
        as.factor(EasterWk) +
        as.factor(MayDayWk) +
        as.factor(SpringBankWk) +
        as.factor(SummerBankWk) +
        as.factor(XmasWk) +
        as.factor(MonthMidWk_2) +
        as.factor(MonthMidWk_3) +
        as.factor(MonthMidWk_4) +
        as.factor(MonthMidWk_5) +
        as.factor(MonthMidWk_6) +
        as.factor(MonthMidWk_7) +
        as.factor(MonthMidWk_8) +
        as.factor(MonthMidWk_9) +
        as.factor(MonthMidWk_10) +
        as.factor(MonthMidWk_11) +
        as.factor(MonthMidWk_12) +
        PropDarkCrime_lag1 +
        splines::ns(T_Dec, df = 5) +
        offset(DifferenceInOffsets) +
        (1 | as.factor(MSOAN112) ),
      family = "beta",
      link = "logit",
      data = spssData %>% dplyr::filter(!is.na(PropDarkCrime_lag1) )
    )
}
prop_pred_GAMM_10 <-
  predict(prop_mod_GAMM_10$gam, newdata = spssData) %>%
  exp() %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>% dplyr::select(WkNoStartFrom1, MSOAN112)) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p10_prop <-
  p_base_prop +
  geom_line(data = prop_pred_GAMM_10 %>% dplyr::filter(!is.na(.)),
            aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle = paste0('Black line is the model including number of lights.\n',
                         'AIC: ', round(AIC(prop_mod_GAMM_10$mer)),
                         '\nBIC: ', round(BIC(prop_mod_GAMM_10$mer)) )
  )
# Interpret coefficient and fit.
print_vals <-
  round(prop_mod_GAMM_10$mer@beta[1:3] %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = text_str_prop, print_vals))
p10_prop # Show plot.
# ----
prop_mod_GAMM_10$gam %>% summary()
prop_mod_GAMM_10$mer %>% summary()
termplot(prop_mod_GAMM_10$gam, all.terms = T)





################################################################################
#####                                                                      #####
#####                    VARIATE: DarknessCrime_sum                        #####
#####                                                                      #####
################################################################################

################################################################################
# Plot base:
# In the following plot, I randomly sample ten MSOAs and show their sum of crimes
# in darkness hours, over the data collection period. This plot will form the
# base of all future plots, which will overlay the expectations from increasingly
# complicated statistical models of the sum of crimes in darkness hours.
#
# I also fit the base Poisson model for the count of crimes in darkness hours,
# fit using a basic GLM, stas::glm().
################################################################################
# ----
# Fit null model.
if (!exists("darknessCount_mod_GLMR_0"))
{
  darknessCount_mod_GLMR_0 <-
    glm(DarknessCrime_sum ~ 1,
        family = "poisson",
        data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
darknessCount_pred_GLMR_0 <-
  predict(darknessCount_mod_GLMR_0) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Plot.
(p_base_darkness <-
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
darknessCount_mod_GLMR_0 %>% summary()


################################################################################
# Model 1:
# Data structure:
# - Grand mean
#
# The most-basic model is the arithmetic average sum of crimes in darkness hours, 
# In the following figure, I will overlay this value over each plot. I will also
# note the fit statistic, AIC.
################################################################################
# ----
# Fit model.
if (!exists("darknessCount_mod_GLMR_1"))
{
  darknessCount_mod_GLMR_1 <-
    glm(DarknessCrime_sum ~ 
          1 +
          DiffLfromMnBy100,
        family = "poisson",
        data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
darknessCount_pred_GLMR_1 <-
  predict(darknessCount_mod_GLMR_1) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)

# Make plot.
p1_darkness <-
    p_base_darkness +
    geom_line(data = darknessCount_pred_GLMR_1, aes(x = WkNoStartFrom1, y = .)) +
    labs(subtitle = paste0('Black line is the model including number of lights.\n',
                           'AIC: ', round(darknessCount_mod_GLMR_1$aic),
                           '\nBIC: ', round(BIC(darknessCount_mod_GLMR_1)) )
         )

# Interpret coefficient and fit.
print_Str <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
print_vals <-
  round(darknessCount_mod_GLMR_1$coefficients %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_Str, print_vals))

print(
  paste0("A log-likelihood ratio test suggests that the inclusion of the ",
         "covariate is statistically significant at any reasonable value:")
)
lmtest::lrtest(darknessCount_mod_GLMR_1)
print(
  paste0("...plotting the expected values from the model show a very poor fit:")
)
p1_darkness # Show plot.
print(
  paste0("The subsequent models represent progressive attempts to improve fit.",
        " I will note the coefficient for `DiffLfromMnBy100` as the models progress.")
)
# ----
darknessCount_mod_GLMR_1 %>% summary()


################################################################################
# Model 2:
# Data structure:
# - Observation-level random effect (ORLE) to handle overdispersion.
#
# The second model will address the concern of overdispersion that is raised
# because the variance of the count variate is many times the mean.
################################################################################
# ----
# Fit model.
if (!exists("darknessCount_mod_GLMR_2"))
{
  darknessCount_mod_GLMR_2 <-
    lme4::glmer(
      DarknessCrime_sum ~
        1 +
        DiffLfromMnBy100 +
        (1 | CaseID),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
darknessCount_pred_GLMR_2 <-
  predict(darknessCount_mod_GLMR_2) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p2_darkness <-
    p_base_darkness +
    geom_line(data = darknessCount_pred_GLMR_2, aes(x = WkNoStartFrom1, y = .)) +
    labs(subtitle =
           paste0('Black line is the MSOA-specific mean.\n',
                  'Random effect: CaseID | Fixed effect: # of lights (centred) \n',
                  'AIC: ', round(AIC(darknessCount_mod_GLMR_2)), ' (previous AIC: ', round(AIC(darknessCount_mod_GLMR_1)), ')',
                  '\nBIC: ', round(BIC(darknessCount_mod_GLMR_2)), ' (previous BIC: ', round(BIC(darknessCount_mod_GLMR_1)), ')'
                  )
         )
  
# Interpret coefficient and fit.
print_Str <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
print_vals <-
  round(darknessCount_mod_GLMR_2@beta %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_Str, print_vals))
print(
  paste0("A log-likelihood ratio test suggests that the inclusion of the ",
         "covariate is statistically significant at any reasonable value:")
)
lmtest::lrtest(darknessCount_mod_GLMR_1,darknessCount_mod_GLMR_2)
print(
  paste0("...plotting the expected values from the model show a very poor fit:")
)
p2_darkness # Show plot.
# ----
darknessCount_mod_GLMR_2 %>% summary()


################################################################################
# Model 3:
# Data structure:
# - Observation-level random effect (ORLE) to handle overdispersion.
# - Random intercept for MSOA.
#
# The third model will respect the data structure by including a random
# intercept for `MSOAN112`. I also include a fixed-effect variable to represent
# each MSOA's count of lamps, centred to the grand mean.
################################################################################
# ----
# Fit model.
if (!exists("darknessCount_mod_GLMR_3"))
{
  darknessCount_mod_GLMR_3 <-
    lme4::glmer(
      DarknessCrime_sum ~
        1 +
        DiffLfromMnBy100 +
        DiffMSOA_MnLfromGrandMnBy100 +
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
darknessCount_pred_GLMR_3 <-
  predict(darknessCount_mod_GLMR_3) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p3_darkness <-
  p_base_darkness +
  geom_line(data = darknessCount_pred_GLMR_3, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle =
         paste0('Black line is the MSOA-specific mean.\n',
                'Random effect: CaseID, MSOA | Fixed effect: # of lights (centred) \n',
                'AIC: ', round(AIC(darknessCount_mod_GLMR_3)), ' (previous AIC: ', round(AIC(darknessCount_mod_GLMR_2)), ')',
                '\nBIC: ', round(BIC(darknessCount_mod_GLMR_3)), ' (previous BIC: ', round(BIC(darknessCount_mod_GLMR_2)), ')'
         )
  )
# Interpret coefficient and fit.
print_Str <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
print_vals <-
  round(darknessCount_mod_GLMR_3@beta %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_Str, print_vals))
print(
  paste0("A log-likelihood ratio test suggests that the inclusion of the ",
         "covariate is statistically significant at any reasonable value:")
)
lmtest::lrtest(darknessCount_mod_GLMR_2,darknessCount_mod_GLMR_3)
print(
  paste0("...plotting the expected values from the model show a very poor fit:")
)
p3_darkness # Show plot.
# ----
darknessCount_mod_GLMR_3 %>% summary()


# ********************
# **** NOT FITTED ****
# ********************
################################################################################
# Model 4:
# Data structure:
# - Observation-level random effect (ORLE) to handle overdispersion.
# - Random intercept for MSOA.
# - Random slope for MSOA.
#
# The fourth model will check to see if a random slope for MSOAN112 also
# improves fit.
################################################################################
# ----
# # Fit model.
# if (!exists("darknessCount_mod_GLMR_4"))
# {
#   darknessCount_mod_GLMR_4 <-
#     lme4::glmer(
#       DarknessCrime_sum ~
#         1 +
#         DiffLfromMnBy100 +
#         (1 | CaseID) +
#         (1 + MSOAN112 | MSOAN112),
#       family = "poisson",
#       data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
#     )
# }
# darknessCount_pred_GLMR_4 <-
#   predict(darknessCount_mod_GLMR_4) %>%
#   as.data.frame() %>%
#   dplyr::bind_cols(spssData %>%
#                      dplyr::select(WkNoStartFrom1, MSOAN112) %>%
#                      dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
#   dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# # Make plot.
# p4_darkness <-
#   p_base_darkness +
#   geom_line(data = darknessCount_pred_GLMR_4, aes(x = WkNoStartFrom1, y = .)) +
#   labs(subtitle =
#          paste0('Black line is the MSOA-specific mean.\n',
#                 'Random effect: CaseID, MSOA | Fixed effect: # of lights (centred) \n',
#                 'AIC: ', round(AIC(darknessCount_mod_GLMR_4)), ' (previous AIC: ', round(AIC(darknessCount_mod_GLMR_3)), ')',
#                 '\nBIC: ', round(BIC(darknessCount_mod_GLMR_4)), ' (previous BIC: ', round(BIC(darknessCount_mod_GLMR_3)), ')'
#          )
#   )
# # Interpret coefficient and fit.
# print_Str <- 
#   paste0("The exponentiated coefficient is %g, which we interpret as saying ",
#          "that the count of crimes committed in the hours of darkness are ",
#          "expected to be %g-times greater with the addition of one light.")
# print_vals <-
#   round(darknessCount_mod_GLMR_4@beta %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
# do.call(sprintf, c(fmt = print_Str, print_vals))
# print(
#   paste0("A log-likelihood ratio test suggests that the inclusion of the ",
#          "covariate is statistically significant at any reasonable value:")
# )
# lmtest::lrtest(darknessCount_mod_GLMR_3,darknessCount_mod_GLMR_4)
# print(
#   paste0("...plotting the expected values from the model show a very poor fit:")
# )
# p4_darkness # Show plot.
# ----
#darknessCount_mod_GLMR_4 %>% summary()


################################################################################
# Model 5:
# Data structure:
# - Observation-level random effect (ORLE) to handle overdispersion.
# - Random intercept for MSOA.
# Splines:
# - Year-Month manifold.
#
# The fifth model will use gamm4::gamm4() to include a joint spline for year
# and month to account for a creeping annual trend and seasonality. I'll need
# to create new variables for this. Both will have to be monotonically-increasing
# integer variables, over the entire study period. Note that when predicting
# using gamm4::gamm4(), you have to predict from the .$gam object in the list.
################################################################################
# ----
# Make additional variables.
if(!"month_idx" %in% colnames(spssData))
{
  spssData$month_idx <-
    spssData %>%
    dplyr::select(starts_with("MonthMidWk_")) %>%
    tidyr::unite(col = unitedMonths, sep = "") %>%
    dplyr::mutate(
      month_idx = dplyr::case_when(
        . == '00000000000' ~ '1',
        . == '10000000000' ~ '2',
        . == '01000000000' ~ '3',
        . == '00100000000' ~ '4',
        . == '00010000000' ~ '5',
        . == '00001000000' ~ '6',
        . == '00000100000' ~ '7',
        . == '00000010000' ~ '8',
        . == '00000001000' ~ '9',
        . == '00000000100' ~ '10',
        . == '00000000010' ~ '11',
        . == '00000000001' ~ '12',
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(month_idx) %>%
    unlist() %>%
    as.numeric()
}

if(!"year_idx" %in% colnames(spssData))
{
  spssData <-
    spssData %>%
    group_by(MSOAN112) %>% 
    dplyr::select(CaseID, month_idx) %>%
    mutate(increment = 1 + cumsum( (month_idx == 12) & (lead(month_idx) == 1) )) %>% 
    dplyr::mutate(
      year_idx = dplyr::lag(increment, n = 1L)
    ) %>% 
    replace_na(list(year_idx=1)) %>%
    dplyr::ungroup() %>%
    dplyr::select(CaseID, year_idx) %>%  
    dplyr::right_join(spssData, by = 'CaseID')
}
  
# Fit model.
if (!exists("darknessCount_mod_GLMR_5"))
{
  darknessCount_mod_GLMR_5 <-
    gamm4::gamm4(
      DarknessCrime_sum ~
        1 +
        DiffLfromMnBy100 +
        DiffMSOA_MnLfromGrandMnBy100 +
        s(year_idx, month_idx),
      random = ~
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
darknessCount_pred_GLMR_5 <-
  predict(darknessCount_mod_GLMR_5$gam) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p5_darkness <-
  p_base_darkness +
  geom_line(data = darknessCount_pred_GLMR_5, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle =
         paste0('Black line is the MSOA-specific mean.\n',
                'Random effect: CaseID, MSOA | Fixed effects: # of lamps (centred) \n',
                'AIC: ', round(AIC(darknessCount_mod_GLMR_5$mer)), ' (previous AIC: ', round(AIC(darknessCount_mod_GLMR_3)), ')',
                '\nBIC: ', round(BIC(darknessCount_mod_GLMR_5$mer)), ' (previous BIC: ', round(BIC(darknessCount_mod_GLMR_3)), ')'
         )
  )
# Interpret coefficient and fit.
print_Str <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
print_vals <-
  round(darknessCount_mod_GLMR_5$mer@beta[1:3] %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_Str, print_vals))
print(
  paste0("A log-likelihood ratio test suggests that the inclusion of the ",
         "covariate is statistically significant at any reasonable value:")
)
lmtest::lrtest(darknessCount_mod_GLMR_3,darknessCount_mod_GLMR_5$mer)
print(
  paste0("...plotting the expected values from the model show a very poor fit:")
)
p5_darkness # Show plot.
# ----
darknessCount_mod_GLMR_5$gam %>% summary()
darknessCount_mod_GLMR_5$mer %>% summary()


################################################################################
# Model 6:
# Data structure:
# - Observation-level random effect (ORLE) to handle overdispersion.
# - Random intercept for MSOA.
# - Random slope for MSOA.
# Splines:
# - Year-Month manifold.
#
# The sixth model will include an additional fixed effect for holiday weeks.
################################################################################
# ----
# Make additional variables.
if(!"holiday_wks" %in% colnames(spssData))
{
  spssData$holiday_wks <-
    spssData %>%
    dplyr::select(NewYrWk, GoodFriWk, EasterWk, MayDayWk,
                  SpringBankWk, SummerBankWk, XmasWk) %>%
    rowSums()
}  

# Fit model.
if (!exists("darknessCount_mod_GLMR_6"))
{
  darknessCount_mod_GLMR_6 <-
    gamm4::gamm4(
      DarknessCrime_sum ~
        1 +
        DiffLfromMnBy100 +
        DiffMSOA_MnLfromGrandMnBy100 +
        holiday_wks +
        s(year_idx, month_idx),
      random = ~
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
darknessCount_pred_GLMR_6 <-
  predict(darknessCount_mod_GLMR_6$gam) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p6_darkness <-
  p_base_darkness +
  geom_line(data = darknessCount_pred_GLMR_6, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle =
         paste0('Black line is the MSOA-specific mean.\n',
                'Random effect: CaseID, MSOA | Fixed effects: # of lamps (centred), Holidays \n',
                'AIC: ', round(AIC(darknessCount_mod_GLMR_6$mer)), ' (previous AIC: ', round(AIC(darknessCount_mod_GLMR_5$mer)), ')',
                '\nBIC: ', round(BIC(darknessCount_mod_GLMR_6$mer)), ' (previous BIC: ', round(BIC(darknessCount_mod_GLMR_5$mer)), ')'
         )
  )
# Interpret coefficient and fit.
print_Str <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
print_vals <-
  round(darknessCount_mod_GLMR_6$mer@beta[1:3] %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_Str, print_vals))
print(
  paste0("A log-likelihood ratio test suggests that the inclusion of the ",
         "covariate is statistically significant at any reasonable value:")
)
lmtest::lrtest(darknessCount_mod_GLMR_5$mer,darknessCount_mod_GLMR_6$mer)
print(
  paste0("...plotting the expected values from the model show a very poor fit:")
)
p6_darkness # Show plot.
# ----
darknessCount_mod_GLMR_6$gam %>% summary()
darknessCount_mod_GLMR_6$mer %>% summary()


################################################################################
# Model 7:
# Data structure:
# - Observation-level random effect (ORLE) to handle overdispersion.
# - Random intercept for MSOA.
# - Random slope for MSOA.
# Splines:
# - Year-Month manifold.
#
# The seventh model will include `DifferenceInOffsets` as a spline.
################################################################################
# ----
# Fit model.
if (!exists("darknessCount_mod_GLMR_7")) # Started at 13:41
{
  darknessCount_mod_GLMR_7 <-
    gamm4::gamm4(
      DarknessCrime_sum ~
        1 +
        DiffLfromMnBy100 +
        DiffMSOA_MnLfromGrandMnBy100 +
        holiday_wks +
        s(year_idx, month_idx) +
        s(DifferenceInOffsets),
      random = ~
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
darknessCount_pred_GLMR_7 <-
  predict(darknessCount_mod_GLMR_7$gam) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p7_darkness <-
  p_base_darkness +
  geom_line(data = darknessCount_pred_GLMR_7, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle =
         paste0('Black line is the MSOA-specific mean.\n',
                'Random effect: CaseID, MSOA | Fixed effects: # of lamps (centred), Holidays \n',
                'AIC: ', round(AIC(darknessCount_mod_GLMR_7$mer)), ' (previous AIC: ', round(AIC(darknessCount_mod_GLMR_6$mer)), ')',
                '\nBIC: ', round(BIC(darknessCount_mod_GLMR_7$mer)), ' (previous BIC: ', round(BIC(darknessCount_mod_GLMR_6$mer)), ')'
         )
  )
# Interpret coefficient and fit.
print_Str <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
print_vals <-
  round(darknessCount_mod_GLMR_7$mer@beta[1:3] %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_Str, print_vals))
print(
  paste0("A log-likelihood ratio test suggests that the inclusion of the ",
         "covariate is statistically significant at any reasonable value:")
)
lmtest::lrtest(darknessCount_mod_GLMR_6$mer,darknessCount_mod_GLMR_7$mer)
print(
  paste0("...plotting the expected values from the model show a very poor fit:")
)
p7_darkness # Show plot.
# ----
darknessCount_mod_GLMR_7$gam %>% summary()
darknessCount_mod_GLMR_7$mer %>% summary()


################################################################################
# Model 8:
# Data structure:
# - Observation-level random effect (ORLE) to handle overdispersion.
# - Random intercept for MSOA.
# Splines:
# - Year-Month manifold.
#
# The eighth model will include four lags of `DarknessCrime_sum`, lagging by
# one, two, three, and four weeks.
################################################################################
# ----
# Make additional data.
if(!"lag1_DarknessCrime_sum" %in% colnames(spssData))
{
  spssData <-
    spssData %>%
    dplyr::group_by(MSOAN112) %>%
    dplyr::mutate(
      lag1_DarknessCrime_sum = lag(DarknessCrime_sum),
      lag2_DarknessCrime_sum = lag(DarknessCrime_sum, n = 2L),
      lag3_DarknessCrime_sum = lag(DarknessCrime_sum, n = 3L),
      lag4_DarknessCrime_sum = lag(DarknessCrime_sum, n = 4L),
    )
}

# Fit model.
if (!exists("darknessCount_mod_GLMR_8"))
{
  darknessCount_mod_GLMR_8 <-
    gamm4::gamm4(
      DarknessCrime_sum ~
        1 +
        DiffLfromMnBy100 +
        DiffMSOA_MnLfromGrandMnBy100 +
        holiday_wks +
        lag1_DarknessCrime_sum +
        lag2_DarknessCrime_sum +
        lag3_DarknessCrime_sum +
        lag4_DarknessCrime_sum +
        s(year_idx, month_idx) +
        s(DifferenceInOffsets),
      random = ~
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
darknessCount_pred_GLMR_8 <-
  predict(darknessCount_mod_GLMR_8$gam) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData[!is.na(spssData$lag4_DarknessCrime_sum),] %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p8_darkness <-
  p_base_darkness +
  geom_line(data = darknessCount_pred_GLMR_8, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle =
         paste0('Black line is the MSOA-specific mean.\n',
                'Random effect: CaseID, MSOA | ',
                'Fixed effects: # of lamps (centred), Holidays, Lagged outcome (order 4) \n',
                'AIC: ', round(AIC(darknessCount_mod_GLMR_8$mer)), ' (previous AIC: ', round(AIC(darknessCount_mod_GLMR_7$mer)), ')',
                '\nBIC: ', round(BIC(darknessCount_mod_GLMR_8$mer)), ' (previous BIC: ', round(BIC(darknessCount_mod_GLMR_7$mer)), ')'
         )
  )
# Interpret coefficient and fit.
print_Str <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
print_vals <-
  round(darknessCount_mod_GLMR_8$mer@beta[1:3] %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_Str, print_vals))
print(
  paste0("A log-likelihood ratio test is not possible because the previous model",
         " was fitted to a dataset without the missing lagged values.")
)
print(
  paste0("Plotting the expected values from the model show a very poor fit:")
)
p8_darkness # Show plot.
# ----
darknessCount_mod_GLMR_8$gam %>% summary()
darknessCount_mod_GLMR_8$mer %>% summary()



############################################################################
# Stuff to calculate the overdispersion.
od.point<-function(modelobject){
  x <- sum(resid(modelobject,type="pearson")^2)
  rdf <- summary(modelobject)$AICtab[5]
  return (c(x, rdf, x/rdf))
}

summary(modelobject)$AICtab

modelobject <- daylightCount_mod_GLMR_3
od.point(modelobject)
x <- sum(resid(daylightCount_mod_GLMR_8$gam,type="pearson")^2); rdf <- summary(daylightCount_mod_GLMR_8$mer)$AICtab[5]; (x/rdf)

# From $gam
darknessCount_mod_GLMR_5$gam$df.residual

# From $mer
############################################################################






################################################################################
#####                                                                      #####
#####                    VARIATE: DaylightCrime_sum                        #####
#####                                                                      #####
################################################################################

################################################################################
# Plot base:
# In the following plot, I randomly sample ten MSOAs and show their sum of crimes
# in darkness hours, over the data collection period. This plot will form the
# base of all future plots, which will overlay the expectations from increasingly
# complicated statistical models of the sum of crimes in darkness hours.
#
# I also fit the base Poisson model for the count of crimes in darkness hours,
# fit using a basic GLM, stas::glm().
################################################################################
# ----
# Fit null model.
if (!exists("daylightCount_mod_GLMR_0"))
{
  daylightCount_mod_GLMR_0 <-
    glm(DaylightCrime_sum ~ 1,
        family = "poisson",
        data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
daylightCount_pred_GLMR_0 <-
  predict(daylightCount_mod_GLMR_0) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Plot.
(p_base_daylight <-
    spssData %>%
    dplyr::select(DaylightCrime_sum, MSOAN112, WkNoStartFrom1) %>%
    dplyr::filter(!MSOAN112 %in% c(5, 111)) %>%
    dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs) %>%
    ggplot() +
    geom_line(aes(x = WkNoStartFrom1, y = DaylightCrime_sum), color = "darkgrey") +
    facet_wrap(vars(MSOAN112), ncol = 5) +
    labs(title = 'Number of crimes in the hours of daylight for a random sample of 10 MSOAs.',
         x = 'Week number', y = 'Weekly count of crimes in darkness hours')
)
# ----
daylightCount_mod_GLMR_0 %>% summary()


################################################################################
# Model 1:
# Data structure:
# - Grand mean
#
# The most-basic model is the arithmetic average sum of crimes in darkness hours, 
# In the following figure, I will overlay this value over each plot. I will also
# note the fit statistic, AIC.
################################################################################
# ----
# Fit model.
if (!exists("daylightCount_mod_GLMR_1"))
{
  daylightCount_mod_GLMR_1 <-
    glm(DaylightCrime_sum ~ 
          1 +
          DiffLfromMnBy100,
        family = "poisson",
        data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
daylightCount_pred_GLMR_1 <-
  predict(daylightCount_mod_GLMR_1) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)

# Make plot.
p1_daylight <-
  p_base_daylight +
  geom_line(data = daylightCount_pred_GLMR_1, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle = paste0('Black line is the model including number of lights.\n',
                         'AIC: ', round(daylightCount_mod_GLMR_1$aic),
                         '\nBIC: ', round(BIC(daylightCount_mod_GLMR_1)) )
  )

# Interpret coefficient and fit.
print_Str <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
print_vals <-
  round(daylightCount_mod_GLMR_1$coefficients %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_Str, print_vals))

print(
  paste0("A log-likelihood ratio test suggests that the inclusion of the ",
         "covariate is statistically significant at any reasonable value:")
)
lmtest::lrtest(daylightCount_mod_GLMR_1)
print(
  paste0("...plotting the expected values from the model show a very poor fit:")
)
p1_daylight # Show plot.
print(
  paste0("The subsequent models represent progressive attempts to improve fit.",
         " I will note the coefficient for `DiffLfromMnBy100` as the models progress.")
)
# ----
daylightCount_mod_GLMR_1 %>% summary()


################################################################################
# Model 2:
# Data structure:
# - Observation-level random effect (ORLE) to handle overdispersion.
#
# The second model will address the concern of overdispersion that is raised
# because the variance of the count variate is many times the mean.
################################################################################
# ----
# Fit model.
if (!exists("daylightCount_mod_GLMR_2"))
{
  daylightCount_mod_GLMR_2 <-
    lme4::glmer(
      DaylightCrime_sum ~
        1 +
        DiffLfromMnBy100 +
        (1 | CaseID),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
daylightCount_pred_GLMR_2 <-
  predict(daylightCount_mod_GLMR_2) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p2_daylight <-
  p_base_daylight +
  geom_line(data = daylightCount_pred_GLMR_2, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle =
         paste0('Black line is the MSOA-specific mean.\n',
                'Random effect: CaseID | Fixed effect: # of lights (centred) \n',
                'AIC: ', round(AIC(daylightCount_mod_GLMR_2)), ' (previous AIC: ', round(AIC(daylightCount_mod_GLMR_1)), ')',
                '\nBIC: ', round(BIC(daylightCount_mod_GLMR_2)), ' (previous BIC: ', round(BIC(daylightCount_mod_GLMR_1)), ')'
         )
  )

# Interpret coefficient and fit.
print_Str <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
print_vals <-
  round(daylightCount_mod_GLMR_2@beta %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_Str, print_vals))
print(
  paste0("A log-likelihood ratio test suggests that the inclusion of the ",
         "covariate is statistically significant at any reasonable value:")
)
lmtest::lrtest(daylightCount_mod_GLMR_1,daylightCount_mod_GLMR_2)
print(
  paste0("...plotting the expected values from the model show a very poor fit:")
)
# Show plot.
p2_daylight

# Fit model for comparing deviance.
if (!exists("daylightCount_mod_GLMR_2_noCovars"))
{
  daylightCount_mod_GLMR_2_noCovars <-
    lme4::glmer(
      DaylightCrime_sum ~
        1 +
        (1 | CaseID),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
# ----
daylightCount_mod_GLMR_2 %>% summary()


################################################################################
# Model 3:
# Data structure:
# - Observation-level random effect (ORLE) to handle overdispersion.
# - Random intercept for MSOA.
#
# The third model will respect the data structure by including a random
# intercept for `MSOAN112`. I also include a fixed-effect variable to represent
# each MSOA's count of lamps, centred to the grand mean.
################################################################################
# ----
# Fit model.
if (!exists("daylightCount_mod_GLMR_3"))
{
  daylightCount_mod_GLMR_3 <-
    lme4::glmer(
      DaylightCrime_sum ~
        1 +
        DiffLfromMnBy100 +
        DiffMSOA_MnLfromGrandMnBy100 +
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
daylightCount_pred_GLMR_3 <-
  predict(daylightCount_mod_GLMR_3) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p3_daylight <-
  p_base_daylight +
  geom_line(data = daylightCount_pred_GLMR_3, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle =
         paste0('Black line is the MSOA-specific mean.\n',
                'Random effect: CaseID, MSOA | Fixed effect: # of lights (centred) \n',
                'AIC: ', round(AIC(daylightCount_mod_GLMR_3)), ' (previous AIC: ', round(AIC(daylightCount_mod_GLMR_2)), ')',
                '\nBIC: ', round(BIC(daylightCount_mod_GLMR_3)), ' (previous BIC: ', round(BIC(daylightCount_mod_GLMR_2)), ')'
         )
  )
# Interpret coefficient and fit.
print_Str <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
print_vals <-
  round(daylightCount_mod_GLMR_3@beta %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_Str, print_vals))
print(
  paste0("A log-likelihood ratio test suggests that the inclusion of the ",
         "covariate is statistically significant at any reasonable value:")
)
lmtest::lrtest(daylightCount_mod_GLMR_2,daylightCount_mod_GLMR_3)
print(
  paste0("...plotting the expected values from the model show a very poor fit:")
)
# Show plot.
p3_daylight

# Fit model for comparing deviance.
if (!exists("daylightCount_mod_GLMR_3_noCovars"))
{
  daylightCount_mod_GLMR_3_noCovars <-
    lme4::glmer(
      DaylightCrime_sum ~
        1 +
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
# ----
daylightCount_mod_GLMR_3 %>% summary()


# ********************
# **** NOT FITTED ****
# ********************
################################################################################
# Model 4:
# Data structure:
# - Observation-level random effect (ORLE) to handle overdispersion.
# - Random intercept for MSOA.
# - Random slope for MSOA.
#
# The fourth model will check to see if a random slope for MSOAN112 also
# improves fit.
################################################################################
# ----
# # Fit model.
# if (!exists("daylightCount_mod_GLMR_4"))
# {
#   daylightCount_mod_GLMR_4 <-
#     lme4::glmer(
#       DaylightCrime_sum ~
#         1 +
#         DiffLfromMnBy100 +
#         (1 | CaseID) +
#         (1 + MSOAN112 | MSOAN112),
#       family = "poisson",
#       data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
#     )
# }
# daylightCount_pred_GLMR_4 <-
#   predict(daylightCount_mod_GLMR_4) %>%
#   as.data.frame() %>%
#   dplyr::bind_cols(spssData %>%
#                      dplyr::select(WkNoStartFrom1, MSOAN112) %>%
#                      dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
#   dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# # Make plot.
# p4_daylight <-
#   p_base_daylight +
#   geom_line(data = daylightCount_pred_GLMR_4, aes(x = WkNoStartFrom1, y = .)) +
#   labs(subtitle =
#          paste0('Black line is the MSOA-specific mean.\n',
#                 'Random effect: CaseID, MSOA | Fixed effect: # of lights (centred) \n',
#                 'AIC: ', round(AIC(daylightCount_mod_GLMR_4)), ' (previous AIC: ', round(AIC(daylightCount_mod_GLMR_3)), ')',
#                 '\nBIC: ', round(BIC(daylightCount_mod_GLMR_4)), ' (previous BIC: ', round(BIC(daylightCount_mod_GLMR_3)), ')'
#          )
#   )
# # Interpret coefficient and fit.
# print_Str <- 
#   paste0("The exponentiated coefficient is %g, which we interpret as saying ",
#          "that the count of crimes committed in the hours of darkness are ",
#          "expected to be %g-times greater with the addition of one light.")
# print_vals <-
#   round(daylightCount_mod_GLMR_4@beta %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
# do.call(sprintf, c(fmt = print_Str, print_vals))
# print(
#   paste0("A log-likelihood ratio test suggests that the inclusion of the ",
#          "covariate is statistically significant at any reasonable value:")
# )
# lmtest::lrtest(daylightCount_mod_GLMR_3,daylightCount_mod_GLMR_4)
# print(
#   paste0("...plotting the expected values from the model show a very poor fit:")
# )
# p4_daylight # Show plot.
# ----
#daylightCount_mod_GLMR_4 %>% summary()


################################################################################
# Model 5:
# Data structure:
# - Observation-level random effect (ORLE) to handle overdispersion.
# - Random intercept for MSOA.
# Splines:
# - Year-Month manifold.
#
# The fifth model will use gamm4::gamm4() to include a joint spline for year
# and month to account for a creeping annual trend and seasonality. I'll need
# to create new variables for this. Both will have to be monotonically-increasing
# integer variables, over the entire study period. Note that when predicting
# using gamm4::gamm4(), you have to predict from the .$gam object in the list.
################################################################################
# ----
# Make additional variables.
if(!"month_idx" %in% colnames(spssData))
{
  spssData$month_idx <-
    spssData %>%
    dplyr::select(starts_with("MonthMidWk_")) %>%
    tidyr::unite(col = unitedMonths, sep = "") %>%
    dplyr::mutate(
      month_idx = dplyr::case_when(
        . == '00000000000' ~ '1',
        . == '10000000000' ~ '2',
        . == '01000000000' ~ '3',
        . == '00100000000' ~ '4',
        . == '00010000000' ~ '5',
        . == '00001000000' ~ '6',
        . == '00000100000' ~ '7',
        . == '00000010000' ~ '8',
        . == '00000001000' ~ '9',
        . == '00000000100' ~ '10',
        . == '00000000010' ~ '11',
        . == '00000000001' ~ '12',
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(month_idx) %>%
    unlist() %>%
    as.numeric()
}

if(!"year_idx" %in% colnames(spssData))
{
  spssData <-
    spssData %>%
    group_by(MSOAN112) %>% 
    dplyr::select(CaseID, month_idx) %>%
    mutate(increment = 1 + cumsum( (month_idx == 12) & (lead(month_idx) == 1) )) %>% 
    dplyr::mutate(
      year_idx = dplyr::lag(increment, n = 1L)
    ) %>% 
    replace_na(list(year_idx=1)) %>%
    dplyr::ungroup() %>%
    dplyr::select(CaseID, year_idx) %>%  
    dplyr::right_join(spssData, by = 'CaseID')
}

# Fit model.
if (!exists("daylightCount_mod_GLMR_5"))
{
  daylightCount_mod_GLMR_5 <-
    gamm4::gamm4(
      DaylightCrime_sum ~
        1 +
        DiffLfromMnBy100 +
        DiffMSOA_MnLfromGrandMnBy100 +
        s(year_idx, month_idx),
      random = ~
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
daylightCount_pred_GLMR_5 <-
  predict(daylightCount_mod_GLMR_5$gam) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p5_daylight <-
  p_base_daylight +
  geom_line(data = daylightCount_pred_GLMR_5, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle =
         paste0('Black line is the MSOA-specific mean.\n',
                'Random effect: CaseID, MSOA | Fixed effects: # of lamps (centred) \n',
                'AIC: ', round(AIC(daylightCount_mod_GLMR_5$mer)), ' (previous AIC: ', round(AIC(daylightCount_mod_GLMR_3)), ')',
                '\nBIC: ', round(BIC(daylightCount_mod_GLMR_5$mer)), ' (previous BIC: ', round(BIC(daylightCount_mod_GLMR_3)), ')'
         )
  )
# Interpret coefficient and fit.
print_Str <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
print_vals <-
  round(daylightCount_mod_GLMR_5$mer@beta[1:3] %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_Str, print_vals))
print(
  paste0("A log-likelihood ratio test suggests that the inclusion of the ",
         "covariate is statistically significant at any reasonable value:")
)
lmtest::lrtest(daylightCount_mod_GLMR_3,daylightCount_mod_GLMR_5$mer)
print(
  paste0("...plotting the expected values from the model show a very poor fit:")
)
# Show plot.
p5_daylight 

# Fit model for comparing deviance.
if (!exists("daylightCount_mod_GLMR_5_noCovars"))
{
  daylightCount_mod_GLMR_5_noCovars <-
    gamm4::gamm4(
      DaylightCrime_sum ~
        1 +
        s(year_idx, month_idx),
      random = ~
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
# ----
daylightCount_mod_GLMR_5$gam %>% summary()
daylightCount_mod_GLMR_5$mer %>% summary()


################################################################################
# Model 6:
# Data structure:
# - Observation-level random effect (ORLE) to handle overdispersion.
# - Random intercept for MSOA.
# - Random slope for MSOA.
# Splines:
# - Year-Month manifold.
#
# The sixth model will include an additional fixed effect for holiday weeks.
################################################################################
# ----
# Make additional variables.
if(!"holiday_wks" %in% colnames(spssData))
{
  spssData$holiday_wks <-
    spssData %>%
    dplyr::select(NewYrWk, GoodFriWk, EasterWk, MayDayWk,
                  SpringBankWk, SummerBankWk, XmasWk) %>%
    rowSums()
}

# Fit model.
if (!exists("daylightCount_mod_GLMR_6"))
{
  daylightCount_mod_GLMR_6 <-
    gamm4::gamm4(
      DaylightCrime_sum ~
        1 +
        DiffLfromMnBy100 +
        DiffMSOA_MnLfromGrandMnBy100 +
        holiday_wks +
        s(year_idx, month_idx),
      random = ~
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
daylightCount_pred_GLMR_6 <-
  predict(daylightCount_mod_GLMR_6$gam) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p6_daylight <-
  p_base_daylight +
  geom_line(data = daylightCount_pred_GLMR_6, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle =
         paste0('Black line is the MSOA-specific mean.\n',
                'Random effect: CaseID, MSOA | Fixed effects: # of lamps (centred), Holidays \n',
                'AIC: ', round(AIC(daylightCount_mod_GLMR_6$mer)), ' (previous AIC: ', round(AIC(daylightCount_mod_GLMR_5$mer)), ')',
                '\nBIC: ', round(BIC(daylightCount_mod_GLMR_6$mer)), ' (previous BIC: ', round(BIC(daylightCount_mod_GLMR_5$mer)), ')'
         )
  )
# Interpret coefficient and fit.
print_Str <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
print_vals <-
  round(daylightCount_mod_GLMR_6$mer@beta[1:3] %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_Str, print_vals))
print(
  paste0("A log-likelihood ratio test suggests that the inclusion of the ",
         "covariate is statistically significant at any reasonable value:")
)
lmtest::lrtest(daylightCount_mod_GLMR_5$mer,daylightCount_mod_GLMR_6$mer)
print(
  paste0("...plotting the expected values from the model show a very poor fit:")
)
p6_daylight # Show plot.
# ----
daylightCount_mod_GLMR_6$gam %>% summary()
daylightCount_mod_GLMR_6$mer %>% summary()


################################################################################
# Model 7:
# Data structure:
# - Observation-level random effect (ORLE) to handle overdispersion.
# - Random intercept for MSOA.
# - Random slope for MSOA.
# Splines:
# - Year-Month manifold.
#
# The seventh model will include `DifferenceInOffsets` as a spline.
################################################################################
# ----
# Fit model.
if (!exists("daylightCount_mod_GLMR_7")) # Started at 13:41
{
  daylightCount_mod_GLMR_7 <-
    gamm4::gamm4(
      DaylightCrime_sum ~
        1 +
        DiffLfromMnBy100 +
        DiffMSOA_MnLfromGrandMnBy100 +
        holiday_wks +
        s(year_idx, month_idx) +
        s(DifferenceInOffsets),
      random = ~
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
daylightCount_pred_GLMR_7 <-
  predict(daylightCount_mod_GLMR_7$gam) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p7_daylight <-
  p_base_daylight +
  geom_line(data = daylightCount_pred_GLMR_7, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle =
         paste0('Black line is the MSOA-specific mean.\n',
                'Random effect: CaseID, MSOA | Fixed effects: # of lamps (centred), Holidays \n',
                'AIC: ', round(AIC(daylightCount_mod_GLMR_7$mer)), ' (previous AIC: ', round(AIC(daylightCount_mod_GLMR_6$mer)), ')',
                '\nBIC: ', round(BIC(daylightCount_mod_GLMR_7$mer)), ' (previous BIC: ', round(BIC(daylightCount_mod_GLMR_6$mer)), ')'
         )
  )
# Interpret coefficient and fit.
print_Str <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
print_vals <-
  round(daylightCount_mod_GLMR_7$mer@beta[1:3] %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_Str, print_vals))
print(
  paste0("A log-likelihood ratio test suggests that the inclusion of the ",
         "covariate is statistically significant at any reasonable value:")
)
lmtest::lrtest(daylightCount_mod_GLMR_6$mer,daylightCount_mod_GLMR_7$mer)
print(
  paste0("...plotting the expected values from the model show a very poor fit:")
)
# Show plot.
p7_daylight 

# Fit model for comparing deviance.
if (!exists("daylightCount_mod_GLMR_7_noCovars"))
{
  daylightCount_mod_GLMR_7_noCovars <-
    gamm4::gamm4(
      DaylightCrime_sum ~
        1 +
        s(year_idx, month_idx) +
        s(DifferenceInOffsets),
      random = ~
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
# ----
daylightCount_mod_GLMR_7$gam %>% summary()
daylightCount_mod_GLMR_7$mer %>% summary()


################################################################################
# Model 8:
# Data structure:
# - Observation-level random effect (ORLE) to handle overdispersion.
# - Random intercept for MSOA.
# Splines:
# - Year-Month manifold.
#
# The eighth model will include four lags of `DaylightCrime_sum`, lagging by
# one, two, three, and four weeks.
################################################################################
# ----
# Make additional data.
if(!"lag1_DaylightCrime_sum" %in% colnames(spssData))
{
  spssData <-
    spssData %>%
    dplyr::group_by(MSOAN112) %>%
    dplyr::mutate(
      lag1_DaylightCrime_sum = lag(DaylightCrime_sum),
      lag2_DaylightCrime_sum = lag(DaylightCrime_sum, n = 2L),
      lag3_DaylightCrime_sum = lag(DaylightCrime_sum, n = 3L),
      lag4_DaylightCrime_sum = lag(DaylightCrime_sum, n = 4L),
    )
}

# Fit model.
if (!exists("daylightCount_mod_GLMR_8"))
{
  daylightCount_mod_GLMR_8 <-
    gamm4::gamm4(
      DaylightCrime_sum ~
        1 +
        DiffLfromMnBy100 +
        DiffMSOA_MnLfromGrandMnBy100 +
        holiday_wks +
        lag1_DaylightCrime_sum +
        lag2_DaylightCrime_sum +
        lag3_DaylightCrime_sum +
        lag4_DaylightCrime_sum +
        s(year_idx, month_idx) +
        s(DifferenceInOffsets),
      random = ~
        (1 | CaseID) +
        (1 | MSOAN112),
      family = "poisson",
      data = spssData %>% dplyr::filter(!MSOAN112 %in% c(5, 111))
    )
}
daylightCount_pred_GLMR_8 <-
  predict(daylightCount_mod_GLMR_8$gam) %>%
  as.data.frame() %>%
  dplyr::bind_cols(spssData[!is.na(spssData$lag4_DaylightCrime_sum),] %>%
                     dplyr::select(WkNoStartFrom1, MSOAN112) %>%
                     dplyr::filter(!MSOAN112 %in% c(5, 111))) %>%
  dplyr::filter(MSOAN112 %in% random_selection_of_MSOAs)
# Make plot.
p8_daylight <-
  p_base_daylight +
  geom_line(data = daylightCount_pred_GLMR_8, aes(x = WkNoStartFrom1, y = .)) +
  labs(subtitle =
         paste0('Black line is the MSOA-specific mean.\n',
                'Random effect: CaseID, MSOA | ',
                'Fixed effects: # of lamps (centred), Holidays, Lagged outcome (order 4) \n',
                'AIC: ', round(AIC(daylightCount_mod_GLMR_8$mer)), ' (previous AIC: ', round(AIC(daylightCount_mod_GLMR_7$mer)), ')',
                '\nBIC: ', round(BIC(daylightCount_mod_GLMR_8$mer)), ' (previous BIC: ', round(BIC(daylightCount_mod_GLMR_7$mer)), ')'
         )
  )
# Interpret coefficient and fit.
print_Str <- 
  paste0("The exponentiated coefficient is %g, which we interpret as saying ",
         "that the count of crimes committed in the hours of darkness are ",
         "expected to be %g-times greater with the addition of one light.")
print_vals <-
  round(daylightCount_mod_GLMR_8$mer@beta[1:3] %>% sum() %>% exp(), 2) %>% rep(,2) %>% as.list()
do.call(sprintf, c(fmt = print_Str, print_vals))
print(
  paste0("A log-likelihood ratio test is not possible because the previous model",
         " was fitted to a dataset without the missing lagged values.")
)
print(
  paste0("Plotting the expected values from the model show a very poor fit:")
)
p8_daylight # Show plot.
# ----
daylightCount_mod_GLMR_8$gam %>% summary()
daylightCount_mod_GLMR_8$mer %>% summary()

