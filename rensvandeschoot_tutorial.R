# This script mimics Rens van de Schoot's tutorial for the `lme4` package
# at https://www.rensvandeschoot.com/tutorials/lme4/
#

###############
# Requisites. #
###############
# ----
if (!("pacman" %in% installed.packages()[,"Package"])) install.packages("pacman")
library(pacman)
pacman::p_load(haven, tidyverse, lme4, RColorBrewer)
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
  mutate(across(everything(), as.vector))

popular2data <- read_sav(file ="https://github.com/MultiLevelAnalysis/Datasets-third-edition-Multilevel-book/blob/master/chapter%202/popularity/SPSS/popular2.sav?raw=true")
popular2data <- select(popular2data, pupil, class, extrav, sex, texp, popular)
# ----

####################
# Initial plotting #
####################
# ----
ggplot(data  = popular2data,
       aes(x = extrav,
           y = popular))+
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter")+# to add some random noise for plotting purposes
  theme_minimal()+
  labs(title = "Popularity vs. Extraversion")

ggplot(data  = spssData,
       aes(x = DiffLfromMnBy100,
           y = PropDarkCrime))+
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter")+
  theme_minimal()+
  labs(title = "Proportion of weekly crimes that occured in the dark hours vs. Number of lamps minus its MSOA mean number of lamps")
# ----

#######################
# Plot with lm() line #
#######################
# ----
ggplot(data  = popular2data,
       aes(x = extrav,
           y = popular))+
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  geom_smooth(method = lm,
              se     = FALSE, 
              col    = "black",
              size   = .5, 
              alpha  = .8)+ # to add regression line
  theme_minimal()+
  labs(title    = "Popularity vs. Extraversion",
       subtitle = "add regression line")

ggplot(data  = spssData,
       aes(x = DiffLfromMnBy100,
           y = PropDarkCrime))+
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter")+
  geom_smooth(method = lm,
              se     = FALSE, 
              col    = "red")+
  theme_minimal()+
  labs(title = "Proportion of weekly crimes that occured in the dark hours vs. Number of lamps minus its MSOA mean number of lamps",
       subtitle = "add LOESS line")
# ----

######################
# Plot separate MSOA # 
######################
# ----
ggplot(data    = popular2data,
       aes(x   = extrav,
           y   = popular,
           col = class))+ #to add the colours for different classes
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_minimal()+
  theme(legend.position = "none")+
  scale_color_gradientn(colours = rainbow(100))+
  labs(title    = "Popularity vs. Extraversion",
       subtitle = "add colours for different classes")

ggplot(data  = spssData,
       aes(x = DiffLfromMnBy100,
           y = PropDarkCrime,
           col = MSOAN112))+
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter")+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_color_gradientn(colours = rainbow(100))+
  labs(title = "Proportion of weekly crimes that occured in the dark hours vs. Number of lamps minus its MSOA mean number of lamps",
       subtitle = "add colurs for different MSOAs")
# ----

#################################
# Plot with LOESS for each MSOA #
#################################
# ----
ggplot(data      = popular2data,
       aes(x     = extrav,
           y     = popular,
           col   = class,
           group = class))+ #to add the colours for different classes
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+ #to add some random noise for plotting purposes
  theme_minimal()+
  theme(legend.position = "none")+
  scale_color_gradientn(colours = rainbow(100))+
  geom_smooth(method = lm,
              se     = FALSE,
              size   = .5, 
              alpha  = .8)+ # to add regression line
  labs(title    = "Popularity vs. Extraversion",
       subtitle = "add colours for different classes and regression lines")

ggplot(data      = spssData,
       aes(x     = DiffLfromMnBy100,
           y     = PropDarkCrime,
           group = MSOAN112))+
  geom_point(size     = 1.2,
             alpha    = .8,
             position = "jitter")+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_color_gradientn(colours = rainbow( length( unique(spssData$MSOAN112) ) ) )+
  geom_smooth(method = lm,
              se     = FALSE,
              size   = .5, 
              alpha  = .8)+
  labs(title = "Proportion of weekly crimes that occured in the dark hours vs. Number of lamps minus its MSOA mean number of lamps",
       subtitle = "add colours for different MSOAs and regression lines")
# ----

###################################
# Fit simplest mixed-effect model #
###################################
# ----
(his_interceptonlymodel <-
  lmer(formula = popular ~ 1 + (1 | class),
       data    = popular2data))

# The `lme4::glmer()` function does not accept quasi-family distriubtions. There
# are a few methods that we can use as alternatives:
# https://www.researchgate.net/post/Why-cant-i-use-the-family-quasi-in-generalised-linear-mixed-model-GLMM
(my_interceptonlymodel <-
  glmer(formula = PropDarkCrime ~ 1 + (1 | MSOAN112),
        family = binomial(link = "logit"),
        data = spssData))
# Fixed intercept is 0.474, i.e. slightly fewer crimes at night.
# For random effects, the residual variance for a crime is 0.080 and the residual
# variance at the MSOA level is 0.002. Not much in either case.
# ----

########################################################################
# Plot different lm() lines if we stratified on fixed-effect variables #
########################################################################
# ----
ggplot(data = popular2data, 
       aes(x   = extrav,
           y   = popular, 
           col = as.factor(sex)))+
  geom_point(size     = 1, 
             alpha    = .7, 
             position = "jitter")+
  geom_smooth(method   = lm,
              se       = T, 
              size     = 1.5, 
              linetype = 1, 
              alpha    = .7)+
  theme_minimal()+
  labs(title    = "Linear Relationship Between Popularity and Extraversion for the 2 Genders", 
       subtitle = "The linear relationship between the two is similar for both genders, with a clear intercept difference")+
  scale_color_manual(name   =" Gender",
                     labels = c("Boys", "Girls"),
                     values = c("lightblue", "pink"))

# NO EQUIVALENT FOR OUR ALL-NUMERIC DATASET.

# ----

#############################################
# Fit mixed-effect model with fixed effects #
#############################################
# ----
(his_model1 <- lmer(formula = popular ~ 1 + sex + extrav + (1|class), 
                data = popular2data))

(my_interceptonlymodel <-
    glmer(formula = PropDarkCrime ~ 1 +
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
           (1 | MSOAN112),
         data    = spssData))
# Fixed intercept is 0.55, i.e. slightly more crimes at night.
# The fixed-effect of interest is that belonging to `DiffLfromMnBy100`, which
# is 0.001, i.e. having an additional 100 lamps more than the MSOA's average was
# associated with 0.1% greater proportion of crimes at night.
# ----
