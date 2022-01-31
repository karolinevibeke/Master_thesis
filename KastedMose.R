################################################################################
#   Karoline Vibeke Dohrmann  //  au540949  //  201507476                      #
#                                                                              #
#   kvdohr@gmail.com  //  au540949@uni.au.dk                                   #
#                                                                              #
#   Master's Thesis in Biology at Aarhus University                            #
#                                                                              #
#   "Spatiotemporal dynamics of large herbivore species Bos taurus and         #
#    Bubalus bubalis in trophic rewilding project Kasted Mose"                 #
#                                                                              #
#   Data available from github repository:                                     #
#   https://github.com/karolinevibeke/Master_thesis.git                        #
#                                                                              #
#   Last revised 30/01/2022                                                    #
################################################################################

# INSTALL AND LOAD PACKAGES ----

# Load package plyr before dplyr to avoid conflicts.
# install.packages("plyr")
library(plyr)

packages <- c("tidyverse",
              "gridExtra", 
              "colourpicker", 
              "egg", 
              "lme4",
              "lmerTest",
              "beanplot", 
              "reshape2", 
              "corrplot", 
              "RColorBrewer", 
              "scales", 
              "ggpubr", 
              "gtable",
              "grid", 
              "ggmap", 
              "RCurl", 
              "rphylopic", 
              "png",
              "car",
              "chisq.posthoc.test",
              "stargazer",
              "ggResidpanel",
              "MuMIn",
              "readr")


# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load packages
lapply(packages, library, character.only = TRUE)


# IMPORT AND INSPECT DATA ----

# Load data files from GitHub:
BOS_TAURUS_SEP <- "https://raw.githubusercontent.com/karolinevibeke/Master_thesis/main/bos_group_sep.csv"
BOS_TAURUS_SEP <- read.csv(url(BOS_TAURUS_SEP))
        
BOS_TAURUS_OCT <- "https://raw.githubusercontent.com/karolinevibeke/Master_thesis/main/bos_group_oct.csv"
BOS_TAURUS_OCT <- read.csv(url(BOS_TAURUS_OCT))

BOS_TAURUS_NOV <- "https://raw.githubusercontent.com/karolinevibeke/Master_thesis/main/bos_group_nov.csv"
BOS_TAURUS_NOV <- read.csv(url(BOS_TAURUS_NOV))

BOS_TAURUS_DEC <- "https://raw.githubusercontent.com/karolinevibeke/Master_thesis/main/bos_group_dec.csv"
BOS_TAURUS_DEC <- read.csv(url(BOS_TAURUS_DEC))

BUBALUS_BUBALIS_SEP <- "https://raw.githubusercontent.com/karolinevibeke/Master_thesis/main/bub_group_sep.csv"
BUBALUS_BUBALIS_SEP <- read.csv(url(BUBALUS_BUBALIS_SEP))

BUBALUS_BUBALIS_OCT <- "https://raw.githubusercontent.com/karolinevibeke/Master_thesis/main/bub_group_oct.csv"
BUBALUS_BUBALIS_OCT <- read.csv(url(BUBALUS_BUBALIS_OCT))

BUBALUS_BUBALIS_NOV <- "https://raw.githubusercontent.com/karolinevibeke/Master_thesis/main/bub_group_nov.csv"
BUBALUS_BUBALIS_NOV <- read.csv(url(BUBALUS_BUBALIS_NOV))

BUBALUS_BUBALIS_DEC <- "https://raw.githubusercontent.com/karolinevibeke/Master_thesis/main/bub_group_dec.csv"
BUBALUS_BUBALIS_DEC <- read.csv(url(BUBALUS_BUBALIS_DEC))

# Data on dominant group behaviour for each cell:
DATA_BEHAVIOUR <- "https://raw.githubusercontent.com/karolinevibeke/Master_thesis/main/dominant_group_behaviour.csv"
DATA_BEHAVIOUR <- read.csv(url(DATA_BEHAVIOUR))

# Data on NDVI, TWI and distance to paths for each cell:
DATA_CELLS <- "https://raw.githubusercontent.com/karolinevibeke/Master_thesis/main/cell_variables.csv"
DATA_CELLS <- read.csv(url(DATA_CELLS))

# Data on temperature and precipitation:
DATA_CLIMATE <- "https://raw.githubusercontent.com/karolinevibeke/Master_thesis/main/weather.csv"
DATA_CLIMATE <- read.csv(url(DATA_CLIMATE))

# Data on hour of each group observation and presence of humans:
DATA_HOUR <- "https://raw.githubusercontent.com/karolinevibeke/Master_thesis/main/hour.csv"
DATA_HOUR <- read.csv(url(DATA_HOUR))

# It is checked that the data has been imported without any mistakes/flaws, i.e.
# looking at column names, number of rows, the types of variables (continuous,
# numeric, categorical, characters) etc.


# WRANGLING - TIDY AND TRANSFORM DATA ----

## CLIMATE DATASET ----

str(DATA_CLIMATE)

# Renaming variables:
DATA_CLIMATE <- DATA_CLIMATE %>%
                rename(MONTH = month,
                       DATE = date,
                       PREC_MM = precip_mm)

# Gathering temperature by key (min, mean, max):
DATA_CLIMATE <- gather(DATA_CLIMATE, TEMP_KEY, TEMP_C,
                       c(min_temp, mean_temp, max_temp))

# Checking variable classes and changing some: 
str(DATA_CLIMATE)

DATA_CLIMATE <- DATA_CLIMATE %>% 
                mutate(across(c(MONTH, TEMP_KEY), 
                                as.factor))

DATA_CLIMATE$DATE <- as.Date(DATA_CLIMATE$DATE, 
                             "%m/%d/%y")

# Checking that variable classes have been changed correctly:
str(DATA_CLIMATE)

# Reordering the variable MONTH: 
DATA_CLIMATE$MONTH <- factor(DATA_CLIMATE$MONTH, 
                             levels = c("sep", "oct", "nov", "dec"),
                             labels = c("sep", "oct", "nov", "dec"))


## FULL DATASET ----

# List of dataframes with observations:
LIST_OBS <- list(BOS_TAURUS_SEP, BOS_TAURUS_OCT, BOS_TAURUS_NOV, BOS_TAURUS_DEC,
                 BUBALUS_BUBALIS_SEP, BUBALUS_BUBALIS_OCT, BUBALUS_BUBALIS_NOV,
                 BUBALUS_BUBALIS_DEC)

# New column in each data frame with percentage observation in all cells: 
PROP_fun <- function(x) { 
            x %>%       
            mutate(x, PROP = COUNT / sum(COUNT))
                   }

LIST_OBS <- LIST_OBS %>%
            lapply(PROP_fun)

# The data frames containing the observational data are merged, to have 
# all the information in one table:
DATA_OBS <- as.data.frame(do.call(rbind, LIST_OBS))

# Some of the column name are changed to ensure consistency:
DATA_CELLS <- rename(DATA_CELLS, sep = NDVI_SEP, oct = NDVI_OCT, 
                     nov = NDVI_NOV, dec = NDVI_DEC, TWI = CELL_TWI)

# The data frame format is changed from wide to long with the tidyr::gather() function.
# Gathering NDVI values by MONTH:
DATA_CELLS <- gather(DATA_CELLS, MONTH, NDVI,
                     c(sep, oct, nov, dec))

# The data frames are joined based on shared attributes.
# By using left_join, all information is kept:
DATA_ALL <- left_join(DATA_OBS, DATA_CELLS, by = c("CELL_ID" = "CELL_ID", 
                                                   "MONTH" = "MONTH"))

# A new column is created for presence/absence data. 
# Return 1 for the cells with observations, and 0 for the ones without:
DATA_ALL <- mutate(DATA_ALL,
              PRESENCE = ifelse(COUNT %in% 0,0,1))

# Inspect the types of variables in the data frame:
str(DATA_ALL)

# Some of the variables has to be transformed into factors (i.e. categorical variables).
# The variables are converted and the original class is overwritten.
# The function dplyr::mutate() is combined with the function across() to modify 
# the columns and the function as.factor() is applied to the selected columns:
DATA_ALL <- DATA_ALL %>%
              mutate(across(c(MONTH, SPECIES, PRESENCE), 
                              as.factor))

# Check that the classes have been changed:
str(DATA_ALL)

# Reordering values in the MONTH variable to appear in a chronological order 
# instead of alphabetical: 
DATA_ALL$MONTH <- factor(DATA_ALL$MONTH, 
                            levels = c("sep", "oct", "nov", "dec"),
                            labels = c("sep", "oct", "nov", "dec"))

## DATA SUBSET - PRESENCE ----

# Create subset only containing data on cells with observations:
DATA_PRESENCE <- filter(DATA_ALL, PRESENCE == 1)

# Reordering values in the MONTH variable: 
DATA_PRESENCE$MONTH <- factor(DATA_PRESENCE$MONTH, 
                         levels = c("sep", "oct", "nov", "dec"),
                         labels = c("sep", "oct", "nov", "dec"))

# Count of number of group observations per species and month:
COUNT_OBS <- DATA_PRESENCE %>% 
                group_by(SPECIES, MONTH) %>% 
                summarise(n_OBS = sum(COUNT))

# Join data frames by shared attributes:
DATA_PRESENCE <- left_join(DATA_PRESENCE, COUNT_OBS, 
                           by = c("SPECIES" = "SPECIES", 
                                  "MONTH" = "MONTH"))

# Count of distinct grid cells with observations and average number of observations per cell:
COUNT_CELLS <- DATA_PRESENCE %>% 
                  group_by(SPECIES, MONTH, n_OBS) %>% 
                  summarise(n_CELLS = n_distinct(CELL_ID)) %>% 
                  mutate(OBS_CELL = n_OBS / n_CELLS)

# Join data frames by shared attributes:
DATA_PRESENCE <- left_join(DATA_PRESENCE, COUNT_CELLS, 
                           by = c("SPECIES" = "SPECIES", 
                                  "MONTH" = "MONTH",
                                  "n_OBS" = "n_OBS"))

# Check classes of behaviour data frame before joining with presence data frame:
str(DATA_BEHAVIOUR)

# Transforming variables into factors:
DATA_BEHAVIOUR <- DATA_BEHAVIOUR %>%
                    mutate(across(c(MONTH, SPECIES, DOM_BEH), 
                                    as.factor))

# Join data frames by shared attributes:
DATA_PRESENCE <- left_join(DATA_PRESENCE, DATA_BEHAVIOUR, 
                           by = c("CELL_ID" = "CELL_ID", 
                                  "MONTH" = "MONTH", 
                                  "SPECIES" = "SPECIES"))

# The labels of behaviour is changed from numbers to words. The original levels 
# are overwritten with new ones. 
# Show factor levels:
levels(DATA_PRESENCE$DOM_BEH)

# Change levels:
levels(DATA_PRESENCE$DOM_BEH) <- c("grazing", "browsing", "standing", "walking", 
                                   "lying", "wallowing", "other")

# Check new levels:
levels(DATA_PRESENCE$DOM_BEH)

# A new variable with overall behavioural classes is created based on the existing 
# variable DOM_BEH. 
# Checking the number of occurrences of each behaviour before compiling them:
(COUNT_BEH <- DATA_PRESENCE %>% 
                group_by(DOM_BEH) %>% 
                summarise(n_BEH = n()))

# "grazing" and "browsing" are compiled in FEEDING.
# "lying" and "wallowing" are compiled in LYING.
# "standing", "walking", and "other" are compiled in UPRIGHT.
# Function dplyr::mutate() is combined with function dplyr::case_when() to create 
# a column that will hold the behavioural class information. The DOM_BEH represents 
# the most dominant group behaviour of all observed in that cell. The grepl() 
# function looks for patterns in the data based on a character string search: 
DATA_PRESENCE <- DATA_PRESENCE %>%
                    mutate(BEH_CLASS = case_when(               
                      grepl("grazing", DOM_BEH) ~ "feeding",
                      grepl("browsing", DOM_BEH) ~ "feeding",
                      grepl("standing", DOM_BEH) ~ "upright",
                      grepl("walking", DOM_BEH) ~ "upright",
                      grepl("lying", DOM_BEH) ~ "lying",
                      grepl("wallowing", DOM_BEH) ~ "lying",
                      grepl("other", DOM_BEH) ~ "upright"))

# Rearranging the order of appearance:
DATA_PRESENCE$BEH_CLASS <- factor(DATA_PRESENCE$BEH_CLASS, 
                                  levels = c("feeding", "lying", "upright"),
                                  labels = c("feeding", "lying", "upright"))

# Counting number of each behavioural class and adding column with distribution
# of behavioural classes among species and months in percentage:
COUNT_BEH <- DATA_PRESENCE %>% 
  group_by(SPECIES, MONTH, BEH_CLASS) %>% 
  summarise(n = sum(COUNT)) %>% 
  mutate(prop = round(n / sum(n) * 100, 3))

# Create subsets for each species with function dplyr::filter():
BOS_TAURUS <- filter(DATA_PRESENCE, SPECIES == "bos_taurus")
BUBALUS_BUBALIS <- filter(DATA_PRESENCE, SPECIES == "bubalus_bubalis")

## DATA SUBSET - HOUR ----

# The data frame DATA_HOUR contains information on every group observation, whereas
# DATA_PRESENCE contains information on the group observations summarized for 
# each grid cell, species and month.

# Checking variable classes and transforming some:
str(DATA_HOUR)

DATA_HOUR <- DATA_HOUR %>%
                mutate(across(c(MONTH, SPECIES, BEHAVIOUR, HUMAN), 
                                as.factor))

DATA_HOUR$DATE <- as.Date(DATA_HOUR$DATE, 
                          "%m/%d/%y")

# Reordering variable MONTH: 
DATA_HOUR$MONTH <- factor(DATA_HOUR$MONTH, 
                          levels = c("sep", "oct", "nov", "dec"),
                          labels = c("sep", "oct", "nov", "dec"))

# Renaming variable with behaviour values:
DATA_HOUR <- rename(DATA_HOUR, DOM_BEH = BEHAVIOUR) 

# Checking levels of behaviour:
levels(DATA_HOUR$DOM_BEH)

# Changing the levels of behaviour from numbers to words:
levels(DATA_HOUR$DOM_BEH) <- c("grazing", "browsing", "standing", "walking", 
                               "lying", "rubbing", "wallowing", "other")

# Checking the new levels:
levels(DATA_HOUR$DOM_BEH)

# As in the previous data frame, a new variable with overall behavioural classes 
# is created based on the existing behaviour variable. In this case, the behaviour
# represents the dominant behaviour of the group for each observation. 
DATA_HOUR <- DATA_HOUR %>%
                mutate(BEH_CLASS = case_when(               
                  grepl("grazing", DOM_BEH) ~ "feeding",
                  grepl("browsing", DOM_BEH) ~ "feeding",
                  grepl("standing", DOM_BEH) ~ "upright",
                  grepl("walking", DOM_BEH) ~ "upright",
                  grepl("lying", DOM_BEH) ~ "lying",
                  grepl("rubbing", DOM_BEH) ~ "upright",
                  grepl("wallowing", DOM_BEH) ~ "lying",
                  grepl("other", DOM_BEH) ~ "upright"))

# Rearranging the order of appearance:
DATA_HOUR$BEH_CLASS <- factor(DATA_HOUR$BEH_CLASS, 
                                  levels = c("feeding", "lying", "upright"),
                                  labels = c("feeding", "lying", "upright"))

# Create subsets for each species with function dplyr::filter():
BOS_TAURUS_HOUR <- filter(DATA_HOUR, SPECIES == "bos_taurus")
BUBALUS_BUBALIS_HOUR <- filter(DATA_HOUR, SPECIES == "bubalus_bubalis")

## DATA SUBSET - HUMAN ----

# Create subset for data with presence of humans:
DATA_HUMAN <- filter(DATA_HOUR, HUMAN == 1)

# Check classes:
str(DATA_HUMAN)

# Reordering month: 
DATA_HUMAN$MONTH <- factor(DATA_HUMAN$MONTH, 
                             levels = c("sep", "oct", "nov", "dec"),
                             labels = c("sep", "oct", "nov", "dec"))


# DATA VISUALISATION ----

# Various plots are created with package ggplot2 to visualise the distribution 
# of the variables' values. The spread of the data might give an idea of the 
# outcome of further analysis.

## DEFINE PLOT FUNCTION AND COLOURS ----

# A ggplot2 function is costumized to make sure that all figures will be consistent 
# and have similar features.
# The new theme function is created using a predefined theme as a base:
plot_theme <- function(){
  theme_bw() +
    theme(axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          plot.title = element_text(size = 24, vjust = 1, hjust = 0.5, face = "bold"),
          legend.position = "bottom", legend.justification = c("right", "bottom"),
          legend.key.width = unit(2, "cm"), legend.key.height = unit(1, "cm"),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 11, face = "italic"),
          legend.box.background = element_rect(color = "grey", size = 0.3))
}


# Define color palettes for variables:

MONTH_col <- terrain.colors(4)

SPECIES_col <- c("rosybrown1", "#deebf7")

BEH_col <- brewer.pal(3, "Pastel2")

NDVI_col <- rev(brewer.pal(4, "Greens"))


## CLIMATE (Done and saved 14-11-2021) ----

# Temperature in Aarhus Municipality during the study period:
(TEMP_plot <- ggplot(DATA_CLIMATE, aes(x = DATE, y = TEMP_C, col = TEMP_KEY)) +
              geom_smooth(size = 1, method = "loess", se = F) +
              plot_theme() +
              scale_x_date(date_breaks = "2 weeks", date_labels = "%m/%d/%y") +
              theme(axis.title.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.text.x = element_blank(),
                    legend.direction = "vertical",
                    legend.position = c(0.99,0.75),
                    legend.title = element_blank()) +
              labs(y = "Temperature [°C]\n") +
              facet_grid(~ MONTH, scales = "free"))

# Precipitation in Aarhus Municipality during the study period:
(PREC_plot <- ggplot(DATA_CLIMATE, aes(x = DATE, y = PREC_MM, fill = PREC_MM)) +
              geom_bar(stat = "identity", size = 0.5) +
              plot_theme() +
              scale_x_date(date_breaks = "2 weeks", date_labels = "%m/%d/%y") +
              theme(legend.position = c(0.99,0.85),
                    legend.direction = "horizontal",
                    legend.key.width = unit(1, "cm"),
                    legend.key.height = unit(0.7, "cm"),
                    legend.title = element_blank(),
                    axis.title.x = element_blank()) +
              scale_fill_distiller(palette = "Blues", direction = 1) +
              labs(y = "Precipitation [mm]\n") +
              facet_grid(~ MONTH, scales = "free"))

# Arrange the elements to be plotted. The inner arrangeGrob() function arranges 
# the plots, the main title, and the global y-axis title. The outer grid.arrange() 
# function arranges and draws the arrangeGrob object and the legend.
(CLIMATE_panel <- grid.arrange(arrangeGrob(
                  TEMP_plot + labs(tag = "(a)"),
                  PREC_plot + labs(tag = "(b)"), 
                  nrow = 2,
                  top = textGrob("Daily climate data partitioned by month\n",  
                      gp = gpar(fontface = "bold", cex = 2.5)),
                  bottom = textGrob("Date [M/D/Y]", 
                      gp = gpar(fontface = "bold", cex = 1.5)))))

# Save plot:
#ggsave(CLIMATE_panel, file = "FIGURES/climate_panel.png", width = 14, height = 14)


## NDVI (Done and saved 14-11-2021) ----

# Frequency distribution of NDVI in Kasted Mose during study period:
(NDVI_freq <- ggplot(DATA_ALL, aes(x = NDVI, y = (stat(count) / sum(count)) * 100)) +
              geom_histogram(aes(fill = ..x..), binwidth = 0.05, color = "black") +
              scale_fill_distiller(palette = "Greens", direction = 1) +
              plot_theme() +
              theme(legend.position = "none") +
              geom_vline(aes(xintercept = mean(NDVI)), 
                         linetype = "dashed", size = 1) +
              labs(y = "Relative frequency [%]"))

# Save plot:
#ggsave(NDVI_freq, file = "FIGURES/ndvi_freq.png", width = 14, height = 7)

# Link of icon from phylopic website:
leaf <- "http://phylopic.org/assets/images/submissions/25e1de47-8f26-4ae0-beda-44bbdc079934.512.png"
# Load logo into R:
leaf_logo <- readPNG(getURLContent(leaf))

# Monthly variation in NDVI in all cells:
(NDVI_box_full <- ggplot(DATA_ALL, aes(x = MONTH, y = NDVI)) +
                  geom_boxplot(aes(fill = MONTH)) +
                  plot_theme() +
                  stat_compare_means(ref.group = ".all.", 
                                     label = "p.signif", 
                                     label.y = c(0.78, 0.7, 0.7, 0.6), 
                                     method = "wilcox.test") +
                  coord_cartesian(ylim = c(0,0.8)) +
                  theme(legend.position = "none",
                        axis.title.x = element_blank()) +
                  geom_hline(aes(yintercept = mean(NDVI)), 
                             linetype = "dashed", size = 1) +
                  scale_fill_manual(values = NDVI_col) +
                  add_phylopic(leaf_logo, alpha = 0.6, x = 4, y = 0.75, 
                               ysize = 0.8, color = "orange"))

# Compare mean NDVI values between the months: 
compare_means(NDVI ~ MONTH, data = DATA_ALL, 
              method = "wilcox.test")


NDVI_mean <- ddply(DATA_PRESENCE, "SPECIES", summarise, grp.mean = mean(NDVI))

# Monthly variation in NDVI in cells with observations:
(NDVI_box_subset <- ggplot(DATA_PRESENCE, aes(x = MONTH, y = NDVI)) +
                    geom_boxplot(aes(fill = SPECIES)) +
                    plot_theme() +
                    stat_compare_means(ref.group = ".all.", 
                                       label = "p.signif",
                                       method = "wilcox.test",
                                       label.y = c(0.78, 0.7, 0.7, 0.6)) +
                    coord_cartesian(ylim = c(0,0.8)) +
                    geom_hline(data = NDVI_mean, aes(yintercept = grp.mean), 
                               linetype = "dashed", color = c("red", "blue"), size = 1) +
                    theme(legend.position = c(0.99,0.85),
                          axis.title.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.text.y = element_blank(),
                          axis.title.x = element_blank()) +
                    scale_fill_manual(values = SPECIES_col, name = "Species"))

# Compare inter-specific mean NDVI values in cells with occurrences:
compare_means(NDVI ~ SPECIES, data = DATA_PRESENCE, 
              group.by = "MONTH", method = "wilcox.test")

# Compare intra-specific mean NDVI values in cells with occurrence between months
# (each month is compared to all the others at once):
compare_means(NDVI ~ MONTH, data = DATA_PRESENCE, 
              group.by = "SPECIES", 
              ref.group = ".all.", 
              method = "wilcox.test")

# Combine plots:
(NDVI_panel <- grid.arrange(arrangeGrob(
                 NDVI_box_full + labs(tag = "(a)"),
                 NDVI_box_subset + labs(tag = "(b)"), 
                 ncol = 2,
                 bottom = textGrob("Month", 
                    gp = gpar(fontface = "bold", cex = 1.5)))))

# Save plot:
#ggsave(NDVI_panel, file = "FIGURES/ndvi_panel.png", width = 14, height = 10)

# mean(DATA_ALL$NDVI) # 0.4524683
# mean(DATA_PRESENCE$NDVI) # 0.4754818
# bos_taurus mean = 0.4929637
# bubalus_bubalis mean = 0.4563848


## TWI (Done and saved 15-11-2021) ----

# Frequency distribution of NDVI in Kasted Mose during study period:
(TWI_freq <- ggplot(DATA_ALL, aes(x = TWI, y = (stat(count) / sum(count)) * 100)) +
             geom_histogram(aes(fill = ..x..), binwidth = 0.25, color = "black") +
             plot_theme() +
             coord_cartesian(xlim = c(2,10)) +
             theme(axis.title.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_blank(),
                   legend.position = "none") +
             scale_fill_distiller(palette = "GnBu", direction = 1,
                                  limits = c(0,10)) +
             geom_vline(aes(xintercept = mean(TWI)), 
                        linetype = "dashed", size = 1) +
             labs(y = "Relative frequency [%]\n"))

# Reordering MONTH values to appear correct on the y-axis in the subsequent barplot:
DATA_PRESENCE$MONTH <- factor(DATA_PRESENCE$MONTH, 
                            levels = rev(c("sep", "oct", "nov", "dec")),
                            labels = rev(c("sep", "oct", "nov", "dec")))

TWI_mean <- ddply(DATA_PRESENCE, "SPECIES", summarise, grp.mean = mean(TWI))

# TWI in cells with observations:
(TWI_box <- ggplot(DATA_PRESENCE, aes(TWI, MONTH)) +
            geom_boxplot(aes(fill = SPECIES)) +
            plot_theme() +
            coord_cartesian(xlim = c(2,10)) +
            scale_fill_manual(values = SPECIES_col, name = "Species") +
            geom_vline(data = TWI_mean, aes(xintercept = grp.mean), 
                       linetype = "dashed", color = c("red", "blue"), size = 1) +
            theme(legend.key.height = unit(0.7, "cm"), 
                  legend.key.width = unit(1.5, "cm"),
                  legend.position = c(0.99,0.82),
                  legend.direction = "vertical",
                  axis.title.x = element_blank()) +
            labs(y = "Month\n"))

# Compare inter-specific mean TWI values in cells with occurrences: 
compare_means(TWI ~ SPECIES, data = DATA_PRESENCE, 
              group.by = "MONTH", 
              method = "wilcox.test")

# Compare intra-specific mean TWI values in cells with occurrences between months:
compare_means(TWI ~ MONTH, data = DATA_PRESENCE, 
              group.by = "SPECIES", 
              method = "wilcox.test")

# Compare intra-specific mean TWI values for each behaviour:
compare_means(TWI ~ BEH_CLASS, data = DATA_PRESENCE, 
              group.by = "SPECIES", 
              ref.groups = ".all.", 
              method = "wilcox.test")

# Combine plots:
(TWI_panel <- grid.arrange(arrangeGrob(
              TWI_freq + labs(tag = "(a)"),
              TWI_box + labs(tag = "(b)"),
              nrow = 2,
              bottom = textGrob("TWI",
                            gp = gpar(fontface = "bold", cex = 1.5)))))

# Save plot:
#ggsave(TWI_panel, file = "FIGURES/twi_panel.png", width = 14, height = 14)

# mean(DATA_ALL$TWI) # 5.265252
# mean(DATA_PRESENCE$TWI) # 5.339236
# bos_taurus mean = 5.362337
# bubalus_bubalis mean = 5.314000

## GRID CELLS (Done and saved 14-11-2021) ----

# Number of occurrences in each cell
N_CELLS <- DATA_ALL %>% 
  group_by(CELL_ID, PRESENCE) %>% 
  summarise(n_OBS = n_distinct(MONTH))

N_CELLS <- N_CELLS %>% 
  group_by(PRESENCE, n_OBS) %>% 
  summarise(n = n())

# Reordering MONTH values to appear correct on the x-axis in the subsequent barplot:
DATA_PRESENCE$MONTH <- factor(DATA_PRESENCE$MONTH, 
                              levels = c("sep", "oct", "nov", "dec"),
                              labels = c("sep", "oct", "nov", "dec"))

(OBS_CELLS <- ggplot(DATA_PRESENCE, aes(x = MONTH, y = OBS_CELL)) + 
              geom_bar(aes(fill = SPECIES), stat = "identity", color = "black", 
                       position = "dodge", width = 0.5) +
              plot_theme() +
              stat_compare_means(ref.group = ".all.",
                                 label = "p.signif",
                                 method = "wilcox.test",
                                 label.y = c(12.5, 10.5, 12.0, 10.0)) +
              scale_fill_manual(values = SPECIES_col, name = "Species") +
              theme(legend.direction = "vertical",
                    legend.position = c(0.99,0.85)) +
              labs(x = "Month", y = "Observations", 
                   title = "Average number of observations per grid cell\n"))

# Save plot:
#ggsave(OBS_CELLS, file = "FIGURES/obs_cells.png", width = 7, height = 7)


## BEHAVIOUR (Done and saved 14-11-2021) ----

# Compare intra-specific proportion performed between behaviours: 
compare_means(prop ~ BEH_CLASS, group.by = "SPECIES", 
              data = COUNT_BEH, method = "wilcox.test")

# Compare inter-specific proportion performed between behaviours:
compare_means(prop ~ SPECIES, group.by = "BEH_CLASS", 
              data = COUNT_BEH, method = "wilcox.test")

# Pie plot of proportional occurrence of behaviours:
(BEH_pie <- ggplot(COUNT_BEH, aes(x = "", y = prop, fill = BEH_CLASS)) +
            geom_col(color = "black", position = "fill") +
            theme_light(base_size = 20) +
            theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  legend.key.height = unit(2, "cm"), 
                  legend.key.width = unit(1, "cm"), 
                  legend.justification = c("right", "bottom"),
                  legend.direction = "vertical",
                  legend.title = element_text(face = "bold")) +
            scale_fill_manual(values = BEH_col, name = "Behaviour") +
            coord_polar(theta = "y") +
            facet_grid(SPECIES ~ MONTH))

# Save plot:
#ggsave(BEH_pie, file = "FIGURES/beh_pie.png", width = 14, height = 7)

# Diurnal behavioural dynamics:
(BEH_violin <- ggplot(DATA_HOUR, aes(x = BEH_CLASS, y = HOUR, fill = SPECIES)) + 
               geom_violin(trim = F) +
               plot_theme() + 
               theme(legend.position = "right",
                     legend.key.height = unit(2, "cm"),
                     legend.key.width = unit(1, "cm"),
                     axis.title.x = element_blank()) +
               scale_fill_manual(values = SPECIES_col, name = "Species") +
               coord_cartesian(ylim = c(5,18)) +
               scale_y_continuous(breaks = seq(6,18,2)) +
               labs(y = "Hour"))

# Save plot:
#ggsave(BEH_violin, file = "FIGURES/beh_violin.png", width = 14, height = 7)


## HUMAN-ANIMAL RELATIONSHIP ----

# Spatial avoidance

PATH_mean <- ddply(DATA_PRESENCE, "SPECIES", summarise, grp.mean = mean(PATH_DIST))

# Plot of distance to walking paths:
(PATH_box <- ggplot(DATA_PRESENCE, aes(x = SPECIES, y = PATH_DIST)) + 
         geom_boxplot(aes(fill = SPECIES), position = "dodge") +
         plot_theme() +
         geom_hline(data = PATH_mean, aes(yintercept = grp.mean), 
               linetype = "dashed", color = c("red", "blue"), size = 1) +
         scale_fill_manual(values = SPECIES_col, name = "Species") +
         theme(legend.direction = "vertical",
               legend.position = c(0.99,0.85),
               axis.title.x = element_blank(),
               axis.text.x = element_blank()) +
         labs(y = "Distance [m]"))

# Mean distance:
# Overall = 115.8225m
# Bos taurus = 124.2369m
# Bubalus bubalis = 106.6308m

# Save plot:
#ggsave(PATH_box, file = "FIGURES/path_box.png", width = 5, height = 10)

# Frequency of visits by hour:
COUNT_HUMAN_HOUR <- DATA_HUMAN %>% 
  group_by(HOUR) %>% 
  summarise(n = n()) %>% 
  mutate(freq = round(n / sum(n) * 100,3))

# Frequency of visits by month:
COUNT_HUMAN_MONTH <- DATA_HUMAN %>% 
  group_by(MONTH) %>% 
  summarise(n = n()) %>% 
  mutate(freq = round(n / sum(n) * 100,3))

# Frequency distribution of human visits
(HUMAN_HOUR <- ggplot(COUNT_HUMAN_HOUR, aes(x = HOUR, y = freq)) +
    geom_col(aes(fill = ..y..), width = 0.8, color = "black", alpha = 0.7) +
    scale_fill_gradient2(low = "black", mid = "yellow", high = "red") +
    plot_theme() +
    theme(legend.position = "none") +
    scale_x_continuous(limits = c(7.5,16.5), breaks = seq(8,16,2)) +
    labs(x = "Hour", y = "Relative frequency [%]"))

# Save plot:
#ggsave(HUMAN_HOUR, file = "O:/Nat_Ecoinformatics/C_Write/_User/KarolineDohrmann_au540949/DATA/DATA_12102021/FIGURES/human_hour.png", width = 10, height = 5)

COUNT_HUMAN_BEH <- DATA_HOUR %>% 
  group_by(SPECIES, HUMAN, BEH_CLASS) %>% 
  summarise(n = n()) %>% 
  mutate(freq = round(n / sum(n) * 100,3))

# Frequency distribution of behaviour
(HUMAN_BEH <- ggplot(DATA_HOUR, aes(x = HUMAN, fill = BEH_CLASS)) + 
             geom_bar(position = "fill", color = "black") +
             plot_theme() +
             scale_fill_manual(values = BEH_col, name = "Behaviour") +
             theme(legend.direction = "vertical",
                   legend.position = "right") +
             labs(x = "Human Presence", y = "Fraction"))

# Save plot:
#ggsave(HUMAN_BEH, file = O:/Nat_Ecoinformatics/C_Write/_User/KarolineDohrmann_au540949/DATA/DATA_12102021/FIGURES/human_beh.png", width = 10, height = 12)

# Combine plots:
(HUMAN_panel <- grid.arrange(arrangeGrob(
                PATH_box + labs(tag = "(a)"),
                HUMAN_BEH + labs(tag = "(b)"),
                nrow = 1)))
                
# Save plot:
#ggsave(HUMAN_panel, file = "FIGURES/human_panel.png", width = 14, height = 10)


# MODELLING ----

## Model 1 - Habitat selection (total area) ----

# The full dataset is tested with the binary variable PRESENCE as response variable. 
# It is tested whether NDVI or TWI has an influence on the choice of habitat. 
# Furthermore, it is tested whether there is an inter-specific difference,
# as well as if there is any interactions. Cell ID is used as a random factor to 
# account for non-independence. A generalized linear mixed-effects model with the 
# family binomial is chosen. It is a logistic regression and binomial is chosen 
# due to the dichotomous nature of the response variable (0/1). Different models 
# are tested and their AIC values are compared to find the best fitted model. 
# The best model from the set of plausible models being considered is the 
# one with the smallest AIC value (the least information loss relative 
# to the true model).

# Distribution of random effect:
(cell_hist_all <- ggplot(DATA_ALL, aes(CELL_ID)) +
                  geom_histogram())

(cell_qq_all <- ggplot(DATA_ALL, aes(sample = CELL_ID)) +
                stat_qq() +
                stat_qq_line(color = "steelblue", lwd = 1) +
                coord_cartesian(ylim = c(0,750)))

(cell_panel_all <- grid.arrange(arrangeGrob(
                   cell_hist_all + labs(tag = "(a)"),
                   cell_qq_all + labs(tag = "(b)"),
                   nrow = 1)))
                   
# The histogram shows a light tail uniform distribution.
# The light tailed distributions yield an s shape depicted in the qq plot.

# Save plot:
#ggsave(cell_panel_all, file = "FIGURES/cell_panel_all.png", width = 10, height = 5)

# Distribution of response variable:
ggplot(DATA_ALL, aes(PRESENCE)) +
       geom_bar(stat = "count")
# Very imbalanced, with a lot more 0's than 1's.

# Full model:
M1_1 <- glmer(PRESENCE ~ TWI*SPECIES + NDVI*SPECIES + (1|CELL_ID), 
              data = DATA_ALL, family = binomial)

summary(M1_1) # AIC 1965.2
drop1(M1_1, test = "Chi")


# Removing non-significant interaction between SPECIES and TWI and refitting model:
M1_2 <- glmer(PRESENCE ~ TWI + NDVI*SPECIES + (1|CELL_ID), 
              data = DATA_ALL, family = binomial)

summary(M1_2) # AIC 1963.7
drop1(M1_2, test = "Chi")

# TWI is nearly non-significant and is removed. Refitting the model:
M1_3 <- glmer(PRESENCE ~ NDVI*SPECIES + (1|CELL_ID), 
              data = DATA_ALL, family = binomial)

summary(M1_3) # AIC 1965.4

# The AIC was lower when TWI was included.

# Comparing models
anova(M1_2, M1_3)


# There is a tendency that TWI has an effect, and the model with TWI
# has a lower AIC - therefore, this model is chosen as the most optimal model.
M1 <- glmer(PRESENCE ~ TWI + NDVI*SPECIES + (1|CELL_ID), 
            data = DATA_ALL, family = binomial)

summary(M1)
(resid_M1 <- resid_panel(M1))

# Save plot:
#ggsave(resid_M1, file = "FIGURES/resid_M1.png", width = 10, height = 10)


# From the Q-Q plot and histogram of residuals it is evident that data is right-skewed.
# Number of each level in PRESENCE:
table(DATA_ALL$PRESENCE)


## Model 2 - Habitat selection (subset area) ----

# Distribution of random effect:
(cell_hist_sub <- ggplot(DATA_PRESENCE, aes(CELL_ID)) +
                  geom_histogram(binwidth = 50))

(cell_qq_sub <- ggplot(DATA_PRESENCE, aes(sample = CELL_ID)) +
                stat_qq() +
                stat_qq_line(color = "steelblue", lwd = 1) +
                coord_cartesian(ylim = c(0,750)))

# Combine plots:
(cell_panel_sub <- grid.arrange(arrangeGrob(
                   cell_hist_sub + labs(tag = "(a)"),
                   cell_qq_sub + labs(tag = "(b)"),
                   nrow = 1)))

# Save plot:
#ggsave(cell_panel_sub, file = "FIGURES/cell_panel_sub.png", width = 10, height = 5)


# The response variable, PROP, represents the relative frequency of 
# observations in a cell:
(prop_hist <- ggplot(DATA_PRESENCE, aes(PROP)) +
   geom_histogram())

(prop_qq <- ggplot(DATA_PRESENCE, aes(sample = PROP)) +
    stat_qq() +
    stat_qq_line(color = "steelblue", lwd = 1))
# Right skewed

# Square-root transformation:
(prop_hist2 <- ggplot(DATA_PRESENCE, aes(sqrt(PROP))) +
    geom_histogram())

(prop_qq2 <- ggplot(DATA_PRESENCE, aes(sample = sqrt(PROP))) +
    stat_qq() +
    stat_qq_line(color = "steelblue", lwd = 1))
# It's better.

# Log-transformation:
(prop_hist3 <- ggplot(DATA_PRESENCE, aes(log(PROP))) +
    geom_histogram())

(prop_qq3 <- ggplot(DATA_PRESENCE, aes(sample = log(PROP))) +
    stat_qq() +
    stat_qq_line(color = "steelblue", lwd = 1))

# Combine plots:
prop_panel <- grid.arrange(prop_hist + labs(tag = "(a)"), prop_qq,
              prop_hist2 + labs(tag = "(b)"), prop_qq2,
              prop_hist3 + labs(tag = "(c)"), prop_qq3,
              ncol = 2, nrow = 3,
              layout_matrix = rbind(c(1,2), c(3,4), c(5,6)))


# Save plot:
#ggsave(prop_panel, file = "FIGURES/prop_panel.png", width = 12, height = 14)

# The data is right-skewed with most cells having few observations and a few cells
# having many observations, as much as 20% of the observations in one cell that month. 
# This is evident for both species and all months.

# Full model fit with ML:
M2_1 <- lmer(sqrt(PROP) ~ NDVI*SPECIES + TWI*SPECIES + (1|CELL_ID), 
             data = DATA_PRESENCE, REML = F)

summary(M2_1) # AIC -885.7
drop1(M2_1, test = "Chi")

# Removes non-significant interaction between SPECIES and NDVI and refit the model:
M2_2 <- lmer(sqrt(PROP) ~ NDVI + TWI*SPECIES + (1|CELL_ID), 
             data = DATA_PRESENCE, REML = F)

summary(M2_2) # -886.6
drop1(M2_2)

# Remove NDVI and refit model:
M2_3 <- lmer(sqrt(PROP) ~ TWI*SPECIES + (1|CELL_ID), 
             data = DATA_PRESENCE, REML = F)

summary(M2_3) # AIC -888.5

# Removing interaction and refitting model:
M2_4 <- lmer(sqrt(PROP) ~ TWI + SPECIES + (1|CELL_ID), 
              data = DATA_PRESENCE, REML = F)
summary(M2_4) # AIC -887.1

#Removing SPECIES and refitting model:
M2_5 <- lmer(sqrt(PROP) ~ TWI + (1|CELL_ID), 
              data = DATA_PRESENCE, REML = F)
summary(M2_5) # AIC -889.1

# Comparing models:
anova(M2_5, M2_3) # M2_5 best
anova(M2_5, M2_4) # M2_5 best

# The inclusion of SPECIES does not lead to a significant improvement of the model.
# The model is refitted with REML
M2 <- lmer(sqrt(PROP) ~ TWI + (1|CELL_ID), 
             data = DATA_PRESENCE, REML = T)
summary(M2)

(resid_M2 <- resid_panel(M2))

# A bit right skewed and heteroscadastic pattern. 

# Save plot:
#ggsave(resid_M2, file = "FIGURES/resid_M2.png", width = 10, height = 10)


## Model 3 - Diurnal activity budget ----

# Using an ANOVA to test if number of cells and behaviour change over time.

# From the summary of the model we can see that relative frequency of presence 
# in a cell varies significantly based on TWI and SPECIES.

# The behaviours are highly significant in predicting where the animals occur. 
# The months and species does not significantly predict where the animals occur. 
# There is a tendency of behaviour depending on month and species in predicting 
# where the animals occur.

# One-way ANOVA:

# Is there a difference in occurrence based on behaviour class?
M3_1 <- aov(sqrt(PROP) ~ BEH_CLASS, data = DATA_PRESENCE)
summary(M3_1)
# P < 0.001 

# Calculate post-hoc test:
TukeyHSD(M3_1)

# There is not a difference in areas selected for FEEDING or UPRIGHT behaviour.
# There is a difference in areas selected for LYING and FEEDING or UPRIGHT, respectively.

# Is there a difference in occurrence based on month of observation?
M3_2 <- aov(sqrt(PROP) ~ MONTH, data = DATA_PRESENCE)
summary(M3_2)
# P = 0.179

# Calculate post-hoc test:
TukeyHSD(M3_2) # None of the interactions are significant.

# Is there a difference on occurrence based on species?
M3_3 <- aov(sqrt(PROP) ~ SPECIES, data = DATA_PRESENCE)
summary(M3_3)
# P = 0.938

# Two-way ANOVA with interaction term:

# Does the effect of behaviour on occurrence depend on month of observation?
M3_4 <- aov(sqrt(PROP) ~ BEH_CLASS*MONTH, data = DATA_PRESENCE)
summary(M3_4)
# BEH P < 0.001
# MONTH P = 0.252
# BEH*MONTH P < 0.1

# Calculate post-hoc test:
TukeyHSD(M3_4)

# Comparing models:
anova(M3_1,M3_4) # The inclusion of month improves the model slightly.

# Does the effect of behaviour on occurrence depend on species?
M3_5 <- aov(sqrt(PROP) ~ BEH_CLASS*SPECIES, data = DATA_PRESENCE)
summary(M3_5)
# BEH P < 0.005
# SPECIES P = 0.426
# BEH*SPECIES P < 0.1

# Calculate post-hoc test:
TukeyHSD(M3_5)

# Comparing models:
anova(M3_1,M3_5) # The inclusion of species does not improve the model.

# When comparing the main effect of Behaviour and the two interaction models,
# the appliance of Month adds a tendency of significant explanatory power
# to the model, while species does not.

# Final model and post-hoc test:
M3 <- aov(sqrt(PROP) ~ BEH_CLASS*MONTH, data = DATA_PRESENCE)
summary(M3)
(tuk <- TukeyHSD(M3))
plot(tuk, col = "red")


## Model 4 - Human-animal relationship (spatial avoidance) ----

# Distribution of the variable PATH_DIST:
(path_hist1 <- ggplot(DATA_PRESENCE, aes(PATH_DIST)) +
               geom_histogram())

(path_qq1 <- ggplot(DATA_PRESENCE, aes(sample = PATH_DIST)) +
    stat_qq() +
    stat_qq_line(color = "steelblue", lwd = 1))

# Square-root transforming as it is right-skewed:
(path_hist2 <- ggplot(DATA_PRESENCE, aes(sqrt(PATH_DIST))) +
    geom_histogram())

(path_qq2 <- ggplot(DATA_PRESENCE, aes(sample = sqrt(PATH_DIST))) +
    stat_qq() +
    stat_qq_line(color = "steelblue", lwd = 1))

# Combine plots:
path_panel <- grid.arrange(path_hist1 + labs(tag = "(a)"), path_qq1,
                           path_hist2 + labs(tag = "(b)"), path_qq2,
                           ncol = 2, nrow = 2,
                           layout_matrix = rbind(c(1,2), c(3,4)))

# Save plot:
ggsave(path_panel, file = "FIGURES/path_panel.png", width = 10, height = 12)

# Full model:
M4_1 <- lmer(sqrt(PROP) ~ sqrt(PATH_DIST)*BEH_CLASS + sqrt(PATH_DIST)*SPECIES + 
                          (1|CELL_ID), data = DATA_PRESENCE, REML = F)

summary(M4_1) # AIC -1049.1
drop1(M4_1)

# Removing interaction between PATH_DIST and SPECIES and refitting model: 
M4_2 <- lmer(sqrt(PROP) ~ sqrt(PATH_DIST)*BEH_CLASS + SPECIES + (1|CELL_ID), 
                          data = DATA_PRESENCE, REML = F)

summary(M4_2) # AIC -1051.1
drop1(M4_2)

# Removing species as main term and refitting modelL:
M4_3 <- lmer(sqrt(PROP) ~ sqrt(PATH_DIST)*BEH_CLASS + (1|CELL_ID), 
             data = DATA_PRESENCE, REML = F)

summary(M4_3) # AIC -1052.0

# Refitting final model with REML:
M4 <- lmer(sqrt(PROP) ~ sqrt(PATH_DIST)*BEH_CLASS + (1|CELL_ID), 
             data = DATA_PRESENCE, REML = T)

summary(M4)

(resid_M4 <- resid_panel(M4))

# Save plot:
#ggsave(resid_M4, file = "FIGURES/resid_M4.png", width = 10, height = 10)


## Model 5 - Human-animal relationship (behaviour) ----

# Frequency table
hum_beh <- table(DATA_HOUR$HUMAN, DATA_HOUR$BEH_CLASS)
hum_beh

# Proportion table by row
prop.table(hum_beh, 1)

# Chi-square test
(xsq <- chisq.test(hum_beh))
xsq$observed
xsq$expected
xsq$residuals
xsq$stdres

# The null hypothesis is rejected as there is a difference between expected and 
# observed frequencies, and behaviour thus is associated with human presence.

chisq.posthoc.test(hum_beh)

# The difference in behaviour class UPRIGHT is significant.
