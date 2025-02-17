library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(lme4)
library(lmtest)
library(lmerTest)
library(sjPlot)
library(readxl)
library(tidyverse)
library(hrbrthemes)
library(forcats)
library(emmeans)
library(multcompView)

##data loading for root architecture
roots <- read_excel("rootdata_analysis_cleaneddata.xlsx", 
                    col_types = c("numeric", "text", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric"))
roots$Length <- (roots$Length)* (25.4/2400) 
roots$AvgDiamMeasured <- (roots$AvgDiamMeasured)*(25.4/2400)/10

##calculations for specific root length -- SRL
roots$SRL <- roots$Length/roots$root_dry_mass

##calculations for BRANCHING INTENSITY -- tip number and forks
roots$branching_intensity_tips <- roots$Tips/roots$Length

#data loading for mycorrhizal growth response  MGR -- just keeping essential columns

Ericoid <- read_excel("Ericoid.xls", col_types = c("numeric", 
                                                   "text", "text", "numeric", "text", "text", 
                                                   "numeric", "text", "text", "numeric", 
                                                   "numeric", "numeric", "numeric", "numeric"))


Ericoid <- Ericoid %>% select(-c("blwg dry BEF", "matches original set up", "blwg wet (Lia)", "notes_1", "notes_2","Dead", "replanted")) %>%
  rename(
    aboveground_dry_biomass = "abvg dry biomass",
    belowground_dry_biomass = "blwg dry (Lia)"
  )

#data loading for fungal information
fungalcomb <- read_excel("fungalcomb.xls", 
                         col_types = c("numeric", "numeric", "text", 
                                       "numeric"))
Ericoid <- Ericoid %>% left_join(fungalcomb, by = "meso_num")

#Replace NA in f_comb column with "STR"
Ericoid$f_comb <- ifelse(Ericoid$f_comb == "NA", NA, Ericoid$f_comb)
Ericoid$f_comb <- ifelse(is.na(Ericoid$f_comb), "STR", Ericoid$f_comb)


# Filter data for STR and calculate average above-ground biomass per species
plant_nam <- c("VAN", "GMI", "GSH", "GPR", "CVU", "VMY", "VCO", "VVI", "PJA", "KLA", "VMA")
STR <- Ericoid%>%
  filter(f_comb == "STR") %>%
  select(plant_sp.x, aboveground_dry_biomass)

avr_STR <- sapply(plant_nam, function(species) {
  mean(STR$aboveground_dry_biomass[STR$plant_sp.x == species], na.rm = TRUE)
})

#table of average above-ground biomass per species
AB_biomass_STR_table <- data.frame(plant_sp.x = plant_nam, AB_biomass_STR = avr_STR)

#adding avg-table to the main data-set: Ericoid 
Ericoid <- Ericoid %>%
  left_join(AB_biomass_STR_table, by = "plant_sp.x")

# Calculate MGR_abvg while avoiding division by zero
Ericoid <- Ericoid %>%
  mutate(MGR_abvg = ifelse(
    is.na(AB_biomass_STR) | AB_biomass_STR == 0, NA,
    (aboveground_dry_biomass - AB_biomass_STR) / AB_biomass_STR
  ))

################################################################################
################################################################################
################################################################################
############################# SPECIFIC ROOT LENGTH #############################
################################################################################
################################################################################
################################################################################

#join root architecture data and MGR data (meso_num and plant_id as common keys)
Ericoid <- Ericoid %>%
  left_join(roots, by = c("plant_id", "meso_num"), relationship = "many-to-many")

# delete useless data
Ericoid <- Ericoid[Ericoid$meso_num %in% c(1:33, 124:156), ]
Ericoid <- Ericoid %>% filter(!is.na(aboveground_dry_biomass))

#Data subset of architecture of sterile plants -- mean of SRL of sterile plants grouped by plant species
sterile_architecture <- subset(Ericoid, f_comb == "STR")

mean_str_arc <- sterile_architecture %>%
  group_by(plant_sp.x) %>% 
  summarise(mean_SRL = mean(SRL, na.rm = TRUE)) 

mean_str_arc <- na.omit(mean_str_arc)

#Data subset of inoculated plants 
inoculated_data <- Ericoid[!is.na(Ericoid$MGR_abvg), ]
inoculated_data <- inoculated_data[inoculated_data$f_comb != "STR", ]

#Group by meso_num and fungal combination (f_comb) to calculate the mean MGR
mean_inoc_mgr_pot <- inoculated_data %>%
  group_by(meso_num, plant_sp.x, f_comb) %>%
  summarise(mean_MGR = mean(MGR_abvg, na.rm = TRUE), .groups = "drop")

#Calculate mean initial height for each meso_num in the Ericoid table
mean_initial_height <- Ericoid %>%
  group_by(meso_num) %>%
  summarise(mean_initial_height = mean(height_i_cm, na.rm = TRUE))

#merge subsets
# Merge the datasets
merged_means_str_arc_inoc_mgr_height_i <- merge(mean_str_arc, mean_inoc_mgr_pot, by = "plant_sp.x")
merged_means_str_arc_inoc_mgr_height_i <- merge(merged_means_str_arc_inoc_mgr_height_i, mean_initial_height, by = "meso_num")

### MODELS
#linear regression model
#lm1 <- lm(mean_MGR ~ mean_SRL * f_comb + mean_initial_height, data = merged_means_str_arc_inoc_mgr_height_i)
#summary(lm1)
#plot(lm1)
#anova(lm1)

# corrected linear regression model 
lm2 <-  lm(mean_MGR ~ mean_SRL + mean_initial_height, data = merged_means_str_arc_inoc_mgr_height_i)
summary(lm2)
anova(lm2)


prediction_data <- data.frame(
  mean_SRL = seq(min(merged_means_str_arc_inoc_mgr_height_i$mean_SRL),
                 max(merged_means_str_arc_inoc_mgr_height_i$mean_SRL),
                 length.out = 100),
  mean_initial_height = mean(merged_means_str_arc_inoc_mgr_height_i$mean_initial_height)
)
prediction_data$predicted_MGR <- predict(lm2, newdata = prediction_data)

dev.off()  
ggplot(data = merged_means_str_arc_inoc_mgr_height_i, aes(x = mean_SRL, y = mean_MGR, color = f_comb)) +
  geom_point(size = 3, alpha = 0.7) +  # Observed data points
  geom_line(data = prediction_data, aes(x = mean_SRL, y = predicted_MGR), color = "grey45", linewidth = 0.8) +  # Regression line
  labs(
    title = "Linear Regression: MGR ~ SRL + Initial Height",
    x = "Specific Root Length",
    y = "Mycorrhizal Growth Response",
    color = "Fungal Species"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           x = unique(merged_means_str_arc_inoc_mgr_height_i$mean_SRL), 
           y = min(merged_means_str_arc_inoc_mgr_height_i$mean_MGR) - 0.5,  # Adjust to position below plot
           label = c(
             "VAN",
             "GSH",
             "GPR",
             "CVU",
             "VMY",
             "VCO",
             "VVI",
             "PJA",
             "KLA",
             "VMA"
           ),
           angle = 90, 
           hjust = 0, 
           size = 3)

################################################################################
################################################################################
################################################################################
############################# BRANCHING INTENSITY###############################
################################################################################
################################################################################
################################################################################


#Data subset of architecture of sterile plants -- mean of BI of sterile plants grouped by plant species
sterile_architecture <- subset(Ericoid, f_comb == "STR")

mean_str_arc <- sterile_architecture %>%
  group_by(plant_sp.x) %>% 
  summarise(mean_BIT = mean(branching_intensity_tips, na.rm = TRUE)) 

mean_str_arc <- na.omit(mean_str_arc)
#view(mean_str_arc)

#Data subset of inoculated plants 
inoculated_data <- Ericoid[!is.na(Ericoid$MGR_abvg), ]
inoculated_data <- inoculated_data[inoculated_data$f_comb != "STR", ]

#Group by meso_num and fungal combination (f_comb) to calculate the mean MGR
mean_inoc_mgr_pot <- inoculated_data %>%
  group_by(meso_num, plant_sp.x, f_comb) %>%
  summarise(mean_MGR = mean(MGR_abvg, na.rm = TRUE), .groups = "drop")

#Calculate mean initial height for each meso_num in the Ericoid table
mean_initial_height <- Ericoid %>%
  group_by(meso_num) %>%
  summarise(mean_initial_height = mean(height_i_cm, na.rm = TRUE))

#merge subsets
merged_means_str_arc_inoc_mgr_height_i <- merge(mean_str_arc, mean_inoc_mgr_pot, by = "plant_sp.x")
merged_means_str_arc_inoc_mgr_height_i <- merge(merged_means_str_arc_inoc_mgr_height_i, mean_initial_height, by = "meso_num")

merged_means_str_arc_inoc_mgr_height_i$mean_BIT <- scale(merged_means_str_arc_inoc_mgr_height_i$mean_BIT)

# Shift the values to make them non-negative
merged_means_str_arc_inoc_mgr_height_i$mean_BIT <- merged_means_str_arc_inoc_mgr_height_i$mean_BIT - min(merged_means_str_arc_inoc_mgr_height_i$mean_BIT)

# Ensure mean_BIT is numeric
merged_means_str_arc_inoc_mgr_height_i$mean_BIT <- as.numeric(as.vector(merged_means_str_arc_inoc_mgr_height_i$mean_BIT))



##MODEL FITTING BI

#lm1_BIT <- lm(mean_MGR ~ mean_BIT * f_comb + mean_initial_height, data = merged_means_str_arc_inoc_mgr_height_i)
#summary(lm1_BIT)
#anova(lm1_BIT)

#lm2_BIT<- lm(mean_MGR ~ mean_BIT + mean_initial_height, data = merged_means_str_arc_inoc_mgr_height_i)
#summary(lm2_BIT)
#anova(lm2_BIT)

lm3_BIT<- lm(mean_MGR ~ f_comb + mean_initial_height, data = merged_means_str_arc_inoc_mgr_height_i)
summary(lm3_BIT)
anova(lm3_BIT)

### GROUPED BARS PLOT
ggplot(merged_means_str_arc_inoc_mgr_height_i, 
       aes(x = plant_sp.x, y = mean_MGR, fill = f_comb)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha= 0.5) +
  labs(
    title = "Impact of Fungal Species on Ericaceous Plants: Mycorrhizal Growth Response",
    x = "Plant Species",
    y = "Mycorrhizal Growth Response",
    fill = "Fungal Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 0.5)
  )

###HEATMAP
ggplot(merged_means_str_arc_inoc_mgr_height_i, 
       aes(x = f_comb, y = plant_sp.x, fill = mean_MGR)) +
  geom_tile() +
  labs(
    title = "Heatmap of Fungal Impact on Ericaceous Plants: Mycorrhizal Growth Response",
    x = "Fungal Species",
    y = "Plant Species",
    fill = "MGR"
  ) +
  scale_fill_gradient2(low = "darkred", mid = "white", high = "palegreen3", 
                       midpoint = 0, na.value = "grey80") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    # Customise background colours
    panel.background = element_rect(fill = "grey79", colour = NA), # Panel (plot area)
    plot.background = element_rect(fill = "white", colour = NA),      # Overall background
    panel.grid = element_blank()                                     # Remove gridlines if needed
  )


emms_BIT <- emmeans(lm3_BIT, ~ f_comb)
cld_results_BIT <- cld(emms_BIT,  Letters = LETTERS)
pairs(emms)
# Convert CLD results to a data frame
df_cld_BIT <- as.data.frame(cld_results_BIT)

ggplot(df_cld_BIT, aes(x = f_comb, y = emmean)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_text(aes(label = .group), vjust = -0.5, hjust = -0.2, size = 5) +
  theme_minimal() +
  labs(y = "Estimated MGR", x = "Fungal Species") +
  ggtitle("Fungal Species Estimated Marginal Means with Grouping Letters") +
  theme(plot.title = element_text(hjust = 0.5))

################################################################################
############################# ROOT DIAMETER ####################################
################################################################################
################################################################################
################################################################################


#Data subset of architecture of sterile plants -- mean of diameter of sterile plants grouped by plant species
sterile_architecture <- subset(Ericoid, f_comb == "STR")

mean_str_arc <- sterile_architecture %>%
  group_by(plant_sp.x) %>% 
  summarise(mean_D = mean(AvgDiamMeasured, na.rm = TRUE)) 


mean_str_arc <- na.omit(mean_str_arc)

#Data subset of inoculated plants 
inoculated_data <- Ericoid[!is.na(Ericoid$MGR_abvg), ]
inoculated_data <- inoculated_data[inoculated_data$f_comb != "STR", ]

#Group by meso_num and fungal combination (f_comb) to calculate the mean MGR
mean_inoc_mgr_pot <- inoculated_data %>%
  group_by(meso_num, plant_sp.x, f_comb) %>%
  summarise(mean_MGR = mean(MGR_abvg, na.rm = TRUE), .groups = "drop")

#Calculate mean initial height for each meso_num in the Ericoid table
mean_initial_height <- Ericoid %>%
  group_by(meso_num) %>%
  summarise(mean_initial_height = mean(height_i_cm, na.rm = TRUE))

#merge subsets
merged_means_str_arc_inoc_mgr_height_i <- merge(mean_str_arc, mean_inoc_mgr_pot, by = "plant_sp.x")
merged_means_str_arc_inoc_mgr_height_i <- merge(merged_means_str_arc_inoc_mgr_height_i, mean_initial_height, by = "meso_num")


merged_means_str_arc_inoc_mgr_height_i$mean_D <- scale(merged_means_str_arc_inoc_mgr_height_i$mean_D)

# Shift the values to make them non-negative
merged_means_str_arc_inoc_mgr_height_i$mean_D <- merged_means_str_arc_inoc_mgr_height_i$mean_D - min(merged_means_str_arc_inoc_mgr_height_i$mean_D)

# Ensure mean_D is numeric
merged_means_str_arc_inoc_mgr_height_i$mean_D <- as.numeric(as.vector(merged_means_str_arc_inoc_mgr_height_i$mean_D))

##MODEL FITTING Diameter
#lm1_D <- lm(mean_MGR ~ mean_D * f_comb + mean_initial_height, data = merged_means_str_arc_inoc_mgr_height_i)
#summary(lm1_D)
#anova(lm1_D)

lm2_D <- lm(mean_MGR ~ mean_D + f_comb + mean_initial_height, data = merged_means_str_arc_inoc_mgr_height_i)
summary(lm2_D)
anova(lm2_D)


emms <- emmeans(lm2_D, ~ f_comb)
cld_results <- cld(emms,  Letters = LETTERS)
pairs(emms)
# Convert CLD results to a data frame
df_cld <- as.data.frame(cld_results)

# Plot with ggplot2
ggplot(df_cld, aes(x = f_comb, y = emmean)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_text(aes(label = .group), vjust = -0.5, hjust = -0.2, size = 5) +
  theme_minimal() +
  labs(y = "Estimated MGR", x = "Fungal Species") +
  ggtitle("Fungal Species Estimated Marginal Means with Grouping Letters") +
  theme(plot.title = element_text(hjust = 0.5))


#similar visualization, but dots are coloured by fungal type
dev.off() 
ggplot(merged_means_str_arc_inoc_mgr_height_i, aes(x = mean_D, y = mean_MGR, color = f_comb)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Linear Regression: MGR ~ Diameter + Fungus + Initial Height",
       x = "Average Root Diameter",
       y = "Mycorrhizal Growth Response (MGR)",
       color = "Fungal Species") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)  # Center the title
  ) +
  # Add plant species names below the plot
  annotate("text", 
           x = unique(merged_means_str_arc_inoc_mgr_height_i$mean_D), 
           y = min(merged_means_str_arc_inoc_mgr_height_i$mean_MGR) - 0.55,  # Adjust to position below plot
           label = c(
             "VAN",
             "GSH",
             "GPR",
             "CVU",
             "VMY",
             "VCO",
             "VVI",
             "PJA",
             "KLA",
             "VMA",
             "GMI"
           ),
           angle = 90, 
           hjust = 0, 
           size = 3)


#models without fungal diversity
lm4_D <- lm(mean_MGR ~ mean_D + mean_initial_height, data = merged_means_str_arc_inoc_mgr_height_i)

# Create a grid of values for prediction
prediction_data <- expand.grid(
  mean_D = seq(min(merged_means_str_arc_inoc_mgr_height_i$mean_D),
               max(merged_means_str_arc_inoc_mgr_height_i$mean_D),
               length.out = 100),
  mean_initial_height = mean(merged_means_str_arc_inoc_mgr_height_i$mean_initial_height)
)


# Predict using the model
prediction_data$predicted_MGR <- predict(lm4_D, newdata = prediction_data)


# Plot the observed data and regression line

ggplot(data = merged_means_str_arc_inoc_mgr_height_i, aes(x = mean_D, y = mean_MGR, color = f_comb)) +
  geom_point(size = 3, alpha = 0.7) +  # Observed data points
  geom_line(data = prediction_data, aes(x = mean_D, y = predicted_MGR), color = "grey45", linewidth = 0.5) +  # Regression line
  labs(
    title = "Linear Regression: MGR ~ D + Initial Height",
    x = "Average Diameter",
    y = "Mycorrhizal Growth Response",
    color = "Fungal Species"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5) ) +
  annotate("text", 
            x = unique(merged_means_str_arc_inoc_mgr_height_i$mean_D), 
            y = min(merged_means_str_arc_inoc_mgr_height_i$mean_MGR) - 0.5,  # Adjust to position below plot
            label = c(
                      "VAN",
                      "GSH",
                      "GPR",
                      "CVU",
                      "VMY",
                      "VCO",
                      "VVI",
                      "PJA",
                      "KLA",
                      "VMA",
                      "GMI"
           ),
           angle = 90, 
           hjust = 0, 
           size = 3)
