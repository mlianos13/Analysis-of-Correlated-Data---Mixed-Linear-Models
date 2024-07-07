# Libraries
library(nlme) #!!!!****
library(ggplot2)
library(emmeans)
library(plyr)
library(dplyr)
library(MASS)
library(car)
library(multcomp)
library(lmerTest)
library(Matrix)
library(lme4)

# Reading in data:
# In a study on the protein percentage in milk, cows were given three types of diets. 37 cows
# were randomly assigned to a diet, and measurements of the protein percentage were taken
# for each cow for 19 consecutive weeks (the first 3 weeks were a settling-in period and are not
# included in the data). The aim of the study was to investigate how the three diets affect the
# milk production over time in terms of the protein percentage.
# The variable [Cow] identifies the cow, [Diet] denotes the assigned diet, [Time] denotes the week
# (4-19) the measurement was taken, and [Protein] denotes the measured protein percentage.
# The data is available in the file assignment5.txt.

# we will treat Time as continous variable

assignment5 <- read.table("assignment5.txt", header = TRUE, sep = "\t")
assignment5$Cow <- as.factor(assignment5$Cow)
assignment5$Diet <- as.factor(assignment5$Diet)
assignment5$Time <- as.numeric(assignment5$Time)
assignment5$Timef <- NULL
summary(assignment5)

# Check for missing values
any(is.na(assignment5))

#Check if Factors are Balanced
assignment5 %>% 
  group_by(Cow, Diet,Time) %>%  # Use actual factor names here
  summarise(count = n(), .groups = "drop") %>% 
  distinct(count) %>% 
  nrow() -> num_unique_counts
if (num_unique_counts == 1) {
  print("The dataset is balanced.")
} else {
  print("The dataset is unbalanced.")
}

#---------------Diagrams--------------# 
ggplot(assignment5, aes(x=Time, y=protein, group=Cow, colour=Diet)) + geom_line()


# Calculate the mean protein level for each combination of Diet and Time
mns <- assignment5 %>% 
  group_by(Diet, Time) %>% 
  summarise(protein = mean(protein, na.rm = TRUE), .groups = "drop")
ggplot(mns, aes(x = Time, y = protein, group = Diet, colour = Diet)) +
  geom_point() +
  geom_line() +
  labs(x = "Time", y = "Protein Percentage") +
  theme_minimal() +
  scale_colour_discrete(name = "Diet")



# Plotting the protein percentage over time for each cow under different diets.
with(assignment5, plot(Time, protein, xlab='Week', ylab='Protein Percentage', pch=20))
unique_diets <- unique(assignment5$Diet)
for(d in unique_diets){
    temp <- assignment5[assignment5$Diet == d,]
    with(temp, lines(Time, protein, type="b", col=which(unique_diets == d), lwd=2, pch=19))
}
legend("topleft", legend=unique_diets, lty=1, col=1:length(unique_diets))


# Histogram for Protein
f <- function(x) {
  dnorm(x, mean = mean(assignment5$protein), sd = sd(assignment5$protein))
}
hist(assignment5$protein, xlab='Protein', probability=T, main = "Distribution of Protein")
curve(f, from = min(assignment5$protein), to = max(assignment5$protein), lwd=3, col="red", add=T)



#----------------------Gaussian model of spatial correlation-----------------------------------#


model <- lme(fixed = protein ~ Diet * Time,  # Interaction between fixed effects
             random = ~ Diet + Time | Cow,  # Random slopes for fixed effects within 'Cow'
             correlation = corGaus(form = ~ Time | Cow), 
             data = assignment5)
summary(model)

#----------------------Model Reduction------------------------#

model2 <- update(model, .~. - (Diet + Time | Cow))
anova(model, model2)
0.5*(1-pchisq(2*(logLik(model)-logLik(model2)),1))

model <- lme(fixed = protein ~ Diet * Time,
                        random = ~ 1 | Cow,  
                        correlation = corGaus(form = ~ Time | Cow),
                        data = assignment5)
summary(model)

model3 <- update(model, correlation = NULL)
anova(model, model3)
0.5*(1-pchisq(2*(logLik(model)-logLik(model3)),1))


model4 <- update(model, .~. - (Diet:Time))
model_ml <- update(model, method = "ML")
model4_ml <- update(model4, method = "ML")
anova(model_ml, model4_ml)


plot(Variogram(model, form= ~ Time | Diet, data=assignment5))
plot(Variogram(model4, form= ~ Time | Diet, data=assignment5))


# Final model

#Gaussian spatial correlation model
model <- lme(fixed = protein ~ Diet + Time,
                        random = ~ 1 | Cow,  
                        correlation = corGaus(form = ~ Time | Cow),
                        data = assignment5)

plot(Variogram(model, form= ~ Time | Diet, data=assignment5))

#Exponential correlation model with the nugget effect
model2 <- lme(protein ~ Diet + Time, 
                    random= ~ 1 | Cow,
correlation=corExp(form= ~ Time | Cow, nugget = TRUE),
data = assignment5)

anova(model,model2)

plot(Variogram(model2, form = ~as.numeric(Time) | Cow, data = assignment5))

model3 <- lme(protein ~ Diet + Time, 
                    random= ~ 1| Cow,
correlation=corExp(form= ~Time | Cow),
data=assignment5)

anova(model2,model3)

#-------------------Model Diagnostics------------------------#

#Exponential correlation model with the nugget effect
model <- lme(protein ~ Diet + Time, 
                    random= ~ 1 | Cow,
correlation=corExp(form= ~ Time | Cow, nugget = TRUE),
data = assignment5)
plot(Variogram(model, form = ~as.numeric(Time) | Cow, data = assignment5))

# Checking residuals of the mixed-effects model
plot(residuals(model))

# Checking normality of residuals
qqnorm(residuals(model))
qqline(residuals(model))

# Plotting residuals vs. fitted values
plot(fitted(model), residuals(model))

plot(augPred(model), addData = TRUE)

# Checking for influential observations
influence_measures <- influence(model)
plot(influence_measures, id = 0.05)  # Plot with 5% cutoff for influential points

#---------------------------Statistical Analysis--------------------------#

model <- lme(protein ~ Diet + Time, 
                    random= ~ 1 | Cow,
correlation=corExp(form= ~ Time | Cow, nugget = TRUE),
data = assignment5)
summary(model)
anova(model)

intervals(model, which = "var-cov")

#---------------------------Post-hoc---------------------------------#

#LS-means:
emm <- emmeans(model, "Diet", by="Time", data=assignment5)
print(summary(emm))

# Pairwise
pairwise_emm <- pairs(emm, adjust = "tukey")
print(pairwise_emm)

emms <- emmeans(model, specs = ~ Diet, at = list(Time = 11.5))
emms_df <- as.data.frame(emms)
emm_plot <- ggplot(emms_df, aes(x = Diet, y = emmean, ymin = lower.CL, ymax = upper.CL)) +
            geom_pointrange() +
            xlab("Diet") +
            ylab("Estimated Protein Percentage") +
            ggtitle("Confidence Intervals for Protein Percentage by Diet at Time = 11.5 weeks") +
            theme_minimal()

print(emm_plot)


#----------------Cross-Validation----------------#

model <- lmer(protein ~ Diet + Time + (1|Cow), data = assignment5)
summary(model)
anova(model)