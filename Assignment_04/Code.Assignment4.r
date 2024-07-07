# Libraries
library(ggplot2)
library(lmerTest)
library(emmeans)
library(dplyr)
library(MASS)
library(car)
library(multcomp)

# Reading in data:
assignment4 <- read.table("assignment4.txt", header = TRUE, sep = "\t")
str(assignment4)
head(assignment4)

# Converting variables to factors
assignment4$loc <- as.factor(assignment4$loc)
assignment4$year <- as.factor(assignment4$year)  
assignment4$nitro <- as.factor(assignment4$nitro)
assignment4$nitro <- as.numeric(as.character(assignment4$nitro))

# Check if Factors are Balanced
assignment4 %>% 
  group_by(loc, year, nitro) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  arrange(count)


#-------------------------------------------------------------------------------------------------------------------------------------------------------------#


#---------Diagrams and Plots--------#

par(mfrow=c(1,1))

# Boxplots for each categorical variable
for (var_name in c('loc', 'year', 'nitro')) {
  k <- length(levels(as.factor(assignment4[,var_name]))) + 1
  boxplot(as.formula(paste("yield ~", var_name)), 
          col = 2:k, main = paste("Yield vs.", var_name), 
          data = assignment4)
}

# Histogram for Yield
f <- function(x) {
  dnorm(x, mean = mean(assignment4$yield), sd = sd(assignment4$yield))
}
hist(assignment4$yield, xlab='Yield', probability=T, main = "Distribution of Yield")
curve(f, from = min(assignment4$yield), to = max(assignment4$yield), lwd=3, col="red", add=T)

par(mfrow=c(1,1))

#------------------------------------------------------------------------------------------------------------------------------------------------------------#

# Fixed effects and interactions
model1 <- lmer(yield ~ loc + nitro + loc:nitro + I(nitro^2):loc 
# Random effects
+ (1 + I(nitro^2) + nitro | year:loc)
+ (1 + nitro | year:loc) 
+ (1 | year:loc)  
+ (1 | year:nitro)              
+ (1 + nitro | year)           
+ (1 | year),                 
  data = assignment4)                  

print(summary(model1), corr = FALSE)
ranova(model1)
anova(model1)


model2 <- update(model1, ~. -(1+I(nitro^2)+nitro|year:loc))
ranova(model2)
anova(model1,model2)
0.5*(1-pchisq(2*(logLik(model1)-logLik(model2)),1))
#p = 0.001318 ** we won't drop (1+I(nitro^2)+nitro|year:loc), 'log Lik.' 1.747123e-06

model2 <- update(model1, ~. -(I(nitro^2):loc))
ranova(model2)
anova(model1,model2)
0.5*(1-pchisq(2*(logLik(model1)-logLik(model2)),1))
#p = 0.0003031 *** we won't drop (I(nitro^2):loc), 'log Lik.' 2.679985e-06

model3 <- update(model1, ~. -(loc:nitro))
ranova(model3)
anova(model1,model3)
0.5*(1-pchisq(2*(logLik(model1)-logLik(model3)),1))
#p = 0.1445 we drop (loc:nitro), 'log Lik.' 0.002328155

model4 <- update(model3, ~. - (1 + nitro | year:loc))
ranova(model4) 
anova(model3, model4)
0.5*(1-pchisq(2*(logLik(model3)-logLik(model4)),1))
#p = 1 so we drop (1 + nitro | year:loc), 'log Lik.' 0.4469559

model5 <- update(model4, ~. -(1 | year:nitro))
ranova(model5)
anova(model4,model5)
0.5*(1-pchisq(2*(logLik(model4)-logLik(model5)),1))
#p = 1 so we drop (1 | year:nitro), 'log Lik.' 0.5

model6 <- update(model5, ~. -(1 + nitro | year))
ranova(model6)
anova(model5,model6)
0.5*(1-pchisq(2*(logLik(model5)-logLik(model6)),1))
#p = 1 so we drop (1 + nitro | year), 'log Lik.' 0.4411685 

model7 <- update(model6, ~. -(1 | year:loc))
ranova(model7)
anova(model6,model7)
0.5*(1-pchisq(2*(logLik(model6)-logLik(model7)),1))
#p = 1 so we drop (1 | year:loc), 'log Lik.' 0.5

model8 <- update(model7, ~. -(1 | year))
ranova(model8)
anova(model7,model8)
0.5*(1-pchisq(2*(logLik(model7)-logLik(model8)),1))
#p = 1 so we drop (1 | year), 'log Lik.' 0.4998249


#-----------------------------------------------------------------------------------------------------------------------------------------------------------#
#Model Diagnostics #1

model1 <- lmer(yield ~ loc + nitro + I(nitro^2):loc + (1 + I(nitro^2) + nitro | year:loc), data = assignment4)

model.lm <- lm(yield ~ loc + nitro + I(nitro^2):loc + nitro:year:loc + I(nitro^2):year:loc, data = assignment4)
par(mfrow=c(1,1))
plot(model.lm, which=1:4)
par(mfrow=c(1,1))
par(mfrow=c(1,1))
stdresid = rstandard(model.lm)
with(assignment4, plot(stdresid ~ nitro, 
                       xlab = "Nitrogen Level", 
                       ylab = "Standardized Residuals", 
                       main = "Standardized Residuals vs Nitrogen Levels"))
with(assignment4, plot(stdresid ~ loc, 
                       xlab = "Location", 
                       ylab = "Standardized Residuals", 
                       main = "Standardized Residuals vs Location"))
with(assignment4, plot(stdresid ~ year, 
                       xlab = "Year", 
                       ylab = "Standardized Residuals", 
                       main = "Standardized Residuals vs Year"))
par(mfrow=c(1,1))

boxcox(model.lm)
par(mfrow = c(1, 1))

#Model Diagnostics #2

model2.lm <- lm((yield^2.6) ~ loc + nitro + I(nitro^2):loc + nitro:year:loc + I(nitro^2):year:loc, data = assignment4)
par(mfrow=c(1,1))
plot(model2.lm,which=1:4)
par(mfrow=c(1,1))
boxcox(model2.lm)
par(mfrow = c(1, 1))

model.lmer <- lmer((yield^2.6) ~ loc + nitro + I(nitro^2):loc + (1 + I(nitro^2) + nitro | year:loc), data = assignment4)
par(mfrow=c(1,1))
# QQ plot for the random intercepts and slopes for nitro at each year:loc combination
qqnorm(unlist(ranef(model.lmer)$`year:loc`), main="Random Effects: year:loc")
qqline(unlist(ranef(model.lmer)$`year:loc`), col="red")
#----------------------------------------------------------------------------------------------------------------------------------------------------#
#Results

model <- lmer((yield^2.6) ~ loc + nitro + I(nitro^2):loc + (1 + I(nitro^2) + nitro | year:loc), data = assignment4)
summary(model)
ranova(model)
anova(model)
#----------------------------------------------------------------------------------------------------------------------------------------------#

#Post-hoc

# Conduct post-hoc analysis for the fixed effect 'loc'
emm_loc <- emmeans(model, specs = pairwise ~ loc)
# Summary of estimated marginal means for 'loc'
print(summary(emm_loc))
# Pairwise comparisons of estimated marginal means for 'loc' with Tukey adjustment
pairwise_loc <- pairs(emm_loc, adjust = "tukey")
print(pairwise_loc)

back_transformed_emm_loc <- summary(emm_loc)$emmeans^(1/2.6)
# Print the back-transformed EMMs
print(back_transformed_emm_loc)
# Plot the back-transformed EMMs
plot(emm_loc, xlab = 'Location', ylab = 'Back-transformed Estimated Marginal Means')

# Here, we use the entire range of 'nitro' values present in the original dataset
nitro_vals <- unique(assignment4$nitro)
locs <- levels(assignment4$loc)
years <- levels(assignment4$year)
new_data <- expand.grid(nitro = nitro_vals, loc = locs, year = years)

# Predict the yield on the transformed scale
new_data$Yield_pred_transformed <- predict(model, newdata = new_data, re.form = NA)

# Back-transform the predictions to the original yield scale
new_data$Yield_pred <- new_data$Yield_pred_transformed^(1/2.6)

# Plot the predictions
ggplot(new_data, aes(x = nitro, y = Yield_pred, color = loc)) +
  geom_line() +
  labs(x = "Nitrogen Level", y = "Predicted Yield", color = "Location") +
  theme_minimal()
