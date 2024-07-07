#Libraries
library(ggplot2)
library(diagram)
library(emmeans)
library(lmerTest)
library(emmeans)
library(dplyr)
library(MASS)
library(car)
library(multcomp)

# Reading in data:
assignment3 <- read.table("assignment3.txt", header = TRUE, sep = "\t")
str(assignment3)
head(assignment3)

#----------------Convert N to Nc---------------------#
assignment3$Nc <- rep(rep(c(0,60,90,120,150,180), each = 4), 3)
#----------------------------------------------------#

# Converting necessary variables to factors
assignment3$Block <- as.factor(assignment3$Block)
assignment3$N <- as.factor(assignment3$N)  
assignment3$Var <- as.factor(assignment3$Var)

#-----------------Explorative-Data-And-Plots--------------------#

# Check if Factors are Balanced
assignment3 %>% 
  group_by(Block, N, Var) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  arrange(count)

#------Checking the effect of the factors and the interactions

#For N:Categorical
interactionN1.lmer <- lmer(Yield ~ N + Var + N:Var + (1|Block) + (1|Block:N), data = assignment3)
interactionN2.lmer <- lmer(Yield ~ N + Var + N:Var + (1|Block), data = assignment3)
0.5*(1-pchisq(2*(logLik(interactionN1.lmer)-logLik(interactionN2.lmer)),1))
interactionN.lm <- lm(Yield ~ N + Var + N:Var, data = assignment3)
0.5*(1-pchisq(2*(logLik(interactionN2.lmer)-logLik(interactionN.lm)),1))

#For N:Numerical
interactionNc1.lmer <- lmer(Yield ~ Nc + Var + Nc:Var + (1|Block) + (1|Block:Nc), data = assignment3)
interactionNc2.lmer <- lmer(Yield ~ Nc + I(Nc^2) + Var + Nc:Var + (1|Block) + (1|Block:Nc), data = assignment3)
anova(interactionNc1.lmer, interactionNc2.lmer)

interactionNc1.lmer <- lmer(Yield ~ Nc + I(Nc^2) + Var + Nc:Var + (1|Block) + (1|Block:Nc), data = assignment3)
interactionNc2.lmer <- lmer(Yield ~ Nc + I(Nc^2) + Var + Nc:Var + (1|Block), data = assignment3)
0.5*(1-pchisq(2*(logLik(interactionNc1.lmer)-logLik(interactionNc2.lmer)),1))
interactionNc.lm <- lm(Yield ~ Nc + I(Nc^2) + Var + Nc:Var, data = assignment3)
0.5*(1-pchisq(2*(logLik(interactionNc2.lmer)-logLik(interactionNc.lm)),1))

#---------Diagramms-and-Plots--------#

par(mfrow=c(1,1))

# Boxplots for each categorical variable
for (var_name in c('Block', 'N', 'Var')) {
  k <- length(levels(as.factor(assignment3[,var_name]))) + 1
  boxplot(as.formula(paste("Yield ~", var_name)), 
          col = 2:k, main = paste("Yield vs.", var_name), 
          data = assignment3)
}

# Histogram for Yield
f <- function(x) {
  dnorm(x, mean = mean(assignment3$Yield), sd = sd(assignment3$Yield))
}
hist(assignment3$Yield, xlab='Yield', probability=T, main = "Distribution of Yield")
curve(f, from = min(assignment3$Yield), to = max(assignment3$Yield), lwd=3, col="red", add=T)

par(mfrow=c(1,1))

# Interaction plots for different pairs of variables
with(assignment3, {
  interaction.plot(N, Var, Yield, legend=FALSE, 
                   bty="n", col=2:6, xtick = TRUE, 
                   main = "Interaction Between N and Var")
  interaction.plot(N, Block, Yield, legend=FALSE, 
                   bty="n", col=2:4, xtick = TRUE,
                   main = "Interaction Between N and Block")
  interaction.plot(Var, Block, Yield, legend=FALSE, 
                   bty="n", col=2:5, xtick = TRUE,
                   main = "Interaction Between Var and Block")
})

par(mfrow=c(1,1))

#Factor Diagram
y.names <- c(expression("[I]" [36]^{72}),
                    expression("N:Var" [10]^{24}),
                    expression("[Block]" [12]^{18}),
                    expression("Var" [3]^{4}),
                    expression("N" [5]^{6}),
                    expression("[0]" [1]^{1}))
M <- matrix(nrow = 6, ncol = 6, byrow = TRUE, data = 0)
M[2,1] <- M[3,1] <- M[5,3] <- M[5,2] <- M[4,2] <- M[6,5] <- M[6,4] <- ""
plotmat(M, pos = c(1,2,2,1), name = y.names, lwd = 2, 
        box.lwd = 2, cex.txt = 1, box.size = 0.12,
        box.type ="square", box.prop = 0.2, arr.type = "triangle", curve = 0)

#----------Models For N:categorical-------------------#

analysis.categ.lm <- lm(Yield ~ N + Var + N:Var, data = assignment3)
anova(analysis.categ.lm)

#----------Models For N:numerical--------------------#

#Fixed Model#
analysis.num.lm <- lm(Yield ~ Nc + Var + Nc:Var, data = assignment3)
anova(analysis.num.lm)


#---------Statistical analysis---------# 

#------------------------------------------N:Categorical---------------------------------------------------#

model1 <- lm(Yield ~ N + Var + N:Var + Block, data = assignment3)

#-----------Model Diagnostic-------------#

par(mfrow=c(2,2))
plot(model1, which=1:4)
par(mfrow=c(1,1))
par(mfrow=c(1,1))
stdresid = rstandard(model1)
par(mfrow=c(1,1))
plot(stdresid ~ fitted(model1), 
     xlab = "Fitted Values", 
     ylab = "Standardized Residuals", 
     main = "Standardized Residuals vs Fitted Values")
with(assignment3, plot(stdresid ~ N, 
                       xlab = "N (Nitrogen Level)", 
                       ylab = "Standardized Residuals", 
                       main = "Standardized Residuals vs N"))
with(assignment3, plot(stdresid ~ Var, 
                       xlab = "Var (Variety)", 
                       ylab = "Standardized Residuals", 
                       main = "Standardized Residuals vs Var"))
par(mfrow=c(1,1))


# Box Cox
boxcox(model1)
par(mfrow = c(1, 1))

#---------Results---------# 
model1 <- lmer(Yield ~ N + Var + N:Var + (1|Block), data = assignment3)
summary(model1)
ranova(model1)
anova(model1)

#-------Post-hoc---------#

# Estimated mean levels and their differences for N
emm_N <- emmeans(model1, pairwise ~ N)
print(emm_N)
tukey_N <- cld(emm_N$emmeans, adjust = "tukey")
print(tukey_N)
plot(emm_N[[2]], xlab = 'Estimate', main = "Pairwise N Comparison")

# Estimated mean levels and their differences for Var
emm_Var <- emmeans(model1, pairwise ~ Var)
print(emm_Var)
tukey_Var <- cld(emm_Var$emmeans, adjust = "tukey")
print(tukey_Var)
plot(emm_Var[[2]], xlab = 'Estimate', main = "Pairwise Var Comparison")


#------------------------------------------N:Numerical---------------------------------------------------#

model2 <- lm(Yield ~ Nc + I(Nc^2) + Var + Nc:Var + Block, data = assignment3)

#-----------Model Diagnostic-------------#
#1st Model Diagnostics
par(mfrow=c(1,1))
plot(model2, which=1:4)
par(mfrow=c(1,1))
plot(as.numeric(assignment3$Nc), rstandard(model2),
     xlab = 'Nc (Nitrogen)', ylab = 'Standardized residuals',
     main = 'Residuals vs. Nitrogen Level')
plot(predict(model2), rstandard(model2), 
     xlab = 'Predicted values', ylab = 'Standardized residuals',
     main = 'Residuals vs. Predicted Values')
plot(as.numeric(assignment3$Var), rstandard(model2), 
     xlab = 'Var (Variety)', ylab = 'Standardized residuals',
     main = 'Residuals vs. Variety')
par(mfrow=c(1,1))
studresid = studres(model2)
par(mfrow=c(1,1))
plot(studresid ~ fitted(model2), 
     xlab = "Fitted Values", 
     ylab = "Studentized Residuals", 
     main = "Studentized Residuals vs Fitted Values")
with(assignment3, plot(studresid ~ Nc, 
                       xlab = "Nc (Nitrogen)", 
                       ylab = "Studentized Residuals", 
                       main = "Studentized Residuals vs Nc"))
with(assignment3, plot(studresid ~ Var, 
                       xlab = "Var (Variety)", 
                       ylab = "Studentized Residuals", 
                       main = "Studentized Residuals vs Var"))
par(mfrow=c(1,1))

# Box Cox
boxcox(model2)
par(mfrow = c(1, 1))

#---------Results---------# 

model2 <- lmer(Yield ~ Nc + I(Nc^2) + Var + Nc:Var + (1|Block), data = assignment3)
summary(model2)
ranova(model2)
anova(model2)

#-------Post-hoc---------#

emm_Var <- emmeans(model2, pairwise ~ Var)
print(emm_Var)
plot(emm_Var[[2]], xlab = 'Estimate')

# Prediction across a range of 'Nc' values for each 'Var'
vars <- levels(assignment3$Var)
nc_vals <- seq(from = min(assignment3$Nc), to = max(assignment3$Nc), length.out = 100)
new_data <- expand.grid(Nc = nc_vals, Var = vars, Block = factor(levels(assignment3$Block)))
# Predict Yield
new_data$Yield_pred <- predict(model2, newdata = new_data, re.form = NA)
new_data$Yield_pred <- new_data$Yield_pred
# Plot the predictions
ggplot(new_data, aes(x = Nc, y = Yield_pred, color = Var)) +
  geom_line(size = 1.5) +
  labs(x = "Nitrogen concentration (Nc)", y = "Predicted Yield", color = "Variety") +
  theme_minimal()
