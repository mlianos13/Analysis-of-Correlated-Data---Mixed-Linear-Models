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

#reading in data:
assignment2 <- read.table("assignment2.txt", header = TRUE, sep = "\t")
str(assignment2)
head(assignment2)
assignment2$Assessor <- as.factor(assignment2$Assessor)
assignment2$TVset <- as.factor(assignment2$TVset)  
assignment2$Picture <- as.factor(assignment2$Picture)
assignment2$Repeat <- as.factor(assignment2$Repeat)

#Check if Factors are Balanced
assignment2 %>% 
  group_by(Assessor, TVset, Picture, Repeat) %>% 
  summarise(count = n()) %>% 
  ungroup() %>% 
  arrange(count)

#---------Diagramms---------#

par(mfrow=c(2,3))
for (var_name in c('Assessor', 'TVset', 'Picture')) {
  k <- length(levels(as.factor(assignment2[,var_name]))) + 1
  boxplot(as.formula(paste("Coloursaturation ~", var_name)), 
          col = 2:k, main = paste("Coloursaturation vs.", var_name), 
          data = assignment2)
}

# Histogram 
f <- function(x) {
  dnorm(x, mean = mean(assignment2$Coloursaturation), sd = sd(assignment2$Coloursaturation))
}
hist(assignment2$Coloursaturation, xlab='Coloursaturation', probability=T, main = "Distribution of Coloursaturation")
curve(f, from = min(assignment2$Coloursaturation), to = max(assignment2$Coloursaturation), lwd=3, col="red", add=T)
rm(f)
par(mfrow=c(1,1))


#Interaction Plot
par(mfrow=c(2,2))
with(assignment2, {
  interaction.plot(TVset, Picture, Coloursaturation, legend=FALSE, 
                   bty="n", col=2:4, xtick = TRUE, 
                   main = "Interaction Between TVset and Picture")
  interaction.plot(TVset, Assessor, Coloursaturation, legend=FALSE, 
                   bty="n", col=2:9, xtick = TRUE,
                   main = "Interaction Between TVset and Assessor")
  interaction.plot(Picture, Assessor, Coloursaturation, legend=FALSE, 
                   bty="n", col=2:9, xtick = TRUE,
                   main = "Interaction Between Picture and Assessor")
})
par(mfrow=c(1,1))

#Check the significance of the interaction

analysis.lmer <- lmer(Coloursaturation ~ TVset + Picture + Repeat + TVset:Picture + TVset:Repeat
                      + Picture:Repeat + (1|Assessor) + (1|Assessor:TVset) + (1|Assessor:Picture) + (1|Assessor:Repeat), data = assignment2)
ranova(analysis.lmer)
anova(analysis.lmer)

#Factor Diagram
y.names <- c(expression("[I]" [114]^{144}),
                    expression("[Assessor:TVset]" [14]^{24}),
                    expression("TVset:Picture" [4]^{9}),
                    expression("[Assessor]" [7]^{8}),
                    expression("TVset" [2]^{3}),
                    expression("Picture" [2]^{3}),
                    expression("[0]" [1]^{1}))
M <- matrix(nrow = 7, ncol = 7, byrow = TRUE, data = 0)
M[2,1] <- M[3,1] <- M[4,2] <- M[5,2] <- M[5,3] <- M[6,3] <- M[7,4] <- M[7,5] <- M[7,6] <- ""
plotmat(M, pos = c(1,2,3,1), name = y.names, lwd = 2, 
        box.lwd = 2, cex.txt = 1, box.size = 0.12,
        box.type ="square", box.prop = 0.2, arr.type = "triangle", curve = 0)

#---------Fixed model---------#
analysis.lm <- lm(Coloursaturation ~ Assessor + TVset + Picture + Assessor:TVset + TVset:Picture, data = assignment2)
summary(analysis.lm)

#---------Mixed model---------#
analysis.lmer <- lmer(Coloursaturation ~ TVset + Picture + TVset:Picture + (1|Assessor) + (1|Assessor:TVset), data = assignment2)
summary(analysis.lmer)

#---------Model Diagnostics---------#

par(mfrow=c(2,2))
plot(analysis.lm, which=1:4)
par(mfrow=c(1,1))
par(mfrow=c(2,2))
plot(as.numeric(assignment2$Assessor), rstandard(analysis.lm),
     xlab = 'Assessor', ylab = 'Standardized residuals')
plot(predict(analysis.lm), rstandard(analysis.lm), 
     xlab = 'Predicted values', ylab = 'Standardized residuals')
plot(as.numeric(assignment2$Picture), rstandard(analysis.lm), 
     xlab = 'Picture', ylab = 'Standardized residuals')
plot(as.numeric(assignment2$TVset), rstandard(analysis.lm), 
     xlab = 'TVset', ylab = 'Standardized residuals')
par(mfrow=c(1,1))
studresid = studres(analysis.lm)
par(mfrow = c(2, 2))
plot(studresid ~ fitted(analysis.lm), 
     xlab = "Fitted Values", 
     ylab = "Studentized Residuals", 
     main = "Studentized Residuals vs Fitted Values")
with(assignment2, plot(studresid ~ Assessor, 
                       xlab = "Assessor", 
                       ylab = "Studentized Residuals", 
                       col = heat.colors(length(unique(Assessor))), 
                       main = "Studentized Residuals vs Assessor"))
with(assignment2, plot(studresid ~ TVset, 
                       xlab = "TVset", 
                       ylab = "Studentized Residuals", 
                       col = rainbow(length(unique(TVset))), 
                       main = "Studentized Residuals vs TVset"))
with(assignment2, plot(studresid ~ Picture, 
                       xlab = "Picture", 
                       ylab = "Studentized Residuals", 
                       col = rainbow(length(unique(Picture))), 
                       main = "Studentized Residuals vs Picture"))
par(mfrow = c(1, 1))

##Check the normality of the random effects
temp<-names(ranef(analysis.lmer))
temp
par(mfrow=c(1,2))
for(i in 1:2){
qqnorm(unlist(ranef(analysis.lmer)[[i]]),main=paste(temp[i]),cex.main=2)
lines((-3):3,sd(unlist(ranef(analysis.lmer)[[i]]))*((-3):3),col="red")
}
par(mfrow=c(1,1))

#Box Cox
boxcox(analysis.lm)
par(mfrow = c(1, 1))

#New model for transformed data
analysis.lmer.tr <- lmer((Coloursaturation)^(3/2) ~ TVset + Picture + TVset:Picture + (1|Assessor) + (1|Assessor:TVset), data = assignment2)
summary(analysis.lmer)

analysis.lm.tr <- lm((Coloursaturation)^(3/2) ~ Assessor + TVset + Picture + Assessor:TVset + TVset:Picture, data = assignment2)
boxcox(analysis.lm.tr)

#---------Model Diagnostics for transformed data---------#
par(mfrow=c(2,2))
plot(analysis.lm.tr, which=1:4)
par(mfrow=c(1,1))
par(mfrow=c(2,2))
plot(as.numeric(assignment2$Assessor), rstandard(analysis.lm.tr),
     xlab = 'Assessor', ylab = 'Standardized residuals')
plot(predict(analysis.lm.tr), rstandard(analysis.lm.tr), 
     xlab = 'Predicted values', ylab = 'Standardized residuals')
plot(as.numeric(assignment2$Picture), rstandard(analysis.lm.tr), 
     xlab = 'Picture', ylab = 'Standardized residuals')
plot(as.numeric(assignment2$TVset), rstandard(analysis.lm.tr), 
     xlab = 'TVset', ylab = 'Standardized residuals')
par(mfrow=c(1,1))
studresid = studres(analysis.lm.tr)
par(mfrow = c(2, 2))
plot(studresid ~ fitted(analysis.lm.tr), 
     xlab = "Fitted Values", 
     ylab = "Studentized Residuals", 
     main = "Studentized Residuals vs Fitted Values")
with(assignment2, plot(studresid ~ Assessor, 
                       xlab = "Assessor", 
                       ylab = "Studentized Residuals", 
                       col = heat.colors(length(unique(Assessor))), 
                       main = "Studentized Residuals vs Assessor"))
with(assignment2, plot(studresid ~ TVset, 
                       xlab = "TVset", 
                       ylab = "Studentized Residuals", 
                       col = rainbow(length(unique(TVset))), 
                       main = "Studentized Residuals vs TVset"))
with(assignment2, plot(studresid ~ Picture, 
                       xlab = "Picture", 
                       ylab = "Studentized Residuals", 
                       col = rainbow(length(unique(Picture))), 
                       main = "Studentized Residuals vs Picture"))
par(mfrow = c(1, 1))

##Check the normality of the random effects
temp<-names(ranef(analysis.lmer.tr))
temp
par(mfrow=c(1,2))
for(i in 1:2){
qqnorm(unlist(ranef(analysis.lmer.tr)[[i]]),main=paste(temp[i]),cex.main=2)
lines((-3):3,sd(unlist(ranef(analysis.lmer.tr)[[i]]))*((-3):3),col="red")
}
par(mfrow=c(1,1))

#---------Statistical analysis---------# 

model <- lmer((Coloursaturation)^(3/2) ~ TVset + Picture + TVset:Picture + (1|Assessor) + (1|Assessor:TVset), data = assignment2, REML = TRUE)
anova(model, test.statistic = "F", type = 3)
ranova(model)
anova(analysis.lmer.tr, analysis.lm.tr)

#Post-hoc

# Estimation of variance parameters
VarCorr(model)

# Profile likelihood-based confidence intervals for the variance parameters
profile_ci <- profile(model, which=1:2, signames=FALSE)
confint(profile_ci)

# Estimated mean levels and their differences for TVset
emm_tvset <- emmeans(model, pairwise ~ TVset)
emm_tvset
plot(emm_tvset[[2]], xlab = 'Estimate')

# Estimated mean levels and their differences for Picture
emm_picture <- emmeans(model, pairwise ~ Picture)
emm_picture
plot(emm_picture[[2]], xlab = 'Estimate')

# Compact letter displays for TVset
tuke_int_tvset <- glht(model, linfct = mcp(TVset="Tukey"))
tuke_int_tvset_2 <- cld(tuke_int_tvset)
tuke_int_tvset_2
plot(tuke_int_tvset_2, col=2:6, main = '')

# Compact letter displays for Picture
tuke_int_picture <- glht(model, linfct = mcp(Picture="Tukey"))
tuke_int_picture_2 <- cld(tuke_int_picture)
tuke_int_picture_2
plot(tuke_int_picture_2, col=2:6, main = '')

# Estimated mean levels and their differences for TVset:Picture
emm_TVPi <- emmeans(model, pairwise ~ TVset:Picture)
emm_TVPi
plot(emm_TVPi[[2]], xlab = 'Estimate')
