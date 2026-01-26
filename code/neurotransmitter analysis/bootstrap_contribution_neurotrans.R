# This script uses bootstrap resampling to assess the relative importance
# of neurotransmitter receptor systems in predicting weighted degree

# Install and load required packages
install.packages("relaimpo")
library("relaimpo")

# Load data
setwd('/data/results/neurotransmitter')
data <- read.table("regression_model.txt",header = F,sep=",")
 

Y = data[,1]  # weighted degree (z-scored)

X = data[,2:9]  # 8 neurotransmitter systems densities (z-scored)

# Regression model
fit_factor <- lm(Y~V2+V3+V4+V5+V6+V7+V8+V9, data = X)

# Bootstrap for assessing relative importance
crf <- calc.relimp(fit_factor,type=c("lmg"))
B <- boot.relimp(fit_factor,b=1000,type=c("lmg"),rank=TRUE,diff=TRUE,rela=TRUE)
booteval.relimp(B)