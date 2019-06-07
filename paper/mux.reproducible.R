# install packages
######################################################################
# Please note that version 1.5.2.2 of PerformanceAnalytics from github
# is needed for chart.QQPlot to function properly, you can check the
# commit history of "braverock/PerformanceAnalytics" for more details.
######################################################################
cran = c("lattice", "sn", "kableExtra")
github = c("chindhanai/skewtInfo", "braverock/PerformanceAnalytics")
newCran <- cran[!(cran %in% installed.packages()[,"Package"])]
if(length(newCran)) {
  install.packages(newCran)
}
devtools::install_github(github)

# load packages
library(sn)
library(kableExtra)
library(PerformanceAnalytics)
library(skewtInfo)
library(lattice)

# Figure 1
tsDRet = ts(Dreturns, start = c(1991, 12), frequency = 12)
plot(tsDRet*100, main = "", col = 4, ylab = "RETURNS(%)")
abline(h = 0, lty = 3)

# Define returns in each time period
returns1 <- Dreturns[1:61]
returns2 <- Dreturns[62:217]
returns3 <- Dreturns[218:289]

# Symmetric-t fit
tFitD1 <- st.mple(y = returns1, symmetr = TRUE)$dp
tFitD2 <- st.mple(y = returns2, symmetr = TRUE)$dp
tFitD3 <- st.mple(y = returns3, symmetr = TRUE)$dp
tFit <- rbind(tFitD1, tFitD2, tFitD3)

# Skew-t fit
stFitD1 <- st.mple(y = returns1, penalty = "Qpenalty")$dp
stFitD2 <- st.mple(y = returns2, penalty = "Qpenalty")$dp
stFitD3 <- st.mple(y = returns3, penalty = "Qpenalty")$dp
stFit <- rbind(stFitD1, stFitD2, stFitD3)

dt <- data.frame(tFit, stFit)
row.names(dt) <- c("D1 returns", "D2 returns", "D3 returns")
colnames(dt) <- c("location", "scale", "dof", "location", "scale", "slant", "dof")
############################################################################
# NOTE: The t-dist. dof estimate is not only large, but also highly unstable
############################################################################

# Table 1
kable(dt, format = "latex", booktabs = T, digits = 3) %>%
  kable_styling() %>%
  add_header_above(c(" " = 1, "Symmetric-t" = 3, "Skew-t" = 4))

cpD1 = dp2cp(stFitD1,family = "ST")
cpD2 = dp2cp(stFitD2,family = "ST")
cpD3 = dp2cp(stFitD3,family = "ST")
mat <- data.frame(round(cbind(cpD1, cpD2, cpD3), 3))
row.names(mat) <- NULL
mat <- cbind(c("Mean", "Volatility", "Skewness", "Excess Kurtosis"), mat)
colnames(mat) <- NULL

# Table 2
kable(mat, format = "latex", booktabs = T, digits = 3) %>%
  kable_styling() %>%
  add_header_above(c("Moments" = 1, "D1 Returns" = 1, "D2 Returns" = 1, "D3 Returns" = 1))

# Figure 2
par(mfrow = c(1, 2))
chart.QQPlot(returns1, envelope = 0.95, distribution = "st",
             distributionParameter = "xi = tFitD1[1],omega = tFitD1[2],alpha = 0,nu = tFitD1[3]",
             pch = 20, col = c(4,2), main = "Symmetric-t Distribution for D1", xlab = "Symmetric-t Quantiles")
chart.QQPlot(returns1, envelope = 0.95, distribution = "st",
             distributionParameter = "xi = stFitD1[1],omega = stFitD1[2],alpha = stFitD1[3],nu = stFitD1[4]",
             pch = 20, col = c(4,2), main = "Skew-t Distribution for D1", xlab = "Skew-t Quantiles")
par(mfrow = c(1, 1))

# Figure 3
par(mfrow = c(1, 2))
chart.QQPlot(returns2, envelope = 0.95, distribution = "st",
             distributionParameter = "xi = tFitD2[1],omega = tFitD2[2],alpha = 0,nu = tFitD2[3]",
             pch = 20, col = c(4,2), main = "Symmetric-t Distribution for D2", xlab = "Symmetric-t Quantiles")
chart.QQPlot(returns2, envelope = 0.95, distribution = "st",
             distributionParameter = "xi = stFitD2[1],omega = stFitD2[2],alpha = stFitD2[3],nu = stFitD2[4]",
             pch = 20, col = c(4,2), main = "Skew-t Distribution for D2", xlab = "Skew-t Quantiles")
par(mfrow = c(1, 1))

# Figure 4
par(mfrow = c(1, 2))
chart.QQPlot(returns3, envelope = 0.95, distribution = "st",
             distributionParameter = "xi = tFitD3[1],omega = tFitD3[2],alpha = 0,nu = tFitD3[3]",
             pch = 20, col = c(4,2), main = "Symmetric-t Distribution for D3", xlab = "Symmetric-t Quantiles")
chart.QQPlot(returns3, envelope = 0.95, distribution = "st",
             distributionParameter = "xi = stFitD3[1],omega = stFitD3[2],alpha = stFitD3[3],nu = stFitD3[4]",
             pch = 20, col = c(4,2), main = "Skew-t Distribution for D3", xlab = "Skew-t Quantiles")
par(mfrow = c(1, 1))

# Figure 5
par(mfrow = c(1, 2))
dens <- density(returns1, bw="SJ")
plot(dens, ylab="Density", xlab="D1 Returns", main="Symmetric-t Density Fit for D1",
     col=3, lwd=1, las = 1, ylim = c(0, 14))
rug(returns1, col = 6)
curve(dst(x, xi = tFitD1[1], omega = tFitD1[2], alpha = 0, nu = tFitD1[3]), add=TRUE, col=4, lwd=1)
legend("topright",legend=c("Sym-t","KDE"), lwd=2, col=c(4,3), bty="n", cex = 0.7)

dens <- density(returns1, bw="SJ")
plot(dens, ylab="Density", xlab="D1 Returns", main="Skew-t Density Fit for D1",
     col=3, lwd=1, las = 1, ylim = c(0, 14))
rug(returns1, col = 6)
curve(dst(x, dp = stFitD1), add=TRUE, col=4, lwd=1)
legend("topright",legend=c("Skew-t","KDE"), lwd=2, col=c(4,3), bty="n", cex = 0.7)
par(mfrow = c(1, 1))


# Figure 6
par(mfrow = c(1, 2))
dens <- density(returns2, bw="SJ")
plot(dens, ylab="Density", xlab="D2 Returns", main="Symmetric-t Density Fit for D2",
     col=3, lwd=1, las = 1, ylim = c(0, 9))
rug(returns2, col = 6)
curve(dst(x, xi = tFitD2[1], omega = tFitD2[2], alpha = 0, nu = tFitD2[3]), add=TRUE, col=4, lwd=1)
legend("topright",legend=c("Sym-t","KDE"), lwd=2, col=c(4,3), bty="n", cex = 0.7)

dens <- density(returns2, bw="SJ")
plot(dens, ylab="Density", xlab="D2 Returns", main="Skew-t Density Fit for D2",
     col=3, lwd=1, las = 1, ylim = c(0, 9))
rug(returns2, col = 6)
curve(dst(x, dp = stFitD2), add=TRUE, col=4, lwd=1)
legend("topright",legend=c("Skew-t","KDE"), lwd=2, col=c(4,3), bty="n", cex = 0.7)
par(mfrow = c(1, 1))


# Figure 7
par(mfrow = c(1, 2))
dens <- density(returns3, bw="SJ")
plot(dens, ylab="Density", xlab="D3 Returns", main="Symmetric-t Density Fit for D3", col=3, lwd=1, las = 1)
rug(returns3, col = 6)
curve(dst(x, xi = tFitD3[1], omega = tFitD3[2], alpha = 0, nu = tFitD3[3]), add=TRUE, col=4, lwd=1)
legend("topright",legend=c("Sym-t","KDE"), lwd=2, col=c(4,3), bty="n", cex = 0.7)

dens <- density(returns3, bw="SJ")
plot(dens, ylab="Density", xlab="D3 Returns", main="Skew-t Density Fit for D3", col=3, lwd=1, las = 1)
rug(returns3, col = 6)
curve(dst(x, dp = stFitD3), add=TRUE, col=4, lwd=1)
legend("topright",legend=c("Skew-t","KDE"), lwd=2, col=c(4,3), bty="n", cex = 0.7)
par(mfrow = c(1,1))

# Figure 8
snFitD1 <- sn.mple(y = returns1, penalty = "Qpenalty")$cp
snFitD1 <- cp2dp(cp = snFitD1, family = "SN")

par(mfrow = c(1, 2))
chart.QQPlot(returns1, envelope = 0.95, distribution = "st",
             distributionParameter = "xi = stFitD1[1],omega = stFitD1[2],alpha = stFitD1[3],nu = stFitD1[4]",
             pch = 20, col = c(4,2), main = "Skew-t Distribution for D1", xlab = "Skew-t Quantiles")
chart.QQPlot(returns1, envelope = 0.95, distribution = "sn",
             distributionParameter = "xi = snFitD1[1], omega = snFitD1[2], alpha = snFitD1[3]",
             pch = 20, col = c(4,2), main = "Skew-normal Distribution for D1", xlab = "Skew-normal Quantiles")
par(mfrow = c(1, 1))

expInfo_D2 <- stInfoMat(dp = stFitD2, type = "expected")$stInfoMat  # stFitD2 computed earlier
kable(expInfo_D2, align = 'c', format = "latex", digits = 3, booktabs = T) %>%
  kable_styling(position = "center")

set.seed(0)
###########################################################################
# n <- 1e5  This value was used for paper
n = 200  #This value is just for a quick computation
# The following computation can take several minutes, and you might want to
# use a smaller n just to see that the code works
###########################################################################
rand_asymp <- rst(n = n, dp = stFitD2)
obsInfo <- stInfoMat(y = rand_asymp, dp = stFitD2, type = "observed")
expInfo <- stInfoMat(dp = stFitD2, type = "expected")
obsInfoMat <- obsInfo$stInfoMat
expInfoMat <- expInfo$stInfoMat
AD <- obsInfoMat - expInfoMat
RD <- (obsInfoMat - expInfoMat)/expInfoMat
SE <- obsInfo$SEMat / sqrt(n)
ratio <- AD/SE
diagonal_elt <- data.frame("Diagonal" = c("I_11", "I_22", "I_33", "I_44"),
                           "Expected" = diag(expInfoMat), "Observed" = diag(obsInfoMat),
                           "AD"=diag(AD), "RD"= diag(RD), "SE"=diag(SE), "tstat.AD" = diag(ratio))

# Table 3
kable(diagonal_elt, format = "latex", digits = 3, align = 'c', booktabs = T) %>%
  kable_styling(position = "center")

offDiag_elt <- data.frame("Off-diagonal" = c("I_12", "I_13", "I_14", "I_23", "I_24", "I_34"),
                          "Expected" = expInfoMat[upper.tri(expInfoMat)],
                          "Observed" = obsInfoMat[upper.tri(obsInfoMat)],
                          "AD" = AD[upper.tri(AD)], "RD" = RD[upper.tri(RD)],
                          "SE" = SE[upper.tri(SE)], "tstat.AD" = ratio[upper.tri(ratio)])

kable(offDiag_elt, format = "latex", digits = 3, align = 'c', booktabs = T) %>%
  kable_styling(position = "center")

expInfo_D1 <- stInfoMat(dp = stFitD1, type = "expected")$stInfoMat
expInfo_D2 <- stInfoMat(dp = stFitD2, type = "expected")$stInfoMat # WAS MISSING
expInfo_D3 <- stInfoMat(dp = stFitD3, type = "expected")$stInfoMat
ev1 <- eigen(expInfo_D1)$values
ev2 <- eigen(expInfo_D2)$values
ev3 <- eigen(expInfo_D3)$values
ev <- rbind(ev1, ev2, ev3)
cd <- c(max(ev1)/min(ev1), max(ev2)/min(ev2), max(ev3)/min(ev3))
dt <- data.frame(ev, cd)
row.names(dt) <- c("D1 returns", "D2 returns", "D3 returns")
colnames(dt) <- c("lambda1", "lambda2", "lambda3", "lambda4", "condition numbers")

# Table 4
kable(dt, align = 'c', format = "latex", digits = 4, booktabs = T) %>%
  kable_styling(position = "center")

# Compute and display asymptotic correlation matrices
kable(cov2cor(solve(expInfo_D1)), align = 'c', format = "latex", digits = 3, booktabs = T) %>%
  kable_styling(position = "center")
kable(cov2cor(solve(expInfo_D2)), align = 'c', format = "latex", digits = 3, booktabs = T) %>%
  kable_styling(position = "center")
kable(cov2cor(solve(expInfo_D3)), align = 'c', format = "latex", digits = 3, booktabs = T) %>%
  kable_styling(position = "center")

# Sample code for Monte Carlo Studies and Huge Dof detection
set.seed(123)
#####################################################################
# Skew-t parameters corresponding to the second period of the returns
# The code below takes about a minute
dp <- stFitD2
# M <- 1000 # This is the value used for the paper
M <- 100  # This value is used for seeing how the code works
#####################################################################
n <- 50
retReps <- matrix(NA, ncol = 50, nrow = M)
dim(retReps)
skewtFits <- matrix(NA, ncol = 4, nrow = M)
for (m in 1:M) {
  retReps[m, ] <- rst(n = n, dp = dp)
}
for (m in 1:M) {
  skewtFits[m, ] <- st.mple(y = retReps[m, ], penalty = "Qpenalty")$dp
}
skewtFitsHugeDof <- skewtFits[(skewtFits[, 4] > 1000),]
nrow(skewtFitsHugeDof)
skewtFitsHugeDof[1, ]

# Table 5
n <- c(50, 100, 200, 400)
threshold <- c(20, 40, 100, 500, 1000)
dat <- data.frame(row.names = c("Fit DoF > 20",
                                "Fit DoF > 40",
                                "Fit DoF > 100",
                                "Fit DoF > 500",
                                "Fit DoF > 1000"),
                  matrix(nrow = 5, ncol = 4))
colnames(dat) <- n

for(i in 1:length(n)) {
  infile <- paste("Data/param", n[i],"_D2",".csv", sep="")
  data <- read.csv(infile, header=TRUE)
  data <- data[1:1000, -1]
  for (j in 1:5) {
    dataHugeDof <- data[(data[, 4] > threshold[j]),]
    dat[j, i] <- paste(nrow(dataHugeDof), "(", nrow(dataHugeDof)/10, "%)", sep = " ")
  }
}

# Table 5
kable(dat, align = 'c', format = "latex", booktabs = T) %>%
  kable_styling(position = "center")

# Computations for Table 6
n <- c(100, 200, 400)
fD1 <- data.frame(row.names = round(stFitD1, 3), matrix(nrow = 4, ncol = 6))
fD2 <- data.frame(row.names = round(stFitD2, 3), matrix(nrow = 4, ncol = 6))
fD3 <- data.frame(row.names = round(stFitD3, 3), matrix(nrow = 4, ncol = 6))
colnames(fD1) <- colnames(fD2) <- colnames(fD3) <- NULL
for(j in 1:3) {
  for (i in 1:length(n)) {
    infile <- paste("Data/param", n[i], "_D", j,".csv", sep="")
    data <- read.csv(infile, header=TRUE)
    data <- data[, -1]
    data <- na.omit(data)
    data <- data[(1:(min(nrow(data), 1000))), ]
    if (j == 1) {
      fD1[, (2*i - 1)] <- round(apply(data, 2, mean), 3)
      SD <- round(apply(data, 2, sd)  / sqrt(1000), 4)
      fD1[, (2*i)] <- round(100*(fD1[, (2*i - 1)] - stFitD1)/ stFitD1, 1)
      fD1[4, (2*i - 1)] <- format(mean(data[, 4]), scientific = TRUE, digits = 3)
      fD1[4, (2*i)] <- format(round(100*(mean(data[, 4]) - stFitD1[4])/stFitD1[4], 3),
                              scientific = TRUE, digits = 0)
      for (k in 1:3) {
        fD1[k, (2*i - 1)] <- paste(fD1[k, (2*i - 1)], " (", SD[k], ")", sep = "")
      }
    }
    if (j == 2) {
      fD2[, (2*i - 1)] <- round(apply(data, 2, mean), 3)
      SD <- round(apply(data, 2, sd) / sqrt(1000), 4)
      fD2[, (2*i)] <- round(100*(fD2[, (2*i - 1)] - stFitD2)/ stFitD2, 1)
      fD2[4, (2*i - 1)] <- format(mean(data[, 4]), scientific = TRUE, digits = 3)
      fD2[4, (2*i)] <- format(round(100*(mean(data[, 4]) - stFitD2[4])/stFitD2[4], 3),
                              scientific = TRUE, digits = 0)
      for (k in 1:3) {
        fD2[k, (2*i - 1)] <- paste(fD2[k, (2*i - 1)], " (", SD[k], ")", sep = "")
      }
    }
    if (j == 3) {
      fD3[, (2*i - 1)] <- round(apply(data, 2, mean), 4)
      SD <- round(apply(data, 2, sd)/ sqrt(1000), 4)
      fD3[, (2*i)] <- round(100*(fD3[, (2*i - 1)] - stFitD3)/ stFitD3, 1)
      fD3[4,(2*i - 1)] <- format(mean(data[, 4]), scientific = TRUE, digits = 3)
      fD3[4,(2*i)] <- format(round(100*(mean(data[, 4]) - stFitD3[4])/stFitD3[4], 3),
                             scientific = TRUE, digits = 0)
      for (k in 1:3) {
        fD3[k, (2*i - 1)] <- paste(fD3[k, (2*i - 1)], " (", SD[k], ")", sep = "")
      }
    }
  }
}

# Table 6
kable(fD1, align = 'c', format = "latex", booktabs = T, digits = 3) %>%
  kable_styling(position = "center") %>%
  add_header_above(c("Theta" = 1, "MC" = 1, "RB%" = 1, "MC" = 1, "RB%" = 1, "MC" = 1, "RB%" = 1)) %>%
  add_header_above(c("D1" = 1, "n = 100" = 2, "n = 200" = 2, "n = 400" = 2))

kable(fD2, align = 'c', format = "latex", booktabs = T, digits = 3) %>%
  kable_styling(position = "center") %>%
  add_header_above(c("Theta" = 1, "MC" = 1, "RB%" = 1, "MC" = 1, "RB%" = 1, "MC" = 1, "RB%" = 1)) %>%
  add_header_above(c("D2" = 1, "n = 100" = 2, "n = 200" = 2, "n = 400" = 2))

kable(fD3, align = 'c', format = "latex", booktabs = T, digits = 3) %>%
  kable_styling(position = "center") %>%
  add_header_above(c("Theta" = 1, "MC" = 1, "RB%" = 1, "MC" = 1, "RB%" = 1, "MC" = 1, "RB%" = 1)) %>%
  add_header_above(c("D3" = 1, "n = 100" = 2, "n = 200" = 2, "n = 400" = 2))

# Computations for Table 7
  n <- c(50, 100, 200, 400)
  threshold = 20
  D1 <- data.frame(row.names = round(stFitD1, 3), matrix(nrow = 4, ncol = 8))
  D2 <- data.frame(row.names = round(stFitD2, 3), matrix(nrow = 4, ncol = 8))
  D3 <- data.frame(row.names = round(stFitD3, 3), matrix(nrow = 4, ncol = 8))
  colnames(D1) <- colnames(D2) <- colnames(D3) <- NULL
  NewD1 <- NewD2 <- NewD3 <- matrix(nrow = 4, ncol = 8)
  for(j in 1:3) {
    for (i in 1:length(n)) {
      infile <- paste("Data/param", n[i], "_D", j,".csv", sep="")
      data <- read.csv(infile, header=TRUE)
      data <- data[, -1]
      data <- na.omit(data)
      data <- data[data[, 4] <= threshold, ]
    data <- data[(1:(min(nrow(data), 1000))), ]
    if (j == 1) {
      D1[, (2*i - 1)] <- round(apply(data, 2, mean), 3)
      NewD1[, (2*i - 1)] <- D1[, (2*i - 1)] - stFitD1
      D1[, (2*i)] <- round(100*(D1[, (2*i - 1)] - stFitD1)/ stFitD1, 1)
      SD <- round(apply(data, 2, sd)/ sqrt(1000), 6)
      for (k in 1:4) {
        D1[k, (2*i - 1)] <- paste(D1[k, (2*i - 1)], " (", round(SD[k], 5), ")", sep="")
      }
    }
    if (j == 2) {
      D2[, (2*i - 1)] <- round(apply(data, 2, mean), 3)
      NewD2[, (2*i - 1)] <- D2[, (2*i - 1)] - stFitD2
      D2[, (2*i)] <- round(100*(D2[, (2*i - 1)] - stFitD2)/ stFitD2, 1)
      SD <- round(apply(data, 2, sd)/ sqrt(1000), 6)
      for (k in 1:4) {
        D2[k, (2*i - 1)] <- paste(D2[k, (2*i - 1)], " (", round(SD[k], 5), ")", sep="")
      }
    }
    if (j == 3) {
      D3[, (2*i - 1)] <- round(apply(data, 2, mean), 3)
      NewD3[, (2*i - 1)] <- D3[, (2*i - 1)] - stFitD3
      D3[, (2*i)] <- round(100*(D3[, (2*i - 1)] - stFitD3)/ stFitD3, 1)
      SD <- round(apply(data, 2, sd)/ sqrt(1000), 6)
      for (k in 1:4) {
        D3[k, (2*i - 1)] <- paste(D3[k, (2*i - 1)], " (", round(SD[k], 5), ")", sep="")
      }
    }
  }
}

# Table 7
kable(D1[, -c(1,2)], align = 'c', format = "latex", booktabs = T, digits = 4) %>%
  kable_styling(position = "center") %>%
  add_header_above(c("Theta" = 1, "MC" = 1, "RB%" = 1, "MC" = 1, "RB%" = 1, "MC" = 1, "RB%" = 1)) %>%
  add_header_above(c("D1" = 1, "n = 100" = 2, "n = 200" = 2, "n = 400" = 2))

kable(D2[, -c(1,2)], align = 'c', format = "latex", booktabs = T, digits = 4) %>%
  kable_styling(position = "center") %>%
  add_header_above(c("Theta" = 1, "MC" = 1, "RB%" = 1, "MC" = 1, "RB%" = 1, "MC" = 1, "RB%" = 1)) %>%
  add_header_above(c("D2" = 1, "n = 100" = 2, "n = 200" = 2, "n = 400" = 2))

kable(D3[, -c(1,2)], align = 'c', format = "latex", booktabs = T, digits = 4) %>%
  kable_styling(position = "center") %>%
  add_header_above(c("Theta" = 1, "MC" = 1, "RB%" = 1, "MC" = 1, "RB%" = 1, "MC" = 1, "RB%" = 1)) %>%
  add_header_above(c("D3" = 1, "n = 100" = 2, "n = 200" = 2, "n = 400" = 2))


# Computations for Table 8
n <- c(50, 100, 200, 400, 800)
threshold <- 20
SED2 <- data.frame(matrix(nrow = 20, ncol = 5))
colnames(SED2) <- c("n", "Theta", "MC", "OBS", "EXP")
SED2$n <- factor(c(rep(50, 4),rep(100, 4),rep(200, 4),rep(400, 4),  rep(800, 4)),
                 levels = c("50", "100", "200", "400", "800"))
SED2$Theta <- factor(rep(round(stFitD2, 3), 5), levels = round(stFitD2, 3))
SED3 <- SED2
SED3$Theta <- factor(rep(round(stFitD3, 3), 5), levels = round(stFitD3, 3))
I_2 <- diag(solve(stInfoMat(dp = stFitD2, type = "expected")$stInfoMat))
I_3 <- diag(solve(stInfoMat(dp = stFitD3, type = "expected")$stInfoMat))
for(j in 2:3) {
  for (i in 1:length(n)) {
    infile <- paste("Data/5000/param", n[i], "_D", j,".csv", sep="")
    data <- read.csv(infile, header=TRUE)
    data <- data[, -1]
    data <- na.omit(data)
    data <- data[data[, 4] <= threshold, ]
    data <- data[1:4000, ]
    if (j == 2) {
      SEavg <- apply(data, 2, sd)
      SED2[(4*i - 3):(4*i), 3] <- SEavg
      SED2[(4*i - 3):(4*i), 4] <- sqrt(diag(st.infoUv(dp = stFitD2, y = returns2)$asyvar.dp)) * sqrt(length(returns2)/n[i])
      SED2[(4*i - 3):(4*i), 5] <- sqrt(I_2/n[i])
    }
    if (j == 3) {
      SEavg <- apply(data, 2, sd)
      SED3[(4*i - 3):(4*i), 3] <- SEavg
      SED3[(4*i - 3):(4*i), 4] <- sqrt(diag(st.infoUv(dp = stFitD3, y = returns3)$asyvar.dp)) * sqrt(length(returns3)/n[i])
      SED3[(4*i - 3):(4*i), 5] <- sqrt(I_3/n[i])
    }
  }
}
##########################################################################
# The above code gives a lot of output of score(dp) and norma(score)
# and also gives a lot of warning messages "dp does not seem to be at MLE"
# But numbers produced above do go into Table 8 below
##########################################################################

# Table 8
kable(cbind(SED2[, -c(1,2)], SED3[, 3:5]), align = 'c', format = "latex", booktabs = T, digits = 4,
      linesep = c("", "", "", "\\addlinespace")) %>%
  kable_styling(position = "center") %>%
  group_rows("n = 50", 1, 4) %>%
  group_rows("n = 100", 5, 8) %>%
  group_rows("n = 200", 9, 12) %>%
  group_rows("n = 400", 13, 16) %>%
  group_rows("n = 800", 17, 20) %>%
  add_header_above(c("D2" = 3, "D3" = 3))

# Figure 9
#####################################################################
# The code below produces some similar output and warnings as for the
# code used to produce Table 8 above
# Does produce the top and bottom figures in Figure 9
#####################################################################
for (i in 1:length(n)) {
  infile <- paste("Data/5000/param", n[i], "_D", 2,".csv", sep="")
  data <- read.csv(infile, header=TRUE)
  data <- data[, -1]
  data <- na.omit(data)
  data <- data[data[, 4] <= threshold, ]
  data <- data[1:4000, ]
  assign(paste("cleanSd", n[i], sep=""), apply(data, 2, sd))
  assign(paste("obsInfo", n[i], sep=""), sqrt(diag(st.infoUv(dp = stFitD2, y = returns2)$asyvar.dp)) * sqrt(length(returns2)/n[i]))
  assign(paste("expInfo", n[i], sep=""), sqrt(I_2/n[i]))
}

panel <- as.data.frame(rbind(MC50 = cleanSd50,
                             obsInfo50, expInfo50,
                             MC100 = cleanSd100,
                             obsInfo100, expInfo100,
                             MC200 = cleanSd200,
                             obsInfo200, expInfo200,
                             MC400 = cleanSd400,
                             obsInfo400, expInfo400,
                             MC800 = cleanSd800,
                             obsInfo800, expInfo800))
exCol <- rep(c("MC", "obsI", "expI"), 20)
ind <- c(rep("xi", 15), rep("omega", 15), rep("alpha", 15), rep("nu", 15))
panel <- stack(panel)
sizes <- rep(c("50", "50", "50",
               "100", "100", "100",
               "200", "200", "200",
               "400", "400", "400",
               "800", "800", "800"), 4)
panel$ind <- factor(ind, levels = c("xi", "omega", "alpha", "nu"))
panel <- cbind(panel, exCol = factor(exCol, levels = c("MC", "obsI", "expI")),
               sizes = factor(sizes, levels = c("50", "100", "200", "400", "800")))

xyplot(values ~ sizes | ind, data = panel, group = exCol, type = "b", xlab ="",
       ylab = "", layout = c(1, 4), as.table=TRUE, lty = 1:3, col = 1,lwd = 2,
       strip.left = TRUE, strip = FALSE,
       scales=list(x=list(alternating=FALSE),
                   y=list(alternating = c(1,1), tck=c(1,0),
                          relation = "free", rot = 0)),
       key=list(space="right",
                lines=list(col=1, lty=1:3, lwd=2),
                text=list(c("MC","obsINFO","expINFO"))),
       par.settings = list(strip.background=list(col="grey")) )

for (i in 1:length(n)) {
  infile <- paste("Data/5000/param", n[i], "_D", 3,".csv", sep="")
  data <- read.csv(infile, header=TRUE)
  data <- data[, -1]
  data <- na.omit(data)
  data <- data[data[, 4] <= threshold, ]
  data <- data[1:4000, ]
  assign(paste("cleanSd", n[i], sep=""), apply(data, 2, sd))
  assign(paste("obsInfo", n[i], sep=""), sqrt(diag(st.infoUv(dp = stFitD3, y = returns3)$asyvar.dp)) * sqrt(length(returns3)/n[i]))
  assign(paste("expInfo", n[i], sep=""), sqrt(I_3/n[i]))
}

panel <- as.data.frame(rbind(MC50 = cleanSd50,
                             obsInfo50, expInfo50,
                             MC100 = cleanSd100,
                             obsInfo100, expInfo100,
                             MC200 = cleanSd200,
                             obsInfo200, expInfo200,
                             MC400 = cleanSd400,
                             obsInfo400, expInfo400,
                             MC800 = cleanSd800,
                             obsInfo800, expInfo800))
exCol <- rep(c("MC", "obsI", "expI"), 20)
ind <- c(rep("xi", 15), rep("omega", 15), rep("alpha", 15), rep("nu", 15))
panel <- stack(panel)
sizes <- rep(c("50", "50", "50",
               "100", "100", "100",
               "200", "200", "200",
               "400", "400", "400",
               "800", "800", "800"), 4)
panel$ind <- factor(ind, levels = c("xi", "omega", "alpha", "nu"))
panel <- cbind(panel, exCol = factor(exCol, levels = c("MC", "obsI", "expI")),
               sizes = factor(sizes, levels = c("50", "100", "200", "400", "800")))

xyplot(values ~ sizes | ind, data = panel, group = exCol, type = "b", xlab ="", ylab = "",
       layout = c(1, 4), as.table=TRUE, lty = 1:3, col = 1,lwd = 2,
       strip.left = TRUE, strip = FALSE,
       scales=list(x=list(alternating=FALSE),
                   y=list(alternating = c(1,1), tck=c(1,0),
                          relation = "free", rot = 0)),
       key=list(space="right",
                lines=list(col=1, lty=1:3, lwd=2),
                text=list(c("MC","obsINFO","expINFO"))),
       par.settings = list(strip.background=list(col="grey")) )

# Code for Appendix C
#####################################################################################
# Gives warning "dp does not seem to be at MLE"
# But table numbers ok (see curious the exact zero fo the xi-xi entry from last line)
#####################################################################################
set.seed(0)
rand_st <- rst(n = 50, dp = stFitD2)
(obsInfo_skewtInfo  <- stInfoMat(y = rand_st, dp = stFitD2, type = "observed")$stInfoMat)
(obsInfo_sn <- st.infoUv(y = rand_st, dp = stFitD2)$info.dp / length(rand_st))
obsInfo_skewtInfo - obsInfo_sn
