# install packages 
cran = c("lattice", "sn", "kableExtra")
github = c("chindhanai/skewtInfo")
newCran <- cran[!(cran %in% installed.packages()[,"Package"])]
if(length(newCran)) {
  install.packages(newCran)
}
devtools::install_github(github)

# load packages
library(sn)
library(kableExtra)
library(skewtInfo)
library(lattice)

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
# But the estimates none-the-less seem good, and go into Table 8 below
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

