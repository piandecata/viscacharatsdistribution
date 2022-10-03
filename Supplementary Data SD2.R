rm(list=ls())
options(digits=3)
setwd()

#Center and standardize
library(Hmisc)

#Historical scenarios LGM CCMS4
LGMCCSM <- read.csv("datasets/LGMccsm.csv")
View(LGMCCSM)
vc1 <- varclus(as.matrix(LGMCCSM[,c(1:23)]))
plot(vc1)

#Historical scenarios LGM MIROC ESM
LGMMIROC<- read.csv("datasets/LGMmiroc-esm.csv")
View(LGMMIROC)
vc2 <- varclus(as.matrix(LGMMIROC[,c(1:22)]))
plot(vc2)

#Historical scenarios Mid Holocene CCMS4
MHCCSM <- read.csv("datasets/MHccsm.csv")
View(MHCCSM)
vc3 <- varclus(as.matrix(MHCCSM[,c(1:22)]))
plot(vc3)

#Historical scenarios Mid Holocene MIROC ESM
MHMIROC <- read.csv("datasets/MHmiroc-esm.csv")
View(MHMIROC)
vc4<- varclus(as.matrix(MHMIROC[,c(1:23)]))
plot(vc4)

#Current scenarios
current <- read.csv("datasets/Current.csv")
View(current)
vc5 <- varclus(as.matrix(current[,c(17,20,32,35,42:60)]))
plot(vc5)

#Future scenarios CCMS4 RCP 26 Year 2070
RCP26CCSM70 <- read.csv("datasets/70ccsm26.csv")
View(RCP26CCSM70)
vc6 <- varclus(as.matrix(RCP26CCSM70[,c(2:23)]))
plot(vc6)

#Future scenarios CCMS4 RCP 85
RCP85CCSM70 <- read.csv("datasets/70ccsm85.csv")
View(RCP85CCSM70)
vc7 <- varclus(as.matrix(RCP85CCSM70[,c(1:19,23,29)]))
plot(vc7)

#Future scenarios MIROC ESM RCP 26
RCP26MIROC70 <- read.csv("datasets/70miroc-esm26.csv")
View(RCP26MIROC70)
vc8 <- varclus(as.matrix(RCP26MIROC70[,c(2:23)]))
plot(vc8)

#Future scenarios MIROC ESM RCP 85
RCP85MIROC70 <- read.csv("datasets/70miroc-esm85.csv")
View(RCP85MIROC70)
vc9 <- varclus(as.matrix(RCP85MIROC70[,c(1:19,23,29)]))
plot(vc9)