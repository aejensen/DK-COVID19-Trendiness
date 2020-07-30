rm(list=ls())

load("DK-COVID19.RData")

library(MASS)

f_post <- total$post[,,1]
df_post <- total$post[,,3]

t0 <- as.Date("2020-03-01")
t_start <- as.numeric(as.Date("2020-07-01") - as.Date("2020-03-01"))
t_end <- length(tPred)

for(i in t_start:t_end) {
  cat(i, "\n")
  
  png(filename = paste("../animation/fig_", sprintf("%03d", i - t_start), ".png", sep=""), width = 1024, height = 1024)
  
  e <- kde2d(f_post[,i], df_post[,i], n = 300, lims=c(0, 10, -2, 2))

  par(bty="n")
  image(e, xlim=c(0, 10), ylim=c(-2, 2),
        xlab="Antal daglige indlæggelser, f(t)", ylab="Hældning for antal daglige indlæggelser, df(t)",
        useRaster=TRUE, col=fields::tim.colors(1024))
  abline(h=0, lty=2)
  title(paste("Dato:", as.Date("2020-03-01") + i))
  
  dev.off()
}
