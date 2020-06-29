load("DK-COVID19.RData")

band <- function(t, l, u, col) {
  polygon(c(t, rev(t)), c(l, rev(u)), col=col, border = NA)
}

makePlot <- function(m) {
  par(mfrow=c(1,3), bty="n", mar = c(2.3, 2.3, 1, 0), mgp=c(1.3,0.4,0))
  plot(m$t, m$y, pch = 19, xlab="Antal dage siden 1. Marts", ylab="f | Y", 
       type="n", ylim=c(0, 100), xaxt="n", xlim=c(0,120))
  axis(1, seq(0, 120, length.out=11))
  #band(tPred, apply(m$post[,,2], 2, quantile, prob = 0.025), 
  #     apply(m$post[,,2], 2, quantile, prob = 0.975), col = "gray85")
  band(tPred, apply(m$post[,,1], 2, quantile, prob = 0.025), 
       apply(m$post[,,1], 2, quantile, prob = 0.975), col = "gray65")
  lines(tPred, apply(m$post[,,1], 2, mean), lwd = 2)
  points(m$t, m$y, pch = 19, cex=0.5)
  #legend("topleft", 
  #       c("Mean", "95% credible interval", "95% posterior prediction interval"), 
  #       col = c("black", "gray65", "gray85"), lwd = 2, bty="n", cex=0.7, 
  #       lty = c(1, NA, NA), pch = c(NA, 15, 15), pt.cex=1.5)
  legend("topleft", 
         c("Gennemsnit", "95% sandsynlighedsinterval"), 
         col = c("black", "gray65"), lwd = 2, bty="n", cex=0.7, 
         lty = c(1, NA), pch = c(NA, 15), pt.cex=1.5)
  title("Antal daglige indlæggelser", font.main=2)
  
  plot(tPred, apply(m$post[,,3], 2, mean), lwd = 2, type="n", yaxt="n", xlim=c(0,120),
       ylim=c(-6, 8), xlab="Antal dage siden 1. Marts", ylab="df | Y", xaxt="n")
  axis(2, seq(-8, 8, 2))
  axis(1, seq(0, 120, length.out=11))
  band(tPred, apply(m$post[,,3], 2, quantile, prob = 0.025), 
       apply(m$post[,,3], 2, quantile, prob = 0.975), col = "gray65")
  lines(tPred, apply(m$post[,,3], 2, mean), lwd = 2)
  abline(h = 0, lty = 2)
  legend("topleft", 
         c("Gennemsnit", "95% sandsynlighedsinterval"), col = c("black", "gray65"), 
         lwd = 2, bty="n", cex=0.7, lty = c(1, NA), pch = c(NA, 15), pt.cex=1.5)
  title("Ændring i antal daglige indlæggelser", font.main=2)
  
  plot(tPred, t(m$post[1,,5])*100, type="l", lty = 1, lwd = 2, xlab="Antal dage siden 1. Marts", 
       ylab="TDI [%]", ylim=c(0,100), xaxt="n", xlim=c(0, 120))
  axis(1, seq(0, 120, length.out=11))
  abline(h = 50, lty = 2)
  title("Trend Direction Index", font.main=2)
}

pdf("covid19DK_fit.pdf", width = 10, height = 4)
makePlot(total)
dev.off()

# Summary statistics
round(total$post[1,,5]*100, 2)
