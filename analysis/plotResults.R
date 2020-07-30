rm(list=ls())

load("DK-COVID19.RData")

band <- function(t, l, u, col) {
  polygon(c(t, rev(t)), c(l, rev(u)), col=col, border = NA)
}

########################################################################
# Figure 1
########################################################################
pdf("fig1.pdf", width = 10, height = 5)

par(mfrow=c(1,2), bty="n", mar = c(2.3, 2.3, 1, 0), mgp=c(1.3,0.4,0))

plot(total$t, total$y, pch = 19, xlab="Antal dage siden 1. marts 2020", 
     ylab="Antal", type="n", ylim=c(0, 100), xaxt="n", xlim=c(0,154))
#axis(1, c(0, 150), label=c("",""))
axis(1, seq(0, 154, 14), cex.axis=0.95)
lines(rep(as.Date("2020-04-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-05-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-06-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-07-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-08-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")

band(tPred, 
     apply(total$post[,,1], 2, quantile, prob = 0.025), 
     apply(total$post[,,1], 2, quantile, prob = 0.975), col = "gray65")
lines(tPred, apply(total$post[,,1], 2, mean), lwd = 2)
points(total$t, total$y, pch = 19, cex=0.5)
legend(82, 98,
       c("Gennemsnit", "95% sandsynlighedsinterval"), 
       col = c("black", "gray65"), lwd = 2, bty="n", cex=0.7, 
       lty = c(1, NA), pch = c(NA, 15), pt.cex=1.5)
title("Antal daglige indlæggelser", font.main=1)

text(as.Date("2020-03-15") - as.Date("2020-03-01"), 96, "Marts", pos=3, cex=0.8)
text(as.Date("2020-04-15") - as.Date("2020-03-01"), 96, "April", pos=3, cex=0.8)
text(as.Date("2020-05-15") - as.Date("2020-03-01"), 96, "Maj", pos=3, cex=0.8)
text(as.Date("2020-06-15") - as.Date("2020-03-01"), 96, "Juni", pos=3, cex=0.8)
text(as.Date("2020-07-15") - as.Date("2020-03-01"), 96, "Juli", pos=3, cex=0.8)

#

plot(tPred, apply(total$post[,,3], 2, mean), lwd = 2, type="n", yaxt="n", xlim=c(0,154),
     ylim=c(-6, 8), xlab="Antal dage siden 1. marts 2020", 
     ylab="Hældning", xaxt="n")

lines(rep(as.Date("2020-04-01") - as.Date("2020-03-01"), 2), c(-6, 8), lty=3, col="gray50")
lines(rep(as.Date("2020-05-01") - as.Date("2020-03-01"), 2), c(-6, 8), lty=3, col="gray50")
lines(rep(as.Date("2020-06-01") - as.Date("2020-03-01"), 2), c(-6, 8), lty=3, col="gray50")
lines(rep(as.Date("2020-07-01") - as.Date("2020-03-01"), 2), c(-6, 8), lty=3, col="gray50")

axis(2, seq(-8, 8, 2))
#axis(1, c(0, 150), label=c("",""))
axis(1, seq(0, 154, 14), cex.axis=0.95)
band(tPred, 
     apply(total$post[,,3], 2, quantile, prob = 0.025), 
     apply(total$post[,,3], 2, quantile, prob = 0.975), col = "gray65")
lines(tPred, apply(total$post[,,3], 2, mean), lwd = 2)
abline(h = 0, lty = 2)
title("Trend for antal indlæggelser", font.main=1)

text(as.Date("2020-03-15") - as.Date("2020-03-01"), 7.5, "Marts", pos=3, cex=0.8)
text(as.Date("2020-04-15") - as.Date("2020-03-01"), 7.5, "April", pos=3, cex=0.8)
text(as.Date("2020-05-15") - as.Date("2020-03-01"), 7.5, "Maj", pos=3, cex=0.8)
text(as.Date("2020-06-15") - as.Date("2020-03-01"), 7.5, "Juni", pos=3, cex=0.8)
text(as.Date("2020-07-15") - as.Date("2020-03-01"), 7.5, "Juli", pos=3, cex=0.8)

dev.off()

########################################################################
# Figure 2
########################################################################
pdf("fig2.pdf", width = 10, height = 5)

par(mfrow=c(1,1), bty="n", mar = c(2.3, 2.3, 1, 0), mgp=c(1.4,0.4,0))
plot(tPred, t(total$post[1,,5])*100, type="n", lty = 1, lwd = 2,
     xlab="Antal dage siden 1. marts 2020", 
     ylab="Sandsynlighed for voksende antal indlæggelser [%]", 
     ylim=c(0,100), xaxt="n", xlim=c(0, 154))

lines(rep(as.Date("2020-03-11") - as.Date("2020-03-01"), 2), c(0, 80), lty=1)
points(as.Date("2020-03-11") - as.Date("2020-03-01"), 80, pch=21, bg="white")
lines(rep(as.Date("2020-04-15") - as.Date("2020-03-01"), 2), c(0, 80), lty=1)
points(as.Date("2020-04-15") - as.Date("2020-03-01"), 80, pch=21, bg="white")
lines(rep(as.Date("2020-04-20") - as.Date("2020-03-01"), 2), c(0, 80), lty=1)
points(as.Date("2020-04-20") - as.Date("2020-03-01"), 80, pch=21, bg="white")
lines(rep(as.Date("2020-05-18") - as.Date("2020-03-01"), 2), c(0, 80), lty=1)
points(as.Date("2020-05-18") - as.Date("2020-03-01"), 80, pch=21, bg="white")
lines(rep(as.Date("2020-05-27") - as.Date("2020-03-01"), 2), c(0, 80), lty=1)
points(as.Date("2020-05-27") - as.Date("2020-03-01"), 80, pch=21, bg="white")

text(as.Date("2020-03-11") - as.Date("2020-03-01"), 80, "Danmark\nlukker", pos = 3, cex=0.7)
text(as.Date("2020-04-15") - as.Date("2020-03-01"), 80, "Skoler op til\n5. klasse åbner", pos = 3, cex=0.7)
text(as.Date("2020-04-20") - as.Date("2020-03-01"), 87, "Frisører og køre-\nskoler åbner", pos = 3, cex=0.7)
text(as.Date("2020-05-18") - as.Date("2020-03-01"), 80, "Restauranter\nåbner", pos = 3, cex=0.7)
text(as.Date("2020-05-27") - as.Date("2020-03-01"), 83, "Gymnasier\nåbner", pos = 3, cex=0.7)

lines(rep(as.Date("2020-04-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-05-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-06-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-07-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")
lines(rep(as.Date("2020-08-01") - as.Date("2020-03-01"), 2), c(0, 100), lty=3, col="gray50")

lines(tPred, t(total$post[1,,5])*100, type="l", lty = 1, lwd = 2)
#axis(1, c(0, 150), label=c("",""))
axis(1, seq(0, 154, 14), cex.axis=0.95)
abline(h = 50, lty = 2)
title("Trend Direction Index", font.main=1)

text(as.Date("2020-03-15") - as.Date("2020-03-01"), 99, "Marts", pos=3, cex=0.8)
text(as.Date("2020-04-15") - as.Date("2020-03-01"), 99, "April", pos=3, cex=0.8)
text(as.Date("2020-05-15") - as.Date("2020-03-01"), 99, "Maj", pos=3, cex=0.8)
text(as.Date("2020-06-15") - as.Date("2020-03-01"), 99, "Juni", pos=3, cex=0.8)
text(as.Date("2020-07-15") - as.Date("2020-03-01"), 99, "Juli", pos=3, cex=0.8)

dev.off()

########################################################################
# Summary statistics
########################################################################
round(total$post[1,,5]*100, 2)

library(rootSolve)
as.Date("2020-03-01") + round(uniroot.all(function(x) approxfun(tPred, total$post[1,,5])(x) - 0.5, c(0, 140)))

#uSeq <- seq(0, 2, length.out=100)
#alarmMat <- sapply(uSeq, function(u) apply(total$post[,,3], 2, function(q) mean(q > u))*100)
#fields::image.plot(tPred[100:200], uSeq, alarmMat[100:200,], col=fields::tim.colors(128))
#fields::image.plot(tPred, uSeq, alarmMat, col=fields::tim.colors(128))
