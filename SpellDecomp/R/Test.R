

# read in HLEdecomp objects.
trans  <- readRDS("/home/tim/git/HLEDecomp/HLEDecomp/Data/Transitions/Collapsed/TR_v06_collapsed.rds")
prev  <- readRDS("/home/tim/git/HLEDecomp/HLEDecomp/Data/Results/mspec06/prev/prev_2.rds")
head(trans)
prev$pi1 + prev$pi2

# select for a single population
previ  <- prev[prev$time==2006 & prev$sex == "f" & prev$edu == "all_edu", ]
transi <- trans[trans$time==2006 & trans$sex == "f" & trans$edu == "all_edu", ]

# 
head(transi)
head(previ)
sum(previ$pi2[-1])
# duration distribution of episodes of state 2 starting at age 52

I     <- diag(31)
U     <- rbind(0, cbind(diag(transi$m22[-1]), 0))
# conditional state-specific remain prob
N2    <- solve(I-U)
# need probability of starting in track to weight columns
# start 2
s2    <- previ$pi1 * transi$m12

E2 <- N2 %*% diag(s2 )
E2[upper.tri(E2,FALSE)] <- NA
matplot(E2,type='l')

# deduct contrib from episodes starting at age 50
d50 <- previ$pi2[1] * c(1,cumprod(transi$m22))
plot(rowSums(E2,na.rm=TRUE),type='l')
lines(previ$pi2[-1] - d50[-c(1,32)], col = "red")
lines(previ$pi2[-1], col = "blue", lty = 2)

sum(previ$pi2) - sum(previ$pi2[1] * c(1,cumprod(transi$m22))) 
sum(E2, na.rm = TRUE)
ind <- prev$sex=="f" & prev$edu=="all_edu"
P1T <- tapply(prev$pi1[ind],prev$time[ind],sum)
P2T <- tapply(prev$pi2[ind],prev$time[ind],sum)

pdf("/home/tim/git/SpellDecomp/SpellDecomp/LabTalk/Figures/Expectancies.pdf", height=8,width=4)
barplot(rbind(P1T,P2T)*2,col=gray(c(.8,.5)),las=1)
mtext(side=2,line=3,"Expected years",cex=1.5)
mtext(side=1,line=3,"Year",cex=1.5)
text(c(.6,.6),c(11,29),c("Healthy","Sick"),srt=90,cex=1.5)
dev.off()

# make object for plotting
entry <- c(col(E2)-1) * 2 + 50
dur   <- c(row(E2) - col(E2)) * 2
pid   <- c(E2)

dE2   <- E2 * c(transi$m24[-1],0)
hE2   <- E2 * c(transi$m21[-1],0)
dE2   <- E2 * transi$m24
hE2   <- E2 * transi$m21

a <- seq(50,110,by=2)
cols <- c("#118925","#ad6b08","#112960","#680862","#025107","#5b0404")
pdf("/home/tim/git/SpellDecomp/SpellDecomp/LabTalk/Figures/TransitionsR.pdf")
matplot(a,transi[,c( "m11","m12","m14","m22","m21","m24")],type='l',
		lwd=2,lty=1,las=1,ylab="Probability",xlab="Age",
		col = cols,cex.lab=1.5,ylim=c(0,1))
dev.off()

pdf("/home/tim/git/SpellDecomp/SpellDecomp/LabTalk/Figures/PrevalenceR.pdf")
plot(a,previ$pi1+previ$pi2,type='l',las=1,cex.lab=1.5,xlab="Age",ylab = "Prevalence")
polygon(c(a,rev(a)),c(rep(0,31),rev(previ$pi1)), border="white",col=gray(.8))
polygon(c(a,rev(a)),c(previ$pi1,rev(previ$pi1+previ$pi2)), border="white",col=gray(.5))
lines(a,previ$pi1+previ$pi2,lwd=2)
dev.off()





# terminal vs other states
T2 <- (row(E2) - col(E2)) * 2 * dE2
I2 <- (row(E2) - col(E2)) * 2 * hE2
sum(T2, na.rm = TRUE) + # terminal
sum(I2, na.rm = TRUE)

# 45\%
sum(T2, na.rm = TRUE) / (sum(T2, na.rm = TRUE) + sum(I2, na.rm = TRUE))

matplot(dE2,type='l')
matplot(t(dE2),type='l')
matplot(hE2/dE2,type='l')
E2D   <- data.frame(entry = entry, dur = dur, Age = seq(50,110,by=2), pi = pid)
E2D   <- E2D[E2D$dur >= 0, ]

lines(tapply(E2D$pi,E2D$Age,sum),col = "magenta")
ageOnly  <- c(0,tapply(E2D$pi,E2D$Age,sum))
durOnly  <- tapply(E2D$pi,E2D$dur,sum)
library(reshape2)
durInAge <- acast(E2D, Age~dur,sum,value.var = "pi",fill=NA)

 age <-seq(50,112,by=2)
 durs <- 0:30 * 2
 library(RColorBrewer)
 blues <- colorRampPalette(brewer.pal(9,"Blues"),space= "Lab")
 rain  <- colorRampPalette(brewer.pal(9,"Spectral"),space= "Lab")
 durInAge[durInAge==0] <- NA
 
plot(age, ageOnly, col = "red",lty=1,lwd=2,type='l')
#matplot(age, durInAge, type = 'l', col = rev(blues(31)),add=TRUE)
durInAgeC <- t(apply(durInAge, 1, cumsum))
durInAgeC <- cbind(0,durInAgeC)
matplot(age, durInAgeC, type = 'l', col = rev(blues(31)), add = TRUE)
lines(age, ageOnly, col = "red", lwd = 2)
cols <- rev(blues(31))

pdf("/home/tim/git/SpellDecomp/SpellDecomp/LabTalk/Figures/PrevalenceDisR.pdf")
plot(a,previ$pi2,type='n',las=1,cex.lab=1.5,xlab="Age",ylab = "Years disabled",xlim=c(50,112))
polygon(c(a,rev(a)),c(rep(0,31),rev(previ$pi2)), border="black",col=gray(.5))
lines(a,previ$pi2,lwd=2)
dev.off()

pdf("/home/tim/git/SpellDecomp/SpellDecomp/LabTalk/Figures/DisByDur.pdf")
plot(age, ageOnly, col = "red",lty=1,lwd=2,type='n',xlab = "Age",
		ylab = "Years disabled",las=1,cex.lab=1.5,xlim=c(50,112))
for (i in 2:31){
	polygon(c(age[-1],rev(age[-1])),c(durInAgeC[,i],rev(durInAgeC[,i+1])),col=cols[i],border="white",lwd=.5)
}
polygon(c(age,rev(age[-1])),c(0,durInAgeC[,1],rev(durInAgeC[,2])),col=cols[1],border="white",lwd=.5)
polygon(c(age,rev(age)),c(ageOnly,0,rev(previ$pi2)),col=gray(.5),border="white",lwd=.5)
lines(age, c(previ$pi2,0), col = "black", lwd = 2)
lines(age, ageOnly, col = "black", lwd = 2)
text(c(60,62,64,66),
		diag(durInAgeC[as.character(c(60,62,64,66)),as.character(c(0,2,4,6))]),
		c(0,2,4,6),pos=1,col="white",cex=1.5,font=2)
dev.off()

d <- seq(0,60,by=2)

pdf("/home/tim/git/SpellDecomp/SpellDecomp/LabTalk/Figures/DisElapsed.pdf")
plot(d,durOnly,type="n",ylab="Years disabled",
		xlab = "Elapsed years in episode",cex.lab=1.5,las=1,xlim=c(0,20))
polygon(c(d,rev(d)),c(durOnly,rep(0,31)),col=gray(.5),lwd=2)
dev.off()



durInAge <- acast(E2D, Age~dur,sum,value.var = "pi")
durInAgeD <- t(apply(durInAge, 2, cumsum))
durInAgeD <- cbind(0,durInAgeD)
plot(age, ageOnly, col = "red",lty=1,lwd=2,type='n',xlab = "Age",
		ylab = "Years sick",las=1,cex.lab=1.5,ylim=c(0,1.4))
for (i in 2:31){
	polygon(c(age[-1],rev(age[-1])),c(durInAgeD[,i],rev(durInAgeD[,i+1])),col=cols[i],border="white",lwd=.5)
}
range(durInAgeD)
par(mfrow=c(1,2))
matplot(durs, t(durInAge), type = 'l', xlim = c(0, 20), 
		col = gray(.7), lty = 1,lwd=1,
		xlab = "duration", ylab = "years")
matplot(durs, t(durInAge[age%%10==0,]), type = 'l', xlim = c(0, 20), 
		col = paste0(rev(rain(7)),"AA"), lty = 1,lwd=3,add=TRUE)

matplot(durs, t(durInAge / rowSums(durInAge, na.rm = TRUE)), type = 'l', xlim = c(0, 20), 
		col = rev(blues(31)), lty = 1,lwd=1)


