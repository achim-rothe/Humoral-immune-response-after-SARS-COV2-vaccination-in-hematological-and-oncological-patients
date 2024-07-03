
library(censReg)
library(plm)
library(emmeans)
library(ggplot2)


foo <- c(min = base::min, max = base::max, median = stats::median, mean = base::mean, sd = stats::sd)

###########################################################################
###                                                                     ###
###                               FUNCTIONS                             ###
###                                                                     ###
###########################################################################


stats <- function(formula, data, fun = c(base::min, base::max, base::mean, stats::sd)){
    funnames <- names(fun)
    res <- aggregate(formula, data = data, FUN = length, drop = FALSE)
    ny <- ncol(res)
    colnames(res)[ncol(res)] <- "N"
    NAs <- aggregate(formula, data = data, FUN = function(x) sum(is.na(x)), drop = FALSE, na.action = na.pass)
    res <- data.frame(res, "NA" = NAs[, ncol(NAs)], check.names=FALSE)
    rest <- sapply(fun, function(f) aggregate(formula, data = data, FUN = f, drop = FALSE)[,ny])
    colnames(rest) <- funnames
    res <- cbind(res, rest)
    return(res)
}


# Boxplots darstellt:
boxPlot <- function(formula, data, add.N = TRUE, ylab, ...){
    dat <- data
    if(missing(ylab)) ylab <- gsub("`", "", rownames(attr(terms(formula), "factors"))[attr(terms(formula), "response")])
    vars <- paste0(gsub("`", "", attr(terms(formula), "term.labels")), collapse="\n")

    boxplot(formula, dat, col=NA, xaxt="n", xlab="", ylab=ylab, xaxs="i", border = NA, yaxt = "n", frame = TRUE, ...)
    abline(h = axTicks(2), col = "lightgrey")
    b <- boxplot(formula, data = dat, sep = "\n", xaxt = "n", xlab = "", ylab = ylab, add = TRUE, ...)
    
    mtext(side = 1, line = 0, at = 1:length(b$names), text = b$names, padj = 1, ...)
    mtext(side = 1, line = 0, at = par("usr")[1], text = vars, padj = 1, adj = 1, ...)
    if(add.N){
        mtext(side = 3, line = 0.1, at = par("usr")[1], text = "N", padj = 0, adj = 1)
        mtext(side = 3, line = 0.1, at = 1:length(b$n), text = b$n, padj = 0, ...)
    }
    return(invisible(dat))
}


###################################################################
###                                                             ###
###  Chapter:                                                   ###
###  Differences in anti-SARS-CoV-2 antibody titers in          ###
###  hematological and oncological patients                     ###
###                                                             ###
###################################################################

load("Analyse 2 Dex Boosterimpfung20220908.RData", verbose=TRUE)
str(dat)
str(dat.long)


stats(`sample titer (BAU/ml)` ~ sample, data = dat.long, fun = foo[1:3])

stats(`2. sample titer (BAU/ml)` ~ disease, data = dat, fun = foo[1:3])

# N
stats(`2. sample titer (BAU/ml)` ~ disease, data = dat, 
       fun=c("titer \u2264 4.8 [%]" = function(x) sum(x == 4.8), 
       "titer (4.8, 2080) [%]" = function(x) sum(x > 4.8 & x < 2080),  
       "titer \u2265 2080 [%]" = function(x) sum(x == 2080)))

# %
stats(`2. sample titer (BAU/ml)` ~ disease, data = dat, 
       fun=c("titer \u2264 4.8 [%]" = function(x) sum(x == 4.8)*100, 
       "titer (4.8, 2080) [%]" = function(x) sum(x > 4.8 & x < 2080)*100,  
       "titer \u2265 2080 [%]" = function(x) sum(x == 2080)*100))
  
# table 1
stats(`2. sample titer (BAU/ml)` ~ `form of therapy`, data=dat, 
       fun=c("titer \u2264 4.8 [%]" = function(x) mean(x == 4.8)*100, 
       "titer (4.8, 2080) [%]" = function(x) mean(x > 4.8 & x < 2080)*100,  
       "titer \u2265 2080 [%]" = function(x) mean(x == 2080)*100))
  


###################################################################
###                                                             ###
###  Chapter:                                                   ###
###  Increase of anti-SARS-CoV-2 IgG following                  ###
###  the booster vaccination                                    ###
###                                                             ###
###################################################################

load("Analyse 2 Dex Boosterimpfung20220908.RData", verbose=TRUE)
str(dat)

stats(`sample titer (BAU/ml)`  ~ sample, data = dat.long, fun = foo[1:3])
boxPlot(`sample titer (BAU/ml)`  ~ sample, data = dat.long)

# disease
stats(`sample titer (BAU/ml)`  ~ disease + sample, data = dat.long, fun = foo[1:3])
boxPlot(`sample titer (BAU/ml)`  ~ disease + sample, data = dat.long)

# gender
stats(`sample titer (BAU/ml)`  ~ gender + sample, data = dat.long, fun = foo[1:3])
boxPlot(`sample titer (BAU/ml)`  ~ gender + sample, data = dat.long)

# glucocorticoids medication
stats(`sample titer (BAU/ml)`  ~ `glucocorticoids medication` + sample, data = dat.long, fun = foo[1:3])
boxPlot(`sample titer (BAU/ml)`  ~ `glucocorticoids medication` + sample, data = dat.long)

# Astrazeneca
stats(`sample titer (BAU/ml)`  ~ Astrazeneca + sample, data = dat.long, fun = foo[1:3])
boxPlot(`sample titer (BAU/ml)`  ~ Astrazeneca + sample, data = dat.long)

# Biontech
stats(`sample titer (BAU/ml)`  ~ Biontech + sample, data = dat.long, fun = foo[1:3])
boxPlot(`sample titer (BAU/ml)`  ~ Biontech + sample, data = dat.long)

# Johnson
stats(`sample titer (BAU/ml)`  ~ Johnson + sample, data = dat.long, fun = foo[1:3])
boxPlot(`sample titer (BAU/ml)`  ~ Johnson + sample, data = dat.long)

# Moderna
stats(`sample titer (BAU/ml)`  ~ Moderna + sample, data = dat.long, fun = foo[1:3])
boxPlot(`sample titer (BAU/ml)`  ~ Moderna + sample, data = dat.long)



# therapy
stats(`sample titer (BAU/ml)` ~ therapy + sample, data=dat.long, fun=foo[1:3])


fun <- function(dat.long, i, v, subset = levels(dat.long[,v]), cex3=1, ...){
    boxplot(`sample titer (BAU/ml)` ~ sample, data=dat.long[which(dat.long[,v] == subset[i]),,drop=FALSE], 
            ylim=range(dat.long$"sample titer (BAU/ml)", na.rm=TRUE), xlab="", ylab="", 
            xaxt = ifelse(i < nlevels(dat.long[,v]), "n", "s"))
    mtext(side=4, line=0, text=subset[i], las=0, cex=cex3, padj =1)
    mtext("sample titer (BAU/ml)", side=2, line=3.3, las=0, cex = 1)
}

op <- par(oma = c(4.5, 5, 0.5, 4), mar = c(0,0,0,0), las=1, mfrow=c(nlevels(dat.long[,"therapy"]),1))
for(i in 1:nlevels(dat.long[,"therapy"])) fun(dat.long, i = i, v ="therapy")
mtext("sample", side=1, line=2, las=0)
par(op)

# form of therapy
stats(`sample titer (BAU/ml)` ~ `form of therapy` + sample, data=dat.long, fun=foo[1:3])

op <- par(oma = c(3.5, 4.3, 0.5, 3), mar = c(0,0,0,0), las=1, mfrow=c(nlevels(dat.long[,"form of therapy"]),1))
for(i in 1:nlevels(dat.long[,"form of therapy"])) fun(dat.long, i = i, v ="form of therapy")
mtext("sample", side=1, line=2, las=0)
par(op)

# diagnose clustered
stats(`sample titer (BAU/ml)` ~ `diagnose clustered` + sample, data=dat.long, fun=foo[1:3])

op <- par(oma = c(3.5, 4.3, 0.5, 3), mar = c(0,0,0,0), las=1, mfrow=c(nlevels(dat.long[,"diagnose clustered"]),1))
for(i in 1:nlevels(dat.long[,"diagnose clustered"])) fun(dat.long, i = i, v ="diagnose clustered")
mtext("sample", side=1, line=2, las=0)
par(op)


###                                                                     ###
###                                   MODEL                             ###
###                                                                     ###

dat.long$pat <- dat.long$"patient-nr."
pdat <- pdata.frame(dat.long, c("pat", "sample" ))

summary(mod <- censReg(sample.titer..BAU.ml. ~ sample  + gender +  
                    glucocorticoids..mg. + form.of.therapy + disease, method = "BHHH", left = 4.8, right = 2080, data = pdat))

namenneu <- function(nam){
    nam <- gsub("..mg.", " (mg)", nam)
    nam <- gsub("form.of.therapy", "form of therapy", nam)
    nam <- gsub("diagnose.clustered", "diagnose clustered", nam)
    return(nam)
}

k <- coef(summary(mod))
rownames(k) <- namenneu(rownames(k))
w <- 1:(nrow(k)-2)


rg <- qdrg(sample.titer..BAU.ml. ~ sample + gender + glucocorticoids..mg. + form.of.therapy + disease, data=pdat, coef=coef(mod)[w],vcov=vcov(mod)[w,w], df=mod$df.residual,
            cov.reduce = function(x) min(x)+c(0,1))
            
joint_tests(rg)

# sample effect
emmeans(rg, revpairwise ~ sample, adjust = "bonferroni")

# therapy effect
emmeans(rg, revpairwise ~ form.of.therapy, adjust = "bonferroni")

# disease effect
emmeans(rg, revpairwise ~ disease, adjust = "bonferroni")


###########################################################################
###                                                                     ###
### Chapter:                                                            ###
### Anti-SARS-CoV-2 IgG titer Covid-19 negative patients                ###
### and breakthrough infections                                         ###
###                                                                     ###
###########################################################################

load("Covid zu No COVID.RData", verbose=TRUE)

stats(`SARS-CoV-2 IgG Titer (BAU/ml)` ~ COVID, data = covid, fun = foo[1:3])
boxplot(`SARS-CoV-2 IgG Titer (BAU/ml)` ~ COVID, data=covid)


stats(`SARS-CoV-2 IgG Titer (BAU/ml)` ~ COVID, data = covidSymp, fun = foo[1:3])
boxplot(`SARS-CoV-2 IgG Titer (BAU/ml)` ~ COVID, data=covidSymp)

###                                                                     ###
###                         MODEL                                       ###
###                                                                     ###

summary(modCovid <- censReg(`SARS-CoV-2 IgG Titer (BAU/ml)` ~ COVID, left = 4.8, right = 2080, data = covid))
mod <- modCovid

k <- coef(summary(mod))
w <- 1:(nrow(k)-1)

rg <- qdrg(`SARS-CoV-2 IgG Titer (BAU/ml)` ~ COVID, data=covid, coef=coef(mod)[w],vcov=vcov(mod)[w,w], df=mod$df.residual)

joint_tests(rg)


###########################################################################
###                                                                     ###
### Chapter:                                                            ###
### Course of anti-SARS-CoV-2 antibody titers before and after          ###
### Covid-19 infection in asymptomatic and symptomatic patients         ###
###                                                                     ###
###########################################################################

load("Covidfinalkorr2.RData", verbose=TRUE)
str(dat.long)

summary(dat.long$`sample titer (BAU/ml)`)

stats(`sample titer (BAU/ml)` ~ sample, data = dat.long, fun=foo[1:3])
boxplot(`sample titer (BAU/ml)` ~ sample, data = dat.long)

stats(`sample titer (BAU/ml)` ~ symptoms, data = dat.long, fun=foo[1:3])
boxplot(`sample titer (BAU/ml)` ~ symptoms, data = dat.long)

## Table 2:
stats(`sample titer (BAU/ml)` ~ symptoms + sample, data = dat.long, fun=foo[1:3])
boxPlot(`sample titer (BAU/ml)` ~ symptoms + sample, data = dat.long)


###                                                                     ###
###                         MODEL                                       ###
###                                                                     ###
dat.long$pat <- dat.long$"patient-nr."
pdat <- pdata.frame(dat.long, c("pat", "sample" ))

summary(mod <- censReg(sample.titer..BAU.ml. ~ sample * (symptoms +  COVID.duration..d.), 
            method = "BHHH", left = 4.8, right = 2080, data = pdat))

k <- coef(summary(mod))
w <- 1:(nrow(k)-2)

rg <- qdrg(sample.titer..BAU.ml. ~ sample * (symptoms + COVID.duration..d.), data=pdat, coef=coef(mod)[w],vcov=vcov(mod)[w,w], df=mod$df.residual,
            cov.reduce = function(x) min(x)+c(0,1))
            
joint_tests(rg)


# effect of sample and symptoms:
emmeans(rg, revpairwise ~ sample * symptoms, by="symptoms", adjust = "bonferroni")
emmeans(rg, revpairwise ~ sample * symptoms, by="sample", adjust = "bonferroni")
