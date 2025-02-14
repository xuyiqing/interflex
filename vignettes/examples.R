### Interflex Demo Code
# v.1.3.0
# please report bugs to Yiqing Xu (yiqingxu@stanford.edu)

## installation
## install.packages('devtools', repos = 'http://cran.us.r-project.org') # if not already installed
## devtools::install_github('xuyiqing/interflex')

library(interflex)
data(interflex)
ls()

## plotting raw data

interflex(estimator = 'raw', Y = "Y", D = "D", X = "X", data = s1, weights = NULL, 
	Ylabel = "Outcome", Dlabel = "Treatment", Xlabel="Moderator", 
	main = "Raw Plot", cex.main = 1.2)

interflex(estimator = 'raw', Y = "Y", D = "D", X = "X", data = s2, 
	Ylabel = "Outcome", Dlabel = "Treatment", Xlabel="Moderator", 
	theme.bw = TRUE, show.grid = FALSE)

interflex(estimator = 'raw', Y = "Y", D = "D", X = "X", data = s3, 
          Ylabel = "Outcome", Dlabel = "Treatment", Xlabel="Moderator")

interflex(estimator = 'raw', Y="Y", D="D", X="X", Z=c("Z1"), data=s2)

# $$\ $$\                                         
# $$ |\__|                                        
# $$ |$$\ $$$$$$$\   $$$$$$\   $$$$$$\   $$$$$$\  
# $$ |$$ |$$  __$$\ $$  __$$\  \____$$\ $$  __$$\ 
# $$ |$$ |$$ |  $$ |$$$$$$$$ | $$$$$$$ |$$ |  \__|
# $$ |$$ |$$ |  $$ |$$   ____|$$  __$$ |$$ |      
# $$ |$$ |$$ |  $$ |\$$$$$$$\ \$$$$$$$ |$$ |      
# \__|\__|\__|  \__| \_______| \_______|\__|      
                                                
                                             

### the binning estimator

out <- interflex(Y = "Y", D = "D", X = "X", Z = "Z1", data = s1, 
	estimator = 'binning', vcov.type = "robust", 
	main = "Marginal Effects", ylim = c(-15, 15))
out$figure

out2 <- interflex(Y = "Y", D = "D", X = "X", Z = "Z1", data = s1, 
	estimator = 'linear', vcov.type = "robust", 
	main = "Marginal Effects", ylim = c(-15, 15))
out2$figure

predict(out2)

plot(out, xlab = "Moderate is X", Xdistr = "none", 
     bin.labs = FALSE, cex.axis = 0.8, cex.lab = 0.8)

out <- interflex(Y = "Y", D = "D", X = "X", Z = "Z1", data = s1, 
	estimator = 'binning', nbins = 4, theme.bw = TRUE, Xunif = TRUE)
out$figure

out <- interflex(Y = "Y", D = "D", X = "X", Z = "Z1", data = s2, 
	estimator = 'binning', Xdistr = "density", bin.labs = FALSE)
out$figure

out <- interflex(Y = "Y", D = "D", X = "X", Z = "Z1", data = s2, 
	estimator = 'binning', cutoffs = c(1,2,4,5))
out$figure

out <- interflex(Y = "Y", D = "D", X = "X", Z = "Z1", 
                 data = s3, estimator = 'binning')
out$tests$p.wald

## the kernel estimator

out <- interflex(Y = "Y", D = "D", X = "X", Z = "Z1", 
                 data = s1, estimator = 'kernel', 
	nboots = 200, parallel = TRUE, cores = 4)

predict(out)

out <- interflex(Y = "Y", D = "D", X = "X", Z = "Z1", 
                 data = s1, estimator = 'kernel',
                 nboots = 200, bw = 1, main = "The Linear Case", 
                 xlim = c(-0.5,6.5), ylim = c(-15, 15))
out$figure

out <- interflex(Y = "Y", D = "D", X = "X", Z="Z1", data = s1, 
                 estimator = 'kernel', nboots = 200, parallel = TRUE, 
                 cores = 4, Xunif = TRUE)
out$figure

out <- interflex(Y = "Y", D = "D", X = "X", Z = "Z1", data = s3, 
                 estimator = 'kernel', theme.bw = TRUE)
out$figure

plot(out, main = "Nonlinear Marginal Effects", ylab = "Coefficients", 
     Xdistr = "density", xlim = c(-3,3), ylim = c(-10,12), CI = FALSE, 
     cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)

predict(out)

### fixed effects

library(ggplot2)
ggplot(s4, aes(x=group, y = Y, colour = as.factor(D))) + 
  geom_point() + guides(colour=FALSE) 

interflex(estimator = "raw", Y = "Y", D = "D", X = "X", 
          data = s4, weights = NULL)

interflex(estimator = "raw", Y = "Y", D = "D", X = "X", Z="Z1", 
          FE=c("unit","year"),data = s4, weights = NULL)

s4.binning <- interflex(Y = "Y", D = "D", X = "X", Z = "Z1", data = s4, 
	estimator = 'binning', FE = NULL, cl = "unit")
s4.binning$figure

s4$wgt <- 1
s4.binning <- interflex(Y = "Y", D = "D", X = "X", Z = "Z1", data = s4, 
	estimator = 'binning', FE = c("unit","year"), cl = "unit", weights = "wgt")
s4.binning$figure

s4.kernel <- interflex(Y = "Y", D = "D", X = "X", Z = "Z1", data = s4, 
	estimator = 'kernel', FE = NULL, cl = "unit", weights = "wgt")
s4.binning$figure

s4.kernel <- interflex(Y = "Y", D = "D", X = "X", Z = "Z1", data = s4, 
	estimator = 'kernel', FE = c("unit","year"), cl = "unit")
s4.binning$figure

s4.kernel <- interflex(Y = "Y", D = "D", X = "X", Z = "Z1", data = s4, 
                       estimator = 'kernel', bw = 0.62, FE = c("unit","year"), 
                       cl = "unit", CI = FALSE, ylim = c(-9, 15), theme.bw = TRUE)
s4.binning$figure

### multiple treatment arms

interflex(estimator = "raw", Y = "Y", D = "D", X = "X", data = s5)

s5.binning <- interflex(Y = "Y", D = "D", X = "X", Z = c("Z1","Z2"), 
                        data = s5, estimator = 'binning',vartype = 'bootstrap')
s4.binning$figure

plot(s5.binning, order=c("C","B"), subtitles=c("Group C","Group B"))

plot(s5.binning, order=c("C","B"), subtitles=c("Control Group","Group C","Group B"),
	pool=TRUE, color = c("Salmon","#7FC97F", "#BEAED4"))

s5.binning2 <- interflex(Y = "Y", D = "D", X = "X",Z = c("Z1","Z2"), data = s5, 
	estimator = 'binning', base='B',vartype = 'bootstrap')
s5.binning2$figure

predict(s5.binning2, order = c('A','B','C'), 
        subtitles = c("Base Group","Group B","Group C"), 
        theme.bw = TRUE, show.grid = FALSE)

s5.kernel <-interflex(Y = "Y", D = "D", X = "X", Z = c("Z1","Z2"),
                      data = s5, estimator = 'kernel')
s5.kernel$figure

predict(s5.kernel, order = c('A','B','C'), subtitles = c("Group A (Base)","Group B","Group C"), 
	theme.bw = TRUE)

predict(s5.kernel, order = c('B','C','A'), subtitles = c("Group B","Group C","Group A"),
	pool = TRUE)


### testing
s5.kernel <-interflex(Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"), 
	data = s5, bw = 0.35, estimator = 'kernel', diff.values = c(-2,0))
s5.kernel$figure

s5.kernel$diff.estimate

dim(s5.kernel$vcov.matrix[[1]]) # why length 2

inter.test(s5.kernel,diff.values=c(0,1))

inter.test(s5.kernel,diff.values=c(-1,0,2))

inter.test(s5.kernel,diff.values=c(-1,0))

#       $$\ $$\                                          $$\               
#       $$ |\__|                                         $$ |              
#  $$$$$$$ |$$\  $$$$$$$\  $$$$$$$\  $$$$$$\   $$$$$$\ $$$$$$\    $$$$$$\  
# $$  __$$ |$$ |$$  _____|$$  _____|$$  __$$\ $$  __$$\\_$$  _|  $$  __$$\ 
# $$ /  $$ |$$ |\$$$$$$\  $$ /      $$ |  \__|$$$$$$$$ | $$ |    $$$$$$$$ |
# $$ |  $$ |$$ | \____$$\ $$ |      $$ |      $$   ____| $$ |$$\ $$   ____|
# \$$$$$$$ |$$ |$$$$$$$  |\$$$$$$$\ $$ |      \$$$$$$$\  \$$$$  |\$$$$$$$\ 
#  \_______|\__|\_______/  \_______|\__|       \_______|  \____/  \_______|
                                                                         

# GLM
s6.linear <- interflex(estimator='linear', method='logit', 
                       Y = "Y", D = "D", X = "X",Z = c("Z1", "Z2"), data = s6)

predict(s6.linear, pool = T, subtitles=c('Group:0','Group:1'))

s6.linear.a <- interflex(estimator='linear',method='logit',data=s6,Y='Y',D='D',diff.values = c(-2,0,2), X='X',Z=c('Z1','Z2'))

plot.interflex(s6.linear.a,diff.values = c(-2,0,2))

s6.linear.d <- interflex(estimator='linear', method='logit', Y = "Y", D = "D", X = "X",Z = c("Z1", "Z2"), data = s6, Z.ref=c(1,1))

plot(s6.linear.d)

s6.linear.f <- interflex(estimator='linear', method='logit', 
                         Y = "Y", D = "D", X = "X",Z = c("Z1", "Z2"), data = s6, full.moderate=TRUE)
plot(s6.linear.f)

s7.linear <- interflex(estimator='linear',method='logit',data=s7,Y='Y',D='D',X='X',Z=c('Z1','Z2'))

plot(s7.linear)

s7.linear.a <- interflex(estimator='linear',method='logit',data=s7,Y='Y',D='D',X='X',Z=c('Z1','Z2'),D.ref=c(0.25,0.5,0.75))

plot(s7.linear.a)

s7.gam <- interflex(estimator='gam',method='logit',data=s7,Y='Y',D='D',X='X',Z=c('Z1','Z2'),CI=FALSE)

s8.linear <- interflex(estimator='linear',data=s8,Y='Y',D='D', X='X',Z=c('Z1','Z2'),method='logit')


## binning estimator

s6.binning <- interflex(estimator='binning',data=s6,Y='Y',D='D', X='X',Z=c('Z1','Z2'),method='logit')

plot(s6.binning)

s8.binning <- interflex(estimator='binning',data=s8,Y='Y',D='D', X='X',Z=c('Z1','Z2'),method='logit')

plot(s8.binning)

predict(s8.binning)

s8.binning$tests

s8.binning.full <- interflex(estimator='binning',data=s8,Y='Y',D='D', X='X',Z=c('Z1','Z2'),method='logit',full.moderate=TRUE)

plot(s8.binning.full)

# GLM Kernel
s8.kernel <- interflex(estimator='kernel',data=s8,Y='Y',D='D', X='X',Z=c('Z1','Z2'),method='logit')

plot(s8.kernel)

s8.kernel.plot <- plot.interflex(s8.kernel,show.all = T,Xdistr = 'none')$`1`

s6.kernel <- interflex(estimator='kernel',data=s6,Y='Y',D='D', 
                       X='X',Z=c('Z1','Z2'),method='logit')

plot.interflex(s6.kernel)

predict(s8.kernel)

s8.kernel.full <- interflex(estimator='kernel',data=s8,Y='Y',D='D', 
           X='X',Z=c('Z1','Z2'),method='logit',full.moderate=TRUE)

plot(s8.kernel.full)

## Multiple Arms(GLM)
s9.linear <- interflex(estimator='linear',data=s9,Y='Y',D='D',X='X',
                       Z=c('Z1','Z2'),method='logit')
plot(s9.linear)

s9.binning <- interflex(estimator='binning',data=s9,Y='Y',D='D', 
                        X='X',Z=c('Z1','Z2'),method='logit')
plot(s9.binning)

s9.kernel <- interflex(estimator='kernel',data=s9,Y='Y',D='D',
                       X='X',Z=c('Z1','Z2'),method='logit')
plot(s9.kernel)

s9.kernel$Avg.estimate

#       $$\               $$\ 
#       $$ |              $$ |
#  $$$$$$$ |$$$$$$\$$$$\  $$ |
# $$  __$$ |$$  _$$  _$$\ $$ |
# $$ /  $$ |$$ / $$ / $$ |$$ |
# $$ |  $$ |$$ | $$ | $$ |$$ |
# \$$$$$$$ |$$ | $$ | $$ |$$ |
#  \_______|\__| \__| \__|\__|


library(ggplot2)
library(ggpubr)
library(interflex)
data(interflex)
ls()

#############################################
## Binary treatment with discrete outcomes
#############################################

s6.DML.nn <- interflex(estimator="DML", data = s6, ml_method="nn",
    Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"), 
	  treat.type = "discrete", base = 0)


## printing the first 10 rows
lapply(s6.DML.nn$est.dml, head, 10)

x <- s6$X
TE <- exp(-1+1+x+1*x)/(1+exp(-1+1+x+1*x))-exp(-1+0+x+0*x)/(1+exp(-1+0+x+0*x))
plot.s6 <- plot(s6.DML.nn, show.all = TRUE, Xdistr = 'none')$`1`
plot.s6 + geom_line(aes(x=x,y=TE),color='red')

#############################################
##Discrete treatment with discrete outcomes
#############################################

### true treatment effects
TE1 <- exp(1+x)/(1+exp(1+x))-exp(-x)/(1+exp(-x))
TE2 <- exp(4-x*x-x)/(1+exp(4-x*x-x))-exp(-x)/(1+exp(-x))

## neural network
s9.DML.nn <- interflex(estimator="DML", data = s9, ml_method="nn",
    Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
    treat.type = "discrete")

plot.s9.nn.1<-plot(s9.DML.nn, show.all = T)$`1`
plot.s9.nn.2<-plot(s9.DML.nn, show.all = T)$`2`

## adding true treatment effects to the plots 
p1.nn<-plot.s9.nn.1 + geom_line(aes(x=x,y=TE1),color='red')
p2.nn<-plot.s9.nn.2 + geom_line(aes(x=x,y=TE2),color='red')
p1.nn
p2.nn

## random forest
s9.DML.rf <- interflex(estimator="DML", data = s9, ml_method="rf",
   Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
   treat.type = "discrete")

plot.s9.rf.1 <- plot(s9.DML.rf, show.all = T)$`1`
plot.s9.rf.2 <- plot(s9.DML.rf, show.all = T)$`2`

p1.rf <- plot.s9.rf.1 + geom_line(aes(x=x,y=TE1),color='red')
p2.rf <- plot.s9.rf.2 + geom_line(aes(x=x,y=TE2),color='red')
p1.rf
p2.rf

## hist gradient boosting
s9.DML.hgb <- interflex(estimator="DML", data = s9, ml_method="hgb",
   Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"),
   treat.type = "discrete")

plot.s9.hgb.1 <- plot(s9.DML.hgb, show.all = T)$`1`
plot.s9.hgb.2 <- plot(s9.DML.hgb, show.all = T)$`2`

p1.hgb <- plot.s9.hgb.1 + geom_line(aes(x=x,y=TE1),color='red')
p2.hgb <- plot.s9.hgb.2 + geom_line(aes(x=x,y=TE2),color='red')
p1.hgb
p2.hgb

p1 <- ggarrange(p1.nn, p2.nn)
p1 <- annotate_figure(p1, top = text_grob("Neural Network",
    face = "bold", size = 13))

p2 <- ggarrange(p1.rf, p2.rf)
p2 <- annotate_figure(p2, top = text_grob("Random Forest",
    face = "bold", size = 13))

p3 <- ggarrange(p1.hgb, p2.hgb)
p3 <- annotate_figure(p3, top = text_grob("Hist Gradient Boosting",
    face = "bold", size = 13))

# put all together
ggarrange(p1, p2, p3, ncol = 1)


s9.DML.nn.2 <- interflex(estimator="DML", data = s9, ml_method="nn",
   Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"), solver = "lbfgs",
   max_iter = 10000, alpha = 1e-5, hidden_layer_sizes = c(3, 3, 2),
   random_state = 1, treat.type = "discrete")

plot(s9.DML.nn.2)

#############################################
## Continuous treatment with discrete outcome
#############################################

library(interflex)
data(interflex)

n <- nrow(s7)
d2 <- s7$D
x_values <- runif(n,min=-3, max = 3) 
d2_value <- median(d2)  # Example value for d2

marginal_effect <- function(D_value, X_value) {
  link_value <- -1 + D_value + X_value + D_value * X_value
  prob_value <- exp(link_value) / (1 + exp(link_value))
  marginal_effect_value <- prob_value * (1 - prob_value) * (1 + X_value)
  return(marginal_effect_value)
}

# Applying the function to the range of x values to calculate true ME
ME <- sapply(x_values, function(x_val) marginal_effect(d2_value, x_val))

s7.DML.nn.1 <- interflex(estimator="DML", data = s7, ml_method="nn",
   Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"), 
   dml_method = "default", treat.type = "continuous")

s7.DML.nn.2 <- interflex(estimator="DML", data = s7, ml_method="nn",
   Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"), 
   dml_method = "regularization", lasso_alpha = 1e-10,
   treat.type = "continuous")

s7.DML.nn.3 <- interflex(estimator="DML", data = s7, ml_method="nn", 
   Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"), 
   dml_method = "polynomial",
   treat.type = "continuous")

s7.DML.nn.4 <- interflex(estimator="DML", data = s7, ml_method="nn", 
   Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"), 
   dml_method = "non-parametric", 
   treat.type = "continuous")


s7.DML.nn.1.p <- plot(s7.DML.nn.1, show.all = TRUE)$`D=0.53 (50%)`
s7.DML.nn.2.p <- plot(s7.DML.nn.2, show.all = TRUE)$`D=0.53 (50%)`
s7.DML.nn.3.p <- plot(s7.DML.nn.3, show.all = TRUE)$`D=0.53 (50%)`
s7.DML.nn.4.p <- plot(s7.DML.nn.4, show.all = TRUE)$`D=0.53 (50%)`


## adding true treatment effects when d = 0.53 (50%) to the plots 
p1.nn<-s7.DML.nn.1.p + geom_line(aes(x=x_values,y=ME),color='red')+
  ylim(-0.5, 0.5) +  ggtitle("dml method: default")+
  theme(plot.title = element_text(size=10))
p1.nn

p2.nn <- s7.DML.nn.2.p +  geom_line(aes(x=x_values,y=ME),color='red')+
  ylim(-0.5, 0.5)+  ggtitle("dml method: polynomial")+
  theme(plot.title = element_text(size=10))
p2.nn

p3.nn <- s7.DML.nn.3.p +geom_line(aes(x=x_values,y=ME),color='red')+
  ylim(-0.5, 0.5) + ggtitle("dml method: regularization")+
  theme(plot.title = element_text(size=10))
p3.nn

p4.nn <- s7.DML.nn.4.p +  geom_line(aes(x=x_values,y=ME),color='red')+
  ylim(-0.5, 0.5) +  ggtitle("dml method: non-parametric")+
  theme(plot.title = element_text(size=10))
p4.nn

ggarrange(p1.nn, p2.nn, p3.nn, p4.nn, ncol = 2, nrow = 2)



s7.DML.nn.lasso1<- interflex(estimator="DML", data = s7, ml_method="nn",
   Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"), 
   dml_method = "regularization", lasso_alpha = 1e-4,
   treat.type = "continuous")

s7.DML.nn.lasso2<- interflex(estimator="DML", data = s7, ml_method="nn",
   Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"), 
   dml_method = "regularization", 
   lasso_alpha = 1e-10, 
                       treat.type = "continuous")

s7.DML.nn.poly1 <- interflex(estimator="DML", data = s7, ml_method="nn", 
   Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"), 
   dml_method = "polynomial", poly_degree = 3,
   treat.type = "continuous")

s7.DML.nn.poly2 <- interflex(estimator="DML", data = s7, ml_method="nn", 
   Y = "Y", D = "D", X = "X", Z = c("Z1", "Z2"), 
   dml_method = "polynomial", poly_degree = 6,
   treat.type = "continuous")


s7.DML.nn.lasso1.p <- plot(s7.DML.nn.lasso1, show.all = TRUE)$`D=0.53 (50%)`
s7.DML.nn.lasso2.p <- plot(s7.DML.nn.lasso2, show.all = TRUE)$`D=0.53 (50%)`
s7.DML.nn.poly1.p <- plot(s7.DML.nn.poly1, show.all = TRUE)$`D=0.53 (50%)`
s7.DML.nn.poly2.p <- plot(s7.DML.nn.poly2, show.all = TRUE)$`D=0.53 (50%)`


## adding true treatment effects when d = 0.53 (50%) to the plots 
p1.nn<-s7.DML.nn.lasso1.p + geom_line(aes(x=x_values,y=ME),color='red')+
  ylim(-0.5, 0.5) + ggtitle("regularization, default alpha")+
  theme(plot.title = element_text(size=10))
p1.nn

p2.nn<-s7.DML.nn.lasso2.p + geom_line(aes(x=x_values,y=ME),color='red')+
  ylim(-0.5, 0.5) + ggtitle(label = "regularization, alpha = 1e-10")+
  theme(plot.title = element_text(size=10))
p2.nn

p3.nn<-s7.DML.nn.poly1.p + geom_line(aes(x=x_values,y=ME),color='red')+
  ylim(-0.5, 0.5) + ggtitle(label = "polynomial, default degree")+
  theme(plot.title = element_text(size=10))
p3.nn

p4.nn<-s7.DML.nn.poly2.p +  geom_line(aes(x=x_values,y=ME),color='red')+
  ggtitle(label = "polynomial, degree = 6") + 
  theme(plot.title = element_text(size=10))
p4.nn

ggarrange(p1.nn, p2.nn, p3.nn, p4.nn, ncol = 2, nrow = 2)

################################################
# Application 1: Huddy, Mason, and Aarøe (2015)
## Binary treatment, continuous outcome
################################################

d <- app_hma2015
head(d)
d$pidstr2_threat<-d$pidstr2 * d$threat
d$issuestr2_threat<-d$issuestr2 * d$threat

Y <- "totangry" 
D <- "threat" 
X <- "pidentity"
Z <- c("issuestr2", "pidstr2", "pidstr2_threat" ,"issuestr2_threat", "knowledge" , "educ" , "male" , "age10" )
Dlabel <- "Threat"
Xlabel <- "Partisan Identity"
Ylabel <- "Anger"
vartype <- "robust"
main <- "Huddy et al. (2015) \n APSR"
cl <- cuts <- cuts2 <- time <- NULL

## neural nets
out.dml.nn <- interflex(estimator='DML', data = d, ml_method="nn",
    Y=Y,D=D,X=X, Z = Z, treat.type = "discrete",  
    Xlabel=Xlabel, Ylabel=Ylabel, Dlabel=Dlabel, ylim=c(-0.18,0.6))

# random forest
out.dml.rf <- interflex(estimator='DML', data = d, ml_method="rf",
    Y=Y,D=D,X=X, Z = Z, treat.type = "discrete",
    Xlabel=Xlabel, Ylabel=Ylabel, Dlabel=Dlabel, ylim=c(-0.18,0.6))

## hist gradient boosting 
out.dml.hgb <- interflex(estimator='DML', data = d, ml_method="hgb",
    Y=Y,D=D,X=X, Z = Z, treat.type = "discrete",
    Xlabel=Xlabel, Ylabel=Ylabel, Dlabel=Dlabel,ylim=c(-0.18,0.6))

## raw data plot
out.raw <- interflex(estimator = "raw", Y=Y, D=D, X=X, data=d, 
    Xlabel=Xlabel, Ylabel=Ylabel, Dlabel=Dlabel)

## linear/binning estimator
out.est1 <- interflex(estimator = "binning",Y=Y,D=D,X=X,Z=Z,data=d,
    Xlabel=Xlabel, Ylabel=Ylabel, Dlabel=Dlabel,
    nbins=3, cutoffs=cuts,  cl=cl, time=time,
    pairwise=TRUE,  Xdistr = "histogram",ylim=c(-0.18,0.6))

## kernel estimator
out.kernel <- interflex(estimator = "kernel", data=d, Y=Y,D=D,X=X,Z=Z,
    cl=NULL, Dlabel=Dlabel, Xlabel=Xlabel, Ylabel=Ylabel,
    bw=0.917, Xdistr = "histogram", ylim=c(-0.18,0.6))

## put all together
ggarrange(out.raw, out.est1$figure, out.kernel$figure,
    out.dml.nn$figure, out.dml.rf$figure, out.dml.hgb$figure, 
    common.legend = TRUE, legend = "none",
    labels = c("Raw", "Binning", "Kernel", "DML: nn","DML: rf", "DML: hgb"))

###########################################
# Application 2: Vernby (2013)
## Continuous treatment, continuous outcome
###########################################

d <- app_vernby2013
Y <- "school_diff" 
D <- "noncitvotsh" 
X <- "noncit15" 
Z <- c("Taxbase2" ,  "Taxbase2_2" , "pop" , "pop_2" ,   "manu" , "manu_2")
Dlabel <- "Share NC"
Xlabel <- "Prop. School-Aged NC"
Ylabel <- "Δ Ed. Services"
vartype <- "robust"
name<-"vernby_2013a"
main<-"Vernby (2013) \n AJPS"
time <- cl <- cuts <- cuts2 <- NULL

head(d)

out.dml.nn<-interflex(estimator='DML', data = d, ml_method="nn", 
    Y=Y,D=D,X=X, Z = Z, treat.type = "continuous", dml_method = "default",  
    Xlabel=Xlabel, Ylabel=Ylabel, Dlabel=Dlabel)

lapply(out.dml.nn$est.dml, head, 20)

## neural nets
out.dml.nn <- interflex(estimator='DML', data = d, ml_method="nn", 
    Y=Y,D=D,X=X, Z = Z, treat.type = "continuous", dml_method = "non-parametric",
    Xlabel=Xlabel, Ylabel=Ylabel, Dlabel=Dlabel)

# random forest
out.dml.rf <- interflex(estimator='DML', data = d, ml_method="rf",
    Y=Y,D=D,X=X, Z = Z, treat.type = "continuous", dml_method = "non-parametric",
    Xlabel=Xlabel, Ylabel=Ylabel, Dlabel=Dlabel)

## hist gradient boosting 
out.dml.hgb<-interflex(estimator='DML', data = d, ml_method="hgb",
   Y=Y,D=D,X=X, Z = Z, treat.type = "continuous", dml_method = "non-parametric",
   Xlabel=Xlabel, Ylabel=Ylabel, Dlabel=Dlabel)

# raw data plot
out.raw <- interflex(estimator = "raw", Y=Y,D=D,X=X,data=d,Xlabel=Xlabel,
   Ylabel=Ylabel, Dlabel=Dlabel, cutoffs=cuts,span=NULL)

## linear/binning estimator
out.est1 <- interflex(estimator = "binning",Y=Y,D=D,X=X,Z=Z,data=d,
  Xlabel=Xlabel, Ylabel=Ylabel, Dlabel=Dlabel,
  cutoffs=cuts,  cl=cl, time=time,
  pairwise=TRUE, Xdistr = "histogram")

## kernel estimator
out.kernel <- interflex(estimator = "kernel", data=d, Y=Y,D=D,X=X,Z=Z,
  cl=cl, Dlabel=Dlabel, Xlabel=Xlabel, Ylabel=Ylabel, Xdistr = "histogram")

## put all together
ggarrange(out.raw, out.est1$figure, out.kernel$figure,
  out.dml.nn$figure, out.dml.rf$figure, out.dml.hgb$figure, 
  common.legend = TRUE, legend = "none",
  labels = c("Raw", "Binning", "Kernel", "DML: nn","DML: rf","DML: hgb"))


