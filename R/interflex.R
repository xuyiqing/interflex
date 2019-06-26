## Jens Hainmueller; Jonathan Mummolo; Yiqing Xu
## This Version: 2018.7.22

######   Interpreting Interaction Models  #######

## 1. inter.raw: first look at the data: D, X, Y
## 2. inter.gam: GAM plots
## 3. inter.binning: estimate
## 4. inter.kernel: kernel smooth plot
## 5. inter.plot: change looks of a plot based on 
##     an interflex object

## Supporting
## 1. coefs: kernel weighted least squares
## 2. cross-validation: cross validate kernel bandwidth
## 3. vcovCluster: clustered standard error


#################################################
coefs <- NULL
crossvalidate <- NULL
inter.binning <- NULL
inter.kernel <- NULL
inter.raw <- NULL
inter.gam <- NULL
inter.plot <- NULL
vcovCluster <- NULL
#################################################

inter.raw<-function(
    data,Y,D,X,weights=NULL,
    nbins=3,
    cutoffs=NULL, 
    span=NULL,
    pos=NULL,
    main = NULL,
    Ylabel=NULL,
    Dlabel= NULL,
    Xlabel=NULL,
    theme.bw = FALSE, 
    show.grid = TRUE,
    cex.main = NULL,
    cex.lab = NULL,
    cex.axis = NULL
    ){ 
    
    ## Y: outcome
    ## D: "treatment" indicator
    ## X: covariate to be interacted with D
    ## nbins: ## of bins
    ## cutoffs: specified cut-points
    ## span: bandwidth for loess

    ## in case data is in tibble format
    if (is.data.frame(data) == FALSE) {
        data <- as.data.frame(data)
    }

    ## check input
    if (is.character(Y) == FALSE) {
        stop("\"Y\" is not a string.")
    } else {
        Y <- Y[1]
    }
    if (is.character(D) == FALSE) {
        stop("\"D\" is not a string.")
    } else {
        D <- D[1]
    }
    if (is.character(X) == FALSE) {
        stop("\"X\" is not a string.")
    } else {
        X <- X[1]
    }
    if (is.null(weights) == FALSE) {
        if (is.character(weights) == FALSE) {
            stop("\"weights\" is not a string.")
        } else {
            weights <- weights[1]
        }   
    }
    if (is.null(nbins) == FALSE) {
        if (nbins%%1 != 0) {
            stop("\"nbins\" is not a positive integer.")
        } else {
            nbins <- nbins[1]
        }
        if (nbins < 1) {
            stop("\"nbins\" is not a positive integer.")
        }
    }
    if (is.null(cutoffs) == FALSE) {
        if (is.numeric(cutoffs) == FALSE) {
            stop("Some element of cutoffs is not numeric.")
        } 
    }
    if (is.null(span) == FALSE) {
        if (is.numeric(span) == FALSE) {
            stop("\"span\" is not numeric.")
        } else {
            span <- span[1]
            if (span <= 0) {
                stop("\"span\" is not a positive number.")
            }
        }
    } 
    if (is.null(Ylabel)==TRUE) {
        Ylabel <- Y
    } else {
        if (is.character(Ylabel) == FALSE) {
            stop("\"Ylabel\" is not a string.")
        } else {
            Ylabel <- Ylabel[1]
        }   
    } 
    if (is.null(Dlabel)==TRUE) {
        Dlabel <- D   
    } else {
        if (is.character(Dlabel) == FALSE) {
            stop("\"Dlabel\" is not a string.")
        } else {
            Dlabel <- Dlabel[1]
        }   
    }
    if (is.null(Xlabel)==TRUE) {
        Xlabel <- X   
    } else {
        if (is.character(Xlabel) == FALSE) {
            stop("\"Xlabel\" is not a string.")
        } else {
            Xlabel <- Xlabel[1]
        }   
    }

    ## font size
    if (is.null(cex.lab)==FALSE) {
        if (is.numeric(cex.lab)==FALSE) {
            stop("\"cex.lab\" is not numeric.")
        }
    }
    if (is.null(cex.main)==FALSE) {
        if (is.numeric(cex.main)==FALSE) {
            stop("\"cex.main\" is not numeric.")
        }
    }
    if (is.null(cex.axis)==FALSE) {
        if (is.numeric(cex.axis)==FALSE) {
            stop("\"cex.axis\" is not numeric.")
        }    
    }  

    ## title
    if (is.null(main)==FALSE) {
        main <- main[1]
    } 

     
    ## load packages
    requireNamespace("ggplot2")

    ## drop missing values
    data <- na.omit(data[,c(Y, D, X)])
     
    ## plotting
    if (length(unique(data[,D]))==2) { ## binary case
        
        ## plotting
        treat.lab<-c(paste(Dlabel,"= 0"),paste(Dlabel,"= 1"))
        data.aug<-data
        data.aug$treat<-factor(data[,D], labels=treat.lab)

        if (is.null(pos)==TRUE) {
            box.pos.co <- min(data[,Y])
            box.pos.tr <- max(data[,Y])
        } else {
            box.pos.co <- pos[1]
            box.pos.tr <- pos[2]
        }

        tr <- which(data[,D]==1)
        co <- which(data[,D]==0)
        qt90.tr <- quantile(data[tr,X],c(0.05,0.95))
        qt90.co <- quantile(data[co,X],c(0.05,0.95))
        qt50.tr <- quantile(data[tr,X],c(0.25,0.75))
        qt50.co <- quantile(data[co,X],c(0.25,0.75))
        med.tr <- median(data[tr,X])
        med.co <- median(data[co,X])
        
        data.aug$qt90 <- NA
        data.aug$qt90[which(data.aug[,D]==1 &
                            data.aug[,X]>=qt90.tr[1] &
                            data.aug[,X]<=qt90.tr[2])]<- box.pos.tr
        data.aug$qt90[which(data.aug[,D]==0 &
                            data.aug[,X]>=qt90.co[1] &
                            data.aug[,X]<=qt90.co[2])]<- box.pos.co
        data.aug$qt50 <- NA
        data.aug$qt50[which(data.aug[,D]==1 &
                            data.aug[,X]>=qt50.tr[1] &
                            data.aug[,X]<=qt50.tr[2])]<- box.pos.tr
        data.aug$qt50[which(data.aug[,D]==0 &
                            data.aug[,X]>=qt50.co[1] &
                            data.aug[,X]<=qt50.co[2])]<- box.pos.co
        data.aug$med <- NA
        data.aug$med[tr[which.min(abs(data.aug[tr,X]-med.tr))]]<- box.pos.tr
        data.aug$med[co[which.min(abs(data.aug[co,X]-med.co))]]<- box.pos.co
        
        
        ## plotting
        if (is.null(weights)==TRUE) {
            p1 <- ggplot(data.aug, aes_string(X, Y))
        } else {
            p1 <- ggplot(data.aug, aes_string(X, Y, weight=weights))
                
        }
        if (theme.bw == TRUE) {
            p1 <- p1 + theme_bw() 
        }
        if (show.grid == FALSE) {
            p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }
        p1 <- p1 + geom_point() +
            geom_smooth(method = "lm", se = F, fullrange = T,
                        colour = "steelblue", size = 1)
        
        if (is.null(span)==TRUE) {
            p1 <- p1+ geom_smooth(method = "loess", formula = y ~ x,
                            se = F, colour="red") 
        } else {
            p1 <- p1 + geom_smooth(method = "loess", formula = y ~ x,
                            se = F, colour="red",span=span)
        }
        p1 <- p1 + xlab(Xlabel) + ylab(Ylabel)

       
        p1 <- p1 + geom_line(aes_string(X,"qt90"),
                             size=1,colour="grey50", na.rm = TRUE)
        p1 <- p1 + geom_line(aes_string(X,"qt50"),
                             size=3,colour="grey50", na.rm = TRUE)
        p1 <- p1 + geom_point(aes_string(X,"med"),size=3,
                              shape=21,fill="white",colour="red", na.rm = TRUE)


        p1 <- p1 + facet_wrap(~treat, ncol=1) 
        #p1 <- p1 + facet_grid(treat ~.)                
        
        
        
    } else { # continuous case
        
        ## grouping by X
        if (is.null(cutoffs)==TRUE) {
            cutoff<-quantile(data[,X],probs=seq(0,1,1/nbins))
            while (length(unique(cutoff))!=nbins+1) {
                nbins<-nbins-1
                cutoff<-quantile(data[,X],probs=seq(0,1,1/nbins))
            } 
        } else {
            cutoffs <-cutoffs[which(cutoffs>min(data[,X]) & cutoffs < max(data[,X]))]
            cutoff<- sort(unique(c(min(data[,X]),cutoffs,max(data[,X]))))
        } 
        groupID<-cut(data[,X],breaks=cutoff, labels = FALSE)
        groupID[which(data[,X]==min(data[,X]))]<-1
        
        ## X labels
        groupID2 <- cut(data[,X],breaks=cutoff)
        gp.lab = paste(Xlabel, ": ", levels(groupID2), sep="")
        gp.lab[1] <- paste(Xlabel, ": [", substring(levels(groupID2)[1],2), sep = "")
        nbins <- length(unique(groupID))
        
        ## if (nbins==2) {
        ##     gp.lab<-c(paste(Xlabel,": low",sep=""),paste(Xlabel,": high",sep="")) 
        ## } else if (nbins==3) {
        ##     gp.lab<-c(paste(Xlabel,": low",sep=""),paste(Xlabel,": medium",sep=""),
        ##               paste(Xlabel,": high",sep="")) 
        ## } else {
        ##     gp.lab<-c();
        ##     for (i in 1:nbins) {
        ##         gp.lab<-c(gp.lab,paste("Grp",i))
        ##     }
        ## }
        
        groupID <- factor(groupID, labels=gp.lab)
        data.aug <- data
        data.aug$groupID<-groupID
         
        ## plotting
        if (is.null(weights)==TRUE) {
            p1 <- ggplot(data.aug, aes_string(D, Y))
        } else {
            p1 <- ggplot(data.aug, aes_string(D, Y,weight=weights))
        }
        if (theme.bw == TRUE) {
            p1 <- p1 + theme_bw() 
        }
        if (show.grid == FALSE) {
            p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        }
        p1 <- p1 + geom_point() + 
            geom_smooth(method = "lm", se = F, fullrange = T,
                        colour = "steelblue", size = 1)
        if (is.null(span)==TRUE) {
            p1 <- p1 +
                geom_smooth(method = "loess", formula = y ~ x,
                            se = F, colour="red") 
        } else {
            p1 <- p1 +
                geom_smooth(method = "loess", formula = y ~ x,
                            se = F, colour="red",span=span) 
        }        
        p1 <- p1 + xlab(Dlabel) + ylab(Ylabel) + facet_grid(.~groupID)
        
    }

    ## axis labels
    if (is.null(cex.lab)==TRUE) {
        cex.lab <- 15
    } else {
        cex.lab <- 15 * cex.lab
    }
    if (is.null(cex.axis)==TRUE) {
        cex.axis <- 15
    } else {
        cex.axis <- 15 * cex.axis
    }
    p1 <- p1 + theme(axis.text = element_text(size=cex.axis), axis.title = element_text(size=cex.lab))

    ## title
    if (is.null(cex.main)==TRUE) {
        cex.main <- 18
    } else {
        cex.main <- 18 * cex.main
    }
    if (is.null(main)==FALSE) {
            p1<-p1 + ggtitle(main) +
            theme(plot.title = element_text(hjust = 0.5, size=cex.main,lineheight=.8, face="bold"))
    } 

    return(list(graph = p1)) 
}


#############################################################
## GAM


inter.gam<-function(data,Y,D,X,
    Z=NULL,
    weights=NULL,
    FE=NULL,
    SE = FALSE,
    k=10,
    angle=c(30, 100,-30,-120),
    Ylabel = NULL, Dlabel = NULL, Xlabel = NULL){

    ## in case data is in tibble format
    if (is.data.frame(data) == FALSE) {
        data <- as.data.frame(data)
    }

    ## check input
    if (is.character(Y) == FALSE) {
        stop("\"Y\" is not a string.")
    } else {
        Y <- Y[1]
    }
    if (is.character(D) == FALSE) {
        stop("\"D\" is not a string.")
    } else {
        D <- D[1]
    }
    if (is.character(X) == FALSE) {
        stop("\"X\" is not a string.")
    } else {
        X <- X[1]
    }
    if (is.null(Z) == FALSE) {
        for (i in 1:length(Z)) {
            if (is.character(Z[i]) == FALSE) {
                stop("Some element in Z is not a string.")
            }
        }
    }
    if (is.null(FE) == FALSE) {
        for (i in 1:length(FE)) {
            if (is.character(FE[i]) == FALSE) {
                stop("Some element in FE is not a string.")
            }
        }
    }
    if (is.null(weights) == FALSE) {
        if (is.character(weights) == FALSE) {
            stop("\"weights\" is not a string.")
        } else {
            weights <- weights[1]
        }   
    }
    if (is.logical(SE) == FALSE & is.numeric(SE)==FALSE) {
        stop("\"SE\" is not a logical flag.")
    }
    if (is.null(k) == FALSE) {
        if (is.numeric(k) == FALSE) {
            stop("\"k\" is not numeric.")
        } else {
            k <- k[1]
            if (k <= 0) {
                stop("\"k\" is not a positive number.")
            }
        }
    }
    if (length(angle)>=5 | length(angle)<1) {
        stop("\"angle\" must have length 1 to 4.")
    } else {
        if (is.numeric(angle)==FALSE) {
            stop("Some element in angle is not numeric.")
        }
    }
    if (is.null(Ylabel)==TRUE) {
        Ylabel <- Y
    } else {
        if (is.character(Ylabel) == FALSE) {
            stop("\"Ylabel\" is not a string.")
        } else {
            Ylabel <- Ylabel[1]
        }   
    } 
    if (is.null(Dlabel)==TRUE) {
        Dlabel <- D   
    } else {
        if (is.character(Dlabel) == FALSE) {
            stop("\"Dlabel\" is not a string.")
        } else {
            Dlabel <- Dlabel[1]
        }   
    }
    if (is.null(Xlabel)==TRUE) {
        Xlabel <- X   
    } else {
        if (is.character(Xlabel) == FALSE) {
            stop("\"Xlabel\" is not a string.")
        } else {
            Xlabel <- Xlabel[1]
        }   
    }
    if (length(unique(data[,D]))==2) {
        warning("The treatment variable is dichotomous.")
    }

    ## drop missing values
    data <- na.omit(data[,c(Y, D, X, Z, FE)])
    
    requireNamespace("mgcv") 
    if (is.null(FE)==FALSE) {
        if (is.null(Z)==TRUE) {Z<-c()}
        for (i in 1:length(FE)) {
            Z<-c(Z,paste("as.factor(",FE[i],")",sep=""))
        } 
    }
    
    if (is.null(Z)==TRUE) { # no controls
        formula<-as.formula(paste(Y,"~","s(",D,",",X,",k=",k,")"))
    } else {
        formula<-as.formula(paste(Y,"~","s(",D,",",X,",k=",k,")+",
                                  paste(Z,collapse="+")))
    }
    if (is.null(weights)==TRUE) {
        model<-gam(formula, data=data)
    } else {
        model<-gam(formula, data=data,weights=data[,weights])
    }
    par(mfrow=c(2,2),mar=c(2,2,0,0))
    if (SE==0) {
        for (i in angle) {
            vis.gam(model, view=c(D,X), type="response",cex.lab=1.5,
                    xlab = Dlabel, ylab = Xlabel, zlab=Ylabel,
                    ticktype="detailed",plot.type="persp",n.grid=40,too.far=0.5,theta=i,phi=20)
        }
    } else {
        for (i in angle) {
            vis.gam(model, view=c(D, X), type="response",cex.lab=1.5,
                    xlab = Dlabel, ylab = Xlabel, zlab=Ylabel,
                    ticktype="detailed",plot.type="persp",se=2,n.grid=40,too.far=0.5,theta=30,phi=20)
        }
    } 
}


####################################################

inter.binning<-function(data,
                        Y, # outcome
                        D, # treatment indicator
                        X, # moderator
                        Z = NULL, # covariates
                        FE = NULL, # fixed effects
                        weights = NULL, # weighting variable
                        na.rm = FALSE,
                        Xunif = FALSE,
                        nbins = 3,  # No. of X bins
                        cutoffs = NULL,
                        vartype = "homoscedastic", # variance type
                        ##  "homoscedastic" (default); "robust"; "cluster", "pcse"
                        cl = NULL, # variable to be clustered on
                        time = NULL, # time variable for pcse
                        pairwise = TRUE, # pcse option
                        wald = TRUE,
                        figure = TRUE,
                        Xdistr = "histogram", # c("density","histogram","none")
                        main = NULL,
                        Ylabel = NULL,
                        Dlabel = NULL,
                        Xlabel = NULL,
                        xlab = NULL,
                        ylab = NULL,
                        xlim = NULL,
                        ylim = NULL,
                        theme.bw = FALSE,
                        show.grid = TRUE,
                        cex.main = NULL,
                        cex.lab = NULL,
                        cex.axis = NULL,
                        bin.labs = TRUE,                        
                        interval = NULL,
                        file = NULL                        
                        ){
    
    x <- NULL
    y <- NULL
    xmin <- NULL
    xmax <- NULL
    ymin <- NULL
    ymax <- NULL

    ## in case data is in tibble format
    if (is.data.frame(data) == FALSE) {
        data <- as.data.frame(data)
    }

    ## check input
    if (is.character(Y) == FALSE) {
        stop("\"Y\" is not a string.")
    } else {
        Y <- Y[1]
    }
    if (is.character(D) == FALSE) {
        stop("\"D\" is not a string.")
    } else {
        D <- D[1]
    }
    if (is.character(X) == FALSE) {
        stop("\"X\" is not a string.")
    } else {
        X <- X[1]
    }
    if (is.null(Z) == FALSE) {
        for (i in 1:length(Z)) {
            if (is.character(Z[i]) == FALSE) {
                stop("Some element in Z is not a string.")
            }
        }
    }
    if (is.null(FE) == FALSE) {
        requireNamespace("lfe")
        for (i in 1:length(FE)) {
            if (is.character(FE[i]) == FALSE) {
                stop("Some element in FE is not a string.")
            }
        }
        if (vartype == "pcse") {
            vartype <- "cluster"
            warning("Fixed-effect models do not allow panel corrected standard errors; changed to clustered standard errors.")
        }
    }
    if (is.null(weights) == FALSE) {
        if (is.character(weights) == FALSE) {
            stop("\"weights\" is not a string.")
        } else {
            weights <- weights[1]
        }   
    }
    if (is.logical(na.rm) == FALSE & is.numeric(na.rm)==FALSE) {
        stop("\"na.rm\" is not a logical flag.")
    }
    if (is.null(nbins) == FALSE) {
        if (nbins%%1 != 0) {
            stop("\"nbins\" is not a positive integer.")
        } else {
            nbins <- nbins[1]
        }
        if (nbins < 1) {
            stop("\"nbins\" is not a positive integer.")
        }
    }
    if (is.null(cutoffs) == FALSE) {
        if (is.numeric(cutoffs) == FALSE) {
            stop("Some element of cutoffs is not numeric.")
        } 
    }
    if (is.null(vartype)==TRUE) {
        vartype <- "homoscedastic"
    }
    if (!vartype %in% c("homoscedastic","robust","cluster","pcse")){
        stop("\"vartype\" must be one of the following: \"homoscedastic\",\"robust\",\"cluster\",\"pcse\".")
    } else if (vartype == "cluster") {
        if (is.null(cl)==TRUE) {
            stop("\"cl\" not specified; set cl = \"varname\".")
        }
    } else if (vartype == "pcse") {
        if (is.null(cl)==TRUE | is.null(time)==TRUE) {
            stop("\"cl\" or \"time\" not specified; set cl = \"varname1\", time = \"varname2\":.")
        }
    }
    if (is.null(cl)==FALSE) {
        if (is.character(cl) == FALSE) {
            stop("cl is not a string.")
        } else {
            cl <- cl[1]
            if (vartype != "pcse") {
                vartype <- "cluster"
            }
        }
    }
    if (is.null(time)==FALSE) {
        if (is.character(time) == FALSE) {
            stop("time is not a string.")
        } else {
            time <- time[1]
        }
    }
    if (is.logical(pairwise) == FALSE & is.numeric(pairwise)==FALSE) {
        stop("pairwise is not a logical flag.")
    } 
    if (is.null(Ylabel)==TRUE) {
        Ylabel <- Y
    } else {
        if (is.character(Ylabel) == FALSE) {
            stop("Ylabel is not a string.")
        } else {
            Ylabel <- Ylabel[1]
        }   
    } 
    if (is.null(Dlabel)==TRUE) {
        Dlabel <- D   
    } else {
        if (is.character(Dlabel) == FALSE) {
            stop("Dlabel is not a string.")
        } else {
            Dlabel <- Dlabel[1]
        }   
    }
    if (is.null(Xlabel)==TRUE) {
        Xlabel <- X   
    } else {
        if (is.character(Xlabel) == FALSE) {
            stop("Xlabel is not a string.")
        } else {
            Xlabel <- Xlabel[1]
        }   
    }
    if (is.null(main)==FALSE) {
        main <- main[1]
    }
    if (!Xdistr %in% c("hist","histogram","density","none")){
        stop("\"Xdistr\" must be \"histogram\", \"density\", or \"none\".")
    }
    if (is.null(xlim)==FALSE) {
        if (is.numeric(xlim)==FALSE) {
            stop("Some element in xlim is not numeric.")
        } else {
            if (length(xlim)!=2) {
                stop("xlim must be of length 2.")
            }
        }
    }
    if (is.null(ylim)==FALSE) {
        if (is.numeric(ylim)==FALSE) {
            stop("Some element in ylim is not numeric.")
        } else {
            if (length(ylim)!=2) {
                stop("ylim must be of length 2.")
            }
        }
    }

    ## axis labels
    if(is.null(xlab)==FALSE){
        if (is.character(xlab) == FALSE) {
            stop("\"xlab\" is not a string.")
        }        
    }
    if(is.null(ylab)==FALSE){
        if (is.character(ylab) == FALSE) {
            stop("\"ylab\" is not a string.")
        }        
    }
    
    ## check missing values
    vars <- c(Y, D, X, Z, FE, cl, weights)
    if (na.rm == TRUE) {        
        data <- na.omit(data[,vars])
    } else {
        if (sum(is.na(data[,vars]))>0) {
            stop("Missing values. Try option na.rm = TRUE\n")
        }
    }
    
    #################################    
    

    n<-dim(data)[1]
    data[,D]<-as.numeric(data[,D])

    ## transform moderator into uniform
    if (Xunif == TRUE) {
        x <- data[, X]
        data[, X] <- rank(x, ties.method = "average")/length(x)*100
        Xlabel <- paste(Xlabel,"(Percentile)")
    }

    ## formula
    mod.f<-paste0(Y,"~",D,"+",X,"+",D,"*",X)
    if (is.null(Z)==FALSE) {
        mod.f <- paste0(mod.f, "+", paste0(Z,collapse="+"))
    }
    if (is.null(FE)==FALSE) {
        mod.f <- paste0(mod.f, "|",paste0(FE, collapse = "+"))
        if (vartype=="cluster") {
            mod.f <- paste0(mod.f, "| 0 |",paste0(cl,collapse = "+"))
        }
    }
    mod.f <- as.formula(mod.f)
    
     
    ## a naive fit
    if (is.null(FE)==TRUE) { #OLS
        if (is.null(weights)==TRUE) {
            mod.naive<-lm(mod.f,data=data)
        } else {
            mod.naive<-lm(mod.f,data=data,weights=data[,weights])
        }
    } else { # FE
        if (is.null(weights)==TRUE) {
            mod.naive<-felm(mod.f,data=data)
        } else {
            mod.naive<-felm(mod.f,data=data,weights=data[,weights])
        }
    }
  
    ## coefficients
    coefs<-summary(mod.naive)$coefficients[,1]
    coef.D<-coefs[D]
    coef.X<-coefs[X]
    coef.DX<-coefs[paste(D,X,sep=":")] #interaction
    
    ## # variance    
    if (is.null(FE)==TRUE) { #OLS
        if (vartype=="homoscedastic") {
            v<-vcov(mod.naive)
        } else if (vartype=="robust") {
            requireNamespace("sandwich")
            v<-vcov(mod.naive,type="HC1") # White with small sample correction
        } else if (vartype=="cluster") {
            v<-vcovCluster(mod.naive,cluster = data[,cl])
        } else if (vartype=="pcse") {
            requireNamespace("pcse")
            v<-pcse(mod.naive,groupN=data[,cl],groupT=data[,time],pairwise=pairwise)$vcov
        }
    } else { # FE
        if (vartype=="homoscedastic") {
            v<-vcov(mod.naive, type = "iid")
        } else if (vartype=="robust") {
            v<-vcov(mod.naive, type="robust") 
        } else if (vartype=="cluster") {
            v<-vcov(mod.naive, type = "cluster") 
        }
    }
    ## get variance
    if (vartype=="pcse") {
        var.D<-v[D,D]
        var.DX<-v[paste(D,X,sep="."),paste(D,X,sep=".")]
        cov.DX<-v[D,paste(D,X,sep=".")]
    } else {
        var.D<-v[D,D]
        var.DX<-v[paste(D,X,sep=":"),paste(D,X,sep=":")]
        cov.DX<-v[D,paste(D,X,sep=":")]
    }
    
    ###  make a vector of the marginal effect of D on Y as X changes 
    X.lvls<-as.numeric(quantile(data[,X], probs=seq(0,1,0.01)))
    marg<-coef.D + coef.DX*X.lvls
    marg
    
    ## the variance is var(B1_D) + X^2*var(B_3) + 2*inst*cov(D, X)
    se<-sqrt(var.D +  X.lvls^2*var.DX + 2*X.lvls*cov.DX)
    df<-mod.naive$df.residual
    crit<-abs(qt(.025, df=df)) # critical values
    
    ##make 95% confidence bands. 
    lb<-marg-crit*se
    ub<-marg+crit*se
    
    ##################################################

    ## grouping by X
    if (is.null(cutoffs)==TRUE) {
        cuts.X<-quantile(data[,X],probs=seq(0,1,1/nbins))
        while (length(unique(cuts.X))!=nbins+1) {
            nbins<-nbins-1
            cuts.X<-quantile(data[,X],probs=seq(0,1,1/nbins))
        } 
    } else {
        cutoffs <-cutoffs[which(cutoffs>min(data[,X]) & cutoffs < max(data[,X]))]
        cuts.X<- sort(unique(c(min(data[,X]),cutoffs,max(data[,X]))))
        
    } 
    groupX<-cut(data[,X],breaks=cuts.X, labels = FALSE)
    groupX[which(data[,X]==min(data[,X]))]<-1
    nbins <- length(unique(groupX))
    
    ## X labels
    groupX2 <- cut(data[,X],breaks=cuts.X)
    gp.lab = paste(Xlabel, ": ", levels(groupX2), sep="")
    gp.lab[1] <- paste(Xlabel, ": [", substring(levels(groupX2)[1],2), sep = "")
    
    
    ############## Discretize X #################
    
    ## mid points
    x0<-rep(NA,nbins)
    for (i in 1:nbins) x0[i]<-median(data[which(groupX==i),X], na.rm=TRUE)
    
    ## create dummies for bins and interactions
    ## G -- a matrix of group dummies
    ## DG -- a matrix of interactions 
    
    G<-DG<-GX<-DGX<-matrix(0,n,nbins)
    for (i in 1:nbins) {
        G[which(groupX==i),i]<-1
        DG[,i]<-data[,D]*G[,i]
        GX[,i]<-G[,i]*(data[,X]-x0[i])
        DGX[,i]<-DG[,i]*(data[,X]-x0[i])
    }
    
    ## formula and esitmation
    Gs<-GXs<-DGs<-DGXs<-c()
    for (i in 1:nbins)  {
        Gs<-c(Gs,paste0("G[,",i,"]"))
        GXs<-c(GXs,paste0("GX[,",i,"]"))
        DGs<-c(DGs,paste0("DG[,",i,"]"))
        DGXs<-c(DGXs,paste0("DGX[,",i,"]"))
    }
    Xf<-paste0(Y,"~ -1+",paste0(DGs,collapse="+"),"+",paste0(DGXs,collapse="+"),
              "+",paste0(Gs,collapse="+"),"+",paste0(GXs,collapse="+"))
    if (is.null(Z)==FALSE) {
        Xf<- paste0(Xf,"+",paste0(Z,collapse="+"))
    }
    if (is.null(FE)==FALSE) {
        Xf <- paste0(Xf, "|",paste0(FE, collapse = "+"))
        if (vartype=="cluster") {
            Xf <- paste0(Xf, "| 0 |",paste0(cl,collapse = "+"))
        }
    }
    mod.Xf<-as.formula(Xf)    
     
    ## fit
    if (is.null(FE)==TRUE) { #OLS
        if (is.null(weights)==TRUE) {
            mod.X<-lm(mod.Xf,data=data)
        } else {
            mod.X<-lm(mod.Xf,data=data,weights=data[,weights])
        }
    } else { # FE
        if (is.null(weights)==TRUE) {
            mod.X<-suppressWarnings(felm(mod.Xf,data=data))
        } else {
            mod.X<-suppressWarnings(felm(mod.Xf,data=data,weights=data[,weights]))
        }
    }

    ## coefficients and CIs
    if (is.null(FE)==TRUE) { #OLS
        if (vartype=="homoscedastic") {
            X.v<-vcov(mod.X)
        } else if (vartype=="robust") {
            X.v<-vcov(mod.X,type="HC1") ## White with small sample correction
        } else if (vartype=="cluster") {
            X.v<-vcovCluster(mod.X,cluster=data[,cl])
        } else if (vartype=="pcse") {
            if (is.null(Z)==FALSE) {
                exclude<-names(which(is.na(mod.X$coefficients)==TRUE))  ## drop colinear variables
                Z.ex<-setdiff(Z,exclude)
                Xf<-paste(Y,"~ -1+",paste(DGs,collapse="+"),"+",paste(DGXs,collapse="+"),
                                "+",paste(Gs,collapse="+"),"+",paste(GXs,collapse="+"),"+",paste(Z.ex,collapse="+"),sep="")
                mod.X<-lm(as.formula(Xf),data=data)
            }
            X.v<-pcse(mod.X,groupN=data[,cl],groupT=data[,time],pairwise=pairwise)$vcov
        }
    } else { # FE
        if (vartype=="homoscedastic") {
            X.v<-vcov(mod.X, type = "iid")
        } else if (vartype=="robust") {
            X.v<-vcov(mod.X, type="robust") 
        } else if (vartype=="cluster") {
            X.v<-vcov(mod.X, type = "cluster") 
        }
    }
    Xcoefs<-mod.X$coefficients[1:nbins]
    X.v<-X.v[1:nbins,1:nbins] 
    X.se<-sqrt(diag(as.matrix(X.v,drop=FALSE)))
    X.se[which(is.na(Xcoefs))]<-NA
    df.X<-mod.X$df.residual
    crit.X<-abs(qt(.025, df=df.X))
    lb.X<-Xcoefs-crit.X*X.se
    ub.X<-Xcoefs+crit.X*X.se

    est.lin<-data.frame(X.lvls, marg, lb, ub)
    est.bin <- data.frame(x0, Xcoefs, X.se, lb.X, ub.X)
    colnames(est.bin) <- c("x0", "coef", "se", "CI_lower", "CI_upper")
    rownames(est.bin) <- gp.lab

    ################### testing  ###############################

    ## variance of treatment in each group 
    varD<-c()
    for (i in 1:nbins) {
        varD<-c(varD,var(data[groupX==i,D]/mean(data[groupX==i,D])))
    }

    ## L kurtosis
    requireNamespace("Lmoments")
    Lkurtosis <- Lmoments(data[,X],returnobject=TRUE)$ratios[4]
    
    ## if the signs are correct
    correctOrder<-ifelse(as.numeric((Xcoefs[1]-Xcoefs[2])*(Xcoefs[2]-Xcoefs[3]))>0,TRUE,FALSE) 
    
    ## p values
    pvalue<-function(i,j){
        stat<-(Xcoefs[i]-Xcoefs[j])/sqrt(X.v[i,i]+X.v[j,j]-2*X.v[i,j])
        p<-(1-pt(abs(stat),df.X))*2
        return(p)
    }
    if (nbins==3) {
        p.twosided<-round(c(pvalue(1,2),pvalue(2,3),pvalue(1,3)),digits=4)
        names(p.twosided)<-c("p.1v2","p.2v3","p.1v3")
        names(Xcoefs)<-c("X_low","X_med","X_high")
    } else if (nbins==2) {
        p.twosided<-round(pvalue(1,2),digits=4)
        names(p.twosided)<-c("p.LvH")
        names(Xcoefs)<-c("X_low","X_high")
    } else if (nbins==4) {
        names(Xcoefs)<-c("X_low","X_med1","X_med2","X_high")
    }

   

    ##############  Ward Test #####################
    
    if (wald == TRUE) { 
        ## formula
        formula0 <- paste(Y,"~",D,"+",X,"+",D,"*",X)
        ## create dummies for bins and interactions
        ## G -- a matrix of group dummies
        ## DG -- a matrix of interactions
        G<-DG<-GX<-DGX<-matrix(0,n,(nbins-1))
        for (i in 1:(nbins-1)) {
            G[which(groupX==(i+1)),i]<-1
            DG[,i]<-data[,D]*G[,i]
            GX[,i]<-data[,X]*G[,i]
            DGX[,i]<-data[,D]*data[,X]*G[,i]
        } 
        ## formula and esitmation
        Gs<-GXs<-DGs<-DGXs<-c()
        for (i in 2:nbins)  {
            Gs<-c(Gs,paste("G",i,sep=""))
            GXs<-c(GXs,paste("GX",i,sep=""))
            DGs<-c(DGs,paste("DG",i,sep=""))
            DGXs<-c(DGXs,paste("DGX",i,sep=""))
        }
        colnames(G) <- Gs
        colnames(DG) <- DGs
        colnames(GX) <- GXs
        colnames(DGX) <- DGXs
        data.aug <- cbind.data.frame(data, G, DG, GX, DGX)
        formula1<-paste(formula0,
            "+",paste(Gs,collapse=" + "),
            "+",paste(GXs,collapse=" + "),
            "+",paste(DGs,collapse=" + "),
            "+",paste(DGXs,collapse=" + "))
        if (is.null(Z)==FALSE) {
            formula0 <- paste0(formula0, "+",paste(Z,collapse=" + "))
            formula1 <- paste0(formula1, "+",paste(Z,collapse=" + "))
        } 
        if (is.null(FE)==FALSE) {
            formula1 <- paste0(formula1, "|",paste0(FE, collapse = "+"))
            if (vartype=="cluster") {
                formula1 <- paste0(formula1, "| 0 |",paste0(cl,collapse = "+"))
            }
        }
        mod.formula0<-as.formula(formula0)
        mod.formula1<-as.formula(formula1)    

        ## fit
        if (is.null(FE)==TRUE) { #OLS
            ## fit
            if (is.null(weights)==TRUE) {
                mod.re<-lm(mod.formula0,data=data.aug)
                mod.un<-lm(mod.formula1,data=data.aug)
            } else {
                mod.re<-lm(mod.formula0,data=data.aug,weights=data.aug[,weights])
                mod.un<-lm(mod.formula1,data=data.aug,weights=data.aug[,weights])
            }
            ## vcov
            if (is.null(vartype)==TRUE) {vartype <- "homoscedastic"}
            if (vartype=="homoscedastic") {
                v<-vcov(mod.un)
            } else if (vartype=="robust") {
                v<-vcov(mod.un,type="HC1") # White with small sample correction
            } else if (vartype=="cluster") {
                v<-vcovCluster(mod.un,cluster = data.aug[,cl])
            } else if (vartype=="pcse") {
                v<-pcse(mod.un,
                    groupN=data.aug[,cl],
                    groupT=data.aug[,time],
                    pairwise=pairwise)$vcov
            }
            ## wald test
            requireNamespace("lmtest")
            wald.out <- tryCatch(
                        p.wald <- round(waldtest(mod.re, mod.un,test="Chisq", vcov=v)[[4]][2],4),
                                        error = function(e){return(NULL)}
                                        )   
            ## warning
            if (is.null(wald.out)==TRUE) {
                p.wald <- NULL
                warning("Var-cov matrix nearly singular in the Wald test.")
            }  
        } else { # FE
            requireNamespace("lfe")
            ## fit
            if (is.null(weights)==TRUE) {
                mod.un<-suppressWarnings(felm(mod.formula1,data=data.aug))
            } else {
                mod.un<-suppressWarnings(felm(mod.formula1,data=data.aug,weights=data.aug[,weights]))
            }
            ## wald test
            constraints <- as.formula(paste0("~",paste0(c(Gs,GXs,DGs,DGXs), collapse = "|")))            
            if (vartype=="homoscedastic") {
                p.wald <- lfe::waldtest(mod.un, constraints, type = "default")[1]
            } else if (vartype=="robust") {
                p.wald <- lfe::waldtest(mod.un, constraints, type = "robust")[1]
            } else {
                p.wald <- lfe::waldtest(mod.un, constraints)[1] # clustered
            }
            names(p.wald) <- NULL            
            p.wald <- round(p.wald,4)        
        }

    } # end of Wald test
     
    ##################################
    ## storage
    ##################################

    ## binary treatment (for plotting)
    binaryD <- FALSE
    if (length(unique(data[,D]))==2) {binaryD <- TRUE}

    ## get densities and histogram
    if (binaryD == TRUE) { ## binary D
        # density
        de <- NULL
        if (is.null(weights)==TRUE) {
            de.co <- density(data[data[,D]==0,X])
            de.tr <- density(data[data[,D]==1,X])
        } else {
            suppressWarnings(de.co<-density(data[data[,D]==0,X],
                weights=data[data[,D]==0,weights]))
            suppressWarnings(de.tr<-density(data[data[,D]==1,X],
                weights=data[data[,D]==1,weights]))
        }
        # histogram
        if (is.null(weights)==TRUE) {
            hist.out<-hist(data[,X],breaks=80,plot=FALSE)
        } else {
            suppressWarnings(hist.out<-hist(data[,X],data[,weights],
                breaks=80,plot=FALSE))
        } 
        n.hist<-length(hist.out$mids)
        ## count the number of treated
        count1<-rep(0,n.hist)
        treat<-which(data[,D]==max(data[,D]))
        for (i in 1:n.hist) {
            count1[i]<-sum(data[treat,X]>=hist.out$breaks[i] &
            data[treat,X]<hist.out$breaks[(i+1)])
        }
        count1[n.hist]<-sum(data[treat,X]>=hist.out$breaks[n.hist] &
            data[treat,X]<=hist.out$breaks[n.hist+1])        
    }  else { ## continuous D
        if (is.null(weights)==TRUE) {
            de <- density(data[,X])
        } else {
            suppressWarnings(de <- density(data[,X],weights=data[,weights]))
        }
        if (is.null(weights)==TRUE) {
            hist.out<-hist(data[,X],breaks=80,plot=FALSE)
        } else {
            suppressWarnings(hist.out<-hist(data[,X],data[,weights],
                breaks=80,plot=FALSE))
        }
        de.co <- de.tr <- NULL 
        count1 <- NULL
    }  

    tests <- list(
        est.binning = round(est.bin, 3),
        bin.size = sprintf(c(table(groupX)/length(groupX)),fmt = '%#.3f'),
        treat.variation.byGroup=sprintf(varD,fmt = '%#.3f'),
        X.Lkurtosis = sprintf(Lkurtosis,fmt = '%#.3f')
        )
    if (nbins==3) {
        tests<-c(tests,list(correctOrder=correctOrder))
    }
    if (nbins%in%c(2,3)) {
        tests<-c(tests,list(p.twosided=sprintf(p.twosided, fmt = '%#.3f'))) 
    }
    if (wald == TRUE) {
        tests <- c(tests, list(p.wald = sprintf(p.wald, fmt = '%#.3f')))
    }    
    
    ## saving
    output<-list(
        type = "binning",
        nbins = nbins,
        est.lin = est.lin,
        est.bin = est.bin,
        binaryD= binaryD,
        Xlabel = Xlabel,
        Dlabel = Dlabel,
        Ylabel = Ylabel,
        de = de, # density
        de.co = de.co, 
        de.tr = de.tr, 
        hist.out = hist.out,
        count.tr = count1,
        tests = tests
        )
    
    
    if (figure==TRUE) {

        class(output) <- "interflex"
        graph <- inter.plot(out = output, CI = TRUE, xlab = xlab, ylab = ylab,
           Ylabel = Ylabel, Dlabel = Dlabel, Xlabel = Xlabel, 
           main = main, xlim = xlim, ylim = ylim, Xdistr = Xdistr,
           file = file, theme.bw = theme.bw, show.grid = show.grid, bin.labs = bin.labs, 
           cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab)  

        output<-c(output,list(graph=graph))
    }
    
    class(output) <- "interflex"
    return(output)
}

###########################

inter.kernel <- function(data,
   Y, D, X,
   Z = NULL,
   weights = NULL, # weighting variable
   FE = NULL,
   na.rm = FALSE,
   Xunif = FALSE, ## transform moderate into a uniform distribution
   CI = TRUE,
   conf.level = 0.95,
   cl = NULL,
   neval = 50,
   nboots = 200,
   parallel = TRUE,
   cores = 4,
   seed = 02139,
   bw = NULL, # bandwidth   
   grid = 20, # either a number of a sequence of numbers
   metric = "MSPE",
   Xdistr = "histogram",
   main = NULL,
   Ylabel = NULL,
   Dlabel = NULL,
   Xlabel = NULL,
   xlab = NULL,
   ylab = NULL,
   xlim = NULL,
   ylim = NULL,                         
   theme.bw = FALSE,
   show.grid = TRUE,
   cex.main = NULL,
   cex.lab = NULL,
   cex.axis = NULL,
   file = NULL                        
   ){

    x <- NULL
    y <- NULL
    xmin <- NULL
    xmax <- NULL
    ymin <- NULL
    ymax <- NULL
    ME <- NULL
    CI_lower <- NULL
    CI_upper <- NULL

    ## in case data is in tibble format
    if (is.data.frame(data) == FALSE) {
        data <- as.data.frame(data)
    }

    ## check input
    if (is.character(Y) == FALSE) {
        stop("\"Y\" is not a string.")
    } else {
        Y <- Y[1]
    }
    if (is.character(D) == FALSE) {
        stop("\"D\" is not a string.")
    } else {
        D <- D[1]
    }
    if (is.character(X) == FALSE) {
        stop("\"X\" is not a string.")
    } else {
        X <- X[1]
    }
    if (is.null(Z) == FALSE) {
        for (i in 1:length(Z)) {
            if (is.character(Z[i]) == FALSE) {
                stop("Some element in Z is not a string.")
            }
        }
    }
    if (is.null(FE) == FALSE) {
        for (i in 1:length(FE)) {
            if (is.character(FE[i]) == FALSE) {
                stop("Some element in FE is not a string.")
            }
        }
    }
    if (is.null(weights) == FALSE) {
        if (is.character(weights) == FALSE) {
            stop("\"weights\" is not a string.")
        } else {
            weights <- weights[1]
        }   
    }
    if (is.logical(na.rm) == FALSE & is.numeric(na.rm)==FALSE) {
        stop("\"na.rm\" is not a logical flag.")
    } 
    if (is.logical(CI) == FALSE & is.numeric(CI)==FALSE) {
        stop("\"CI\" is not a logical flag.")
    }
    if (is.null(conf.level)==FALSE) {
        if (is.numeric(conf.level)==FALSE) {
            stop("\"conf.level\" should be a number between 0.5 and 1.")
        } else {
            if (conf.level<=0.5 | conf.level>1) {
                stop("\"conf.level\" should be a number between 0.5 and 1.")
            }
        } 
    }
    if (is.null(cl)==FALSE) {
        if (is.character(cl) == FALSE) {
            stop("\"cl\" is not a string.")
        } else {
            cl <- cl[1] 
        }
    }
    if (is.null(neval)==FALSE) {
        if (is.numeric(neval)==FALSE) {
            stop("\"neval\" is not a positive integer.")
        } else {
            neval <- neval[1]
            if (neval%%1!=0 | neval<=0) {
                stop("\"neval\" is not a positive integer.")
            }  
        } 
    } 
    if (is.null(nboots) == FALSE) {
        if (is.numeric(nboots)==FALSE) {
            stop("\"nboots\" is not a positive integer.")
        } else {
            nboots <- nboots[1]
            if (nboots%%1 != 0 | nboots < 1) {
                stop("\"nboots\" is not a positive number.")
            }
        } 
    }
    if (is.logical(parallel) == FALSE & is.numeric(parallel)==FALSE) {
        stop("\"paralell\" is not a logical flag.")
    }
    if (is.numeric(cores)==FALSE) {
        stop("\"cores\" is not a positive integer.")
    } else {
        cores <- cores[1]
        if (cores%%1!= 0 | cores<=0) {
            stop("\"cores\" is not a positive integer.")
        }
    }
    if (is.null(bw)==FALSE) {
        if (is.numeric(bw)==FALSE) {
            stop("\"bw\" should be a positive number.")
        } else {
            bw <- bw[1]
        } 
        if (bw<=0) {
            stop("\"bw\" should be a positive number.")
        }
    }
    if (is.numeric(seed)==FALSE) {
        stop("\"seed\" should be a number.")
    }
    if (is.numeric(grid) == FALSE) {
        stop("\"grid\" should be numeric.")
    } else {
        if (length(grid)==1) {
            if (grid%%1 != 0 | grid<1) {
                stop("\"grid\" is not a positive integer.")
            }
        } else {
            grid <- grid[which(grid>0)]
        }
    }
    if (!metric%in%c("MSPE","MAPE")) {
        stop("\"metric\" should be either \"MSPE\" or \"MAPE\".")
    } 
    if (is.null(Ylabel)==TRUE) {
        Ylabel <- Y
    } else {
        if (is.character(Ylabel) == FALSE) {
            stop("\"Ylabel\" is not a string.")
        } else {
            Ylabel <- Ylabel[1]
        }   
    } 
    if (is.null(Dlabel)==TRUE) {
        Dlabel <- D   
    } else {
        if (is.character(Dlabel) == FALSE) {
            stop("\"Dlabel\" is not a string.")
        } else {
            Dlabel <- Dlabel[1]
        }   
    }
    if (is.null(Xlabel)==TRUE) {
        Xlabel <- X   
    } else {
        if (is.character(Xlabel) == FALSE) {
            stop("\"Xlabel\" is not a string.")
        } else {
            Xlabel <- Xlabel[1]
        }   
    }
    if (is.null(main)==FALSE) {
        main <- main[1]
    } 
    if (!Xdistr %in% c("hist","histogram","density","none")){
        stop("\"Xdistr\" must be \"histogram\", \"density\", or \"none\".")
    }
    if (is.null(xlim)==FALSE) {
        if (is.numeric(xlim)==FALSE) {
            stop("Some element in xlim is not numeric.")
        } else {
            if (length(xlim)!=2) {
                stop("\"xlim\" must be of length 2.")
            }
        }
    }
    if (is.null(ylim)==FALSE) {
        if (is.numeric(ylim)==FALSE) {
            stop("Some element in ylim is not numeric.")
        } else {
            if (length(ylim)!=2) {
                stop("\"ylim\" must be of length 2.")
            }
        }
    }

    ## axis labels
    if(is.null(xlab)==FALSE){
        if (is.character(xlab) == FALSE) {
            stop("\"xlab\" is not a string.")
        }        
    }
    if(is.null(ylab)==FALSE){
        if (is.character(ylab) == FALSE) {
            stop("\"ylab\" is not a string.")
        }        
    }
    
    
    ## set seed 
    if (is.null(seed)==FALSE) {
        set.seed(seed)
    }

    ## check variable name
    M <- c(Y,D,X,Z,FE,cl,weights)
    for (var in M){
        if((var %in% names(data))==FALSE )
            stop("Wrong variable name.")
    }
    
    ## check missing values
    vars <- unique(c(Y, D, X, Z, FE, cl, weights))
    if (na.rm == TRUE) {        
        data <- na.omit(data[,vars])
    } else {
        if (sum(is.na(data[,vars]))>0) {
            stop("Missing values. Try option na.rm = TRUE\n")
        }
    }
    n<-dim(data)[1]

    ## transform moderator into uniform
    if (Xunif == TRUE) {
        x <- data[, X]
        data[, X] <- rank(x, ties.method = "average")/length(x)*100
        Xlabel <- paste(Xlabel,"(Percentile)")
    }

    ## check cluster
    if (is.null(cl)==TRUE & is.null(FE)==FALSE) {
        warnings("Fixed effects model assumed. Clustering standard errors highly recommended.")
    }

    ## change fixed effect variable to factor
    if (is.null(FE)==FALSE) {
        if (length(FE) == 1) {
            data[, FE] <- as.numeric(as.factor(data[, FE]))
        } else {
            data[, FE] <- sapply(data[,FE],function(vec){as.numeric(as.factor(vec))})
        }
    }
    
    ## labels
    if (is.null(Xlabel)==TRUE) {Xlabel = X}
    if (is.null(Ylabel)==TRUE) {Ylabel = Y}
    if (is.null(Dlabel)==TRUE) {Dlabel = D}

    ## X more than 5 values
    if (length(unique(data[ ,X]))< 5) {
        warning("Moderator has less than 5 values; consider a fully saturated model.")
    } 

     
    ## preparation for parallel computing
    if (parallel == TRUE & (CI == TRUE|is.null(bw))) {
        requireNamespace("doParallel")
        ## require(iterators)
        maxcores <- detectCores()
        cores <- min(maxcores, cores)
        pcl<-makeCluster(cores)  
        doParallel::registerDoParallel(pcl)
        cat("Parallel computing with", cores,"cores...\n") 
    }
    
    ## evaluation points
    X.eval <- seq(min(data[,X]), max(data[,X]), length.out = neval)
    
    ## bandwidth selection
    if(is.null(bw) == TRUE){
        CV <- 1
        if (length(grid) == 1) {
            rangeX <- max(data[, X]) - min(data[, X])
            grid <- exp(seq(log(rangeX/50), log(rangeX), length.out = grid))
        }
        cv.out <- crossvalidate(data = data, X.eval = X.eval,
                                Y = Y, D = D, X = X, Z = Z,
                                FE = FE, cl = cl,
                                weights = weights,
                                grid = grid, metric = metric,
                                parallel = parallel)
        bw <- cv.out$opt.bw
    }  else {
        CV <- 0
    }

    ## Estimates
    if (CI == FALSE) {
        est <- coefs(data = data, bw = bw, Y = Y, X = X, D = D, Z = Z,
                     FE = FE, X.eval = X.eval, weights = weights)[,c(1,3)] 
    } else { ## bootstrap
        ## point estimate
        coef <- coefs(data = data, bw = bw, Y = Y, X = X, D = D, Z = Z,
                      FE = FE, X.eval = X.eval, weights = weights)[,3]
        ## boostrap SEs
        if (is.null(cl)==FALSE) { ## find clusters
            clusters<-unique(data[,cl])
            id.list<-split(1:n,data[,cl])
        } 
        oneboot <- function() {
            if (is.null(cl)==TRUE) {
                smp<-sample(1:n,n,replace=TRUE)
            } else { ## block bootstrap
                cluster.boot<-sample(clusters,length(clusters),replace=TRUE)
                smp<-unlist(id.list[match(cluster.boot,clusters)])
            }   
            s<-data[smp,]
            out <- coefs(data = s, bw = bw, Y = Y, X = X, D = D, Z = Z,
                         FE = FE, X.eval = X.eval, weights = weights)[,3]
            return(out)
        } 
        cat("Bootstrapping ...")
        if (parallel==TRUE) {
            suppressWarnings(
                bootout <- foreach (i=1:nboots, .combine=cbind,
                                    .export=c("oneboot","coefs"),
                                    .inorder=FALSE) %dopar% {oneboot()}
            ) 
            cat("\r")
        } else {
            bootout<-matrix(NA,length(X.eval),nboots)
            for (i in 1:nboots) { 
                bootout[,i] <- oneboot()
                if (i%%50==0) cat(i) else cat(".")
            }
            cat("\r")        
        } 
        ## summary
        CI.lvl <- c((1-conf.level)/2, (1-(1-conf.level)/2))
        ci<-t(apply(bootout, 1, quantile, CI.lvl))
        est<-data.frame(cbind("X" = X.eval, "ME" = coef,
                              "SE"=apply(bootout,1,sd),
                              "CI_lower"=ci[,1], "CI_upper"=ci[,2]))
        
    }

    if (parallel == TRUE & (CI == TRUE|CV==1)) {
        suppressWarnings(stopCluster(pcl))
        cat("\n") 
    }

    ## binary treatment
    binaryD <- FALSE
    if (length(unique(data[,D]))==2) {binaryD <- TRUE}

    ## get densities and histogram
    if (binaryD == TRUE) { ## binary D
        # density
        de <- NULL
        de.co <- density(data[data[,D]==0,X])
        de.tr <- density(data[data[,D]==1,X])
        # histogram
        hist.out<-hist(data[,X],breaks=80,plot=FALSE)
        n.hist<-length(hist.out$mids)        
        ## count the number of treated
        count1<-rep(0,n.hist)
        treat<-which(data[,D]==max(data[,D]))
        for (i in 1:n.hist) {
            count1[i]<-sum(data[treat,X]>=hist.out$breaks[i] &
            data[treat,X]<hist.out$breaks[(i+1)])
        }
        count1[n.hist]<-sum(data[treat,X]>=hist.out$breaks[n.hist] &
            data[treat,X]<=hist.out$breaks[n.hist+1])        
    }  else { ## continuous D
        de <- density(data[,X]) 
        de.co <- de.tr <- NULL
        hist.out<-hist(data[,X],breaks=80,plot=FALSE)
        count1 <- NULL
    }  

    ######################################################
    
    output<-list(
        type = "kernel",
        bw = bw,
        est = est,                 
        binaryD = binaryD, # binary treatment
        Xlabel = Xlabel,
        Dlabel = Dlabel,
        Ylabel = Ylabel,
        de = de, # density
        de.co = de.co, 
        de.tr = de.tr, 
        hist.out = hist.out,
        count.tr = count1)

    if (CV == 1) {
        output <- c(output, list(CV.out = cv.out$CV.out))
    }    
    
    #################
    ## plotting
    #################
    
    class(output) <- "interflex"
    graph <- inter.plot(out = output, CI = CI, xlab = xlab, ylab = ylab,
     Ylabel = Ylabel, Dlabel = Dlabel, Xlabel = Xlabel, 
     main = main, xlim = xlim, ylim = ylim, Xdistr = Xdistr,
     file = file, theme.bw = theme.bw, show.grid = show.grid,
     cex.main = cex.main, cex.axis = cex.axis, cex.lab = cex.lab)

    output <- c(output, list(graph = graph))
    
    class(output) <- "interflex"
    return (output)
}

#### Plot kernel results #######

inter.plot <- function(
    out,
    CI = TRUE,
    Xdistr = "histogram",
    main = NULL,
    Ylabel = NULL,
    Dlabel = NULL,
    Xlabel = NULL,
    xlab = NULL,
    ylab = NULL,
    xlim = NULL,
    ylim = NULL,
    theme.bw = FALSE,
    show.grid = TRUE,     
    cex.main = NULL,
    cex.lab = NULL,
    cex.axis = NULL,
    bin.labs = TRUE, # bin labels    
    interval = NULL, # interval in replicated papers
    file = NULL
   ){

    x <- NULL
    y <- NULL
    xmin <- NULL
    xmax <- NULL
    ymin <- NULL
    ymax <- NULL
    count1 <- NULL
    X <- NULL
    ME <- NULL
    CI_lower <- NULL
    CI_upper <- NULL
    marg <- NULL
    lb <- NULL
    ub <- NULL
    x0 <- NULL
    coef <- NULL

    if (!class(out) %in% c("interflex")) {
        stop("Not an \"interflex\" object.")
    }     
    type <- out$type
    

    requireNamespace("ggplot2")

    binaryD <- out$binaryD
    Xlabel <- out$Xlabel
    Dlabel <- out$Dlabel
    Ylabel <- out$Ylabel
    de <- out$de
    de.co <- out$de.co
    de.tr <- out$de.tr
    hist.out <- out$hist.out
    count.tr <- out$count.tr

    if (is.null(xlim)==FALSE) {
        if (is.numeric(xlim)==FALSE) {
            stop("Some element in xlim is not numeric.")
        } else {
            if (length(xlim)!=2) {
                stop("\"xlim\" must be of length 2.")
            }
        }
    }
    if (is.null(ylim)==FALSE) {
        if (is.numeric(ylim)==FALSE) {
            stop("Some element in ylim is not numeric.")
        } else {
            if (length(ylim)!=2) {
                stop("\"ylim\" must be of length 2.")
            }
        }
    }

    ## font size
    if (is.null(cex.lab)==FALSE) {
        if (is.numeric(cex.lab)==FALSE) {
            stop("\"cex.lab\" is not numeric.")
        }
    }
    if (is.null(cex.main)==FALSE) {
        if (is.numeric(cex.main)==FALSE) {
            stop("\"cex.main\" is not numeric.")
        }
    }
    if (is.null(cex.axis)==FALSE) {
        if (is.numeric(cex.axis)==FALSE) {
            stop("\"cex.axis\" is not numeric.")
        }    
    }   


    ## axis labels
    if(is.null(xlab)==TRUE){
        xlab<-c(paste("Moderator: ", Xlabel, sep=""))
    } else {
        if (is.character(xlab) == FALSE) {
            stop("\"xlab\" is not a string.")
        }        
    }
    if(is.null(ylab)==TRUE){
        ylab<-c(paste("Marginal Effect of ",Dlabel," on ",Ylabel,sep=""))
    } else {
        if (is.character(ylab) == FALSE) {
            stop("\"ylab\" is not a string.")
        }        
    }
    if (is.null(main)==FALSE) {
        main <- main[1]
    }
    if (!Xdistr %in% c("hist","histogram","density","none")){
        stop("\"Xdistr\" must be \"histogram\", \"density\", or \"none\".")
    }
    if (is.null(Xdistr) == TRUE) {
        Xdistr <- "density"
    } else if (!Xdistr %in% c("density","histogram","hist","none")) {
        Xdistr <- "density"
    }

    #########################################

    if (type == "binning") {

        # est.lin<-data.frame(X.lvls,marg,lb,ub)
        # est.bin<-data.frame(x0,Xcoefs,se,lb.X,ub.X)
        nbins <- out$nbins
        est.lin <- out$est.lin
        est.bin<-out$est.bin
        est.bin2<-est.bin[which(is.na(est.bin[,2])==FALSE),] ## non missing part
        est.bin3<-est.bin[which(is.na(est.bin[,2])==TRUE),]  ## missing part
        X.lvls <- est.lin[,1]
        errorbar.width<-(max(X.lvls)-min(X.lvls))/20        
        
        # y-range
        yrange<-na.omit(unlist(c(est.lin[,c(3,4)],est.bin[,c(4,5)])))
        if (is.null(ylim)==FALSE) {yrange<-c(ylim[2],ylim[1]+(ylim[2]-ylim[1])*1/8)}
        maxdiff<-(max(yrange)-min(yrange))
        pos<-max(yrange)-maxdiff/20

    } else if (type == "kernel") {

        est <- out$est 
        # Checking options
        if (is.logical(CI) == FALSE & is.numeric(CI)==FALSE) {
            stop("\"CI\" is not a logical flag.")
        }
        ## y-range
        if (CI == TRUE) {
            yrange<-na.omit(c(est$CI_lower, est$CI_upper))
        }  else {
            yrange<-na.omit(c(est$ME))
        }
        if (is.null(ylim)==FALSE) {yrange<-c(ylim[2],ylim[1]+(ylim[2]-ylim[1])*1/8)}
        maxdiff<-(max(yrange)-min(yrange)) 

    }

    ## Start Plotting ###
    p1 <- ggplot()

    ## black white theme and mark zero
    if (theme.bw == FALSE) {
        p1 <- p1 + geom_hline(yintercept=0,colour="white",size=2)
    } else {
        p1 <- p1 + theme_bw() + geom_hline(yintercept=0,colour="#AAAAAA50",size=2)
    }

    if (show.grid == FALSE) {
        p1 <- p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }

    ## plotting moderator distribution      
    if (Xdistr == "density") { # density plot
        if (binaryD == TRUE) { ## binary D
            ## put in data frames
            deX.ymin <- min(yrange)-maxdiff/5
            deX.co <- data.frame(x = de.co$x,
                y = de.co$y/max(de.co$y) * maxdiff/5 + min(yrange) - maxdiff/5)
            deX.tr <- data.frame(x = de.tr$x,
                y = de.tr$y/max(de.tr$y) * maxdiff/5 + min(yrange) - maxdiff/5)
            ## color
            feed.col<-col2rgb("gray50")
            col.co<-rgb(feed.col[1]/1000, feed.col[2]/1000,feed.col[3]/1000)
            col.tr<-rgb(red=1, blue=0, green=0)

            ## plotting
            p1 <- p1 + geom_ribbon(data = deX.co, aes(x = x, ymax = y, ymin = deX.ymin),
             fill = col.co, alpha = 0.2) +
            geom_ribbon(data = deX.tr, aes(x = x, ymax = y, ymin = deX.ymin),
                fill = col.tr, alpha = 0.2)          

        } else { ## continuous D
            # put in data frames
            deX.ymin <- min(yrange)-maxdiff/5
            deX <- data.frame(x = de$x,
                y = de$y/max(de$y) * maxdiff/5 + min(yrange) - maxdiff/5) 

            ## color
            feed.col<-col2rgb("gray50")
            col<-rgb(feed.col[1]/1000, feed.col[2]/1000,feed.col[3]/1000)

            ## plotting
            p1 <- p1 + geom_ribbon(data = deX, aes(x = x, ymax = y, ymin = deX.ymin),
                fill = col, alpha = 0.2)
        }
    } else if (Xdistr %in% c("histogram","hist")) { # histogram plot
        if (binaryD == TRUE) { ## binary D
            n.hist<-length(hist.out$mids)
            dist<-hist.out$mids[2]-hist.out$mids[1]
            hist.max<-max(hist.out$counts)
            ## put in a data frame
            histX<-data.frame(ymin=rep(min(yrange)-maxdiff/5,n.hist),
              ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
              xmin=hist.out$mids-dist/2,
              xmax=hist.out$mids+dist/2,
              count1=count.tr/hist.max*maxdiff/5+min(yrange)-maxdiff/5)     
            p1 <- p1 + geom_rect(data=histX,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                colour="gray50",alpha=0,size=0.5) + # control
            geom_rect(data=histX,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=count1),
              fill="red",colour="grey50",alpha=0.3,size=0.5) # treated

        }  else { ## continuous D

            n.hist<-length(hist.out$mids)
            dist<-hist.out$mids[2]-hist.out$mids[1]
            hist.max<-max(hist.out$counts)            
            histX<-data.frame(ymin=rep(min(yrange)-maxdiff/5,n.hist),
              ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
              xmin=hist.out$mids-dist/2,
              xmax=hist.out$mids+dist/2)

            p1 <- p1 + geom_rect(data=histX,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
               colour="gray50",alpha=0,size=0.5)
        }
    } # end of plotting X distribution  


    # kernel estimates
    if (type == "kernel") {  
        ## point estimates
        p1 <-  p1 + geom_line(data=est,aes(X,ME))
        ## confidence interval
        if (CI == TRUE) {
            p1 <- p1 + geom_ribbon(data=est, aes(x=X,ymin=CI_lower,ymax=CI_upper),alpha=0.2)
        }         
    } # end of kernel-specific part

    ## binning estimates and bin labels
    if (type == "binning") {    
        ## linear plot 
        p1<-p1 + geom_line(data=est.lin,aes(X.lvls,marg))+
            geom_ribbon(data=est.lin, aes(x=X.lvls,ymin=lb,ymax=ub),alpha=0.2)
        ## bin estimates
        p1<-p1+ geom_errorbar(data=est.bin2, aes(x=x0, ymin=CI_lower, ymax=CI_upper),colour="red",
                              width= errorbar.width)+
                                  geom_point(data=est.bin2,aes(x0,coef),size=4,shape=21,fill="white",colour="red") 

        ## in case there's non-overlap
        p1<-p1+annotate(geom="text", x=est.bin3[,1], y=rep(0,dim(est.bin3)[1]),
                        label="NaN",colour="red")  
        ## labels: L, M, H and so on 
        if (bin.labs == TRUE) {
            if (nbins==3) {
                p1<-p1 + annotate(geom="text", x=est.bin[1,1], y=pos,
                  label="L",colour="gray50",size=10) +
                annotate(geom="text", x=est.bin[2,1], y=pos,
                   label="M",colour="gray50",size=10) +
                annotate(geom="text", x=est.bin[3,1], y=pos,
                   label="H",colour="gray50",size=10)
            } else if (nbins==4) {
                p1<-p1 + annotate(geom="text", x=est.bin[1,1], y=pos,
                  label="L",colour="gray50",size=10) +
                annotate(geom="text", x=est.bin[2,1], y=pos,
                   label="M1",colour="gray50",size=10) +
                annotate(geom="text", x=est.bin[3,1], y=pos,
                   label="M2",colour="gray50",size=10) +
                annotate(geom="text", x=est.bin[4,1], y=pos,
                   label="H",colour="gray50",size=10)
            } 
        }        
    } # end of binning-specific case

    ## mark the original interval (in replicated papers)
    if (is.null(interval)==FALSE) {
        p1<- p1 + geom_vline(xintercept=interval,colour="steelblue", linetype=2,size=1.5)
    }

    ## Other universal options
    ## axis labels
    if (is.null(cex.lab)==TRUE) {
        cex.lab <- 15
    } else {
        cex.lab <- 15 * cex.lab
    }
    if (is.null(cex.axis)==TRUE) {
        cex.axis <- 15
    } else {
        cex.axis <- 15 * cex.axis
    }
    p1 <- p1 + xlab(xlab) + ylab(ylab) + 
        theme(axis.text = element_text(size=cex.axis), axis.title = element_text(size=cex.lab))

    ## title
    if (is.null(cex.main)==TRUE) {
        cex.main <- 18
    } else {
        cex.main <- 18 * cex.main
    }
    if (is.null(main)==FALSE) {
            p1<-p1 + ggtitle(main) +
            theme(plot.title = element_text(hjust = 0.5, size=cex.main,
                lineheight=.8, face="bold"))
    } 

    ## xlim and ylim
    if (is.null(ylim)==FALSE) {
        ylim2 = c(ylim[1]-(ylim[2]-ylim[1])*0.25/6, ylim[2]+(ylim[2]-ylim[1])*0.4/6)
    }
    if (is.null(xlim)==FALSE & is.null(ylim)==FALSE) {
        p1<-p1+coord_cartesian(xlim = xlim, ylim = ylim2)
    }
    if (is.null(xlim)==TRUE & is.null(ylim)==FALSE) {
        p1<-p1+coord_cartesian(ylim = ylim2)
    }
    if (is.null(xlim)==FALSE & is.null(ylim)==TRUE) {
        p1<-p1+coord_cartesian(xlim = xlim)
    } 
    
    ## save to file
    if (is.null(file)==FALSE) {
        ggsave(file, plot = p1)         
    } 
    ######################################################
    
    graph <- p1     
    return (graph)
}



######################################
## Weighted Least-Squares
######################################

coefs <- function(data,bw,Y,X,D,
                  Z=NULL,
                  FE=NULL,
                  weights,
                  X.eval = NULL,
                  neval = 50){

    
    
    ## evaluation points
    if (is.null(X.eval) == TRUE) {
        X.eval <- seq(min(data$x), max(data$x), length.out = neval)
    }
    ## conditions
    noZ <- is.null(Z)
    noFE <- is.null(FE)
    
    ## survey weights
    n <- dim(data)[1]
    if (is.null(weights)==TRUE) {
        weights <- rep(1, n)
    } else {
        weights <- data[,weights]
    }
    
    ## storage
    result <- matrix(NA, length(X.eval), (5 + length(Z)))
    colnames(result)<-c("X","Intercept","ME","x","Dx",Z)
    
    ## main algorithm
    if (noFE) { # no fixed effects, wls
        ## formula
        if (noZ) { ## no covariates
            formula <- as.formula("y ~ d + xx + d_xx")
        } else { ## with covariates
            formula<- as.formula(paste("y ~ d + xx + d_xx + ",paste(Z,collapse="+")))               
        } 
        ## weighted least squares without fixed effect
        wls<-function(x, dat, weights){
            dat1<-dat
            dat1$xx<-dat1$x-x
            dat1$d_xx <- dat1$d * dat1$xx
            dat1$w <- dnorm(dat1$xx/bw) * weights 
            reg <- lm(formula, data=dat1, weights = w)
            result <- c(x, reg$coef)
            result[which(is.na(result))] <- 0 
            return(result)   
        } ## end of WLS function
        
        ## data to be used
        dat<-data[, c(Y, D, X, Z)]
        colnames(dat)[1:3] <- c("y","d","x")
        ## estimation (a loop)
        for (i in 1:length(X.eval)) {
            result[i,] <- wls(X.eval[i], dat, weights = weights)
        }
        result <- data.frame(result)    
    } else{ # with fixed effects
        ## data to be used
        dat <- data[, c(Y, D, X, Z)] 
        dat.FE <- as.matrix(data[,FE],drop = FALSE)
        for (i in 1:length(X.eval)) {
            xx <- dat[,X]-X.eval[i]
            dat1 <- as.matrix(cbind(dat[,Y],dat[,D],
                                    xx, dat[,D] * xx, dat[,Z]))
            w<-dnorm(xx/bw) * weights
            estcoef <- fastplm(data = dat1, FE = dat.FE,
                               weight = w, FEcoefs = 0)$coefficients
            estcoef[which(is.nan(estcoef))] <- 0
            result[i,] <- c(X.eval[i],0,estcoef)
        }
        result <- data.frame(result) 
    } 
    return(result)
}


######################################
## Cross-validation
######################################

crossvalidate <- function(data, Y, D, X, Z = NULL, weights = NULL,
                          FE = NULL, cl = NULL, 
                          X.eval, grid, nfold = 5, parallel = FALSE, bw,
                          metric = "MSPE"){

   
    ## %dopar% <- NULL
    
    
    ## calculate error for testing set
    getError <- function(train, test, bw, Y, X, D, Z, FE, weights, X.eval){
 
        coef<-coefs(data=train,bw=bw,Y=Y,X=X,D=D,Z=Z,FE=FE,
                    weights = weights, X.eval= X.eval)
        coef[is.na(coef)] <- 0
        n2<-length(Z)
        esCoef<-function(x){
            Xnew<-abs(X.eval-x)
            d1<-min(Xnew)     ## distance between x[i] and the first nearest x in training set
            label1<-which.min(Xnew)
            Xnew[label1]<-Inf
            d2<-min(Xnew)     ## distance between x[i] and the second nearest x in training set
            label2<-which.min(Xnew)
            if(d1==0){
                func <- coef[label1,-c(1,4,5)] # X.eval (1), intercept (2), d (3), xx (4), d:xx (5), z
            }  else if(d2==0){
                func <- coef[label2,-c(1,4,5)]  
            } else{ ## weighted average 
                func <- (coef[label1,-c(1,4,5)] * d2 + coef[label2,-c(1,4,5)] * d1)/(d1 + d2) 
            }
            return(func)
        } 
        Knn<-t(sapply(test[,X],esCoef)) ## coefficients for test  class==matrix
         
        ## predicting
        test.Y <- test[,Y]
        test.X <- as.matrix(cbind(rep(1,dim(test)[1]),test[,c(D,Z)])) 
        sumOfEst<-matrix(lapply(1:length(test.X), function(i){test.X[i]*Knn[[i]]}),
                         nrow=nrow(test.X), ncol=ncol(test.X))
        error<-test.Y - rowSums(matrix(unlist(sumOfEst),length(test.Y)))  
        return(c(mean(abs(error)),mean(error^2)))
    } 
    ## grid search and 5 fold cross validation
    cv<-function(bw){
        mse<-matrix(NA,nfold,2)
        for(j in 1:nfold){ # 5-fold CV
            testid <- which(fold==j)
            mse[j,] <- getError(train= data[-testid,], test = data[testid, ],
                          bw = bw, Y=Y, X=X, D=D, Z=Z, FE=FE, weights = weights, X.eval= X.eval)
        }
        return(c(bw, apply(mse,2,mean)))
    }

    
    cat("Cross-validating bandwidth ... ")
    ## generate 5 random folds
    ## if is.null(cluster)==false then fold by cluster
    n<-dim(data)[1]
    fold<-rep(0,n)
    
    if(is.null(cl)==TRUE){        
        fold <- c(0:(n-1))%%nfold + 1
        fold<-sample(fold, n, replace = FALSE)
    } else{
        clusters<-unique(data[,cl])
        m <- length(clusters)
        id.list<-split(1:n,data[,cl])
        cl.fold <- c(0:(m-1))%%nfold + 1
        cl.fold <- sample(cl.fold, m, replace = FALSE)
        for (i in 1:m) {
            id.list[[i]] <- rep(cl.fold[i],length(id.list[[i]]))
        }
        fold <- unlist(id.list)
    }
    
    ## preprocessing if there's fixed effects
    if (is.null(FE)==FALSE) { # with FE, first demean test data
        for (var in c(Y,D,Z)) {
            data[,var] <- fastplm(data = as.matrix(data[, var]),
                                  FE = as.matrix(data[, FE],drop = FALSE),
                                  weight = rep(1, n),
                                  FEcoefs = 0)$residuals   
        } 
    }
    ## calculation
    if (parallel == TRUE) {
        Error<-suppressWarnings(foreach(bw = grid, .combine = rbind,
                                        .export = c("coefs","getError"),
                      .inorder = FALSE) %dopar% {cv(bw)})
    } else {
        Error <- matrix(NA, length(grid), 3)
        for (i in 1:length(grid)) {
            Error[i, ] <- cv(grid[i])
            cat(".")
        } 
    } 
    colnames(Error) <- c("bandwidth","MAPE","MSPE")
    rownames(Error) <- NULL

    if (metric=="MAPE") {
        opt.bw <- grid[which.min(Error[,2])]
    } else {
        opt.bw <- grid[which.min(Error[,3])]
    } 
    output <- list(CV.out = round(Error,3),
                   opt.bw = opt.bw)
    cat(paste("Bandwidth =", round(output$opt.bw,3),"\n"))
    return(output)
}


#############

## vcovCluster.r 
## function to compute var-cov matrix using clustered robust standard errors
## inputs:
## model = model object from call to lm or glm
## cluster = vector with cluster ID indicators
## output:
## cluster robust var-cov matrix
## to call this for a model directly use:
## coeftest(model,vcov = vcovCluster(model, cluster))
## formula is similar to Stata's , cluster command

vcovCluster <- function(
    model,
    cluster
)
{
    requireNamespace("sandwich")
    ## require(lmtest)
    
    if(nrow(model.matrix(model))!=length(cluster)){
        stop("check your data: cluster variable has different N than model")
    }
    M <- length(unique(cluster))
    N <- length(cluster)           
    K <- model$rank   
    ## if(M<50){
    ##     warning("Fewer than 50 clusters, variances may be unreliable (could try block bootstrap instead).")
    ## }
    dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
    uj  <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
    rcse.cov <- dfc * sandwich(model, meat. = crossprod(uj)/N)
    ##colnames(rcse.cov)<-rownames(rcse.cov)<-names(model$coefficients)
    return(rcse.cov)
}

