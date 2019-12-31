plot.interflex <- function(
  out,
  order = NULL,
  subtitles = NULL,
  show.subtitles = NULL,
  CI = NULL,
  diff.values = NULL,
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
  cex.sub = NULL,
  cex.lab = NULL,
  cex.axis = NULL,
  bin.labs = TRUE, # bin labels    
  interval = NULL, # interval in replicated papers
  file = NULL,
  ncols = NULL,
  #pool plot
  pool = FALSE,
  legend.title = NULL,
  color = NULL,
  jitter = FALSE
){
  if (!class(out) %in% c("interflex")) {
    stop("Not an \"interflex\" object.")
  } 


  if(out$type=='linear'){
	out$type <- 'binning'
	out$nbins <- 1
  }
  
  if(out$type=='kernel'){
	if(is.null(CI)==TRUE){
		CI <- out$CI
	}
  }
  
  if(out$type=='binning'){
	if(is.null(CI)==TRUE){
		CI <- TRUE
	}
  }
  
  treat.type <- out$treat.type
  if(is.null(order)==FALSE){
	order <- as.character(order)
  }
  
  show.subtitle <- show.subtitles
  subtitle <- subtitles
  
  if(pool==TRUE & treat.type=='discrete'){
	graph <- inter.plot.pool(out=out,
	order = order,
	subtitle = subtitle,
	CI = CI,
	Xdistr = Xdistr,
	main = main,
	Ylabel = Ylabel,
	Dlabel = Dlabel,
	Xlabel = Xlabel,
	xlab = xlab,
	ylab = ylab,
	xlim = xlim,
	ylim = ylim,
	theme.bw = theme.bw,
	show.grid = show.grid,     
	cex.main = cex.main,
	cex.lab = cex.lab,
	cex.axis = cex.axis,
	bin.labs = bin.labs, # bin labels    
	interval = interval, # interval in replicated papers
	legend.title = legend.title,
	color = color,
	file = file,
	jitter = jitter)
	
	return(graph)
  }
  
  
  
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
  
  if(is.null(Ylabel)==T){
	Ylabel <- out$Ylabel
  }
  if(is.null(Xlabel)==T){
	Xlabel <- out$Xlabel
  }
  if(is.null(Dlabel)==T){
	Dlabel <- out$Dlabel
  }


  de <- out$de
  base <- out$base
  treat.type <- out$treat.type
  all_treat <- out$treatlevels
  
  de.tr <- out$de.tr
  hist.out <- out$hist.out
  count.tr <- out$count.tr
  
  if(treat.type=='discrete') {
    other_treat <- sort(all_treat[which(all_treat!=base)])
    if(is.null(order)==F){
    if(length(order)!=length(other_treat)){
      stop("\"order\" should contain all kinds of treatments except for the base group.")
    }
    if(length(order)!=length(unique(order))){
        stop("\"order\" should not contain repeated values.")
    }
      
    if(sum(!is.element(order,other_treat))!=0 | sum(!is.element(other_treat,order))!=0){
      stop("\"order\" should contain all kinds of treatments except for the base group.")
    }
    other_treat <- order
    }
	
    if(is.null(show.subtitle)==T){
    if(length(other_treat)==1){
      show.subtitle <- FALSE
    }
    else{
      show.subtitle <- TRUE
    }
      
    if(is.null(subtitle)==F){
      show.subtitle <- TRUE
    }
    }
    

    if(is.null(subtitle)==F){
      if(length(subtitle)!=length(other_treat)){
        stop("\"subtitle\" had a wrong length.")
      }
    }
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
  if (is.null(cex.sub)==FALSE) {
    if (is.numeric(cex.sub)==FALSE) {
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
  
  if (is.null(ncols) == FALSE) {
    if (ncols%%1 != 0) {
      stop("\"ncols\" is not a positive integer.")
    } else {
      ncols <- ncols[1]
    }
    if (ncols < 1) {
      stop("\"ncols\" is not a positive integer.")
    }
  } else{
  if(treat.type=="discrete"){
	ncols <- length(other_treat)
  }
  if(treat.type=="continuous"){
	ncols <- 1
  }
  }
  treat_sc <- max(1,ncols-1)
  
  if(is.null(diff.values)==FALSE){
	if(is.numeric(diff.values)==FALSE){
		stop("\"diff.values\" is not numeric.")
	}
	if(length(diff.values)!=2){
		stop("\"diff.values\" must be of length 2.")
	}
	
	if(treat.type=='discrete' & type=='binning'){
		tempxx <- out$est.lin[[other_treat[1]]][,'X.lvls']
	}
	if(treat.type=='discrete' & type=='kernel'){
		tempxx <- out$est[[other_treat[1]]][,'X']
	}
	if(treat.type=='continuous' & type=='binning'){
		tempxx <- out$est.lin[,'X.lvls']
	}
	if(treat.type=='continuous' & type=='kernel'){
		tempxx <- out$est[,'X']
	}
	min.XX <- min(tempxx)
	max.XX <- max(tempxx)
	for(a in diff.values){
		if(a<min.XX|a>max.XX){
			stop("Elements in \"diff.values\" should be larger than the minimum of moderator and less than the maximum of it.")
		}
	}
  }
  

  #########################################
  
  if (type == "binning") {
    nbins <- out$nbins
    if(nbins>1){
    if (treat.type=='discrete'){
      est.lin <- out$est.lin
      est.bin <- out$est.bin
      est.bin2 <- list() ## non missing part
      est.bin3 <- list() ## missing part
      yrange <- c()
      
      for(char in other_treat) {
        est.bin2[[char]] <- est.bin[[char]][which(is.na(est.bin[[char]][,2])==FALSE),] 
        est.bin3[[char]] <- est.bin[[char]][which(is.na(est.bin[[char]][,2])==TRUE),] 
        yrange <- c(yrange,na.omit(unlist(c(est.lin[[char]][,c(4,5)],est.bin[[char]][,c(4,5)]))))
      }
      
      if (is.null(ylim)==FALSE) {yrange<-c(ylim[2],ylim[1]+(ylim[2]-ylim[1])*1/8)}
      X.lvls <- est.lin[[other_treat[1]]][,1]
      errorbar.width<-(max(X.lvls)-min(X.lvls))/20
      maxdiff<-(max(yrange)-min(yrange))
      pos<-max(yrange)-maxdiff/20
    }
    
    if (treat.type=='continuous'){
      nbins <- out$nbins
      est.lin <- out$est.lin
      est.bin<-out$est.bin
      est.bin2<-est.bin[which(is.na(est.bin[,2])==FALSE),] ## non missing part
      est.bin3<-est.bin[which(is.na(est.bin[,2])==TRUE),]  ## missing part
      X.lvls <- est.lin[,1]
      errorbar.width<-(max(X.lvls)-min(X.lvls))/20
      yrange<-na.omit(unlist(c(est.lin[,c(3,4)],est.bin[,c(4,5)])))
      if (is.null(ylim)==FALSE) {yrange<-c(ylim[2],ylim[1]+(ylim[2]-ylim[1])*1/8)}
      maxdiff<-(max(yrange)-min(yrange))
      pos<-max(yrange)-maxdiff/20
    }
    }
    else if(nbins==1){
      if (treat.type=='discrete'){
        est.lin <- out$est.lin
        yrange <- c()
        
        for(char in other_treat) {
          yrange <- c(yrange,na.omit(unlist(est.lin[[char]][,c(4,5)])))
        }
        
        if (is.null(ylim)==FALSE) {yrange<-c(ylim[2],ylim[1]+(ylim[2]-ylim[1])*1/8)}
        X.lvls <- est.lin[[other_treat[1]]][,1]
        errorbar.width<-(max(X.lvls)-min(X.lvls))/20
        maxdiff<-(max(yrange)-min(yrange))
        pos<-max(yrange)-maxdiff/20
      }
      
      if (treat.type=='continuous'){
        est.lin <- out$est.lin
        X.lvls <- est.lin[,1]
        errorbar.width<-(max(X.lvls)-min(X.lvls))/20
        yrange<-na.omit(unlist(est.lin[,c(3,4)]))
        if (is.null(ylim)==FALSE) {yrange<-c(ylim[2],ylim[1]+(ylim[2]-ylim[1])*1/8)}
        maxdiff<-(max(yrange)-min(yrange))
        pos<-max(yrange)-maxdiff/20
      }
    }
  } 
  
  else if (type == "kernel") {
    if(treat.type=='continuous'){
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
    
    if(treat.type=='discrete'){
      est <- out$est
      # Checking options
      if (is.logical(CI) == FALSE & is.numeric(CI)==FALSE) {
        stop("\"CI\" is not a logical flag.")
      }
      yrange <- c()
      for(char in other_treat) {
        tempest <- est[[char]]
        
        if(CI==TRUE) {
          yrange <- c(yrange,na.omit(c(tempest$CI_lower,tempest$CI_upper)))
        }
        else {
          yrange<-c(yrange,na.omit(tempest$ME))
        }
      }
      
      if (is.null(ylim)==FALSE) {yrange<-c(ylim[2],ylim[1]+(ylim[2]-ylim[1])*1/8)}
      maxdiff<-(max(yrange)-min(yrange)) 
    }
  }
  
  
  # initialization
  if(treat.type=='continuous'){
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
  }
  
  if(treat.type=='discrete') {
    p.group <- list()
    for (char in other_treat) {
      
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
      p.group[[char]] <- p1
    }
  }
  
  ## plotting moderator distribution      
  if (Xdistr == "density") { # density plot
    
    if (treat.type == 'discrete') { ## binary D
      ## put in data frames
      deX.ymin <- min(yrange)-maxdiff/5
      
      deX.co <- data.frame(x = de.tr[[base]]$x,
                           y = de.tr[[base]]$y/max(de.tr[[base]]$y) * maxdiff/5 + min(yrange) - maxdiff/5)
      
      ## color
      feed.col<-col2rgb("gray50")
      col.co<-rgb(feed.col[1]/1000, feed.col[2]/1000,feed.col[3]/1000)
      col.tr<-rgb(red=1, blue=0, green=0)
    
      for(char in other_treat){
        deX.tr <- data.frame(x = de.tr[[char]]$x,
                             y = de.tr[[char]]$y/max(de.tr[[char]]$y) * maxdiff/5 + min(yrange) - maxdiff/5)
        
        p1 <- p.group[[char]] + geom_ribbon(data = deX.co, aes(x = x, ymax = y, ymin = deX.ymin),
                               fill = col.co, alpha = 0.2) +
          geom_ribbon(data = deX.tr, aes(x = x, ymax = y, ymin = deX.ymin),
                      fill = col.tr, alpha = 0.2)
        
        p.group[[char]] <- p1
      }
    } 
    
    if(treat.type=='continuous') { ## continuous D
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
  }
  
  else if (Xdistr %in% c("histogram","hist")) { # histogram plot
    if (treat.type=='discrete') { ## binary D
      n.hist<-length(hist.out$mids)
      dist<-hist.out$mids[2]-hist.out$mids[1]
      hist.max<-max(hist.out$counts)
      
      hist.col<-data.frame(ymin=rep(min(yrange)-maxdiff/5,n.hist),
                        #ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
                        xmin=hist.out$mids-dist/2,
                        xmax=hist.out$mids+dist/2,
                        count1=count.tr[[base]]/hist.max*maxdiff/5+min(yrange)-maxdiff/5) 
      
      for (char in other_treat){
        
        hist.treat<-data.frame(ymin=hist.col[,'count1'],
                          #ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
                          xmin=hist.out$mids-dist/2,
                          xmax=hist.out$mids+dist/2,
                          count1=count.tr[[char]]/hist.max*maxdiff/5+hist.col[,'count1']) 
        
        p1 <- p.group[[char]] + geom_rect(data=hist.col,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=count1),fill="gray50",
                             colour="gray50",alpha=0.3,size=0.3) + # control
          geom_rect(data=hist.treat,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=count1),
                    fill="red",colour="gray50",alpha=0.3,size=0.3)
        
       
        p.group[[char]] <- p1
      }
    } 
    
    if(treat.type=='continuous') { ## continuous D
      
      n.hist<-length(hist.out$mids)
      dist<-hist.out$mids[2]-hist.out$mids[1]
      hist.max<-max(hist.out$counts)            
      histX<-data.frame(ymin=rep(min(yrange)-maxdiff/5,n.hist),
                        ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
                        xmin=hist.out$mids-dist/2,
                        xmax=hist.out$mids+dist/2)
      
      p1 <- p1 + geom_rect(data=histX,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                           colour="gray50",alpha=0.3,size=0.5)
    }
  } # end of plotting X distribution  
  
  


  # kernel estimates
  if (type == "kernel") { 
    if(treat.type=='continuous'){
    ## point estimates
    p1 <-  p1 + geom_line(data=est,aes(X,ME))
    ## confidence interval
    if (CI == TRUE) {
      p1 <- p1 + geom_ribbon(data=est, aes(x=X,ymin=CI_lower,ymax=CI_upper),alpha=0.2)
    }
	if(is.null(diff.values)==FALSE){
		for(target.value in diff.values){
			Xnew<-abs(est[,'X']-target.value)
			d1<-min(Xnew)     
			label1<-which.min(Xnew)
			Xnew[label1]<-Inf
			d2<-min(Xnew)     
			label2<-which.min(Xnew)
			if(d1==0){
				est.mark <- est[label1,"ME"]
				if(CI==TRUE){
					lb.mark <- est[label1,"CI_lower"]
					ub.mark <- est[label1,"CI_upper"]
				}
			}  
			else if(d2==0){
				est.mark <- est[label2,"ME"]
				if(CI==TRUE){
					lb.mark <- est[label2,"CI_lower"]
					ub.mark <- est[label2,"CI_upper"]
				}
			} 
			else{ ## weighted average
				est.mark1 <- est[label1,"ME"]
				est.mark2 <- est[label2,"ME"]
				est.mark <- ((est.mark1 * d2 + est.mark2 * d1)/(d1 + d2))
				if(CI==TRUE){
					lb.mark1 <- est[label1,"CI_lower"]
					ub.mark1 <- est[label1,"CI_upper"]
					lb.mark2 <- est[label2,"CI_lower"]
					ub.mark2 <- est[label2,"CI_upper"]
					lb.mark <- ((lb.mark1 * d2 + lb.mark2 * d1)/(d1 + d2))
					ub.mark <- ((ub.mark1 * d2 + ub.mark2 * d1)/(d1 + d2))
				}
			}
			p1 <- p1 + annotate("point",x=target.value,y=est.mark,size=1,colour='red')
			if(CI==TRUE){
				p1 <- p1+ annotate("errorbar",x=target.value,ymin=lb.mark,ymax=ub.mark,colour='red',size=0.5,width= (max(tempxx)-min(tempxx))/30)
			}
			}
		}
    }
    if(treat.type=='discrete'){
      for(char in other_treat) {
        p1 <- p.group[[char]]
        tempest <- est[[char]]
        p1 <-  p1 + geom_line(data=tempest,aes(X,ME))
        if (CI == TRUE) {
          p1 <- p1 + geom_ribbon(data=tempest, aes(x=X,ymin=CI_lower,ymax=CI_upper),alpha=0.2)
        }
        ymin=min(yrange)-maxdiff/5
		
		if(is.null(diff.values)==FALSE){
			for(target.value in diff.values){
				Xnew<-abs(tempest[,'X']-target.value)
				d1<-min(Xnew)     
				label1<-which.min(Xnew)
				Xnew[label1]<-Inf
				d2<-min(Xnew)     
				label2<-which.min(Xnew)
				if(d1==0){
					est.mark <- tempest[label1,"ME"]
					if(CI==TRUE){
						lb.mark <- tempest[label1,"CI_lower"]
						ub.mark <- tempest[label1,"CI_upper"]
					}
				}  
				else if(d2==0){
					est.mark <- tempest[label2,"ME"]
					if(CI==TRUE){
						lb.mark <- tempest[label2,"CI_lower"]
						ub.mark <- tempest[label2,"CI_upper"]
					}
				} 
				else{ ## weighted average
					est.mark1 <- tempest[label1,"ME"]
					est.mark2 <- tempest[label2,"ME"]
					est.mark <- ((est.mark1 * d2 + est.mark2 * d1)/(d1 + d2))
					if(CI==TRUE){
						lb.mark1 <- tempest[label1,"CI_lower"]
						ub.mark1 <- tempest[label1,"CI_upper"]
						lb.mark2 <- tempest[label2,"CI_lower"]
						ub.mark2 <- tempest[label2,"CI_upper"]
						lb.mark <- ((lb.mark1 * d2 + lb.mark2 * d1)/(d1 + d2))
						ub.mark <- ((ub.mark1 * d2 + ub.mark2 * d1)/(d1 + d2))
					}
				}
				p1 <- p1 + annotate("point",x=target.value,y=est.mark,size=1,colour='red')
				if(CI==TRUE){
					p1 <- p1+ annotate("errorbar",x=target.value,ymin=lb.mark,ymax=ub.mark,colour='red',size=0.5,width= (max(tempxx)-min(tempxx))/30)
				}
			}
		}
        p.group[[char]] <- p1
      }
  } # end of kernel-specific part
  }
  
  
  ## binning estimates and bin labels
  if (type == "binning") {
    if(nbins>1){
    if(treat.type == 'continuous') {
    
    ## linear plot 
    p1<-p1 + geom_line(data=est.lin,aes(X.lvls,marg))
	
	if(CI==T){
      p1 <- p1 + geom_ribbon(data=est.lin, aes(x=X.lvls,ymin=lb,ymax=ub),alpha=0.2)
	  }
    ## bin estimates
    p1<-p1+ geom_errorbar(data=est.bin2, aes(x=x0, ymin=CI_lower, ymax=CI_upper),colour="red",
                          width= errorbar.width)+
      geom_point(data=est.bin2,aes(x0,coef),size=4/treat_sc,shape=21,fill="white",colour="red") 
    
    ## in case there's non-overlap
    p1<-p1+annotate(geom="text", x=est.bin3[,1], y=rep(0,dim(est.bin3)[1]),
                    label="NaN",colour="red")  
    ## labels: L, M, H and so on 
    if (bin.labs == TRUE) {
      if (nbins==3) {
        p1<-p1 + annotate(geom="text", x=est.bin[1,1], y=pos,
                          label="L",colour="gray50",size=10/treat_sc) +
          annotate(geom="text", x=est.bin[2,1], y=pos,
                   label="M",colour="gray50",size=10/treat_sc) +
          annotate(geom="text", x=est.bin[3,1], y=pos,
                   label="H",colour="gray50",size=10/treat_sc)
      } else if (nbins==4) {
        p1<-p1 + annotate(geom="text", x=est.bin[1,1], y=pos,
                          label="L",colour="gray50",size=10/treat_sc) +
          annotate(geom="text", x=est.bin[2,1], y=pos,
                   label="M1",colour="gray50",size=10/treat_sc) +
          annotate(geom="text", x=est.bin[3,1], y=pos,
                   label="M2",colour="gray50",size=10/treat_sc) +
          annotate(geom="text", x=est.bin[4,1], y=pos,
                   label="H",colour="gray50",size=10/treat_sc)
      } 
      else if (nbins==2) {
        p1<-p1 + annotate(geom="text", x=est.bin[1,1], y=pos,
                          label="L",colour="gray50",size=10/treat_sc) +
          annotate(geom="text", x=est.bin[2,1], y=pos,
                   label="H",colour="gray50",size=10/treat_sc)
      }
    }
    }
    
    if (treat.type == 'discrete'){
      
      for(char in other_treat) {
        
        p1 <- p.group[[char]]
        
        p1<-p1 + geom_line(data=est.lin[[char]],aes(X.lvls,marg))
		
		if(CI==TRUE){
         p1 <- p1 + geom_ribbon(data=est.lin[[char]], aes(x=X.lvls,ymin=lb,ymax=ub),alpha=0.2)
		 }
        ## bin estimates
        p1<-p1+ geom_errorbar(data=est.bin2[[char]], aes(x=x0, ymin=CI_lower, ymax=CI_upper),colour="red",
                              width= errorbar.width)+
          geom_point(data=est.bin2[[char]],aes(x0,coef),size=2/treat_sc,shape=21,fill="white",colour="red") 
        
        ## in case there's non-overlap
        p1<-p1+annotate(geom="text", x=est.bin3[[char]][,1], y=rep(0,dim(est.bin3[[char]])[1]),
                        label="NaN",colour="red")  
        ## labels: L, M, H and so on 
        if (bin.labs == TRUE) {
          if (nbins==3) {
            p1<-p1 + annotate(geom="text", x=est.bin[[char]][1,1], y=pos,
                              label="L",colour="gray50",size=10/treat_sc) +
              annotate(geom="text", x=est.bin[[char]][2,1], y=pos,
                       label="M",colour="gray50",size=10/treat_sc) +
              annotate(geom="text", x=est.bin[[char]][3,1], y=pos,
                       label="H",colour="gray50",size=10/treat_sc)
          } 
          else if (nbins==4) {
            p1<-p1 + annotate(geom="text", x=est.bin[[char]][1,1], y=pos,
                              label="L",colour="gray50",size=10/treat_sc) +
              annotate(geom="text", x=est.bin[[char]][2,1], y=pos,
                       label="M1",colour="gray50",size=10/treat_sc) +
              annotate(geom="text", x=est.bin[[char]][3,1], y=pos,
                       label="M2",colour="gray50",size=10/treat_sc) +
              annotate(geom="text", x=est.bin[[char]][4,1], y=pos,
                       label="H",colour="gray50",size=10/treat_sc)
          }
          else if (nbins==2) {
            p1<-p1 + annotate(geom="text", x=est.bin[[char]][1,1], y=pos,
                              label="L",colour="gray50",size=10/treat_sc) +
              annotate(geom="text", x=est.bin[[char]][2,1], y=pos,
                       label="H",colour="gray50",size=10/treat_sc)
          }
        }
        ymin=min(yrange)-maxdiff/5
        #p1 <- p1 + annotate(geom='text',x=median(X.lvls),y=ymin+maxdiff/20,label=char,colour='red',size=5/treat_sc)
        #p1 <- p1 + annotate(geom='text',x=median(X.lvls),y=ymin-maxdiff/20,label=base,colour='blue',size=5/treat_sc)
        p.group[[char]] <- p1
      }
    }
    }
    
	if(nbins==1){
      if(treat.type == 'continuous') {
        ## linear plot 
        p1<-p1 + geom_line(data=est.lin,aes(X.lvls,marg))
		
		if(CI==TRUE){
		 p1 <- p1+geom_ribbon(data=est.lin, aes(x=X.lvls,ymin=lb,ymax=ub),alpha=0.2)
		}
		tempest <- est.lin
		if(is.null(diff.values)==FALSE){
			for(target.value in diff.values){
				Xnew<-abs(tempest[,'X.lvls']-target.value)
				d1<-min(Xnew)     
				label1<-which.min(Xnew)
				Xnew[label1]<-Inf
				d2<-min(Xnew)     
				label2<-which.min(Xnew)
				if(d1==0){
					est.mark <- tempest[label1,"marg"]
					if(CI==TRUE){
						lb.mark <- tempest[label1,"lb"]
						ub.mark <- tempest[label1,"ub"]
					}
				}  
				else if(d2==0){
					est.mark <- tempest[label2,"marg"]
					if(CI==TRUE){
						lb.mark <- tempest[label2,"lb"]
						ub.mark <- tempest[label2,"ub"]
					}
				} 
				else{ ## weighted average
					est.mark1 <- tempest[label1,"marg"]
					est.mark2 <- tempest[label2,"marg"]
					est.mark <- ((est.mark1 * d2 + est.mark2 * d1)/(d1 + d2))
					if(CI==TRUE){
						lb.mark1 <- tempest[label1,"lb"]
						ub.mark1 <- tempest[label1,"ub"]
						lb.mark2 <- tempest[label2,"lb"]
						ub.mark2 <- tempest[label2,"ub"]
						lb.mark <- ((lb.mark1 * d2 + lb.mark2 * d1)/(d1 + d2))
						ub.mark <- ((ub.mark1 * d2 + ub.mark2 * d1)/(d1 + d2))
					}
				}
				p1 <- p1 + annotate("point",x=target.value,y=est.mark,size=1,colour='red')
				if(CI==TRUE){
					p1 <- p1+ annotate("errorbar",x=target.value,ymin=lb.mark,ymax=ub.mark,colour='red',size=0.5,width= (max(tempxx)-min(tempxx))/30)
				}
			}
		}
      }
      
      if (treat.type == 'discrete'){
        
        for(char in other_treat) {
          
          p1 <- p.group[[char]]
          
          p1<-p1 + geom_line(data=est.lin[[char]],aes(X.lvls,marg))
		  
		  if(CI==T){
            p1 <- p1+geom_ribbon(data=est.lin[[char]], aes(x=X.lvls,ymin=lb,ymax=ub),alpha=0.2)
          }
          ymin=min(yrange)-maxdiff/5
		  tempest <- est.lin[[char]]
		  if(is.null(diff.values)==FALSE){
			for(target.value in diff.values){
				Xnew<-abs(tempest[,'X.lvls']-target.value)
				d1<-min(Xnew)     
				label1<-which.min(Xnew)
				Xnew[label1]<-Inf
				d2<-min(Xnew)     
				label2<-which.min(Xnew)
				if(d1==0){
					est.mark <- tempest[label1,"marg"]
					if(CI==TRUE){
						lb.mark <- tempest[label1,"lb"]
						ub.mark <- tempest[label1,"ub"]
					}
				}  
				else if(d2==0){
					est.mark <- tempest[label2,"marg"]
					if(CI==TRUE){
						lb.mark <- tempest[label2,"lb"]
						ub.mark <- tempest[label2,"ub"]
					}
				} 
				else{ ## weighted average
					est.mark1 <- tempest[label1,"marg"]
					est.mark2 <- tempest[label2,"marg"]
					est.mark <- ((est.mark1 * d2 + est.mark2 * d1)/(d1 + d2))
					if(CI==TRUE){
						lb.mark1 <- tempest[label1,"lb"]
						ub.mark1 <- tempest[label1,"ub"]
						lb.mark2 <- tempest[label2,"lb"]
						ub.mark2 <- tempest[label2,"ub"]
						lb.mark <- ((lb.mark1 * d2 + lb.mark2 * d1)/(d1 + d2))
						ub.mark <- ((ub.mark1 * d2 + ub.mark2 * d1)/(d1 + d2))
					}
				}
				p1 <- p1 + annotate("point",x=target.value,y=est.mark,size=1,colour='red')
				if(CI==TRUE){
					p1 <- p1+ annotate("errorbar",x=target.value,ymin=lb.mark,ymax=ub.mark,colour='red',size=0.5,width= (max(tempxx)-min(tempxx))/30)
				}
			}
		}
			
        p.group[[char]] <- p1
        }
      }
    }
  } # end of binning-specific case
  
  
  if(treat.type=='continuous'){
  
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
    cex.axis <- 15/treat_sc
  } else {
    cex.axis <- 15 * cex.axis/treat_sc
  }
  p1 <- p1 + xlab(xlab) + ylab(ylab) + 
    theme(axis.text = element_text(size=cex.axis), axis.title = element_text(size=cex.lab))
  
  ## title
  if (is.null(cex.main)==TRUE) {
    cex.main <- 18
  } else {
    cex.main <- 18 * cex.main
  }
  
  if (is.null(cex.sub)==TRUE) {
    cex.sub <- 12
  } else {
    cex.sub <- 12 * cex.sub
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
    ggsave(file, plot = p1,scale=1.2)         
  } 
  ######################################################
  
  graph <- p1
  }
  
  if(treat.type=='discrete')
  {
    ## axis labels
    if (is.null(cex.lab)==TRUE) {
      cex.lab <- 15
    } else {
      cex.lab <- 15 * cex.lab
    }
    if (is.null(cex.axis)==TRUE) {
      cex.axis <- 15/treat_sc
    } else {
      cex.axis <- 15 * cex.axis/treat_sc
    }
    ## title
    if (is.null(cex.main)==TRUE) {
      cex.main <- 18
    } else {
      cex.main <- 18 * cex.main
    }
    
    if (is.null(cex.sub)==TRUE) {
      cex.sub <- 12
    } else {
      cex.sub <- 12 * cex.sub
    }
    
    ## xlim and ylim
    if (is.null(ylim)==FALSE) {
      ylim2 = c(ylim[1]-(ylim[2]-ylim[1])*0.25/6, ylim[2]+(ylim[2]-ylim[1])*0.4/6)
    }
    
    k <- 1
    for(char in other_treat) {
      p1 <- p.group[[char]]
      ## mark the original interval (in replicated papers)
      if (is.null(interval)==FALSE) {
        p1<- p1 + geom_vline(xintercept=interval,colour="steelblue", linetype=2,size=1.5)
      }
      
      ## Other universal options
      
      p1 <- p1 + xlab(NULL) + ylab(NULL) + 
        theme(axis.text = element_text(size=cex.axis))

      if(show.subtitle==T){
      if(is.null(subtitle)==T){
          subtitle.temp <- paste0("Treat Group = ",char,", Base Group = ",base)
	        p1 <- p1 + labs(subtitle = subtitle.temp) + theme(plot.subtitle = element_text(hjust = 0.5, size=cex.sub,
	                                                                                 lineheight=.8))
      }
      if(is.null(subtitle)==F){
        subtitle.temp<- subtitle[k]
        p1 <- p1 + labs(subtitle = subtitle.temp) + theme(plot.subtitle = element_text(hjust = 0.5, size=cex.sub,
                                                                                  lineheight=.8))
      }
      }
      #p1 <- p1 + labs(subtitle = subtitle)

      if (is.null(xlim)==FALSE & is.null(ylim)==FALSE) {
        p1<-p1+coord_cartesian(xlim = xlim, ylim = ylim2)
      }
      if (is.null(xlim)==TRUE & is.null(ylim)==FALSE) {
        p1<-p1+coord_cartesian(ylim = ylim2)
      }
      if (is.null(xlim)==FALSE & is.null(ylim)==TRUE) {
        p1<-p1+coord_cartesian(xlim = xlim)
      } 
      

      ######################################################
      p.group[[char]] <- p1
      k <- k+1
    }
    
    ## ylim
    for(char in other_treat){
      y.limits <- layer_scales(p.group[[char]])$y$range$range
      if(char == other_treat[1]){
        ymaxmax <- y.limits[2]
        yminmin <- y.limits[1]
      }
      else{
        ymaxmax <- max(ymaxmax,y.limits[2])
        yminmin <- min(yminmin,y.limits[1])
      }
    }
    for(char in other_treat){
      p.group[[char]] <- p.group[[char]]+ylim(c(yminmin,ymaxmax))
    }
    
    requireNamespace("gridExtra")
    requireNamespace("grid")
	requireNamespace("ggplotify")
    suppressMessages(
      graph <- arrangeGrob(grobs=p.group,ncol=ncols,align = 'h',
                            top = textGrob(main,gp=gpar(fontsize=cex.main,face="bold")),
                            left = textGrob(ylab, rot = 90, vjust=1,gp=gpar(fontsize=cex.lab)),
                            bottom = textGrob(xlab, gp=gpar(fontsize=cex.lab))
                            )
      )
    
  ## save to file
    if (is.null(file)==FALSE) {
      ggsave(file, graph,scale = 1.2)         
    } 
	graph <- as.ggplot(graph)
  }
  
  return(graph)
}




