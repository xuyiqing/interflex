inter.plot.pool <- function( # only for discrete treatments 
  out,
  order = NULL,
  subtitle = NULL,
  legend.title = NULL,
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
  color = NULL,
  file = NULL,
  jitter = F
) {
  
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
  other_treat <- sort(all_treat[which(all_treat!=base)])
  ntreat <- length(other_treat)
  
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
   
   if(is.null(subtitle)==F){
      if(length(subtitle)!=length(other_treat)){
        stop("\"subtitle\" had a wrong length.")
      }
   }else{
	subtitle <- other_treat
   }
  
  #print(other_treat)
  
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
  
  if(is.null(legend.title)==F){
    if (is.character(legend.title) == FALSE) {
      stop("\"legend.title\" is not a string.")
    }        
  }
  
  if (!Xdistr %in% c("hist","histogram","density","none")){
    stop("\"Xdistr\" must be \"histogram\", \"density\", or \"none\".")
  }
  if (is.null(Xdistr) == TRUE) {
    Xdistr <- "density"
  } else if (!Xdistr %in% c("density","histogram","hist","none")) {
    Xdistr <- "density"
  }
  
  if (type == "binning") {
      nbins <- out$nbins
      if(nbins>1){
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
      else if(nbins==1){
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
  }
  
  else if (type == "kernel") {
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
  
  # color
  requireNamespace("RColorBrewer")
  if(ntreat>=8) {
    stop("Too many kinds of treatments. Pool plot will be too messy.")
  }
  platte <- brewer.pal(n=max(3,ntreat), "Set2")
  if(ntreat<3){
    platte <- platte[c(1:ntreat)]
  }
  
  if(is.null(color)==FALSE){
	platte <- c(color,platte)
  }
  
  # initialize
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
  

  
  # kernel estimates
  if (type == "kernel"){
      
      for(char in other_treat) {
        if(char==other_treat[1]){
          tempest <- est[[char]]
        }else{
          tempest <- rbind(tempest,est[[char]])
        }
      }
		
	  tempest$Treatment <- factor(tempest$Treatment, levels = other_treat)
	  
      p1 <- p1 + geom_line(data=tempest,aes(x = X,y = ME,color = Treatment))
	  #p1 <- p1 + scale_color_discrete(name = Dlabel, )
	  #p1 <- p1 + scale_x_discrete(limits = other_treat,labels=NULL)
      p1 <- p1 + scale_color_manual(values = platte[1:length(other_treat)],labels = subtitle)
      if (CI == TRUE) {
          p1 <- p1 + geom_ribbon(data=tempest, aes(x=X,ymin=CI_lower,ymax=CI_upper,fill = Treatment),alpha=0.2)
          p1 <- p1 + scale_fill_manual(values = platte[1:length(other_treat)],labels = subtitle)
      }

      ymin=min(yrange)-maxdiff/5
  }
  # end of kernel-specific part
  
  
  
  if (type == "binning") {
    
    for(char in other_treat) {
      if(char==other_treat[1]){
        tempest <- est.lin[[char]]
      }else{
        tempest <- rbind(tempest,est.lin[[char]])
      }
    }
	
    tempest$Treatment <- factor(tempest$Treatment, levels = other_treat)
	#print(tempest)
    p1 <- p1 + geom_line(data=tempest,aes(x = X.lvls,y = marg,color = Treatment))
    p1 <- p1 + scale_color_manual(values = platte[1:length(other_treat)], labels = subtitle)
    if (CI == TRUE) {
    p1 <- p1 + geom_ribbon(data=tempest, aes(x=X.lvls,ymin=lb,ymax=ub,fill = Treatment),alpha=0.2)
    p1 <- p1 + scale_fill_manual(values = platte[1:length(other_treat)],labels = subtitle)
    }

  ## bin estimates
    
  if(nbins>1){
  k <- 1
  for(char in other_treat) {
    tempest2 <- est.bin2[[char]]
    tempest3 <- est.bin3[[char]]
       
    if(jitter==T){
       tempest2[,'x0'] <- tempest2[,'x0'] + rnorm(dim(tempest2)[1],0,sd(tempest[,'X.lvls'])/20)
       tempest3[,'x0'] <- tempest3[,'x0'] + rnorm(dim(tempest3)[1],0,sd(tempest[,'X.lvls'])/20)
       }
    
      p1 <- p1+ geom_errorbar(data=tempest2, aes(x=x0, ymin=CI_lower, ymax=CI_upper),color = platte[k],
                            width= errorbar.width/3)+
          geom_point(data=tempest2,aes(x=x0,y=coef),size=3,shape=21,fill = platte[k],color = platte[k])
      
  
      ## in case there's non-overlap
      
      if(dim(tempest3)[1]!=0){
      p1 <- p1+geom_text(data=tempest3,aes(x=x0,y=0),label="NaN",color = platte[k])
      }
      k <- k+1
     }


    ## labels: L, M, H and so on 
      if (bin.labs == TRUE) {
        char0 <- other_treat[1]
        if (nbins==3) {
            p1 <- p1 + annotate(geom="text", x=est.bin[[char0]][1,1], y=pos,
                              label="L",colour="gray50",size=10) +
              annotate(geom="text", x=est.bin[[char0]][2,1], y=pos,
                       label="M",colour="gray50",size=10) +
              annotate(geom="text", x=est.bin[[char0]][3,1], y=pos,
                       label="H",colour="gray50",size=10)
          } 
        else if (nbins==4) {
            p1 <- p1 + annotate(geom="text", x=est.bin[[char0]][1,1], y=pos,
                              label="L",colour="gray50",size=10) +
              annotate(geom="text", x=est.bin[[char0]][2,1], y=pos,
                       label="M1",colour="gray50",size=10) +
              annotate(geom="text", x=est.bin[[char0]][3,1], y=pos,
                       label="M2",colour="gray50",size=10) +
              annotate(geom="text", x=est.bin[[char0]][4,1], y=pos,
                       label="H",colour="gray50",size=10)
        } 
        else if (nbins==2) {
          p1 <- p1 + annotate(geom="text", x=est.bin[[char0]][1,1], y=pos,
                              label="L",colour="gray50",size=10) +
            annotate(geom="text", x=est.bin[[char0]][2,1], y=pos,
                     label="H",colour="gray50",size=10)
        } 
      }
  }
     ymin=min(yrange)-maxdiff/5
  }
  # end of binning-specific part
  

  
  # plot moderate distribution
  if (Xdistr == "density") { # density plot
    ## put in data frames
    dist<-hist.out$mids[2]-hist.out$mids[1]
    deX.ymin <- min(yrange)-maxdiff/5
    deX.co <- data.frame(x = de.tr[[base]]$x,
                         y = de.tr[[base]]$y/max(de.tr[[base]]$y) * maxdiff/10 + min(yrange) - maxdiff/5)
    ## color
    feed.col<-col2rgb("gray50")
    col.co<-rgb(feed.col[1]/1000, feed.col[2]/1000,feed.col[3]/1000)
    p1 <- p1 + geom_ribbon(data = deX.co, aes(x = x, ymax = y, ymin = deX.ymin),color='gray50',
                           fill = col.co, alpha = 0.0, size=0.3)
    k <- 1
    char0 <- other_treat[1]
    start_level <- rep(deX.ymin,length(de.tr[[char0]]$x))
    for(char in other_treat){
      
      dex.tr.plot <- data.frame(x = de.tr[[char]]$x,
                           start_level = start_level,
                           end_level = de.tr[[char]]$y/max(de.tr[[char0]]$y)*maxdiff/10+start_level)
      
      p1 <- p1+geom_ribbon(data = dex.tr.plot, aes(x = x, ymax = end_level, ymin = start_level), color=platte[k],
                           alpha = 0.0,fill = platte[k],size=0.3)
      
      
      k <- k+1
    }
    p1 <- p1+geom_line(data = dex.tr.plot, aes(x = x, y = ymin), color='gray50',size=0.3)
    
    #p1 <- p1 + annotate(geom='text',x=min(hist.out$mids)+dist*6,y=ymin-maxdiff/20,label=paste0("Base Group: ",base),
    #                   colour='gray50',size=2.2)
    
}
  
  
  
  if (Xdistr %in% c("histogram","hist")) { # histogram plot
      n.hist<-length(hist.out$mids)
      dist<-hist.out$mids[2]-hist.out$mids[1]
      hist.max<-max(hist.out$counts)
      hist.col<-data.frame(ymin=rep(min(yrange)-maxdiff/5,n.hist),
                           ymax=hist.out$counts/hist.max*maxdiff/5+min(yrange)-maxdiff/5,
                           xmin=hist.out$mids-dist/2,
                           xmax=hist.out$mids+dist/2,
                           count1=count.tr[[base]]/hist.max*maxdiff/5+min(yrange)-maxdiff/5) 
      
      p1 <- p1 + geom_rect(data=hist.col,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=count1),fill='gray50',color='gray50',
                           alpha=0.3,size=0.5) 
      
      #p1 <- p1 + annotate(geom='text',x=min(hist.out$mids)+dist*6,y=ymin-maxdiff/20,label=paste0("Base Group: ",base),
       #                  colour='gray50',size=2.2)
      
      k <- 1
      start_level <- count.tr[[base]]/hist.max*maxdiff/5+min(yrange)-maxdiff/5
      for (char in other_treat){
        hist.treat<-data.frame(ymin=start_level,
                               ymax=count.tr[[char]]/hist.max*maxdiff/5+start_level,
                               xmin=hist.out$mids-dist/2,
                               xmax=hist.out$mids+dist/2)
        
        start_level <- count.tr[[char]]/hist.max*maxdiff/5+start_level
        
        p1 <- p1 + geom_rect(data=hist.treat,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=platte[k],color='gray50',
                             alpha=0.5,size=0.5)
        
        
        k <- k + 1
        }
      }
   # end of plotting X distribution  
  
 
  ## other properties
  
  if(is.null(legend.title)==F){
    p1 <- p1 + labs(fill = legend.title,color = legend.title)
  }
 
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
    ggsave(file, p1,scale = 1.5)         
  } 

  return(p1)
}


