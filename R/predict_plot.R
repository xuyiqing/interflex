predict.interflex <- function(
  out,
  order = NULL,
  subtitles = NULL,
  show.subtitles = NULL,
  #plot
  CI = NULL,
  pool = FALSE,
  main = NULL,
  Ylabel = NULL,
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
  color = NULL,
  file = NULL,
  interval = NULL,
  legend.title = NULL,
  ncols = NULL
){

  if (!class(out) %in% c("interflex")) {
    stop("Not an \"interflex\" object.")
  }
  
  if(out$predict==FALSE){
	return(0)
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
  
  if(is.null(Xlabel)==T){
    X <- out$X
  } else{X <- Xlabel}

  if(is.null(Ylabel)==T){
    Y <- out$Y
  } else{Y <- Ylabel}
  
  treat.type <- out$treat.type
  all_treat <- out$all.treat
  labelname <- out$labelname
  D.ref <- out$D.ref
  est_predict <- out$est.predict
  
  if(is.null(order)==FALSE){
	order <- as.character(order)
  }
  
  show.subtitle <- show.subtitles
  subtitle <- subtitles
  
  if(treat.type=='discrete') {
    #all_treat <- sort(all_treat[which(all_treat!=base)])
    if(is.null(order)==F){
    if(length(order)!=length(all_treat)){
      stop("\"order\" should contain all kinds of treatments.")
    }
    if(length(order)!=length(unique(order))){
        stop("\"order\" should not contain repeated values.")
    }
      
    if(sum(!is.element(order,all_treat))!=0 | sum(!is.element(all_treat,order))!=0){
      stop("\"order\" should contain all kinds of treatments.")
    }
    all_treat <- order
    }
	

    if(is.null(show.subtitle)==T){
      show.subtitle <- TRUE
   }
   
    if(is.null(subtitle)==F){
      if(length(subtitle)!=length(all_treat)){
        stop("\"subtitle\" had a wrong length.")
      }
    }
  }
  
  if(treat.type=='continuous') {
  
	if(is.null(order)==F & is.null(D.ref)==F){
    if(length(order)!=length(D.ref)){
      stop("\"order\" should contain all kinds of treatments.")
    }
    if(length(order)!=length(unique(order))){
        stop("\"order\" should not contain repeated values.")
    }
      
    if(sum(!is.element(order,D.ref))!=0 | sum(!is.element(D.ref,order))!=0){
      stop("\"order\" should contain all kinds of treatments.")
    }
	
	labelname <- c()
    for (targetD in order){
      labelname <- c(labelname,paste0("D=",round(targetD,2)))
    }
	
    all_treat <- labelname
    }
	
    if(is.null(show.subtitle)==T){
      show.subtitle <- TRUE
   }
   
    if(is.null(subtitle)==F){
      if(length(subtitle)!=length(all_treat)){
        stop("\"subtitle\" had a wrong length.")
      }
    } 
  }
  
  if (is.logical(CI) == FALSE & is.numeric(CI)==FALSE) {
    stop("\"CI\" is not a logical flag.")
  }
  
  if (is.logical(pool) == FALSE & is.numeric(pool)==FALSE) {
    stop("\"pool\" is not a logical flag.")
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
      stop("\"cex.sub\" is not numeric.")
    }
  }
  if (is.null(cex.axis)==FALSE) {
    if (is.numeric(cex.axis)==FALSE) {
      stop("\"cex.axis\" is not numeric.")
    }    
  }   
  
  ## axis labels
  if(is.null(xlab)==TRUE){
    xlab<-c(paste("Moderator: ", X, sep=""))
  } else {
    if (is.character(xlab) == FALSE) {
      stop("\"xlab\" is not a string.")
    }        
  }
  if(is.null(ylab)==TRUE){
    ylab<-c(paste("Expected Value of",Y,sep=" "))
	
  } else {
    if (is.character(ylab) == FALSE) {
      stop("\"ylab\" is not a string.")
    }        
  }
  if (is.null(main)==FALSE) {
    main <- main[1]
  }

 
 ntreat <- length(all_treat)
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
  ncols <- ntreat
 }
  

  ## axis labels
  if (is.null(cex.lab)==TRUE) {
    cex.lab <- 12
  } else {
    cex.lab <- 12 * cex.lab
  }
  if (is.null(cex.axis)==TRUE) {
    cex.axis <- 12
  } else {
    cex.axis <- 12 * cex.axis
  }
  ## title
  if (is.null(cex.main)==TRUE) {
    cex.main <- 18
  } else {
    cex.main <- 18 * cex.main
  }
  
    ## title
  if (is.null(cex.sub)==TRUE) {
    cex.sub <- 12
  } else {
    cex.sub <- 12 * cex.sub
  }
  
  ## xlim and ylim
  if (is.null(ylim)==FALSE) {
    ylim2 = c(ylim[1]-(ylim[2]-ylim[1])*0.25/6, ylim[2]+(ylim[2]-ylim[1])*0.4/6)
  }
  

  ## color  
  requireNamespace("RColorBrewer")
  platte <- rep('#999999',ntreat)
  if(ntreat>=9) {
    stop("Too many kinds of treatments. Pool plot will be too messy.")
  }
  platte <- c(brewer.pal(n=min(8,max(3,ntreat)), "Dark2"),platte)
  if(is.null(color)==FALSE){
	platte <- c(color,platte)
  }

if(pool==F){
  plot_list <- list()
  k <- 1
  for(char in all_treat){
    p <- ggplot()
    if(theme.bw==T){
      p <- p + theme_bw()
    }
	if (show.grid == FALSE) {
		p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	  }
    p <- p + geom_line(data=est_predict[[char]],aes(x=X,y=EY),size=1/ntreat,color=platte[k])
    if(CI==T){
      p <- p + geom_ribbon(data=est_predict[[char]],aes(x=X,ymin=CI_lower,ymax=CI_upper),
                           alpha=0.3,fill=platte[k])
    }
	

    if(treat.type=='discrete'){
	  if(show.subtitle==TRUE){
	  
	  if(is.null(subtitle)==TRUE){
		p <- p + labs(subtitle = paste0("Group:",char)) + theme(plot.subtitle = element_text(hjust = 0.5, size=cex.sub,
	                                                                                 lineheight=.8))
	  }
	  
	  if(is.null(subtitle)==FALSE){
		p <- p + labs(subtitle = subtitle[k]) + theme(plot.subtitle = element_text(hjust = 0.5, size=cex.sub,
	                                                                                 lineheight=.8))
	  }
	  
	  }
    }
	
    if(treat.type=='continuous'){
	if(show.subtitle==TRUE){
	if(is.null(subtitle)==TRUE){
      p <- p + labs(subtitle = paste0("Group:",labelname[k])) + theme(plot.subtitle = element_text(hjust = 0.5, size=cex.sub,
                                                                              lineheight=.8))
	}
	if(is.null(subtitle)==FALSE){
		p <- p + labs(subtitle = subtitle[k]) + theme(plot.subtitle = element_text(hjust = 0.5, size=cex.sub,
	                                                                                 lineheight=.8))
	}
      #p <- p + scale_fill_discrete(labels=labelname)
    }
	}
    k <- k+1
    plot_list[[char]] <- p
  }
}


if(pool==T){
  for(char in all_treat){
    if(char==all_treat[1]){
      tograph <- est_predict[[char]]
      tograph[,'Treatment'] <- char
    }
    else{
      tograph1 <- est_predict[[char]]
      tograph1[,'Treatment'] <- char
      tograph<- rbind(tograph,tograph1)
    }
  }
  p <- ggplot()
  if(theme.bw==T){
    p <- p + theme_bw()
  }
  if (show.grid == FALSE) {
		p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	}
	
  
  
  
  if(treat.type=='discrete'){
	#print(all_treat)
	tograph$Treatment <- factor(tograph$Treatment, levels = all_treat)
	p <- p + geom_line(data=tograph,aes(x=X,y=EY,color=Treatment),size=0.5)
    #p <- p + scale_color_manual(values = platte[1:ntreat],labels=subtitle)
	if(is.null(subtitle)==TRUE){
      p <- p + scale_color_manual(values = platte[1:ntreat])
	} else{
	  p <- p + scale_color_manual(values = platte[1:ntreat],labels=subtitle)
	}
  }
  if(treat.type=='continuous'){
	tograph$Treatment <- factor(tograph$Treatment, levels = all_treat)
	p <- p + geom_line(data=tograph,aes(x=X,y=EY,color=Treatment),size=0.5)
    #p <- p + scale_color_manual(values = platte[1:ntreat],labels=labelname)
	if(is.null(subtitle)==TRUE){
      p <- p + scale_color_manual(values = platte[1:ntreat],labels=labelname)
	} else{
	  p <- p + scale_color_manual(values = platte[1:ntreat],labels=subtitle)
	}

  }
  

  

  
  if(CI==T){
    p <- p + geom_ribbon(data=tograph,aes(x=X,ymin=CI_lower,ymax=CI_upper,fill=Treatment),
                         alpha=0.2,show.legend = T,size=0)
    if(treat.type=='discrete'){
	
	if(is.null(subtitle)==TRUE){
      p <- p + scale_fill_manual(values = platte[1:ntreat])
	} else{
	  p <- p + scale_fill_manual(values = platte[1:ntreat],labels=subtitle)
	}

    }
    if(treat.type=='continuous'){
	  if(is.null(subtitle)==TRUE){
      p <- p + scale_fill_manual(values = platte[1:ntreat],labels=labelname)
	  } else{
	  p <- p + scale_fill_manual(values = platte[1:ntreat],labels=subtitle)
	  }
    }
  }
  
  if(is.null(legend.title)==F){
    p <- p + labs(fill = legend.title,color = legend.title)
  }
  
  plot_list <- p
}
  

# estimate plot_list
  
  if(pool==T){
    ncols <- 1
  }
  
  if(pool==F){
    for(char in all_treat) {
      p1 <- plot_list[[char]]
      ## mark the original interval (in replicated papers)
      if (is.null(interval)==FALSE) {
        p1<- p1 + geom_vline(xintercept=interval,colour="steelblue", linetype=2,size=1.5)
      }
      
      ## Other universal options
      p1 <- p1 + xlab(NULL) + ylab(NULL) + 
        theme(axis.title = element_text(size=cex.lab))
      

      
      if (is.null(xlim)==FALSE & is.null(ylim)==FALSE) {
        p1<-p1+coord_cartesian(xlim = xlim, ylim = ylim2)
      }
      if (is.null(xlim)==TRUE & is.null(ylim)==FALSE) {
        p1<-p1+coord_cartesian(ylim = ylim2)
      }
      if (is.null(xlim)==FALSE & is.null(ylim)==TRUE) {
        p1<-p1+coord_cartesian(xlim = xlim)
      } 
      plot_list[[char]] <- p1
    }
    
	
	
    ## ylim
    for(char in all_treat){
      y.limits <- layer_scales(plot_list[[char]])$y$range$range
      if(char == all_treat[1]){
        ymaxmax <- y.limits[2]
        yminmin <- y.limits[1]
      }
      else{
        ymaxmax <- max(ymaxmax,y.limits[2])
        yminmin <- min(yminmin,y.limits[1])
      }
    }
    for(char in all_treat){
      plot_list[[char]] <- plot_list[[char]]+ylim(c(yminmin,ymaxmax))
    }
    
    requireNamespace("gridExtra")
    
    suppressMessages(   
      graph <- arrangeGrob(grobs=plot_list,ncol=ncols,align = 'h',
                          top = textGrob(main,gp=gpar(fontsize=cex.main,face="bold")),
                          left = textGrob(ylab, rot = 90, vjust=1,gp=gpar(fontsize=cex.lab)),
                          bottom = textGrob(xlab, gp=gpar(fontsize=cex.lab)))
    )
  }
  if(pool==T){
    p1 <- plot_list
    if (is.null(interval)==FALSE) {
      p1<- p1 + geom_vline(xintercept=interval,colour="steelblue", linetype=2,size=1.5)
    }
    
    ## Other universal options
    p1 <- p1 + xlab(xlab) + ylab(ylab) + 
      theme(axis.text = element_text(size=cex.axis), axis.title = element_text(size=cex.lab))
    
    if (is.null(main)==FALSE) {
      p1<-p1 + ggtitle(main)+ 
      theme(
        plot.title = element_text(hjust = 0.5, size=cex.main,lineheight=.8, face="bold")
      )
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
    graph <- p1
  }
  
  if (is.null(file)==FALSE) {
    ggsave(file, graph,scale =1.5)         
  }
  requireNamespace("ggplotify")
  graph <- as.ggplot(graph)
  return(graph)
  
}