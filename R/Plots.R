# Plot node paths
# path_Z = a list of membership solutions (assumed to have a continuous meaning)
# lambdas = threshold parameters, which is assumed to be related to path_Z
# colors = node colors (one for each )
# linetypes = one for each node
plot_nodepaths <- function(comdet_path, lambdas, colors = NULL, linetypes = NULL) {
  
  path_z <- lapply(comdet_path$Z_path, function(x) x$Zraw)
  N <- nrow(path_z[[1]])
  K <- ncol(path_z[[1]])
  if(is.null(colors)) {
    if(K <= 3) {
      colfunc1 <- colorRampPalette(c("red", "orange"))
      colfunc2 <- colorRampPalette(c("blue", "lightblue"))
      colfunc3 <- colorRampPalette(c("darkgreen", "green"))
      memberships <- binarize_one(path_z[[1]])
      colors = rep(0, N)
      #order(table(memberships))
      colors[memberships == 1] <- colfunc1(sum(memberships ==1))
      colors[memberships == 2] <- colfunc2(sum(memberships ==2))
      colors[memberships == 3] <- colfunc3(sum(memberships ==3)) 
    }else{
      colors = terrain.colors(N)
    }
  }
    
  if(is.null(linetypes))
    linetypes <- rep("solid", N)
  
  require(ggplot2)
  df <- lapply(1:length(path_z), function(i) cbind(Reduce(rbind, path_z[[i]]), mapply(paste, rep(1:N, K), kronecker(1:K, rep(1,N))), 
                                                   rep(1:N, K), kronecker(1:K, rep(1,N)), lambdas[i] ))
  df <- Reduce(rbind, df)
  df <- as.data.frame(df)
  names(df) <- c("Membership", "NodeK", "Node", "Community", "Lambda")
  df$Membership <- as.numeric(as.character(df$Membership))
  df$Node <- factor(df$Node, levels = 1:N)
  df$Lambda <- as.numeric(as.character(df$Lambda))
  ggplot(df, aes(x = Lambda, y = Membership, col = Node)) + geom_line(aes(linetype = Node)) + 
    facet_wrap(~Community)+  scale_colour_manual(values = colors) +theme_bw() +
    theme(legend.position="none") + scale_linetype_manual(values = linetypes) +
    ggtitle("Overlapping memberships by community")
  
  # Base plot function
  #maxs = c(0,sapply(1:K,function(i) max(abs(sapply(path_z, function(x) x[,i])))))
  #plot(lambdas,abs(sapply(path_z, function(x) x[1,]))[1,],
  #     ylim = c(0,sum(maxs)*1.1),type = 'l',
  #      col = colors[1],
  #      ylab=paste("Scaled Z"),xlab = "lambda")
  # for(i in 1:K) {
  #   lim_min = 0
  #   lim_max = max(abs(sapply(path_z, function(x) x[,i])))
  #   for(j in 1:N) {
  #     lines(lambdas,(Reduce(sum,maxs[1:i]))*1.1+abs(sapply(path_z, function(x) x[j,]))[i,],
  #           col = colors[j])
  #   }
  #   abline(h = Reduce(sum,maxs[1:i])*1.1)
  # }
}

#########################################################
# Network diagram with pie plots colored by community membership
# A = network
# Z = membership matrix (rows normalized with L1 norm)
# colors = community colors
# layoutpos = positions for layout
membership_pie_plots <- function(A, Z,  colors = NULL, layoutpos = NULL,
                                 nodesizes = NULL) {
  net <- graph.adjacency(A, mode = "undirected", diag = F)
  K <- ncol(Z)
  if(is.null(colors))
    colors = rainbow(K)
  if(is.null(layoutpos))
    layoutpos <- layout.fruchterman.reingold(net)
  if(is.null(nodesizes)) {
    degree <- apply(A,1,sum)
    nodesizes <- 3*(log(degree)+1)
  }
    
  values = lapply(seq_len(nrow(Z)), function(i) Z[i,])
  shapes = ifelse(sapply(values, function(x) sum(x!=0))==1, "circle", "pie")
  cols = colors[apply(Z, 1,which.max)]
  par(mar=c(0.5,0.5,0.5,0.5))
  
  V(net)$pie.color=list(colors)
  plot.igraph(net, layout = layoutpos,
              vertex.color = cols,
              vertex.shape = shapes,
              vertex.label=NA,
              vertex.size=nodesizes, 
              vertex.pie=values)
              #xlim=c(-0.2,0.2)*3, ylim=c(-0.2,0.2)*3) #manually adjusted

}

gif.membership_pie_plots <- function(A, Zlist, filename = "gifpieplot.gif", colors = NULL, layoutpos = NULL,
                                     nodesizes = NULL) {
  net <- graph.adjacency(A, mode = "undirected", diag = F)
  K <- ncol(Zlist[[1]])
  if(is.null(colors))
    colors = rainbow(K)
  if(is.null(layoutpos))
    layoutpos <- layout.fruchterman.reingold(net)
  if(is.null(nodesizes)) {
    degree <- apply(A,1,sum)
    nodesizes <- 3*(log(degree)+1)
  }
  
  plotpie = function(i)  membership_pie_plots(A, Zlist[[i]], colors, layoutpos, nodesizes)
  
  # requires 
  dir.create("gif_temp")
  setwd("gif_temp")
  #https://www.r-bloggers.com/animate-gif-images-in-r-imagemagick/
  png(file="gifpieplot%02d.png", width=1000, height=1000)
  for (i in 1:length(Zlist)){
    plotpie(i)
  }
  dev.off()
  # convert the .png files to one .gif file using ImageMagick. 
  # The system() function executes the command as if it was done
  # in the terminal. the -delay flag sets the time between showing
  # the frames, i.e. the speed of the animation.
  system(paste("magick -delay 10 *.png ", filename))
  # to not leave the directory with the single jpeg files
  # I remove them.
  file.remove(list.files(pattern=".png"))
  setwd("..")
}
