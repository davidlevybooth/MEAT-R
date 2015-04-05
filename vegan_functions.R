require("vegan")

### multivariate analysis function for Shiny App

anal<-function(anal_type, biplot, data, constrained) {
    par(mfrow=c(1,1))
    if(anal_type == "PCA")
    model1 <- rda(data)
    if(anal_type == "RDA") {
       model1<-rda(data, constrained)}
    if(anal_type == "PCoA") {
      model1<-capscale(data~1, distance="bray") }
    if(anal_type == "db-RDA") {
      form <- formula(paste("data ~ ", paste(colnames(constrained), collapse= "+")))
      model1<-capscale(form, data = constrained, dist="bray") }
    eig<- eigenvals(model1); prop<-eig/sum(eig) 
    mod_plot<-plot(model1, type="n", xlab=paste(anal_type,"1",100*round(prop[1],3),"%"),ylab=paste(anal_type,"2",100*round(prop[2],3),"%"))
    mod_plot
    points(model1, display = "sites", cex=2, col="grey40",bg = "deepskyblue3", pch = 21)
    if(biplot == TRUE) { 
      text(model1, display = "species", cex=0.9, col="grey20")
      arrows(0,0, as.matrix((mod_plot$species)[,1])/1.2,as.matrix((mod_plot$species)[,2])/1.2,length=0.04, lty=1,col="grey20") 
    }   
    
    if(anal_type == "RDA") { 
      r2a_p<-round(r2p(model1), 3)
      usr <- par( "usr" )
      r2_x<-usr[ 2 ]
      r2_y<-usr[ 3 ]
      text(r2_x, r2_y, bquote(R^{2}-adj: .(r2a_p)), cex=0.9, adj = c( 1.2, -0.3 ))}   
    if(biplot == TRUE && is.null(constrained) == FALSE) {
      if(anal_type == "RDA" || anal_type == "db-RDA")
        text(model1, display = "cn", cex=0.9, col="deepskyblue4")}  
}

anal_group<-function(anal_type, biplot, groups, pal, legend, data, constrained) {
  if(anal_type == "PCA")
    model1 <- rda(data)
  if(anal_type == "RDA") {
    model1<-rda(data, constrained) }
  if(anal_type == "PCoA") {
    model1<-capscale(data~1, distance="bray") }
  if(anal_type == "db-RDA") {
    form <- formula(paste("data ~ ", paste(colnames(constrained), collapse= "+")))
    model1<-capscale(form, data = constrained, dist="bray") }
  eig<- eigenvals(model1); prop<-eig/sum(eig) 
  mod_plot<-plot(model1, type="n", xlab=paste(anal_type,"1",100*round(prop[1],3),"%"),ylab=paste(anal_type,"2",100*round(prop[2],3),"%"))
  mod_plot
  if(is.null(pal) == FALSE)
  points(model1, display = "sites", cex=2, col="grey40",bg = pal[groups], pch = 21) 
if(biplot == TRUE) { 
  text(model1, display = "species", cex=0.9, col="grey20")
  arrows(0,0, as.matrix((mod_plot$species)[,1])/1.2,as.matrix((mod_plot$species)[,2])/1.2,length=0.04, lty=1,col="grey20") }   
if(anal_type == "RDA") { 
  r2a_p<-round(r2p(model1), 3)
  usr <- par( "usr" )
  r2_x<-usr[ 2 ]
  r2_y<-usr[ 3 ]
  text(r2_x, r2_y, bquote(R^{2}-adj: .(r2a_p)), cex=0.9, adj = c( 1.2, -0.3 ))}   
if(biplot == TRUE && is.null(constrained) == FALSE) {
  if(anal_type == "RDA" || anal_type == "db-RDA")
    text(model1, display = "cn", cex=0.9, col="deepskyblue4")}
if(legend == TRUE) {
    groups <- subset(groups, duplicated(groups) != TRUE)
    legend( "topleft"
          , cex = 0.8,
          , pt.cex = 1.4,
          , y.intersp = 0.9,  
          , bty = "n", 
          , legend = groups,
          , text.col = "black",
          , col = "grey40", 
          , pt.bg = pal[groups],
          , pch = 21 
          ) }
}

anal_group_overlay<-function(anal_type, biplot, groups, pal, legend, overlay, data, constrained) {
  if(anal_type == "PCA")
    model1 <- rda(data)
  if(anal_type == "RDA") {
    model1<-rda(data, constrained) }
  if(anal_type == "PCoA") {
    model1<-capscale(data~1, distance="bray") }
  if(anal_type == "db-RDA") {
    form <- formula(paste("data ~ ", paste(colnames(constrained), collapse= "+")))
    model1<-capscale(form, data = constrained, dist="bray") }
  eig<- eigenvals(model1); prop<-eig/sum(eig) 
  mod_plot<-plot(model1, type="n", xlab=paste(anal_type,"1",100*round(prop[1],3),"%"),ylab=paste(anal_type,"2",100*round(prop[2],3),"%"))
  mod_plot
  if(is.null(pal) == FALSE)
  points(model1, display = "sites", cex=2, col="grey40",bg = pal[groups], pch = 21) 
    var_fit <- envfit(model1, overlay, add=TRUE)
    plot(var_fit, col = "deepskyblue3", add=TRUE)
  if(biplot == TRUE) { 
    text(model1, display = "species", cex=0.9, col="grey20")
    arrows(0,0, as.matrix((mod_plot$species)[,1])/1.2,as.matrix((mod_plot$species)[,2])/1.2,length=0.04, lty=1,col="grey20") }   
  if(anal_type == "RDA") { 
    r2a_p<-round(r2p(model1), 3)
    usr <- par( "usr" )
    r2_x<-usr[ 2 ]
    r2_y<-usr[ 3 ]
    text(r2_x, r2_y, bquote(R^{2}-adj: .(r2a_p)), cex=0.9, adj = c( 1.2, -0.3 ))}
  if(biplot == TRUE && is.null(constrained) == FALSE) {
    if(anal_type == "RDA" || anal_type == "db-RDA")
      text(model1, display = "cn", cex=0.9, col="deepskyblue4")}
  if(legend == TRUE) {
    groups <- subset(groups, duplicated(groups) != TRUE)
    legend( "topleft"
            , cex = 0.8,
            , pt.cex = 1.4,
            , y.intersp = 0.9,  
            , bty = "n", 
            , legend = groups,
            , text.col = "black",
            , col = "grey40", 
            , pt.bg = pal[groups],
            , pch = 21 
    ) }
}

anal_overlay<-function(anal_type, biplot, overlay, data, constrained) {
  if(anal_type == "PCA")
    model1 <- rda(data)
  if(anal_type == "RDA") {
    model1<-rda(data, constrained) }
  if(anal_type == "PCoA") {
    model1<-capscale(data~1, distance="bray") }
  if(anal_type == "db-RDA") {
    form <- formula(paste("data ~ ", paste(colnames(constrained), collapse= "+")))
    model1<-capscale(form, data = constrained, dist="bray") }
  eig<- eigenvals(model1); prop<-eig/sum(eig) 
  mod_plot<-plot(model1, type="n", xlab=paste(anal_type,"1",100*round(prop[1],3),"%"),ylab=paste(anal_type,"2",100*round(prop[2],3),"%"))
  mod_plot
  points(model1, display = "sites", cex=2, col="grey40",bg = "deepskyblue3", pch = 21) 
    var_fit <- envfit(model1, overlay, add=TRUE)
    plot(var_fit, col = "deepskyblue3", add=TRUE)
  if(biplot == TRUE) { 
    text(model1, display = "species", cex=0.9, col="grey20")
    arrows(0,0, as.matrix((mod_plot$species)[,1])/1.2,as.matrix((mod_plot$species)[,2])/1.2,length=0.04, lty=1,col="grey20") }   
  if(anal_type == "RDA") { 
    r2a_p<-round(r2p(model1), 3)
    usr <- par( "usr" )
    r2_x<-usr[ 2 ]
    r2_y<-usr[ 3 ]
    text(r2_x, r2_y, bquote(R^{2}-adj: .(r2a_p)), cex=0.9, adj = c( 1.2, -0.3 ))}
  if(biplot == TRUE && is.null(constrained) == FALSE) {
    if(anal_type == "RDA" || anal_type == "db-RDA")
      text(model1, display = "cn", cex=0.9, col="deepskyblue4")}  
}


## Scores function  for results table
scores_2<-function(anal_type, data) {
  #if(anal_type == "PCA")
  model1 <- rda(data)
   #if(anal_type == "RDA")
    #    model1<-rda(data ~ ., constrained)
  plot(model1)
}

# ## function for RDA P-symbol generation 
# ## currently not working due to timeout issues wiht anova.cca
# p_sym <- function(model){ p_val<-anova.cca(model);
#                        if (p_val <= 0.001) { x<-"***" } 
#                        else if (p_val <= 0.01) {x<-"**" } 
#                        else if (p_val <=0.05) { x<-"*" } 
#                        else x<-"NS"
#                        return(x)} 

## Get character string for R2a + p-value symbol
r2p<- function(model) {
  r2a = RsquareAdj(model)$adj.r.squared
  return(r2a)
}
