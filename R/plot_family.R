#' @title  plotJuxtaposition
#' @author Dieter Henrik Heiland
#' @description plotJuxtaposition
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export


plotJuxtaposition <- function(dist, points=T, span=0.2,rescale=F, quantile.cor=T, q = 0.1){
  
  if(quantile.cor==T){
    dist <- NFCN2::getCleaned(dist, feat = "aligned.target",q = q)
  }
  
  if(rescale==T){
    dist$aligned.target <- scales::rescale(dist$aligned.target, c(0,1))
  }
  
  p=ggplot2::ggplot()+ ggplot2::theme_classic()
  if(points==T){
    p=p+ggplot2::geom_point(data=dist, mapping=aes(x=dist, y=aligned.target))
  }
  p=p+geom_smooth(data=dist, mapping=aes(x=dist, y=aligned.target), se=F, span=span)
  return(p)
  
}




