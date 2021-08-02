#' @title  plotJuxtaposition
#' @author Dieter Henrik Heiland
#' @description plotJuxtaposition
#' @inherit 
#' @return 
#' @examples 
#' 
#' @export


plotJuxtaposition <- function(dist, points=T, span){
  
  p=ggplot2::ggplot()+ ggplot2::theme_classic()
  if(points==T){
    p=p+ggplot2::geom_point(data=dist, mapping=aes(x=dist, y=aligned.target))
  }
  p=p+geom_smooth(data=dist, mapping=aes(x=dist, y=aligned.target), se=F, span=span)
  return(p)
  
}




