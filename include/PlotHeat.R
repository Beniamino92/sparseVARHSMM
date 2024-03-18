#' Alejandra Avalos Pacheco's function to splot sparsity 

#' Generate heatmap of matrices
#' @param Matrix matrix to plot
#' @param limit limit for the values of the plot
#' @export
#' @example
#' library(ggplot2)
#' sigma = matrix(rbinom(10,1,.30), 5, 2)
#' plot.heat(sigma,limit=c(-1,1))
plot.heat <- function(Matrix,Xlab="",Ylab="",Main="",limit=c(-2,2)){
  Matrix = as.matrix(Matrix)
  colnames(Matrix)<-NULL
  rownames(Matrix)<-NULL
  x = reshape2::melt(data.frame(Matrix,ind=c(nrow(Matrix):1)),id="ind")
  colnames(x)<-c("X1","X2","value")
  p_X_heat = ggplot(data = x, aes(x=X2, y=X1, fill=value)) +
    theme_bw() +
    #geom_tile(show.legend = F) +
    geom_tile() +
    xlab(Xlab) +
    ylab(Ylab) +
    ggtitle(Main) +
    theme(axis.title=element_text(size=14,face="bold"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    scale_fill_gradient2(limits=limit) + 
    theme(legend.position="bottom")
  return(p_X_heat)
}
