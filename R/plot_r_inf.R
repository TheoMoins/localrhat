
plot_local_r <- function(chaines, r_values, theoretical_r = NULL, xlabels = NULL,
                         plot_infinite_val = TRUE, ylim, xlim=NULL, title ="",
                         plot_legend = TRUE, plot_theoretical_r = FALSE, threshold = 1.02,
                         col=c(rgb(0.1,0.1,0.7,0.7), rgb(0.7,0.1,0.1,0.7))){
  par(mar=c(4,6.5,2,2))
  if (is.null(xlim)){
    xlim = c(min(chaines)-0.01, max(chaines)+0.01)
  }
  plot(sort(chaines), r_values,
       type = "l", main = title, xlab="x", ylab="",
       xaxs="i", yaxs="i", bty = "n",
       xaxt = "n", yaxt = "n",
       cex.lab = 2, cex.main = 2,
       lwd=5, col=col[1],
       xlim = xlim, ylim=ylim)
  if (is.null(xlabels)){
    axis(2, las= 1, cex.axis=2, mgp = c(1,0.5,0), lwd=2)
  } else {
    axis(2, las= 1, at = xlabels, labels = xlabels,
         cex.axis=2, mgp = c(1,0.5,0), lwd=3)
  }
  axis(1, cex.axis=2, lwd=3, mgp = c(1,1,0))

  if(plot_infinite_val){
    abline(h=threshold, col="black", lty = 4, lwd = 4)
    abline(h=max(r_values, na.rm = TRUE), col=col[1], lty = 4, lwd = 4)
    axis(2, at = max(r_values, na.rm = TRUE), labels = expression(italic(hat(R)[infinity])),
         cex.axis = 2.5, las=1, lwd=3, mgp = c(1,3.5,0),
         col.axis=col[1], col.ticks = col[1])
  } else {
    abline(h=400, col="black", lty = 4, lwd = 4)
    abline(h=min(r_values, na.rm = TRUE), col=col[1], lty = 4, lwd = 4)
  }

  if (!is.null(theoretical_r)){
    lines(theoretical_r[,1], theoretical_r[,2], type="l", col = col[2], lwd=5)
    if (plot_legend){
      legend("topright", c(expression(italic(R(x))), expression(italic(hat(R)(x)))),
             col = c(col[2], col[1]), lty=1, cex=2, lwd = 5)
    }
    abline(h=max(theoretical_r[,2], na.rm = TRUE), col=col[2], lty = 4, lwd = 4)
    if (plot_theoretical_r){
      axis(2, at = max(theoretical_r[,2], na.rm = TRUE),
           labels =expression(italic(R[infinity])),
           cex.axis = 2.5, las=1, lwd=3, mgp = c(1,3.5,0),
           col.axis=col[2], col.ticks = col[2])
    }
  }
}

plot_hist <- function(R_matrix, colors = c("red", "blue", "green", "yellow"),
                      vaxis_pos = 0.999, bin_size = 0.1, lim_y_axis = NULL,
                      plot_legend=TRUE, xlabels = NULL, nlabels = 6,
                      plot_threshold = T, threshold = 1.02){
  par(mar=c(4,5.5,2,2))
  r_names = colnames(R_matrix)

  if (is.null(lim_y_axis)){
    lim_y_axis = length(R_matrix[,1])
  }
  if (is.null(xlabels)){
    xlim <- c(min(R_matrix), max(R_matrix))
  } else {
    xlim <- c(min(xlabels), max(xlabels))
  }

  means = round(seq(from = min(R_matrix), to = max(R_matrix), length.out = nlabels), digits = 2)

  for (i in 1:length(r_names)){
    bool <- (i!=1)
    h <- ceiling((max(R_matrix[,i])-min(R_matrix[,i])) / bin_size)+1

    b <- seq(min(R_matrix[,i]),  max(R_matrix[,i]), length.out = h)

    hist(R_matrix[,r_names[i]], col = colors[i], breaks = b,
         xlim=xlim, ylim = c(0,lim_y_axis), add=bool,
         xlab="", yaxt="n", xaxt="n", ylab="", main="", prob=FALSE)
    means <- c(means, round(mean(R_matrix[,r_names[i]]), digits=2))
  }

  axis(2, cex.axis = 2, las=1, pos = vaxis_pos, lwd = 4)
  if (is.null(xlabels)){
    axis(1, cex.axis = 2, pos=0, padj=0.5, lwd = 4)
  } else {
    axis(1, at = xlabels, labels = xlabels,
         cex.axis = 2, pos=0, padj=0.5, lwd = 4)
  }

  if (plot_threshold){
    abline(v=1.01, lty=2, lwd=4)
  }
  if (threshold > 1.01){
    abline(v=threshold, lty=2, lwd=4, col = colors[length(colors)])
    # abline(v=threshold, lty=2, lwd=4)
  }

  r_exp = c()
  for (i in 1:length(r_names)){
    if (!is.expression(r_names[i])){
      if (r_names[i] == "R-hat"){
        r_exp = c(r_exp, expression(italic(hat(R))))
      }
    }
    if (!is.expression(r_names[i])) {
      if (r_names[i] == "Rank-R-hat"){
        r_exp = c(r_exp, expression(italic(Rank-hat(R))))
      }
    }
    if (!is.expression(r_names[i])) {
      if (r_names[i] == "Brooks Multivariate R-hat"){
        r_exp = c(r_exp, expression(italic(hat(R))))
      }
    }
    if (!is.expression(r_names[i])) {
      if (r_names[i] == "R-hat-infinity"){
        r_exp = c(r_exp, expression(italic(hat(R)[infinity])))
      }
    }
    if (!is.expression(r_names[i])) {
      if (r_names[i] == "max-R-hat"){
        r_exp = c(r_exp, expression(italic(max-hat(R)[infinity])))
      }
    }
  }

  if (plot_legend){
    if (is.null(r_exp)){
      legend("topright", r_names, fill = colors, cex = 2)
    } else {
      legend("topright", r_exp, fill = colors, cex = 2)
    }
  }
}

