#tiny modification of awesome logihist plotting function by M. de la Cruz Rot 2005
#i just replaced the breaks and labels of the first y axis on the left to save some ink
my_logihist<-function (x, y, scale.hist = 5, breaks = "Sturges", counts = TRUE, 
          intervalo = 0, ylab2 = "Frequency", fillb = 1, colob = 1, 
          sizeb = 1, pglm = FALSE, se = FALSE, sizeglm = 1, colglm = 1, 
          y1breaks=c(0, 0.5, 1), y1labs=c("0.0", "0.5", "1.0")) 
{
  labx <- deparse(substitute(x))
  laby <- deparse(substitute(y))
  id <- NULL
  if (sum(class(x) == "glm") > 0) {
    laby <- x$terms[[2]]
    labx <- x$terms[[3]]
    y <- x$data[, names(x$data) == laby]
    x <- x$data[, names(x$data) == labx]
    if (pglm == FALSE) 
      pglm <- TRUE
  }
  a <- data.frame(x, y)
  names(a) <- c(labx, laby)
  if (is.null(fillb)) 
    fillb <- NA
  if (!is.null(fillb)) {
    if (length(fillb) > 1) {
      fillb0 <- fillb[1]
      fillb1 <- fillb[2]
    }
    if (length(fillb) < 2) {
      fillb0 <- fillb1 <- fillb
    }
  }
  if (!is.null(colob)) {
    if (length(colob) > 1) {
      colob0 <- colob[1]
      colob1 <- colob[2]
    }
    if (length(colob) < 2) {
      colob0 <- colob1 <- colob
    }
  }
  if (!is.null(sizeb)) {
    if (length(sizeb) > 1) {
      sizeb0 <- sizeb[1]
      sizeb1 <- sizeb[2]
    }
    if (length(sizeb) < 2) {
      sizeb0 <- sizeb1 <- sizeb
    }
  }
  h.br <- hist(x, breaks = breaks, plot = FALSE)$br
  if (intervalo > 0) 
    h.br <- seq(from = range(h.br)[1], to = range(h.br)[2], 
                by = intervalo)
  h.x <- hist(x[y == 0], breaks = h.br, plot = FALSE)$mid
  h.y0 <- hist(x[y == 0], breaks = h.br, plot = FALSE)$counts
  h.y1 <- hist(x[y == 1], breaks = h.br, plot = FALSE)$counts
  h.y0n <- h.y0/(max(c(h.y0, h.y1)) * scale.hist)
  h.y1n <- 1 - h.y1/(max(c(h.y0, h.y1)) * scale.hist)
  datospol <- NULL
  for (i in 1:length(h.y0n)) {
    if (h.y0n[i] > 0) 
      datospol <- rbind(datospol, cbind(rep(0, 4), rep(i, 
                                                       4), c(rep(h.br[i], 2), rep(h.br[i + 1], 2)), 
                                        c(0, rep(h.y0n[i], 2), 0)))
  }
  for (i in 1:length(h.y1n)) {
    if (h.y1n[i] < 1) 
      datospol <- rbind(datospol, cbind(rep(1, 4), rep(i, 
                                                       4), c(rep(h.br[i], 2), rep(h.br[i + 1], 2)), 
                                        c(h.y1n[i], 1, 1, h.y1n[i])))
  }
  axis.hist <- function(h.y0, h.y1, scale.hist) {
    tope <- max(c(h.y0, h.y1))
    label.down <- c(0, (ceiling(tope/10)) * 5, (ceiling(tope/10)) * 
                      10)
    label.up <- c((ceiling(tope/10)) * 10, (ceiling(tope/10)) * 
                    5, 0)
    at.down <- label.down/(tope * scale.hist)
    at.up <- 1 - (label.up/(tope * scale.hist))
    at.hist <- c(at.down, at.up)
    label.hist <- c(label.down, label.up)
    return(list(at = at.hist, labels = label.hist))
  }
  datos.ax <- axis.hist(h.y0, h.y1, scale.hist)
  datospol <- data.frame(datospol)
  names(datospol) <- c("value", "id", "x", 
                       "y")
  p <- ggplot(a, aes(x, y))
  p <- p + geom_polygon(data = datospol[datospol$value == 1, 
  ], aes(x = x, y = y, group = id), fill = fillb1, colour = colob1, 
  size = sizeb1) + geom_polygon(data = datospol[datospol$value == 
                                                  0, ], aes(x = x, y = y, group = id), fill = fillb0, colour = colob0, 
                                size = sizeb0) + scale_y_continuous(breaks=y1breaks, labels = y1labs,
                                                                    sec.axis = sec_axis(trans = ~., 
                                                                                        breaks = datos.ax$at, labels = datos.ax$labels, name = ylab2)) + 
    guides(fill = FALSE) + ylab(laby) + xlab(labx)
  if (pglm == TRUE) 
    p <- p + stat_smooth(method = "glm", method.args = list(family = "binomial"), 
                         se = se, size = sizeglm, colour = colglm)
  p
}
