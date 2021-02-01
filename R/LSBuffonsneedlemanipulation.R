#
# Copyright (C) 2019 University of Amsterdam
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

LSBuffonsneedlemanipulation   <- function(jaspResults, dataset, options, state = NULL){

  # input values
  crosses <- options[["k"]]
  observations <- options[["n"]]
  a <- options[["a"]]
  b <- options[["b"]]
  lengthRatio <- options[["length"]]

  # example d for computation
  d <- 5
  l <- lengthRatio * d / 100

  # some warnings
  if(crosses > observations)
    .quitAnalysis(gettext("Warning: Number of crosses needs to be smaller than or equal to number of observations!"))


  ## Summary Table
  summaryTable <- createJaspTable(title = gettext("Summary Table"))

  summaryTable$dependOn(c("crosses", "observations", "a", "b", "lengthRatio"))
  summaryTable$addCitation("JASP Team (2018). JASP (Version 0.9.2) [Computer software].")

  summaryTable$addColumnInfo(name = "lowerCI", title = gettext("95% CI (Lower)"), type = "string")
  summaryTable$addColumnInfo(name = "Median", title = gettext("Median"),   type = "string")
  summaryTable$addColumnInfo(name = "upperCI", title = gettext("95% CI (upper)"), type = "string")


  # fill in the table
  CI95lower <- 2 * l / (qbeta(.025, crosses, observations - crosses, lower.tail = FALSE) * d)
  med <- 2 * l / (qbeta(.5, crosses, observations - crosses, lower.tail = FALSE) * d)
  CI95upper <- 2 * l / (qbeta(.975, crosses, observations - crosses, lower.tail = FALSE) * d)
  summaryTable$addRows(list(lowerCI = CI95lower, Median = med,   upperCI = CI95upper))



  ## Distribution Plot
  distPlot <- createJaspPlot(title = "Prior and Posterior distribution",  width = 480, height = 320)
  distPlot$dependOn(c("crosses", "observations", "a", "b", "lengthRatio"))
  distPlot$addCitation("JASP Team (2018). JASP (Version 0.9.2) [Computer software].")

  # values
  x <- seq(2,4,0.01)
  yPost <- 2 * l / (x^2 * d) * dbeta((2 * l / (x * d)), a + crosses, b + observations - crosses)
  yPrior <- 2 * l / (x^2 * d) * dbeta((2 * l / (x * d)), a, b)
  # axis specification

  df <- data.frame(
    x = c(x, rep(pi, 2), x),
    y = c(yPost, c(0, 1.2*max(yPost)), yPrior),
    g = factor(rep(c("Posterior for Pi", "Pi", "Implied Prior for Pi"), c(length(yPost), 2, length(yPrior))),
               levels = c("Posterior for Pi", "Pi", "Implied Prior for Pi"))
  )

  values <- factor(c(
    "Pi"                   = "red",
    "Posterior for Pi"     = "black",
    "Implied Prior for Pi" = "blue"
  ))

  distPlot$plotObject <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, y = y, group = g, color = g)) +
    ggplot2::geom_line() +
    ggplot2::scale_color_manual(name = NULL, values = levels(values)) +
    ggplot2::coord_cartesian(xlim = c(2, 4), ylim = c(0, 1.2*max(yPost))) +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw(legend.position = "right")

  jaspResults[["summaryTable"]] <- summaryTable
  jaspResults[["distPlot"]] <- distPlot

  return()
}
