#' Prepare necessary components of the IAR model for use in \code{stan}.
#'
#' The underlying `stan` code can work without explicitly calling or building the adjacency
#' matrix, but to do so requires some processing, which is done here.
#'
#' @param adjBUGS A vector of neighbor indices.
#' @param numBUGS A vector containing the number of neighbors for each point.
#' @return Returns a list where \code{J} is the number predictors, \code{J_edges} is the number
#' of edges in connected graph representing the spatial data, and \code{node1, node2} are required
#' to fit the model without explicitly specifying an adjacency matrix.
#' @note  This function is derivative from Mitzi Morris's code \insertCite{Morris:2017,Morris:2019}{ssnet}.
#' The "original" is \code{mungeCARdata4stan()}, but in that function it was assumed that the number of
#' edges was easily calculated. Some spatial settings do not adhere to this setup, and errors arise. This
#' versions, thus far, appears to work as intended in less standard settings.
#' @importFrom Rdpack reprompt
#' @references
#'
#' \insertRef{Morris:2017}{ssnet}
#'
#' \insertRef{Morris:2019}{ssnet}
#' @export
mungeCARdata4stan_irregular = function(adjBUGS,numBUGS) {
  J <- length(numBUGS)
  nn <- numBUGS
  node1 <- c()
  node2 <- c()
  iAdj <- 0
  iEdge <- 0
  for (i in 1:J) {
    for (j in 1:nn[i]) {
      iAdj <- iAdj + 1
      if (i < adjBUGS[iAdj]) {
        iEdge <- iEdge + 1
        node1[iEdge] <- i
        node2[iEdge] <- adjBUGS[iAdj]
      }
    }
  }
  return (list("J" = J, "J_edges" = length(node1),
               "node1" = node1, "node2" = node2))
}
