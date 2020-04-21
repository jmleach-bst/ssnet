#' Prepare necessary components of the IAR model for use in \code{stan}.
#'
#' The underlying `stan` code can work without explicitly calling or building the adjacency
#' matrix, but to do so requires some processing, which is done here.
#'
#' @param adjBUGS A vector of neighbor indices.
#' @param numBUGS A vector containing the number of neighbors for each point.
#' @return Returns a list where \code{J} is the number predictors, \code{J_edges} is the number
#' of edges in connected graph representing the spatial data, and \code{node1, node2} are required
#' to fit the model without explicitly specifying an adjacency matrix. This function is essentially
#' unchanged from Mitzi Morris's code \insertCite{Morris:2017,Morris:2019}{ssnet}.
#' @importFrom Rdpack reprompt
#' @references
#'
#' \insertRef{Morris:2017}{ssnet}
#'
#' \insertRef{Morris:2019}{ssnet}
#' @export
mungeCARdata4stan = function(adjBUGS,numBUGS) {
  J = length(numBUGS);
  nn = numBUGS;
  J_edges = length(adjBUGS) / 2;
  node1 = vector(mode="numeric", length=J_edges);
  node2 = vector(mode="numeric", length=J_edges);
  iAdj = 0;
  iEdge = 0;
  for (i in 1:J) {
    for (j in 1:nn[i]) {
      iAdj = iAdj + 1;
      if (i < adjBUGS[iAdj]) {
        iEdge = iEdge + 1;
        node1[iEdge] = i;
        node2[iEdge] = adjBUGS[iAdj];
      }
    }
  }
  return (list("J"=J,"J_edges"=J_edges,"node1"=node1,"node2"=node2));
}
