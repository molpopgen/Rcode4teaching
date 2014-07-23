#' Rcode4teaching.
#'
#' @name Rcode4teaching
#' @docType package
#' @description Stuff I use to generate examples, etc., for lectures
#' @author Kevin R. Thornton \email{krthornt@@uci.edu}
#' @references \url{https://github.com/molpopgen/Rcode4teaching}
#' @exportPattern "^[[:alpha:]]+"
NULL

#' coalSample
#'
#' @name coalSample
#' @docType class
#' @description A class representing the output of a coalescent simulation
class('coalSample')

#' coalTree
#'
#' @name coalTree
#' @docType class
#' @description A bifurcating genealogy representing the history of a sample of size n
#' @details a coalTree is a list with the following elements:
#' nodes: An array of 2n-1 integers.  nodes[i] is the ancestor of the i-th node.  
#'	 The root is labeled with a -1
#' times: An array of 2n-1 times, corresponding to the values in nodes.
#' tmrca: The tmrca of the sample, equivalent to times[2n-1]
#' ttot:  The total time on the tree
class('coalTree')
