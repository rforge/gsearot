addSignMatrix <-
function( graph ){
  #
  #   Adding the sign matrix to the graph object. A sign matrix is
  #   a matrix that indicate the action between nodes (+1 or -1, 
  #   corresponding to activation and inhibition) in a network. A 0 means
  #   the two corresponding nodes are independent (not connected or
  #   activations/inhibitions cancel each other).
  #
  #   Originally by Lars Snipen
  #   Biostatistics, Norwegian University of Life Sciences
  require( igraph )
  
  inhib.idx <- grep( "inhibition", get.edge.attribute( graph, "Action" ) )
  
  sign.mat <- matrix( rep( 1, vcount( graph )^2 ), ncol=vcount( graph ) )
  if( length( inhib.idx ) > 0 ){
    for( i in 1:length( inhib.idx ) ){
      ixx <- get.edge( graph, (inhib.idx[i]-1) )
      sign.mat[ixx[1]+1,ixx[2]+1] <- -1
      sign.mat[ixx[2]+1,ixx[1]+1] <- -1
    }
  }
  sign.mat <- .signProp( graph, get.adjacency( graph ) * sign.mat )
  graph <- set.graph.attribute( graph, "Sign.matrix", sign.mat )
  return( graph )
}
