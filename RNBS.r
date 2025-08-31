library(igraph)
library(Matrix)
library(nnls)
library(parallel)

network_inf_kNN_glap <- function(network, gamma = 0.01, kn = 11, verbose = TRUE) {
    if (verbose) {
        glap_inv_starttime <- Sys.time()
    }

    network_nodes <- V(network)$name
    if (is.null(network_nodes)) {
        network_nodes <- as.character(1:vcount(network))
    }

    L_arr <- laplician_matrix(network, sparse = FALSE)

    L_vandin <- L_arr + gamma * diag(nrow(L_arr))

    L_inv_arr <- solve(L_vandin)
    rownames(L_inv_arr) <- network_nodes
    colnames(L_inv_arr) <- network_nodes    

    if (verbose) {
      cat("Graph influence matrix calculated:", 
          as.numeric(Sys.time() - glap_inv_starttime), "seconds\n")
      KNN_starttime <- Sys.time()
    }

    KNN_graph <- make_empty_graph(directed = FALSE)
    KNN_graph <- add_vertices(KNN_graph, length(network_nodes), name = network_nodes)

    for (gene in network_nodes) {
      gene_influences <- L_inv_arr[gene, ]
      gene_knn <- names(sort(gene_influences, decreasing = TRUE)[1:kn])

      for (neighbor in gene_knn) {
        if (L_inv_arr[gene, neighbor] > 0 && gene != neighbor) {
          KNN_graph <- add_edges(KNN_graph, c(gene, neighbor))
        }
      }
    }

    KNN_nodes <- V(KNN_graph)$name
    knnGlap <- laplacian_matrix(KNN_graph, sparse = FALSE)
    rownames(knnGlap) <- KNN_nodes
    colnames(knnGlap) <- KNN_nodes
    
    if (verbose) {
      cat("KNN Laplacian constructed:", 
          as.numeric(Sys.time() - KNN_starttime), "seconds\n")
    }
    
    return(knnGlap)
}