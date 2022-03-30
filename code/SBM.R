library(igraph)
library(tidygraph)
library(ggraph)
library(control)
library(netcontrol)
library(Matrix)
library(matrixcalc)
library(JuliaCall) 
julia_setup()
julia_library("ControlSystems")
julia_library("LinearAlgebra")
julia_library("SparseArrays")

MAXIMUM_GRAPH_SIXE <- 40 # Let's not get too ambitious

split_into_k <- function(n, k) {
    # Split {1,..., n} into k even-ish groups
    # split_into_k(100,3) -> 33, 33, 34
    first <- rep(round(n/k), k-1)
    remainder <- n - sum(first)
    return(c(first, n - sum(first)))
}

n <- 100
k <- 3
prob.within_group <- 0.25
prob.between_group <- 0.01

M <- matrix(0, k,k)
diag(M) <- rep(prob.within_group, k)
M[upper.tri(M)] <- prob.between_group
M[lower.tri(M)] <- prob.between_group


# SAMPLING FUNCTIONS
# We are going to sample graphs in the following way:
#   1. Sample a directed SBM 
#   2. for every connection, get the value of the weight by sampling 
#      it from a distribution
#   3. 
# NOTE: There is a paper "Adapting the Stochastic Block Model to Edge-Weighted Networks by Aicher et Al"
# which discusses estimation on these models, could be useful for later.

# Generic SBM Sampler
undirected_asymetric_sbm <- function(n, pref.matrix, block.sizes, sampler) {
    # Sampler is a function that generates the value of the weight.
    A <- sample_sbm(n, pref.matrix, block.sizes, directed=TRUE) %>% as_adj
    A[which(A==1)] <- sampler(length(which(A==1)))
    # TODO reweight
    colSums(A)
    # TODO add 1 if a certain condition holds
    # TODO Remove nodes that have no connections
    return(A)
}

create_sparse <- function(p, n, k) {
    # I think this is equivalent to Erdos-Renyi.
    # p: probability of any edge connecting
    # n: number of nodes
    # k: number of clusters
    sbm <- matrix(p, k, k)
    A <- undirected_asymetric_sbm(n, sbm, split_into_k(n, k), sampler=runif)
    return(A)
}

create_one_big_cluster_graph( ) {
    ...
}

create_clusters <- function(n, k) {
    # TODO figure out how we want to do the clusters
    # p:
    # k: number of clusters
    # n: number of nodes 
    ...
}

graph_ <- create_sparse(0.05,n,k)
graph_ %>% ggraph() +
    geom_edge_link(color="grey") +
    geom_node_point(size=5, color="firebrick")


minimum_driver_nodes <- function(A) {
    """
    Only works with symmetric adjacency matries
    """
    counter <- table(eigen(A)$values)
    return(max(t))
}

# node_level_metrics <- function(A, node_index) {
#     n <- dim(A)[1]
#     e <- rep(0, n )
#     e[node_index] <- 1
#     Ci <- ctrb(A, e) # Controllability matrix for node
#     control_centrality <- rankMatrix(Ci)[[1]]
#     Wi <- Ci %*% t(Ci) # TODO use the network control code to calculate this (doens't work directly)
#                        # Would just be a matter of copy and pasting code I think

#     node_to_network_centrality <- sum(diag(Wi))

#     network_to_node_centrality <- NA # TODO
    
#     return(list(
#                 control_centrality=control_centrality,
                
#                 ))
# }

metrics <- function(graph) {
    A <- as_adjacency_matrix(graph)
    # n_D <- minimum_driver_nodes(A) <- would need function to find maximum of GEOMETRIC multiplicity

}


