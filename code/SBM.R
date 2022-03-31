library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(control)
library(netcontrol)
library(Matrix)
library(matrixcalc)
library(JuliaCall) 
julia_setup()
julia_source("metrics.jl")

MAXIMUM_GRAPH_SIZE <- 50 # Let's not get too ambitious

## UTILITY FUNCTIONs
split_into_k <- function(n, k) {
    # Split {1,..., n} into k even-ish groups
    # split_into_k(100,3) -> 33, 33, 34
    first <- rep(round(n/k), k-1)
    remainder <- n - sum(first)
    return(c(first, n - sum(first)))
}

# SAMPLING FUNCTIONS
# We are going to sample graphs in the following way:
#   1. Sample a directed SBM 
#   2. for every connection, get the value of the weight by sampling 
#      it from a distribution
#   3. 
# NOTE: There is a paper "Adapting the Stochastic Block Model to Edge-Weighted Networks by Aicher et Al"
# which discusses estimation on these models, could be useful for later.

# Generic SBM Sampler
undirected_asymetric_sbm <- function(n, pref.matrix, block.sizes, sampler=rlnorm) {
    # Sampler is a function that generates the value of the weight.
    graph <- sample_sbm(n, pref.matrix, block.sizes, directed=TRUE)
    A <- as_adjacency_matrix(graph)
    A[which(A==1)] <- sampler(length(which(A==1))) # Sample from distribution for those edges that are connected.
    A = matrix(A, dim(A)[1], dim(A)[1]) # Convert to non-sparse
    return(as_tbl_graph(A))
}

# preferential_attachment TODO

create_ER <- function(n, p) {
    # I think this is equivalent to Erdos-Renyi.
    # p: probability of any edge connecting
    # n: number of nodes
    # k: number of clusters
    sbm <- matrix(p, 1, 1)
    graph <- undirected_asymetric_sbm(n, sbm, split_into_k(n, 1))
    return(graph)
}

create_one_big_cluster_graph <- function(n, ...) {
    cluster_self_prob <- 0.7
    # from big cluster out
    cluster_out_prob <- 0.2
    # from big cluster to in
    cluster_in_prob <- 0.15
    non_cluster_between <- 0.05
    
    sbm <- matrix( c( 
                     cluster_self_prob, cluster_out_prob, 
                     cluster_in_prob, non_cluster_between), 
                  nrow=2, ncol=2 )
    undirected_asymetric_sbm(n, sbm, c(15, 35), runif)
}

create_two_big_cluster_graph <- function(n, ...) {
    first_cluster_self_prob <- 0.75
    # from big cluster out
    cluster_between_prob <- 0.1
    # from big cluster to in
    cluster_between_prob <- 0.1
    second_cluster_between <- 0.8
    
    sbm <- matrix( c( 
                     first_cluster_self_prob, cluster_between_prob, 
                     cluster_between_prob, second_cluster_between), 
                  nrow=2, ncol=2 )
    undirected_asymetric_sbm(n, sbm, c(25, 25), runif)
}

# TESTING
graph <- create_two_big_cluster_graph(MAXIMUM_GRAPH_SIZE, sampler=runif)
graph <- graph %>% activate(nodes) %>% mutate(id = row_number()) %>% mutate(cluster = ifelse( id <= 25, "dense", "sparse"))
ggraph(graph) + 
    geom_edge_link(color="black", arrow= arrow(length = unit(4, 'mm'))) + 
    geom_node_point(aes(color=factor(cluster)), size=4) +
    theme_void()
# geom_node_text(aes_string(label=label), size=7, col="black") +
# scale_x_reverse() +
# ggtitle(message) +
# labs(color="groups", title=message) +
# theme_void() +
# theme(plot.title = element_text(lineheight=1.5, size=20, face="bold.italic"),
#       legend.position = "none"
#       ) +
# scale_color_tableau("Superfishel Stone")

create_clusters <- function(n, k) {
    # TODO figure out how we want to do the clusters
    # p:
    # k: number of clusters
    # n: number of nodes 
    ...
}

regularize_graph <- function(graph) {
    n <- gorder(graph)
    A <- matrix(as_adjacency_matrix(graph), n, n)
    A <- julia_call("regularize", A)
    as_tbl_graph(A)
}

calc_control_metrics <- function(graph) {
    A <- matrix(as_adjacency_matrix(graph), gorder(graph), gorder(graph))
    control_metrics <- julia_call("get_metrics", A)
    return(control_metrics)
}


simulate <- function() {
    result <- NA
    result_df <- NA
    N <- MAXIMUM_GRAPH_SIZE

    # ATTENTION: every function must take two arguments
    graph_simulators <- list( 
                             "ER" = function(n, p) create_ER(n, p),
                             "1C" = function(n, sampler) create_one_big_cluster_graph(n, runif),
                             "2C" = function(n, sampler) create_two_big_cluster_graph(n, runif)
                             )

    param_grid <- list(
                       "ER"=seq(0, 1, 0.05)[-1],
                       "1C"=list(0), # ignore
                       "2C"=list(0) # ignore
                     )

    simulation_number <- list(
                        "ER"=1,
                        "1C"=50,
                        "2C"=50
                     )

    i <- 1
    for (name in names(graph_simulators)) {
        for (p in param_grid[[name]]) {
            for (k in 1:simulation_number[[name]]) {
                graph <- graph_simulators[[name]](N, p) 
                graph <- regularize_graph(graph)
                control_metrics <- calc_control_metrics(graph)
                graph_index <- i
                result <- tibble(
                    graph_index=i,
                    graph_simulator=paste(name, p, sep="_"),
                    node_index=1:gorder(graph),
                    eigenvector_centrality=eigen_centrality(graph)$vector,
                    meandegree=mean(degree(graph)),
                    edges=gsize(graph),
                    density=graph.density(graph),
                    betweeness=betweenness(graph),
                    node_to_network=control_metrics[[1]], # tr(W)
                    network_to_node=control_metrics[[2]], # tr(M)
                    ratio = node_to_network / network_to_node
                )
                print(name)
                if (i == 1) {
                    result_df <- result
                }
                else {
                    result_df <- bind_rows(result_df, result)
                }
                i <- i + 1
            }
        }
    }
    return(result_df)
}

result_df <- simulate()

save(result_df, file="data.RData")
