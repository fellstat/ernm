# Tests for model fitting that may reveal problems with the some part of the fitting procedure
context("model tests")

library(testthat)
library(ergm)
library(network)
library(ernm)
data(sampson)

test_that("models", {
    
    # make undirected for ease:
    adj_matrix <- as.matrix(samplike, matrix.type = "adjacency")
    undirected_adj_matrix <- adj_matrix | t(adj_matrix)
    undirected_net <- network(undirected_adj_matrix, directed = FALSE)
    set.vertex.attribute(undirected_net, "cloisterville ", samplike %v% 'cloisterville ')
    set.vertex.attribute(undirected_net, "group", samplike %v% 'group')
    samplike_undir <- undirected_net
    
    # Display the undirected network
    samplike <- as.network(samplike_undir, directed = FALSE)
    
    # Test MRF version of ERNM
    MRF <- ernm(samplike_undir ~ edges + homophily("group") + logisticNeighbors('group','group','Loyal') | group,verbose = FALSE)
    
    # Test ERGM verision of ERNM
    ERGM <- ernm(samplike_undir ~ edges + gwesp(0.5) + gwdegree(0.5) + homophily("group") + logisticNeighbors('group','group','Loyal'),verbose = FALSE)
    
    # Test ERNM
    t_1 <- proc.time()[3]
    ERNM <- ernm(samplike_undir ~ edges + gwesp(0.5) + gwdegree(0.5) + homophily("group") + logisticNeighbors('group','group','Loyal') | group,verbose = FALSE)
    t_1 <- proc.time()[3] - t_1
    
    # Test tapered ERNM:
    ERNM_formula <- as.formula("samplike_undir ~edges + gwesp(0.5) + gwdegree(0.5) + homophily('group') + logisticNeighbors('group','group','Loyal') | group")
    stats <- ernm::calculateStatistics(ERNM_formula)
    t_2 <- proc.time()[3]
    ERNM_tapered_1 <- ernm(ERNM_formula,
                           modelArgs = list(tau = 1/(2*stats),
                                            centers = stats,
                                            modelClass = 'TaperedModel'),
                           verbose = TRUE)
    t_2 <- proc.time()[3] - t_2
    
    # Test tapered ERNM:
    # more tapering needed here
    ERNM_formula <- as.formula("samplike_undir ~ edges + triangles() + star(2) + homophily('group') + logisticNeighbors('group','group','Loyal') | group")
    stats <- ernm::calculateStatistics(ERNM_formula)
    ERNM_tapered_2 <- ernm(ERNM_formula,
                           modelArgs = list(modelClass = 'TaperedModel',
                                            tau = 1/(2*stats),
                                            centers = stats),
                           verbose = FALSE)
    
    
    # ERNM taper should be quicker due to not having to sample from bad theta suggestions:
    testthat::expect_true(t_2 < t_1)
    # All models should converge
    testthat::expect_true(ERGM$converged)
    testthat::expect_true(MRF$converged)
    testthat::expect_true(ERNM$converged)
    testthat::expect_true(ERNM_tapered_1$converged)
    testthat::expect_true(ERNM_tapered_2$converged)
}
)
    
    