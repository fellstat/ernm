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
    t_1 <- proc.time()
    MRF <- ernm(samplike_undir ~ edges + homophily("group") + logisticNeighbors('group','group','Loyal') | group,
                tapered = FALSE,
                verbose = TRUE)
    t_1 <- proc.time() - t_1
    
    # Test ERGM verision of ERNM
    t_2 <- proc.time()
    ERGM <- ernm(samplike_undir ~ edges + gwesp(0.5) + gwdegree(0.5) + homophily("group") + logisticNeighbors('group','group','Loyal'),
                 tapered = FALSE,
                 verbose = FALSE)
    t_2 <- proc.time() - t_2
    
    # Test ERNM
    t_3 <- proc.time()[3]
    ERNM <- ernm(samplike_undir ~ edges + gwesp(0.5) + gwdegree(0.5) + homophily("group") + logisticNeighbors('group','group','Loyal') | group,
                 tapered = FALSE,
                verbose = FALSE)
    t_3 <- proc.time()[3] - t_3
    
    # Test tapered ERNM:
    ERNM_formula <- as.formula("samplike_undir ~edges + gwesp(0.5) + gwdegree(0.5) + homophily('group') + logisticNeighbors('group','group','Loyal') | group")
    stats <- ernm::calculateStatistics(ERNM_formula)
    t_4 <- proc.time()[3]
    ERNM_tapered_1 <- ernm(ERNM_formula,
                           tapered = TRUE,
                           modelArgs = list(tau = 1 / (3^2 * (stats + 5)),
                                            centers = stats,
                                            modelClass = 'TaperedModel'),
                           verbose = FALSE)
    t_4 <- proc.time()[3] - t_4
    
    # Test tapered ERNM:
    # more tapering needed here
    ERNM_formula <- as.formula("samplike_undir ~ edges + triangles() + star(2) + homophily('group') + logisticNeighbors('group','group','Loyal') | group")
    stats <- ernm::calculateStatistics(ERNM_formula)
    t_5 <- proc.time()[3]
    ERNM_tapered_2 <- ernm(ERNM_formula,
                           tapered = TRUE,
                           modelArgs = list(tau = 1 / (2 * (stats + 5)),
                                            centers = stats,
                                            modelClass = 'TaperedModel'),
                           verbose = FALSE)
    t_5 <- proc.time()[3] - t_5
    
    # Print the time it took to  fit the models vs expected time:
    print(as.data.frame(cbind(c("MRF", "ERGM", "ERNM", "ERNM_tapered_1", "ERNM_tapered_2"), c(t_1, t_2, t_3, t_4, t_5))))
    
    
    # All models should converge
    testthat::expect_true(ERGM$converged)
    testthat::expect_true(MRF$converged)
    testthat::expect_true(ERNM$converged)
    testthat::expect_true(ERNM_tapered_1$converged)
    testthat::expect_true(ERNM_tapered_2$converged)
}
)
    
    