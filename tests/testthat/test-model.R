# Tests for model fitting that may reveal problems with the some part of the fitting procedure
context("model tests")

library(testthat)
library(network)
library(ernm)
data("samplike")

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
    MRF <- ernm(samplike_undir ~ edges + homophily("group") + logisticNeighbors('group','group','Loyal') | group,
                tapered = FALSE,
                verbose = FALSE)
    
    # Test ERGM verision of ERNM
    ERGM <- ernm(samplike_undir ~ edges + gwesp(0.5) + gwdegree(0.5) + homophily("group") + logisticNeighbors('group','group','Loyal'),
                 tapered = FALSE,
                 verbose = FALSE)
    
    # Test ERNM
    t_1 <- proc.time()[3]
    ERNM <- ernm(samplike_undir ~ edges + gwesp(0.5) + gwdegree(0.5) + homophily("group") + logisticNeighbors('group','group','Loyal') | group,
                 tapered = FALSE,
                 verbose = FALSE)
    t_1 <- proc.time()[3] - t_1
    
    # Test tapered ERNM:
    ERNM_formula <- as.formula("samplike_undir ~edges + gwesp(0.5) + gwdegree(0.5) + homophily('group') + logisticNeighbors('group','group','Loyal') | group")
    stats <- ernm::calculateStatistics(ERNM_formula)
    t_2 <- proc.time()[3]
    ERNM_tapered_1 <- ernm(ERNM_formula,
                           tapered = TRUE,
                           modelArgs = list(tau = 1 / (3^2 * (stats + 5)),
                                            centers = stats,
                                            modelClass = 'TaperedModel'),
                           verbose = FALSE)
    t_2 <- proc.time()[3] - t_2
    
    # Test tapered ERNM:
    # more tapering needed here
    ERNM_formula <- as.formula("samplike_undir ~ edges + triangles() + star(2) + homophily('group') + logisticNeighbors('group','group','Loyal') | group")
    stats <- ernm::calculateStatistics(ERNM_formula)
    ERNM_tapered_2 <- ernm(ERNM_formula,
                           tapered = TRUE,
                           modelArgs = list(tau = 1 / (2 * (stats + 5)),
                                            centers = stats,
                                            modelClass = 'TaperedModel'),
                           verbose = FALSE)
    
    # Test bulk dyad updates vs single dyad updates
    tails <- sample(1:network.size(samplike_undir), 100, replace = TRUE)
    heads <- sample(1:network.size(samplike_undir), 100, replace = TRUE)
    # make tails less than heads:
    old_tails <- tails
    tails <- pmin(tails, heads)
    heads <- pmax(old_tails, heads)
    drop <- which(tails == heads)
    if(length(drop) > 0){
      tails <- tails[-drop]
      heads <- heads[-drop]
    }
    ERNM <- ernm(samplike_undir ~ edges + gwesp(0.5) + gwdegree(0.5) | group,
                 tapered = FALSE,
                 verbose = FALSE)
    model <- ERNM$m$sampler$getModel()
    t_3 <- proc.time()[3]
    change_stats_1 <- mapply(tails,heads,FUN = function(tail, head) {
      old <- model$statistics()
      model$dyadUpdate(tail,head)
      new <- model$statistics()
      model$dyadUpdate(tail,head)
      return(new-old)
    },SIMPLIFY = F)
    t_3 <- proc.time()[3] - t_3
    
    
    t_4 <- proc.time()[3]
    change_stats_2 <- model$computeChangeStats(tails, heads)
    t_4 <- proc.time()[3] - t_4
    
    r_t <- microbenchmark::microbenchmark( 
      change_stats_1 <- mapply(tails,heads,FUN = function(tail, head) {
        old <- model$statistics()
        model$dyadUpdate(tail,head)
        new <- model$statistics()
        model$dyadUpdate(tail,head)
        return(new-old)
    },SIMPLIFY = F),times = 100)
    r_t
    
    cpp_t <- microbenchmark::microbenchmark(change_stats_2 <- model$computeChangeStats(tails, heads),
                                            times = 100)
    cpp_t
    
    # do it repeadtly:
    tmp <- model$computeChangeStats(rep(tails,each = 5), rep(heads,each =5))
    change_stats_rep <- mapply(rep(tails,each=5),rep(heads,each=5),FUN = function(tail, head) {
      old <- model$statistics()
      model$dyadUpdate(tail,head)
      new <- model$statistics()
      model$dyadUpdate(tail,head)
      return(new-old)
    },SIMPLIFY = F)
    change_stats_rep <- do.call(rbind,change_stats_rep)
    cbind(tmp,change_stats_rep)
    # GWESP SEEMS UNSTABLE but maybe not wrong?
      
    # test that all rows are the same:
    bulk_change_stat_test <- sapply(1:length(tails),function(i){
      all(change_stats_1[[i]] == change_stats_2[i,])
    })
    bulk_change_stat_test_1 <- all(bulk_change_stat_test)
    bulk_change_stat_test_2 <- t_4<t_3
    
    
    
    # All models should converge
    testthat::expect_true(ERGM$converged)
    testthat::expect_true(MRF$converged)
    testthat::expect_true(ERNM$converged)
    testthat::expect_true(ERNM_tapered_1$converged)
    testthat::expect_true(ERNM_tapered_2$converged)
    testthat::expect_true(bulk_change_stat_test)
    testthat::expect_true(bulk_change_stat_test_1)
    testthat::expect_true(bulk_change_stat_test_2)
}
)
    
    