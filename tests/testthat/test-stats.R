context("stat tests")

library(testthat)
library(ergm)
library(network)

test_that("Stats", {
    
    # ================
    # Setup
    # ================
    set.seed(1)
    # make 1 100 node network with some variables:
    add_treated_neighs <- function(net,treatment_var){
        tmp <- as.numeric(get.vertex.attribute(net,treatment_var))
        tmp <- tmp == max(tmp)
        set.vertex.attribute(net,paste(treatment_var,"_neighbors",sep=""),as.character(sapply(1:(net%n%'n'),function(i){sum(tmp[get.neighborhood(net,i)])})))
        return(net)
    }
    
    make_net = function(...){
        tmp <- matrix(rnorm(10000)>2,nrow = 100)
        net <- as.network(tmp,directed =F)
        set.vertex.attribute(net,"var_1",as.character(apply(rmultinom(100,1,rep(1/3,3)),2,FUN = function(x){which(x==1)})))
        set.vertex.attribute(net,"var_2",as.character(apply(rmultinom(100,1,rep(1/2,2)),2,FUN = function(x){which(x==1)})))
        set.vertex.attribute(net,"var_3",as.character(apply(rmultinom(100,1,rep(1/2,2)),2,FUN = function(x){which(x==1)})))
        summary(as.factor(get.vertex.attribute(net,"var_1")))
        
        net <- add_treated_neighs(net,"var_1")
        net <- add_treated_neighs(net,"var_2")
        net <- add_treated_neighs(net,"var_3")
    }
    
    net = make_net()
    nets = lapply(1:100,FUN = make_net)
    
    
    # ================
    # Hamming Distance
    # ================
    
    edges = as.edgelist(net)*1
    edges = cbind(as.double(edges[,1]),as.double(edges[,2]))
    edges = cbind(edges[,1],edges[,2])
    head(edges)
    
    hamming_calc = function(edges,net){
        
        e_list = as.edgelist(net)*1.0
        e_list = cbind(as.double(e_list[,1]),as.double(e_list[,2]))
        
        tmp = rbind(e_list,edges)
        shared = sum(duplicated(tmp))
        
        dist = (dim(e_list)[1] - shared) + (dim(edges)[1] - shared)
        return(dist)
    }
    
    hamming_test_1 <- (
        ernm::calculateStatistics(net ~ hamming(edges,100)) == hamming_calc(edges,net)
    )
    
    edges[1,] = c(1,11)
    hamming_test_2 <- (
        ernm::calculateStatistics(net ~ hamming(edges,100)) == hamming_calc(edges,net)
    )
    
    edges[2,] = c(1,13)
    hamming_test_3 <- (
        ernm::calculateStatistics(net ~ hamming(edges,100)) == hamming_calc(edges,net)
    )
    
    
    # check dyad update function
    edges = as.edgelist(net)*1
    edges = cbind(as.double(edges[,1]),as.double(edges[,2]))
    
    # update included edge that network has
    model = ernm(net ~ hamming(edges,100),maxIter = 2,mcmcBurnIn = 100,mcmcInterval = 10, mcmcSampleSize = 100)
    model = model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    model$calculate()
    
    hamming_calc(edges,net)
    net_tmp = net
    delete.edges(net_tmp,get.edgeIDs(net_tmp,1,61))
    r_stat_1 = hamming_calc(edges,net_tmp)
    
    model$statistics()
    model$dyadUpdate(1,61)
    cpp_stat_1 =  model$statistics()
    
    # update included edge that network has
    model = ernm(net ~ hamming(edges,100),maxIter = 2,mcmcBurnIn = 100,mcmcInterval = 10, mcmcSampleSize = 100)
    model = model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    model$calculate()
    
    hamming_calc(edges,net)
    net_tmp = net
    delete.edges(net_tmp,get.edgeIDs(net_tmp,61,1))
    r_stat_2 = hamming_calc(edges,net_tmp)
    
    model$statistics()
    model$dyadUpdate(61,1)
    cpp_stat_2 =  model$statistics()
    
    # update included edge that network does not have
    edges = matrix(c(1,2),nrow = 1)
    model = ernm(net ~ hamming(edges,100),maxIter = 2,mcmcBurnIn = 100,mcmcInterval = 10, mcmcSampleSize = 100)
    model = model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    model$calculate()
    model$statistics()
    
    hamming_calc(edges,net)
    net_tmp = net
    add.edges(net_tmp,1,2)
    r_stat_3 = hamming_calc(edges,net_tmp)
    
    model$statistics()
    model$dyadUpdate(1,2)
    cpp_stat_3 = model$statistics()
    
    # update not included edge that network has
    edges = matrix(c(1,2),nrow = 1)
    model = ernm(net ~ hamming(edges,100),maxIter = 2,mcmcBurnIn = 100,mcmcInterval = 10, mcmcSampleSize = 100)
    model = model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    
    hamming_calc(edges,net)
    net_tmp = net
    delete.edges(net_tmp,get.edgeIDs(net_tmp,1,61))
    r_stat_4 = hamming_calc(edges,net_tmp)
    
    model$statistics()
    model$dyadUpdate(1,61)
    cpp_stat_4 = model$statistics()
    
    # update not included edge that network does not have
    edges = as.edgelist(net)*1.0
    edges = cbind(as.double(edges[,1]),as.double(edges[,2]))
    head(edges)
    
    model = ernm(net ~ hamming(edges,100),maxIter = 2,mcmcBurnIn = 100,mcmcInterval = 10, mcmcSampleSize = 100)
    model = model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    model$calculate()

    hamming_calc(edges,net)
    net_tmp = net
    add.edges(net_tmp,1,2)
    r_stat_5 = hamming_calc(edges,net_tmp)
    
    model$statistics()
    model$dyadUpdate(1,2)
    cpp_stat_5 = model$statistics()
    
    hamming_test_4 = (r_stat_1 == cpp_stat_1)
    hamming_test_5 = (r_stat_2 == cpp_stat_2)
    hamming_test_6 = (r_stat_3 == cpp_stat_3)
    hamming_test_7 = (r_stat_4 == cpp_stat_4)
    hamming_test_7 = (r_stat_5 == cpp_stat_5)

    testthat::expect_true(hamming_test_1)
    testthat::expect_true(hamming_test_2)
    testthat::expect_true(hamming_test_3)
    testthat::expect_true(hamming_test_4)
    testthat::expect_true(hamming_test_5)
    testthat::expect_true(hamming_test_6)
    testthat::expect_true(hamming_test_7)
}
)