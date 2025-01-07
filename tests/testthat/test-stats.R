context("stat tests")

library(testthat)
library(ergm)
library(network)
library(ernm)

test_that("Stats", {
    
    # ================
    # Setup
    # ================
    set.seed(1)
    N = 100
    # make 1 100 node network with some variables:
    add_treated_neighs <- function(net,treatment_var){
        tmp <- as.numeric(get.vertex.attribute(net,treatment_var))
        tmp <- tmp == max(tmp)
        set.vertex.attribute(net,paste(treatment_var,"_neighbors",sep=""),as.character(sapply(1:(net%n%'n'),function(i){sum(tmp[get.neighborhood(net,i)])})))
        return(net)
    }
    
    make_net <- function(...){
        tmp <- matrix(rnorm(10000)>2,nrow = N)
        net <- as.network(tmp,directed =F)
        set.vertex.attribute(net,"var_1",as.character(apply(rmultinom(N,1,rep(1/3,3)),2,FUN = function(x){which(x==1)})))
        set.vertex.attribute(net,"var_2",as.character(apply(rmultinom(N,1,rep(1/2,2)),2,FUN = function(x){which(x==1)})))
        set.vertex.attribute(net,"var_3",as.character(apply(rmultinom(N,1,rep(1/2,2)),2,FUN = function(x){which(x==1)})))
        set.vertex.attribute(net,"var_4",rnorm(N))
        set.vertex.attribute(net,"var_5",rnorm(N))
        summary(as.factor(get.vertex.attribute(net,"var_1")))
        
        net <- add_treated_neighs(net,"var_1")
        net <- add_treated_neighs(net,"var_2")
        net <- add_treated_neighs(net,"var_3")
    }
    
    net<-make_net()
    nets<-lapply(1:100,FUN = make_net)
    
    # ================
    # Node Count
    # ================
    v1 = (net %v% "var_2")
    v2 = (net %v% "var_3")
    
    r_stat_1 = sum((v1=="2"))
    cpp_stat_1 = as.numeric(ernm::calculateStatistics(net ~ nodeCount("var_2")))
    
    model <- ernm(net ~ nodeCount("var_2","1")| var_2 ,
                  tapered = FALSE,
                  maxIter = 2,
                  mcmcBurnIn = 100,
                  mcmcInterval = 10,
                  mcmcSampleSize = 100,
                  verbose = FALSE)
    model <- model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    model$calculate()
    model$statistics()
    
    # Flip var 2 from 0 to 1 when var 3 is 0
    y_pos = which(net %v% "var_2"== "2")
    vert = y_pos[1]
    model$discreteVertexUpdate(vert,"var_2",1)
    cpp_stat_2 <- model$statistics()
    r_stat_2 <- r_stat_1 - 1
    #reset:
    model$calculate()
    
    y_neg = which(net %v% "var_2"== "1")
    vert = y_neg[1]
    model$discreteVertexUpdate(vert, "var_2",2)
    cpp_stat_3 <- model$statistics()
    r_stat_3 <- r_stat_1 + 1
    #reset:
    model$calculate()
    
    # check a dyad update makes not difference:
    r_stat_4 = r_stat_1
    model$dyadUpdate(61,1)
    cpp_stat_4 <- model$statistics()
    #reset:
    model$calculate()
    
    # check a dyad update makes not difference:
    r_stat_5 = r_stat_1
    model$dyadUpdate(1,61)
    cpp_stat_5 <- model$statistics()
    #reset:
    model$calculate()
    
    nodeCount_test_1 <- (r_stat_1 == cpp_stat_1)
    nodeCount_test_2 <- (r_stat_2 == cpp_stat_2)
    nodeCount_test_3 <- (r_stat_3 == cpp_stat_3)
    nodeCount_test_4 <- (r_stat_4 == cpp_stat_4)
    nodeCount_test_5 <- (r_stat_5 == cpp_stat_5)
    
    testthat::expect_true(nodeCount_test_1)
    testthat::expect_true(nodeCount_test_2)
    testthat::expect_true(nodeCount_test_3)
    testthat::expect_true(nodeCount_test_4)
    testthat::expect_true(nodeCount_test_5)
        
    
    # ================
    # Logistic Regression
    # ================
    
    v1 = (net %v% "var_2")
    v2 = (net %v% "var_3")
    
    r_stat_1 = sum((v1=="2")[v2=="2"])
    r_stat_2 = r_stat_1
    r_stat_3 = sum((v1=="2")[v2=="1"])
    
    cpp_stat_1 = as.numeric(ernm::calculateStatistics(net ~ logistic("var_2","var_3","1") | var_2))
    cpp_stat_2 = as.numeric(ernm::calculateStatistics(net ~ logistic("var_2","var_3") | var_2))
    cpp_stat_3 = as.numeric(ernm::calculateStatistics(net ~ logistic("var_2","var_3","2") | var_2))

    # Check the discreteVarUpdate function works

    y_pos = which(net %v% "var_2"== "2")
    y_neg = which(net %v% "var_2"== "1")
    x_pos = which(net %v% "var_3"== "2")
    x_neg = which(net %v% "var_3"== "1")
    
    model <- ernm(net ~ logistic("var_2","var_3","1") | var_2 ,
                  tapered = FALSE,
                  maxIter = 2,
                  mcmcBurnIn = 100,
                  mcmcInterval = 10,
                  mcmcSampleSize = 100,
                  verbose = FALSE)
    model <- model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    model$calculate()
    model$statistics()
    
    # Note updates take in the factor level that the 
    
    # Flip var 2 from 0 to 1 when var 3 is 0
    vert = intersect(y_neg,x_neg)[1]
    model$discreteVertexUpdate(vert, "var_2",2)
    cpp_stat_4 <- model$statistics()
    r_stat_4 <- sum((v1=="2")[v2=="2"])
    #reset:
    model$calculate()
    
    # Flip var 2 from 0 to 1 when var 3 is 1
    vert = intersect(y_neg,x_pos)[1]
    model$discreteVertexUpdate(vert, "var_2",2)
    cpp_stat_5 <- model$statistics()
    r_stat_5 <- sum((v1=="2")[v2=="2"]) + 1
    #reset:
    model$calculate()
    
    # Flip var 2 from 1 to 0 when var 3 is 0
    vert = intersect(y_pos,x_neg)[1]
    model$discreteVertexUpdate(vert, "var_2",1)
    cpp_stat_6 <- model$statistics()
    r_stat_6 <- sum((v1=="2")[v2=="2"])
    #reset:
    model$calculate()
    
    # Flip var 2 from 1 to 0 when var 3 is 1
    vert = intersect(y_pos,x_pos)[1]
    model$discreteVertexUpdate(vert, "var_2",1)
    cpp_stat_7 <- model$statistics()
    r_stat_7 <- sum((v1=="2")[v2=="2"]) -1
    #reset:
    model$calculate()
    
    logistic_test_1 <- (r_stat_1 == cpp_stat_1)
    logistic_test_2 <- (r_stat_2 == cpp_stat_2)
    logistic_test_3 <- (r_stat_3 == cpp_stat_3)
    logistic_test_4 <- (r_stat_4 == cpp_stat_4)
    logistic_test_5 <- (r_stat_5 == cpp_stat_5)
    logistic_test_6 <- (r_stat_6 == cpp_stat_6)
    logistic_test_7 <- (r_stat_7 == cpp_stat_7)
    
    testthat::expect_true(logistic_test_1)
    testthat::expect_true(logistic_test_2)
    testthat::expect_true(logistic_test_3)
    testthat::expect_true(logistic_test_4)
    testthat::expect_true(logistic_test_5)
    testthat::expect_true(logistic_test_6)
    testthat::expect_true(logistic_test_7)
    
    # ================
    # Logistic Neighbors
    # ================
    
    # Calculate for base level = "2"
    r_stat_1 <- sum(as.numeric(net %v% "var_2_neighbors")[which(net %v% "var_2" == "1")])
    cpp_stat_1 <- ernm::calculateStatistics(net ~ logisticNeighbors("var_2","var_2","2") | var_2)
    
    # Calculate for base level = "1"
    r_stat_2 <- sum(as.numeric(net %v% "var_2_neighbors")[which(net %v% "var_2" == "2")])
    cpp_stat_2 <- ernm::calculateStatistics(net ~ logisticNeighbors("var_2","var_2","1") | var_2)
    
    # Test over all the nets to be sure:
    r_stat_3 <- TRUE
    cpp_stat_3 <- TRUE
    for(i in 1:length(nets)){
        method_1 <- sum(as.numeric(nets[[i]] %v% "var_2_neighbors")[which(nets[[i]] %v% "var_2" == "1")])
        method_2 <- as.numeric(ernm::calculateStatistics(nets[[i]] ~ logisticNeighbors("var_2","var_2","2") | var_2))
        if(method_1 != method_2){
            print(paste("failed on net ",i))
            r_stat_3 = FALSE
        }
    }
    
    # test the dyadUpdate function :
    which(net %v% "var_2"== "2")
    net_2 <- net
    tmp_model <- ernm(net ~ logisticNeighbors("var_2","var_2","2") | var_2,
                      tapered = FALSE,
                      maxIter = 2,
                      mcmcSampleSize = 1000,
                      mcmcBurnIn = 100,
                      verbose = FALSE)
    r_stat_4 <- TRUE
    cpp_stat_4 <- TRUE
    for(i in which(net %v% "var_2"== "2")[-1]){
        net_2 <- net
        net_2[i,2] <- 1 - net_2[i,2]
        
        change_1 <- ernm::calculateStatistics(net_2 ~ logisticNeighbors("var_2","var_2","2") | var_2) -
            ernm::calculateStatistics(net ~ logisticNeighbors("var_2","var_2","2") | var_2)
        
        model <- tmp_model
        
        model <- model$m$sample$getModel()
        model$setNetwork(as.BinaryNet(net))
        model$calculate()
        
        stat1 <- model$statistics()
        model$dyadUpdate(2,i)
        stat2 <- model$statistics()
        change_2 <- stat2-stat1
        
        if(change_1 != change_2){
            r_stat_4 <- FALSE
        }
    }
    
    # test for changing a node value to
    r_stat_5 <- TRUE
    cpp_stat_5 <- TRUE
    for(i in 1:(net%n% "n")){
        net_2 <- net
        if((net_2 %v% "var_2")[i] == "1"){
            new_value <- "2"
        }else{
            new_value <- "1"
        }
        set.vertex.attribute(net_2,"var_2",new_value,i)
        
        change_1 <- ernm::calculateStatistics(net_2 ~ logisticNeighbors("var_2","var_2","2") | var_2) - ernm::calculateStatistics(net ~ logisticNeighbors("var_2","var_2","2") | var_2)
        
        model <- tmp_model
        
        model <- model$m$sample$getModel()
        model$setNetwork(as.BinaryNet(net))
        model$calculate()
        
        stat1 <- model$statistics()
        model$discreteVertexUpdate(i,"var_2",as.numeric(new_value))
        stat2 <- model$statistics()
        change_2 <- stat2-stat1
        
        if(change_1 != change_2){
            r_stat_5 <- FALSE
        }
    }
    
    logistic_neighbot_test_1 <- (r_stat_1 == cpp_stat_1)
    logistic_neighbot_test_2 <- (r_stat_2 == cpp_stat_2)
    logistic_neighbot_test_3 <- (r_stat_3 == cpp_stat_3)
    logistic_neighbot_test_4 <- (r_stat_4 == cpp_stat_4)
    logistic_neighbot_test_5 <- (r_stat_5 == cpp_stat_5)
    
    testthat::expect_true(logistic_neighbot_test_1)
    testthat::expect_true(logistic_neighbot_test_2)
    testthat::expect_true(logistic_neighbot_test_3)
    testthat::expect_true(logistic_neighbot_test_4)
    testthat::expect_true(logistic_neighbot_test_5)

    # ================
    # Hamming Distance
    # ================
    
    edges <- as.edgelist(net)*1
    edges <- cbind(as.double(edges[,1]),as.double(edges[,2]))
    edges <- cbind(edges[,1],edges[,2])
    head(edges)
    
    hamming_calc <- function(edges,net){
        
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
    
    edges[1,] <-  c(1,11)
    hamming_test_2 <- (
        ernm::calculateStatistics(net ~ hamming(edges,100)) == hamming_calc(edges,net)
    )
    
    edges[2,] <-  c(1,13)
    hamming_test_3 <- (
        ernm::calculateStatistics(net ~ hamming(edges,100)) == hamming_calc(edges,net)
    )
    
    
    # check dyad update function
    edges <- as.edgelist(net)*1
    edges <- cbind(as.double(edges[,1]),as.double(edges[,2]))
    
    # update included edge that network has
    model <- ernm(net ~ hamming(edges,100),
                  tapered = FALSE,
                  maxIter = 2,
                  mcmcBurnIn = 100,
                  mcmcInterval = 10,
                  mcmcSampleSize = 100,
                  verbose = FALSE)
    model <- model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    model$calculate()
    
    hamming_calc(edges,net)
    net_tmp <- net
    delete.edges(net_tmp,get.edgeIDs(net_tmp,1,61))
    r_stat_1 <- hamming_calc(edges,net_tmp)
    
    model$statistics()
    model$dyadUpdate(1,61)
    cpp_stat_1 <- model$statistics()
    
    # update included edge that network has
    model <- ernm(net ~ hamming(edges,100),
                  tapered = FALSE,
                  maxIter = 2,
                  mcmcBurnIn = 100,
                  mcmcInterval = 10,
                  mcmcSampleSize = 100,
                  verbose = FALSE)
    model <- model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    model$calculate()
    
    hamming_calc(edges,net)
    net_tmp <- net
    delete.edges(net_tmp,get.edgeIDs(net_tmp,61,1))
    r_stat_2 <- hamming_calc(edges,net_tmp)
    
    model$dyadUpdate(61,1)
    cpp_stat_2 <- model$statistics()
    
    # update included edge that network does not have
    edges <- matrix(c(1,2),nrow = 1)
    model <- ernm(net ~ hamming(edges,100),
                  tapered = FALSE,
                  maxIter = 2,
                  mcmcBurnIn = 100,
                  mcmcInterval = 10,
                  mcmcSampleSize = 100,
                  verbose = FALSE)
    model <- model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    model$calculate()
    
    hamming_calc(edges,net)
    net_tmp <- net
    add.edges(net_tmp,1,2)
    r_stat_3 <- hamming_calc(edges,net_tmp)
    
    model$dyadUpdate(1,2)
    cpp_stat_3 <- model$statistics()
    
    # update not included edge that network has
    edges <- matrix(c(1,2),nrow = 1)
    model <- ernm(net ~ hamming(edges,100),
                  tapered = FALSE,
                  maxIter = 2,
                  mcmcBurnIn = 100,
                  mcmcInterval = 10,
                  mcmcSampleSize = 100,
                  verbose = FALSE)
    model <- model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    
    hamming_calc(edges,net)
    net_tmp <- net
    delete.edges(net_tmp,get.edgeIDs(net_tmp,1,61))
    r_stat_4 <- hamming_calc(edges,net_tmp)
    
    model$statistics()
    model$dyadUpdate(1,61)
    cpp_stat_4 <- model$statistics()
    
    # update not included edge that network does not have
    edges <- as.edgelist(net)*1.0
    edges <- cbind(as.double(edges[,1]),as.double(edges[,2]))
    head(edges)
    
    model <- ernm(net ~ hamming(edges,100),
                  tapered = FALSE,
                  maxIter = 2,
                  mcmcBurnIn = 100,
                  mcmcInterval = 10,
                  mcmcSampleSize = 100,
                  verbose = FALSE)
    model <- model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    model$calculate()

    hamming_calc(edges,net)
    net_tmp <- net
    add.edges(net_tmp,1,2)
    r_stat_5 <- hamming_calc(edges,net_tmp)
    
    model$statistics()
    model$dyadUpdate(1,2)
    cpp_stat_5 <- model$statistics()
    
    hamming_test_4 <- (r_stat_1 == cpp_stat_1)
    hamming_test_5 <- (r_stat_2 == cpp_stat_2)
    hamming_test_6 <- (r_stat_3 == cpp_stat_3)
    hamming_test_7 <- (r_stat_4 == cpp_stat_4)
    hamming_test_7 <- (r_stat_5 == cpp_stat_5)

    testthat::expect_true(hamming_test_1)
    testthat::expect_true(hamming_test_2)
    testthat::expect_true(hamming_test_3)
    testthat::expect_true(hamming_test_4)
    testthat::expect_true(hamming_test_5)
    testthat::expect_true(hamming_test_6)
    testthat::expect_true(hamming_test_7)
    
    # ========================
    # gaussRegression
    # ========================
    # Show 1) Linear regression Y ~ X gives the same results as ernm(net ~ GaussRegression(Y,X) | Y)
    y = (net %v% "var_4")
    X = (net %v% "var_5")
    
    r_stat_1 <- sum(y**2)
    r_stat_2 <- sum(y*X)
    
    cpp_stat_1 <- as.numeric(ernm::calculateStatistics(net ~ gaussRegression("var_4","var_5") | var_4))
    cpp_stat_2 <- cpp_stat_1[1]
    cpp_stat_1 <- cpp_stat_1[2]
    regression_test_1 <- abs(r_stat_1 - cpp_stat_1) <= 10**-10
    regression_test_2 <- abs(r_stat_2 - cpp_stat_2) <= 10**-10

    # Check the continVarUpdate function works change first y = 10
    k = 10
    r_stat_3 <- sum(y*X) - y[k]*X[k] + 1*X[k]
    r_stat_4 <- sum(y**2) - y[k]**2 + 1**2
    r_stat_5 <- sum(y*X) - y[k]*X[k] + y[k]*1
    
    model <- ernm(net ~ gaussRegression("var_4","var_5") | var_4,
                  tapered = FALSE,
                  maxIter = 2,
                  mcmcBurnIn = 100,
                  mcmcInterval = 10,
                  mcmcSampleSize = 100,
                  verbose = FALSE)
    model <- model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    model$calculate()
    model$continVertexUpdate(k, "var_4",1)
    
    regression_test_3 <- (abs(r_stat_3 - model$statistics()[1]) <= 10**-10)
    regression_test_4 <- (abs(r_stat_4 - model$statistics()[2]) <= 10**-10)
    
    model$calculate()
    model$continVertexUpdate(k, "var_5",1)
    regression_test_5 <- (abs(r_stat_5 - model$statistics()[1]) <= 10**-10)
    
    testthat::expect_true(regression_test_1)
    testthat::expect_true(regression_test_2)
    testthat::expect_true(regression_test_3)
    testthat::expect_true(regression_test_4)
    testthat::expect_true(regression_test_5)
    
    # ========================
    # regressNeighbors
    # ========================
    # Test : 
    # - Change in Y var
    # - Change in X var
    # - dyad update off
    # - dyad update on
    k <- 10
    neighs <- get.neighborhood(net,k)
    Y <- (net %v% "var_4")
    X <- (net %v% "var_5")
    r_stat_1 <- sum(sapply(1:(net %n% 'n'), function(i){
        neighs <- get.neighborhood(net,i)
        sum(Y[i]*X[neighs])
        }
    ))
    
    model <- ernm(net ~ regressNeighbors("var_4","var_5") | var_4,
                  tapered = FALSE,
                  maxIter = 10,
                  mcmcBurnIn = 100,
                  mcmcInterval = 10,
                  mcmcSampleSize = 100,
                  nodeSamplingPercentage = 1,
                  verbose = FALSE)
    model <- model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    model$calculate()
    model$statistics()
    regress_neigh_test_1 <- (abs(r_stat_1 -  model$statistics()) <= 10**-10)
    
    k = 10
    r_stat_2 <- sum(sapply(1:(net %n% 'n'), function(i){
        neighs <- get.neighborhood(net,i)
        if(i==k){
            sum(10*X[neighs])
        }else{
            sum(Y[i]*X[neighs])
        }
    }))
    
    X_tmp <- X
    X_tmp[k] <- 10
    r_stat_3 <- sum(sapply(1:(net %n% 'n'), function(i){
        neighs <- get.neighborhood(net,i)
        sum(Y[i]*X_tmp[neighs])
    }))
    
    neighs <- get.neighborhood(net,k)
    delete.edges(net,get.edgeIDs(net,k,neighs[1]))
    old_neighs <- neighs
    neighs <- get.neighborhood(net,k)
    r_stat_4 <- sum(sapply(1:(net %n% 'n'), function(i){
        neighs <- get.neighborhood(net,i)
        sum(Y[i]*X[neighs])
    }))
    add.edges(net,k,old_neighs[1])
    neighs_4 <- old_neighs
    
    neighs <- get.neighborhood(net,k)
    add.edges(net,k,(1:N)[-c(neighs,k)][1])
    old_neighs <- neighs
    neighs <- get.neighborhood(net,k)
    r_stat_5 <- sum(sapply(1:(net %n% 'n'), function(i){
        neighs <- get.neighborhood(net,i)
        sum(Y[i]*X[neighs])
    }))
    delete.edges(net,get.edgeIDs(net,k,(1:N)[-c(old_neighs,k)][1]))
    neighs_5 <- old_neighs
    
    
    model$continVertexUpdate(k, "var_4",10)
    regress_neigh_test_2 <- (abs(r_stat_2 -  model$statistics()) <= 10**-10)
    model$calculate()
    model$continVertexUpdate(k, "var_5",10)
    regress_neigh_test_3 <- (abs(r_stat_3 -  model$statistics()) <= 10**-10)
    model$calculate()
    model$dyadUpdate(k,neighs_4[1])
    regress_neigh_test_4 <- (abs(r_stat_4 -  model$statistics()) <= 10**-10)
    model$calculate()
    model$dyadUpdate(k,(1:N)[-c(neighs_5,k)][1])
    regress_neigh_test_5 <- (abs(r_stat_5 -  model$statistics()) <= 10**-10)
    
    
    testthat::expect_true(regress_neigh_test_1)
    testthat::expect_true(regress_neigh_test_2)
    testthat::expect_true(regress_neigh_test_3)
    testthat::expect_true(regress_neigh_test_4)
    testthat::expect_true(regress_neigh_test_5)
    
    # ========================
    # absDiff
    # ========================
    # Test : 
    # - Change in Y var
    # - dyad update off
    # - dyad update on
    k <- 10
    Y <- (net %v% "var_4")
    tmp <- as.edgelist(net)
    r_stat_1 <- sum(abs(Y[tmp[,1]] - Y[tmp[,2]]))
    
    model <- ernm(net ~ absDiff("var_4") | var_4,
                  tapered = FALSE,
                  maxIter = 10,
                  mcmcBurnIn = 100,
                  mcmcInterval = 10,
                  mcmcSampleSize = 100,
                  nodeSamplingPercentage = 1,
                  verbose = FALSE)
    model <- model$m$sampler$getModel()
    model$setNetwork(ernm::as.BinaryNet(net))
    model$calculate()
    model$statistics()
    absDiff_test_1 <- (abs(r_stat_1 -  model$statistics()) <= 10**-10)
    
    neighs <- get.neighborhood(net,k)
    r_stat_2 <- r_stat_1 - abs(Y[k] - Y[neighs[1]])
    neighs_2 <- neighs
    
    neighs <- get.neighborhood(net,k)
    r_stat_3 <- r_stat_1 + abs(Y[k] - Y[(1:N)[-c(neighs,k)][1]])
    neighs_3 <- neighs
    
    model$calculate()
    model$statistics()
    model$dyadUpdate(k,neighs_2[1])
    model$statistics()
    absDiff_test_2 <- (abs(r_stat_2 -  model$statistics()) <= 10**-10)
    model$calculate()
    model$dyadUpdate(k,(1:N)[-c(neighs_3,k)][1])
    absDiff_test_3 <- (abs(r_stat_3 -  model$statistics()) <= 10**-10)
    
    
    testthat::expect_true(absDiff_test_1)
    testthat::expect_true(absDiff_test_2)
    testthat::expect_true(absDiff_test_3)
}
)
