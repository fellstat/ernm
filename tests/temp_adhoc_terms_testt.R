# This script carries out testing for the newly coded up ERNM terms:

library(statnet)
library(ernm)

# =====================================================================
# Setup
# =====================================================================
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

#=====================================================================
# nodeCountTopLevel
#=====================================================================
# Check what the logistic term actually does:
for(i in 1:length(nets)){
  print(paste("net ",i ))
  method_1 = sum(as.numeric(nets[[i]] %v% "var_2" == "2"))
  method_2 = as.numeric(ernm::calculateStatistics(nets[[i]] ~ nodeCountTopLevel("var_2")))
  if(method_1 != method_2){
    print("FAILED!!!")
    print(method_1)
    print(method_2)
  }
}
#=====================================================================
# Logistic
#=====================================================================
net = nets[[2]]
v1 = (net %v% "var_2")
v2 = (net %v% "var_3")

sum((v1=="1")[v2=="1"])
sum((v1=="1")[v2=="2"])
sum((v1=="2")[v2=="1"])
sum((v1=="2")[v2=="2"])

as.numeric(ernm::calculateStatistics(net ~ logistic("var_2","var_3") | var_2))
as.numeric(ernm::calculateStatistics(net ~ logisticTopLevel("var_2","var_3") | var_2))

debug(createCppModel)
ernm::createCppModel(net ~ logisticTopLevel("var_2","var_3"))

# Check what the logistic term actually does:
for(i in 1:length(nets)){
  print(paste("net ",i ))
  method_1 = sum(as.numeric(nets[[i]] %v% "var_2" == "2")[which(nets[[i]] %v% "var_3" == "2")])
  method_2 = as.numeric(ernm::calculateStatistics(nets[[i]] ~ logistic("var_2","var_3") | var_2))
  if(method_1 != method_2){
    print("FAILED!!!")
    print(method_1)
    print(method_2)
  }
}

# Conclusion "Logistic counts the number of occurences of the topLevel of the LHS var with the bottomLevel of the RHS var


# =====================================================================
# Logistic Neighbors
# =====================================================================
sum(as.numeric(net %v% "var_2_neighbors")[which(net %v% "var_2" == "1")])
ernm::calculateStatistics(net ~ logisticNeighbors("var_2","var_2") | var_2)

for(i in 1:length(nets)){
  print(paste("net ",i ))
  method_1 = sum(as.numeric(nets[[i]] %v% "var_2_neighbors")[which(nets[[i]] %v% "var_2" == "1")])
  method_2 = as.numeric(ernm::calculateStatistics(nets[[i]] ~ logisticNeighbors("var_2","var_2") | var_2))
  if(method_1 != method_2){
    print("FAILED!!!")
  }
}

# test the dyadUpdate function :
which(net %v% "var_2"== "2")
net_2 <- net
tmp_model <- ernm(net ~ logisticNeighbors("var_2","var_2") | var_2,
                  maxIter = 2,
                  mcmcSampleSize = 1000,
                  mcmcBurnIn = 100)
for(i in which(net %v% "var_2"== "2")[-1]){
  net_2 <- net
  net_2[i,2] <- 1 - net_2[i,2] 
  
  change_1 <- ernm::calculateStatistics(net_2 ~ logisticNeighbors("var_2","var_2") | var_2) - ernm::calculateStatistics(net ~ logisticNeighbors("var_2","var_2") | var_2)
  
  model <- tmp_model
  
  model <- model$m$sample$getModel()
  model$setNetwork(as.BinaryNet(net))
  model$calculate()
  
  stat1 <- model$statistics()
  model$dyadUpdate(2,i)
  stat2 <- model$statistics()
  change_2 <- stat2-stat1
  
  if(change_1 != change_2){
    print("PROBLEM")
    print(i)
    print(change_1)
    print(change_2)
  }
}
 # looks okay 

# test for changing a node value to 

for(i in 1:(net%n% "n")){
  net_2 <- net
  if((net_2 %v% "var_2")[i] == "1"){
    new_value <- "2"
  }else{
    new_value <- "1"
  }
  set.vertex.attribute(net_2,"var_2",new_value,i)
  
  change_1 <- ernm::calculateStatistics(net_2 ~ logisticNeighbors("var_2","var_2") | var_2) - ernm::calculateStatistics(net ~ logisticNeighbors("var_2","var_2") | var_2)
  
  model <- tmp_model
  
  model <- model$m$sample$getModel()
  model$setNetwork(as.BinaryNet(net))
  model$calculate()
  
  stat1 <- model$statistics()
  model$discreteVertexUpdate(i,"var_2",as.numeric(new_value))
  stat2 <- model$statistics()
  change_2 <- stat2-stat1
  
  if(change_1 != change_2){
    print("PROBLEM")
    print(i)
    print(change_1)
    print(change_2)
  }
}

# looks fine ! 

# =====================================================================
# Logistic Neighbors Top Level
# =====================================================================
sum(as.numeric(net %v% "var_2_neighbors")[which(net %v% "var_2" == "2")])
ernm::calculateStatistics(net ~ logisticNeighborsTopLevel("var_2","var_2") | var_2)

# test the dyadUpdate function :
which(net %v% "var_2"== "2")
net_2 <- net
tmp_model <- ernm(net ~ logisticNeighborsTopLevel("var_2","var_2") | var_2,
                  maxIter = 2,
                  mcmcSampleSize = 1000,
                  mcmcBurnIn = 100)
for(i in which(net %v% "var_2"== "2")[-1]){
  net_2 <- net
  net_2[i,2] <- 1 - net_2[i,2] 

  change_1 <- ernm::calculateStatistics(net_2 ~ logisticNeighborsTopLevel("var_2","var_2") | var_2) - ernm::calculateStatistics(net ~ logisticNeighborsTopLevel("var_2","var_2") | var_2)
  
  model <- tmp_model
  
  model <- model$m$sample$getModel()
  model$setNetwork(as.BinaryNet(net))
  model$calculate()
  
  stat1 <- model$statistics()
  model$dyadUpdate(2,i)
  stat2 <- model$statistics()
  change_2 <- stat2-stat1
  
  if(change_1 != change_2){
    print("PROBLEM")
    print(i)
    print(change_1)
    print(change_2)
  }
}

for(i in 1:(net%n% "n")){
  net_2 <- net
  if((net_2 %v% "var_2")[i] == "1"){
    new_value <- "2"
  }else{
    new_value <- "1"
  }
  set.vertex.attribute(net_2,"var_2",new_value,i)
  
  change_1 <- ernm::calculateStatistics(net_2 ~ logisticNeighborsTopLevel("var_2","var_2") | var_2) - ernm::calculateStatistics(net ~ logisticNeighborsTopLevel("var_2","var_2") | var_2)
  
  model <- tmp_model
  
  model <- model$m$sample$getModel()
  model$setNetwork(as.BinaryNet(net))
  model$calculate()
  
  stat1 <- model$statistics()
  model$discreteVertexUpdate(i,"var_2",as.numeric(new_value))
  stat2 <- model$statistics()
  change_2 <- stat2-stat1
  
  if(change_1 != change_2){
    print("PROBLEM")
    print(i)
    print(change_1)
    print(change_2)
  }
}


# =====================================================================
# Differential Homophilly
# =====================================================================
# Already includeed ! 
ernm::calculateStatistics(net ~ homophily(name = "var_1", collapse = F, mix = F) | var_2)
ergm::summary_formula(net ~ nodematch("var_1",diff = T))

# =====================================================================
# Homophillous ESP
# =====================================================================
tmp <- sapply(unique(net %v% "var_1"),function(i){
  remove <- which((net %v% "var_1") != i) 
  net_tmp <-net
  delete.vertices(net_tmp,remove)
  return(ergm::summary_formula(net_tmp ~ gwesp(0.5,fixed = T)))
})
names(tmp) <- unique(net %v% "var_1")

# ERNM Calcs
ernm::calculateStatistics(net ~ gwesp(alpha = 0.5,homogenous = T,variableName = "var_1"))

# Check that regular ESP is still working:
ergm::summary_formula(net ~ gwesp(0.5,fixed = T))
ernm::calculateStatistics(net ~ gwesp(0.5))

# Check the dyad update works :
tmp_model <- ernm(net ~ gwesp(alpha = 0.5,homogenous = T,variableName = "var_1"),
                  maxIter = 2,
                  mcmcSampleSize = 1000,
                  mcmcBurnIn = 100)

for(i in 1:(net%n%"n")){
  for(j in 1:(net%n%"n")){
    if(i==j){next}
    net_2 <- net
    net_2[i,j] <- 1 - net_2[i,j] 
    change_1 <- ernm::calculateStatistics(net_2 ~ gwesp(alpha = 0.5,homogenous = T,variableName = "var_1")) - 
      ernm::calculateStatistics(net ~ gwesp(alpha = 0.5,homogenous = T,variableName = "var_1"))
  
    model <- tmp_model
    model <- model$m$sample$getModel()
    model$setNetwork(as.BinaryNet(net))
    model$calculate()
    
    stat1 <- model$statistics()
    model$dyadUpdate(i,j)
    stat2 <- model$statistics()
    change_2 <- stat2-stat1
    
    if(abs(change_1 -  change_2) > 10**-3){
      print("PROBLEM")
      print(i)
      print(j)
      print(change_1)
      print(change_2)
      Sys.sleep(1)
    }
  }
}

# Check reuglar gwesp is still working
tmp_model <- ernm(net ~ nodeMatch("var_3"))
tmp_model <- ernm(net ~ gwesp(alpha = 0.5,homogenous = T,variableName = "var_2"),
                  maxIter = 2,
                  mcmcSampleSize = 1000,
                  mcmcBurnIn = 100)

m <- tmp_model$m$sampler$getModel()
for(i in 1:1000){
  print(i)
  tmp <- as.BinaryNet(net)
}
# after enough calls something in c++ breaks.. 


delete.vertex.attribute(net,"na")
for(i in which(net %v% "var_2"== "2")){
  for(j in which(net %v% "var_2"== "2")){
    if(i==j){next}
    print(i)
    print(j)
    net_2 <- net
    net_2[i,j] <- 1 - net_2[i,j] 
    change_1 <- tryCatch({ernm::calculateStatistics(net_2 ~ gwesp(alpha = 0.5,homogenous = T,variableName = "var_2")) - 
      ernm::calculateStatistics(net ~ gwesp(alpha = 0.5,homogenous = T,variableName = "var_2"))},
      error = function(e){return(10**10)})
    
    model <- tmp_model$m$sampler$getModel()
    
    model$setNetwork(as.BinaryNet(net))
    model$calculate()
    stat1 <- tryCatch({model$statistics()},error = function(e){10**5})
    model$dyadUpdate(i,j)
    stat2 <- tryCatch({model$statistics()},error = function(e){10**7})
    change_2 <- stat2-stat1
    if(abs(change_1 -  change_2) > 10**-3){
      print("PROBLEM")
      print(i)
      print(j)
      print(change_1)
      print(change_2)
      Sys.sleep(1)
    }
  }
}

# =====================================================================
# hamming Function
# =====================================================================

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


ernm::calculateStatistics(net ~ hamming(edges,100))
hamming_calc(edges,net)

edges[1,] = c(1,11)
ernm::calculateStatistics(net ~ hamming(edges,100))
hamming_calc(edges,net)

edges[2,] = c(1,13)
ernm::calculateStatistics(net ~ hamming(edges,100))
hamming_calc(edges,net)


# check dyad update function
edges = as.edgelist(net)*1
edges = cbind(as.double(edges[,1]),as.double(edges[,2]))
head(edges)

# update included edge that network has
model = ernm(net ~ hamming(edges,100),maxIter = 2,mcmcBurnIn = 100,mcmcInterval = 10, mcmcSampleSize = 100)
model = model$m$sampler$getModel()
model$setNetwork(ernm::as.BinaryNet(net))
model$calculate()
model$statistics()

print("R calc")
hamming_calc(edges,net)
net_tmp = net
delete.edges(net_tmp,get.edgeIDs(net_tmp,1,61))
hamming_calc(edges,net_tmp)

print("Cpp calc")
model$statistics()
model$dyadUpdate(1,61)
model$statistics()

print("R calc")
hamming_calc(edges,net)
net_tmp = net
delete.edges(net_tmp,get.edgeIDs(net_tmp,61,1))
hamming_calc(edges,net_tmp)

print("Cpp calc")
model$statistics()
model$dyadUpdate(61,1)
model$statistics()



# update included edge that network does not have
edges = matrix(c(1,2),nrow = 1)
model = ernm(net ~ hamming(edges),maxIter = 2,mcmcBurnIn = 100,mcmcInterval = 10, mcmcSampleSize = 100)
model = model$m$sampler$getModel()
model$setNetwork(ernm::as.BinaryNet(net))
model$calculate()
model$statistics()

hamming_calc(edges,net)
net_tmp = net
add.edges(net_tmp,1,2)
hamming_calc(edges,net_tmp)

model$statistics()
model$dyadUpdate(1,2)
model$statistics()

# update not included edge that network has
edges = matrix(c(1,2),nrow = 1)
model = ernm(net ~ hamming(edges),maxIter = 2,mcmcBurnIn = 100,mcmcInterval = 10, mcmcSampleSize = 100)
model = model$m$sampler$getModel()
model$setNetwork(ernm::as.BinaryNet(net))
model$calculate()
model$statistics()

hamming_calc(edges,net)
net_tmp = net
delete.edges(net_tmp,get.edgeIDs(net_tmp,1,61))
hamming_calc(edges,net_tmp)

model$statistics()
model$dyadUpdate(1,61)
model$statistics()

# update not included edge that network does not have
edges = as.edgelist(net)*1.0
edges = cbind(as.double(edges[,1]),as.double(edges[,2]))
head(edges)

model = ernm(net ~ hamming(edges),maxIter = 2,mcmcBurnIn = 100,mcmcInterval = 10, mcmcSampleSize = 100)
model = model$m$sampler$getModel()
model$setNetwork(ernm::as.BinaryNet(net))
model$calculate()
model$statistics()

hamming_calc(edges,net)
net_tmp = net
add.edges(net_tmp,1,2)
hamming_calc(edges,net_tmp)

model$statistics()
model$dyadUpdate(1,2)
model$statistics()