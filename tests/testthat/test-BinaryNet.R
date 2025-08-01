# 
# Author: ianfellows
###############################################################################


library(testthat)
library(ernm)
library(network)

context("BinaryNet")

# test_that("DirectedNet",{
# 
# 			data(samplike)
# 			nw <- samplike
# 			el <- as.matrix(nw,matrix.type="edgelist")
# 			nVerts <- 18
# 			net <- as.BinaryNet(nw)
# 			expect_true(all(el[order(el[,1], el[,2]),]==net$edges()))
# 
# 			expect_true(all(list.vertex.attributes(nw) %in% net$variableNames(TRUE)))
# 
# 			net1 <- as.network.DirectedNet(net)
# 			expect_identical(network.vertex.names(nw),network.vertex.names(net1))
# 			net2 <- new(DirectedNet,el,nVerts)
# 			expect_true(all(as.matrix(net1,matrix.type="edgelist")==el[order(el[,1],el[,2]),]))
# 
# 			expect_true(all(nw[1:10,1:5]==net1[1:10,1:5]))
# 
# 			expect_true(all(rowSums(as.matrix(nw))==net$outDegree(1:nVerts)))
# 			expect_true(all(colSums(as.matrix(nw))==net$inDegree(1:nVerts)))
# 
# 			net$setDyads(c(1,2,2),c(4,5,6),c(TRUE,FALSE,NA))
# 			expect_identical(net$getDyads(c(1,2,2),c(4,5,6)),c(TRUE,FALSE,NA))
# 
# 			net1[1:5,1:2] <- matrix(rep(1,10),ncol=2)
# 			expect_identical(as.logical(net1[1:5,1:2]),
# 					structure(c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE),
# 					.Dim = c(1L, 10L)))
# 
# 			expect_true(is.null(net[["dsd"]]))
# 			f <- nw %v% "group"
# 			net[["group"]] <- f
# 			expect_true(identical(net[["group"]],as.factor(f)))
# 			net[["group"]] <- f
# 
# 			f <- as.factor(f)
# 			f[1:floor(length(f)/2)] <- NA
# 			net[["f"]] <- f
# 			expect_true(identical(net[["f"]],f))
# 
# 			rn <- rnorm(nVerts)
# 			net[["norm"]] <- rn
# 			expect_true(all(net[["norm"]] == rn))
# 			net[["norm"]] <- NULL
# 			expect_true(is.null(net[["norm"]]))
# 
# 			expect_identical(net$nMissing(1:18), c(0L, 1L, 0L, 0L, 0L, 0L, 0L,
# 							0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
# 							0L, 0L, 0L))
# 
# 			net[1,1:18] <- rep(NA,18)
# 			expect_identical(net$edges(),structure(c(2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L,
# 				4L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L, 7L, 7L, 7L, 7L, 7L, 7L,
# 				8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 9L, 9L, 10L, 10L, 10L, 10L, 10L,
# 				11L, 11L, 11L, 11L, 11L, 12L, 12L, 12L, 12L, 12L, 12L, 13L, 13L,
# 				13L, 13L, 13L, 14L, 14L, 14L, 14L, 14L, 15L, 15L, 15L, 16L, 16L,
# 				16L, 16L, 16L, 17L, 17L, 17L, 17L, 18L, 18L, 18L, 18L, 18L, 18L,
# 				1L, 3L, 9L, 15L, 1L, 2L, 7L, 8L, 14L, 1L, 2L, 3L, 5L, 9L, 1L,
# 				2L, 4L, 6L, 2L, 5L, 7L, 11L, 14L, 1L, 3L, 4L, 5L, 9L, 10L, 9L,
# 				10L, 11L, 13L, 2L, 7L, 8L, 10L, 11L, 16L, 8L, 9L, 11L, 12L, 13L,
# 				6L, 8L, 9L, 10L, 12L, 2L, 8L, 9L, 10L, 11L, 13L, 7L, 8L, 9L,
# 				10L, 14L, 8L, 9L, 10L, 12L, 13L, 2L, 13L, 18L, 9L, 10L, 15L,
# 				17L, 18L, 10L, 15L, 16L, 18L, 9L, 10L, 13L, 15L, 16L, 17L), .Dim = c(83L,
# 				2L))
# 			)
# 			expect_true(length(net$outNeighbors(c(1))[[1]])==0)
# 			net[1:18,1] <- rep(NA,18)
# 			expect_true(all(net$outNeighbors(c(2))[[1]]==c(3,9,15)))
# 			expect_true(all(net$outDegree(c(1,2))==c(0,3)))
# 
# 			a <- 1:nVerts
# 			attr(a,"lowerBound") <- 1
# 			attr(a,"upperBound") <- nVerts
# 			net[["a"]] <- a
# 			expect_true(attr(net[["a"]],"lowerBound")==1)
# 			expect_true(attr(net[["a"]],"upperBound")==17)
# 		})
# 
# 
# test_that("UndirectedNet",{
# 
# 			data(flo)
# 			nw<-network(flo,directed=FALSE)
# 			el <- as.matrix(nw,matrix.type="edgelist")
# 			nVerts <- 16
# 
# 			net <- as.BinaryNet(nw)
# 			expect_true(all(el[order(el[,2]),2:1]==net$edges()))
# 
# 			expect_true(all(list.vertex.attributes(nw) %in% net$variableNames(TRUE)))
# 
# 			net1 <- as.network(net)
# 			expect_identical(network.vertex.names(nw),network.vertex.names(net1))
# 			net2 <- new(UndirectedNet,el,nVerts)
# 			expect_true(all(as.matrix(net1,matrix.type="edgelist")==el[order(el[,2]),2:1]))
# 
# 			expect_true(all(nw[1:10,1:5]==net[1:10,1:5]))
# 
# 			expect_true(all(rowSums(as.matrix(nw))==net$degree(1:nVerts)))
# 
# 			net$setDyads(c(1,2,2),c(4,5,6),c(TRUE,FALSE,NA))
# 			expect_identical(net$getDyads(c(1,2,2),c(4,5,6)),c(TRUE,FALSE,NA))
# 
# 			net[1:5,1:2] <- matrix(rep(1,10),ncol=2)
# 			expect_identical(net[1:5,1:2],
# 					structure(c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE,
# 									TRUE, TRUE), .Dim = c(5L, 2L)))
# 
# 			expect_true(is.null(net[["dsd"]]))
# 			nw %v% "group" <- rep(c("a","b","c","d"),4)
# 			f <- nw %v% "group"
# 			net[["group"]] <- f
# 			expect_true(identical(net[["group"]],as.factor(f)))
# 
# 			f <- as.factor(f)
# 			f[1:floor(length(f)/2)] <- NA
# 			net[["f"]] <- f
# 			expect_true(identical(net[["f"]],f))
# 
# 			rn <- rnorm(nVerts)
# 			net[["norm"]] <- rn
# 			expect_true(all(net[["norm"]] == rn))
# 			net[["norm"]] <- NULL
# 			expect_true(is.null(net[["norm"]]))
# 
# 			expect_identical(net$nMissing(1:nVerts), c(0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
# 							0L))
# 
# 			net[1,1:nVerts] <- rep(NA,nVerts)
# 			expect_true(length(net$neighbors(c(1))[[1]])==0)
# 			expect_true(all(net$degree(c(1))==c(0)))
# 
# 			a <- 1:nVerts
# 			attr(a,"lowerBound") <- 1
# 			attr(a,"upperBound") <- nVerts
# 			net[["a"]] <- a
# 			expect_true(attr(net[["a"]],"lowerBound")==1)
# 			expect_true(attr(net[["a"]],"upperBound")==17)
# 		})
		
