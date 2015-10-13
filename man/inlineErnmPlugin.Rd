\name{inlineErnmPlugin}
\alias{inlineErnmPlugin}
\title{An ernm plug-in for easy C++ prototyping and access}
\usage{
  inlineErnmPlugin(...)
}
\arguments{
  \item{...}{plug-in arguments}
}
\description{
  An ernm plug-in for easy C++ prototyping and access

  The inline plug-in for ernm
}
\examples{
library(inline)
registerPlugin("ernm",inlineErnmPlugin)
src <- "
		Rcpp::IntegerMatrix tmp(0,2);
		ernm::BinaryNet<ernm::Directed> net(tmp,Rcpp::as<int>(n));
		return net;
		"
emptyNetwork <- cxxfunction(signature(n="integer"), src, plugin="ernm")
net <- emptyNetwork(10)
net[1:10,1:10]
}

