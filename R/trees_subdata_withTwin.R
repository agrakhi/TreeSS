#' TreeSS predictions when using Twinning from each stratum
#'
#' @description The function first creates a best partition and then find the expected sample size from each node. 
#' Post that, it samples from each stratum using \code{twinning}. Finally, it builds a random forest on this
#' sample and find final predicted values.
#'
#'
#' @param pop a full dataset with the first column being the response and remaining p columns as the predictors
#'
#' @param pop.val a test dataset with the first column being the response and remaining p columns as the predictors
#'
#' @param n the desired sample size.
#'
#' @param NumCores the number of cores the computer has to make use of parallel computing
#'
#' @param t_s the size of the observations on which the tree used for partition should be built. By default, 
#' it is equal to sqrt(\code{N})
#'
#' @param t_d the maximum depth of the tree for the partition. The \code{t_n} trees from depth 3 to \code{t_d} are 
#' built. By default, it is equal to log(sqrt(\code{N})), where log is base 2.
#'
#' @param thhold the minimum number of observations to be taken from each stratum. By default, 
#' it is set to 5.
#'
#' @param t_n the number of trees to be built for each depth. By default, it is set to 10.
#'
#' @param t_p the number of variables to be tried on each node while creating the partition. By default, 
#' it is set to p.
#'
#' @param ntree the number of trees to be built when random froest is built on each sample for the final predictions.
#' By default, it is set to 100.
#'
#' @param ptry the number of variables to be tried at each split when random froest is built on each sample for the
#' final predictions. By default, it is set at p if p is at most 5, and at p/3 otherwise. 
#'
#' @return A list with three objects. First object is indices indicating the selected TreeSS sample. The second object 
#' is the mean squared prediction error on the test dataset, whereas the third
#' object is the predicted values on the test dataset.
#'
#' @source Singh, R. and Stufken, J. (2023+). Model-free tree-based subdata selection, 1--31. 
#' \doi{10.48550/arXiv.2205.13497}
#' 
#' @seealso \code{\link{trees_subdata_withUni}}
#'
#' @examples
#' ## 1. A simulted dataset of 10000 observations and 2 variables with response 
#' ## as Michalewicz function plus some error
#' NumCores = min(1,parallel :: detectCores() - 4)
#' X1 = rnorm(n=10000, mean=0, sd=1)
#' X2 = rnorm(n=10000, mean=0, sd=1)
#' Y1 = sin(pi*X1)*(sin(pi*X1^2))^20 + sin(pi*X2)*(sin(2*pi*X2^2))^20  
#' Y = Y1+rnorm(n=10000, mean = 0, sd = (max(Y1)-min(Y1))/10^3)
#' data = cbind(Y, X1, X2)
#' 
#' # creating train-test split using twinning of 80-20
#' twin_indices = twinning :: twin(data, r=5) 
#' pop.val = data[twin_indices,]
#' pop = data[-twin_indices,]
#' n = 1000
#' 
#' N = dim(pop)[1]
#' p = dim(pop)[2]-1
#' trees_subdata = trees_subdata_withTwin(pop, pop.val, n, NumCores = NumCores) 
#' print(trees_subdata$valTreeSS)
#' 
#' 
#'
#' @export

trees_subdata_withTwin <- function(pop, pop.val, n, NumCores = 8, t_s = c(), t_d = c(), thhold = 5, t_n = 10, t_p = c(), ntree = 100, ptry = c()){
  #tthold is the same as t_h in the paper
  sz = dim(pop)
  N = sz[1]
  p = sz[2] - 1 # number of predictors
  
  `%dopar%` <- foreach::`%dopar%`
  
  if(is.null(t_s)){
	t_s = round(sqrt(N))
  }
   
  if(is.null(t_d)){
	t_d = round(log2(sqrt(N)))
  }

  if(is.null(t_p)){
	t_p = p
  }
       
  # depth.vec is a vector of (unique) choices for depth of each tree
  depth.vec = 3: min(t_d, floor(log2(n)))
  
  # depth.all.vec repeats each depth choice 10 times and hence eventually we build 
  # 10 trees of each depth choice
  depth.all.vec = rep(depth.vec,each=t_n)
  
  #number of built trees would be the same as the length of the vector depth.all.vec
  treech = length(depth.all.vec)
  
  #setup parallel backend to use many processors
   cl <- parallel :: makeCluster(NumCores) #do not overload your computer
   doParallel ::registerDoParallel(cl)
   finalMatrix<- foreach :: foreach (temploop = 1:treech, .combine=cbind, .packages=c("ranger"))%dopar% {
    temploop= 1
	depthch = depth.all.vec[temploop]
	
    if (t_p==1){
      Sel.vars = seq(2,p+1)
    } else{
      Sel.vars = sample(seq(2,p+1),t_p, replace=F)
    }  
    randsample = sample(seq(1,N),t_s, replace=F)  
	
    tree.result<- ranger :: ranger(x= pop[randsample,Sel.vars,drop=F], y=pop[randsample,1],num.trees = 1,replace=F, sample.fraction=1,max.depth = depthch, mtry =  length(Sel.vars))
    
    ##### Predicting terminal Nodes for each tree
    predID = as.integer(stats :: predict(tree.result, type="terminalNodes", data=pop[, Sel.vars,drop=F])$predictions)
   
    Ni = table(predID)
    
    # finding S_i for each node, if there is only one observation in a node then we assign S_i = 0  
    Si=(stats :: aggregate(pop[,1], by=list(Category=predID), FUN=stats :: sd)[,2])
    Si[is.na(Si)]=0
    MSPEall = sum(((Ni-1)*Si^2))
   gc()
   
   # storing how many nodes in each tree
   leng = length(Ni)
   c(list(MSPEall=MSPEall, leng=leng,Ni = Ni,Si = Si, predID= predID)) 
   }
 
  parallel :: stopCluster(cl)
 
  # extracting the information from partition
  MSPEall = finalMatrix[1,]
  leng = finalMatrix[2,]
  #finding the best tree 
  besttree = find_tree(MSPEall, leng, floor(n/thhold)) 
  stratsize = unlist(finalMatrix[3,besttree])
  stratsd= unlist(finalMatrix[4,besttree])
  predID = unlist(finalMatrix[5,besttree])
  
  valTreeSS = numeric(1)

  # the expected sample size from each node
  expected.SampleSize = SampleSizesNodes(stratsd, n, stratsize, thhold)
  
 
   predY.val = matrix(0, nrow= dim(pop.val)[1], ncol=1)

   cl <- parallel :: makeCluster(NumCores) #do not overload your computer
   doParallel ::registerDoParallel(cl)
 
   # twin sample from each strata
   finalMatrix2<- foreach :: foreach (temploop = 1:length(expected.SampleSize), .combine=cbind, .packages=c("ranger","twinning"))%dopar% {
    SpecificNode = as.numeric(names(expected.SampleSize)[temploop])
    r = expected.SampleSize[as.numeric(names(expected.SampleSize))==SpecificNode]
    rel_indices = which(predID == SpecificNode)
    poptempYXperP = pop[rel_indices, , drop=F]
    
    ratio_req1 = floor(length(rel_indices)/r) #round(1/0.2)
    if(ratio_req1 >= 2 & ratio_req1 <=  dim(poptempYXperP)[1]/2) {#} & ratio_req1)
    twinsamp = twinning::twin(poptempYXperP, r=ratio_req1)
    if (length(twinsamp) >  r){
    twinsampTemp = sample(1:length(twinsamp), r, replace=F)
    twinsamp = twinsamp[twinsampTemp]
    }
    } else {
      twinsamp = sample(1:length(rel_indices), r, replace=F)
    }
   c(list(indices = rel_indices[twinsamp]))
  }
	parallel :: stopCluster(cl)
	
  if (is.null(ptry)){  
	if (p <= 5){
		ptry = p
	} else {
		ptry = ceiling(p/3) 
	}
  }
   
	TreeSS.sample = unlist(finalMatrix2)
	
	# prediction using random forests
	forest.twin<- ranger :: ranger(x= pop[TreeSS.sample, seq(2,p+1),drop=F], y=pop[TreeSS.sample,1],num.trees = ntree,mtry = ptry)
    predY.val = stats :: predict(forest.twin, data=pop.val[, -1,drop=F])$predictions
    valTreeSS = base :: mean((pop.val[,1]-predY.val)^2)
   
  return(list("TreeSS.sample" =  TreeSS.sample, "valTreeSS" = valTreeSS, "predY.val" = predY.val)) 
 } 