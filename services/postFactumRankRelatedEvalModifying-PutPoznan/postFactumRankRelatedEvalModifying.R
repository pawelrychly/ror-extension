##########################################
# Usage:
# R --slave --vanilla --args "[inDirectory]" "[outDirectory]" < postFactumRankRelatedEvalModyfying.R
# Example: 
# R --slave --vanilla --args "${PWD}/in" "${PWD}/out" < postFactumRankRelatedEvalModyfying.R
##########################################

inDirectory <- commandArgs()[5] #"D:/projects/praca-magisterska/ROR-ranking/services/postFactumRankRelatedEvalModifying-PutPoznan/tests/in2" # 
outDirectory <- commandArgs()[6] # "D:/projects/praca-magisterska/ROR-ranking/services/postFactumRankRelatedEvalModifying-PutPoznan/tests/out2" # 

##########################################
# Set the working directory as the "in" directory
##########################################

setwd(inDirectory)

errFile<-NULL
errData<-NULL
errCalc<-NULL
execFlag<-FALSE
M_BIG_NUMBER <- 2 
initLibrary <- function() {
  options( java.parameters = "-Xmx2g" )
  library(ror)
  library(RXMCDA)
  library(rJava)
  assignInNamespace("MINEPS", ns="ror", value=0.0001)
  M_BIG_NUMBER <- 2 
}  

############################################################################
############################ Helper functions ##############################
############################################################################

SolveModel <- function(perf, model) {
  
  objective <- ror:::buildObjectiveFunction(perf)
  objective <- GetNormalizedMatrix(matrix=objective, width=ncol(model$lhs))
  obj <- L_objective(objective)
  roiConst <- L_constraint(model$lhs, model$dir, model$rhs)
  lp <- OP(objective=obj, constraints=roiConst, maximum=TRUE)
  ROI_solve(lp, ror:::.solver)
}

CheckConstraintsConsistency <- function(perf, model) {
  #check constraints consistency  
  #  model: structure of constraints with the following elements:
  #    model$lhs: matrix - left side of constraints
  #    model$dir: list of operators
  #    model$rhs: matrix - right side of constraints
  
  
  ret <- SolveModel(perf, model)
  return(ret$status$code == 0 && ret$objval >= ror:::MINEPS)
}


ShowSolutionForRanks <- function(perf,  model = NULL, number.of.binary.variables, nums.of.characteristic.points = NULL, ret=NULL) {
  
  if (is.null(ret) && is.null(model)) {
    stop("Konieczne jest podanie modelu problemu (model) lub jego rozwiązania (ret)")
  }
  if (is.null(ret)) {
    ret <- ret <- PfaRanksSolveModel(perf, model, number.of.binary.variables)
  }
  if (ret$status$code != 0) {
    print("Model jest nieosiągalny przy podanych ograniczeniach")
  }
  par(mfrow = c(3,2))
  #print(perf)
  levels <- c()
  
  
  #print(model)
  
  levels <- ror:::getLevels(perf)
  solution <- ret$solution
  i <- 1
  if (!is.null(nums.of.characteristic.points)) {
    i <- sum(as.numeric(lapply(levels, length))) + 2
    
    levels = GetCharacteristicPoints(perf=perf, nums.of.characteristic.points = nums.of.characteristic.points)
  } else {
    i <- 1
  }
  
  #print(solution)
  
  for (values in levels) {
    labels <- c()
    v <- c()
    for (value in values) {
      labels <- c(labels, value)
      v <- c(v, solution[i])
      i = i+1
    }
    
    print("Punkty charakterystyczne:")
    print(labels)
    print("Wartości w punktach charakterystycznych")
    print(v)
    #plot(labels, v, type="b", col="red")
  }
  
}


# model <- jeśli ret == NULL szukane jest rozwiązanie problemu, na podstawie , którego prezentowany jest wynik
# ret <- wynik ( rozwiązanie problemu) - jeśli jest różny od null wartość model jest nie potrzebna
ShowSolution <- function(perf, model = NULL, nums.of.characteristic.points = NULL, ret=NULL) {
  
  if (is.null(ret) && is.null(model)) {
    stop("Konieczne jest podanie modelu problemu (model) lub jego rozwiązania (ret)")
  }
  if (is.null(ret)) {
    ret <- SolveModel(perf, model)
  } 
  if (ret$status$code != 0) {
    print("Model jest nieosiągalny przy podanych ograniczeniach")
  }
  par(mfrow = c(3,2))
  #print(perf)
  levels <- c()
  
  #print(model)
  #print(ret$solution)
  solution <- ret$solution
  
  levels <- ror:::getLevels(perf)
  
  mapping.all = list()
  #jesli odcinkami liniowe - mapowanie wartosci uzytecznosci na level
  j = 1;
  k = 1;
  for (values in levels) {
    mapping = list()
    for (value in values) {
      mapping[[toString(value)]] = solution[j]
      j = j+1
    }
    mapping.all[[toString(k)]] = mapping
    k = k + 1
  }
  
  num_of = 1
  #print("OCENY WARIANTÓW NA KRYTERIACH")
  #for(variant_index in seq(nrow(perf))) {
  #  print(sprintf("Oceny wariantu %d:", variant_index))
  #  for(value_index in seq(ncol(perf))) {
  #    v <- perf[[value_index]][[variant_index]]
  #    
  #    print(sprintf("Kryterium %d: %f", value_index , mapping.all[[toString(value_index)]][[toString(v)]]  ))
  #  }
  #  
  #}
  print(ret)
  i <- 1
  if (!is.null(nums.of.characteristic.points)) {
    i <- sum(as.numeric(lapply(levels, length))) + 2
    print(i)
    levels = GetCharacteristicPoints(perf=perf, nums.of.characteristic.points = nums.of.characteristic.points)
  } else {
    i <- 1
  }
  
  #print(solution)
  
  for (values in levels) {
    labels <- c()
    v <- c()
    for (value in values) {
      labels <- c(labels, value)
      v <- c(v, solution[i])
      i = i+1
    }
    
    print("Punkty charakterystyczne:")
    print(labels)
    print("Wartości w punktach charakterystycznych")
    print(v)
    plot(labels, v, type="b", col="red")
  }
  
}

EqualMatrix <- function(mat1, mat2) {
  
  if (ncol(mat1) != ncol(mat2) || (nrow(mat1) != nrow(mat2))) {
    return(FALSE)
  }
  width <- ncol(mat1)
  height <- nrow(mat1)
  
  for (i in seq(height)) {
    for (j in seq(width)) {
      if (mat1[i,j] != mat2[i,j]) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

GetNumberOfVariables <- function(perf, numbers.of.characteristic.points){
  # Function returns number of variables used in problem
  #
  # perf - matrix of options and their ratings on criteria
  # numbers.of.characteristic.points - sequence of numbers. Position i of this array means number of characteristic
  #                        points, on i - criteria
  # return - number of variables
  levels <- ror:::getLevels(perf);
  num.of.values = ror:::getNrVars(levels)
  num.of.characteristic.points = sum(numbers.of.characteristic.points);
  num.of.variables = num.of.values + num.of.characteristic.points;
  return(num.of.variables);
}

GetNumberOfValues <- function(perf) {
  levels <- ror:::getLevels(perf);
  nr.vars <- ror:::getNrVars(levels);
  return(nr.vars)
}

GetNormalizedMatrix <- function(matrix, width, right=TRUE) {
  #Function adds empty columns if width of matrix is to small
  # matrix - matrix
  # witdth - expected width of matrix
  # return normalized matrix
  if (!is.matrix(matrix)) {
    matrix = matrix(matrix, nrow=1)  
  }
  if (length(matrix) == 0) {
    return(NULL)
  }
  
  
  if (ncol(matrix) < width) {
    
    numOfAdditionalCols <- width - ncol(matrix);
    tempMatrix <- matrix(0.0, ncol=numOfAdditionalCols, nrow(matrix));
    if (right == TRUE) {
      matrix <- cbind(matrix, tempMatrix);
    } else {
      matrix <- cbind(tempMatrix, matrix);
    }
  }
  return(matrix);
}

GetLastElement <- function(vector) {
  #Function returns last element of vector
  tail(vector, n=1)
}

GetIntervalData <- function(characteristic.points, criterion.index, value) {
  i = 2;
  while((characteristic.points[[criterion.index]][i] < value) & (i <= length(characteristic.points[[criterion.index]]))) {
    i = i + 1;
  }
  
  bottom = 0
  top = 0
  if (i <= length(characteristic.points[[criterion.index]])) {
    bottom = characteristic.points[[criterion.index]][i-1];
    top = characteristic.points[[criterion.index]][i];
  }
  interval = list(index=i, value = c(bottom, top));
  return(interval)
} 

BuildPiecewiseLinearMonotonousConstraints <- function(perf, strict.vf=FALSE,
                                                      nums.of.characteristic.points) {
  
  stopifnot(is.logical(strict.vf));
  levels.list <- ror:::getLevels(perf);
  characteristic.points = GetCharacteristicPoints(perf, nums.of.characteristic.points);
  num.of.values = ror:::getNrVars(levels.list)
  num.of.characteristic.points = ror:::getNrVars(characteristic.points) -1;
  
  offsets.points = ror:::getOffsets(characteristic.points);
  offsets.points = offsets.points + num.of.values
  num.of.variables = num.of.values + num.of.characteristic.points;
  c1 <- BuildMonotonousConstraintsForCharacteristicPoints(
    perf, characteristic.points, offsets.points, num.of.variables, strict.vf)
  c2 <- BuildConstraintsForValuesEvaluation(
    characteristic.points, offsets.points, levels.list, num.of.variables)
  monotonousConst <- ror:::combineConstraints(c1, c2)
  return(monotonousConst)
}

GetCharacteristicPoints <- function(perf, nums.of.characteristic.points) {
  # Function return list of lists of characteristic points for each criterion
  # 
  # perf: the performance matrix
  # nums.of.characteristic.points: list of nums of characteristic points for each criterion 
  levels.list <- ror:::getLevels(perf);
  list.of.characteristic.points = list()
  nrAlts <- nrow(perf)
  nrCrit <- ncol(perf)
  for (i in c(1:nrCrit)){  #dla każdego kryterium
    list.of.characteristic.points[[i]] = vector(mode="numeric", length=0)
    num.of.characteristic.points = nums.of.characteristic.points[i]
    if (!is.numeric(num.of.characteristic.points)) {
      num.of.characteristic.points = 0
    }
    if (num.of.characteristic.points > 1) {
      for (j in c(1:num.of.characteristic.points)) {
        last.level = GetLastElement(levels.list[[i]])      
        range.of.performance = last.level - levels.list[[i]][1]
        characteristic.point = levels.list[[i]][1] + (range.of.performance)*(j-1)/(num.of.characteristic.points-1) 
        list.of.characteristic.points[[i]] = c(list.of.characteristic.points[[i]], characteristic.point)
      }  
    } else {
      list.of.characteristic.points[[i]] = list()
    }
  }
  return(list.of.characteristic.points);
}

BuildConstraintsForValuesEvaluation <- function(characteristic.points, 
                                                offsets.for.characteristic.points, levels.list, numOfVariables) {
  offsets.levels = ror:::getOffsets(levels.list);
  left.side.of.constraints = c()
  for (i in seq(1:length(levels.list))) {
    if (length(characteristic.points[[i]] > 0)) {  
      for (j in seq(1:length(levels.list[[i]]))) {
        interval <- GetIntervalData(characteristic.points, i, levels.list[[i]][j])
        proc = (levels.list[[i]][j] - interval$value[[1]]) / (interval$value[[2]] - interval$value[[1]])
        variables = vector(mode="numeric", length=numOfVariables)
        index.levels <- offsets.levels[i] + j - 1;
        index.interval <- offsets.for.characteristic.points[i] + interval$index - 1; 
        variables[index.levels] = -1;
        variables[index.interval-1] = 1-proc;
        variables[index.interval] = proc;
        left.side.of.constraints <- rbind(left.side.of.constraints, variables)
      }
    }
  }
  if (!is.null(left.side.of.constraints)) {
    return(list(lhs=left.side.of.constraints, dir=rep("==", nrow(left.side.of.constraints)), rhs=rep(0.0, nrow(left.side.of.constraints))))
  }
  return(list())
}

BuildMonotonousConstraints <- function(perf, strict.vf=FALSE, nums.of.characteristic.points) {
  
  stopifnot(is.logical(strict.vf));
  levels.list <- ror:::getLevels(perf);
  characteristic.points = GetCharacteristicPoints(perf, nums.of.characteristic.points);
  num.of.values = ror:::getNrVars(levels.list)
  num.of.characteristic.points = ror:::getNrVars(characteristic.points) -1;
  offsets <- ror:::getOffsets(levels.list);
  offsets.for.characteristic.points <- ror:::getOffsets(characteristic.points);
  offsets.for.characteristic.points = offsets.for.characteristic.points + num.of.values
  num.of.variables = num.of.values + num.of.characteristic.points;
  c1 <- BuildMonotonousConstraintsForCharacteristicPoints(
    perf, characteristic.points, offsets,  offsets.for.characteristic.points, num.of.variables, strict.vf)
  c2 <- BuildConstraintsForValuesEvaluation(
    characteristic.points, offsets.for.characteristic.points, levels.list, num.of.variables)
  monotonousConst <- ror:::combineConstraints(c1, c2)
  return(monotonousConst)
}

BuildPairwiseComparisonConstraints <- function(perf, strong.prefs = NULL,
                                               weak.prefs = NULL, indif.prefs = NULL) {
  #Function builds constraints for pairwise comparison 
  alt.vars <- ror:::buildAltVariableMatrix(perf)
  all.constraints <- NULL
  
  if (is.matrix(strong.prefs)) {
    for (i in 1:nrow(strong.prefs)) {
      pref.constraints <- ror:::buildStrongPreferenceConstraint(strong.prefs[i,1], strong.prefs[i,2], alt.vars)
      all.constraints <- ror:::combineConstraints(all.constraints, pref.constraints);
      row.names(all.constraints$lhs)[nrow(all.constraints$lhs)] <- c(paste(strong.prefs[i,1], ">", strong.prefs[i,2]))
    }
  }
  if (is.matrix(weak.prefs)) {
    for (i in 1:nrow(weak.prefs)) {
      pref.constraints <- ror:::buildWeakPreferenceConstraint(weak.prefs[i,1], weak.prefs[i,2], alt.vars)
      all.constraints <- ror:::combineConstraints(all.constraints, pref.constraints);
      row.names(all.constraints$lhs)[nrow(all.constraints$lhs)] <- c(paste(weak.prefs[i,1], ">=", weak.prefs[i,2]))
    }
  }
  if (is.matrix(indif.prefs)) {
    for (i in 1:nrow(indif.prefs)) {
      pref.constraints <- ror:::buildIndifPreferenceConstraint(indif.prefs[i,1], indif.prefs[i,2], alt.vars)
      all.constraints <- ror:::combineConstraints(all.constraints, pref.constraints);
      row.names(all.constraints$lhs)[nrow(all.constraints$lhs)] <- c(paste(indif.prefs[i,1], "==", indif.prefs[i,2]))
    }
  }
  return(all.constraints)  
}  
GenerateAllSubsetsOfConstraints <- function(all.constraints) {
  
  indexes <- lapply(1:nrow(all.constraints$lhs), function(x) combn(nrow(all.constraints$lhs),x))
  all.subsets.of.constraints = list()
  k = 0
  for (n in seq(length(indexes))) {
    for(j in seq(ncol(indexes[[n]]))) {
      constraints <- NULL
      for(i in seq(nrow(indexes[[n]]))) {
        constraints$lhs <- rbind(constraints$lhs, all.constraints$lhs[indexes[[n]][i,j], ])
        constraints$dir <- cbind(constraints$dir, all.constraints$dir[indexes[[n]][i,j]])
        constraints$rhs <- cbind(constraints$rhs, all.constraints$rhs[indexes[[n]][i,j]])
        row.names(constraints$lhs)[i] <- row.names(all.constraints$lhs)[indexes[[n]][i,j]] 
      }
      k <- k + 1 
      
      all.subsets.of.constraints[[k]] <- constraints
      
    }
  }
  
  return(all.subsets.of.constraints)
}

BuildMonotonousConstraintsForCharacteristicPoints <- function(perf,
                                                              characteristic.points, offsets,  offsets.for.characteristic.points, num.of.variables, strict.vf = FALSE) {
  levels.list <- ror:::getLevels(perf)  
  left.side.of.constraints = c()
  for (i in seq(1:length(characteristic.points))){  
    if (length(characteristic.points[[i]] > 0)) {
      for (j in seq(1:(length(characteristic.points[[i]])-1))) {
        index <- offsets.for.characteristic.points[i] + j - 1
        #list of variables, first numOfValues variables means g(x) and epsilon. Next numOfCharacteristicPoints
        #values means g(q_j) and epsilon
        piecewise.linear = vector(mode = "numeric", length = num.of.variables)
        piecewise.linear[index] <- 1
        piecewise.linear[index+1] <- -1
        if (strict.vf == TRUE) {
          nrVars <- ror:::getNrVars(levels.list);
          piecewise.linear[nrVars] = 1;
        }
        left.side.of.constraints <- rbind(left.side.of.constraints, piecewise.linear)
      }
    } else {
      res <- c()
      for (j in seq(1:(length(levels.list[[i]])-1))) {
        index <- offsets[i] + j - 1
        general <- array(0, dim=num.of.variables)
        general[index] <- 1
        general[index+1] <- -1
        if (strict.vf == TRUE) {
          nrVars <- ror:::getNrVars(levels.list);
          general[nrVars] = 1;
        }
        res <- rbind(res, general)
      }
      left.side.of.constraints <- rbind(left.side.of.constraints, res)
    }
  }
  result <- list(lhs=left.side.of.constraints, dir=rep("<=", nrow(left.side.of.constraints)),
                 rhs=rep(0.0, nrow(left.side.of.constraints)))
  return(result)
}

BuildBaseConstraints <- function( perf, num.of.variables, strict.vf = FALSE,  nums.of.characteristic.points = NULL) {
  
  # perf - performances matrix
  # num.of.variables - number of variables used in problem
  # strict.vf - TRUE -> value functions strictly increasing (instead of monotonous increasing)
  # is.piecewise.linear - Boolean variable determines whether functions should 
  #    be piecewise linear
  # nums.of.characteristicPoints - list of values - each value means num of characteristic points on criterium
  c = list(); #list of constraints
  c[[1]] <- BuildMonotonousConstraints(perf, strict.vf, nums.of.characteristic.points)  
  c[[2]] <- ror:::buildFirstLevelZeroConstraints(perf)
  c[[3]] <- ror:::buildBestLevelsAddToUnityConstraint(perf)
  c[[4]] <- ror:::buildAllVariablesLessThan1Constraint(perf)
  c[[5]] <- ror:::buildEpsilonStrictlyPositiveConstraint(perf)
  
  normalized.c = list();
  for (constraint in c) {
    constraint$lhs <- GetNormalizedMatrix(constraint$lhs, num.of.variables)
    normalized.c[length(normalized.c)+1] <- list(constraint)
  }
  
  base.constraints <- ror:::combineConstraints(normalized.c[[1]], normalized.c[[2]], normalized.c[[3]], normalized.c[[4]], normalized.c[[5]])
  return(base.constraints)
}

BuildBaseLPModel <- function(perf, strict.vf, strong.prefs = NULL,
                             weak.prefs = NULL, indif.prefs = NULL, nums.of.characteristic.points=NULL) {
  #Function build and return constraints for base LP model
  #  perf - performances matrix
  #  strict.vf - TRUE -> value functions strictly increasing (instead of monotonous increasing)
  #  *.prefs - an n x 2 matrix, where each row (a, b) means that a is 
  #    [strongly or weakly preferred, or indifferent] to b.
  #  nums.of.characteristicPoints - list of values - each value means num of characteristic points on criterium
  
  if (is.null(nums.of.characteristic.points)) {
    nums.of.characteristic.points = rep(0, ncol(perf))
  }
  
  num.of.variables = GetNumberOfVariables(perf, nums.of.characteristic.points)
  all.constraints <- BuildBaseConstraints(perf, num.of.variables, strict.vf, nums.of.characteristic.points)
  pairwise.comparison.constraints <- BuildPairwiseComparisonConstraints(perf,
                                                                        strong.prefs, weak.prefs, indif.prefs)
  
  if (!is.null(pairwise.comparison.constraints)) {
    pairwise.comparison.constraints$lhs <- GetNormalizedMatrix(pairwise.comparison.constraints$lhs, num.of.variables);  
    all.constraints <- ror:::combineConstraints(all.constraints, pairwise.comparison.constraints);    
  }
  
  return(all.constraints)
  
}



#############################################################################################
########################### Post Factum Analysis ############################################
#############################################################################################

PfaImproveVariant <- function(perf, a, q, which.attributes = NULL){
  
  if (is.null(which.attributes)){
    which.attributes = rep(1, ncol(perf))
  }
  for (i in seq(ncol(perf))) { 
    
    if (which.attributes[i] == 1) {
      max.val <- max(perf[,i])
      min.val <- min(perf[,i])
      perf[a, i] <- perf[a, i] * q
      
      if (perf[a, i] > max.val) {
        perf[a, i] <- max.val
      }
      if (perf[a, i] < min.val) {
        perf[a, i] <- min.val
      }
    }
  }
  return(perf)
}



BinarySearchForPossibleRanks <- function(perf, a, k, strict.vf, strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL,
                                         nums.of.characteristic.points=NULL, which.attributes = NULL, precision=0.00005, start.q) {
  print("BINARY SEARCH")
  number.of.binary.variables <- nrow(perf)
  q <- start.q
  diff <- q/2
  
  new.perf <- perf
  repeat {
    
    old.perf <- new.perf
    new.perf <- PfaImproveVariant(perf, a, q, which.attributes)
    
    constraints <- PfaRanksBuildConstraintsForMin(perf = new.perf, a =a, k =k, strict.vf = strict.vf,
                                                  strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                                  nums.of.characteristic.points=nums.of.characteristic.points)
    
    
    if (diff < precision) {
      
      new.perf <- PfaImproveVariant(perf, a, q + precision, which.attributes)
      
      constraints <-  PfaRanksBuildConstraintsForMin(perf = new.perf, a =a, k =k, strict.vf = strict.vf,
                                                     strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                                     nums.of.characteristic.points=nums.of.characteristic.points)
      
      print("NUMS OF CHARACTERISTIC POINTS")
      print(nums.of.characteristic.points)
      ShowSolutionForRanks(new.perf, constraints, number.of.binary.variables=number.of.binary.variables,  nums.of.characteristic.points = nums.of.characteristic.points)
      
      
      return(q + precision)
    }
    print(q)
    if (PfaRanksCheckConstraintsConsistency(new.perf, constraints, number.of.binary.variables)){       
      
      new.q <- q - (diff / 2)
      diff <- q - new.q
    } else { #infeasible - jeśli nieosiągalne, to musimy zwiększyć wartość 
      new.q <- q + (diff / 2)
      diff <- new.q - q
    }
    q <- new.q
  }
  
}

BinarySearchForNecessaryRanks <- function(perf, a, k, strict.vf, strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL,
                                          nums.of.characteristic.points=NULL, which.attributes = NULL, precision=0.00005, start.q) {
  print("BINARY SEARCH")
  number.of.binary.variables <- nrow(perf)
  q <- start.q
  diff <- q/2
  
  new.perf <- perf
  repeat {
    
    old.perf <- new.perf
    new.perf <- PfaImproveVariant(perf, a, q, which.attributes)
    
    constraints <- PfaRanksBuildConstraintsForMax(perf = new.perf, a =a, k =k, strict.vf = strict.vf,
                                                  strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                                  nums.of.characteristic.points=nums.of.characteristic.points)
    
    if (diff < precision) {
      
      new.perf <- PfaImproveVariant(perf, a, q - precision, which.attributes)
      
      constraints <- PfaRanksBuildConstraintsForMax(perf = new.perf, a =a, k =k, strict.vf = strict.vf,
                                                    strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                                    nums.of.characteristic.points=nums.of.characteristic.points)
      print("NUMS OF CHARACTERISTIC POINTS")
      print(nums.of.characteristic.points)
      ShowSolutionForRanks(new.perf, constraints, number.of.binary.variables= number.of.binary.variables,  nums.of.characteristic.points = nums.of.characteristic.points)
      
      
      return(q - precision)
    }
    
    if (!PfaRanksCheckConstraintsConsistency(new.perf, constraints, number.of.binary.variables)){       
      new.q <- q - (diff / 2)
      diff <- q - new.q
    } else { #feasible - jeśli osiągalne, to musimy zwiększyć wartość Necessary
      new.q <- q + (diff / 2)
      diff <- new.q - q
    }
    q <- new.q
  }
  
}



PfaRanksBuildConstraintsForMin<- function(perf, a, k, strict.vf = TRUE,
                                          strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL, 
                                          nums.of.characteristic.points=NULL) {
  #a - index of variant a
  #b - index of variant b
  
  alt.vars <- ror:::buildAltVariableMatrix(perf)
  base.model <- BuildBaseLPModel(perf, strict.vf=strict.vf, strong.prefs = strong.prefs,
                                 weak.prefs = weak.prefs, indif.prefs = indif.prefs, nums.of.characteristic.points=nums.of.characteristic.points)
  
  #if (!CheckConstraintsConsistency(perf, base.model)) {
  #  stop("Model infeasible - Podstawowe ograniczenia są nieosiągalne.")
  #} 
  pref.constraints <- list()
  number.of.variables <- nrow(perf) 
  all.binaries <- vector()
  
  binaries.sum.constraint <- list()
  binaries.sum.constraint$lhs <- matrix(ncol= (number.of.variables), nrow =  1, data=0)
  binaries.sum.constraint$dir <- "<="
  binaries.sum.constraint$rhs <- k - 1
  
  for (b in seq(number.of.variables)) {
    if (b != a) {
      a.pref.b <- list()
      binaries <- matrix(ncol= (number.of.variables), nrow =  1, data=0)
      binaries[b] <- M_BIG_NUMBER
      all.binaries <- rbind(all.binaries, binaries)
      binaries.sum.constraint$lhs[b] <- 1
      # Weak
      a.pref.b <- ror:::buildWeakPreferenceConstraint(a, b, alt.vars) 
      pref.constraints <- ror:::combineConstraints(pref.constraints, a.pref.b)
    }
  }
  
  if (ncol(base.model$lhs) >= ncol(pref.constraints$lhs)) {
    pref.constraints$lhs <- GetNormalizedMatrix(pref.constraints$lhs, ncol(base.model$lhs))  
  } else {
    base.model$lhs <- GetNormalizedMatrix(base.model$lhs, ncol(pref.constraints$lhs))
  }
  pref.constraints$lhs <- cbind(pref.constraints$lhs, all.binaries)
  pref.constraints$rhs[1:length(pref.constraints$rhs)] <- 0
  
  base.model$lhs <- GetNormalizedMatrix(base.model$lhs, ncol(pref.constraints$lhs))
  binaries.sum.constraint$lhs <- GetNormalizedMatrix(binaries.sum.constraint$lhs, ncol(pref.constraints$lhs), right=FALSE)
  all.constraints <-ror:::combineConstraints(base.model, pref.constraints, binaries.sum.constraint)
  return(all.constraints)
}


PfaRanksBuildConstraintsForMax <- function(perf, a, k, strict.vf = TRUE,
                                           strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL, nums.of.characteristic.points=NULL) {
  #a - index of variant a
  #b - index of variant b
  
  alt.vars <- ror:::buildAltVariableMatrix(perf)
  base.model <- BuildBaseLPModel(perf, strict.vf, strong.prefs = strong.prefs,
                                 weak.prefs = weak.prefs, indif.prefs = indif.prefs,
                                 nums.of.characteristic.points=nums.of.characteristic.points)
  
  #if (!CheckConstraintsConsistency(perf, base.model)) {
  #  stop("Model infeasible")
  #}
  pref.constraints <- list()
  number.of.variables <- nrow(perf) 
  all.binaries <- vector()
  
  binaries.sum.constraint <- list()
  binaries.sum.constraint$lhs = matrix(ncol= (number.of.variables), nrow =  1, data=0)
  binaries.sum.constraint$dir = ">="
  binaries.sum.constraint$rhs = k
  
  for (b in seq(number.of.variables)) {
    if (b != a) {
      a.pref.b <- list()
      binaries <- matrix(ncol= (number.of.variables), nrow =  1, data=0)
      binaries[b] <- M_BIG_NUMBER
      all.binaries <- rbind(all.binaries, binaries)
      binaries.sum.constraint$lhs[b] <- 1
      a.pref.b <- ror:::buildWeakPreferenceConstraint(a, b, alt.vars) 
      a.pref.b$lhs[length(a.pref.b$lhs)] <- 1
      a.pref.b$dir <- "<="
      
      pref.constraints <- ror:::combineConstraints(pref.constraints, a.pref.b)
    }
  }
  
  if (ncol(base.model$lhs) >= ncol(pref.constraints$lhs)) {
    pref.constraints$lhs <- GetNormalizedMatrix(pref.constraints$lhs, ncol(base.model$lhs))  
  } else {
    base.model$lhs <- GetNormalizedMatrix(base.model$lhs, ncol(pref.constraints$lhs))
  }
  pref.constraints$lhs <- cbind(pref.constraints$lhs, all.binaries)
  pref.constraints$rhs[1:length(pref.constraints$rhs)] <- M_BIG_NUMBER
  base.model$lhs <- GetNormalizedMatrix(base.model$lhs, ncol(pref.constraints$lhs))
  binaries.sum.constraint$lhs <- GetNormalizedMatrix(binaries.sum.constraint$lhs, ncol(pref.constraints$lhs), right=FALSE)
  all.constraints <-ror:::combineConstraints(base.model, pref.constraints, binaries.sum.constraint)
  return(all.constraints)
}



PfaFindQPossibleRankImprovement <- function(perf, a, k, strict.vf,
                                            strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL,
                                            nums.of.characteristic.points=NULL, which.attributes = NULL, precision=0.00005) {
  number.of.binary.variables <- nrow(perf)
  new.perf = perf
  q <- 1
  stop.searching = FALSE
  
  repeat {    
    old.perf = new.perf
    
    new.perf <- PfaImproveVariant(perf, a, q, which.attributes = which.attributes)
    if (EqualMatrix(new.perf, old.perf) && q > 1) {
      stop.searching = TRUE  
    }
    
    constraints <- PfaRanksBuildConstraintsForMin(perf = new.perf, a =a, k =k, strict.vf = strict.vf,
                                                  strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                                  nums.of.characteristic.points=nums.of.characteristic.points)
    
    
    
    result <- PfaRanksCheckConstraintsConsistency(new.perf, constraints, number.of.binary.variables)
    
    if (result) {
      if (q == 1) {
        constraints <- PfaRanksBuildConstraintsForMin(perf = new.perf, a =a, k =k, strict.vf = strict.vf,
                                                      strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                                      nums.of.characteristic.points=nums.of.characteristic.points)
        
        
        ShowSolutionForRanks(new.perf, constraints, number.of.binary.variables=number.of.binary.variables,  nums.of.characteristic.points = nums.of.characteristic.points)
        
        return(1)
      }
      q <- BinarySearchForPossibleRanks(perf = perf, a = a, k = k, strict.vf = strict.vf,
                                        strong.prefs = strong.prefs, weak.prefs = weak.prefs, indif.prefs = indif.prefs,
                                        nums.of.characteristic.points = nums.of.characteristic.points,
                                        which.attributes = which.attributes, precision=precision, start.q = q)
      return(q);#binary search
    } else {
      if (stop.searching) {
        
        print("STOP")
        stop("model is infeasible")
      }
    }
    
    
    q <- q * 2
  }
}


PfaFindQPossibleRankDeterioration <- function(perf, a, k, strict.vf,
                                              strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL,
                                              nums.of.characteristic.points=NULL, which.attributes = NULL, precision=0.00005) {
  number.of.binary.variables <- nrow(perf)
  new.perf = perf
  q <- 1
  stop.searching = FALSE
  
  repeat { 
    print(q)
    old.perf = new.perf
    new.perf <- PfaImproveVariant(perf, a, q, which.attributes = which.attributes)
    if (EqualMatrix(new.perf, old.perf) && q < 1) {
      stop.searching = TRUE  
    }
    constraints <- PfaRanksBuildConstraintsForMin(perf = new.perf, a =a, k =k, strict.vf = strict.vf,
                                                  strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                                  nums.of.characteristic.points=nums.of.characteristic.points)
    
    
    result <- PfaRanksCheckConstraintsConsistency(new.perf, constraints, number.of.binary.variables)
    if (!result) {
      if (q == 1) {
        
        stop("model infeasible ")
      }
      q <- BinarySearchForPossibleRanks(perf = perf, a = a, k = k, strict.vf = strict.vf,
                                        strong.prefs = strong.prefs, weak.prefs = weak.prefs, indif.prefs = indif.prefs,
                                        nums.of.characteristic.points = nums.of.characteristic.points,
                                        which.attributes = which.attributes, precision=precision, start.q = q)
      return(q);#binary search
    } else {
      
      if (stop.searching) {
        constraints <-  PfaRanksBuildConstraintsForMin(perf = new.perf, a =a, k =k, strict.vf = strict.vf,
                                                       strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                                       nums.of.characteristic.points=nums.of.characteristic.points)
        
        print("NUMS OF CHARACTERISTIC POINTS")
        print(nums.of.characteristic.points)
        ShowSolutionForRanks(new.perf, constraints, number.of.binary.variables = number.of.binary.variables, nums.of.characteristic.points = nums.of.characteristic.points)
        
        
        return(
          FindQLimesRanks( perf, a, k=k,  start.q = q, is.possible=TRUE, strict.vf=strict.vf, which.attributes = which.attributes,
                           direction.up = 1, precision=precision, nums.of.characteristic.points=nums.of.characteristic.points, max=FALSE)) ##search limes
      }
      print("STOP-FALSE")
    }
    
    
    q <- q * 0.5
  }
}

FindQLimesRanks <- function(perf, a, k,  start.q, is.possible=FALSE, strict.vf,
                            strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL,
                            nums.of.characteristic.points=NULL, which.attributes = NULL, direction.up = 1, precision=0.00005, max=FALSE){
  
  q <- start.q
  old.perf = perf
  new.perf = perf
  diff <- q
  #diff <- q / 2
  #
  # Stara wersja
  #q <- q + (direction.up * diff/2)
  #diff <- diff /2
  #Nowa wersja
  q <- q + (direction.up * diff)
  diff <- diff /2
  
  new.perf <- PfaImproveVariant(perf, a, q, which.attributes = which.attributes)
  old.perf <- new.perf
  repeat {
    print(q)
    print(diff)
    print("diff:")
    
    #print(new.perf)
    if (diff < precision) {
      constraints = list()
      if (max == TRUE) {
        constraints <- PfaRanksBuildConstraintsForMax(perf = new.perf, a =a, k =k, strict.vf = strict.vf,
                                                      strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                                      nums.of.characteristic.points=nums.of.characteristic.points)
      } else {
        constraints <- PfaRanksBuildConstraintsForMax(perf = new.perf, a =a, k =k, strict.vf = strict.vf,
                                                      strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                                      nums.of.characteristic.points=nums.of.characteristic.points)
        
      }
      
      ShowSolution(new.perf, constraints, nums.of.characteristic.points = nums.of.characteristic.points)
      
      return(q)
    }
    
    
    if (EqualMatrix(new.perf, old.perf)) {
      new.q <- q + (direction.up * diff/2)
      diff <- diff /2
    } else {
      
      new.q <- q - (direction.up * diff/2)
      diff <- diff /2
    }
    q <- new.q
    old.perf <- new.perf
    new.perf <- PfaImproveVariant(perf, a, q, which.attributes = which.attributes)
    # print(new.perf)
    
  }
} 


PfaFindQNecessaryRankImprovement <- function(perf, a, k, strict.vf,
                                             strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL,
                                             nums.of.characteristic.points=NULL, which.attributes = NULL, precision=0.00005) {
  
  number.of.binary.variables <- nrow(perf)
  new.perf = perf
  q <- 1
  stop.searching = FALSE
  
  repeat {    
    old.perf = new.perf
    new.perf <- PfaImproveVariant(perf, a, q, which.attributes = which.attributes)
    if (EqualMatrix(new.perf, old.perf) && q > 1) {
      stop.searching = TRUE  
    }
    
    constraints <- PfaRanksBuildConstraintsForMax(perf = new.perf, a =a, k =k, strict.vf = strict.vf,
                                                  strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                                  nums.of.characteristic.points=nums.of.characteristic.points)
    
    result <- PfaRanksCheckConstraintsConsistency(new.perf, constraints, number.of.binary.variables)
    print(q)
    if (!result) {
      if (q == 1) {
        stop("infeasible")
      }
      q <- BinarySearchForNecessaryRanks(perf = perf, a = a, k = k, strict.vf = strict.vf,
                                         strong.prefs = strong.prefs, weak.prefs = weak.prefs, indif.prefs = indif.prefs,
                                         nums.of.characteristic.points = nums.of.characteristic.points,
                                         which.attributes = which.attributes, precision=precision, start.q = q)
      return(q);#binary search
    } else {
      if (stop.searching) {
        #wydaje mi się, że tutaj powinno być infeasible, gdyż a nigdy nie osiągnie szukanego stanu
        stop("infeasible")
        
        # lub jeśli mimo wszystko powinno się szukać granicy. 
        #return(
        #  FindQLimes( perf, a, start.q = q, is.possible=TRUE, which.attributes = which.attributes,
        #              direction.up = -1, precision=precision)) ##search limes
      }
    }
    q <- q * 2
  }
}


PfaFindQNecessaryRankDeterioration <- function(perf, a, k, strict.vf,
                                               strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL,
                                               nums.of.characteristic.points=NULL, which.attributes = NULL, precision=0.00005) {
  number.of.binary.variables <- nrow(perf)
  new.perf = perf
  q <- 1
  stop.searching = FALSE
  
  repeat {    
    old.perf = new.perf
    new.perf <- PfaImproveVariant(perf, a, q, which.attributes = which.attributes)
    if (EqualMatrix(new.perf, old.perf) && q < 1) {
      
      stop.searching = TRUE  
    }
    
    constraints <- PfaRanksBuildConstraintsForMax(perf = new.perf, a =a, k =k, strict.vf = strict.vf,
                                                  strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                                  nums.of.characteristic.points=nums.of.characteristic.points)
    
    
    result <- PfaRanksCheckConstraintsConsistency(new.perf, constraints, number.of.binary.variables)
    if (result) {
      if (q == 1) {
        return(1)
      }
      q <- BinarySearchForNecessaryRanks(perf = perf, a = a, k = k, strict.vf = strict.vf,
                                         strong.prefs = strong.prefs, weak.prefs = weak.prefs, indif.prefs = indif.prefs,
                                         nums.of.characteristic.points = nums.of.characteristic.points,
                                         which.attributes = which.attributes, precision=precision, start.q = q)
      return(q);#binary search
    } else {
      
      if (stop.searching) {
        #wydaje mi się, że powinno być infeasible
        #stop("infeasible. Równianie nigdy nie będzie spełnione.")
        #JEśli jednak szukać bounda.
        constraints <- PfaRanksBuildConstraintsForMax(perf = new.perf, a =a, k =k, strict.vf = strict.vf,
                                                      strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                                      nums.of.characteristic.points=nums.of.characteristic.points)
        
        ShowSolutionForRanks(new.perf, constraints, number.of.binary.variables= number.of.binary.variables,  nums.of.characteristic.points = nums.of.characteristic.points)
        
        
        return(
          FindQLimesRanks( perf, a, k=k, start.q = q, is.possible=FALSE, strict.vf=strict.vf, which.attributes = which.attributes,
                           direction.up = 1, precision=precision, nums.of.characteristic.points=nums.of.characteristic.points ,max = TRUE)) ##search limes
      }
    }
    
    q <- q * 0.5
  }
}


###--------


PfaRanksBuildConstraints <- function(perf, a, k, is.possible=TRUE, strict.vf = FALSE,
                                     strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL, is.piecewise.linear = FALSE,
                                     nums.of.characteristic.points=NULL) {
  #a - index of variant a
  #b - index of variant b
  
  alt.vars <- ror:::buildAltVariableMatrix(perf)
  base.model <- BuildBaseLPModel(perf, strict.vf, strong.prefs = strong.prefs,
                                 weak.prefs = weak.prefs, indif.prefs = indif.prefs, is.piecewise.linear = is.piecewise.linear,
                                 nums.of.characteristic.points=nums.of.characteristic.points)
  
  if (!CheckConstraintsConsistency(perf, base.model)) {
    stop("Model infeasible")
  }
  pref.constraints <- list()
  number.of.variables <- nrow(perf) 
  all.binaries <- vector()
  
  binaries.sum.constraint <- list()
  binaries.sum.constraint$lhs = matrix(ncol= (number.of.variables), nrow =  1, data=0)
  
  if (is.possible) {
    binaries.sum.constraint$dir = "<="
    binaries.sum.constraint$rhs = k - 1
  } else {
    binaries.sum.constraint$dir = ">="
    binaries.sum.constraint$rhs = k
  }
  
  for (b in seq(number.of.variables)) {
    if (b != a) {
      a.pref.b <- list()
      binaries <- matrix(ncol= (number.of.variables), nrow =  1, data=0)
      binaries[b] <- M_BIG_NUMBER
      all.binaries <- rbind(all.binaries, binaries)
      binaries.sum.constraint$lhs[b] <- 1
      # Week
      a.pref.b <- ror:::buildStrongPreferenceConstraint(a, b, alt.vars) 
      a.pref.b$dir <- "<="
      
      pref.constraints <- ror:::combineConstraints(pref.constraints, a.pref.b)
    }
  }
  
  if (ncol(base.model$lhs) >= ncol(pref.constraints$lhs)) {
    pref.constraints$lhs <- GetNormalizedMatrix(pref.constraints$lhs, ncol(base.model$lhs))  
  } else {
    base.model$lhs <- GetNormalizedMatrix(base.model$lhs, ncol(pref.constraints$lhs))
  }
  pref.constraints$lhs <- cbind(pref.constraints$lhs, all.binaries)
  pref.constraints$rhs[1:length(pref.constraints$rhs)] <- M_BIG_NUMBER
  base.model$lhs <- GetNormalizedMatrix(base.model$lhs, ncol(pref.constraints$lhs))
  binaries.sum.constraint$lhs <- GetNormalizedMatrix(binaries.sum.constraint$lhs, ncol(pref.constraints$lhs), right=FALSE)
  all.constraints <-ror:::combineConstraints(base.model, pref.constraints, binaries.sum.constraint)
  return(all.constraints)
}


PfaRanksIsBestVariant <- function(perf, a) {
  is.best <- TRUE
  for(i in seq(ncol(perf))) {
    if (perf[a, i] < max(perf[,i])) {
      return (FALSE)
    }
  }
  return(TRUE)
  
}

PfaRanksFindQ <- function(perf, a, k, is.possible=FALSE, strict.vf,
                          strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL, is.piecewise.linear = FALSE,
                          nums.of.characteristic.points=NULL) {
  
  q <- 0.5;
  number.of.binary.variables <- nrow(perf)
  repeat{
    q <- 2*q
    new.perf <- PfaImproveVariant(perf, a, q)
    
    
    
    constraints <- PfaRanksBuildConstraints(perf = new.perf, a = a, k = k, is.possible = is.possible, strict.vf = strict.vf,
                                            strong.prefs = strong.prefs, weak.prefs = weak.prefs, indif.prefs = indif.prefs, is.piecewise.linear = is.piecewise.linear,
                                            nums.of.characteristic.points = nums.of.characteristic.points)
    
    result <- PfaRanksCheckConstraintsConsistency(new.perf, constraints, number.of.binary.variables)
    if ((is.possible) && (result)){
      break
    } 
    if ((!is.possible) &&(!result)) {
      if (q == 1) {
        stop("problem is infeasible")
      }
      break
    }
    if (PfaRanksIsBestVariant(new.perf, a)) {
      stop("Model infeasible")
    }
    
  }
  return(q)
}


#Rank related improvement - possible 
PfaRanksPossibleComprehensiveImprovement <- function(perf, a, k, strict.vf,
                                                     strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL, is.piecewise.linear = FALSE,
                                                     nums.of.characteristic.points=NULL, precision=0.05, greater.than.one=TRUE) {
  
  is.possible=TRUE
  
  q <- 1
  if (greater.than.one) {
    q <- PfaRanksFindQ(perf = perf, a = a, k = k,  is.possible=is.possible, strict.vf=strict.vf,
                       strong.prefs=strong.prefs, weak.prefs=weak.prefs, indif.prefs = indif.prefs, is.piecewise.linear = is.piecewise.linear,
                       nums.of.characteristic.points=nums.of.characteristic.points)
  }
  ## TODO
  number.of.binary.variables <- nrow(perf)
  diff <- q
  
  repeat{
    new.perf <- PfaImproveVariant(perf, a, q)
    constraints <- PfaRanksBuildConstraints(perf = new.perf, a = a, k = k, is.possible = is.possible, strict.vf = strict.vf,
                                            strong.prefs = strong.prefs, weak.prefs = weak.prefs, indif.prefs = indif.prefs, is.piecewise.linear = is.piecewise.linear,
                                            nums.of.characteristic.points = nums.of.characteristic.points)
    
    
    if (diff < precision) {
      return(q)
    }
    if (PfaRanksCheckConstraintsConsistency(new.perf, constraints,  number.of.binary.variables)){ 
      new.q <- q - (diff / 2)
      diff <- q - new.q
    } else { #infeasible 
      new.q <- q + (diff / 2)
      diff <- new.q - q
    }
    q <- new.q
    
    if (greater.than.one) {
      if (q < 1) {
        return(1)
      }
    } else {
      if (q > 1) {
        stop("ifnfeasible: q is greater than one. ")
      }
    }
  }
}

#Rank related improvement - necessary
PfaRanksNecessaryComprehensiveImprovement <- function(perf, a, k, strict.vf,
                                                      strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL, is.piecewise.linear = FALSE,
                                                      nums.of.characteristic.points=NULL, precision=0.05, greater.than.one = TRUE) {
  
  is.possible=FALSE
  q <- 0.5
  if (greater.than.one) {
    q <- PfaRanksFindQ(perf = perf, a = a, k = k,  is.possible=is.possible, strict.vf=strict.vf,
                       strong.prefs=strong.prefs, weak.prefs=weak.prefs, indif.prefs = indif.prefs, is.piecewise.linear = is.piecewise.linear,
                       nums.of.characteristic.points=nums.of.characteristic.points)
  }
  ## TODO
  number.of.binary.variables <- nrow(perf)
  diff <- q
  
  repeat{
    new.perf <- PfaImproveVariant(perf, a, q)
    constraints <- PfaRanksBuildConstraints(perf = new.perf, a = a, k = k, is.possible = is.possible, strict.vf = strict.vf,
                                            strong.prefs = strong.prefs, weak.prefs = weak.prefs, indif.prefs = indif.prefs, is.piecewise.linear = is.piecewise.linear,
                                            nums.of.characteristic.points = nums.of.characteristic.points)
    
    
    if (diff < precision) {
      return(q)
    }
    if (PfaRanksCheckConstraintsConsistency(new.perf, constraints,  number.of.binary.variables)){ 
      new.q <- q + (diff / 2)
      diff <- new.q - q
    } else { #infeasible 
      new.q <- q - (diff / 2)
      diff <- q - new.q
    }
    q <- new.q
    
    if (greater.than.one) {
      if (q < 1) {
        stop("Infeasible: q is smaller than one")
      }
    } else {
      if (q > 1) {
        return(1)
      }
    }
  }
}

PfaRanksSolveModel <- function(performances, model, number.of.binary.variables) {
  number.of.variables <- ncol(model$lhs)
  position.of.first.binary.variable <- number.of.variables - number.of.binary.variables + 1
  types <- rep('C', number.of.variables)
  types[position.of.first.binary.variable : number.of.variables] <- 'B'
  
  objective <- ror:::buildObjectiveFunction(performances)
  objective <- GetNormalizedMatrix(matrix=objective, width=ncol(model$lhs))
  obj <- L_objective(objective)
  roiConst <- L_constraint(L = model$lhs, dir =model$dir, rhs=model$rhs)
  lp <- OP(objective=obj, constraints=roiConst, maximum=TRUE, types=types)
  res <- ROI_solve(lp, ror:::.solver)
  
  return(res)
}


PfaRanksCheckConstraintsConsistency <- function(perf, model, number.of.binary.variables = 0) {
  #check constraints consistency  
  #  model: structure of constraints with the following elements:
  #    model$lhs: matrix - left side of constraints
  #    model$dir: list of operators
  #    model$rhs: matrix - right side of constraints
  
  #baseModel <- BuildBaseLPModel(perf, strict.vf, strong.prefs = NULL,
  #                           weak.prefs = NULL, indif.prefs = NULL, is.piecewise.linear = FALSE,
  #                           nums.of.characteristicPoints=NULL)
  
  
  ret <- PfaRanksSolveModel(perf, model, number.of.binary.variables)
  
  return(ret$status$code == 0 && ret$objval >= ror:::MINEPS)
}






#PfaRanksBuildConstraintsForUmissing <- function(perf, a, b, is.possible=TRUE, strict.vf = FALSE,
#                                           strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL, is.piecewise.linear = FALSE,
#                                           nums.of.characteristic.points=NULL, improvement = TRUE) {
#a: index of variant a
#improvement: bool  TRUE/FALSE(deterioration)

#  alt.vars <- ror:::buildAltVariableMatrix(perf)
#  base.model <- BuildBaseLPModel(perf, strict.vf, strong.prefs = strong.prefs,
#                                 weak.prefs = weak.prefs, indif.prefs = indif.prefs, is.piecewise.linear = is.piecewise.linear,
#                                 nums.of.characteristic.points=nums.of.characteristic.points)

#  if (!CheckConstraintsConsistency(perf, base.model)) {
#    stop("Model infeasible")
#  }

#  pref.constraints <- list()
#  number.of.variables <- nrow(perf) 
#  all.binaries <- vector()

#  binaries.sum.constraint <- list()
#  binaries.sum.constraint$lhs = matrix(ncol= (number.of.variables), nrow =  1, data=0)
#  
#  if (is.possible) {
#    binaries.sum.constraint$dir = "<="
#    binaries.sum.constraint$rhs = k - 1
#  } else {
#    binaries.sum.constraint$dir = ">="
#    binaries.sum.constraint$rhs = k
#  }
#  
#  for (b in seq(number.of.variables)) {
#    if (b != a) {
#      a.pref.b <- list()
#      binaries <- matrix(ncol= (number.of.variables), nrow =  1, data=0)
#      binaries[b] <- M_BIG_NUMBER
#      all.binaries <- rbind(all.binaries, binaries)
#      binaries.sum.constraint$lhs[b] <- 1
#      a.pref.b <- ror:::buildWeakPreferenceConstraint(a, b, alt.vars) 
#      a.pref.b$dir <- "<="
#      
#      if (ncol(base.model$lhs) >= ncol(a.pref.b$lhs)) {
#        a.pref.b$lhs <- GetNormalizedMatrix(a.pref.b$lhs, ncol(base.model$lhs))  
#      } else {
#        base.model$lhs <- GetNormalizedMatrix(base.model$lhs, ncol(a.pref.b$lhs))
#      }

#      if (((is.possible) && (improvement)) || ((!is.possible) && (!improvement))) {
#        a.pref.b$lhs <- cbind(a.pref.b$lhs, 1)  
#      } else {
#        a.pref.b$lhs <- cbind(a.pref.b$lhs, -1)
#      }
#      


#      pref.constraints <- ror:::combineConstraints(pref.constraints, a.pref.b)
#    }
#  }



#  base.model$lhs <- GetNormalizedMatrix(base.model$lhs, ncol(a.pref.b$lhs))
#  greater.than.zero <- list(lhs = rep(0.0, ncol(base.model$lhs)), dir=">=", rhs=0)
#  greater.than.zero$lhs[ncol(base.model$lhs)] <- 1
#  
#  all.constraints <-ror:::combineConstraints(base.model, a.pref.b, greater.than.zero)
#  
#  return(all.constraints)
#}



PfaRanksFindUmissingForPossiblyOrNecessaryPrefference <- function(perf, a, k, strict.vf, is.possibly.preffered,
                                                                  strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL, is.piecewise.linear = FALSE,
                                                                  nums.of.characteristic.points=NULL, precision=0.05, improvement=TRUE) {
  
  constraints <- PfaRanksBuildConstraintsForUmissing(perf = perf, a = a, b = k, is.possible =  is.possibly.preffered, strict.vf = strict.vf,
                                                     strong.prefs = strong.prefs, weak.prefs = weak.prefs, indif.prefs = indif.prefs, is.piecewise.linear = is.piecewise.linear,
                                                     nums.of.characteristic.points = nums.of.characteristic.points, improvement=improvement)
  
  maximum = TRUE
  if (((is.possibly.preffered) && (improvement)) || ((!is.possibly.preffered) && (!improvement))) { 
    maximum = FALSE
  } else {
    maximum = TRUE
  }
  
  ret <- PfaUmissingSolveModel(perf, constraints,  maximum = maximum)
  
  
  return(ret$objval)
}



#based on MIN
PfaRanksBuildConstraintsForPossibly <- function(perf, a, k, strict.vf = TRUE,
                                                strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL, 
                                                nums.of.characteristic.points=NULL, improvement = FALSE) {
  #a - index of variant a
  #b - index of variant b
  
  alt.vars <- ror:::buildAltVariableMatrix(perf)
  base.model <- BuildBaseLPModel(perf, strict.vf=strict.vf, strong.prefs = strong.prefs,
                                 weak.prefs = weak.prefs, indif.prefs = indif.prefs, nums.of.characteristic.points=nums.of.characteristic.points)
  
  #if (!CheckConstraintsConsistency(perf, base.model)) {
  #  stop("Model infeasible - Podstawowe ograniczenia są nieosiągalne.")
  #} 
  pref.constraints <- list()
  number.of.variables <- nrow(perf) 
  all.binaries <- vector()
  
  binaries.sum.constraint <- list()
  binaries.sum.constraint$lhs <- matrix(ncol= (number.of.variables), nrow =  1, data=0)
  binaries.sum.constraint$dir <- "<="
  binaries.sum.constraint$rhs <- k - 1
  u <- 1
  if (improvement) {
    u <- 1
  } else {
    u <- -1
  }
  position.of.u = 0
  for (b in seq(number.of.variables)) {
    if (b != a) {
      a.pref.b <- list()
      binaries <- matrix(ncol= (number.of.variables), nrow =  1, data=0)
      binaries[b] <- M_BIG_NUMBER
      all.binaries <- rbind(all.binaries, binaries)
      binaries.sum.constraint$lhs[b] <- 1
      # Weak
      a.pref.b <- ror:::buildWeakPreferenceConstraint(a, b, alt.vars) 
      
      pref.constraints <- ror:::combineConstraints(pref.constraints, a.pref.b)
    }
  }
  
  if (ncol(base.model$lhs) >= ncol(pref.constraints$lhs)) {
    pref.constraints$lhs <- GetNormalizedMatrix(pref.constraints$lhs, ncol(base.model$lhs))  
  } else {
    base.model$lhs <- GetNormalizedMatrix(base.model$lhs, ncol(pref.constraints$lhs))
  }
  pref.constraints$lhs <- cbind(pref.constraints$lhs, u) 
  position.of.u = ncol(pref.constraints$lhs)
  pref.constraints$lhs <- cbind(pref.constraints$lhs, all.binaries)
  
  #constraints u >= 0
  u.greater <- list(
    lhs = rep(0, ncol(pref.constraints$lhs)),
    dir = ">=",
    rhs = 0
  )
  u.greater$lhs[position.of.u] = 1
  pref.constraints$rhs[1:length(pref.constraints$rhs)] <- 0
  pref.constraints <- ror:::combineConstraints(pref.constraints, u.greater)
  base.model$lhs <- GetNormalizedMatrix(base.model$lhs, ncol(pref.constraints$lhs))
  binaries.sum.constraint$lhs <- GetNormalizedMatrix(binaries.sum.constraint$lhs, ncol(pref.constraints$lhs), right=FALSE)
  all.constraints <-ror:::combineConstraints(base.model, pref.constraints, binaries.sum.constraint)
  return(all.constraints)
}


#based on MAX
PfaRanksBuildConstraintsForNecessary <- function(perf, a, k, strict.vf = TRUE,
                                                 strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL, nums.of.characteristic.points=NULL, improvement = FALSE) {
  #a - index of variant a
  #b - index of variant b
  
  alt.vars <- ror:::buildAltVariableMatrix(perf)
  base.model <- BuildBaseLPModel(perf, strict.vf, strong.prefs = strong.prefs,
                                 weak.prefs = weak.prefs, indif.prefs = indif.prefs,
                                 nums.of.characteristic.points=nums.of.characteristic.points)
  
  #if (!CheckConstraintsConsistency(perf, base.model)) {
  #  stop("Model infeasible")
  #}
  pref.constraints <- list()
  number.of.variables <- nrow(perf) 
  all.binaries <- vector()
  
  binaries.sum.constraint <- list()
  binaries.sum.constraint$lhs = matrix(ncol= (number.of.variables), nrow =  1, data=0)
  binaries.sum.constraint$dir = ">="
  binaries.sum.constraint$rhs = k
  
  u <- 1
  if (improvement) {
    u <- 1
  } else {
    u <- -1
  }
  position.of.u = 0
  
  for (b in seq(number.of.variables)) {
    if (b != a) {
      a.pref.b <- list()
      binaries <- matrix(ncol= (number.of.variables), nrow =  1, data=0)
      binaries[b] <- M_BIG_NUMBER
      all.binaries <- rbind(all.binaries, binaries)
      binaries.sum.constraint$lhs[b] <- 1
      a.pref.b <- ror:::buildWeakPreferenceConstraint(a, b, alt.vars) 
      a.pref.b$lhs[length(a.pref.b$lhs)] <- 1
      a.pref.b$dir <- "<="
      
      pref.constraints <- ror:::combineConstraints(pref.constraints, a.pref.b)
    }
  }
  
  if (ncol(base.model$lhs) >= ncol(pref.constraints$lhs)) {
    pref.constraints$lhs <- GetNormalizedMatrix(pref.constraints$lhs, ncol(base.model$lhs))  
  } else {
    base.model$lhs <- GetNormalizedMatrix(base.model$lhs, ncol(pref.constraints$lhs))
  }
  
  pref.constraints$lhs <- cbind(pref.constraints$lhs, u) 
  position.of.u = ncol(pref.constraints$lhs)
  
  
  pref.constraints$lhs <- cbind(pref.constraints$lhs, all.binaries)
  pref.constraints$rhs[1:length(pref.constraints$rhs)] <- M_BIG_NUMBER
  
  #constraints u >= 0
  u.greater <- list(
    lhs = rep(0, ncol(pref.constraints$lhs)),
    dir = ">=",
    rhs = 0
  )
  u.greater$lhs[position.of.u] = 1
  pref.constraints <- ror:::combineConstraints(pref.constraints, u.greater)
  base.model$lhs <- GetNormalizedMatrix(base.model$lhs, ncol(pref.constraints$lhs))
  
  
  binaries.sum.constraint$lhs <- GetNormalizedMatrix(binaries.sum.constraint$lhs, ncol(pref.constraints$lhs), right=FALSE)
  all.constraints <-ror:::combineConstraints(base.model, pref.constraints, binaries.sum.constraint)
  return(all.constraints)
}

#PfaRanksU <- function(perf, model, number.of.binary.variables = 0) {
#check constraints consistency  
#  model: structure of constraints with the following elements:
#    model$lhs: matrix - left side of constraints
#    model$dir: list of operators
#    model$rhs: matrix - right side of constraints

#baseModel <- BuildBaseLPModel(perf, strict.vf, strong.prefs = NULL,
#                           weak.prefs = NULL, indif.prefs = NULL, is.piecewise.linear = FALSE,
#                           nums.of.characteristicPoints=NULL)


#  ret <- PfaRanksSolveModel(perf, model, number.of.binary.variables)
#  print(ret$solution)
#  return(ret$status$code == 0 && ret$objval >= ror:::MINEPS)
#}

PfaRanksSolveModelForU <- function(performances, model, maximum) {
  number.of.binary.variables <- nrow(performances)
  number.of.variables <- ncol(model$lhs)
  position.of.first.binary.variable <- number.of.variables - number.of.binary.variables + 1
  position.of.u <- number.of.variables - number.of.binary.variables
  types <- rep('C', number.of.variables)
  types[position.of.first.binary.variable : number.of.variables] <- 'B'
  
  objective <- rep(0, ncol(model$lhs))
  objective[position.of.u] <- 1
  objective <- GetNormalizedMatrix(matrix=objective, width=ncol(model$lhs))
  obj <- L_objective(objective)
  roiConst <- L_constraint(L = model$lhs, dir =model$dir, rhs=model$rhs)
  lp <- OP(objective=obj, constraints=roiConst, maximum=maximum, types=types)
  res <- ROI_solve(lp, ror:::.solver)
  #print("MAXIMUM")
  #print("OBJECTIVE:")
  #print(obj)
  return(res)
}





#############################################################################################
########################### Rank Related Eval Modifying ###############################
#############################################################################################



# RANK RELATED IMPROVEMENT - POST FACTUM ANALYSIS - Powiązane funkcje, znajdują się w pliku postFactumAnalysusRank.R. 
#
#  Funkcja znajduje odpowiedzi na pytania:
#   1.) Jak mało muszę polepszyć a by osiągnało ono co najmniej k tą pozycję w rankingu, dla przynajmniej jednej funkcji (greater.than.one = TRUE - q >= 1 )
#   2.) Jak dużo mogę pogorszyć a by a ciągle posiadało możliwą pozycję w rankingu "co najmniej k" (greather.than.one = FALSE - q <= 1)
#  Przykład:  q<-PfaPossibleComprehensiveRankImprovement(perf = perf, a =1, k=2, strict.vf = TRUE, precision=0.00005, which.attributes = c(1,1), greater.than.one = TRUE)
#
PfaPossibleComprehensiveRankImprovement <- function(perf, a, k, strict.vf, strong.prefs=NULL,
                                                    weak.prefs=NULL, indif.prefs = NULL, nums.of.characteristic.points=NULL, precision=0.00005, 
                                                    which.attributes = NULL, greater.than.one = TRUE) {
  # greater than one  = TRUE q >= 1 else q < 1  
  
  if (greater.than.one) {
    q <- PfaFindQPossibleRankImprovement(perf = perf, a = a, k = k, strict.vf = strict.vf,
                                         strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                         nums.of.characteristic.points=nums.of.characteristic.points, which.attributes = which.attributes,
                                         precision=precision)
  } else {
    q <- PfaFindQPossibleRankDeterioration(perf = perf, a = a, k = k, strict.vf = strict.vf,
                                           strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                           nums.of.characteristic.points=nums.of.characteristic.points, which.attributes = which.attributes,
                                           precision=precision)
  }
  return(q)
}


#  NECESSARY RANK IMPROVEMENT
#  Funkcja znajduje odpowiedzi na pytania:
#   1.) Jak mało muszę polepszyć a by osiągnało ono co najmniej k tą pozycję w rankingu, dla Każdej funkcji (greater.than.one = TRUE - q >= 1 )
#   2.) Jak dużo mogę pogorszyć a by a ciągle posiadało konieczną pozycję w rankingu "co najmniej k" (greather.than.one = FALSE - q <= 1)
#  Przykład:  q<-PfaNecessaryComprehensiveRankImprovement(perf = performances, a =1, k=2, strict.vf = TRUE, precision=0.00005, which.attributes = c(1,1), greater.than.one = TRUE)

PfaNecessaryComprehensiveRankImprovement <- function(perf, a, k, strict.vf, strong.prefs=NULL,
                                                     weak.prefs=NULL, indif.prefs = NULL, nums.of.characteristic.points=NULL, precision=0.00005, 
                                                     which.attributes = NULL, greater.than.one = TRUE) {
  # greater than one  = TRUE q >= 1 else q < 1  
  
  if (greater.than.one) {
    q <- PfaFindQNecessaryRankImprovement(perf = perf, a = a, k = k, strict.vf = strict.vf,
                                          strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                          nums.of.characteristic.points=nums.of.characteristic.points, which.attributes = which.attributes,
                                          precision=precision)
  } else {
    q <- PfaFindQNecessaryRankDeterioration(perf = perf, a = a, k = k, strict.vf = strict.vf,
                                            strong.prefs=strong.prefs,weak.prefs=weak.prefs, indif.prefs = indif.prefs,
                                            nums.of.characteristic.points=nums.of.characteristic.points, which.attributes = which.attributes,
                                            precision=precision)
  }
  return(q)
}


initLibrary()

##############################################################################
####################### Reading Data #########################################
##############################################################################


getPerformances <- function(args) {
  
  treeAlternatives <- NULL
  treeCriteria <- NULL
  treePerformanceTable <- NULL
  
  tmpErr <- try(
    treeAlternatives<-xmlTreeParse("alternatives.xml",useInternalNodes=TRUE)
  )
  if (inherits(tmpErr, 'try-error')) {
    print("try")
    errFile <<- "Cannot read alternatives file."
    print(errFile)
  }
  if (is.null(errFile)) {
    if (checkXSD(treeAlternatives) == 0) {
      print(treeAlternatives)
      errFile <<- " Alternatives file is not XMCDA valid."
      print(errFile)
    }
  }
  
  tmpErr <- try(
    treeCriteria<-xmlTreeParse("criteria.xml",useInternalNodes=TRUE)
  )
  if (inherits(tmpErr, 'try-error')) {
    errFile <<- "Cannot read criteria file."
  }
  if (is.null(errFile)) {
    if (checkXSD(treeCriteria) == 0) {
      errFile <<- " Criteria file is not XMCDA valid."
    }
  }
  
  tmpErr <- try(
    treePerformanceTable<-xmlTreeParse("performances.xml",useInternalNodes=TRUE)
  )
  if (inherits(tmpErr, 'try-error')) {
    errFile <<- "Cannot read performance file."
  }
  if (is.null(errFile)) {
    if (checkXSD(treePerformanceTable) == 0) {
      errFile <<- " Performances file is not XMCDA valid."
    }
  }
  
  altIDsFrame <- NULL
  critIDsFrame <- NULL
  performancesFrame <- NULL
  critIDs <- NULL
  altIDs <- NULL
  performances <- NULL
  if (is.null(errFile)) {
    print("OK")
    flag <- TRUE
    critIDsFrame <- getCriteriaIDs(treeCriteria)
    if (critIDsFrame$status == "OK") {
      print("OK")
      critIDs <- critIDsFrame[[1]] 
      altIDsFrame <- getAlternativesIDs(treeAlternatives)
    } else {
      errData <<- critIDsFrame$status
      flag <- FALSE
    }
    if (flag) {
      if (altIDsFrame$status == "OK") {
        altIDs <- altIDsFrame[[1]]
        performancesFrame <- getPerformanceTables(treePerformanceTable,altIDs = altIDs, critIDs = critIDs)    
      } else {
        errData <<- altIDsFrame$status
        flag <- FALSE
      }
    }
    
    if (flag) {
      if (performancesFrame$status == "OK") {
        performances <- performancesFrame[[1]]
      } else {
        errData <<- performancesFrame$status
      }
    }
    
  } 
  
  return(performances)  
}


convertLabelsToIds <- function(mat, rowmap) {
  # mat = matrix(mat)
  #print("AAAA")
  if ((nrow(mat) < 1) || (ncol(mat) < 1)) {
    return(NULL)
  } 
  #print(mat)
  m = matrix(0, nrow = nrow(mat), ncol=ncol(mat))
  #print(m)
  for (x in 1:nrow(mat)) {
    for(y in 1:ncol(mat)) {
      m[x,y] <- rowmap[[mat[x,y]]]
    }
  }
  
  return(m)
}

getParametersData <- function() {
  tree <- NULL
  if (is.null(errFile)) {
    tmpErr<- try(
{
  tree<-xmlTreeParse("parameters.xml",useInternalNodes=TRUE)  
}
    )
  }
  if (is.null(errFile)) {
    if (!is.null(tree)) {
      if (checkXSD(tree)==0) {
        errFile<<-"File with parameters values is not XMCDA valid."
      }
    }
  }
  
  parameters <- list()
  parameters[["strict"]] <- FALSE
  parameters[["precision"]] <- 0.00005
  parameters[["is-improvement"]] <- TRUE
  parameters[["is-possible-comprehensive-modifying"]] <- TRUE
  if (is.null(errFile)) {
    parameters <- getParameters(tree)
    if (parameters[["strict"]] == 1) parameters[["strict"]] = TRUE else parameters[["strict"]] = FALSE 
    if (parameters[["is-improvement"]] == 1) parameters[["is-improvement"]] = TRUE else parameters[["is-improvement"]] = FALSE
    if (parameters[["is-possible-comprehensive-modifying"]] == 1) parameters[["is-possible-comprehensive-modifying"]] = TRUE else parameters[["is-possible-comprehensive-modifying"]] = FALSE  
  }
  return(parameters)
}


getTarget <- function(args, performances) {
  treeTarget <- NULL
  
  if (is.null(errFile)) {
    tmpErr<-try(
{
  treeTarget<-xmlTreeParse("target-rank.xml",useInternalNodes=TRUE)
}
    ) 
    if (inherits(tmpErr, 'try-error')) {
      errFile <<- "Cannot read target-rank file."
      print(errFile)
    }
  }
  if (is.null(errFile)) {
    if (!is.null(treeTarget)) {
      if (checkXSD(treeTarget)==0)
      {
        errFile<<-"Target-rank file is not XMCDA valid."        
      } 
    }
  }
  
  if (is.null(errFile)){
    errData <- NULL
    flag <- TRUE
    cmp <- NULL
    
    cmpFrame <- getAlternativesValues(treeTarget, rownames(performances))#$preferences#[[1]]
    
    print("targetRank")
    if (cmpFrame$status == "OK") {
      print("OK")
      cmp <- cmpFrame[["target-rank"]]
      i = 1
      rowmap = list()
      for (element in rownames(performances)) {
        rowmap[element] = i
        i = i + 1
      }
      print(cmp)
      target = c(rowmap[[cmp[1,1]]], cmp[1,2])
      
      return(target)
    } else {
      errData <<- cmpFrame$status
      flag <- FALSE
    }
  } 
  return(NULL)
  
}


getPreferences <- function(args, performances) {
  
  
  treeAltComparison <- NULL
  if (is.null(errFile)) {
    tmpErr<-try(
{
  treeAltComparison<-xmlTreeParse("preferences.xml",useInternalNodes=TRUE)
}
    )
  }
  if (is.null(errFile)) {
    if (!is.null(treeAltComparison)) {
      if (checkXSD(treeAltComparison)==0)
      {
        errFile<<-"Preferences file is not XMCDA valid."        
      } 
    }
  }
  print("1.3")
  if (is.null(errFile)){
    print("1.4")
    errData <- NULL
    flag <- TRUE
    cmp <- NULL
    
    cmpFrame <- getAlternativesComparisonsLabels(treeAltComparison, rownames(performances))#$preferences#[[1]]
    print("1.5")
    if (cmpFrame$status == "OK") {
      cmp <- cmpFrame$preferences
      i = 1
      rowmap = list()
      for (element in rownames(performances)) {
        rowmap[element] = i
        i = i + 1
      }
      print("1.6")
      print(rowmap)
      preferences = list()
      preferences[['indif']] <- convertLabelsToIds(matrix(cmp[cmp[,3] == 0,1:2], ncol = 2), rowmap)
      print("1.7")
      preferences[['weak']] <- convertLabelsToIds(matrix(cmp[cmp[,3] == 1,1:2], ncol = 2), rowmap)
      preferences[['strong']] <- convertLabelsToIds(matrix(cmp[cmp[,3] == 2,1:2], ncol = 2), rowmap)
      print("1.7")
      print("preferences")
      print(cmp)
      print(preferences)
      print(rowmap)
      return(preferences)
    } else {
      errData <<- cmpFrame$status
      flag <- FALSE
    }
  } 
  pref = list()
  pref[['strong']] = NULL
  pref[['weak']] = NULL
  pref[['indif']] = NULL
  return(pref)
}

getCharacteristicPoints <- function(args, performances){
  treePoints <- NULL
  if (is.null(errFile)) {
    tmpErr<-try(
{
  treePoints<-xmlTreeParse("characteristic-points.xml",useInternalNodes=TRUE)
}
    )
  }
  if (is.null(errFile)){  
    if ((!is.null(treePoints))){
      if (checkXSD(treePoints)==0)
      {
        errFile<-"Characteristic Points File is not XMCDA valid."        
      }
    }        
  }
  if (is.null(errFile)) {
    
    
    print(colnames(performances))
    values <- getCriteriaValues(treePoints, colnames(performances))#[[1]]
    if (values$status == "OK") {
      characteristicPoints <- rep(0, length(colnames(performances)))
      
      for (row_num in 1:nrow(values$characteristicPoints)) {
        row <- values$characteristicPoints[row_num,]
        characteristicPoints[row[1]] <- row[2] 
      } 
      return(characteristicPoints)
    } else {
      errData <<- values$status
    }
  } 
  return(NULL)
}

getListOfAttributes <- function(args, performances){
  treeAttributes <- NULL
  if (is.null(errFile)) {
    tmpErr<-try(
{
  treeAttributes<-xmlTreeParse("which-criteria.xml",useInternalNodes=TRUE)
}
    )
  }
  if (is.null(errFile)){  
    if ((!is.null(treeAttributes))){
      if (checkXSD(treeAttributes)==0)
      {
        errFile<-"which-criteria File is not XMCDA valid."        
      }
    }        
  }
  if (is.null(errFile)) {
    
    
    print(colnames(performances))
    values <- getCriteriaValues(treeAttributes, colnames(performances))#[[1]]
    if (values$status == "OK") {
      attributes <- rep(0, length(colnames(performances)))
      
      for (row_num in 1:nrow(values$criteriaIds)) {
        row <- values$criteriaIds[row_num,]
        if (row[2] != 0) {
          row[2] = 1
        }
        attributes[row[1]] <- row[2] 
      } 
      
      return(attributes)
    } else {
      errData <<- values$status
    }
  } 
  return(NULL)
}




#args <- commandArgs(trailingOnly = TRUE)
#args[1] <- "./data/test/alt.xml"
#args[2] <- "./data/test/cri.xml"
#args[3] <- "./data/test/per.xml"
#args[4] <- "./data/test/alt_comparison.xml"
#args[5] <- "./data/test/points.xml"
#args[6] <- "./data/test/parameters.xml"
#args[7] <- "./data/test/result-worst.xml"
#args[8] <- "./data/test/result-best.xml"
#args[9] <- "./data/test/messages.xml"


resultFile <- "evaluation-modifying-value.xml"
resultFileMessages <- "messages.xml"



#mandatory
is_proper_data = TRUE
performances = getPerformances(args)
print(errFile)
print(performances)
targetRelation <- NULL
if (is.null(performances)) {
  is_proper_data = FALSE
} else {
  targetRelation <- getTarget(args, performances)
  if (is.null(targetRelation)) {
    is_proper_data = FALSE
  }
}

#optional
pref = NULL
cp = NULL
parameters <- NULL
which.attributes = NULL
if (is_proper_data) {
  pref <- getPreferences(args, performances)
  cp <- getCharacteristicPoints(args, performances)
  parameters <- getParametersData()
  which.attributes <- getListOfAttributes(args, performances)
}

result <-NULL

if (is.null(errFile) && is_proper_data){
  tmpErr<- try (
{
  print(parameters)
  print(which.attributes)
  print(targetRelation)
  if (parameters[["is-possible-comprehensive-modifying"]]) {
    result <- PfaPossibleComprehensiveRankImprovement(perf=performances, a = targetRelation[1], k = targetRelation[2], strict.vf = parameters[["strict"]], strong.prefs = pref$strong, weak.prefs=pref$weak, indif.prefs=pref$indif, nums.of.characteristic.points=cp, precision=parameters[["precision"]], 
                                                  which.attributes = which.attributes, greater.than.one = parameters[["is-improvement"]])
    
  } else {
    result <- PfaNecessaryComprehensiveRankImprovement(perf=performances, a = targetRelation[1], k = targetRelation[2], strict.vf = parameters[["strict"]], strong.prefs = pref$strong, weak.prefs=pref$weak, indif.prefs=pref$indif, nums.of.characteristic.points=cp, precision=parameters[["precision"]], 
                                                   which.attributes = which.attributes, greater.than.one = parameters[["is-improvement"]])  
    
  }
  
}
  )  
  if (inherits(tmpErr, 'try-error')){
    errCalc<<-"Cannot find value."
  } else {
    execFlag<-TRUE
  }
}

setwd(outDirectory)

if (execFlag) {
  outTree = newXMLDoc()
  
  newXMLNode("xmcda:XMCDA",
             attrs=c("xsi:schemaLocation" = "http://www.decision-deck.org/2012/XMCDA-2.2.0 http://www.decision-deck.org/xmcda/_downloads/XMCDA-2.2.0.xsd"),
             suppressNamespaceWarning=TRUE,
             namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2012/XMCDA-2.2.0"),
             parent=outTree)
  
  
  
  ids <- rownames(performances)
  print(ids[targetRelation[1]])
  
  putAlternativeValue(outTree, result, alternativesIDs = ids[targetRelation[1]], mcdaConcept="evaluation-modifying")
  saveXML(outTree, file=resultFile)  
}

if (!is.null(errCalc)){
  print("errCalc")
  # something went wrong at the calculation step
  
  outTreeMessage = newXMLDoc()
  
  newXMLNode("xmcda:XMCDA", 
             attrs=c("xsi:schemaLocation" = "http://www.decision-deck.org/2012/XMCDA-2.2.0 http://www.decision-deck.org/xmcda/_downloads/XMCDA-2.2.0.xsd"),
             suppressNamespaceWarning=TRUE, 
             namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2012/XMCDA-2.2.0"), 
             parent=outTreeMessage)
  
  status<-putErrorMessage(outTreeMessage, errCalc, name="Error")
  
  status<-saveXML(outTreeMessage, file=resultFileMessages)
  
}

if (!is.null(errData)){
  print("errData")
  # something went wrong at the calculation step
  
  outTreeMessage = newXMLDoc()
  
  newXMLNode("xmcda:XMCDA", 
             attrs=c("xsi:schemaLocation" = "http://www.decision-deck.org/2012/XMCDA-2.2.0 http://www.decision-deck.org/xmcda/_downloads/XMCDA-2.2.0.xsd"),
             suppressNamespaceWarning=TRUE, 
             namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2012/XMCDA-2.2.0"), 
             parent=outTreeMessage)
  
  status<-putErrorMessage(outTreeMessage, errData, name="Error")
  
  status<-saveXML(outTreeMessage, file=resultFileMessages)
  
}

if (!is.null(errFile)){
  print("errFile")
  # something went wrong at the calculation step
  
  outTreeMessage = newXMLDoc()
  
  newXMLNode("xmcda:XMCDA", 
             attrs=c("xsi:schemaLocation" = "http://www.decision-deck.org/2012/XMCDA-2.2.0 http://www.decision-deck.org/xmcda/_downloads/XMCDA-2.2.0.xsd"),
             suppressNamespaceWarning=TRUE, 
             namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2012/XMCDA-2.2.0"), 
             parent=outTreeMessage)
  
  status<-putErrorMessage(outTreeMessage, errFile, name="Error")
  
  status<-saveXML(outTreeMessage, file=resultFileMessages)
  
}
#analogicznie errData i errFile
#https://github.com/sbigaret/DecisionDeck-Web-Services/blob/master/XMCDAWebServicesData/dummyRXMCDAService/dummy.R
