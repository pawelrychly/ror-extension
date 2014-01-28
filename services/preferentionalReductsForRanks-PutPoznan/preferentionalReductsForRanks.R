##########################################
# Usage:
# R --slave --vanilla --args "[inDirectory]" "[outDirectory]" < preferentionalReductsForRanks.R
# Example: 
# R --slave --vanilla --args "${PWD}/in" "${PWD}/out" < preferentionalReductsForRanks.R
##########################################

inDirectory <- commandArgs()[5] # "D:/projects/praca-magisterska/ROR-ranking/services/preferentionalReductsForRanks-PutPoznan/tests/in1" #
outDirectory <- commandArgs()[6] # "D:/projects/praca-magisterska/ROR-ranking/services/preferentionalReductsForRanks-PutPoznan/tests/out1" #

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
    
    #print("Punkty charakterystyczne:")
    #print(labels)
    #print("Wartości w punktach charakterystycznych")
    #print(v)
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
  print("OCENY WARIANTÓW NA KRYTERIACH")
  for(variant_index in seq(nrow(perf))) {
    print(sprintf("Oceny wariantu %d:", variant_index))
    for(value_index in seq(ncol(perf))) {
      v <- perf[[value_index]][[variant_index]]
      
      print(sprintf("Kryterium %d: %f", value_index , mapping.all[[toString(value_index)]][[toString(v)]]  ))
    }
    
  }
  #print(ret)
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

######################################################################################
############################### ERA ##################################################
######################################################################################



EraSolveModel <- function(performances, model) {
  number.of.variables <- ncol(model$lhs)
  number.of.binary.variables <- nrow(performances)
  position.of.first.binary.variable <- number.of.variables - number.of.binary.variables + 1
  types <- rep('C', number.of.variables)
  types[position.of.first.binary.variable : number.of.variables] <- 'B'
  obj <- L_objective(EraBuildObjective(performances, number.of.variables))
  
  #obj <- L_objective(c(EraBuildObjective(performances, number.of.variables), c(0,0,0,0)))
  
  #vars = model$lhs[15:18, ]
  #vars2 = matrix(c(-1, 0, 0, 0, 0, -1, 0, 0, 0,0, -1, 0, 0, 0, 0, -1), nrow=4)
  #vars = cbind(vars, vars2)
  #model$lhs = GetNormalizedMatrix(matrix=model$lhs, ncol(vars))
  #model$dir = c(model$dir, c("==", "==", "==", "=="))
  #model$rhs = c(model$rhs, c(0,0,0,0))
  #model$lhs = rbind(model$lhs, vars)
  ##model$dir = cbind(model$dir, matrix(data=c("==", "==", "==", "=="), nrow= 1))
  #model$rhs = cbind(model$rhs, matrix(data=c(0, 0, 0, 0), nrow= 1))
  #vars = rbind(model$lhs, vars)
  #types <-c(types, c("C","C","C","C"))
  
  roiConst <- L_constraint(L = model$lhs, dir =model$dir, rhs=model$rhs)
  lp <- OP(objective=obj, constraints=roiConst, maximum=FALSE, types=types)
  res <- ROI_solve(lp, ror:::.solver)
  return(res)
}

EraBuildObjective <- function(performances, number.of.all.variables) {
  
  objective <- rep(0.0, number.of.all.variables)
  number.of.binary.variables <- nrow(performances)
  position.of.first.binary.variable <- number.of.all.variables - number.of.binary.variables + 1
  objective[position.of.first.binary.variable:number.of.all.variables] <- 1  
  return(objective)
}

EraBuildVariantsDiffMatrix <- function(index.of.variant.1, index.of.variant.2, perf) {
  alt.vars <- ror:::buildAltVariableMatrix(perf)
  mat <- alt.vars[index.of.variant.1,]
  mat <- mat - alt.vars[index.of.variant.2,]
  mat <- matrix(data = mat, nrow = 1)
  row.names(mat) <- paste(index.of.variant.1, "-", index.of.variant.2)
  return(mat)
}

EraBuildConstraints <- function(rank.constraints, performances, is.rank.max = TRUE){
  rank.constraints.length = ncol(rank.constraints$lhs)
  number.of.binary.variables = nrow(performances) 
  variants <- seq(number.of.binary.variables)
  constraints.list <- list()
  for (i in variants){
    differences = c()
    additional.variables.list = c()
    for (j in variants) {
      if (i != j) {        
        if (is.rank.max) {
          difference <- EraBuildVariantsDiffMatrix(i,j, performances)
        } else {
          difference <- EraBuildVariantsDiffMatrix(j,i, performances)
          difference[GetNumberOfValues(performances)] = -1
        }
        
        additional.variables = vector(mode="numeric",
                                      length=number.of.binary.variables)
        additional.variables[j] = M_BIG_NUMBER
        additional.variables.list = rbind(additional.variables.list,
                                          additional.variables)
        differences <- rbind(differences, difference)
      }
    }
    differences <- GetNormalizedMatrix(differences, rank.constraints.length)
    constraints <- list(
      "lhs" = cbind(differences, additional.variables.list),
      "dir" = rep(">=", nrow(additional.variables.list)),
      "rhs" = rep( 0, nrow(additional.variables.list))
    )
    all.constraints.length <- ncol(constraints$lhs)
    rank.constraints$lhs <- GetNormalizedMatrix(rank.constraints$lhs, all.constraints.length)
    allConstraints <- ror:::combineConstraints(rank.constraints, constraints)  
    constraints.list[i] <- list(allConstraints)
    
  } 
  return(constraints.list)
}


##########################################################################################
########################## UTA GMS Extension #############################################
##########################################################################################
UTA_GMS_SolveModel <- function(perf, model) {
  
  objective <- ror:::buildObjectiveFunction(perf)
  objective <- GetNormalizedMatrix(matrix=objective, width=ncol(model$lhs))
  obj <- L_objective(objective)
  roiConst <- L_constraint(model$lhs, model$dir, model$rhs)
  lp <- OP(objective=obj, constraints=roiConst, maximum=TRUE)
  ROI_solve(lp, ror:::.solver)
}

CheckNeccesseryRelation <- function(perf, a, b, strict.vf = FALSE,
                                    strong.prefs = NULL, weak.prefs = NULL, indif.prefs = NULL, nums.of.characteristic.points = NULL) {
  
  ## check vars
  stopifnot(is.logical(strict.vf))
  if (a == b) {
    return(TRUE)
  }
  num.of.variables = GetNumberOfVariables(perf, nums.of.characteristic.points)
  alt.vars <- ror:::buildAltVariableMatrix(perf)  
  base.model <- BuildBaseLPModel(perf, strict.vf=strict.vf, 
                                 strong.prefs=strong.prefs, weak.prefs=weak.prefs,
                                 indif.prefs=indif.prefs, nums.of.characteristic.points)
  
  add.const <- c()
  add.const <- ror:::buildStrongPreferenceConstraint(b, a, alt.vars) 
  add.const$lhs <- GetNormalizedMatrix(add.const$lhs, num.of.variables);
  all.const <- ror:::combineConstraints(base.model, add.const)
  
  ret <- SolveModel(perf, all.const)
  
  return(ret$status$code != 0 || ret$objval < ror:::MINEPS)
}

#Przetestowane
UtaGMSFindNeccesseryRelations <- function(perf, strong.prefs = NULL,
                                          weak.prefs = NULL, indif.prefs = NULL, strict.vf=FALSE, nums.of.characteristic.points = NULL) {
  
  if (is.null(nums.of.characteristic.points)) {
    nums.of.characteristic.points = rep(0, ncol(perf))
  }
  
  rel <- matrix(nrow=nrow(perf), ncol=nrow(perf))
  
  base.model <- BuildBaseLPModel(perf, strict.vf, strong.prefs = strong.prefs,
                                 weak.prefs = weak.prefs, indif.prefs = indif.prefs,
                                 nums.of.characteristic.points = nums.of.characteristic.points)
  
  if (!CheckConstraintsConsistency(perf, base.model)) {
    stop("Model infeasible")
  }
  
  necessary.relations = matrix(nrow=0, ncol=2)
  for (i in 1:nrow(rel)) {
    for(j in 1:nrow(rel)) {
      
      rel[i,j] = CheckNeccesseryRelation(perf, i, j, strict.vf=strict.vf,
                                         strong.prefs=strong.prefs, weak.prefs=weak.prefs, indif.prefs=indif.prefs, nums.of.characteristic.points)
      if ((rel[i,j] == TRUE) && (i != j)) {
        necessary.relations <- rbind(necessary.relations, c(i,j))
      }
    }
  }
  #print(rel)
  return(necessary.relations)
}

##########################################################################################
########################## Preferentional Reducts ########################################
##########################################################################################



PrFindAllPreferentionalReducts <- function(perf, e.ts, e.pi.rank, name, number.of.binary.variables = 0) {
  prs <- list()
  
  if (!PrCheckConstraintsConsistency(perf, e.ts, number.of.binary.variables)) {
    return(" EMPTY SET ")
    
  }
  
  ss <- GenerateAllSubsetsOfConstraints(e.pi.rank)
  
  i <- 0
  while ( i < length(ss)) {
    i = i + 1
    
    #print(paste("i: ", i))
    constraints <- ss[[i]]
    ts <- e.ts
    ts <- ror:::combineConstraints(ts, constraints)
    if (!PrCheckConstraintsConsistency(perf, ts, number.of.binary.variables)) {
      
      prs <- append(prs, rownames(constraints$lhs))
      #prs <-append(prs, constraints)
      ##removing all supersets
      j <- i
      
      while(j < length(ss)) {
        j= j+1
        a <- ss[[j]]
        if (IsASubset(constraints, ss[[j]])) {
          
          ss[[j]] <- NULL
          j = j - 1  #Usunieto element tablicy więc trzeba cofnąć wskaźnik, żeby nie ominąć kolejnego elementu
        } 
        
      }
    }
  }
  
  return(prs)
}


IsASubset <- function(candidate.set, superset) {
  
  is.identical = TRUE
  for (i in seq(nrow(candidate.set$lhs))) {
    candidate.row <- list()
    candidate.row$lhs <- candidate.set$lhs[i, ]
    candidate.row$dir <- candidate.set$dir[i]
    candidate.row$rhs <- candidate.set$rhs[i]
    is.in.superset = FALSE
    for (j in seq(nrow(superset$lhs))) {
      superset.row <-list()
      superset.row$lhs <- superset$lhs[j, ]
      superset.row$dir <- superset$dir[j]
      superset.row$rhs <- superset$rhs[j]
      if (identical(candidate.row$lhs, superset.row$lhs) &&
            identical(candidate.row$dir, superset.row$dir) &&
            identical(candidate.row$rhs, superset.row$rhs)) {
        
        is.in.superset = TRUE;
        break
      }
    }
    if (is.in.superset == FALSE) {
      
      return(FALSE)
    }  
  }
  
  return(TRUE)
}



PrBuildConstraintsForRanksPreferentionalReducts <- function(base.constraints, id.of.variant, perf, ranks, alt.vars ) {
  
  base.constraints.width <- ncol(base.constraints$lhs)
  a.num <- nrow(perf)
  all.constraints <- list()
  number.of.binary.variables = nrow(perf) 
  differences.best = vector()
  differences.worst = vector()
  additional.variables.list.best = c()
  additional.variables.list.worst = c()
  best <- ranks[id.of.variant,2];
  worst <- ranks[id.of.variant,1]
  all.binaries.sum = matrix(ncol= (2 * number.of.binary.variables) + 2, nrow=3, data=0)
  number.of.values = GetNumberOfValues(perf)  
  
  for ( i in seq(nrow(perf))) {
    if (i != id.of.variant) {
      if (best > 1) {
        difference <- EraBuildVariantsDiffMatrix(id.of.variant, i, perf)
        difference[number.of.values] = -1
        differences.best <- rbind(differences.best, difference)
        
        additional.variables.best = vector(mode="numeric", length=number.of.binary.variables * 2)
        additional.variables.best[i] = M_BIG_NUMBER
        additional.variables.list.best = rbind(additional.variables.list.best, additional.variables.best)
        all.binaries.sum[1, i] <- 1
        
      }
      if (worst < a.num) {
        difference <- EraBuildVariantsDiffMatrix(i, id.of.variant, perf)
        difference[number.of.values] = -1
        additional.variables.worst = vector(mode="numeric", length=number.of.binary.variables * 2)
        additional.variables.worst[number.of.binary.variables + i] = M_BIG_NUMBER
        
        additional.variables.list.worst = rbind(additional.variables.list.worst,
                                                additional.variables.worst)
        differences.worst <- rbind(differences.worst, difference)
        all.binaries.sum[2, number.of.binary.variables + i] <- 1
        
      }
    } 
  }
  
  binaries.constraints <- list(
    lhs = all.binaries.sum,
    dir = c("<=", "<=", "=="),
    rhs = c(0, 0, 0)
  )
  if (best > 1) {
    binaries.constraints$lhs[1, 2 * number.of.binary.variables + 1] <- -a.num
    binaries.constraints$lhs[3, 2 * number.of.binary.variables + 1] <- 1
    binaries.constraints$rhs[1]= best - 2
  }
  if (worst < a.num) {
    binaries.constraints$lhs[2, 2 * number.of.binary.variables + 2] <- -a.num
    binaries.constraints$lhs[3, 2 * number.of.binary.variables + 2] <- 1
    binaries.constraints$rhs[2]= a.num - worst - 1
  }
  binaries.constraints$rhs[3]= 1
  
  if (!is.null(nrow(differences.best))) {
    #Normalizacja best
    differences.best <- GetNormalizedMatrix(
      matrix = differences.best,
      width = base.constraints.width
    )
    differences.best <- cbind(differences.best, additional.variables.list.best)
    
  }
  if (!is.null(nrow(differences.worst))) {
    #Normalizacja worst
    differences.worst <- GetNormalizedMatrix(
      matrix = differences.worst,
      width = base.constraints.width
    )
    
    differences.worst <- cbind(differences.worst, additional.variables.list.worst) 
  }
  
  differences.constraint.best <- list()
  differences.constraint.worst <- list()
  if (!is.null(nrow(differences.best))) {
    differences.constraint.best <- list(
      lhs = differences.best,  
      dir = rep(">=", nrow(differences.best)),
      rhs = rep( 0, nrow(differences.best))
    )
  }
  if (!is.null(nrow(differences.worst))) {
    differences.constraint.worst <- list(
      lhs = differences.worst,  
      dir = rep(">=", nrow(differences.worst)),
      rhs = rep( 0, nrow(differences.worst))
    )
  }
  
  differences.constraints <- list()
  if ((!is.null(nrow(differences.best))) && (!is.null(nrow(differences.worst)))) {
    differences.constraints <- ror:::combineConstraints(differences.constraint.worst, differences.constraint.best )
  } else {
    if (!is.null(nrow(differences.best))) {
      differences.constraints <- differences.constraint.best 
    }
    if (!is.null(nrow(differences.worst))) {
      differences.constraints <- differences.constraint.worst     
    }    
  }
  differences.constraints$lhs <- GetNormalizedMatrix(
    matrix = differences.constraints$lhs,
    width = ncol(differences.constraints$lhs) + 2
  )
  
  binaries.constraints$lhs <- GetNormalizedMatrix(
    matrix = binaries.constraints$lhs,
    width = ncol(differences.constraints$lhs),
    right = FALSE
  )
  #binaries.constraints <- 
  base.constraints$lhs <- GetNormalizedMatrix(
    matrix = base.constraints$lhs,
    width = ncol(differences.constraints$lhs)
  )
  all.constraints <- ror:::combineConstraints(base.constraints, differences.constraints, binaries.constraints)
  
  return(all.constraints)
}


PrSolveModel <- function(performances, model, number.of.binary.variables) {
  number.of.variables <- ncol(model$lhs)
  position.of.first.binary.variable <- number.of.variables - number.of.binary.variables + 1
  types <- rep('C', number.of.variables)
  types[position.of.first.binary.variable : number.of.variables] <- 'B'
  obj <- L_objective(EraBuildObjective(performances, number.of.variables))
  roiConst <- L_constraint(L = model$lhs, dir =model$dir, rhs=model$rhs)
  
  lp <- OP(objective=obj, constraints=roiConst, maximum=TRUE, types=types)
  #print("--l_objective----")
  #print(obj)
  res <- ROI_solve(lp, ror:::.solver)
  #print("--objval----")
  #print(res$objval)
  #print("---solution---")
  #print(res$solution)
  return(res)
}


PrCheckConstraintsConsistency <- function(perf, model, number.of.binary.variables = 0) {
  #check constraints consistency  
  #  model: structure of constraints with the following elements:
  #    model$lhs: matrix - left side of constraints
  #    model$dir: list of operators
  #    model$rhs: matrix - right side of constraints
  
  #baseModel <- BuildBaseLPModel(perf, strict.vf, strong.prefs = NULL,
  #                           weak.prefs = NULL, indif.prefs = NULL, is.piecewise.linear = FALSE,
  #                           nums.of.characteristicPoints=NULL)
  
  if (number.of.binary.variables == 0) {
    ret <- UTA_GMS_SolveModel(perf, model)
  } else {
    ret <- PrSolveModel(perf, model, number.of.binary.variables)
  }
  return(ret$status$code == 0 && ret$objval >= ror:::MINEPS)
}

#####################################################################################
############################# Preferentional reducts for ranks ######################
#####################################################################################


# Funkcja znajduje redukty preferencyjne dla każdego wariantu takie, że jego najlepsza pozycja w rankingu nie będzie
# lepsza od P^*, oraz jego najgorsza pozycja w rankingu nie będzie gorsza od P_*.
# ranks - wynik metody ExtremeRankingAnalysis, macierz n x 2, w której pierwsza kolumna oznacza P_*, 2 kolumna to P^*
# Przykład uruchomienia metody:
#   pr <- PrFindAllRanksPreferentionalReducts(perf=performances, ranks = reducts, strong.prefs = str, strict.vf=FALSE)
#
PrFindAllRanksPreferentionalReducts <- function(perf, ranks, strong.prefs = NULL, weak.prefs = NULL,
                                                indif.prefs = NULL, strict.vf=FALSE, nums.of.characteristic.points = NULL) {
  
  #Function finds all ranks preferentional reducts
  # ranks:  an n x 2 matrix, where each row describes worst and best rank for variant 
  if (is.null(nums.of.characteristic.points)) {
    nums.of.characteristic.points = rep(2, ncol(perf))
  }
  
  num.of.variables = GetNumberOfVariables(perf, nums.of.characteristic.points)  
  base.constraints <- BuildBaseConstraints(perf, num.of.variables, strict.vf, nums.of.characteristic.points)
  pairwise.comparison.constraints <- BuildPairwiseComparisonConstraints(perf, strong.prefs, weak.prefs, indif.prefs)
  pairwise.comparison.constraints$lhs <- GetNormalizedMatrix(pairwise.comparison.constraints$lhs, num.of.variables);
  
  reducts <- list();
  alt.vars <- ror:::buildAltVariableMatrix(perf)  
  for(i in seq(nrow(perf))) {
    
    constraints <- PrBuildConstraintsForRanksPreferentionalReducts( base.constraints, i, perf, ranks, alt.vars )
    pairwise.comparison.constraints$lhs <- GetNormalizedMatrix(pairwise.comparison.constraints$lhs, ncol(constraints$lhs));
    
    #print(pairwise.comparison.constraints)
    num.of.bin <- 2 * nrow(perf) + 2
    
    result <-PrFindAllPreferentionalReducts(perf, 
                                            e.ts = constraints , e.pi.rank = pairwise.comparison.constraints, name = paste(i), num.of.bin)
    
    if (length(result) == 0) {
      result <- "EMPTY SET"
    }
    #reducts<- append(reducts, list(name = paste(i,"-", ranks[i,1], ":", ranks[i,2]), reduct = result))
    reducts[[paste(i, ":[", ranks[i,1], ":", ranks[i,2], "]")]] <- result
  }
  
  #for (alternative in names(reducts)) {
  #  print(paste(alternative, ":"))
  #  elements = list()
  #  
  #  for (element in reducts[alternative] ){
  #    element <- paste("[", element , "]")
      #print(element)
  #    elements <- paste(elements, " ", element)
  #  }
  #  print(elements)
  #}
  
  return(reducts)
  #PrFindAllPreferentionalReducts(e.ts = , e.pi.rank = peirwise.comparison.constraints)
}



initLibrary()

##############################################################################
####################### RXMCDA Extension ######################################
##############################################################################

putAlternativesComparisonsReductsForRanks <- function (tree, alternativesComparisons, mcdaConcept = NULL) 
{
  out <- list()
  err1 <- NULL
  err2 <- NULL
  racine <- NULL
  tmpErr <- try({
    racine <- xmlRoot(tree)
  })
  if (inherits(tmpErr, "try-error")) {
    err1 <- "No <xmcda:XMCDA> found."
  }
  if (length(racine) != 0) {
    if (!is.null(mcdaConcept)) {
      altVals <- newXMLNode("alternativesComparisons", 
                            attrs = c(mcdaConcept = mcdaConcept), parent = racine, 
                            namespace = c())
    }
    else {
      altVals <- newXMLNode("alternativesComparisons", 
                            parent = racine, namespace = c())
    }
    pairs <- newXMLNode("pairs", parent = altVals, namespace = c())
    for (i in 1:dim(alternativesComparisons)[1]) {
      tmpErr <- try({
        pair <- newXMLNode("pair", parent = pairs, namespace = c())
        initial <- newXMLNode("initial", parent = pair, 
                              namespace = c())
        newXMLNode("alternativeID", alternativesComparisons[i,1], parent = initial, namespace = c())
        terminal <- newXMLNode("terminal", parent = pair, 
                               namespace = c())
        newXMLNode("alternativeID", alternativesComparisons[i,2], parent = terminal, namespace = c())
        
        name_value <- NULL
        if (dim(alternativesComparisons)[2] > 3) {
          name_value <- alternativesComparisons[i,4]
        }
        val <- newXMLNode("value", attrs=c(name = name_value ), parent = pair, namespace = c())
        if (is.na(alternativesComparisons[i, 3])) {
          newXMLNode("NA", alternativesComparisons[i,3], parent = val, namespace = c())
        }
        else newXMLNode("label", alternativesComparisons[i, 3], parent = val, namespace = c())
        
      })
      if (inherits(tmpErr, "try-error")) {
        err2 <- "Impossible to put (a) value(s) in a <alternativesComparisons>."
      }
    }
  }
  if (!is.null(err1) | (!is.null(err2))) {
    out <- c(out, list(status = c(err1, err2)))
  }
  else {
    out <- c(out, list(status = "OK"))
  }
  return(out)
}


##############################################################################
####################### Reading Data #########################################
##############################################################################
getRanks <- function(performances) {
  treeWorstRanks <- NULL
  treeBestRanks <- NULL
  worstRanks <- NULL
  bestRanks <- NULL
  ranks <- NULL
  tmpErr <- try(
    treeWorstRanks<-xmlTreeParse("worst-ranking.xml",useInternalNodes=TRUE)
  )
  if (inherits(tmpErr, 'try-error')) {
    errFile <<- "Cannot read worst-ranking file."
    print(errFile)
  }
  if (is.null(errFile)) {
    if (checkXSD(treeWorstRanks) == 0) {
      errFile <<- " Worst ranking file is not XMCDA valid."
      print(errFile)
    }
  }
  tmpErr <- try(
    treeBestRanks<-xmlTreeParse("best-ranking.xml",useInternalNodes=TRUE)
  )
  if (inherits(tmpErr, 'try-error')) {
    errFile <<- "Cannot read best-ranking file."
    print(errFile)
  }
  if (is.null(errFile)) {
    if (checkXSD(treeWorstRanks) == 0) {
      errFile <<- " Best ranking file is not XMCDA valid."
      print(errFile)
    }
  }
  
  
  
  if (is.null(errFile)) {
    
    alternativesIds = rownames(performances)
    flag <- TRUE
    worstRanksFrame <- getAlternativesValues(treeWorstRanks, alternativesIDs=alternativesIds)
    bestRanksFrame <- NULL
    if (worstRanksFrame$status == "OK") {
      worstRanks <- worstRanksFrame[[1]] 
      bestRanksFrame <- getAlternativesValues(treeBestRanks, alternativesIDs=alternativesIds)
    } else {
      errData <<- worstRanksFrame$status
      flag <- FALSE
    }
    if (flag) {
      if (bestRanksFrame$status == "OK") {
        bestRanks <- bestRanksFrame[[1]]
        ranks <- cbind(worstRanks[,2], bestRanks[,2])
        
        
      } else {
        errData <<- bestRanksFrame$status
      }
    }
  }
  
  return(ranks)
}


getPerformances <- function(args) {
  
  treeAlternatives <- NULL
  treeCriteria <- NULL
  treePerformanceTable <- NULL
  
  tmpErr <- try(
    treeAlternatives<-xmlTreeParse("alternatives.xml",useInternalNodes=TRUE)
  )
  if (inherits(tmpErr, 'try-error')) {
    #print("try")
    errFile <<- "Cannot read alternatives file."
    print(errFile)
  }
  if (is.null(errFile)) {
    if (checkXSD(treeAlternatives) == 0) {
      #print(treeAlternatives)
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
    
    flag <- TRUE
    critIDsFrame <- getCriteriaIDs(treeCriteria)
    if (critIDsFrame$status == "OK") {
      
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

getStrictParameter <- function(args) {
  tree <- NULL
  if (is.null(errFile)) {
    tmpErr<- try(
{
  tree<-xmlTreeParse("strict.xml",useInternalNodes=TRUE)  
}
    )
  }
  if (is.null(errFile)) {
    if (!is.null(tree)) {
      if (checkXSD(tree)==0) {
        errFile<<-"File with parameter value is not XMCDA valid."
      }
    }
  }
  if (is.null(errFile)) {
    parameters <- getParameters(tree, "strict")
    strict <- FALSE
    if (parameters$strict) return(TRUE) else return(FALSE)    
  }
  return(FALSE)
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
  
  if (is.null(errFile)){
    
    errData <- NULL
    flag <- TRUE
    cmp <- NULL
    
    cmpFrame <- getAlternativesComparisonsLabels(treeAltComparison, rownames(performances))#$preferences#[[1]]
    
    if (cmpFrame$status == "OK") {
      cmp <- cmpFrame$preferences
      i = 1
      rowmap = list()
      for (element in rownames(performances)) {
        rowmap[element] = i
        i = i + 1
      }
      preferences = list()
      preferences[['indif']] <- convertLabelsToIds(matrix(cmp[cmp[,3] == 0,1:2], ncol = 2), rowmap)
      preferences[['weak']] <- convertLabelsToIds(matrix(cmp[cmp[,3] == 1,1:2], ncol = 2), rowmap)
      preferences[['strong']] <- convertLabelsToIds(matrix(cmp[cmp[,3] == 2,1:2], ncol = 2), rowmap)
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
    
    
    #print(colnames(performances))
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


reductsFile <- "reducts-by-alternatives.xml"
resultFileMessages <- "messages.xml"



#mandatory
is_proper_data = TRUE
ranks = NULL
performances = getPerformances(args)
#print(errFile)
#print(performances)
if (is.null(performances)) {
  is_proper_data = FALSE
} else {
  ranks <- getRanks(performances)
  if (is.null(ranks)) {
    is_proper_data = FALSE
  }
}



#optional
pref = NULL
cp = NULL
strict = FALSE
if (is_proper_data) {
  pref <- getPreferences(args, performances)
  cp <- getCharacteristicPoints(args, performances)
  strict <- getStrictParameter(args)
}
results <-NULL
#print(performances)

if (is.null(errFile) && is_proper_data){
  tmpErr<- try (
{    
  #print(performances)
  #print(pref)
  results <- PrFindAllRanksPreferentionalReducts(perf=performances, ranks=ranks, strong.prefs=pref$strong, weak.prefs=pref$weak, indif.prefs=pref$indif, strict.vf=strict, nums.of.characteristic.points=cp)
  print("results:")
  #print(results)
}
  )  
  if (inherits(tmpErr, 'try-error')){
    errCalc<<-"Cannot find reducts"
  } else {
    execFlag<-TRUE
  }
}

setwd(outDirectory)

if (execFlag) {
  outTreeReducts = newXMLDoc()
  
  newXMLNode("xmcda:XMCDA",
             attrs=c("xsi:schemaLocation" = "http://www.decision-deck.org/2012/XMCDA-2.2.0 http://www.decision-deck.org/xmcda/_downloads/XMCDA-2.2.0.xsd"),
             suppressNamespaceWarning=TRUE,
             namespace = c("xsi" = "http://www.w3.org/2001/XMLSchema-instance", "xmcda" = "http://www.decision-deck.org/2012/XMCDA-2.2.0"),
             parent=outTreeReducts)
  
  alternativeIds <- rownames(performances)
  comparisons_matrix = matrix(nrow = 0, ncol = 4)
  for (key in names(results)) {
    alternativeNumber <- strsplit(key, " :[", TRUE)[[1]]
    alternativeId <- alternativeIds[as.numeric(alternativeNumber[1])]
    comparisons <- results[[key]]
    for (comparison in comparisons) {
      if (comparison == "EMPTY SET") {
        next
      }
      comparison_elements <- strsplit(comparison, " ", TRUE)[[1]]
      operator = "indifference"
      if (comparison_elements[2] == ">=") {
        operator = "greater-or-equal"  
      }
      if (comparison_elements[2] == ">") {
        operator = "greater-than"  
      }
      comparisons_matrix <- rbind(comparisons_matrix, c(alternativeIds[as.numeric(comparison_elements[1])], alternativeIds[as.numeric(comparison_elements[3])], alternativeId,  operator))
      #print(comparison_elements)
    }    
  }
  print(comparisons_matrix) 
  putAlternativesComparisonsReductsForRanks(outTreeReducts, alternativesComparisons=comparisons_matrix)
  #print(relations)
  #print(comparisons_matrix)
  #saveXML(outTreeRelations, file=necessaryRelationsFile)
  saveXML(outTreeReducts, file=reductsFile)
  
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





