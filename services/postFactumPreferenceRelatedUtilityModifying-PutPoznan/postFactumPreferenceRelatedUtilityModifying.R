##########################################
# Usage:
# R --slave --vanilla --args "[inDirectory]" "[outDirectory]" < postFactumPreferenceRelatedUtilityModyfying.R
# Example: 
# R --slave --vanilla --args "${PWD}/in" "${PWD}/out" < postFactumPreferenceRelatedUtilityModyfying.R
##########################################

inDirectory <- commandArgs()[5]  #"D:/projects/praca-magisterska/ROR-ranking/services/postFactumPreferenceRelatedUtilityModifying-PutPoznan/tests/in3" #commandArgs()[5]  #  
outDirectory <- commandArgs()[6] #"D:/projects/praca-magisterska/ROR-ranking/services/postFactumPreferenceRelatedUtilityModifying-PutPoznan/tests/out3" #commandArgs()[6]  #

##########################################
# Set the working directory as the "in" directory
##########################################

setwd(inDirectory)

errFile<-NULL
errData<-NULL
errCalc<-NULL
execFlag<-FALSE
M_BIG_NUMBER <- 1000
initLibrary <- function() {
  options( java.parameters = "-Xmx2g" )
  library(ror)
  library(RXMCDA)
  library(rJava)
  assignInNamespace("MINEPS", ns="ror", value=0.0001)
  
  M_BIG_NUMBER <- 1000 
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
  #par(mfrow = c(3,2))
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
  #print(perf)
  #print(seq(nrow(perf)))
  #print(seq(ncol(perf)))
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
  #print(ret)
  #mozliwa jest sytuacja, że częsc atrybutow jest odcinkami liniowa
  i_attributes <- 1 #tutaj zaczynają sie atrybuty
  i_points <- sum(as.numeric(lapply(levels, length))) + 2   #tutaj zaczynaja sie punkty charakt.
  orginal.or.characteristic = rep(0, ncol(perf))
  
  
  
  if (!is.null(nums.of.characteristic.points)) {
    #i <- sum(as.numeric(lapply(levels, length))) + 2
    levels_characteristic = GetCharacteristicPoints(perf=perf, nums.of.characteristic.points = nums.of.characteristic.points)
    index <- 1
    for ( num in nums.of.characteristic.points) {
      if (num > 1) {
        orginal.or.characteristic[index] = 1
        levels[[index]] <- levels_characteristic[[index]]
      }
      index = index + 1
    }
  } 
  
  #print(solution)
  
  index <- 1
  for (values in levels) {
    
    if (orginal.or.characteristic[index] == 1) {
      labels <- c()
      v <- c()
      for (value in values) {
        labels <- c(labels, value)
        v <- c(v, solution[i_points])
        i_points = i_points+1
      }
      i_attributes = i_attributes + ncol(perf)
      
      print("Punkty charakterystyczne:")
      print(labels)
      print("Wartości w punktach charakterystycznych")
      print(v)
      #plot(labels, v, type="b", col="red")
    } else {
      labels <- c()
      v <- c()
      for (value in values) {
        labels <- c(labels, value)
        v <- c(v, solution[i_attributes])
        i_attributes = i_attributes+1
      }
      
      print("wartosci atrybutow:")
      print(labels)
      print("Wartości dla atrybutow")
      print(v)
      #plot(labels, v, type="b", col="red")
      i_points = i_points + nums.of.characteristic.points[[index]]
    }
    index <- index + 1
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

PfaUmissingSolveModel <- function(performances, model, maximum=FALSE) {
  number.of.variables <- ncol(model$lhs)
  types <- rep('C', number.of.variables)
  
  obj <- L_objective(PfaUmissingBuildObjective(performances, number.of.variables))
  roiConst <- L_constraint(L = model$lhs, dir = model$dir, rhs=model$rhs)
  lp <- OP(objective=obj, constraints=roiConst, maximum=maximum, types=types)
  
  res <- ROI_solve(lp, ror:::.solver)
  return(res)
}

PfaUmissingBuildObjective <- function(performances, number.of.all.variables) {
  
  objective <- rep(0.0, number.of.all.variables)
  objective[number.of.all.variables] <- 1  
  return(objective)
}




PfaBuildConstraintsForUmissing <- function(perf, a, b, is.possible=TRUE, strict.vf = FALSE,
                                           strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL,
                                           nums.of.characteristic.points=NULL, improvement = TRUE) {
  #a: index of variant a
  #b: index of variant b
  #improvement: bool  TRUE/FALSE(deterioration)
  alt.vars <- ror:::buildAltVariableMatrix(perf)
  base.model <- BuildBaseLPModel(perf, strict.vf, strong.prefs = strong.prefs,
                                 weak.prefs = weak.prefs, indif.prefs = indif.prefs,
                                 nums.of.characteristic.points=nums.of.characteristic.points)
  
  if (!CheckConstraintsConsistency(perf, base.model)) {
    stop("Model infeasible")
  }
  
  a.pref.b = vector()
  if (is.possible) {
    a.pref.b <- ror:::buildWeakPreferenceConstraint(a, b, alt.vars) 
  } else {  #necessary
    a.pref.b <- ror:::buildWeakPreferenceConstraint(b, a, alt.vars)
  }
  a.pref.b$lhs <- matrix(data = a.pref.b$lhs, nrow=1)
  
  if (ncol(base.model$lhs) >= ncol(a.pref.b$lhs)) {
    a.pref.b$lhs <- GetNormalizedMatrix(a.pref.b$lhs, ncol(base.model$lhs))  
  } else {
    base.model$lhs <- GetNormalizedMatrix(base.model$lhs, ncol(a.pref.b$lhs))
  }
  
  #TODO
  if (((is.possible) && (improvement)) || ((!is.possible) && (!improvement))) {
    #if (improvement) {
    a.pref.b$lhs <- cbind(a.pref.b$lhs, 1)  
  } else {
    a.pref.b$lhs <- cbind(a.pref.b$lhs, -1)
  }
  base.model$lhs <- GetNormalizedMatrix(base.model$lhs, ncol(a.pref.b$lhs))
  greater.than.zero <- list(lhs = rep(0.0, ncol(base.model$lhs)), dir=">=", rhs=0)
  greater.than.zero$lhs[ncol(base.model$lhs)] <- 1
  
  all.constraints <-ror:::combineConstraints(base.model, a.pref.b, greater.than.zero)
  
  return(all.constraints)
}



#############################################################################################
########################### Preference Related Utility Modifying ###############################
#############################################################################################



# Funkcja znajduje odpowiedzi na pytania:
#   1.) Jak mało muszę się polepszyć by coś możliwie zachodziło ( improvement == TRUE and is.possibly.preffered == TRUE ) 
#   2.) Jak dużo mogę się pogorszyć, by coś wciąż możliwie zachodziło ( improvement == FALSE and is.possibly.preffered == TRUE )
#   3.) Jak dużo muszę się poprawiać by możliwie zachodziło przeciwne do tego co sprawdzamy ( improvement == TRUE and is.possibly.preffered == FALSE )
#   4.) Jak mało muszę się pogorszyć by możliwie zachodziło przeciwne do tego co sprawdzamy ( improvement == FALSE and is.possibly.preffered == FALSE )
#   
#   Przykład wywołania: Umiss <-PfaUmissingNecessaryOrPossibleComprehensiveImprovement(perf=performances, strong.prefs=str, a = 4, b = 3, strict.vf = TRUE, is.possibly.preffered = TRUE, nums.of.characteristic.points=c(4,4), precision=0.000005, improvement=FALSE)
#


PfaUmissingNecessaryOrPossibleComprehensiveImprovement <- function(perf, a, b, strict.vf, is.possibly.preffered,
                                                                   strong.prefs=NULL,weak.prefs=NULL, indif.prefs = NULL,
                                                                   nums.of.characteristic.points=NULL, improvement=TRUE) {
  
  #print("PRECISION")
  
  constraints <- PfaBuildConstraintsForUmissing(perf = perf, a = a, b = b, is.possible =  is.possibly.preffered, strict.vf = strict.vf,
                                                strong.prefs = strong.prefs, weak.prefs = weak.prefs, indif.prefs = indif.prefs,
                                                nums.of.characteristic.points = nums.of.characteristic.points, improvement=improvement)
  
  maximum = TRUE
  if (((is.possibly.preffered) && (improvement)) || ((!is.possibly.preffered) && (!improvement))) { 
    maximum = FALSE
  } else {
    maximum = TRUE
  }
  
  ret <- PfaUmissingSolveModel(perf, constraints,  maximum = maximum)
  ShowSolution(perf=perf, nums.of.characteristic.points=nums.of.characteristic.points, ret=ret)
  return(ret$objval)
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


getTargetRelation <- function(args, performances) {
  treeTarget <- NULL
  
  
  
  if (is.null(errFile)) {
    tmpErr<-try(
{
  treeTarget<-xmlTreeParse("target-relation.xml",useInternalNodes=TRUE)
}
    ) 
    if (inherits(tmpErr, 'try-error')) {
      errFile <<- "Cannot read target-relation file."
      print(errFile)
    }
  }
  if (is.null(errFile)) {
    if (!is.null(treeTarget)) {
      if (checkXSD(treeTarget)==0)
      {
        errFile<<-"Target-relation file is not XMCDA valid."        
      } 
    }
  }
  
  if (is.null(errFile)){
    errData <- NULL
    flag <- TRUE
    cmp <- NULL
    
    cmpFrame <- getAlternativesComparisonsLabels(treeTarget, rownames(performances))#$preferences#[[1]]
    
    
    if (cmpFrame$status == "OK") {
      cmp <- cmpFrame[["target-relation"]]
      i = 1
      rowmap = list()
      for (element in rownames(performances)) {
        rowmap[element] = i
        i = i + 1
      }
      print(cmp)
      target = c(rowmap[[cmp[1,1]]], rowmap[[cmp[1,2]]])
      
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


resultFile <- "utility-modifying-value.xml"
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
  targetRelation <- getTargetRelation(args, performances)
  if (is.null(targetRelation)) {
    is_proper_data = FALSE
  }
}

#optional
pref = NULL
cp = NULL
parameters <- NULL
#which.attributes = NULL
if (is_proper_data) {
  pref <- getPreferences(args, performances)
  cp <- getCharacteristicPoints(args, performances)
  parameters <- getParametersData()
  #which.attributes <- getListOfAttributes(args, performances)
}

result <-NULL

if (is.null(errFile) && is_proper_data){
  tmpErr<- try (
    {
      print(pref)
      print(parameters)
      result <- PfaUmissingNecessaryOrPossibleComprehensiveImprovement(perf=performances, a = targetRelation[1], b =  targetRelation[2], strict.vf = parameters[["strict"]], is.possibly.preffered = parameters[["is-possible-comprehensive-modifying"]],
                                                                       strong.prefs = pref$strong, weak.prefs=pref$weak, indif.prefs=pref$indif, nums.of.characteristic.points=cp, improvement=parameters[["is-improvement"]]) 
      print(result)
      
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
  #print(ids[targetRelation[1]])
  
  putAlternativeValue(outTree, result, alternativesIDs = ids[targetRelation[1]], mcdaConcept="utility-modifying")
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
