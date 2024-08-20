checkTestResults <- function(testresults) {
    if (1 %in% testresults) {
        nTests <- length(testresults)
        passedTests <- nTests - length(which(1 == testresults))
        cat(passedTests, "out of ", nTests, "checks passed. \n")
        stop()
    }
    else {
        cat(length(testresults), "out of", length(testresults), "checks passed. \n")
    }
    
}

myTestX <- function(x) {
    if (x == 3) {
        cat("yay the check passed! \n")
        return(0)
    }
    else {
        cat("Something about x is wrong! \n")
        return(1)
    }
}

myTestY <- function(y) {
    if (y == 3) {
        cat("yay the check passed! \n")
        return(0)
    }
    else {
        cat("Something about Y is wrong! \n")
        return(1)
    }
}

myResultX <- 3
myResultY <- 3

myTestResults <- c(myTestX(myResultX), myTestY(myResultY))
checkTestResults(myTestResults)