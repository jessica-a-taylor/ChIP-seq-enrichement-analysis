# Function to add a column for expression level.
expressionColumn <- function(dataToUse, set) {
  if (nrow(dataToUse) >= 1) {
    dataToUse <- cbind(dataToUse, data.frame(GeneSet = rep(set, times = nrow(dataToUse))))
  }
  else dataToUse <- dataToUse
  return(dataToUse)
}