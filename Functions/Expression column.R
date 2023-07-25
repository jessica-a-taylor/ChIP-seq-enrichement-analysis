# Function to add a column for expression level.
expressionColumn <- function(dataToUse, level) {
  if (nrow(dataToUse) >= 1) {
    dataToUse <- cbind(dataToUse, data.frame(Expression = rep(level, times = nrow(dataToUse))))
  }
  else dataToUse <- dataToUse
  return(dataToUse)
}