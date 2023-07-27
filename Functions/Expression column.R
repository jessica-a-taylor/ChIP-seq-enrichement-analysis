# Function to add a column for expression level.
expressionColumn <- function(dataToUse, level) {
  if (nrow(dataToUse) >= 1) {
    dataToUse <- cbind(dataToUse, data.frame(Expression = rep(str_match(level, "^.*_(.*)$")[,-1], times = nrow(dataToUse))))
  }
  else dataToUse <- dataToUse
  return(dataToUse)
}