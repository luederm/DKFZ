# We will try to estimate petal width (column 4), from sepal length, sepal width, and petal length

# First we must standardize our data
irisData = scale(iris[,1:4], center = T, scale = T)
N = nrow(irisData)

# Check if we standardized properly
(1/N) * sum(irisData[,1]) # Should be 0 if properly centered
(1/N) * sum(irisData[,4]) # Should be 0 if properly centered (Response Variable)
(1/N) * sum(irisData[,1] ^ 2) # Should equal 1

lambda = 0
beta = c(0, 0, 0)

sumOfSquares = 0
for (i in 1:N)
{
  sumOfSquares = sumOfSquares + ( irisData[i,4] - sum(irisData[i,1:3] * beta) )^2
}
solution = (1/(2*N)) * sumOfSquares + lambda * sum(abs(beta))


# ***** Simple case: Single predictor *****
predictor = iris[,3]
response = iris[,4]

# Before standardization
plot(predictor, response)

# Standardize
predictor = scale(predictor, center = T, scale = T)
response = scale(response, center = T, scale = T)

# Check if we standardized properly
N = nrow(predictor)
(1/N) * sum(predictor) # Should be 0 if properly centered
(1/N) * sum(response) # Should be 0 if properly centered (Response Variable)
(1/N) * sum(predictor ^ 2) # Should equal 1

# After standardization
plot(predictor, response)

# Solution
lambda = 0.2
beta = ((1/N)*(as.vector(predictor) %*% as.vector(response)))[1,1]
minimizationFunc = (1/(2*N)) * sum( (response - (predictor*beta)) ^ 2 ) + lambda * abs(beta)

# Apply soft-thresholding
if (beta > lambda)
{
  beta = beta - lambda
  
} else if (beta < (0-lambda))
{
  beta = beta + lambda
  
} else
{
  beta = 0
}

# Plot what we got in blue and the general least squares solution in red
plot(predictor, response)
abline(0, beta, col = "blue")
abline(lm(predictor~response), col = "red")

