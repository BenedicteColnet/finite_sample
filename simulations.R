library(ggplot2)

# simulation set up according to Wager & Nie
generate_simulation <- function(n = 1000, p = 12, setup = "D"){
  
  # set-ups
  if (setup == "D"){
    X = matrix(rnorm(n*p), n, p)
    b = (pmax(X[,2] + X[,3] + X[,4], 0) + pmax(X[,5] + X[,6], 0)) / 2
    e = 1/(1 + exp(-X[,1]) + exp(-X[,2] + exp(-X[,3])))
    tau = pmax(X[,2] + X[,3] + X[,4], 0) - pmax(X[,5] + X[,6], 0)
    
  } else if (setup == "A") {
    X = matrix(runif(n*p, min=0, max=1), n, p)
    b = sin(pi * X[,2] * X[,3]) + 2 * (X[,4] - 0.5)^2 + X[,5] + 0.5 * X[,6]
    eta = 0.1
    e = pmax(eta, pmin(sin(pi * X[,1] * X[,2] * X[,3]), 1-eta))
    tau = (X[,3] + X[,4]) / 2
  } else {
    print("Wrong setup.")
    break
  }
  
  
  # complete potential outcomes, treatment, and observed outcome
  simulation <- data.frame(X = X, b = b, tau = tau, e = e)
  simulation$Y_0 <- simulation$b - 0.5*simulation$tau + rnorm(n, mean = 0, sd = 0.1)
  simulation$Y_1 <- simulation$b + 0.5*simulation$tau + rnorm(n, mean = 0, sd = 0.1)
  simulation$A <- rbinom(n, size = 1, prob = simulation$e)
  simulation$Y <- ifelse(simulation$A == 1, simulation$Y_1, simulation$Y_0)
  
  return(simulation)
}

# toy simulation
a_simulation <- generate_simulation()

ggplot(a_simulation, aes(x = X.1, y = Y)) +
  geom_point()
ggsave("./data/test.pdf")

write.csv(x=a_simulation, file="./data/tmp.csv")