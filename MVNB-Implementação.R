k <- 11  # Inclui o intercepto como parte dos betas
n <- nrow(Y)  # Número de observações
d <- ncol(Y)  # Número de variáveis de resposta

empirical_means <- colMeans(Y)
empirical_variances <- apply(Y, 2, var)

# Valores iniciais para os parâmetros, agora com intercepto
init_beta <- matrix(NA, nrow = d, ncol = k + 1)  # Inclui coluna para o intercepto

for (t in 1:d) {
  nb_model <- glm.nb(Y[, t] ~ X)  # Ajuste do modelo NB com intercepto
  init_beta[t, ] <- coef(nb_model)  # Armazena os coeficientes, incluindo o intercepto
}

init_lambda <- rep(0, d * (d - 1) / 2)  # Inicializando lambda

init_m <- numeric(d)
for (t in 1:d) {
  theta_t <- 1 - empirical_means[t] / empirical_variances[t]  # Estimativa de theta_t
  init_m[t] <- empirical_means[t] * (1 - theta_t) / theta_t   # Estimativa de m_t
}

init_params <- c(as.vector(init_beta), init_lambda, init_m)

iteration_counter <<- 0  # Variável global para contagem de iterações

# Definindo a função de verossimilhança
log_likelihood <- function(params, y, X, d, k, c) {
  iteration_counter <<- iteration_counter + 1
  
  # Extraindo beta, lambda e m dos parâmetros
  beta <- matrix(params[1:(d * (k + 1))], nrow = d, ncol = k + 1)
  lambda_values <- params[(d * (k + 1) + 1):(d * (k + 1) + d * (d - 1) / 2)]
  m <- params[(d * (k + 1) + d * (d - 1) / 2 + 1):length(params)]
  
  # Construindo a matriz simétrica de lambda
  lambda <- matrix(0, d, d)
  idx <- 1
  for (l in 1:(d - 1)) {
    for (nu in (l + 1):d) {
      lambda[l, nu] <- lambda_values[idx]
      lambda[nu, l] <- lambda_values[idx]
      idx <- idx + 1
    }
  }
  
  # Calculando mu com o intercepto
  X_intercept <- cbind(1, X)  # Adicionando a coluna de intercepto
  mu <- exp(X_intercept %*% t(beta))
  
  # Inicializando a verossimilhança
  logL <- 0
  
  # Computando a verossimilhança
  for (i in 1:n) {
    term1 <- sum(y[i, ] * log(m * mu[i, ] / (1 + m * mu[i, ])))
    term2 <- -sum(1 / m * log(1 + m * mu[i, ]))
    term3 <- sum(lgamma(1 / m + y[i, ]) - lgamma(1 / m) - lgamma(y[i, ] + 1))
    
    interaction_term <- 0
    for (l in 1:(d - 1)) {
      for (nu in (l + 1):d) {
        interaction_term <- interaction_term + lambda[l, nu] * 
          (exp(-y[i, l]) - c[l]) * (exp(-y[i, nu]) - c[nu])
      }
    }
    
    logL <- logL + term1 + term2 + term3 + log(1 + interaction_term)
  }
  
  cat("Iteration:", iteration_counter, "\n")
  cat("Current Parameters (Beta, Lambda, m):\n")
  print(list(beta = beta, lambda_values = lambda_values, m = m))
  cat("Current Log-Likelihood:", -logL, "\n\n")
  
  return(-logL)
}

# Otimização com o modelo ajustado para intercepto
result <- optim(init_params, log_likelihood, y = Y, X = X, d = d, k = k, c = c, 
                method = "BFGS", control = list(maxit = 1000, trace = TRUE))

# Extraindo parâmetros otimizados
beta_hat <- matrix(result$par[1:(d * (k + 1))], nrow = d, ncol = k + 1)
lambda_values_hat <- result$par[(d * (k + 1) + 1):(d * (k + 1) + d * (d - 1) / 2)]
m_hat <- result$par[(d * (k + 1) + d * (d - 1) / 2 + 1):length(result$par)]

# Reconstruindo a matriz simétrica lambda
lambda_hat <- matrix(0, d, d)
idx <- 1
for (l in 1:(d - 1)) {
  for (nu in (l + 1):d) {
    lambda_hat[l, nu] <- lambda_values_hat[idx]
    lambda_hat[nu, l] <- lambda_values_hat[idx]
    idx <- idx + 1
  }
}

# Exibindo os resultados
list(beta = beta_hat, lambda = lambda_hat, m = m_hat)


#calculando estimadores

beta_hat_2 <- beta_hat
lambda_values_hat_2 <- lambda_hat
m_hat_2 <- m_hat

# Reconstituir a matriz simétrica lambda a partir dos valores otimizados
lambda_hat_2 <- lambda_hat

mu_hat_2 <- exp(beta_hat_2[, 1] + X %*% t(beta_hat_2[, 2:12]))

par(mfrow = c(2, 2))  # Divide the plotting area into 2x2 for 4 plots

```{r Betas por derivada, beeeem melhor}
inv_x <- ginv(cbind(1,X))
beta_hat_der <- solve(t(cbind(1,X)) %*% cbind(1,X)) %*% t(cbind(1,X)) %*% log(Y + 0.01)

```

