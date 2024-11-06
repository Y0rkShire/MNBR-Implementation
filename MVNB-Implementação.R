inv_x <- ginv(cbind(1,X))
beta_hat_der <- solve(t(cbind(1,X)) %*% cbind(1,X)) %*% t(cbind(1,X)) %*% log(Y + 0.001)


k <- ncol(X)  # Inclui o intercepto como parte dos betas
n <- nrow(Y)  # Número de observações
d <- ncol(Y)  # Número de variáveis de resposta

c <- colMeans(exp(-Y))

empirical_means <- colMeans(Y)
empirical_variances <- apply(Y, 2, var)

# Valores iniciais para os parâmetros, agora com intercepto
init_beta <- matrix(NA, nrow = d, ncol = k + 1)  # Inclui coluna para o intercepto

for (t in 1:d) {
  nb_model <- glm.nb(Y[, t] ~ X)  # Ajuste do modelo NB com intercepto
  init_beta[t, ] <- coef(nb_model)  # Armazena os coeficientes, incluindo o intercepto
}

init_lambda <- rep(1e-5, d * (d - 1) / 2)  # Inicializando lambda

init_m <- numeric(d)
for (t in 1:d) {
  theta_t <- 1 - empirical_means[t] / empirical_variances[t]  # Estimativa de theta_t
  init_m[t] <- empirical_means[t] * (1 - theta_t) / theta_t   # Estimativa de m_t
}

init_params <- c(init_lambda, init_m)

iteration_counter <<- 0  # Variável global para contagem de iterações

# Definindo a função de verossimilhança
log_likelihood <- function(params, y, X, d, k, c, beta_hat) {
  iteration_counter <<- iteration_counter + 1
  
  # Extract lambda and m from parameters
  lambda_values <- params[1:(d * (d - 1) / 2)]
  m <- params[(d * (d - 1) / 2 + 1):length(params)]
  
  # Build symmetric lambda matrix
  lambda <- matrix(0, d, d)
  idx <- 1
  for (l in 1:(d - 1)) {
    for (nu in (l + 1):d) {
      lambda[l, nu] <- lambda_values[idx]
      lambda[nu, l] <- lambda_values[idx]
      idx <- idx + 1
    }
  }
  
  # Print lambda and m values at the start of each iteration
  cat("Iteration:", iteration_counter, "\n")
  cat("Lambda matrix:\n")
  print(lambda)
  cat("m values:\n")
  print(m)
  
  # Calculate mu using fixed beta
  X_intercept <- cbind(1, X)
  mu <- exp(X_intercept %*% beta_hat)
  
  # Initialize log-likelihood
  logL <- 0
  
  # Compute log-likelihood
  for (i in 1:n) {
    term1 <- sum(y[i, ] * log(m * mu[i, ] / (1 + m * mu[i, ])))
    term2 <- -sum(1 / m * log(1 + m * mu[i, ]))
    term3 <- sum(lgamma(1 / m + y[i, ]) - lgamma(1 / m) - lgamma(y[i, ] + 1))
    
    interaction_term <- 0
    for (l in 1:(d - 1)) {
      for (nu in (l + 1):d) {
        interaction_term <- interaction_term + lambda[l, nu] * (exp(-y[i, l]) - c[l]) * (exp(-y[i, nu]) - c[nu])
      }
    }
    if(interaction_term < -1){
      interaction_term <- -1 + 1e-5
    }
    
    logL <- logL + term1 + term2 + term3 + log(1 + interaction_term)
  }
  
  cat("- Log-Likelihood:", -logL, "\n")
  
  return(-logL)
}

# Parâmetros iniciais ajustados (excluindo beta)
init_params <- c(init_lambda, init_m)

# Limites inferiores ajustados (excluindo beta)
lower_bounds <- rep(1e-5, length(init_params))

# Rodar otimização apenas para lambda e m
result <- optim(
  init_params, log_likelihood, y = Y, X = X, d = d, k = k, c = c, beta_hat = beta_hat_der,
  method = "L-BFGS-B", control = list(maxit = 1000, trace = TRUE),
  lower = lower_bounds
)

# Extrair parâmetros otimizados
lambda_values_hat <- result$par[1:(d * (d - 1) / 2)]
m_hat <- result$par[(d * (d - 1) / 2 + 1):length(result$par)]

# Reconstruir matriz simétrica lambda
lambda_hat <- matrix(0, d, d)
idx <- 1
for (l in 1:(d - 1)) {
  for (nu in (l + 1):d) {
    lambda_hat[l, nu] <- lambda_values_hat[idx]
    lambda_hat[nu, l] <- lambda_values_hat[idx]
    idx <- idx + 1
  }
}

# Imprimir resultados finais
cat("Optimization completed\n")
cat("Final Iteration Count:", iteration_counter, "\n")
cat("Final Negative Log-Likelihood:", result$value, "\n\n")

# Exibir parâmetros otimizados (beta fixo)
list(beta = beta_hat_der, lambda = lambda_hat, m = m_hat)


#calculando estimadores

beta_hat_2 <- beta_hat
lambda_values_hat_2 <- lambda_hat
m_hat_2 <- m_hat

# Reconstituir a matriz simétrica lambda a partir dos valores otimizados
lambda_hat_2 <- lambda_hat

mu_hat_2 <- exp(beta_hat_2[, 1] + X %*% t(beta_hat_2[, 2:12]))



AIC_params <- c(result$par)

verossimilhancas <- c()
X <- X_orig[,-5]
for(K in 1:k - 1){
  iteration_counter <<- 0
  inv_x <- ginv(cbind(1,X[,-K]))
  beta_hat_AIC <- solve(t(cbind(1,X[,-K])) %*% cbind(1,X[,-K])) %*% t(cbind(1,X[,-K])) %*% log(Y + 0.001)
  verossimilhancas <- c(verossimilhancas,
                        optim(
    AIC_params, log_likelihood, y = Y, X = X[,-K], d = d, k = k - 2, c = c, beta_hat = beta_hat_AIC,
    method = "L-BFGS-B", control = list(maxit = 1000, trace = TRUE),
    lower = lower_bounds
  )$value*2 + 20)
}

```{r Betas por derivada, beeeem melhor}
inv_x <- ginv(cbind(1,X))
beta_hat_der <- solve(t(cbind(1,X)) %*% cbind(1,X)) %*% t(cbind(1,X)) %*% log(Y + 0.01)

```

