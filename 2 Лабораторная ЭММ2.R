if (!require(tseries)) 
{
  install.packages("tseries")
}
library(tseries)
set.seed(123)

# Функция для моделирования процесса GARCH(1,0)
garch_process_1_0 <- function(a0, a1, num_obs) 
{
  if (a0 <= 0) stop("a0 должно быть положительным")
  if (a1 <= 0 || a1 >= 1) stop("a1 должно быть в диапазоне (0, 1)")
  h <- numeric(num_obs)
  process_values <- numeric(num_obs)
  epsilon <- rnorm(num_obs)
  h[1] <- a0 / (1 - a1)  # Начальная волатильность
  for (t in 2:num_obs) 
  {
    h[t] <- a0 + a1 * epsilon[t - 1]^2 * h[t - 1]
    process_values[t] <- sqrt(h[t]) * epsilon[t]
  }
  par(mfrow = c(2, 1))
  plot(process_values, type = "l", col = "seagreen", main = "Процесс GARCH(1,0)", ylab = "h_n", xlab = "Время")
  plot(sqrt(h), type = "l", col = "blue", main = "Волатильность GARCH(1,0)", ylab = "σ_n", xlab = "Время")
  return(list(process_values = process_values, volatility = sqrt(h)))
}

# Функция для моделирования процесса GARCH(3,0)
garch_process_3_0 <- function(a0, a1, a2, a3, num_obs) 
{
  if (a0 <= 0) stop("a0 должно быть положительным")
  if (a1 <= 0 || a2 <= 0 || a3 <= 0 || (a1 + a2 + a3) >= 1) stop("Сумма a1, a2, a3 должна быть меньше 1")
  h <- numeric(num_obs)
  process_values <- numeric(num_obs)
  epsilon <- rnorm(num_obs)
  h[1:3] <- a0 / (1 - a1 - a2 - a3)  # Начальная волатильность
  for (t in 4:num_obs) 
  {
    h[t] <- a0 + a1 * (process_values[t - 1])^2 + a2 * (process_values[t - 2])^2 + a3 * (process_values[t - 3])^2
    process_values[t] <- sqrt(h[t]) * epsilon[t]
  }
  par(mfrow = c(2, 1))
  plot(process_values, type = "l", col = "seagreen", main = "Процесс GARCH(3,0)", ylab = "h_n", xlab = "Время")
  plot(sqrt(h), type = "l", col = "blue", main = "Волатильность GARCH(3,0)", ylab = "σ_n", xlab = "Время")
  return(list(process_values = process_values, volatility = sqrt(h)))
}

# Функция для моделирования процесса GARCH(1,1)
garch_process_1_1 <- function(a0, a1, b1, num_obs) 
{
  if (a0 <= 0) stop("a0 должно быть положительным")
  if (a1 <= 0 || b1 <= 0 || (a1 + b1) >= 1) stop("a1 и b1 должны быть положительными и a1 + b1 < 1")
  h <- numeric(num_obs)
  process_values <- numeric(num_obs)
  epsilon <- rnorm(num_obs)
  h[1] <- a0 / (1 - a1 - b1)  # Начальная волатильность
  for (t in 2:num_obs) 
  {
    h[t] <- a0 + a1 * (epsilon[t - 1])^2 + b1 * h[t - 1]
    process_values[t] <- sqrt(h[t]) * epsilon[t]
  }
  par(mfrow = c(2, 1))
  plot(process_values, type = "l", col = "seagreen", main = "Процесс GARCH(1,1)", ylab = "h_n", xlab = "Время")
  plot(sqrt(h), type = "l", col = "blue", main = "Волатильность GARCH(1,1)", ylab = "σ_n", xlab = "Время")
  return(list(process_values = process_values, volatility = sqrt(h)))
}

# Параметры для GARCH(1,0), GARCH(3,0) и GARCH(1,1)
n <- 1100
a0 <- 0.1
a1 <- 0.4

# Генерация данных для GARCH(1,0)
garch_1_0 <- garch_process_1_0(a0, a1, n)

# Параметры для GARCH(3,0)
a1_3_0 <- 0.3
a2_3_0 <- 0.2
a3_3_0 <- 0.1

# Генерация данных для GARCH(3,0)
garch_3_0 <- garch_process_3_0(a0, a1_3_0, a2_3_0, a3_3_0, n)

# Параметры для GARCH(1,1)
b1 <- 0.5

# Генерация данных для GARCH(1,1)
garch_1_1 <- garch_process_1_1(a0, a1, b1, n)

# Прогноз для GARCH(3,0) ----
train_size <- 1000
train_data <- garch_3_0$process_values[1:train_size]
test_data <- garch_3_0$process_values[(train_size + 1):n]

# Оценка параметров GARCH(3,0) на обучающей выборке
garch_fit_3_0 <- garch(train_data, order = c(0, 3), trace = FALSE)
cat("Оценка параметров GARCH(3,0):\n")
print(coef(garch_fit_3_0))

# Функция для последовательного прогнозирования GARCH(3,0)
sequential_garch_forecast <- function(garch_fit, data, n_ahead) {
  params <- coef(garch_fit)  # Оцененные параметры
  a0 <- params["a0"]
  a1 <- params["a1"]
  a2 <- params["a2"]
  a3 <- params["a3"]
  
  n <- length(data)
  h <- numeric(n_ahead)  # Вектор для прогнозов волатильности
  forecast_values <- numeric(n_ahead)  # Вектор для прогнозов значений процесса
  epsilon <- rnorm(n_ahead)  # Генерация случайных ошибок для прогноза
  
  # Инициализация последними значениями из обучающей выборки
  last_values <- tail(data, 3)^2  # Квадраты последних трех наблюдений
  h[1] <- a0 + a1 * last_values[3] + a2 * last_values[2] + a3 * last_values[1]
  forecast_values[1] <- sqrt(h[1]) * epsilon[1]
  
  # Последовательное прогнозирование
  for (t in 2:n_ahead) {
    if (t == 2) {
      h[t] <- a0 + a1 * forecast_values[t-1]^2 + a2 * last_values[3] + a3 * last_values[2]
    } else if (t == 3) {
      h[t] <- a0 + a1 * forecast_values[t-1]^2 + a2 * forecast_values[t-2]^2 + a3 * last_values[3]
    } else {
      h[t] <- a0 + a1 * forecast_values[t-1]^2 + a2 * forecast_values[t-2]^2 + a3 * forecast_values[t-3]^2
    }
    forecast_values[t] <- sqrt(h[t]) * epsilon[t]
  }
  
  return(list(forecast_volatility = sqrt(h), forecast_values = forecast_values))
}

# Выполнение прогноза на 100 шагов вперед
n_ahead <- 100
forecast <- sequential_garch_forecast(garch_fit_3_0, train_data, n_ahead)

# График квадратов наблюдений и прогнозов волатильности
par(mfrow = c(1, 1))
plot(test_data[1:n_ahead]^2, type = "l", col = "blue", 
     main = "Сравнение квадратов наблюдений и прогнозов GARCH(3,0)", 
     ylab = "Квадраты значений / Волатильность", xlab = "Время")
lines(forecast$forecast_volatility^2, col = "red", lty = 2)
legend("topright", legend = c("Квадраты наблюдений", "Прогноз волатильности"), 
       col = c("blue", "red"), lty = c(1, 2))