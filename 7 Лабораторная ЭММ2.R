# Параметры
set.seed(123)
lambda <- 2  # Средняя частота страховых случаев
t <- 50      # Период времени
num_simulations <- 1000

# Моделирование количества страховых случаев
N_t <- rpois(num_simulations, lambda * t)

# Построение гистограммы
hist(N_t, breaks=20, probability=TRUE, main="Гистограмма N_t",
     xlab="Количество страховых случаев", ylab="Частота")

# Добавим теоретическое распределение Пуассона
curve(dpois(x, lambda * t), add=TRUE, col="red", lwd=2)
legend("topright", legend=c("Эмпирическое", "Теоретическое"), col=c("black", "red"), lty=c(1, 1))
# Параметры для модели капитала компании
U0 <- 50       # Начальный капитал
c <- 1         # Доход от премий в единицу времени
lambda <- 0.3  # Параметр распределения Пуассона
mu <- 3        # Средняя выплата
t_max <- 100   # Время для моделирования

# Реализация процесса капитала для заданного времени
simulate_ruin_process <- function(U0, c, lambda, mu, t_max) {
  U <- U0
  times <- c(0)
  capital <- c(U0)
  
  while (TRUE) {
    tau <- rexp(1, rate=lambda)
    t <- sum(times) + tau
    
    if (t > t_max) {
      # ДОБАВЛЯЕМ капитал до конца периода, даже если убытка уже не будет
      U <- U + c * (t_max - sum(times))  # Доход до конца периода
      times <- c(times, t_max)
      capital <- c(capital, U)
      break
    }
    
    X <- rexp(1, rate=1/mu)
    U <- U + c * tau - X
    times <- c(times, t)
    capital <- c(capital, U)
    
    if (U < 0) break
  }
  
  list(times=times, capital=capital)
}


# Случай 1: условие разорения выполняется
res1 <- simulate_ruin_process(U0, c, 0.3, 3, t_max)
plot(res1$times, res1$capital, type="l", main="Процесс капитала компании (условие разорения выполняется)",
     xlab="Время", ylab="Капитал")

# Случай 2: условие разорения не выполняется
res2 <- simulate_ruin_process(U0, c, 0.1, 3, t_max)
plot(res2$times, res2$capital, type="l", main="Процесс капитала компании (условие разорения не выполняется)",
     xlab="Время", ylab="Капитал")
# Параметры для оценки вероятности разорения
num_simulations <- 1000
U0 <- 100
c <- 1
lambda <- 0.3
mu <- 3
t_max <- 1000

# Функция для проверки разорения
check_ruin <- function() {
  res <- simulate_ruin_process(U0, c, lambda, mu, t_max)
  return(any(res$capital < 0))
}

# Проведение симуляций и расчет вероятности разорения
ruin_results <- replicate(num_simulations, check_ruin())
ruin_probability <- mean(ruin_results)

cat("Выборочная вероятность разорения:", ruin_probability, "\n")
# Расчет теоретической вероятности разорения по условию Лундберга
rho <- c / (lambda * mu) - 1
psi_theoretical <- exp(-U0 * (rho / (1 + rho)) / mu)

cat("Теоретическая вероятность разорения по условию Лундберга:", psi_theoretical, "\n")
