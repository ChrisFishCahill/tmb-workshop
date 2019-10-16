#data:
adults    <- c(2.4,3.2,3.9,5.7,6,7.4,8.2,10,10.1,10.4,11.3,12.8,18,24)
juveniles <- c(11.6,7.1,14.3,19.1,12.4,19.7,31.5,18.5,22.1,26.9,19.2,21,18.1,26.8)
temps     <- c(12.00,15.04,13.45,12.03,13.00,13.48,9.24,15.33,12.08,9.50,17.10,12.16,22.00,13.64)
data      <- data.frame(adults=adults, juveniles=juveniles, temps=temps)
Nobs      <- nrow(data)
plot(data)
