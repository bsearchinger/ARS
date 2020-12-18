library(dplyr)
library(ggplot2)
source('./R/modular_functions.R')

# replicate gilks & wild simulation
x1 <- c(-.5, -1, -2, -5, -10, -9, -8, -7, -6)

x2 <- c(.5, 1, 2, 5, 10, 1, 2, 3, 4)

df <- data.frame(x1 = double(), x2 = double(), x_abs_added = integer(),
                     total_iters = integer(), samples_rejected = integer())
set.seed(1)

for(p in 1:length(x1)){
  for(q in 1:500){
    out <- ars(1, c(x1[p], x2[p]), dnorm, c(mean = 0, sd = 1), c(-Inf,Inf))
    k <- unname(out[2])
    i <- unname(out[3])
    m <- unname(out[4])
    df <- rbind(df, c(x1 = x1[p], x2 = x2[p], x_abs_added = k,
                      total_iters = i, samples_rejected = m))
    q <- q + 1
  }
}

res <- df %>% group_by(x1, x2) %>% summarise(mean_iters = mean(total_iters),
                                             max_iters = max(total_iters),
                                             mean_rej = mean(samples_rejected),
                                             max_rej = max(samples_rejected),
                                             mean_added = mean(x_abs_added),
                                             max_added = max(x_abs_added))

write.table(res, "./tests/table.csv")
