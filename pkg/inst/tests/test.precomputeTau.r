source("../utils.r")

M <- 7
N <- 5
times <- c(1,2,3,4,5,6,7)
sen <- c(1,3,3,1,2,5,2)
rec <- c(3,1,1,3,5,1,4)

tau <- precomputeTau(cbind(times,sen,rec),N)
expect_that(all(tau[1,,] == 0), is_true())
expect_that(tau[2,1,3],equals(0))
expect_that(tau[4,1,3],equals(2))
expect_that(tau[6,1,4],equals(3))
expect_that(tau[6,1,4],equals(3))
expect_that(tau[6,1,4],equals(3))
expect_that(tau[7,2,3],equals(4))
expect_that(tau[7,2,4],equals(4))
expect_that(tau[7,2,1],equals(5))
