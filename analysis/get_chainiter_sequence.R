##  Script to generate random sequence of chains and iterations for simulations

nchains <- 3
niters <- 1000
nsims <- 2500

randchain <- sample(x = c(1:nchains), size = nsims, replace = TRUE)
randiter <- sample(x = c(1:niters), size = nsims, replace = TRUE)

randchainiters <- data.frame(chain=randchain, iter=randiter)
saveRDS(randchainiters, "random_chain_iter_sequence.RDS")
