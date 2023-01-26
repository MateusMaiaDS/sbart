rm(list=ls())
Rcpp::sourceCpp("src/sbart.cpp")
source("R/wrap_bart.R")
source("R/other_functions.R")
n_ <- 100
x <- matrix(seq(-pi,pi,length.out = n_))
x_new <- matrix(seq(-pi,pi,length.out = n_*10))
y <- sin(x) + rnorm(n = n_,sd = 0.1)
colnames(x) <- "x"
colnames(x_new) <- "x"
# Testing the GP-BART
bart_test <- rbart(x_train = x,y = y,x_test = x_new,n_tree = 10,n_mcmc = 2000,alpha = 0.95,beta = 2,
                  n_burn = 500)


# unit_test_grow(x_train = x,y_train = y)

# dim(bart_test[[2]])


plot(x,y)
points(x_new,apply(bart_test[[2]],1,mean), pch = 20, col = "red")
points(x,apply(bart_test[[1]],1,mean),pch=20)


bartmod <- dbarts::bart(x.train = x,y.train = y,ntree = 20,x.test = x_new)
# points(x,bartmod$yhat.train %>% colMeans(), col = "red",pch=20)
points(x_new,bartmod$yhat.test %>% colMeans(), col = " blue", pch = 20)

# Replicating the same plot on ggplot
ggplot()+
     geom_point(data = data.frame(x = x, y= y ), mapping = aes(x = x, y =y ))+
     geom_point(data = data.frame(x = x, y = apply(bart_test[[1]],1,mean)), mapping = aes(x = x,  y = y), col = "darkblue", alpha = 0.7, pch= 3)+
     geom_line(data = data.frame(x = x_new, y = apply(bart_test[[2]],1,mean)), mapping = aes(x = x, y = y) , col = "blue") +
     geom_line(data = data.frame(x = x_new, y = bartmod$yhat.test.mean), mapping = aes(x = x, y = y), col ="red")+
     theme_bw()

plot(bart_test$tau_post, type = "l")
# microbenchmark::microbenchmark(bart(x_train = x,y_train = y,x_test = x,n_tree = 200,n_mcmc = 5000,
#                                     n_burn = 0,tau = 1,mu = 1,
#                                     tau_mu = 4*4*200,naive_sigma = 1,alpha = 0.95,
#                                     beta = 2,a_tau = 1,d_tau = 1,nsigma = 1),
#                                dbarts::bart(x.train = x,y.train = y,x.test = ,ntree = 200,ndpost = 5000,nskip = 0),times = 10)
