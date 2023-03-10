rm(list=ls())
library(tidyverse)

Rcpp::sourceCpp("src/sbart.cpp")
source("R/wrap_bart.R")
source("R/other_functions.R")
n_ <- 100
x <- matrix(seq(-pi,pi,length.out = n_))
x_new <- matrix(seq(-pi,pi,length.out = n_*10))
y <- sin(x) + rnorm(n = n_,sd = 0.1)
y[x<0] <- y[x<0] + 2
y[x>0] <- y[x>0] - 2

colnames(x) <- "x"
colnames(x_new) <- "x"

# Testing over the motorbike data
library(boot)
data("motor")
x <- motor$times %>% as.matrix
y <- motor$accel %>% as.matrix()
x_new <- seq(min(x),max(x),length.out = 500) %>% as.matrix()
colnames(x) <- "x"
colnames(x_new) <- "x"

# Testing the GP-BART
bart_test <- rbart(x_train = x,y = y,x_test = x_new,n_tree = 10,n_mcmc = 2000,
                   alpha = 0.95,beta = 2,
                   n_burn = 500,scale_bool = TRUE)



bartmod <- dbarts::bart(x.train = x,y.train = y,ntree = 20,x.test = x_new)

# All ll trees prediction plot
all_tree_posterior_mean <- Reduce("+",bart_test$all_tree_post)/length(bart_test$all_tree_post)
if(!(ncol(all_tree_posterior_mean) %>% is_null())){
        colnames(all_tree_posterior_mean) <- paste0("tree.",1:ncol(all_tree_posterior_mean))
} else {
        all_tree_posterior_mean <- all_tree_posterior_mean %>% as.matrix()
        colnames(all_tree_posterior_mean) <- "tree.1"
}
all_tree_posterior_mean <- all_tree_posterior_mean %>% as.data.frame %>% add_column(x) %>% pivot_longer(starts_with("tree"))


# Getting quantiles for ribon
pi_ <- bart_test$y_hat_test %>% apply(1,function(x){quantile(x,probs = c(0.025,0.975))}) %>% t %>% cbind(x_new)
colnames(pi_) <- c("lower","upper","x")
pi_ <- pi_ %>% as.data.frame()
# Replicating the same plot on ggplot
ggplot()+
     geom_ribbon(data = pi_,mapping = aes(x = x,ymin = lower, ymax = upper), col = NA, alpha = 0.1, lty = "dashed",fill = "blue")+
     geom_point(data = data.frame(x = x, y= y ), mapping = aes(x = x, y =y ))+
     geom_point(data = data.frame(x = x, y = apply(bart_test[[1]],1,mean)), mapping = aes(x = x,  y = y), col = "darkblue", alpha = 0.7, pch= 3)+
     geom_line(data = data.frame(x = x_new, y = apply(bart_test[[2]],1,mean)), mapping = aes(x = x, y = y) , col = "blue") +
     geom_line(data = data.frame(x = x_new, y = bartmod$yhat.test.mean), mapping = aes(x = x, y = y), col ="red")+
     geom_line(data = all_tree_posterior_mean,
               mapping = aes(x = x, y = value, col = name), alpha = 0.5,show.legend = FALSE)+
     ylim(y = range(y)*1.4)+
     theme_bw()

# plot(bart_test$tau_post, type = "l")
# microbenchmark::microbenchmark(bart(x_train = x,y_train = y,x_test = x,n_tree = 200,n_mcmc = 5000,
#                                     n_burn = 0,tau = 1,mu = 1,
#                                     tau_mu = 4*4*200,naive_sigma = 1,alpha = 0.95,
#                                     beta = 2,a_tau = 1,d_tau = 1,nsigma = 1),
#                                dbarts::bart(x.train = x,y.train = y,x.test = ,ntree = 200,ndpost = 5000,nskip = 0),times = 10)
plot(bart_test$tau_b_post[-length(bart_test$tau_b_post)], type = "l")
