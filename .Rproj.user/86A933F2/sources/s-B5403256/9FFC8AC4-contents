

# source("C:/Users/ericc/mlike_approx/algo/setup.R")
source("C:/Users/ericc/mlike_approx/covariance/HIW_helper.R")
# source("C:/Users/ericc/mlike_approx/EP/hybml.R")

testG  = matrix(c(1,1,0,0,0,
                  1,1,1,1,1,
                  0,1,1,1,1,
                  0,1,1,1,1,
                  0,1,1,1,1), 5, 5)

testG = matrix(c(1,1,1,0,0,
                 1,1,1,0,0,
                 1,1,1,1,1,
                 0,0,1,1,1,
                 0,0,1,1,1), 5,5)

# Given graph testG
testG = matrix(c(1,1,0,0,1,0,0,0,0,
                 1,1,1,1,1,0,0,0,0,
                 0,1,1,1,0,0,0,0,0,
                 0,1,1,1,1,1,1,0,0,
                 1,1,0,1,1,1,0,0,0,
                 0,0,0,1,1,1,1,1,1,
                 0,0,0,1,0,1,1,1,1,
                 0,0,0,0,0,1,1,1,1,
                 0,0,0,0,0,1,1,1,1), 9, 9)

a = c(1, 3, 2, 5, 4, 6, 7, 8, 9)
testG = testG[a, a]

D = nrow(testG)
b = 3          # prior degrees of freedom
V = diag(1, D) # prior scale matrix

D_0 = 0.5 * D * (D + 1) # num entries on diagonal and upper diagonal
J = 1000

log_density = function(u, data) {
    -psi(u, data)
}

# hme_exp_term = function(u) {
#    Lt = matrix(0, D, D)     # (D x D) lower triangular matrix
#    Lt_vec_0 = numeric(D_0)  # (D_0 x 1) vector to fill upper triangular, Lt
#    Lt_vec_0[edgeInd] = u
#    Lt[upper.tri(Lt, diag = T)] = Lt_vec_0   # populate lower triangular terms
#
#    - N * log_det(Lt) + 0.5 * matrix.trace(t(Lt) %*% Lt %*% S)
# }

# logical vector determining existence of edges between vertices
edgeInd = testG[upper.tri(testG, diag = TRUE)] %>% as.logical
upperInd = testG[upper.tri(testG)] %>% as.logical
D_u = sum(edgeInd)

# Specify true value of Sigma
set.seed(1)
true_params = HIWsim(testG, b, V)
Sigma_G = true_params$Sigma
Omega_G = true_params$Omega # precision matrix -- this is the one we work with

# chol(Omega_G)

# Generate data Y based on Sigma_G
N = 100
Y = matrix(0, N, D)
for (i in 1:N) {
    Y[i, ] = t(t(chol(Sigma_G)) %*% rnorm(D, 0, 1)) # (500 x D)
}

S = t(Y) %*% Y

params = list(N = N, D = D, D_0 = D_0, testG = testG, edgeInd = edgeInd,
              upperInd = upperInd, S = S, V = V, b = b)

J = 5000
postIW = sampleHIW(J, D_u, D_0, testG, b, N, V, S, edgeInd)
post_samps = postIW$post_samps                 # (J x D_u)

u_df = preprocess(post_samps, D_u, params)     # J x (D_u + 1)

(LIL = logmarginal(Y, testG, b, V, S))

lambda = function(u, params) { pracma::grad(psi, u, params = params) }
hess   = function(u, params) { pracma::hessian(psi, u, params = params) }

hybml_const(u_df)$zhat
hybml(u_df, params, grad = lambda, hess = hess)
(LIL = logmarginal(Y, testG, b, V, S))
- 0.5 * D * N * log(2 * pi) +
    gnorm(testG, b + N, V + S, iter = 1000) - gnorm(testG, b, V, iter = 1000)



gnorm_approx = numeric(20)
approx = numeric(20)
for (i in 1:20) {
    print(i)
    postIW = sampleHIW(J, D_u, D_0, testG, b, N, V, S, edgeInd)
    post_samps = postIW$post_samps                 # (J x D_u)

    u_df = preprocess(post_samps, D_u, params)     # J x (D_u + 1)
    gnorm_approx[i] = - 0.5 * D * N * log(2 * pi) +
        gnorm(testG, b + N, V + S, iter = 5000) - gnorm(testG, b, V, iter = 5000)
    approx[i] = hyb(u_df, psi, params)

}
(gnorm_approx[i] = - 0.5 * D * N * log(2 * pi) +
        gnorm(testG, b + N, V + S, iter = 1000) - gnorm(testG, b, V, iter = 1000))

# hybrid = hybrid_ml(D_u, u_df, J, params)
# hybrid$zhat

approx[i] = hyb(u_df, psi, params)

mean(abs(LIL - gnorm_approx))
mean(abs(LIL - approx))

sqrt(mean((LIL - gnorm_approx)^2))




# u_samp = as.matrix(post_samps)
# colnames(u_samp) = names(u_df)[1:D_u]
# # prepare bridge_sampler input()
# lb = rep(-Inf, D_u)
# ub = rep(Inf, D_u)
# names(lb) <- names(ub) <- colnames(u_samp)
#
# bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
#                                                log_posterior = log_density,
#                                                data = params,
#                                                lb = lb, ub = ub,
#                                                silent = TRUE)
# bridge_result$logml



hme_approx = function(u_df, params, J, D, N) {

    N   = params$N
    D   = params$D
    D_0 = params$D_0
    S   = params$S

    tmp_J = apply(u_df[,1:D_u], 1, hme_exp_term) # term in the exponential
    log(J) - 0.5 * N * D * log(2 * pi) - log_sum_exp(tmp_J)
}

hme_approx(u_df, params, J, D, N)


B = 100 # number of replications
hyb = numeric(B)
bridge = numeric(B)
hme = numeric(B)
set.seed(123)
J = 25
i = 1


while (i <= B) {
    i = 11

    postIW = sampleHIW(J, D_u, D_0, testG, b, N, V, S, edgeInd)
    post_samps = postIW$post_samps                 # (J x D_u)

    u_df = preprocess(post_samps, D_u, params)     # J x (D_u + 1)

    hybrid = hybrid_ml(D_u, u_df, J, param_list)

    if (any(is.na(hybrid))) {print(paste("error in iteration", i)); next;}

    hyb[i] = hybrid$zhat

    u_samp = as.matrix(post_samps)
    colnames(u_samp) = names(u_df)[1:D_u]
    lb = rep(-Inf, D_u)
    ub = rep(Inf, D_u)
    names(lb) <- names(ub) <- colnames(u_samp)
    bridge_result = bridgesampling::bridge_sampler(samples = u_samp,
                                                   log_posterior = log_density,
                                                   data = params,
                                                   lb = lb, ub = ub,
                                                   silent = TRUE,
                                                   method = 'warp3')

    bridge[i] = bridge_result$logml
    hme[i] = hme_approx(u_df, params, J, D, N)

    avg_hyb = mean(hyb[hyb!=0])
    avg_bridge = mean(bridge[bridge!=0])
    avg_hme = mean(hme[hme!=0])

    print(paste("iter ", i, ': ',
                "hybrid = ", round(avg_hyb, 3),
                '; ', "ae = ", round(mean((LIL - hyb[hyb!=0])), 4),
                ' // ',
                "bridge = ", round(avg_bridge, 3), '; ',
                "ae = ", round(mean((LIL - bridge[bridge!=0])), 4),
                ' // ',
                "hme = ", round(avg_hme, 3), '; ',
                "ae = ", round(mean((LIL - hme[hme!=0])), 4),
                sep = ''))
    i = i + 1
}


set.seed(1)
gnorm_approx = numeric(B)
for (i in 1:B) {
    gnorm_approx[i] = - 0.5 * D * N * log(2 * pi) +
        gnorm(testG, b + N, V + S, iter = 13) -
        gnorm(testG, b, V, iter = 12)
}


approx = data.frame(LIL, hyb = hyb[1:B], bridge = bridge[1:B],
                    hme = hme[1:B], gnorm_approx = gnorm_approx)

names(approx) = c("LIL", "HybE", "WBSE", "HME", "GNORM")


(error = data.frame(approx = colMeans(approx), approx_sd = apply(approx, 2, sd),
            mae = colMeans(abs(LIL - approx)),
            ae = colMeans(LIL - approx),
            rmse = sqrt(colMeans((LIL - approx)^2))) %>% round(3))


delta_df = (LIL - approx) %>% dplyr::select(-c('LIL')) %>% melt()

delta_df %>% head

### figure dimensions for saving: (954 x 488)
ggplot(delta_df, aes(x = variable, y = value)) + geom_boxplot() +
    coord_flip() +
    labs(y = expression(paste(Delta, ' ', ln, ' ', p(y))), x = '') +
    theme_bw() +
    theme(axis.text  = element_text(size=25),
          axis.title = element_text(size=25,face="bold")) +
    scale_y_continuous(breaks = seq(-60, 0, 10))


fig_loc = 'C:/Users/ericc/mlike_approx/final_sims/'
saveRDS(list(J = J, D = D, D_u = D_u, N = N, approx = approx, error = error),
        file = paste(fig_loc, 'hiw_d5_j25.RData', sep = ''))
hiw_d5 = readRDS(paste(fig_loc, 'hiw_d5_j25_final.RData', sep = ''))

approx = hiw_d5$approx

hiw_d5$error

delta_df = (LIL - approx) %>% dplyr::select(-c('LIL')) %>% melt()

delta_df %>% head

### figure dimensions for saving: (954 x 488)
ggplot(delta_df, aes(x = variable, y = value)) + geom_boxplot() +
    geom_hline(yintercept = 0, col = 'red', size = 1, linetype = 'dashed') +
    coord_flip() +
    labs(y = expression(paste(Delta, ' ', ln, ' ', p(y))), x = '') +
    theme_bw() +
    theme(axis.text  = element_text(size=25),
          axis.title = element_text(size=25,face="bold")) +
    scale_y_continuous(breaks = seq(-60, 0, 10))




hiw_5e3 = readRDS(paste(fig_loc, 'hiw_d5_j5e3.RData', sep = ''))
hiw_5e3$error

approx_5e3 = hiw_5e3$approx_df
approx_5e3 %>% head
names(approx_5e3) = c("LIL", "HybE", "BSE", "HME")


delta_df = (LIL - approx_5e3) %>% dplyr::select(-c('LIL')) %>% melt()

delta_df %>% head

### figure dimensions for saving: (954 x 488)
ggplot(delta_df, aes(x = variable, y = value)) + geom_boxplot() +
    coord_flip() +
    labs(y = expression(paste(Delta, ' ', ln, ' ', p(y))), x = '') +
    theme_bw() +
    theme(axis.text  = element_text(size=25),
          axis.title = element_text(size=25,face="bold")) +
    scale_y_continuous(breaks = seq(-60, 0, 10))






