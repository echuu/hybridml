
# logML_approx.R


#### hml_const() ---------------------------------------------------------------
#### this function is called from the candidatePartition.R file to fit the
#### first stage partition
#
hybrid_ml_call = function(u_df) {

    # const_vec  = numeric(N_approx) # store constant approximation
    D = ncol(u_df) - 1
    J = nrow(u_df)

    ## (2) fit the regression tree via rpart()
    u_rpart = rpart::rpart(psi_u ~ ., u_df)

    ## (3) process the fitted tree

    # (3.1) obtain the (data-defined) support for each of the parameters
    param_support = extractSupport(u_df, D) #

    # (3.2) obtain the partition
    u_partition = extractPartition(u_rpart, param_support)


    # param_out = u_star_cand(u_rpart, u_df, u_partition, D) # partition.R
    param_out = u_star(u_rpart, u_df, u_partition, D) # partition.R
    # opt_part = param_out$optimal_part

    # ----------------------------------------------------------------------
    n_partitions = nrow(u_partition) # number of partitions

    # ----------------------------------------------------------------------

    psi_partition = param_out %>%
        dplyr::select(-c('leaf_id', 'psi_choice', 'logQ_cstar', 'n_obs'))

    bounds = psi_partition %>% dplyr::select(-c("psi_star"))

    log_vol_vec = (bounds[seq(2, 2 * D, 2)] - bounds[seq(1, 2 * D, 2)]) %>%
        log() %>% rowSums()

    zhat = (-psi_partition$psi_star + log_vol_vec) %>% log_sum_exp


    return(list(n_partitions  = n_partitions,
                param_out     = param_out,
                u_rpart       = u_rpart,
                param_support = param_support,
                zhat          = zhat))

}
# end of hybrid_ml_call() function --------------------------------------------------


#### logml():
#### Wrapper function that calls the main logml_call() function -- we include
#### this here to catch potential errors during simulations so that rare/
#### problematic draws don't cause the entire simulation to stop running.
#### Replications that result in an error will return NA, which is something we
#### check for in the final simulation results
hybml_const <- function(u_df) {
    out <- tryCatch(
        {
            hybrid_ml_call(u_df)
            # The return value of `readLines()` is the actual value
            # that will be returned in case there is no condition
            # (e.g. warning or error).
            # You don't need to state the return value via `return()` as code
            # in the "try" part is not wrapped insided a function (unlike that
            # for the condition handlers for warnings and error below)
        },
        error=function(cond) {
            message(paste("hybrid computation error"))
            message(cond)
            return(NA)
        },
        warning=function(cond) {
            # message("Here's the original warning message:")
            # message(cond)
            return(NULL)
        },
        finally={
        }
    )
    return(out)
} # end of logml() function ----------------------------------------------------



bridge_approx <- function(samples, log_density, prior, lb, ub, method = 'normal') {
    out <- tryCatch(
        {

            bridge_result <- bridgesampling::bridge_sampler(samples = samples,
                                                            log_posterior = log_density,
                                                            data = prior,
                                                            lb = lb, ub = ub,
                                                            silent = TRUE,
                                                            method = method)
            bridge_result$logml
            # The return value of `readLines()` is the actual value
            # that will be returned in case there is no condition
            # (e.g. warning or error).
            # You don't need to state the return value via `return()` as code
            # in the "try" part is not wrapped insided a function (unlike that
            # for the condition handlers for warnings and error below)
        },
        error=function(cond) {
            message(paste("bridge error"))
            message(cond)
            return(NA)
        },
        warning=function(cond) {
            # message("Here's the original warning message:")
            message(cond)
            return(out)
        },
        finally={
        }
    )
    return(out)
} # end of logml() function ----------------------------------------------------



#### old, unvectorized version of constant hybrid approximation
#### hml_const() ---------------------------------------------------------------
#
#
# hml_const = function(N_approx, D, u_df_full, J, prior) {

#     const_vec  = numeric(N_approx) # store constant approximation

#     # compute approximation to LIL N_approx times
#     for (t in 1:N_approx) {

#         ## (1) subset out rows in u_df_full to be used in the t-th approximation
#         row_id = J * (t - 1) + 1
#         u_df = u_df_full[row_id:(row_id+J-1),]

#         ## (2) fit the regression tree via rpart()
#         u_rpart = rpart(psi_u ~ ., u_df)

#         ## (3) process the fitted tree

#         # (3.1) obtain the (data-defined) support for each of the parameters
#         param_support = extractSupport(u_df, D) #

#         # (3.2) obtain the partition
#         u_partition = extractPartition(u_rpart, param_support)

#         # organize all data into single data frame --> ready for approximation
#         param_out = u_star(u_rpart, u_df, u_partition, D) # partition.R

#         # ----------------------------------------------------------------------
#         n_partitions = nrow(u_partition) # number of partitions

#         # ----------------------------------------------------------------------

#         K = nrow(u_partition)

#         # new declarations here: additional storage for all 3 approximations
#         const_approx   = numeric(K)       # store approx that uses 1-term taylor

#         # declare terms that will be used in the log-sum-exp trick
#         eta_k = numeric(K) # log of the area of each partition A_k

#         ck_1 = numeric(K)

#         # ----------------------------------------------------------------------

#         # (4) compute closed form integral over each partition
#         for (k in 1:n_partitions) {

#             # print(k)

#             # compute the following for log-sum-exp trick
#             # ck_1[k] = -psi(u, prior)

#             ck_1[k] = -param_out$psi_star[k]
#             # ------------------------------------------------------------------
#             for (d in 1:D) {

#                 # find column id of the first lower bound
#                 col_id_lb = grep("u1_lb", names(param_out)) + 2 * (d - 1)
#                 col_id_ub = col_id_lb + 1

#                 # limits of integration, length of the interval for param d
#                 upper = param_out[k, col_id_ub]
#                 lower = param_out[k, col_id_lb]

#                 eta_k[k] = eta_k[k] + log(upper - lower)

#             } # end of loop computing each of 1-dim integrals

#             const_approx[k]  = ck_1[k] + eta_k[k]

#         } # end of for loop over the K partitions

#         # update approximations
#         const_vec[t]  = log_sum_exp(const_approx)


#         psi_star_df = param_out %>% dplyr::select(leaf_id, psi_star)

#         u_df = u_df %>% mutate(leaf_id = u_rpart$where)
#         u_df = merge(u_df, psi_star_df, by = 'leaf_id')

#     } # end of N_approx outer loop


#     return(list(const_vec     = const_vec,
#                 const_approx  = const_approx,   # used to determine logML approx
#                 n_partitions  = n_partitions,
#                 u_df_fit      = u_df,
#                 param_out     = param_out,
#                 u_rpart       = u_rpart,
#                 param_support = param_support))

# }
# end of hml_const() function --------------------------------------------------







# end of logML_approx.R


