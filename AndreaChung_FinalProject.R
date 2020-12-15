#Create data tables from supplementary data
#Observed Data

#Type 1
df_same_time1 = matrix(data = c(66, 87, 25, 22, 4,
                                13, 14, 15, 9, 4,
                                NA, 4, 4, 9, 1,
                                NA, NA, 4, 3, 1,
                                NA, NA, NA, 1, 1,
                                NA, NA, NA, NA, 0), nrow = 6,
                       ncol = 5, byrow = TRUE)
df_same_time1
rownames(df_same_time1) <- 0:5
colnames(df_same_time1) <- 1:5

df_same_time2 = matrix(data = c(44, 62, 47, 38, 9,
                                10, 13, 8, 11, 5,
                                NA, 9, 2, 7, 3,
                                NA, NA, 3, 5, 1,
                                NA, NA, NA, 1, 0,
                                NA, NA, NA, NA, 1), nrow = 6,
                       ncol = 5, byrow = TRUE)
df_same_time2
rownames(df_same_time2) <- 0:5
colnames(df_same_time2) <- 1:5

#Type 2
df_diff_A = matrix(data = c(9, 12, 18, 9, 4,
                            1, 6, 6, 4, 3,
                            NA, 2, 3, 4, 0,
                            NA, NA, 1, 3, 2,
                            NA, NA, NA, 0, 0,
                            NA,NA,NA,NA,0), nrow = 6,
                   ncol = 5, byrow= T)

rownames(df_diff_A) <- 0:5
colnames(df_diff_A) <- 1:5

df_diff_B = matrix(data = c(15, 12 ,4,
                            11, 17, 4,
                            NA, 21, 4,
                            NA, NA, 5), nrow = 4,
                   ncol = 3, byrow = T)
df_diff_B
rownames(df_diff_B) <- 0:3
colnames(df_diff_B) <- 1:3

#===========================================================================
#abc_sample_generate
abc_sample_generate = function(obs_data1, obs_data2, sum_stat, eps) {
  #randomly chosen parameters -- qc and qh
  q_values1 = runif(n = 2, min = 0, max = 1)
  q_values2 = runif(n = 2, min = 0, max = 1)
  s1 = ncol(obs_data1)
  s2 = ncol(obs_data2)
  y1 = data_gen_fun(q_values1, s1, obs_data1)
  y2 = data_gen_fun(q_values2, s2, obs_data2)
  data = list("y1" = y1, "y2" = y2)
  sum_stat = distance(obs_data1, obs_data2, data)
  sample = if(sum_stat < eps) {
    return(sample)
  }
}

#===========================================================================
#probability matrix
prob_model_matrix = function(q_values, s, data) {
  w = matrix(0, nrow = 6, ncol = s)
  
  #calculating w0s
  w[1,] = sapply(1:s, function(n, theta) {theta^n}, q_values[2])
  
  #calculation of first column
  w[2,1] = 1 - w[1,1] 
  
  #calculation for the rest
  for(i in 2:s) {
    for(j in 2:nrow(w)) {
      if(i >= j) {
        w[j,i] = choose(i, j-1)*w[j, j-1]*(q_values[2] * q_values[1]^(j-1))^(
          i-j+1)
      }else{
        w[j,i] <- 1 - sum(w[,i])
        break
      }
    }
  }
  init = colSumns(data)[1]
  prob_matrix = rmultinom(n = 1, size = init, prob = w[,1])
  for(x in 2:s) {
    next_i = colSums(data)[x]
    prob_calc = rmultinom(n = 1, size = next_i, prob = w[,x])
    prob_matrix <- cbind(prob_matrix, prob_calc)
  }
  return(prob_matrix)
}

#===========================================================================

#data generating function
data_gen_fun = function(q_values, s, data) {
  y = prob_model_matrix(q_values, s, data)
  return(y)
}

y1 = data_gen_fun(q_values1, s1, data1)
y2 = data_gen_fun(q_values2, s2, data2)
gen_data = list("y1" = y1, "y2" = y2)

#===========================================================================
#Summary Statistics Function

distance = function(data1, data2, gen_data) {
  return(0.5 * (norm(data1 - gen_data$y1, type = "F") + norm(data2 - gen_data$y2, type = "F")))
}

#===========================================================================

#posterior samples
posterior_dist = replicate(n = 10000, 
                           abc_sample_generate(obs_data1, 
                                               obs_data2, sum_stat, eps))
