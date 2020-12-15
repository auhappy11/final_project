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
abc_sample_generate = function(obs_data1, obs_data2, eps) {
  #randomly chosen parameters -- qc and qh
  q_values1 <- runif(n = 2, min = 0, max = 1)
  q_values2 <- runif(n = 2, min = 0, max = 1)
  y1 <- data_gen_fun(q_values1, ncol(obs_data1), obs_data1)
  y2 <- data_gen_fun(q_values2, ncol(obs_data2), obs_data2)
  data <- list("y1" = y1, "y2" = y2)
  sum_stat = distance(obs_data1, obs_data2, data)
  sample = while(TRUE) {
    if(sum_stat <= eps) {
    return(list(round(q_values1[1], 5), round(q_values1[2], 5), 
                round(q_values2[1], 5), round(q_values2[2], 5)))
    }
  }
  return(sample)
}


#testing

abc_sample_generate(df_same_time1, df_same_time2, eps = 100)
abc_sample_generate(df_diff_A, df_diff_B, eps = 50)

#===========================================================================
#probability matrix
prob_model_matrix = function(q_values, obs_data) {
  obs_data[is.na(obs_data) <- 0]
  s = ncol(obs_data)
  w = matrix(0, nrow = nrow(obs_data), ncol = s)
  
  #calculating w0s
  w[1,] <- sapply(1:s, function(n) (q_values[2])^n)
  
  #calculation of first column
  w[2,1] <- 1 - w[1,1] 
  #print(w)
  #calculation for the rest
  for(i in 2:s) {
    for(j in 2:nrow(w)) {
      if(j <= i) {
        w[j,i] = choose(i, j-1) * (w[j, j-1]) * ((q_values[2]) * 
                                             (q_values[1])^(j-1))^(i-j+1)
      }else{
        w[j,i] <- 1 - sum(w[,i])
        break
      }
    }
  } 
  init = colSums(obs_data, na.rm = T)[1]
  prob_matrix = rmultinom(n = 1, size = init, prob = w[,1])
  for(k in 2:s) {
    next_i = colSums(obs_data, na.rm = T)[k]
    prob_calc = rmultinom(n = 1, size = next_i, prob = w[,k])
    prob_matrix <- cbind(prob_matrix, prob_calc)
  }
  return(prob_matrix)
}

#testing
prob_model_matrix(runif(2), df_same_time1)
prob_model_matrix(runif(2), df_same_time2)

#===========================================================================

#data generating function
data_gen_fun = function(q_values, obs_data) {
  y = prob_model_matrix(q_values, obs_data)
  return(y)
}

y1 = data_gen_fun(q_values, obs_data)
y2 = data_gen_fun(q_values, obs_data)
gen_data = list("y1" = y1, "y2" = y2)

#testing
y_test1 = data_gen_fun(runif(2), 5, df_same_time1)
y_test2 = data_gen_fun(runif(2), 5, df_same_time2)
test = list("y1" = y_test1, "y2" = y_test2)
test

#===========================================================================
#Summary Statistics Function

distance = function(obs_data1, obs_data2, gen_data) {
  obs_data1[is.na(obs_data1)] <- 0
  obs_data2[is.na(obs_data2)] <- 0
  return(0.5 * (norm(obs_data1 - gen_data$y1, type = "F") + 
                  norm(obs_data2 - gen_data$y2, type = "F")))
}


#testing
distance(df_same_time1, df_same_time2, test)

#===========================================================================

#posterior samples
posterior_samples = replicate(n = 10000, 
                           abc_sample_generate(obs_data1, 
                                               obs_data2, sum_stat, eps))

#===========================================================================

model1 <- replicate(n = 10000, 
                   abc_sample_generate(df_same_time1, df_same_time2,eps = 200))

model1 <- as.data.frame(t(model1))
rownames(model1) <- NULL
colnames(model1) <- c("qh1", "qc1", "qh2", "qc2")
model1

plot(model1[,1], model1[,2], col = c("red"), xlim = c(0,1), ylim = c(0,1), xlab = "qh", ylab = "qc",
     points(model1[,3], model1[,4], col = c("blue")))

model2 = replicate(n = 10, 
                   abc_sample_generate(df_diff_A, df_diff_B,eps = 50))
model2 <- as.data.frame(t(model2))   
rownames(model2) <- NULL
colnames(model2) <- c("qh1", "qc1", "qh2", "qc2")

plot(model2[,1], model2[,2], col = c("red"), xlim = c(0,1), ylim = c(0,1), xlab = "qh", ylab = "qc",
     points(model2[,3], model2[,4], col = c("blue")))
