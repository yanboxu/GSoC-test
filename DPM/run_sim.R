setwd("~/Dropbox/Yanbo/LDTR/experiments/C_sampling/DPM/")
source("sim.R")

dimI = 400

################
sim_data = sim_dpm(dimI)
################

K = sim_data$K

dimI = sim_data$dimI
dimP = sim_data$dimP

dataX = sim_data$dataX

paramZ = sim_data$Z
paramMu = sim_data$Mu
paramS = sim_data$Sigma

#######################
# plot data
plot(1, type = "n", xlab = "", ylab = "", 
     xlim = c(min(dataX[,1]), max(dataX[,1])),
     ylim = c(min(dataX[,2]), max(dataX[,2])))

for(i in 1:dimI){
  points(x = dataX[i,1], y = dataX[i,2], 
       col = paramZ[i], cex = .8, pch = 20,
       xlab = "", ylab = "")
}

for(k in 1:K){
  points(x = paramMu[1,k], y = paramMu[2,k], col = k, cex = 1.5, pch = 8, lwd = 2)
}
#######################

source("DPM.R")
tol_iter = 1000
thin = 10
thined_tol_iter = tol_iter/thin
model_fit = dpm_fit(tol_iter, thin, K, sim_data)

fitZ = model_fit$Z[,thined_tol_iter-1] + 1
fitMu = model_fit$Mu[,,thined_tol_iter-1]
fitS =  model_fit$Sigma[,,,thined_tol_iter-1]
#######################
# plot data
plot(1, type = "n", xlab = "", ylab = "", 
     xlim = c(min(dataX[,1]), max(dataX[,1])),
     ylim = c(min(dataX[,2]), max(dataX[,2])))

for(i in 1:dimI){
  points(x = dataX[i,1], y = dataX[i,2], 
         col = fitZ[i], cex = .8, pch = 20,
         xlab = "", ylab = "")
}

for(k in 1:K){
  points(x = fitMu[1,k], y = fitMu[2,k], col = k, cex = 1.5, pch = 8, lwd = 2)
}
#######################
