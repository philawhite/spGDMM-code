library(splines)
library(fields)
library(geosphere)
library(splines2)
library(nimble)
library(vegan)
library(TruncatedNormal)
library(emulator)
############## map for families
# setwd("C:/Users/philaw/Box/Research/Turnover/code")
#setwd("C:/Users/philaw/Box/Research/ERT/data")
rm(list = ls())

# K = 10
# 
# mod_out = mclapply(1:K,function(qqq){

dat_all = read.csv("../../data/sa_family_data.csv")

location_mat = dat_all[,1:3] 
envr_use = dat_all[,4:12] 
species_mat = dat_all[,-(1:12)] 
ns = nrow(location_mat)

# set.seed(1)
# tmp = sample(1:ns)
# ind_loc_hold = lapply(1:K,function(i){
#   idx_keep = which((1:ns)%% K == (i-1))
#   sort(tmp[idx_keep])
# })

dist_use = as.matrix(vegdist(species_mat,"bray"))
mean(dist_use[upper.tri(dist_use)] == 0)
mean(dist_use[upper.tri(dist_use)] == 1)


# cord.dec = SpatialPoints(cbind(location_mat$longitude, location_mat$latitude), 
#                          proj4string = CRS("+proj=longlat"))
# 
# Transforming coordinate to UTM using EPSG=32748 for WGS=84, UTM Zone=48M,
# Southern Hemisphere)
# cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=34S +datum=WGS84"))
# location_mat$easting = cord.UTM@coords[,1] - cord.UTM@bbox[1,1]
# location_mat$northing = cord.UTM@coords[,2] - cord.UTM@bbox[2,1]


dist_mat= distm(cbind(location_mat$longitude,location_mat$latitude))/1e3


X = envr_use[,c("gmap","RFL_CONC","Elevation30m","HeatLoadIndex30m","tmean13c",
                "SoilConductivitymSm","SoilTotalNPercent")]

# cor_melt = melt(cor(X[,-c(3,4)]))


# ggplot(cor_melt) + geom_tile(aes(x = Var1,y = Var2,fill = value))+
#   scale_fill_viridis()

deg = 3
knots = 2
df_use = deg + knots

formula_use = as.formula(paste("~ 0 +",paste(
  paste("iSpline(",colnames(X),",degree=",deg - 1 ,",df = ",df_use, 
        " ,intercept = TRUE)"),collapse = "+")))

y_GDM = dist_use[upper.tri(dist_use)]


I_spline_bases = model.matrix(formula_use,data = X)

vec_distance = dist_mat[upper.tri(dist_mat)]


X_GDM = cbind(sapply(1:ncol(I_spline_bases),function(i){
  
  dist_temp = rdist(I_spline_bases[,i])
  vec_dist = dist_temp[upper.tri(dist_temp)]
  vec_dist
  
}),
iSpline(vec_distance,degree = deg -1,
        df = df_use,intercept = TRUE)
)


colnames(X_GDM) = c(
  paste(rep(colnames(X),each = df_use ),"I",rep(1:df_use,times = ncol(X)),sep = ""),
  paste("dist","I",1:df_use,sep = "")
)
# dat = as.data.frame(X_GDM)
y = y_GDM
# dat$Dissimilarity = y_GDM
# dat$log_dissimilarity = log(dat$y)
dist = vec_distance


tmp = matrix(rep(1:nrow(dist_use),each = nrow(dist_use)),nrow = nrow(dist_use))
col_ind = tmp[upper.tri(tmp)]
tmp = matrix(rep(1:nrow(dist_use),times = nrow(dist_use)),nrow = nrow(dist_use))
row_ind = tmp[upper.tri(tmp)]

# attach(dat)
# idx_hold = sort(which(col_ind %in% ind_loc_hold[[qqq]] | row_ind %in% ind_loc_hold[[qqq]]))
# y_hold = y[idx_hold]
# y[idx_hold] = NA
# n_miss = length(idx_hold)
Y_is_one = which(y == 1)
# Y_is_missing = idx_hold
Y_is_not_one = which(y != 1)

n1 = length(Y_is_one)
N = length(y)
p = ncol(X_GDM)

# lm_mod= lm(as.formula(paste("log(y)~",paste(colnames(X_GDM),collapse = "+",sep = ""),sep = "")),
#            data = dat)
lm_out = optim(c(.3, rnorm(p),rnorm(ns)) ,function(par){
  sum((log(y)  - par[1] - X_GDM %*% exp(par[2:(p + 1)]) -
         abs(par[p+1+row_ind] - par[p+1 + col_ind]) )^2)
},method = "BFGS")

# saveRDS(lm_out,"init_par.rds")
# lm_out = readRDS("init_par.rds")
# log_beta_init = ifelse(coef(lm_mod)[-1]> 0,log(coef(lm_mod)[-1]), -10)
# alpha_init = coef(lm_mod)[1]
# eta_init = rep(0,ns)
# summary(lm_mod$p)

reps = 10e3
thin = 10
burn = 100e3
tune = 100
# reps = 1e2
# thin = 1
# burn = 2e2
# tune = 100
# cov_tune = 1e3
# cov_tune_start = cov_tune + 1e3
### tuning stuff

cand_sd_eta = .01
count_eta = 0

count_beta = rep(0,p+1)
cand_sd_beta = rep(.01,p+1)

# count_alpha = 0
# cand_sd_alpha = .01

# X_sigma = cbind(1,X_GDM[,grep("dist",colnames(X_GDM))])
X_sigma = cbind(1,poly(vec_distance,degree = 3))

p_sigma = ncol(X_sigma)

count_sigma = rep(0,p_sigma+1)
cand_sd_sigma = rep(.01,p_sigma+1)

#### parameter

rho_fix = max(dist_mat)/10
R_spat = exp(-dist_mat/rho_fix)
chol_R = t(chol(R_spat))
R_inv = solve(R_spat)

sig2_eta_save = numeric(reps) ; sig2_eta_now = 1
log_beta_sigma = matrix(NA,reps,p_sigma); log_beta_sigma_now =  rep(0,p_sigma)
sig2_V_now = exp(c(X_sigma %*% log_beta_sigma_now))
log_v_now = numeric(N)
eta_save = matrix(NA,reps, ns) ; eta_now = lm_out$par[p + 1 + 1:ns]
log_beta_save = matrix(NA,reps,p); log_beta_now = c(lm_out$par[2:(p+1)])
alpha_save = numeric(reps); alpha_now = lm_out$par[1]
# y_pred_save = matrix(0,reps,n_miss); y_pCred_now = y_hold
log_lik = numeric(reps)

# log_beta_tune = matrix(NA,burn - cov_tune,p+1) 
# log_warp_tune = array(0,c(burn - cov_tune,n_w-1))

xb = c(alpha_now + X_GDM %*% exp(log_beta_now) + abs(eta_now[row_ind]  - eta_now[col_ind]))
# var(log(dat$y) - xb) /var(log(dat$y))

log_v_now[Y_is_not_one] = log(y[Y_is_not_one])
log_v_now[Y_is_one] = rtnorm(n1,mu = xb[Y_is_one],sd = sqrt(sig2_V_now[Y_is_one]),lb=0,ub= Inf)
# log_v_now[Y_is_missing] = rnorm(n_miss,mean = xb[Y_is_missing],sd = sqrt(sig2_V_now[Y_is_missing]))

like_now = sum(dnorm(log_v_now,xb,sqrt(sig2_V_now),log = TRUE))

st = proc.time()

for(i in 1:(reps*thin + burn)){
  
  ### Update folded normal parameters
  
  eta_cand = c(sqrt(1 - cand_sd_eta^2) * eta_now +
                 cand_sd_eta * sqrt(sig2_eta_now) * chol_R %*% rnorm(ns))
  
  xb_cand = c(alpha_now + X_GDM %*% exp(log_beta_now) +
                abs(eta_cand[row_ind]  - eta_cand[col_ind]))
  
  like_cand = sum(dnorm(log_v_now,xb_cand,sqrt(sig2_V_now),log = TRUE))
  
  r = like_cand - like_now
  if(!is.na(r)){
    if(r  > log(runif(1))){
      xb = xb_cand
      eta_now = eta_cand
      like_now = like_cand
      # prior_now = prior_cand
      count_eta = count_eta + 1
    }
  }
  ### Update variance of eta
  
  a_eta = 1 + ns/2
  b_eta = 1 + quad.form(R_inv,eta_now)/2 
  sig2_eta_now = 1/rgamma(1,a_eta,b_eta)
  
  ### Update log_beta and alpha parameters
  lbeta_now = c(alpha_now, log_beta_now)
  lbeta_cand = lbeta_now
  
  for(j in 1:(p+1)){
    
    lbeta_cand[j] =rnorm(1,lbeta_now[j],cand_sd_beta[j])
    
    xb_cand = c(lbeta_cand[1] + X_GDM %*% exp(lbeta_cand[-1]) +
                  abs(eta_now[row_ind] - eta_now[col_ind]))
    
    like_cand = sum(dnorm(log_v_now,xb_cand,sqrt(sig2_V_now),log = TRUE))
    
    r = like_cand - like_now+ 
      dnorm(lbeta_cand[j],0,sqrt(10),log = TRUE) - 
      dnorm(lbeta_now[j],0,sqrt(10),log = TRUE) 
    if(!is.na(r)){
      if(r  > log(runif(1))){
        xb = xb_cand
        lbeta_now = lbeta_cand
        
        like_now = like_cand
        # prior_now = prior_cand
        count_beta[j] = count_beta[j] + 1
      }
    }
  }
  
  alpha_now = lbeta_now[1]
  log_beta_now = lbeta_now[-1]
  
  
  # lbeta_cand = c(c(alpha_now,log_beta_now) + chol_lbeta %*% rnorm(p + 1))
  # xb_cand = c(lbeta_cand[1] + X_GDM %*% exp(lbeta_cand[-1]) +
  #               abs(eta_now[row_ind]  - eta_now[col_ind]))
  # 
  # like_cand = sum(dnorm(log_v_now,xb_cand,sqrt(sig2_V_now),log = TRUE))
  # 
  # r = like_cand - like_now
  # 
  # if(r  > log(runif(1))){
  #   xb = xb_cand
  #   alpha_now = lbeta_cand[1]
  #   log_beta_now = lbeta_cand[-1]
  #   
  #   like_now = like_cand
  #   # prior_now = prior_cand
  #   count_beta = count_beta + 1
  # }
  
  ### Update sigV
  
  log_beta_sigma_cand = log_beta_sigma_now
  
  for(j in 1:p_sigma){
    
    log_beta_sigma_cand[j] =rnorm(1,log_beta_sigma_now[j],cand_sd_sigma[j])
    
    sig2_V_cand = exp(c(X_sigma %*% log_beta_sigma_cand))
    
    like_cand = sum(dnorm(log_v_now,xb,sqrt(sig2_V_cand),log = TRUE))
    
    r = like_cand - like_now + 
      dnorm(log_beta_sigma_cand[j],0,sqrt(10),log = TRUE) - 
      dnorm(log_beta_sigma_now[j],0,sqrt(10),log = TRUE) 
    if(!is.na(r)){
      if(r  > log(runif(1))){
        sig2_V_now = sig2_V_cand
        log_beta_sigma_now = log_beta_sigma_cand
        
        like_now = like_cand
        # prior_now = prior_cand
        count_sigma[j] = count_sigma[j] + 1
      }
    }
  }
  
  
  ### Update log of V
  
  # log_v_now[Y_is_one] = rtruncnorm(1,a=0,b= Inf,mean = xb[Y_is_one],sd = sqrt(sig2_V_now))
  log_v_now[Y_is_one] = rtnorm(n1,mu = xb[Y_is_one],sd = sqrt(sig2_V_now[Y_is_one]),lb=0,ub= Inf)
  # log_v_now[Y_is_missing] = rnorm(n_miss,mean = xb[Y_is_missing],sd = sqrt(sig2_V_now[Y_is_missing]))
  ####### Update Likelihood
  
  like_now = sum(dnorm(log_v_now,xb,sqrt(sig2_V_now),log = TRUE))
  
  ####### Predict Y
  
  # exp_y_pred = exp(log_v_now[Y_is_missing])
  # 
  # y_pred_now = ifelse(exp_y_pred > 1, 1,exp_y_pred)
  
  
  if(i %% 100 == 0){
    time_its <- (proc.time() - st)[3] / (i)
    time_used <- round((proc.time() - st)[3]/(60),digits=4)
    time_left <- round(time_its * (reps*thin + burn - i )/(60),digits=4)
    cat("\r", i, " of ", reps*thin + burn," |||| like = ",like_now,"|| count_beta",
        count_beta,"|| count_eta = ",count_eta,"||| Time left: ",
        round(floor(time_left/60),1), " hours",time_left%%60," minutes")
    flush.console()
    
  }
  
  
  if(i < burn){
    
    
    if(i %% tune == 0){
      
      acc_beta = count_beta / tune; count_beta = rep(0,p+1)
      # acc_beta = count_beta / tune; count_beta = 0
      
      acc_eta = count_eta / tune; count_eta = 0
      
      acc_sigma = count_sigma/ tune; count_sigma = rep(0,p_sigma+1)
      
      
      cand_sd_sigma = ifelse(acc_sigma > 0.6, cand_sd_sigma*2,
                             ifelse(acc_sigma < 0.2, cand_sd_sigma/1.5, cand_sd_sigma) )
      
      cand_sd_beta = ifelse(acc_beta > 0.6, cand_sd_beta*2,
                            ifelse(acc_beta < 0.2, cand_sd_beta/1.5, cand_sd_beta) )
      
      cand_sd_eta = ifelse(acc_eta > 0.45, cand_sd_eta*2,
                           ifelse(acc_eta < 0.15, cand_sd_eta/1.5, cand_sd_eta) )
      
      
      # if(i > cov_tune_start){
      #   
      # 
      #   scale_par_beta = ifelse(acc_beta > 0.45, scale_par_beta*2,
      #                           ifelse(acc_beta < 0.15, scale_par_beta/1.5, scale_par_beta) )
      #   
      #   cov_alp_bet = cov(unique(log_beta_tune[1:(i-cov_tune),])) + 1e-12 * diag(p+1)
      #   cand_bet = scale_par_beta^2 * cov_alp_bet/(p+1)
      #   chol_lbeta = t(chol(cov_alp_bet))
      #   
      # }
      
      
    } 
    
    
    
    
  }
  
  if(i > burn & i %% thin == 0){
    sig2_eta_save[(i- burn)/thin] = sig2_eta_now
    log_beta_sigma[(i- burn)/thin,] = log_beta_sigma_now
    eta_save[(i- burn)/thin,] = eta_now
    log_beta_save[(i- burn)/thin,] = log_beta_now
    alpha_save[(i- burn)/thin] = alpha_now
    # y_pred_save[(i- burn)/thin,] = y_pred_now
    log_lik[(i- burn)/thin] = like_now
    
  }
  
}

# out = apply(reduce(mod_out,rbind),2,mean)
rm(list=ls())
