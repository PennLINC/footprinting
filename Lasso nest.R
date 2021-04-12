require("mgcv")
require("glmnet")
setwd("/storage/xia_mobile/beiwe_output_043020/Results/Group")
#load("./fc_accgps_lasso_nested.RData")
edge_vars = readRDS("./edge_vars.RDS")
edge_data = readRDS("./fpt_fc_edge.RDS")
lasso_edge_perm = function(edge_data ,edge_vars, lambdas) {
  lasso_out = data.frame()
  for (i in 1:dim(edge_data)[1]) {
    print(i)
    df_train = edge_data[-i,]
    df_test = edge_data[i,]
    
    df_train_fit = gam(accgps  ~ s(age) + sex + gps_days + acc_days  + motion , data = df_train, method = "REML")
    df_train$accgps_res = df_train_fit$residuals
    
    x_train = as.matrix(df_train[,edge_vars])
    y_train = df_train$accgps_res
    y_train = sample(y_train)
    
    x_test = as.matrix(df_test[,edge_vars])
    y_test = df_test$accgps - predict(df_train_fit,newdata = df_test)
    
    lasso_reg = glmnet(x_train, y_train, lambda = lambdas, alpha = 1, standardize = TRUE)
    predictions_test <- predict(lasso_reg, newx = x_test)
    
    lasso_out[i,"real"] = y_test
    lasso_out[i,colnames(predictions_test)] = as.numeric(predictions_test)
  }
  lasso_cor = apply(lasso_out,2,function(col) cor(col,lasso_out$real))[-1]
  return(lasso_cor)
}


lasso_perm = list()
for (j in 1:100){
  lambdas = rev(seq(0.01,0.10,0.01))
  lasso_tune_l = list()
  lasso_out = data.frame()
  for (i in 1:dim(edge_data)){
    df_train = edge_data[-i,]
    df_test = edge_data[i,]
    
    lasso_tune = lasso_edge_perm(df_train,edge_vars,lambdas)
    lasso_tune_l[[i]] = lasso_tune
    best_lambda = lambdas[which.max(lasso_tune)]
    
    print(lasso_tune)
    print(paste(i,"-----", "best lambda is", best_lambda))
    
    df_train_fit = gam(accgps  ~ s(age) + sex + gps_days + acc_days  + motion , data = df_train, method = "REML")
    df_train$accgps_res = df_train_fit$residuals
    
    x_train = as.matrix(df_train[,edge_vars])
    y_train = df_train$accgps_res
    y_train = sample(y_train)
    
    x_test = as.matrix(df_test[,edge_vars])
    y_test = df_test$accgps - predict(df_train_fit,newdata = df_test)
    
    lasso_reg = glmnet(x_train, y_train, lambda = best_lambda, alpha = 1, standardize = TRUE)
    
    predictions_test <- predict(lasso_reg, s = best_lambda, newx = x_test)
    
    lasso_out[i,"predict"] = as.numeric(predictions_test)
    lasso_out[i,"real"] = y_test
  }
  lasso_perm[[j]] = cor(lasso_out$predict,lasso_out$real)
  saveRDS(lasso_perm,file = "./lasso_perm_1.RDS")
}

screen -S nest10
R

setwd("/storage/xia_mobile/beiwe_output_043020/Results/Group")
perm_all = c()
for (i in 1:10){
  perm_all = c(perm_all,unlist(readRDS(paste0("./lasso_perm_",i,".RDS"))))
}
length(which(perm_all>0.2887))/length(perm_all)
length(perm_all)

