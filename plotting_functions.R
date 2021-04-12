### GPS plotting functions ###
minMiss_histplot<-function(data, bins, title="", tag ="", percent = T){
  total_days = dim(data)[1]
  if (percent == F) {
    p<-ggplot(data, aes(x=MinsMissing)) + 
      geom_histogram(fill="lightgrey",bins = bins) +
      geom_vline(data=data, aes(xintercept=1440, color = paste0("Missing all data (1440 min, n= ", length(which(MinsMissing == 1440))," days)")),linetype="dashed") +
      geom_vline(data=data, aes(xintercept=1296, color = paste0("Missing no data (1296 min, n= ", length(which(MinsMissing <= 1296))," days)")),linetype="dashed")
    
  } else {
    mean = round(mean(data$MinsMissing),0)
    percents = round(quantile(data$MinsMissing,prob = c(0.75, 0.8, 0.85, 0.9), na.rm=TRUE),0)
    num_percent = sapply(percents, FUN = function(percent) length(which(data$MinsMissing <= percent)))
    legends = c()
    for (i in 1:length(percents)){
      legends[i] = paste0(names(percents[i])," (",percents[i]," min, n= ",num_percent[i]," days)")
    }
    
    p<-ggplot(data, aes(x=MinsMissing)) + 
      geom_histogram(fill="lightgrey",bins = bins) +
      geom_vline(data=data, aes(xintercept=1440, color = paste0("Missing all data (1440 min, n= ", length(which(MinsMissing == 1440))," days)")),linetype="dashed") +
      geom_vline(data=data, aes(xintercept=1296, color = paste0("Missing no data (1296 min, n= ", length(which(MinsMissing <= 1296))," days)")),linetype="dashed") +
      geom_vline(data=data, aes(xintercept=mean, color = paste0("mean (", mean, " min)")),linetype="dashed") +
      geom_vline(data=data, aes(xintercept=as.numeric(percents[1]), color = legends[1]),linetype="dashed") +
      geom_vline(data=data, aes(xintercept=as.numeric(percents[2]), color = legends[2]),linetype="dashed") +
      geom_vline(data=data, aes(xintercept=as.numeric(percents[3]), color = legends[3]),linetype="dashed") +
      geom_vline(data=data, aes(xintercept=as.numeric(percents[4]), color = legends[4]),linetype="dashed") 
  }
  p +  ggtitle(title) + theme_cowplot() + labs(color = paste0("n=",total_days," days")) +
    xlab("Minutes Missing in a Day") + theme(plot.title = element_text(hjust = 0.5), legend.title.align = 0.5) + labs(tag = tag)
}


### cor plot enhanced ###
rquery.cormat <-function(x, type=c('lower', 'upper', 'full', 'flatten'),
         graph=TRUE, graphType=c("correlogram", "heatmap"),
         col=NULL, ...)
{
  library(corrplot)
  # Helper functions
  #+++++++++++++++++
  # Compute the matrix of correlation p-values
  cor.pmat <- function(x, ...) {
    mat <- as.matrix(x)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], na.action = na.omit,...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  # Get lower triangle of the matrix
  getLower.tri<-function(mat){
    upper<-mat
    upper[upper.tri(mat)]<-""
    mat<-as.data.frame(upper)
    mat
  }
  # Get upper triangle of the matrix
  getUpper.tri<-function(mat){
    lt<-mat
    lt[lower.tri(mat)]<-""
    mat<-as.data.frame(lt)
    mat
  }
  # Get flatten matrix
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  # Define color
  if (is.null(col)) {
    col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                              "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                              "#4393C3", "#2166AC", "#053061"))(200)
    col<-rev(col)
  }
  
  # Correlation matrix
  cormat<-signif(cor(x, use = "pairwise.complete.obs", ...),2)
  pmat<-signif(cor.pmat(x, ...),2)
  # Reorder correlation matrix
  #ord<-corrMatOrder(cormat, order="alphabet")
  #cormat<-cormat[ord, ord]
  #pmat<-pmat[ord, ord]
  # Replace correlation coeff by symbols
  sym<-symnum(cormat, abbr.colnames=FALSE)
  # Correlogram
  if(graph & graphType[1]=="correlogram"){
    corrplot(cormat, method = "color", type=ifelse(type[1]=="flatten", "lower", type[1]),
             tl.col="black", tl.srt=45,col=col,...)
  }
  else if(graphType[1]=="heatmap")
    heatmap(cormat, col=col, symm=TRUE)
  # Get lower/upper triangle
  if(type[1]=="lower"){
    cormat<-getLower.tri(cormat)
    pmat<-getLower.tri(pmat)
  }
  else if(type[1]=="upper"){
    cormat<-getUpper.tri(cormat)
    pmat<-getUpper.tri(pmat)
    sym=t(sym)
  }
  else if(type[1]=="flatten"){
    cormat<-flattenCorrMatrix(cormat, pmat)
    pmat=NULL
    sym=NULL
  }
  list(r=cormat, p=pmat, sym=sym)
}


### nice histogram ##
hist_chx = function(vector, bins=30, title="",xaxis="",yaxis="", fill="lightgrey"){
  df = data.frame(val = vector)
  ggplot(df, aes(x=val)) + 
    geom_histogram(bins = bins, fill = fill) + 
    theme_cowplot() +
    labs(title=title,
         x =xaxis, y = yaxis) + 
    theme(plot.title = element_text(hjust = 0.5))
}

### nice scatterplot ###
scatter_chx = function(vector){
  df = data.frame(val = vector)
  ggplot(df, aes(x=val)) + 
    gemo
}


### match functions ###
calc_match_vector = function(subj_mat_1,subj_mat_2,method) {
  part_times = length(subj_mat_1[[1]])
  database = list() 
  for (time in 1:part_times) {
      database[[time]] = lapply(subj_mat_2, function(subjmat) subjmat[[time]])
  }
  
    # match target to database
    match_cor = list()
    for (subj1 in subject_seq){ #loop through each subj
      # create a list of target across partitions
      target_list = lapply(subj_mat_1[[subj1]], function(part) part)
      # create a match list
      for (time in 1:part_times){
        target_subj_time = target_list[[time]] #loop through each partition
        for (subj2 in names(subj_mat_2)){ #loop everyone in 2nd half
          data_subj_time = subj_mat_2[[subj2]][[time]]
          if (method == "cor") {
          match_cor[[subj1]][[as.character(time)]][[subj2]] = cor(target_subj_time,data_subj_time,use = "na.or.complete")
          } 
          else if (method == "rmse") {
            match_cor[[subj1]][[as.character(time)]][[subj2]] = rmse(target_subj_time,data_subj_time)
          }
        }
      }
    }
  return(match_cor)
}


calc_acc_time=function(match_cor,method = "max") {
  acc_time = array()
  part_times = length(match_cor[[1]])
  for (time in 1:part_times){
    acc_time[time] = 0
    for (subj in subject_seq){
      if (method == "max"){
        position = which.max(unlist(match_cor[[subj]][[as.character(time)]]))
      } 
      else if (method == "min"){
        position = which.min(unlist(match_cor[[subj]][[as.character(time)]]))
      }
      
      predicted_subj = subject_seq[position]
      if (predicted_subj == subj) {
        acc_time[time] = acc_time[time] + 1
      }
    }
  }
  acc_time = acc_time/length(subject_seq)
  return(acc_time)
}

calc_acc_subj = function(match_cor, method = "max"){
  acc_subj = array()
  part_times = length(match_cor[[1]])
  for (subj in subject_seq){
    acc_subj[subj] = 0
    for (time in 1:part_times){
      if (method == "max"){
        position = which.max(unlist(match_cor[[subj]][[as.character(time)]]))
      } 
      else if (method == "min"){
        position = which.min(unlist(match_cor[[subj]][[as.character(time)]]))
      }
      predicted_subj = subject_seq[position]
      if (predicted_subj == subj) {
        acc_subj[subj] = acc_subj[subj] + 1
      }
    }
  }
  acc_subj = acc_subj/part_times
  acc_subj = acc_subj[-1]
  return(acc_subj)
}

make_feature_matrix = function(gps_df, subj_seq, range) {
  subj_mat_1 = list()
  subj_mat_2 = list()
  for (subj in names(subj_seq)){
    print(subj)
    subj_data = subset(gps_df, IID==subj)
    subj_mat_1[[subj]] = lapply(subj_seq[[subj]], function(list) rquery.cormat(subj_data[list,range], type = "flatten", graph = F)$r)
    subj_mat_2[[subj]] = lapply(subj_seq[[subj]], function(list) rquery.cormat(subj_data[-list,range], type = "flatten", graph = F)$r)
  }
  return(list(subj_mat_1 = subj_mat_1, subj_mat_2 = subj_mat_2 ))
}

make_subj_seq = function(gps_df_clean2,part_times) {
  subj_seq = list()
  for (subj in unique(gps_df_clean2$IID)){
    subj_data = subset(gps_df_clean2, IID==subj)
    if (dim(subj_data)[1] >5) {
    subj_seq[[subj]] <-createDataPartition(subj_data$IID,times = part_times, p =0.5)
    }
  }
  return(subj_seq)
}

combine_perm_subj = function(subject_seq, perm_acc_subj) {
  perm_subj = list()
  for (subj in subject_seq) {
    perm_subj$val[[subj]] = sapply(perm_acc_subj, function(perm) perm[which(names(perm) == subj)])
    perm_subj$hist[[subj]] = hist_chx(perm_subj[[subj]]$val, bins = 8, title = paste(subj,": Accuracy Across \n",perm_time,"Permutations"), xaxis = "Prediction Accuracy", yaxis = "Count")
    perm_subj$acc[[subj]] = sum(perm_subj$val[[subj]])/perm_time
  }
  return(perm_subj)
}


subj_scatter_perm = function(gps_df,acc_subj,perm_subj, method){
    subj_df <- data.frame(x=unique(gps_df$IID))
  subj_df$y = acc_subj
  subj_df$y_perm = unlist(perm_subj$acc)
  subj_df = subj_df[order(subj_df$y),]
  p = ggplot(subj_df, aes(x = reorder(x, y), y = value)) + 
    geom_point(aes(y = y, col = "subject data")) + 
    geom_point(aes(y = y_perm, col = "permutation")) +
    theme_cowplot() + 
    labs(title = paste(method,"Features \n Subject Level Accuracy Across",part_times,"Data Partitions"), 
         x = "Subjects", y = "Prediction Accuracy") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), plot.title = element_text(hjust = 0.5))
  return(p)
}

conf_scatter_plot = function(conf_var, acc_subj, col){
  if (!identical(conf_var$IID,names(acc_subj))) {
    stop("subj names do not match")
  }
  conf_df = data.frame(conf = unlist(conf_var[,col]), acc_subj = acc_subj)
  conf_scatter = ggscatter(conf_df, y = "acc_subj", x = "conf",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   ) + stat_cor(method = "pearson") +
    theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) +
    labs(title = paste(col), y = "Footprint Distinctiveness", x = col)
  return(conf_scatter)
  }

get_data_mean = function(subj){
  acc_data = readRDS(file.path("~/Documents/xia_gps/beiwe_output_043020/Results/Group/accelerometer/",subj,"/accelerometer_raw.rds"))
  mean(sapply(unique(acc_data$date), function(i) dim(subset(acc_data,date == i))[1]))
}

get_question_ans = function(question_of_interest) {
  question_idx = lapply(names(survey_ft), function(subj) which(grepl(question_of_interest, names(survey_ft[[subj]]), fixed = TRUE) == TRUE))
  names(question_idx) = names(survey_ft)
  
  question_answer = lapply(names(survey_ft), function(subj) {q_answer = rbindlist(survey_ft[[subj]][question_idx[[subj]]]);
                                                             q_answer = q_answer[order(q_answer$date,q_answer$hour)]})
  names(question_answer) = names(survey_ft)
  
  return(question_answer)
}

get_question_variance = function(question_of_interest) {
  question_idx = lapply(names(survey_ft), function(subj) which(grepl(question_of_interest, names(survey_ft[[subj]]), fixed = TRUE) == TRUE))
  names(question_idx) = names(survey_ft)
  
  question_answer = lapply(names(survey_ft), function(subj) {q_answer = rbindlist(survey_ft[[subj]][question_idx[[subj]]]);
                                                             q_answer = q_answer[order(q_answer$date,q_answer$hour)]})
  names(question_answer) = names(survey_ft)
  
  q_var = unlist(lapply(subj_list, function(subj) rmssd(question_answer[[subj]]$ans_num, na.rm = T)))
  return(q_var)
}

get_question_mean = function(question_of_interest) {
  question_idx = lapply(names(survey_ft), function(subj) which(grepl(question_of_interest, names(survey_ft[[subj]]), fixed = TRUE) == TRUE))
  names(question_idx) = names(survey_ft)
  
  question_answer = lapply(names(survey_ft), function(subj) {q_answer = rbindlist(survey_ft[[subj]][question_idx[[subj]]]);
                                                             q_answer = q_answer[order(q_answer$date,q_answer$hour)]})
  names(question_answer) = names(survey_ft)
  
  q_mean = unlist(lapply(subj_list, function(subj) mean(question_answer[[subj]]$ans_num, na.rm = T)))
  return(q_mean)
}

get_q_var_pvals = function(question_of_interest){
  print(question_of_interest)
  q_var_lm_df = acc_demo_df
  q_var_lm_df$mood_var = q_variable_list[,question_of_interest]
  if (substr(question_of_interest,0,5) == "rmssd") {
    mean_qoi = paste0("mean",substr(question_of_interest,6,nchar(question_of_interest)))
    q_var_lm_df$mood_mean = q_variable_list[,mean_qoi]
    q_var_accgps_fit = gam(accgps ~ mood_var + gps_days + acc_days + s(admin_age) + admin_sex + mood_mean, data = q_var_lm_df,  method = "REML")
  } else {
    q_var_accgps_fit = gam(accgps ~ mood_var + gps_days + acc_days + s(admin_age) + admin_sex , data = q_var_lm_df,  method = "REML")
  }
  q_var_accgps_p = summary(q_var_accgps_fit)[[4]]['mood_var']
  return(data.frame(p_accgps = q_var_accgps_p))
}

vector_to_mat = function(fc_vector, method = "half"){
  
  if (method == "full"){
    num_node = sqrt(length(fc_vector))
    matrix_fc = matrix(NA, num_node, num_node)
    matrix_fc = matrix(fc_vector, num_node, num_node)
  } else if (method == "half"){
    num_node = ceiling(sqrt(length(fc_vector)*2))
    matrix_fc = matrix(NA, num_node, num_node)
    matrix_fc[lower.tri(matrix_fc, diag = F)] = fc_vector
  } else if (method == "half+d") {
    num_node = floor(sqrt(length(fc_vector)*2))
    matrix_fc = matrix(NA, num_node, num_node)
    matrix_fc[lower.tri(matrix_fc, diag = T)] = fc_vector
  }
  
  matrix_fc = symmetrize(matrix_fc, rule = "lower")
  return(matrix_fc)
}

get_fc_com = function(atlas_now, subj_net) {
    com_nums = read.table(file.path(atlas_path,atlas_now,paste0(atlas_now,"CommunityAffiliation.1D")))$V1
    com_names = read.table(file.path(atlas_path,atlas_now,paste0(atlas_now,"CommunityNames.txt")))$V1
    com_names = substr(com_names,0,3)
    if (atlas_now == "power264") {com_names[c(1,2)] = c("smH","smM")}
    fc_mat_now_all = list()
    com_name_vec = c()
    for (net_1 in 1:max(com_nums)){
      node_now_1 = which(com_nums == net_1)
      net_name_1 = com_names[net_1]
      for (net_2 in 1:max(com_nums)){
        node_now_2 = which(com_nums == net_2)
        net_name_2 = com_names[net_2]
        network_now = paste(net_name_1,net_name_2,sep="_")
        com_name_vec = c(com_name_vec,network_now)
        for (subj in subjects) {
        fc_mat = vector_to_mat(subj_net[[atlas_now]][[subj]]$V1)
        fc_mat_now = fc_mat[node_now_1,node_now_2]
        fc_mat_now = fc_mat_now[lower.tri(fc_mat_now)]
        fc_mat_now_all[[network_now]][[subj]]$V1 = as.vector(fc_mat_now)
        }
      }
    }
    com_name_mat = vector_to_mat(com_name_vec, "full")
    com_name_vec_short = com_name_mat[lower.tri(com_name_mat, diag = T)]
    fc_mat_now_all_short = fc_mat_now_all[names(fc_mat_now_all) %in% com_name_vec_short]   
    return(fc_mat_now_all_short)
  }

  get_fc_nodal = function(atlas_now, subj_net) {
    com_nums = read.table(file.path(atlas_path,atlas_now,paste0(atlas_now,"CommunityAffiliation.1D")))$V1
    com_names = read.table(file.path(atlas_path,atlas_now,paste0(atlas_now,"CommunityNames.txt")))$V1
    com_names = substr(com_names,0,3)
  
  if (atlas_now == "power264") {com_names[c(1,2)] = c("smH","smM")}
  fc_node_now_all = list()

  for (subj in subjects) {
    fc_mat = vector_to_mat(subj_net[[atlas_now]][[subj]]$V1)
    fc_node = apply(fc_mat,2,function(x) sum(x, na.rm = T))
    names(fc_node) = paste(com_names[com_nums], 1:length(com_nums),sep = "_")
    fc_node_now_all[[subj]]$V1 = fc_node
  }
  
  return(fc_node_now_all)
}

get_fc_mat = function(atlas_now, subj_net) {
  fc_mat_all = list()
  for (subj in subjects) {
    fc_mat = vector_to_mat(subj_net[[atlas_now]][[subj]]$V1)
    fc_mat_all[[subj]]$fc_matrix = fc_mat
    fc_mat_all[[subj]]$covar = data.frame(motion = motion_df[subj],  subset(acc_subj_mood_df, BBLID == subj) )
  }
  return(fc_mat_all)
}

get_fc_edge = function(atlas_now, subj_net) {
    com_nums = read.table(file.path(atlas_path,atlas_now,paste0(atlas_now,"CommunityAffiliation.1D")))$V1
    com_names = read.table(file.path(atlas_path,atlas_now,paste0(atlas_now,"CommunityNames.txt")))$V1
    com_names = substr(com_names,0,3)
  
  if (atlas_now == "power264") {com_names[c(1,2)] = c("smH","smM")}
  fc_edge_all = list()
  
  mat_names = matrix(NA, nrow = length(com_nums), ncol = length(com_nums))
  for (node_i in 1:length(com_nums)){
    for (node_j in 1:length(com_nums)){
      mat_names[node_i,node_j] = paste(paste0(com_names[com_nums],1:length(com_nums),"")[node_i],
                                       paste0(com_names[com_nums],1:length(com_nums),"")[node_j],sep = "_")
    }
  }
  
  for (subj in subjects) {
    fc_mat = vector_to_mat(subj_net[[atlas_now]][[subj]]$V1)
    
    edge_names = mat_names[lower.tri(mat_names, diag = F)]
    fc_edge = subj_net[[atlas_now]][[subj]]$V1
    names(fc_edge) = edge_names
    
    fc_edge_all[[subj]] = fc_edge
  }
  return(fc_edge_all)
}

lasso_edge_perm = function(edge_data ,edge_vars, lambdas) {
  lasso_out = data.frame()
  for (i in 1:dim(edge_data)[1]) {
    print(i)
    df_train = edge_data[-i,]
    df_test = edge_data[i,]
    
    df_train_fit = gam(accgps  ~ s(admin_age) + admin_sex + gps_days + acc_days  + motion , data = df_train, method = "REML")
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

lasso_edge = function(edge_data ,edge_vars, lambdas) {
  lasso_out = data.frame()
  for (i in 1:dim(edge_data)[1]) {
    print(i)
    df_train = edge_data[-i,]
    df_test = edge_data[i,]
    
    df_train_fit = gam(accgps  ~ s(admin_age) + admin_sex + gps_days + acc_days  + motion , data = df_train, method = "REML")
    df_train$accgps_res = df_train_fit$residuals
    
    x_train = as.matrix(df_train[,edge_vars])
    y_train = df_train$accgps_res
    
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

get_final_lass = function(edge_data,edge_vars, sub_seq){
  edge_data = edge_data[sub_seq,]
  edge_data$accgps_res = gam(accgps  ~ s(admin_age) + admin_sex + gps_days + acc_days  + motion , data = edge_data, method = "REML")$residuals
  x_all = as.matrix(edge_data[,edge_vars])
  y_all = edge_data$accgps_res
  lasso_reg = glmnet(x_all, y_all, lambda = 0.07, alpha = 1, standardize = TRUE)
}

calc_match_cor = function(subj_mat_1,subj_mat_2) {
  part_times = length(subj_mat_1[[1]])
  database = list() 
  for (time in 1:part_times) {
      database[[time]] = lapply(subj_mat_2, function(subjmat) subjmat[[time]]$cor)
  }
  
    # match target to database
    match_cor = list()
    for (subj1 in names(subj_mat_1)){ #loop through each subj
      # create a list of target across partitions
      target_list = lapply(subj_mat_1[[subj1]], function(part) part$cor)
      # create a match list
      for (time in 1:part_times){
        target_subj_time = target_list[[time]] #loop through each partition
        for (subj2 in names(subj_mat_2)){ #loop everyone in 2nd half
          data_subj_time = subj_mat_2[[subj2]][[time]]$cor
          match_cor[[subj1]][[as.character(time)]][[subj2]] = cor(target_subj_time,data_subj_time,use = "na.or.complete")
        }
      }
    }
  return(match_cor)
}