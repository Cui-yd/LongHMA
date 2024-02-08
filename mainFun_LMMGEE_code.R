longFun = function(Y, X, M, COV = NULL, id, wave, topN = NULL,
                   method = c("lme", "gee"), family, verbose = F){
  ## longitudinal method, based on linear mixed model or generalized estimating equation
  # input:
  # y:      longitudinal outcome, with id and waves
  # X:      longitudinal form of exposure
  # M:      longitudinal form of mediators
  # COV:    covariates
  # topN:   number of top mediators in SIS
  
  
  ##### FIRST STEP: INTERSECT SIS #####
  
  if(verbose) message("Step 1: Screening...", "     (", Sys.time(), ")")
  
  if(is.null(topN)) topN = 2*ceiling(n/log(n))
  
  SIS.beta = longBetaMarginal(Y = Y, M = M, X = X, COV = COV, id = id, wave = wave, p = p, method = method, family = family)
  SIS.alpha = longAlphaMarginal(X = X, M = M, COV = COV, p = p, wave = wave)
  SIS.beta.rank = sort(abs(SIS.beta[1, ]), decreasing = T)
  SIS.alpha.rank = sort(abs(SIS.alpha[1, ]), decreasing = T)
  top.beta = names(head(SIS.beta.rank, topN))
  top.alpha = names(head(SIS.alpha.rank, topN))
  subID_SIS = intersect(top.beta, top.alpha)
  
  if (length(subID_SIS) == 0) {
    s = 0
    while (s >= 0) {
      s = s + 1
      topN = s * ceiling(n/log(n))
      top.beta = names(head(SIS.beta.rank, topN))
      top.alpha = names(head(SIS.alpha.rank, topN))
      subID_SIS = intersect(top.beta, top.alpha)
      if (length(subID_SIS) > 0) {
        break
      }
    }
  }
  
  subM_SIS = data.frame(M[, subID_SIS])
  colnames(subM_SIS) = subID_SIS
  subdat_SIS = data.frame(id, Y, X, COV, subM_SIS)
  
  if(verbose) message("        Top ", length(subID_SIS), " mediators selected (by marginal correlation): ", paste(subID_SIS, sep = " "), "\n")
  
  
  
  ##### SECOND STEP: INDIRECT EFFECT TEST #####
  
  if(verbose) message("Step 2: Indirect effect testing ...", "     (", Sys.time(), ")")
  
  # beta estimation
  switch (method,
          "lme" = {
            lme.fit = suppressMessages(lmer(Y ~ . -id + (1+Z1|id), data = subdat_SIS))
            lme.coef = as.data.frame(summary(lme.fit)$coef)
            lme.M = lme.coef[grep("M", rownames(lme.coef), value=TRUE),]
            beta_coef = lme.M[,1]
            beta_var = lme.M[,2]^2
            beta_pval = lme.M[,5]
          },
          "gee" = {
            gee.fit = geeglm(Y ~ . -id, id=id, data=subdat_SIS, corstr="exchangeable", family = family)
            gee.coef = summary(gee.fit)$coef
            gee.M = gee.coef[grep("M", rownames(gee.coef), value=TRUE),]
            beta_coef = gee.M[, 1]
            beta_var = gee.M[, 2]^2
            beta_pval = gee.M[, 4]
          }
  )
  
  
  # alpha estimation
  alpha.fit = longAlphaMarginal(X = X, M = subM_SIS, COV = COV, 
                                p = ncol(subM_SIS), wave = wave)
  alpha_coef = alpha.fit[1, subID_SIS]
  alpha_var = alpha.fit[2, subID_SIS]
  alpha_pval = alpha.fit[3, subID_SIS]
  
  
  # p value calculation
  ab_coef = alpha_coef * beta_coef
  ab_var = (alpha_coef^2) * (beta_var) + (beta_coef^2) * (alpha_var)
  # confidence interval
  conf_low = ab_coef - qnorm(0.975)*sqrt(ab_var);  conf_up = ab_coef + qnorm(0.975)*sqrt(ab_var)
  
  
  # sobel test for alpha and beta
  s.test = (abs(ab_coef))/(sqrt(ab_var))                             # z-score of sobel test
  sob_pval = 2 * (1 - pnorm(s.test))                                     # p-value of sobel test
  sob_pval_fdr = p.adjust(sob_pval, "fdr", length(subID_SIS))
  sob_pval_bon = p.adjust(sob_pval, "bonferroni", length(subID_SIS))
  sob_pval_by = p.adjust(sob_pval, "BY", length(subID_SIS))
  
  if(verbose){
    message("        Significant FDR(BH) Sobel mediator(s): ", names(sob_pval_fdr)[which(sob_pval_fdr < 0.05)])
    message("        Significant bonforroni Sobel mediator(s): ", names(sob_pval_bon)[which(sob_pval_bon < 0.05)])
    message("        Significant FDR(BY) Sobel mediator(s): ", names(sob_pval_by)[which(sob_pval_by < 0.05)], "\n")
  }
  
  # joint test for alpha and beta
  alpha_pval_fdr = p.adjust(alpha_pval, "fdr", length(subID_SIS))
  beta_pval_fdr = p.adjust(beta_pval, "fdr", length(subID_SIS))
  pval_bind_fdr = rbind(alpha_pval_fdr, beta_pval_fdr)
  join_pval_fdr = apply(pval_bind_fdr, 2, max)
  
  alpha_pval_bon = p.adjust(alpha_pval, "bonferroni", length(subID_SIS))
  beta_pval_bon = p.adjust(beta_pval, "bonferroni", length(subID_SIS))
  pval_bind_bon = rbind(alpha_pval_bon, beta_pval_bon)
  join_pval_bon = apply(pval_bind_bon, 2, max)
  
  alpha_pval_by = p.adjust(alpha_pval, "BY", length(subID_SIS))
  beta_pval_by = p.adjust(beta_pval, "BY", length(subID_SIS))
  pval_bind_by = rbind(alpha_pval_by, beta_pval_by)
  join_pval_by = apply(pval_bind_by, 2, max)
  
  
  if(verbose){
    message("        Significant FDR(BH) Joint mediator(s): ", names(join_pval_fdr)[which(join_pval_fdr < 0.05)])
    message("        Significant bonforroni Joint mediator(s): ", names(join_pval_bon)[which(join_pval_bon < 0.05)])
    message("        Significant FDR(BY) Joint mediator(s): ", names(join_pval_by)[which(join_pval_by < 0.05)])
  }
  
  result = data.frame(sob_pval_fdr, sob_pval_bon, sob_pval_by, 
                      join_pval_fdr, join_pval_bon, join_pval_by,
                      ab_coef, ab_var, conf_low, conf_up, sob_pval, 
                      alpha_coef, alpha_var, alpha_pval, 
                      beta_coef, beta_var, beta_pval)
  
  return(list(result = result, subID_SIS = subID_SIS, subM_SIS = subM_SIS, 
              SIS.beta = SIS.beta, SIS.alpha = SIS.alpha,
              id = id, wave = wave, X = X, Y = Y, COV = COV))
}



longAlphaMarginal = function(X, M, COV, p, wave){
  ### linear model for X => M
  # input:
  # X: exposure
  # M: mediators
  # COV: covariates
  # p: dimension of M
  # wave: visit wave of all sample
  
  Est_alpha = matrix(0, 3, p)
  colnames(Est_alpha) = colnames(M)
  X = X[which(wave == min(wave))]
  M = data.frame(M[which(wave == min(wave)),])
  COV = COV[which(wave == min(wave)),]
  for(j in 1:p){
    MX = data.frame(M = M[,j], X = X, COV = COV)
    fit = lm(M ~ ., data = MX)
    Est_alpha[1,j] = summary(fit)$coef[2,1]          # coefficients of alpha
    Est_alpha[2,j] = summary(fit)$coef[2,2]^2        # variance of alpha
    Est_alpha[3,j] = summary(fit)$coef[2,4]          # p-value of alpha
  }
  return(Est_alpha=Est_alpha)
}



longBetaMarginal = function(Y, M, X, COV, id, wave, p, method = c("lme", "gee"), family){
  ### longitudinal model for M => Y
  # input:
  # Y: outcome
  # M: mediators
  # COV: covariates
  # id: id of all sample (longitudinal)
  # wave: visit wave of all sample
  # p: dimension of M
  # method: estimate method
  
  Est_beta = matrix(0, 3, p)
  colnames(Est_beta) = colnames(M)
  for(j in 1:p){
    MY = data.frame(M = M[, j], Y = Y, X = X, COV = COV, id = id, wave = wave)
    switch (method,
            "lme" = {
              lme.fit = suppressMessages(lmer(Y ~ . -id + (1+COV.Z1|id), data = MY))
              lme.coef = summary(lme.fit)$coef
              Est_beta[1,j] = lme.coef[2,1]                   # coefficients of beta
              Est_beta[2,j] = lme.coef[2,2]^2                 # variance of beta
              Est_beta[3,j] = lme.coef[2,5]                   # p-value of beta
            },
            "gee" = {
              gee.fit = geeglm(formula = Y ~ . -id, family = family, data = MY, id = id, corstr="exchangeable")
              gee.coef = summary(gee.fit)$coef
              Est_beta[1,j] = gee.coef[2,1]                   # coefficients of beta
              Est_beta[2,j] = gee.coef[2,2]^2                 # variance of beta
              Est_beta[3,j] = gee.coef[2,4]                   # p-value of beta
            }
    )             
  }
  return(Est_beta=Est_beta)
}






lineFun = function(Y, X, M, COV = NULL, id, wave, topN = NULL, verbose = F){
  ## the linear method, ignore correlation of longitudinal data
  # input:
  # y:      longitudinal outcome, with id and times
  # X:      longitudinal form of exposure
  # M:      longitudinal form of mediators
  # COV:    covariates
  # topN:   number of top mediators in SIS
  
  
  
  ##### FIRST STEP: SIS #####
  
  if(verbose) message("Step 1: Screening...", "     (", Sys.time(), ")")
  
  if(is.null(topN)) topN = 2*ceiling(n/log(n))
  
  SIS.beta = lineBetaMarginal(Y = Y, M = M, COV = COV, p = p)
  SIS.alpha = lineAlphaMarginal(X = X, M = M, COV = COV, p = p, wave = wave)
  SIS.beta.rank = sort(abs(SIS.beta[1, ]), decreasing = T)
  SIS.alpha.rank = sort(abs(SIS.alpha[1, ]), decreasing = T)
  top.beta = names(head(SIS.beta.rank, topN))
  top.alpha = names(head(SIS.alpha.rank, topN))
  subID_SIS = intersect(top.beta, top.alpha)
  
  if (length(subID_SIS) == 0) {
    s = 0
    while (s >= 0) {
      s = s + 1
      topN = s * ceiling(n/log(n))
      top.beta = names(head(SIS.beta.rank, topN))
      top.alpha = names(head(SIS.alpha.rank, topN))
      subID_SIS = intersect(top.beta, top.alpha)
      if (length(subID_SIS) > 0) {
        break
      }
    }
  }
  
  subM_SIS = data.frame(M[, subID_SIS])
  colnames(subM_SIS) = subID_SIS
  subdat_SIS = data.frame(id, Y, X, COV, subM_SIS, wave)
  
  
  if(verbose) message("        Top ", length(subID_SIS), 
                      " mediators selected (ranked by marginal p-value): ", subID_SIS, "\n")  
  
  
  ##### SECOND STEP: SIGNIFICANT TEST #####
  
  if(verbose) message("Step 2: Joint significance test ...", "     (", Sys.time(), ")")
  
  # beta estimate for comparison
  lm_form = "Y ~ ."
  lm.fit = lm(formula = lm_form, data = subdat_SIS)
  lm.coef = as.data.frame(summary(lm.fit)$coef)
  lm.M = lm.coef[grep("M", rownames(lm.coef), value=TRUE),]
  beta_coef = lm.M[,1]
  beta_var = (lm.M[,2])^2
  beta_pval = lm.M[,4]
  
  # alpha estimated from linear regression
  alpha.fit = lineAlphaMarginal(X = X, M = subM_SIS, COV = COV, p = length(subID_SIS), wave = wave)
  alpha_coef = alpha.fit[1,subID_SIS]
  alpha_var = alpha.fit[2,subID_SIS]
  alpha_pval = alpha.fit[3,subID_SIS]
  
  # p value calculation
  ab_coef = alpha_coef * beta_coef                               # estimation of alpha * beta
  ab_var = (alpha_coef^2) * (beta_var) + (beta_coef^2) * (alpha_var)
  # confidence interval
  conf_low = ab_coef - 1.96 * sqrt(ab_var);  conf_up = ab_coef + 1.96 * sqrt(ab_var)
  
  # sobel test for alpha and beta
  s.test = (abs(ab_coef))/(sqrt(ab_var))                            # z-score for sobel test
  sob_pval = 2 * (1-pnorm(s.test))                                   # p-value of sobel test
  sob_pval_fdr = p.adjust(sob_pval, 'fdr', length(subID_SIS))      
  sob_pval_fdr[sob_pval_fdr > 1] = 1
  sob_pval_bon = p.adjust(sob_pval, 'bonferroni', length(subID_SIS))           
  sob_pval_bon[sob_pval_bon > 1] = 1
  sob_pval_by = p.adjust(sob_pval, 'BY', length(subID_SIS))      
  sob_pval_by[sob_pval_by > 1] = 1
  
  
  if(verbose){
    message("        Significant FDR(BH) Sobel mediator(s): ", names(sob_pval_fdr)[which(sob_pval_fdr < 0.05)])
    message("        Significant bonforroni Sobel mediator(s): ", names(sob_pval_bon)[which(sob_pval_bon < 0.05)])
    message("        Significant FDR(BY) Sobel mediator(s): ", names(sob_pval_by)[which(sob_pval_by < 0.05)], "\n")
  }
  
  
  # joint test for alpha and beta
  alpha_pval_fdr = p.adjust(alpha_pval, 'fdr', length(subID_SIS))   
  beta_pval_fdr = p.adjust(beta_pval, 'fdr', length(subID_SIS))
  pval_bind_fdr = rbind(alpha_pval_fdr, beta_pval_fdr)
  join_pval_fdr = apply(pval_bind_fdr, 2, max)
  
  alpha_pval_bon = p.adjust(alpha_pval, 'bonferroni', length(subID_SIS))     
  beta_pval_bon = p.adjust(beta_pval, 'bonferroni', length(subID_SIS))
  pval_bind_bon = rbind(alpha_pval_bon, beta_pval_bon)
  join_pval_bon = apply(pval_bind_bon, 2, max)
  
  alpha_pval_by = p.adjust(alpha_pval, 'BY', length(subID_SIS))   
  beta_pval_by = p.adjust(beta_pval, 'BY', length(subID_SIS))
  pval_bind_by = rbind(alpha_pval_by, beta_pval_by)
  join_pval_by = apply(pval_bind_by, 2, max)
  
  if(verbose){
    message("        Significant FDR(BH) Joint mediator(s): ", names(join_pval_fdr)[which(join_pval_fdr < 0.05)])
    message("        Significant bonforroni Joint mediator(s): ", names(join_pval_bon)[which(join_pval_bon < 0.05)])
    message("        Significant FDR(BY) Joint mediator(s): ", names(join_pval_by)[which(join_pval_by < 0.05)])
  }
  
  result = data.frame(sob_pval_fdr, sob_pval_bon, sob_pval_by, 
                      join_pval_fdr, join_pval_bon, join_pval_by,
                      ab_coef, ab_var, conf_low, conf_up, sob_pval, 
                      alpha_coef, alpha_var, alpha_pval, 
                      beta_coef, beta_var, beta_pval)
  
  return(list(result = result, subID_SIS = subID_SIS, subM_SIS = subM_SIS, 
              SIS.beta = SIS.beta, SIS.alpha = SIS.alpha,
              id = id, wave = wave, X = X, Y = Y, COV = COV))
}



lineAlphaMarginal = function(X, M, COV, p, wave){
  # input:
  # X: exposure
  # M: mediators
  # COV: covariates
  # p: dimension of M
  
  Est_alpha = matrix(0, 3, p)
  colnames(Est_alpha) = colnames(M)
  X = X[which(wave == min(wave))]
  M = data.frame(M[which(wave == min(wave)),])
  COV = COV[which(wave == min(wave)),]
  for(j in 1:p){
    MX = data.frame(M = M[, j], X = X, COV = COV)
    fit = lm(M ~ ., data = MX)
    Est_alpha[1,j] = summary(fit)$coef[2,1]               # coefficients of alpha
    Est_alpha[2,j] = summary(fit)$coef[2,2]^2             # variance of alpha
    Est_alpha[3,j] = summary(fit)$coef[2,4]               # p-value of alpha
  }
  return(Est_alpha=Est_alpha)
}



lineBetaMarginal = function(Y, M, COV, p){
  # input:
  # Y: outcome
  # M: mediators
  # COV: covariates
  # p: dimension of M
  
  Est_beta = matrix(0, 3, p)
  colnames(Est_beta) = colnames(M)[1:p]
  for(j in 1:p){
    MY = data.frame(M = M[, j], Y = Y, COV = COV)
    fit = lm(Y ~ ., data = MY)
    Est_beta[1,j] = summary(fit)$coef[2,1]                # coefficients of beta
    Est_beta[2,j] = summary(fit)$coef[2,2]^2              # variance of beta
    Est_beta[3,j] = summary(fit)$coef[2,4]                # p-value of beta
  }
  return(Est_beta=Est_beta)
}
