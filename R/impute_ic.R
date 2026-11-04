#' Time Imputation for Interval-Censored Data
#' 
#' @description Time Imputation for Interval-Censored Data
#'
#' @param formula Formula of the type Surv(left, right) ~ treatment (or ~ 1) 
#' @param data Data frame with the original data
#' @param method A character string specifying the method or model to be used. 
#'   Supported options are:
#'   \itemize{
#'     \item \code{"midpoint"}, \code{"left"}: Simple imputation methods.
#'     \item \code{"lognormal"}, \code{"weibull"}, \code{"exponential"}, 
#'     \code{"gaussian"}: Parametric survival models.
#'     \item \code{"xgboost"}: eXtreme Gradient Boosting.
#'     \item \code{"RSF"}: Random Survival Forests.
#'     \item \code{"survivalsvm"}: Support Vector Regression optimized with a 
#'     linear kernel to predict survival times.
#'     \item \code{"glmboost"}: Boosting-based Cox regression with 
#'     component-wise linear base-learners.
#'     \item \code{"gamboost"}: Boosting-based Cox regression automatically 
#'     parsing smooth (\code{bbs}) and linear (\code{bols}) terms.
#'   }
#' @param n_rep Number of replicates for the lognorm, weibull, exponential, 
#' gaussian, xgboost methods (default 1000)
#' @param nrounds Number of boosting iterations for the \code{"xgboost"} method 
#' (default 100)
#' @param ntree Number of trees to grow for the \code{"RSF"} method 
#' (default 100)
#' @return The original data frame with two extra columns: 'time_imp' and 
#' 'status_imp'
#' 
#' @examples
#' \dontrun{

#' data(breastCancer)
#' head(breastCancer)
#' dim(breastCancer)
#' 
#' breastCancer[1:6, c("rectime", "censrec", "survtime", "censdead", "left_rec", 
#'                     "right_rec", "left", "right", "age", "menopause")]
#' 
#' imputed_mid <- impute_ic(Surv(left, right) ~ 1, 
#'                data = breastCancer, 
#'                method = "midpoint")
#' head(imputed_mid[, c("left", "right", "time_imp", "status_imp")])
#' 
#' imputed_lnorm<- impute_ic(Surv(left, right) ~ age + menopause + 
#'                 hormone + size + grade + nodes + prog_recp + 
#'                 estrg_recp, data = breastCancer, method = "lognormal", 
#'                 n_rep = 1000)
#' head(imputed_lnorm[, c("left", "right", "time_imp", "status_imp")])
#' 
#' imputed_xgb <- impute_ic(Surv(left, right_rec) ~ age + menopause + 
#'                hormone + size + grade + nodes + prog_recp + estrg_recp, 
#'                data = breastCancer, method = "xgboost", n_rep = 1000, 
#'                nrounds = 50) 
#'                
#' imputed_rsf <- impute_ic(Surv(left, right_rec) ~ age + menopause + 
#'                hormone + size + grade + nodes + prog_recp + estrg_recp, 
#'                data = breastCancer, method = "RSF", n_rep = 1000,
#'                ntree = 100) 
#' 
#' imputed_svm <- impute_ic(Surv(left, right_rec) ~ age + menopause +
#'                hormone + size + grade + nodes + prog_recp + estrg_recp, 
#'                data = breastCancer, method = "survivalsvm",  
#'                n_rep = 50)
#' 
#' fit_pm<- fit_ic(imputed_mid$time_imp, imputed_mid$status_imp)
#' 
#' summary(fit_pm, times = c(365, 730, 1095, 1460))
#' 
#' fit_xgb<- fit_ic(imputed_wei$time_imp, imputed_wei$status_imp)
#' 
#' summary(fit_xgb, times = c(365, 730, 1095, 1460))
#' 
#' ints <- get_intervals(breastCancer$left, breastCancer$right_rec\)
#' times <- seq(0, max(breastCancer$survtime), 1)
#' tb_surv <- survival_ic(ints, times = times, conf = TRUE) 
#' 
#' plot_ic(tb_surv[1:2600,], conf = T, showlegend = F)
#' lines(fit_pm, col=1, lty=1, conf.int = F) 
#' lines(fit_lnorm, col=2, lty=2, conf.int = F) 
#' lines(fit_wei, col=3, lty=2, conf.int = F) 
#' lines(fit_exp, col=4, lty=2, conf.int = F) 
#' lines(fit_gau, col=5, lty=2, conf.int = F)
#' legend("topright", legend = c("TB", "Midpoint", "Log-normal", "Weibull",
#'        "Exponential", "Gauss"), lwd = 1.2, 
#'        col = c(2, 1, 2, 3,  4, 5), cex = 0.7, lty = c(1,1,2,2, 2, 2))
#' }

impute_ic<- function(formula, data, method = "midpoint", n_rep = 1000, 
                          nrounds = 100, ntree = 100){
  
  termo_lhs <- formula[[2]] 
  
  left_name  <- as.character(termo_lhs[[2]])
  
  right_name <- as.character(termo_lhs[[3]])
  
  if(method %in% c("midpoint", "left", "right")){
    
    left  <- as.numeric(data[[left_name]])
    right <- as.numeric(data[[right_name]])
    
    right_rec<-ifelse(right=='Inf', NA, right)
    
    is_censored <- is.infinite(right) | is.na(right)
    
  }else{
    
    dat.ini<-data
    
    dat.na<-data[is.na(data[[right_name]]),]  #data[[right_name]]  sao if e nao NA  
    
    data<-data[!is.na(data[[right_name]]),]
    
    dim(data) #686 --- todos os registos 
    
    left  <- as.numeric(data[[left_name]])
    
    right <- as.numeric(data[[right_name]])
    
    right_rec<-ifelse(right=='Inf', NA, right)
    
    is_censored <- is.infinite(right) | is.na(right)
    
  }

  if (method == "midpoint") {
    #time_imp2 <- ifelse(!is_censored, (left + right_rec) / 2, left)
    time_imp <- ifelse(!is_censored, (left + right) / 2, left)
    status_imp <- ifelse(!is_censored, 1, 0)
    
  }else{ 
    
    if(method=='left'){
      
      time_imp <- ifelse(!is_censored,left, left)
      status_imp <- ifelse(!is_censored, 1, 0)
    }
    
    if(method=='right'){
      
      time_imp <- ifelse(!is_censored, right , left)
      status_imp <- ifelse(!is_censored, 1, 0)
      
    }
    
    if(method %in% c("lognormal", "weibull", "exponential", "gaussian")){
      
      #set.seed(123)
      
      right_calc <- ifelse(is_censored, 1e9, right)
      
      rhs <- if(length(formula) == 3) as.character(formula[3]) else "1" 
      
      bs_fits <- lapply(1:n_rep, function(k) {
        
        cat('Replica: ', k, '\n')
        
        data$times_sim <- runif(nrow(data), min = left, max = right_calc)
        
        data$status_sim <- rep(1, nrow(data))
        
        rhs <- if(length(formula) == 3) as.character(formula[3]) else "1"
        
        f_sim <- as.formula(paste("survival::Surv(times_sim, status_sim) ~", 
                                  rhs))
        
        survival::survreg(f_sim, data = data, dist = method) 
        
      })
      
      time_imp <- sapply(1:nrow(data), function(j) {
        preds <- sapply(bs_fits, function(fit) predict(fit)[j])
        t_min <- min(preds); t_max <- max(preds)
        if(t_max == t_min) return(left[j])
        val <- left[j] + (preds - t_min) / (t_max - t_min) * (right_calc[j] - 
                                                                left[j])
        return(median(val, na.rm = TRUE))
      })
      
      data<-rbind(data, dat.na)
      
      #time_imp 
      
      status_imp <- ifelse(is_censored, 0, 1)
      
      time_imp<-c(time_imp, dat.na$left)
      
      status_imp<-c(status_imp, rep(0, length(dat.na$left)))
      
      time_imp<-ifelse(is_censored, data$left, time_imp)
      
    }
    
    if(method == "xgboost"){
      
      #set.seed(123)
      
      right_calc <- ifelse(is_censored, 1e9, right)
      
      rhs_formula <- if(length(formula) == 3) formula[c(1, 3)] else ~ 1
      
      X <- model.matrix(rhs_formula, data = data)
      
      all_preds <- matrix(NA, nrow = nrow(data), ncol = n_rep)
      
      for (k in 1:n_rep) {
      
        cat('Replica: ', k, '\n')
      
        times_sim <- sapply(1:nrow(data), function(i) runif(1, left[i], 
                                                            right_calc[i]))
        
        dtrain <- xgb.DMatrix(data = X)
        setinfo(dtrain, "label_lower_bound", times_sim)
        setinfo(dtrain, "label_upper_bound", times_sim)
        
        params <- list(
          objective = "survival:aft",
          eval_metric = "aft-nloglik",
          aft_loss_distribution = "normal", # normal, logistic, ou extreme
          aft_loss_distribution_scale = 1.0,
          tree_method = "hist"
        )
        
        bst <- xgb.train(params = params, data = dtrain, nrounds = nrounds, 
                         verbose = 0)
        all_preds[, k] <- predict(bst, newdata = dtrain)
        
      }
      
      time_imp <- sapply(1:nrow(data), function(j) {
        
        preds_j <- all_preds[j, ]
        t_min <- min(preds_j); t_max <- max(preds_j)
        
        if(t_max == t_min) return(left[j])
  
        val <- left[j] + (preds_j - t_min) / (t_max - t_min) * (right_calc[j] - 
                                                                  left[j])
        return(median(val, na.rm = TRUE))
      })
      
      data<-rbind(data, dat.na)
      
      #time_imp 
      
      status_imp <- ifelse(is_censored, 0, 1)
      
      time_imp<-c(time_imp, dat.na$left)
      
      status_imp<-c(status_imp, rep(0, length(dat.na$left)))
      
    }
    
    if(method=='RSF'){
      
      right_calc <- ifelse(is_censored, 1e9, right)
      
      if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
        stop("The 'randomForestSRC' package is required for the RSF method.")
      }
      
      #set.seed(123)
      
      all_preds <- matrix(NA, nrow = nrow(data), ncol = n_rep)
      
      rhs_vars <- all.vars(formula[[3]])
      
      for (k in 1:n_rep) {
        
        cat('Replica: ', k, '\n')
        
        df_rsf <- data.frame(
          Y1 = runif(nrow(data), min = left, max = right_calc),
          Y2 = rep(1, nrow(data))
        )
        
        df_rsf <- cbind(df_rsf, data[, rhs_vars, drop = FALSE])
        

        f_str <- paste(
          "Surv(Y1, Y2) ~",
          paste(rhs_vars, collapse = "+")
        )
        
        fit_rsf <- randomForestSRC::rfsrc(
          formula = as.formula(f_str), 
          data = df_rsf, 
          ntree = ntree,
          perf.type = "none"
        )
        
        all_preds[, k] <- as.numeric(fit_rsf$predicted)
      }
      
      time_imp <- sapply(1:nrow(data), function(j) {
        preds_j <- all_preds[j, ]
        t_min <- min(preds_j); t_max <- max(preds_j)
        if(t_max == t_min) return(left[j])
        val <- left[j] + (preds_j - t_min) / (t_max - t_min) * (right_calc[j] - 
                                                                  left[j])
        return(median(val, na.rm = TRUE))
      })
      
      data<-rbind(data, dat.na)
      
      #time_imp 
      
      status_imp <- ifelse(is_censored, 0, 1)
      
      time_imp<-c(time_imp, dat.na$left)
      
      status_imp<-c(status_imp, rep(0, length(dat.na$left)))
    
    }
    
    if(method == "survivalsvm"){
      
      # install.packages("survivalsvm")   # se ainda nao estiver instalado
      if (!requireNamespace("survivalsvm", quietly = TRUE)) {
        stop("The 'survivalsvm' package is required for the survivalsvm 
             method.")
      }
      
      right_calc <- ifelse(is_censored, 1e9, right)
      
      rhs_vars <- all.vars(formula[[3]])
      
      all_preds <- matrix(NA, nrow = nrow(data), ncol = n_rep)
      
      for (k in 1:n_rep) {
        
        #k<-1
        
        cat('Replica: ', k, '\n')
        
        df_svm <- data.frame(
          Y1 = runif(nrow(data), min = left, max = right_calc),
          Y2 = rep(1, nrow(data))
        )
        df_svm <- cbind(df_svm, data[, rhs_vars, drop = FALSE])
        
        # Remover eventuais variaveis constantes na replica (mesma protecao do 
        #gamboost)
        keep <- rhs_vars[vapply(rhs_vars, function(v) {
          x <- df_svm[[v]]
          if (is.factor(x)) length(levels(droplevels(x))) >= 2
          else length(unique(x[!is.na(x)])) >= 2
        }, logical(1))]
        
        f_str <- paste("Surv(Y1, Y2) ~", paste(keep, collapse = " + "))
        
        fit_svm <- tryCatch(
          survivalsvm::survivalsvm(
            formula     = as.formula(f_str),
            data        = df_svm,
            type        = "regression", # produz previsoes na escala do tempo
            gamma.mu    = 0.1,          # parametro de regularizacao 
            opt.meth    = "quadprog",   # alternativa: "ipop" (mais robusto, 
                                        #mais lento)
            kernel      = "lin_kernel"  # alternativas: "add_kernel", 
                                        # "rbf_kernel", "poly_kernel"
          ),
          error = function(e) {
            message("Replica failed ", k, ": ", conditionMessage(e))
            NULL
          }
        )
        
        if (is.null(fit_svm)) next
        
        pred_svm <- predict(fit_svm, newdata = df_svm)
        # predict.survivalsvm devolve um objeto S3; o vetor de previsoes esta 
        #em $predicted
        all_preds[, k] <- as.numeric(pred_svm$predicted)
      }
      
      time_imp <- sapply(1:nrow(data), function(j) {
        preds_j <- all_preds[j, ]
        preds_j <- preds_j[!is.na(preds_j)]
        if (length(preds_j) == 0) return(left[j])
        t_min <- min(preds_j); t_max <- max(preds_j)
        if (t_max == t_min) return(left[j])
        val <- left[j] + (preds_j - t_min) / (t_max - t_min) * (right_calc[j] - 
                                                                  left[j])
        return(median(val, na.rm = TRUE))
      })
      
      data <- rbind(data, dat.na)
      
      status_imp <- ifelse(is_censored, 0, 1)
      time_imp   <- c(time_imp, dat.na$left)
      status_imp <- c(status_imp, rep(0, length(dat.na$left)))
    }
    
    if(method == "glmboost"){

      right_calc <- ifelse(is_censored, 1e9, right)
      
      bs_fits <- lapply(1:n_rep, function(k) {
        
        #k<-1
        
        cat('Replica: ', k, '\n')
        
        data$times_sim <- runif(nrow(data), min = left, max = right_calc)
        
        data$status_sim <- rep(1, nrow(data))
        
        rhs <- if(length(formula) == 3) as.character(formula[3]) else "1"
        
        #f_sim <- as.formula(paste("survival::Surv(times_sim, status_sim) ~", 
        #rhs))
        
        #survival::survreg(f_sim, data = data, dist = method) 
        
        f_sim <- as.formula(paste("Surv(times_sim, status_sim) ~", rhs))
        
        glmboost(f_sim, data=data, family = CoxPH(), 
                 control = boost_control(mstop = 100, nu = 0.1))
        
      })
      
      time_imp <- sapply(1:n_rep, function(k) {
        fit <- bs_fits[[k]]
        sf <- survFit(fit, newdata = data)   # data com todas as linhas
        # sf$surv é matriz: linhas = tempos, colunas = individuos
        apply(sf$surv, 2, function(s) {
          idx <- which(s <= 0.5)[1]
          if (is.na(idx)) return(max(sf$time))
          sf$time[idx]
        })
      })
      # time_imp fica matriz: linhas = indivíduos, colunas = replicas
      # (se sapply devolver por colunas, faz t() conforme necessario)
      
      # Agora aplica a normalizacao por individuo entre replicas:
      time_imp_final <- sapply(1:nrow(data), function(j) {
        preds <- time_imp[j, ]
        t_min <- min(preds); t_max <- max(preds)
        if (t_max == t_min) return(left[j])
        val <- left[j] + (preds - t_min) / (t_max - t_min) * (right_calc[j] - 
                                                                left[j])
        return(median(val, na.rm = TRUE))
      })
      
      length(time_imp_final)
      
      data<-rbind(data, dat.na)
      
      dim(data)
      
      #time_imp 
      
      status_imp <- ifelse(is_censored, 0, 1)
      
      length(status_imp)
      
      time_imp<-c(time_imp_final, dat.na$left)
      
      length(time_imp)
      
      status_imp<-c(status_imp, rep(0, length(dat.na$left)))
      
      time_imp<-ifelse(is_censored, data$left, time_imp)  #acrescentei
      
    }    
    
    if(method == "gamboost"){
      
      # ---- 1. Parse da formula do utilizador ----------------------------
      # Objetivo: identificar termos bbs(...) / bols(...) e termos "nus"
      
      rhs_terms <- attr(terms(formula), "term.labels")
      # ex.: c("bbs(age)", "bbs(nodes)", "rx", "sex", "obstruct", ...)
      
      # Termos com wrapper de base learner do mboost
      is_smooth <- grepl("^bbs\\(",  rhs_terms)
      is_linear <- grepl("^bols\\(", rhs_terms)
      is_wrapped <- is_smooth | is_linear
      
      # Extrai o nome da variavel dentro do wrapper (ex.: "bbs(age)" -> "age")
      #extract_var <- function(term) {
      #  sub("^[a-zA-Z]+\\(\\s*([^,\\)\\s]+).*", "\\1", term)
      #}
      
      
      #Correcao simples: use uma regex mais limpa, sem \s dentro da classe 
      #negada:
        
        
      #  extract_var <- function(term) {
      #    # remove o wrapper "nomeFuncao(" no inicio e tudo a partir de 
      #vírgula ou ")"
      #    x <- sub("^[a-zA-Z]+\\(\\s*", "", term)   # tira "bbs(  "
      #    x <- sub("[,\\)].*$", "", x)              # tira ", ..." ou ")"
      #    trimws(x)
      #  }
      
      
      #Ou, ainda mais robusto (usando all.vars do proprio R, sem regex):"
      extract_var <- function(term) {
        v <- all.vars(parse(text = term))
        v[1]   # primeira variável dentro do wrapper
      }
      
      smooth_vars <- vapply(rhs_terms[is_smooth], extract_var, character(1))
      linear_vars <- vapply(rhs_terms[is_linear], extract_var, character(1))
      factor_vars <- rhs_terms[!is_wrapped]   # variaveis "nuas" -> factores
      
      # Variaveis continuas = as que aparecem em bbs() ou bols() 
      #com var numerica
      numeric_vars <- unname(c(smooth_vars, linear_vars))
      
      # ---- 2. Coercao automatica das colunas ---------------------------
      for (v in numeric_vars) {
        if (v %in% names(data)) data[[v]] <- as.numeric(data[[v]])
      }
      
      for (v in factor_vars) {
        if (v %in% names(data)) data[[v]] <- as.factor(data[[v]])
      }
      # Garantir também que left/right são numericos
      data[[left_name]]  <- as.numeric(data[[left_name]])
      data[[right_name]] <- as.numeric(data[[right_name]])
      
      # ---- 3. Reconstrucao da formula "expandida" ----------------------
      # Utilizador escreveu bbs(age) -> expandimos para bbs(age, knots=10, df=3)
      # Se ele proprio já passou argumentos (ex.: bbs(age, df=4)), respeitamos.
      
      rebuild_smooth <- function(term) {
        if (grepl("^bbs\\([^,]+\\)\\s*$", term)) {
          var <- extract_var(term)
          # heurística: menos knots se poucos valores únicos
          nu <- length(unique(data[[var]]))
          k  <- if (nu < 10) max(3, nu - 2) else 10
          sprintf("bbs(%s, knots = %d, df = 3)", var, k)
        } else {
          term  # já tem argumentos definidos pelo utilizador
        }
      }
      
      smooth_expanded <- vapply(rhs_terms[is_smooth], 
                                rebuild_smooth, character(1))
      linear_expanded <- rhs_terms[is_linear]                   # inalterados
      factor_expanded <- sprintf("bols(%s, df = 1)", factor_vars) # nus -> bols
      
      rhs_full <- paste(c(smooth_expanded, linear_expanded, factor_expanded),
                        collapse = " + ")
      
      # Formula final que sera usada dentro do loop bootstrap
      formula_full <- as.formula(
        paste("Surv(times_sim, status_sim) ~", rhs_full)
      )
      
      # ---- 4. Preparacao de left/right/censura -------------------------
      left  <- as.numeric(data[[left_name]])
      right <- as.numeric(data[[right_name]])
      
      is_censored <- is.infinite(right) | is.na(right)
      right_calc  <- ifelse(is_censored, 1e9, right)
      
      # ---- 5. Loop bootstrap: ajuste de gamboost por replica -----------
      bs_fits <- lapply(1:n_rep, function(k) {
        cat("Replica: ", k, "\n")
        
        data$times_sim  <- runif(nrow(data), min = left, max = right_calc)
        data$status_sim <- rep(1, nrow(data))
        
        suppressWarnings(
        
        tryCatch(
          mboost::gamboost(formula_full, data = data,
                           family = CoxPH(),
                           control = boost_control(mstop = 100, nu = 0.1)),
          error = function(e) {
            message("Replica failed ", k, ": ", conditionMessage(e))
            NULL
          }
        )
        
        )
      })
      
      # ---- 6. Predicao do tempo mediano por replica --------------------
      time_imp <- sapply(1:n_rep, function(k) {
        fit <- bs_fits[[k]]
        if (is.null(fit)) return(rep(NA_real_, nrow(data)))
        sf <- survFit(fit, newdata = data)
        apply(sf$surv, 2, function(s) {
          idx <- which(s <= 0.5)[1]
          if (is.na(idx)) return(max(sf$time))
          sf$time[idx]
        })
      })
      # time_imp: linhas = individuos, colunas = replicas
      
      # ---- 7. Normalizacao por individuo entre replicas ----------------
      time_imp_final <- sapply(1:nrow(data), function(j) {
        preds <- time_imp[j, ]
        t_min <- min(preds, na.rm = TRUE)
        t_max <- max(preds, na.rm = TRUE)
        if (!is.finite(t_min) || !is.finite(t_max) || t_max == t_min) 
          return(left[j])
        val <- left[j] + (preds - t_min) / (t_max - t_min) * (right_calc[j] - 
                                                                left[j])
        return(median(val, na.rm = TRUE))
      })
      
      # ---- 8. Reincorporacao das linhas com right = NA -----------------
      data <- rbind(data, dat.na)
      
      status_imp <- ifelse(is_censored, 0, 1)
      time_imp   <- c(time_imp_final, dat.na[[left_name]])
      status_imp <- c(status_imp, rep(0, nrow(dat.na)))
      time_imp<-ifelse(is_censored, data$left, time_imp)  #acrescentei
    }
  
  }
  
  res_data <- data
  res_data$time_imp <- as.numeric(time_imp)
  res_data$status_imp <- as.numeric(status_imp)
  
  return(res_data)
  
}