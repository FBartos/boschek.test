#' Boschek's test of reliability
#'
#' @param data A matrix or a dataframe containing observations. 
#' It can be either in long format - with columns corresponding to a raters and rows to subjects rated, or,
#' in an aggregated manner, with one column labeled "frequency" containing counts of combinations of rating 
#' which are specified in remaining columns. 
#' @param Q 
#' @param model 
#' @param exact_prob 
#' @param print.all 
#'
#' @return
#' @export
#'
#' @examples
usethis::use_package('nleqslv')
boschek.test <- function(data, Q, model = "A", exact_prob = TRUE, print.all = FALSE){
  
  # check input
  if(any(apply(Q,1,sum) != 1) | any(apply(Q,2,sum) != 1))stop("Q misspecified - doesn't sum to 1")
  if(any(diag(Q) <= .5))stop("Assumptions not met - diagonal values of Q must be larger than .5")
  if(model != "A" & model != "B")stop("Select model type")
  if(!is.data.frame(data) & !is.matrix(data))stop("Incompatible data type")
  if(any(is.na(data)))stop("NA in data is not allowed")
  
  # reshape data
  if(any(colnames(data) == "frequency")){
    
    for(i in c(1:ncol(data))[colnames(data) != "frequency"]){
      data <- data[order(data[,i]),]
    }
    
    N_flatt      <- cbind.data.frame(data[,colnames(data) != "frequency"],"frequency" = data[,"frequency"])
    combinations <- apply(as.matrix(data[,colnames(data) != "frequency"]),2,as.numeric)
    
    NR <- ncol(data) - 1
    NK <- max(unlist(combinations))
    
    if(nrow(data) != NK^NR)stop("Supply frequency for all combinations of possible raters evaluations")
  }else{
    N_flatt      <- as.data.frame(table(as.data.frame(data)), responseName = "frequency") 
    combinations <- as.matrix(do.call(expand.grid,lapply(dim(table(as.data.frame(data))),seq)))
    
    NR <- ncol(data)
    NK <- max(dim(table(data)))
  }
  
  
  # solve T
  X <- nleqslv::nleqslv(rep(1/NK,NK-1),solve_T, N_flatt = N_flatt, combinations = combinations, Q = Q, NK = NK, NR = NR,
                        control = list(xtol = 1e-15, ftol = 1e-15))
  
  
  # compute residuals
  N_flatt <- cbind(N_flatt,"expected_frequency" = solve_P(X$x, combinations = combinations, Q = Q, NK = NK, NR = NR)*sum(N_flatt$frequency))
  N_flatt$standardized_residuals <- with(N_flatt,(frequency - expected_frequency)/sqrt(expected_frequency))
  
  N <- sum(N_flatt$frequency)
  
  # data description
  data_description <- list(
    "observations" = N,
    "number_of_categories" = NK,
    "number_of_raters" = NR
  )
  
  
  # model B
  combinations_help <- apply(combinations,1,function(x){paste(table(c(x,1:NK))-1,sep = "", collapse = "|")})
  
  N_flattB <- NULL
  symetry_chisq <- NULL
  for(com in unique(combinations_help)){
    temp_c <- c(1:nrow(N_flatt))[combinations_help == com]
    N_flattB <- rbind.data.frame(N_flattB,
                                 cbind.data.frame(com,
                                                  "frequency" = sum(N_flatt[temp_c,"frequency"]),
                                                  "expected_frequency" = sum(N_flatt[temp_c,"expected_frequency"])))
    
    # test of symetry
    if(!grepl(NK, com)){
      symetry_chisq <- c(symetry_chisq,(N_flatt[temp_c,"frequency"]-mean(N_flatt[temp_c,"frequency"]))^2/mean(N_flatt[temp_c,"frequency"]))
    }
  }
  colnames(N_flattB)[1] <- paste(c("evaluation (",
                                   paste(1:NK,sep = "", collapse = "|"),
                                   ")"),sep = "", collapse = "") 
  N_flattB$standardized_residuals <- with(N_flattB,(frequency - expected_frequency)/sqrt(expected_frequency))
  
  # test of symetry
  test_of_symmetry <- list(
    "X^2" = sum(symetry_chisq),
    "df"  = NK^NR-(factorial(NK+NR-1)/(factorial(NR)*factorial(NK-1))),
    "p.value" = pchisq(sum(symetry_chisq),NK^NR-(factorial(NK+NR-1)/(factorial(NR)*factorial(NK-1))),lower.tail = F)
  )
  
  ### solve exactly in case of two categories
  probabilities_estimates <- NULL
  if(NK == 2 & exact_prob){
    
    TQ <- nleqslv::nleqslv(c(
      N_flattB$frequency[1]/(N_flattB$frequency[1] + N_flattB$frequency[nrow(N_flattB)]),
      (sum(N_flattB$frequency) - N_flattB$frequency[1] - N_flattB$frequency[nrow(N_flattB)])/(NR*sum(N_flattB$frequency))),
      solve_TQ, N_flatt = N_flattB, NR = NR, control = list(xtol = 1e-15, ftol = 1e-15))$x
    
    # check that the solution satisfies boundaries
    if((TQ[1] >= 1 | TQ[1] <= 0) | TQ[2] >= .5 | TQ[1] <= 0){
      # if not, try alternative starting values
      TQ_seq <- sapply(seq(0.01,0.49,0.01),function(q){
        nleqslv::nleqslv(c(
          N_flattB$frequency[1]/(N_flattB$frequency[1] + N_flattB$frequency[nrow(N_flattB)]),q),
          solve_TQ, N_flatt = N_flattB, NR = NR, control = list(xtol = 1e-15, ftol = 1e-15))$x
      })
      TQ_seq <- TQ_seq[,!(TQ_seq[1,] >= 1 | TQ_seq[1,] <= 0) & !(TQ_seq[2,] >= .5 | TQ_seq[2,] <= 0)]
      if(ncol(TQ_seq) == 0)stop("Solution to routine B not found.")
      TQ_seq <- TQ_seq[,!duplicated(round(TQ_seq[2,],3))]
      if(is.matrix(TQ_seq))stop("Multiple solutions to routine B found.")
    }
    
    t1 <- TQ[1]
    t2 <- 1 - TQ[1]
    # compute SE
    t1_se <- compute_SE(t1,TQ[2], N)
    t2_se <- compute_SE(t2,TQ[2], N)
    
    probabilities_estimates <- list(
      "T1" = list(
        "estimate" = t1,
        "se" = t1_se,
        "95% CI" = c(ifelse(t1-1.96*t1_se < 0, 0, t1-1.96*t1_se),
                     ifelse(t1+1.96*t1_se > 1, 1, t1+1.96*t1_se))
      ),
      "T2" = list(
        "estimate" = t2,
        "se" = t2_se,
        "95% CI" = c(ifelse(t2-1.96*t2_se < 0, 0, t2-1.96*t2_se),
                     ifelse(t2+1.96*t2_se > 1, 1, t2+1.96*t2_se))
      )
    )
  }
  
  # chi squared test
  chi.sq <- list(
    "X^2" = sum(with(N_flatt,(frequency - expected_frequency)^2/expected_frequency)),
    "df"  = NK^NR - NK,
    "p.value" =   pchisq(sum(with(N_flatt,(frequency - expected_frequency)^2/expected_frequency)),NK^NR - NK, lower.tail = F)
  )
  
  chi.sqB <- list(
    "X^2" = sum(with(N_flattB,(frequency - expected_frequency)^2/expected_frequency)),
    "df"  = factorial(NK+NR-1)/(factorial(NR)*factorial(NK-1)) - NK,
    "p.value" =   pchisq(sum(with(N_flattB,(frequency - expected_frequency)^2/expected_frequency)),
                         factorial(NK+NR-1)/(factorial(NR)*factorial(NK-1)) - NK, lower.tail = F)
  )  
  
  
  # console print code
  print_probabilities_estimates <- NULL
  print_test_of_symmetry        <- NULL
  print_data_summary_add        <- NULL
  print_data_summary_add1       <- NULL
  if(model == "A"){
    model_info <- "test .... - přidat další popis který se má zobrazit uživatelům"
    N_flatt_print <- N_flatt
    N_flatt_print[,c(NR+2):ncol(N_flatt_print)] <- apply(N_flatt_print[c(NR+2):ncol(N_flatt_print)],2,function(x)as.numeric(sprintf("%.3f", x)))
    
    colnames(N_flatt_print)[1:NR] <- sapply(1:NR,function(x)paste(c("R",x),collapse = ""))
    colnames(N_flatt_print)[(NR+1):ncol(N_flatt_print)] <- c("Ob Frq","Exp Frq","Std Res")
    
    if(print.all){
      print_data_summary_rows <- 1:nrow(N_flatt_print)
    }else{
      print_data_summary_rows <- 1:ifelse(nrow(N_flatt_print)<10,nrow(N_flatt_print),10)
      if(NR > 7){
        print_data_summary_add1 <- paste(c(' and ',ncol(N_flatt_print)-10,' columns'),collapse = "")
        N_flatt_print <- N_flatt_print[,c(1:8,(ncol(N_flatt_print)-2):ncol(N_flatt_print))]
        N_flatt_print[,8] <- "..."
        colnames(N_flatt_print)[8] <- "..."
      }
      if(nrow(N_flatt_print) - 10 > 0)print_data_summary_add  <- paste(c('... (shortened by ',nrow(N_flatt_print) - 10,' rows',print_data_summary_add1,')','\n'),collapse = "")
    }
    print_data_summary <- c(
      c(paste(colnames(N_flatt_print),collapse = '\t'),'\n'),
      apply(N_flatt_print[print_data_summary_rows,],1,function(x)c(paste(x,collapse = '\t'),'\n')),
      print_data_summary_add)
    
    print_model_fit <- paste(c(
      'X^2(',chi.sq$df,') = ',sprintf("%.2f", chi.sq$`X^2`),', ',
      ifelse(chi.sq$p.value < .001, 'p < .001', paste(c('p = ', sprintf("%.3f", chi.sq$p.value)),collapse = "")),
      '\n'),sep = "",collapse = "")
  }else if(model == "B"){
    model_info <- "test assuming symmetry - přidat další popis který se má zobrazit uživatelům"
    N_flattB_print <- N_flattB
    N_flattB_print[,3:4] <- apply(N_flattB_print[,3:4],2,function(x)as.numeric(sprintf("%.3f", x)))
    
    colnames(N_flattB_print)[1] <- "N Cat"
    colnames(N_flattB_print)[2:ncol(N_flattB_print)] <- c("Ob Frq","Exp Frq","Std Res")
    
    if(print.all){
      print_data_summary_rows <- 1:nrow(N_flattB_print)
    }else{
      print_data_summary_rows <- 1:ifelse(nrow(N_flattB_print)<10,nrow(N_flattB_print),10)
      if(nrow(N_flattB_print) - 10 > 0)print_data_summary_add  <- paste(c('... (shortened by ',nrow(N_flattB_print) - 10,' rows)','\n'),collapse = "")
    }
    print_data_summary <- c(
      c(paste(colnames(N_flattB_print),collapse = '\t'),'\n'),
      apply(N_flattB_print[print_data_summary_rows,],1,function(x)c(paste(x,collapse = '\t'),'\n')),
      print_data_summary_add)
    
    if(NK == 2 & exact_prob){
      ep <- t(sapply(probabilities_estimates,function(x)sprintf("%.3f", unlist(x))))
      print_probabilities_estimates <- paste(c(
        'Estimated probabilities:','\n',
        ' ','\t','Est','\t','SE','\t','95% CI','\n',
        ' T1','\t',ep[1,1],'\t',ep[1,2],'\t','[',ep[1,3],', ',ep[1,4],']','\n',
        ' T2','\t',ep[2,1],'\t',ep[2,2],'\t','[',ep[2,3],', ',ep[2,4],']','\n',
        '\n'
      ),sep = "",collapse = "")
    }
    print_test_of_symmetry <- paste(c(
      '\n',
      ' Test of symmetry:','\n',
      ' X^2(',test_of_symmetry$df,') = ',sprintf("%.2f", test_of_symmetry$`X^2`),', ',
      ifelse(test_of_symmetry$p.value < .001, 'p < .001', paste(c('p = ', sprintf("%.3f", test_of_symmetry$p.value)),collapse = "")),
      '\n'),sep = "",collapse = "")
    
    print_model_fit <- paste(c(
      'X^2(',chi.sqB$df,') = ',sprintf("%.2f", chi.sqB$`X^2`),', ',
      ifelse(chi.sqB$p.value < .001, 'p < .001', paste(c('p = ', sprintf("%.3f", chi.sqB$p.value)),collapse = "")),
      '\n'),sep = "",collapse = "")
  }
  
  
  # print to console
  output <- c(
    "Boschek's test of reliability",'\n',
    'Model',model,": ",model_info,'\n',
    '\n',
    'Data description:','\n',
    'Observations: ',data_description$observations,'\n',
    'Categories:   ',data_description$number_of_categories,'\n',
    'Raters:       ',data_description$number_of_raters,'\n',
    '\n',
    'Data summary:','\n',
    print_data_summary,
    '\n',
    print_probabilities_estimates,
    'Model fit:','\n',
    print_model_fit,
    print_test_of_symmetry
  )
  cat(output)
  
  # return silently
  if(model == "A"){
    return(invisible(list(
      "summary" = output,
      "model" = model,
      "data_description" = data_description,
      "table" = N_flatt,
      "allowed_classifications" = Q,
      "test_of_model_fit" = chi.sq
    )))
  }else if(model == "B"){
    return(invisible(list(
      "summary" = output,
      "model" = model,
      "data_description" = data_description,
      "table" = N_flattB,
      "allowed_classifications" = Q,
      "probabilities_estimates" = probabilities_estimates,
      "test_of_model_fit" = chi.sqB,
      "test_of_symmetry" = test_of_symmetry
    )))
  }
}
