#' Boschek's test of reliability
#' @description Implementation of the Boschek's test of reliability as described in "Pape name".
#' @references Boschek, P., Vranka, M. A., Bartoš, F., Name of the Paper, 2019
#' @import nleqslv
#' @param data A matrix or a dataframe containing observations. 
#' It can be either in long format - with columns corresponding to a raters and rows to subjects rated, or,
#' in an aggregated manner, with one column labeled "frequency" containing counts of combinations of ratings 
#' which are specified in remaining columns. 
#' @param Q Classification matrix against which are the data tested.
#' @param exact_prob Whether the exact probability estimates should be computed. Only possible when number
#' of categories is 2.
#'
#'
#' @examples ### 2 raters, 2 categories and data in long format
#' # load data1
#' data(d1)
#' 
#' # create classification matrix
#' Q1 <- matrix(c(.9, .1,
#'                .1, .9),byrow = TRUE, nrow = 2, ncol = 2)
#'                
#' # fit the model
#' m1 <- boschek.test(data = d1,Q = Q1)
#' 
#' # quickly inspect the model
#' m1
#' 
#' # for more details use summary function
#' summary(m1)
#' 
#' # to print residuals and aggregated data matrix use residuals function
#' residuals(m1)
#' 
#' # in order to print all columns and rows set argument print.all = TRUE
#' residuals(m1, print.all = TRUE)
#' 
#' # or access the residuals directly from the fitted object
#' m1$model_A$table
#' m1$model_B$table
#' 
#' ### ------------------------------------------------- ###
#' ### 3 raters, 3 categories and data in long format
#' # load data1
#' data(d2)
#' 
#' # note that column containing the number of rattings is called frequency
#' head(d2)
#' 
#' # create classification matrix
#' Q2 <- matrix(c(0.85,  0.1, 0.05,
#'                0.1 ,  0.8, 0.1 ,
#'                0.05,  0.1, 0.85 ),byrow = TRUE, nrow = 3, ncol = 3)
#' 
#' # fit the model
#' m2 <- boschek.test(data = d2,Q = Q2)
#' 
# check input
#' @export boschek.test
boschek.test <- function(data, Q, exact_prob = TRUE){
  
  call <- match.call()
  
  if(any(apply(Q,1,sum) != 1) | any(apply(Q,2,sum) != 1))stop("Q misspecified - doesn't sum to 1")
  if(any(diag(Q) <= .5))stop("Assumptions not met - diagonal values of Q must be larger than .5")
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
  X <- nleqslv::nleqslv(rep(1/NK,NK-1),.solve_T, N_flatt = N_flatt, combinations = combinations, Q = Q, NK = NK, NR = NR,
                        control = list(xtol = 1e-15, ftol = 1e-15))
  
  
  # compute residuals
  N_flatt <- cbind(N_flatt,"expected_frequency" = .solve_P(X$x, combinations = combinations, Q = Q, NK = NK, NR = NR)*sum(N_flatt$frequency))
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
    "p.value" = stats::pchisq(sum(symetry_chisq),NK^NR-(factorial(NK+NR-1)/(factorial(NR)*factorial(NK-1))),lower.tail = F)
  )
  
  ### solve exactly in case of two categories
  probabilities_estimates <- NULL
  if(NK == 2 & exact_prob){
    
    TQ <- nleqslv::nleqslv(c(
      N_flattB$frequency[1]/(N_flattB$frequency[1] + N_flattB$frequency[nrow(N_flattB)]),
      (sum(N_flattB$frequency) - N_flattB$frequency[1] - N_flattB$frequency[nrow(N_flattB)])/(NR*sum(N_flattB$frequency))),
      .solve_TQ, N_flatt = N_flattB, NR = NR, control = list(xtol = 1e-15, ftol = 1e-15))$x
    
    # check that the solution satisfies boundaries
    if((TQ[1] >= 1 | TQ[1] <= 0) | TQ[2] >= .5 | TQ[1] <= 0){
      # if not, try alternative starting values
      TQ_seq <- sapply(seq(0.01,0.49,0.01),function(q){
        nleqslv::nleqslv(c(
          N_flattB$frequency[1]/(N_flattB$frequency[1] + N_flattB$frequency[nrow(N_flattB)]),q),
          .solve_TQ, N_flatt = N_flattB, NR = NR, control = list(xtol = 1e-15, ftol = 1e-15))$x
      })
      TQ_seq <- TQ_seq[,!(TQ_seq[1,] >= 1 | TQ_seq[1,] <= 0) & !(TQ_seq[2,] >= .5 | TQ_seq[2,] <= 0)]
      if(ncol(TQ_seq) == 0)stop("Solution to routine B not found.")
      TQ_seq <- TQ_seq[,!duplicated(round(TQ_seq[2,],3))]
      if(is.matrix(TQ_seq))stop("Multiple solutions to routine B found.")
    }
    
    t1 <- TQ[1]
    t2 <- 1 - TQ[1]
    # compute SE
    t1_se <- .compute_SE(t1,TQ[2], N)
    t2_se <- .compute_SE(t2,TQ[2], N)
    
    probabilities_estimates <- rbind(
      "T1" = c(
        "Estimate" = t1,
        "SE"       = t1_se,
        "95% l.CI" = ifelse(t1-1.96*t1_se < 0, 0, t1-1.96*t1_se),
        "95% u.CI" = ifelse(t1+1.96*t1_se > 1, 1, t1+1.96*t1_se)
      ),
      "T2" = c(
        "Estimate" = t2,
        "SE"       = t2_se,
        "95% l.CI" = ifelse(t2-1.96*t2_se < 0, 0, t2-1.96*t2_se),
        "95% u.CI" = ifelse(t2+1.96*t2_se > 1, 1, t2+1.96*t2_se)
      )
    )
  }
  
  # chi squared test
  chi.sq <- list(
    "X^2" = sum(with(N_flatt,(frequency - expected_frequency)^2/expected_frequency)),
    "df"  = NK^NR - NK,
    "p.value" =   stats::pchisq(sum(with(N_flatt,(frequency - expected_frequency)^2/expected_frequency)),NK^NR - NK, lower.tail = F)
  )
  
  chi.sqB <- list(
    "X^2" = sum(with(N_flattB,(frequency - expected_frequency)^2/expected_frequency)),
    "df"  = factorial(NK+NR-1)/(factorial(NR)*factorial(NK-1)) - NK,
    "p.value" =   stats::pchisq(sum(with(N_flattB,(frequency - expected_frequency)^2/expected_frequency)),
                         factorial(NK+NR-1)/(factorial(NR)*factorial(NK-1)) - NK, lower.tail = F)
  )
  
  colnames(N_flatt)[1:NR] <- sapply(1:NR,function(x)paste(c("R",x),collapse = ""))
  colnames(N_flatt)[(NR+1):ncol(N_flatt)] <- c("Ob Freq","Exp Freq","Std Res")
  
  colnames(N_flattB)[1] <- "N Cat"
  colnames(N_flattB)[2:ncol(N_flattB)]    <- c("Ob Freq","Exp Freq","Std Res")
  
ret         <- list(
  "call"             = call,
  "data_description" = data_description,
  "allowed_classifications" = Q,
  "exact_prob"       = exact_prob,
  "model_A"          = list(
    "table"                   = N_flatt,
    "test_of_model_fit"       = chi.sq
  ),
  "model_B"          = list(
    "table"                   = N_flattB,
    "probabilities_estimates" = probabilities_estimates,
    "test_of_model_fit"       = chi.sqB,
    "test_of_symmetry"        = test_of_symmetry
  )
)
class(ret)  <- "boschek.test"
return(ret)
}

# methods
print.boschek.test         <- function(x){
  print_model_fit_A <- paste(c(
    'X^2(',x$model_A$test_of_model_fit$df,') = ',sprintf("%.2f", x$model_A$test_of_model_fit$`X^2`),', ',
    ifelse(x$model_A$test_of_model_fit$p.value < .001, 'p < .001', paste(c('p = ', sprintf("%.3f", x$model_A$test_of_model_fit$p.value)),collapse = "")),
    '\n'),sep = "",collapse = "")
  
  print_test_of_symmetry <- paste(c(
    'X^2(',x$model_B$test_of_symmetry$df,') = ',sprintf("%.2f", x$model_B$test_of_symmetry$`X^2`),', ',
    ifelse(x$model_B$test_of_symmetry$p.value < .001, 'p < .001', paste(c('p = ', sprintf("%.3f", x$model_B$test_of_symmetry$p.value)),collapse = "")),
    '\n'),sep = "",collapse = "")
  
  print_model_fit_B <- paste(c(
    'X^2(',x$model_B$test_of_model_fit$df,') = ',sprintf("%.2f", x$model_B$test_of_model_fit$`X^2`),', ',
    ifelse(x$model_B$test_of_model_fit$p.value < .001, 'p < .001', paste(c('p = ', sprintf("%.3f", x$model_B$test_of_model_fit$p.value)),collapse = "")),
    '\n'),sep = "",collapse = "")

  
  cat("Call:\n")
  print(x$call)
  
  cat("\nModel fit:\n")
  cat(print_model_fit_A)
  
  cat("\nModel fit assuming symetry:\n")
  cat(print_model_fit_B)
  
  if(x$data_description$number_of_categories == 2 & x$exact_prob){
    cat("\nEstimated probabilities:\n")
    print(x$model_B$probabilities_estimates[,1])
  }
  
  cat("\nTest of symmetry:\n")
  cat(print_test_of_symmetry)
}
summary.boschek.test       <- function(x){
  print_model_fit_A <- paste(c(
    'X^2(',x$model_A$test_of_model_fit$df,') = ',sprintf("%.2f", x$model_A$test_of_model_fit$`X^2`),', ',
    ifelse(x$model_A$test_of_model_fit$p.value < .001, 'p < .001', paste(c('p = ', sprintf("%.3f", x$model_A$test_of_model_fit$p.value)),collapse = "")),
    '\n'),sep = "",collapse = "")
  
  print_test_of_symmetry <- paste(c(
    'X^2(',x$model_B$test_of_symmetry$df,') = ',sprintf("%.2f", x$model_B$test_of_symmetry$`X^2`),', ',
    ifelse(x$model_B$test_of_symmetry$p.value < .001, 'p < .001', paste(c('p = ', sprintf("%.3f", x$model_B$test_of_symmetry$p.value)),collapse = "")),
    '\n'),sep = "",collapse = "")
  
  print_model_fit_B <- paste(c(
    'X^2(',x$model_B$test_of_model_fit$df,') = ',sprintf("%.2f", x$model_B$test_of_model_fit$`X^2`),', ',
    ifelse(x$model_B$test_of_model_fit$p.value < .001, 'p < .001', paste(c('p = ', sprintf("%.3f", x$model_B$test_of_model_fit$p.value)),collapse = "")),
    '\n'),sep = "",collapse = "")
  
  
  cat("Call:\n")
  print(x$call)
  
  cat(paste(c(
    '\n','Data description:','\n',
    'Observations: ',x$data_description$observations,'\n',
    'Categories:   ',x$data_description$number_of_categories,'\n',
    'Raters:       ',x$data_description$number_of_raters,'\n'
  ), collapse = "")
  )
  
  cat("\nModel fit:\n")
  cat(print_model_fit_A)
  
  cat("\nModel fit assuming symetry:\n")
  cat(print_model_fit_B)
  
  if(x$data_description$number_of_categories == 2 & x$exact_prob){
    cat("\nEstimated probabilities:\n")
    stats::printCoefmat(x$model_B$probabilities_estimates, digits = 2,
                 cs.ind  = c(1:4), tst.ind = integer(), zap.ind = integer())
  }
  
  cat("\nTest of symmetry:\n")
  cat(print_test_of_symmetry)
}
residuals.boschek.test     <- function(x, print.all = F){
  
  table_A <- x$model_A$table
  print_data_summary_add_A  <- ""
  if(!print.all){
    print_data_summary_rows <- 1:ifelse(nrow(table_A)<10,nrow(table_A),10)
    if(x$data_description$number_of_raters > 7){
      print_data_summary_add1 <- paste(c(' and ',ncol(table_A)-10,' columns'),collapse = "")
      table_A <- table_A[,c(1:8,(ncol(table_A)-2):ncol(table_A))]
      table_A[,8] <- "..."
      colnames(table_A)[8] <- "..."
    }else{
      print_data_summary_add1 <- NULL
    }
    if(nrow(table_A) - 10 > 0){
      print_data_summary_add_A  <- paste(c('... (shortened by ',nrow(table_A) - 10,' rows',print_data_summary_add1,')','\n'),collapse = "")
    }
    table_A <- table_A[print_data_summary_rows,]
  }
  
  
  table_B <- x$model_B$table
  print_data_summary_add_B <- ""
  if(!print.all){
    print_data_summary_rows <- 1:ifelse(nrow(table_B)<10,nrow(table_B),10)
    if(nrow(table_B) - 10 > 0){
      print_data_summary_add_B  <- paste(c('... (shortened by ',nrow(table_B) - 10,' rows)','\n'),collapse = "")
    }
    table_B <- table_B[print_data_summary_rows,]
  }
  
  cat("Model data matrix:\n")
  print(table_A)
  cat(print_data_summary_add_A)
  
  cat("\n")
  
  cat("Model data matrix assuming symetry:\n")
  print(table_B)
  cat(print_data_summary_add_B)
}

# data documentation
#' d1
#'
#' test dataset for ilustration
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{V1}{ratings by rater 1}
#'   \item{V2}{ratings by rater 2}
#'   \item{V3}{ratings by rater 3}
#'   \item{V4}{ratings by rater 4}
#' }
"d1"


# data documentation
#' d2
#'
#' test dataset for ilustration
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{V1}{rating by rater 1}
#'   \item{V2}{rating by rater 2}
#'   \item{V3}{rating by rater 3}
#'   \item{frequency}{number of rating with the particular combination of ratings}
#' }
"d2"