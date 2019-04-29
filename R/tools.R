solve_T <- function(x, N_flatt, combinations, Q, NK, NR){
  y <- rep(0,length =  NK-1)
  
  P <- solve_P(x = x, combinations = combinations, Q = Q, NK = NK, NR = NR)
  
  for(k in 1:(NK-1)){
    y[k] <- sum(sapply(1:nrow(combinations),function(c){
      (N_flatt[c,"frequency"]/P[c]) * 
        (exp(sum(log(sapply(1:length(combinations[c,]),function(j)Q[k,combinations[c,j]])))) - exp(sum(log(sapply(1:length(combinations[c,]),function(j)Q[NK,combinations[c,j]])))))
      
    })) 
  }
  y
}
solve_P <- function(x, combinations, Q, NK, NR){
  P <- NULL
  for(i in 1:nrow(combinations)){
    comb <- combinations[i,]
    
    P <- c(P,
           sum(sapply(1:(NK-1),function(k){
             x[k]*(exp(sum(log(sapply(1:length(comb),function(j)Q[k,comb[j]])))) - exp(sum(log(sapply(1:length(comb),function(j)Q[NK,comb[j]]))))) 
           })) + exp(sum(log(sapply(1:length(comb),function(j)Q[NK,comb[j]]))))
    )
  }
  return(P)
}
solve_TQ <- function(x, N_flattB, NR){
  y <- numeric(2)
  
  A <- NULL
  B <- NULL
  C <- NULL
  
  for(i in 0:NR){
    A <- c(A,x[1]*(x[2]^i*(1-x[2])^(NR-i)-x[2]^(NR-i)*(1-x[2])^i)+x[2]^(NR-i)*(1-x[2])^i)
    B <- c(B,i*x[2]^(i-1)*(1-x[2])^(NR-1)-(NR-i)*x[2]^i*(1-x[2])^(NR-i-1))
    C <- c(C,(NR-i)*x[2]^(NR-i-1)*(1-x[2])^i-i*x[2]^(NR-i)*(1-x[2])^(i-1))
  }
  
  y[1] <- sum(sapply(0:NR,function(i){
    N_flattB[i+1,"frequency"]*(x[1]*(B[i+1]-C[i+1])+C[i+1])/A[i+1]
  }))
  y[2] <- sum(sapply(0:NR,function(i){
    N_flattB[i+1,"frequency"]*(x[2]^i*(1-x[2])^(NR-i)-x[2]^(NR-i)*(1-x[2])^i)/A[i+1]
  }))
  
  return(y)
}
compute_SE <- function(t, Q, N){
  se <- ((4*N*(Q-t))/(t-2*t*Q+Q^2) + (2*N*(1-2*Q)^2)/(Q-Q^2) + (4*N*(Q+t+1)^2)/((1-Q)^2-t+2*t*Q))^(-1/2)
}