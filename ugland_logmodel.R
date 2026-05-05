#' this fits the semilog model from Ugland et al
#' which takes a log-log regression
#' on the species frequencies table
#' and then uses expected counts
#' after hypergeometric sampling
#' to predict species
fit_slog_spec <- function(freq_table,
                          m,
                          cutoff=20,
                          show_lm=FALSE,
                          out="test1"){
#' Step 1: log-log fit of  species count
#' we take cutoff before the first 0
counts1 <- unname(freq_table[as.character(1:cutoff)])
counts1[is.na(counts1)] <- 0 
counts1 <- counts1[1:(which.min(counts1)-1)]
countlm <- lm(log(counts1) ~ log(seq(along=counts1))) 
if (show_lm){
pdf(paste0(out,".pdf"))  
plot(log(1:length(counts1)),log(counts1),
     main=paste("Spec abundance freqs,",out,
                "R2=",
                round(summary(countlm)$r.squared,
                      2)),
     xlab="log(spec count)",
     ylab="log(spec count abundance)")
abline(countlm,col="red",lty=2)
dev.off()}
cat("R2 loglog freqtable:",
    summary(countlm)$r.squared,"\n")
count_coeff <- unname(countlm$coefficients)
#' we now construct the semilog model
#' as in Ugland et al Appendix 2
#' with their notation
#' A is n in our paper, sample size of already observed 
#' S is species observed in sample (size A=n)
#' ksum is the denominator
est <- rep(-1,length(m))
names(est) <- m
A <- c(freq_table %*% as.numeric(names(freq_table)))
S <- sum(freq_table)
ksum <- sum((1:A)^count_coeff[2])
#' We need to approx L(a) from Appendix 2 of U. et al.
#' a is subsample sizes 
a <- c(seq(1,A,5),A)
L2 <- function(a1){
  L1 <- function(k){prod(1 - (a1/(A:(A-(k-1)))))*k^(count_coeff[2])}
  out1 <- sum(sapply(1:(A-a1),L1))}
Llm_data <- data.frame(a,L=sapply(a,L2)) 
logLlm <- lm(L~log(a),data=Llm_data)
cat("R2 semilog L:",summary(logLlm)$r.squared,"\n")
Llm_coeff <- unname(logLlm$coefficients)
for (m1 in m){
a <- A + m1
est[as.character(m1)] <- S*(1-ksum^(-1)*(Llm_coeff[2]*log(a)+Llm_coeff[1]))
}
return(est)
}
