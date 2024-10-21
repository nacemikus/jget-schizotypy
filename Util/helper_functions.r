
logit <- function(x) log(x/(1-x))
inv_logit <- function(x) 1/(1+exp(-x))
pw = function(x,delta) x^delta/(x^delta + (1-x)^delta) 
sf <- function(x, pub = 0, prob_vect = c(0.5,0.025,0.975), dec_no = 3) { 
  y <- quantile(x, probs=prob_vect)%>% round(dec_no) 
  y[[4]] <- mean(x<0) %>% round(3)
  names(y)[4] <- "p"
  p_val = y[[4]];
  if (y[[4]] > 0.5) p_val = 1 - p_val
  
  if (y[[1]] < 0) {
    y_text = paste("b = ", y[[1]], ", 95% CrI [",y[[2]], ", ",y[[3]], "], P(b>0) = ",p_val, sep ="")
  } else   y_text = paste("b = ", y[[1]], ", 95% CrI [",y[[2]], ", ",y[[3]], "], P(b<0) = ",p_val, sep ="")
  
  if (pub == 0) {
    return(y)
  } else {
    
    return(y_text)  
  } }
wo <- function(x, d = 3, remove.na = FALSE) {
  if (remove.na) {
    x <- x[abs(x - mean(x, na.rm =TRUE)) <d*sd(x, na.rm =TRUE)]
  } else  {
    x[abs(x - mean(x, na.rm =TRUE)) > d*sd(x, na.rm =TRUE)] <- NA
    
  }
  return(x)
}
