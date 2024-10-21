ICC3k <-  function(x) 
{
  
  if (is.matrix(x)) 
    x <- data.frame(x)
  items <- colnames(x)
  n.items <- NCOL(x)
  
  
  n.obs <- dim(x)[1]
  nj <- dim(x)[2]
  x.s <- stack(x)
  x.df <- data.frame(x.s, subs = rep(paste("S", 1:n.obs, sep = ""), 
                                     nj))
  
  
  
  colnames(x.df) <- c("values", "items", "id")
  # mod.lmer <- lme4::lmer(values ~ 1 + (1 | id) + (1 | items), 
  #                        data = x.df, na.action = na.omit)
  mod.lmer <- lme4::lmer(values ~ 1 + (1 | id) + items, 
                         data = x.df, na.action = na.omit)
  
  vc <- lme4::VarCorr(mod.lmer)
  MS_id <- vc$id[1, 1]
  MS_items <- vc$items[1, 1]
  MSE <- error <- MS_resid <- (attributes(vc)$sc)^2 
  MSB <- nj * MS_id + error
  ICC_coef <- (MSB - MSE)/MSB
  return(ICC_coef)
}

