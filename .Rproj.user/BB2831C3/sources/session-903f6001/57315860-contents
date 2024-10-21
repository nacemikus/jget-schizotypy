# 
rm(list = ls())

## load packages  --------------------------------------------------------------------
library(tidyverse)
library(PMA)
#library(gdata)
#library(lme4)
# library(rmatio)
library(ggridges)
library(wesanderson)
# library(smooth)
# library(ggplot2)
# library(reshape2)
library(psych)
# library(gtools)
# library(plyr)
# library(hrbrthemes)
# library(dplyr)
# library(tidyr)
# library(viridis)
library(gridExtra)
library(cowplot)
#library(psych)
#library(tidyr)
#library(ggpubr)
#library(grid)
# library(ggthemes)
#library(ggforce)
#library(stargazer)
#library(tables)
# library(knitr)
#library(xtable)
#library(lmerTest)
#library(ez)
#library(gtools)
#library(perm)
#library(caret)
#library(pracma)
#library(ggridges)
# library(nlme)
library(lme4)
library(lmerTest)
library(brms)
# library(MASS)
# library(tidyverse)  # ggplot, dplyr, and friends
# library(ggridges)   # Ridge plots
# library(ggstance)   # Horizontal pointranges and bars
# library(patchwork)  # Lay out multiple ggplot plots; install from https://github.com/thomasp85/patchwork
# library(scales)     # Nicer formatting for numbers
# library(broom)
# library(loo)
library(rstan)

library(bayesplot)
# library(rstanarm)
## load helper functions ####
source("Util/theme_functions_trt.r")
source("Util/ICC3k.r")
source("Util/helper_functions.r")
# pdi_colour = c("#005828", "#C5E0B4") #
# c( "#CCE0CC", "#006000")
pdi_colour = c(  "#4472C4", "grey50") 
############################## #
## Load the TRT data #####
############################## #
data_trial <- readRDS(file = "Data/data_trial_preprocessed.rds")
data_all <- readRDS(file = "Data/data_sub_preprocessed.rds")

data_all_sess <- data_all %>% group_by(ID) %>% summarize(ses_no = length(unique(session)))
data_all <- merge(data_all, data_all_sess, by = "ID")


## save lables for plots #### 
my_labels = c( expression(paste("Tonic belief volatility (\U1D714)")),
               expression(paste("Variance belief volatility (", "\U1D714"["\U1D6FC"], ")   ")) ,
               expression("Environmental volatility (\U1D717)"))

my_greeks = c( "\U1D714",
               expression(paste("\U1D714"["\U1D6FC"])) ,
               "\U1D717",
               expression("\U1D707"["\U1D6FC"]) ,
               expression("\U1D707"["\U1D463"]) ,
               expression("\U1D713"["1"]) ,
               "\U1D707",
               "\U1D6FC",
               "\U1D713")


Bootstrapped_summary_new <- readRDS("Data/ICC_data_plot.rds")
x_txt_size <- 20
txt_size = 15
g_trt <- 
  Bootstrapped_summary_new %>% filter(ntrials == 241, par_name != "noise", analysis == "pooling") %>% 
  ggplot(aes(x = par_name, y = est)) + 
 
  geom_errorbar(aes( ymin = Q025, ymax =Q975), size = 2, width= 0)+
  geom_point(aes(x = par_name, y = est),size = 5, shape = 21, colour = "black", fill =  "#C27D38")  + ylim( c(0,1))+
  # geom_boxplot(aes(x = par_name, lower = Q25, upper = Q75, middle = est, ymin = Q025, ymax =Q975),
  #              stat = "identity", position= position_dodge(0.5), fill =  "#C27D38", width = 0.5)  +
  # 
  theme_Publication(base_size = txt_size) + 
  ylab("ICC (means with 95% CI)")  +  xlab("") + # + xlab("HGF, Belief Volatility, \U03C9") , fill =  "#C27D38"  expression("\U03C9"[a]) xlab("\U1D717")
  theme(axis.line = element_line(size = 1.5),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 20)) + #facet_wrap(~analysis) 
  # scale_fill_manual(values =   c("#C27D38", "#CEAB07"))+
  scale_x_discrete(labels = c(my_greeks[[1]], my_greeks[[3]], my_greeks[[2]] ))#c("\U1D714", "\U1D717", expression("\U1D714"["\U1D6FC"])))
g_trt 

g_trt <- 
  Bootstrapped_summary_new %>% filter(ntrials == 241, par_name != "noise", analysis == "pooling") %>% 
  ggplot(aes(x = par_name, y = est)) + 
  
  geom_errorbar(aes( ymin = Q025, ymax =Q975), size = 2, width= 0)+
  geom_point(aes(x = par_name, y = est),size = 5, shape = 21, colour = "black", fill =  "#C27D38")  + ylim( c(0,1))+
  # geom_boxplot(aes(x = par_name, lower = Q25, upper = Q75, middle = est, ymin = Q025, ymax =Q975),
  #              stat = "identity", position= position_dodge(0.5), fill =  "#C27D38", width = 0.5)  +
  # 
  theme_Publication(base_size = txt_size) + 
  ylab("ICC (means with 95% CI)")  +  xlab("") + # + xlab("HGF, Belief Volatility, \U03C9") , fill =  "#C27D38"  expression("\U03C9"[a]) xlab("\U1D717")
  theme(axis.line = element_line(size = 1.5),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 20)) + #facet_wrap(~analysis) 
  # scale_fill_manual(values =   c("#C27D38", "#CEAB07"))+
  scale_x_discrete(limits = c("om", "om2", "oma"), labels = c(my_greeks[[1]], my_greeks[[3]], my_greeks[[2]] ))#c("\U1D714", "\U1D717", expression("\U1D714"["\U1D6FC"])))
g_trt 
# trt from PE to conf ######
sub2use = data_trial %>% filter(!is.na(pdi_t), session == 4  ) %>% pull(ID) %>% unique
sess_no = 4


data_trial$signed_update_s =  data_trial$signed_update %>% ave(FUN = scale)
data_trial$signed_update_next_trial_s  = data_trial$signed_update_s %>% lead()
data_trial$signed_PE_s = data_trial$signed_PE %>% ave(FUN = scale)
data_trial = data_trial %>% mutate( log_lr = ((data_trial$update+1)/(data_trial$PE+1)) %>% log %>% ave(FUN = scale))

betas_conf_PE <- matrix(NA, ncol =sess_no, nrow = 45)
# betas_conf_sd <- matrix(NA, ncol =sess_no, nrow = 45)
betas_update_PE<- matrix(NA, ncol =sess_no, nrow = 45)
betas_conf_update<- matrix(NA, ncol =sess_no, nrow = 45)
betas_lr_sd <- matrix(NA, ncol =sess_no, nrow = 45)
betas_conf_sd <- matrix(NA, ncol =sess_no, nrow = 45)
betas_conf_PE2<- matrix(NA, ncol =sess_no, nrow = 45)
for (s in c(1:sess_no)) { 
  print(s)
  model1 <- lmer(data = data_trial %>% filter(ID %in% sub2use, session == s  ), confidence_next_trial ~  (PE_t)  + (PE_t| ID))
  
  model2 <- lmer(data = data_trial %>% filter(ID %in% sub2use, session == s  ), signed_update_next_trial_s ~ (signed_PE_s)  + (signed_PE_s| ID))
  model3 <- lmer(data = data_trial %>% filter(ID %in% sub2use, session == s  ), log_lr ~ (sd)  + (sd| ID))
  
  
  betas_conf_PE[,s] = coef(model1)$ID[,2]
  betas_update_PE[,s] = coef(model2)$ID[,2]
  betas_lr_sd[,s] = coef(model3)$ID[,2]
}




betas_conf_PE %>% ICC()
betas_conf_PE %>% cor

betas_update_PE %>% ICC()
betas_lr_sd %>% ICC()
betas_lr_sd %>% cor

# filter all subjects
data_all2 = data_all %>% filter(ID %in% sub2use)

data_all2$betas_conf_PE = as.vector(t(betas_conf_PE) )
data_all2$betas_update_PE = as.vector(t(betas_update_PE) )
data_all2$betas_conf_update = as.vector(t(betas_conf_update) )  
data_all2$betas_lr_sd = as.vector(t(betas_lr_sd) )  

## model checks (fig 3 a - c) ####

get_total_d = function(model) { 
  get_sum = summary(model)
  total_sig = sqrt(get_sum$varcor$ID %>% diag() %>% sum + get_sum$sigma^2)
  fixef_temp = model %>% fixef() /total_sig
  print(fixef_temp)
}

model_mua <- lmer(data = data_trial , mua_h ~ (sd)  + (sd| ID))
model_mua %>% summary()
get_total_d(model_mua)

model_mux2 <- lmer(data = data_trial , mux2_h ~ (sd)  + (sd| ID))
model_mux2 %>% summary()
get_total_d(model_mux2)


model_sd <- lmer(data = data_trial , log_lr ~ (sd)  + (sd| ID))
model_sd %>% summary()
get_total_d(model_sd)

# predicting the behavioral pattern from parameters  ####
# scale parameters
data_all2$om_h_s = data_all2$om_h %>% ave(FUN = scale)
data_all2$oma_h_s = data_all2$oma_h %>% ave(FUN = scale)
data_all2$th_h_s = data_all2$th_h %>% ave(FUN = scale)



lm(data = data_all2 ,betas_lr_sd %>% ave(FUN = scale, centre  =F )  ~   om_h_s + th_h_s + oma_h_s) %>% summary() 

lm(data = data_all2 ,betas_conf_PE %>% ave(FUN = scale, centre  =F )  ~  om_h_s + th_h_s + oma_h_s) %>% summary() 

lm(data = data_all2 ,betas_update_PE %>% ave(FUN = scale, centre  =F )  ~  om_h_s + th_h_s + oma_h_s) %>% summary() 


g_om_update_PE = data_all2  %>% ggplot(aes(x = om_h, y = betas_update_PE)) + geom_point(alpha = 0.5) + geom_smooth(method = "lm", se = FALSE, colour = "black") + theme_Publication() +
  ylab(expression(paste("Update ~ PE (\U1D6A9"["U"], ")"))) + xlab(my_greeks[[1]])+ theme(axis.title.y = element_text(size = 15)) + ylim(c(0.35,1.2))
g_oma_lr_sd =   data_all2 %>% ggplot(aes(x = oma_h, y = betas_lr_sd)) + geom_point(alpha = 0.5) + geom_smooth(method = "lm", se = FALSE, colour = "black") + theme_Publication() +
  ylab(expression(paste("Learning Rate ~ Variance (\U1D6A9"["L"], ")")))  + xlab(my_greeks[[2]])+ theme(axis.title.y = element_text(size = 15))
g_th_conf_PE =  data_all2  %>% ggplot(aes(x = th_h, y = betas_conf_PE)) + geom_point(alpha = 0.5) + geom_smooth(method = "lm", se = FALSE, colour = "black") + theme_Publication() +
  ylab(expression(paste("Confidence ~ PE (\U1D6A9"["C"], ")")))   + xlab(my_greeks[[3]]) + theme(axis.title.y = element_text(size = 15))
g_2step_task =  data_all %>% filter(!is.na(ts_score), session == 1) %>% ggplot(aes(x=th_h, y = ts_score)) + geom_point(alpha = 0.5) + geom_smooth(method = "lm", colour = "black", se = F) + theme_Publication()+
  xlab(my_greeks[[3]]) + ylab("Points Earned") + labs(subtitle = "Two-step task")

g_betas_oms = plot_grid(g_om_update_PE, g_oma_lr_sd, g_th_conf_PE, g_2step_task, nrow = 1)
g_betas_oms


# plot icc lrs ####
get_me_those_ICC= function(par_name, data = data_all2) {
data_lr_h_trt <- data_all2  %>% select(ID, !!sym(par_name), session) %>% spread(session, !!sym(par_name)) 
data_lr_h_trt <- data_lr_h_trt %>% ungroup %>% select(!ID) 
colnames(data_lr_h_trt) = c("session1", "session2", "session3", "session4")
return(data_lr_h_trt)
}

get_me_those_ICC("mean_mua")
get_me_those_ICC("mean_mux2")
get_me_those_ICC("mean_loglr1")
get_me_those_ICC("mean_loglr2_h")

data_trial$index = data_trial$ID%>% as.numeric() + 60*(data_trial$session%>% as.numeric() -1)
data_trial2 = data_trial %>% filter(!is.na(pdi_total_cat))


plot_noise_data = function(yvar = "lr", data_temp =data_trial2) {
  g <- data_temp %>% group_by(ID, sd) %>%   summarize(
    yvar_mean_id = mean(!!sym(yvar), na.rm  = T)) %>% group_by(sd) %>% 
    summarize(N=n(),
              mean_var = mean(yvar_mean_id),
              se_var = sd(yvar_mean_id)/sqrt(N-1)) %>% 
    ggplot(aes(x=sd, y = mean_var)) + 
    geom_errorbar(aes(ymin= mean_var - se_var, ymax = mean_var + se_var), size = 1, colour = "black", width = 0.2, position = position_dodge(width = 0.5))+
    geom_bar(stat = "identity", width = 0.5 , position = position_dodge2(), colour = "black")  + theme_Publication() +
    theme(axis.line = element_line(size = 1.5),
          panel.grid = element_blank(),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15)) + 
    scale_fill_manual(values =   rev(pdi_colour))+
    scale_colour_manual(values = rev(pdi_colour))+
    scale_x_discrete(labels = c("Low", "High"))+
    xlab("Variance") + ylab("Unsigned PE (log)") 
  return(g)
}


## iccs mus ####
# data_trt =rbind(get_me_those_ICC("mean_mua_h"), get_me_those_ICC("mean_mux2_h"))
# get_me_those_ICC("mean_loglr1_h")

icc_data = tibble(ICC_mean = ICC(get_me_those_ICC("mean_loglr1_h"))$results$ICC[6],
       ICC_lci = ICC(get_me_those_ICC("mean_loglr1_h"))$results[6,7],
       ICC_uci = ICC(get_me_those_ICC("mean_loglr1_h"))$results[6,8],
       par = "mean_loglr1_h")
icc_data_temp =  tibble(ICC_mean = ICC(get_me_those_ICC("mean_mua_h"))$results$ICC[6],
                        ICC_lci = ICC(get_me_those_ICC("mean_mua_h"))$results[6,7],
                        ICC_uci = ICC(get_me_those_ICC("mean_mua_h"))$results[6,8],
                        par = "mean_mua_h")

icc_data = rbind(icc_data, icc_data_temp)
icc_data_temp =  tibble(ICC_mean = ICC(get_me_those_ICC("mean_mux2_h"))$results$ICC[6],
                        ICC_lci = ICC(get_me_those_ICC("mean_mux2_h"))$results[6,7],
                        ICC_uci = ICC(get_me_those_ICC("mean_mux2_h"))$results[6,8],
                        par = "mean_mux2_h")

icc_data_mus = rbind(icc_data, icc_data_temp)

#ICC3

icc_data = tibble(ICC_mean = ICC(get_me_those_ICC("mean_loglr1_h"))$results$ICC[3],
                  ICC_lci = ICC(get_me_those_ICC("mean_loglr1_h"))$results[3,7],
                  ICC_uci = ICC(get_me_those_ICC("mean_loglr1_h"))$results[3,8],
                  par = "mean_loglr1_h")
icc_data_temp =  tibble(ICC_mean = ICC(get_me_those_ICC("mean_mua_h"))$results$ICC[3],
                        ICC_lci = ICC(get_me_those_ICC("mean_mua_h"))$results[3,7],
                        ICC_uci = ICC(get_me_those_ICC("mean_mua_h"))$results[3,8],
                        par = "mean_mua_h")

icc_data = rbind(icc_data, icc_data_temp)
icc_data_temp =  tibble(ICC_mean = ICC(get_me_those_ICC("mean_mux2_h"))$results$ICC[3],
                        ICC_lci = ICC(get_me_those_ICC("mean_mux2_h"))$results[3,7],
                        ICC_uci = ICC(get_me_those_ICC("mean_mux2_h"))$results[3,8],
                        par = "mean_mux2_h")

icc_data_mus_ICC3 = rbind(icc_data, icc_data_temp)

x_txt_size <- 20
txt_size = 15
dodge_width = 0.9
g_trt_mus <- 
  icc_data_mus %>% 
  ggplot(aes(x = par, y = ICC_mean)) + 
  
  geom_errorbar(aes( ymin = ICC_lci, ymax =ICC_uci), size = 2, width= 0 )+
  geom_point(size = 5, shape = 21, colour = "black", fill =  "#C27D38" )  + ylim( c(0,1))+
  geom_errorbar(data = icc_data_mus_ICC3, aes(ymin = ICC_lci, ymax =ICC_uci), size = 2, width= 0 )+
  geom_point(data = icc_data_mus_ICC3,size = 5, shape = 21, colour = "black", fill =  "#CEAB07" )  + ylim( c(0,1))+
  # geom_boxplot(aes(x = par_name, lower = Q25, upper = Q75, middle = est, ymin = Q025, ymax =Q975),
  #              stat = "identity", position= position_dodge(0.5), fill =  "#C27D38", width = 0.5)  +
  # 
  theme_Publication(base_size = 15) + 
  ylab("ICC (means with 95% CI)")  +  xlab("") + # + xlab("HGF, Belief Volatility, \U03C9") , fill =  "#C27D38"  expression("\U03C9"[a]) xlab("\U1D717")
  theme(axis.line = element_line(size = 1.5),
        panel.grid = element_blank(),
        legend.position = "none",  
        axis.text.x = element_text(size = 15,angle = 0))+ #+ #facet_wrap(~analysis) 
  # scale_fill_manual(values =   c("#C27D38", "#CEAB07"))+ PH1
  scale_x_discrete(labels = c(my_greeks[[6]],
                              my_greeks[[4]],
                              my_greeks[[5]]))
  # scale_x_discrete(labels = c(expression(paste("Precision-weight \U1D713"[1])),
  #                             expression(paste("Mean belief about variance")),
  #                             expression(paste("Mean belief about volatility"))))
g_trt_mus 


# plot icc data ####
#ICC 3k
icc_data = tibble(ICC_mean = ICC(betas_lr_sd)$results$ICC[6],
                  ICC_lci = ICC(betas_lr_sd)$results[6,7],
                  ICC_uci = ICC(betas_lr_sd)$results[6,8],
                  par = "betas_conf_PE")
icc_data_temp =  tibble(ICC_mean = ICC(betas_update_PE)$results$ICC[6],
                        ICC_lci = ICC(betas_update_PE)$results[6,7],
                        ICC_uci = ICC(betas_update_PE)$results[6,8],
                        par = "betas_update_PE")

icc_data = rbind(icc_data, icc_data_temp)
icc_data_temp =  tibble(ICC_mean = ICC(betas_conf_PE)$results$ICC[6],
                        ICC_lci = ICC(betas_conf_PE)$results[6,7],
                        ICC_uci = ICC(betas_conf_PE)$results[6,8],
                        par = "betas_lr_sd")

icc_data_trt = rbind(icc_data, icc_data_temp)
# ICC 3
icc_data = tibble(ICC_mean = ICC(betas_lr_sd)$results$ICC[3],
                  ICC_lci = ICC(betas_lr_sd)$results[3,7],
                  ICC_uci = ICC(betas_lr_sd)$results[3,8],
                  par = "betas_conf_PE")
icc_data_temp =  tibble(ICC_mean = ICC(betas_update_PE)$results$ICC[3],
                        ICC_lci = ICC(betas_update_PE)$results[3,7],
                        ICC_uci = ICC(betas_update_PE)$results[3,8],
                        par = "betas_update_PE")

icc_data = rbind(icc_data, icc_data_temp)
icc_data_temp =  tibble(ICC_mean = ICC(betas_conf_PE)$results$ICC[3],
                        ICC_lci = ICC(betas_conf_PE)$results[3,7],
                        ICC_uci = ICC(betas_conf_PE)$results[3,8],
                        par = "betas_lr_sd")

icc_data_trt_ICC3 = rbind(icc_data, icc_data_temp)

x_txt_size <- 20
txt_size = 15

g_trt_beh <- 
  icc_data_trt %>% 
  ggplot(aes(x = par, y = ICC_mean)) + 
  
  geom_errorbar(aes( ymin = ICC_lci, ymax =ICC_uci), size = 2, width= 0)+
  geom_point(size = 5, shape = 21, colour = "black", fill =  "#C27D38")  + ylim( c(0,1))+
  geom_errorbar(data = icc_data_trt_ICC3, aes(ymin = ICC_lci, ymax =ICC_uci), size = 2, width= 0)+
  geom_point(data = icc_data_trt_ICC3,size = 5, shape = 21, colour = "black", fill =  "#CEAB07")  + ylim( c(0,1))+
  # geom_boxplot(aes(x = par_name, lower = Q25, upper = Q75, middle = est, ymin = Q025, ymax =Q975),
  #              stat = "identity", position= position_dodge(0.5), fill =  "#C27D38", width = 0.5)  +
  # 
  theme_Publication(base_size = 15) + 
  ylab("ICC (means with 95% CI)")  +  xlab("") + # + xlab("HGF, Belief Volatility, \U03C9") , fill =  "#C27D38"  expression("\U03C9"[a]) xlab("\U1D717")
  theme(axis.line = element_line(size = 1.5),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 15,angle = 0))+ #+ #facet_wrap(~analysis) 
  # scale_fill_manual(values =   c("#C27D38", "#CEAB07"))+
  scale_x_discrete(labels = c(expression(paste("\U1D6A9"["C"])),
                              expression(paste("\U1D6A9"["U"])),
                              expression(paste("\U1D6A9"["L"]))))
g_trt_beh

## plot Figure 3 ####
g_fig3_part1 =   plot_grid(plot_noise_data("log_lr", data_trial) + ylab("Learning rate (log)"),
                    plot_noise_data("mua_h_sc",  data_trial)  + ylab(paste("Belief about variance")),
                    plot_noise_data("mux2_h_sc",  data_trial) + ylab(paste("Belief about volatility")),
                    g_trt, g_trt_mus,g_trt_beh,
                    nrow = 2, align = "v")                
g_fig3 = plot_grid(
  g_fig3_part1,
  g_betas_oms, ncol=1, rel_heights = c(1,0.5))
g_fig3

ggsave("g_fig3.png", plot = g_fig3, width = 10, height = 12, dpi = 300)

## associations of hgf parameters with pdi ####

model_name= "hgfxPDI.rds"
if (!file.exists(paste("Brms_models/", model_name,sep = ""))) {
  print(model_name)
  data_sum <- data_all %>% group_by(ID) %>% summarize(Mean_PE = mean(mean_PE, na.rm = TRUE), 
                                                      pdi_t = pdi_t[1], 
                                                      psc_total_s = psc_total_s[1],
                                                      pdi_total = pdi_total[1], 
                                                      mean_om_h = mean(om_h, na.rm = TRUE),
                                                      mean_om2_h = mean(th_h, na.rm = TRUE),
                                                      mean_oma_h = mean(oma_h, na.rm = TRUE),
                                                      mean_om = mean(om, na.rm = TRUE),
                                                      mean_om2 = mean(om2, na.rm = TRUE),
                                                      mean_oma = mean(oma, na.rm = TRUE)) %>% 
    mutate( mean_om_h_sc = ave(mean_om_h, FUN = scale),
            mean_om2_h_sc = ave(mean_om2_h, FUN = scale),
            mean_oma_h_sc = ave(mean_oma_h, FUN = scale))
  
  bay.model <- brm(pdi_t ~  mean_om_h_sc + mean_om2_h_sc + mean_oma_h_sc,
                   data = data_sum,
                   prior =c(set_prior("normal(0,1)", class = "b")),
                   warmup = 500, iter =2000, chains =4)
  saveRDS(bay.model,file =  paste("Brms_models/", model_name,sep = ""))
}
hgfxPDI <- readRDS("Brms_models/hgfxPDI.rds")
hgfxPDI %>% summary()

post_scz <- hgfxPDI  %>% posterior_samples()
post_scz <- post_scz/post_scz$sigma
(post_scz$b_mean_om2_h_sc) %>% sf(1)
(post_scz$b_mean_om_h_sc) %>% sf(1)
(post_scz$b_mean_oma_h_sc) %>% sf(1)

d_om2 <- (post_scz$b_mean_om2_h_sc)
d_om2 <- d_om2 %>% as.data.frame()  %>%  `colnames<-`("om2")   %>% mutate(type = "NA", type_name = NA)

d_om <- (post_scz$b_mean_om_h_sc)
d_om <- d_om %>% as.data.frame()  %>%  `colnames<-`("om")   %>% mutate(type = "NA", type_name = NA)

d_oma <- (post_scz$b_mean_oma_h_sc)
d_oma <- d_oma %>% as.data.frame()  %>%  `colnames<-`("oma")   %>% mutate(type = "NA", type_name = NA)


# d_om2 %>% glimpse()
d_om2$par_name = "om2"
d_om$par_name = "om"
d_oma$par_name = "oma"
colnames(d_om2)[1] <- "par"
colnames(d_om)[1] <- "par"
colnames(d_oma)[1] <- "par"

d_all <- rbind(d_om, d_om2, d_oma)

d_int <- d_all %>% group_by(par_name, type) %>% summarise(mean =quantile(par, prob = .5),
                                                          ll   = quantile(par, prob = .025),
                                                          ul   = quantile(par, prob = .975),
                                                          lls   = quantile(par, prob = .10),
                                                          uls   = quantile(par, prob = .90))   # since the `key` variable is really two variables in one, here we split them up
# pdi plot (Figure 4a) ####
g_pdi <- ggplot(data = d_int) +
  geom_density_ridges(data = d_all %>% filter(type == "NA"), aes(x = par, y = as.factor(par_name),fill =  as.factor(par_name), height = ..density..), alpha = 0.5, scale = 0.8) + 
  geom_errorbar(data = d_int %>% filter(type == "NA"), aes(x = mean, y = par_name, color = type, group = type,xmin = lls, xmax = uls), size = 1.5, width = 0, position = position_dodge(0.2))+
  geom_vline(xintercept = 0,  linetype = "dashed") +
  geom_pointrange(data = d_int %>% filter(type == "NA"), aes(x = mean, xmin = ll, xmax = ul, y = par_name, color = type, group = type), position = position_dodge(0.2)) + # c("firebrick", "black", "red")
  xlab("Effect sizes (95% and 80% CrI)")+
  labs(title = "Association with PDI")+ 
  theme_Publication(base_size = 15) +
  theme(legend.position = "none",
        axis.line = element_line(size = 1.5),
        # 
        # panel.grid = element_blank(),
        # axis.title.y = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=20),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "transparent"))  +
  scale_color_manual(values  = c("firebrick")) + 
  scale_fill_manual(values  = c("#CCCCCC", "#4472C4", "#CCCCCC")) +
  scale_y_discrete(limits = c("oma", "om2", "om"), 
                   labels = c(expression("\U1D714"["\U1D6FC"]), "\U1D717", "\U1D714")) #labels = c(expression(paste("Noise volatility (", "\U03C9"[a], ")")),"Environmental Volatility (\U1D717)",  "Belief Volatility (\U03C9)")  )
g_pdi

## plot hgf params against PDI
  data_sum <- data_all %>% group_by(ID) %>% summarize(Mean_PE = mean(mean_PE, na.rm = TRUE), 
                                                      pdi_t = pdi_t[1], 
                                                      psc_total_s = psc_total_s[1],
                                                      pdi_total = pdi_total[1], 
                                                      mean_om_h = mean(om_h, na.rm = TRUE),
                                                      mean_om2_h = mean(th_h, na.rm = TRUE),
                                                      mean_oma_h = mean(oma_h, na.rm = TRUE),
                                                      mean_noise_h = mean(noise_h, na.rm = TRUE),
                                                      mean_om = mean(om, na.rm = TRUE),
                                                      mean_om2 = mean(om2, na.rm = TRUE),
                                                      mean_oma = mean(oma, na.rm = TRUE),
                                                      mean_noise = mean(noise, na.rm = TRUE),
                                                      ts_score = ts_score[1]) 
  
  
 
  
  g_om2_mean <- data_sum  %>% 
    ggplot(aes(x = pdi_t, y = mean_om2_h)) + 
    geom_point()+
    stat_smooth(method = "lm", se = F, colour = "black") + 
    # geom_ribbon(aes(ymin = par_lb, ymax = par_up), alpha = 0.5, fill = "grey70") +  
    
    theme_Publication(base_size = 15) +
    # coord_cartesian(ylim = c(-10,-3.5))+
    # annotate("text", x= 0.4, y = 21, label = "\u03B2 = 0.33, 95% CI [0.14, 0.53]", size = 8)+
    theme(legend.position = "none",
          axis.line = element_line(size = 1.5),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5)) + xlab("PDI (scaled)") + ylab("Environmental volatility (\U1D717)")
  
  
  g_om2_mean 
  
  g_oma_mean <- data_sum  %>% 
    ggplot(aes(x = pdi_t, y = mean_oma_h)) + 
    geom_point()+
    stat_smooth(method = "lm", se = F, colour = "black") + 
    # geom_ribbon(aes(ymin = par_lb, ymax = par_up), alpha = 0.5, fill = "grey70") +  
    
    theme_Publication(base_size = 15) +
    coord_cartesian(ylim = c(-10,-3.5))+
    # annotate("text", x= 0.4, y = 21, label = "\u03B2 = 0.33, 95% CI [0.14, 0.53]", size = 8)+
    theme(legend.position = "none",
          axis.line = element_line(size = 1.5),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5)) + xlab("PDI (scaled)") + ylab("Environmental volatility (\U1D717)")
  
  
  g_oma_mean 
  
  g_om_mean <- data_sum  %>% 
    ggplot(aes(x = pdi_t, y = mean_om_h)) + 
    geom_point()+
    stat_smooth(method = "lm", se = F, colour = "black") + 
    # geom_ribbon(aes(ymin = par_lb, ymax = par_up), alpha = 0.5, fill = "grey70") +  
    
    theme_Publication(base_size = 15) +
    # coord_cartesian(ylim = c(-10,-3.5))+
    # annotate("text", x= 0.4, y = 21, label = "\u03B2 = 0.33, 95% CI [0.14, 0.53]", size = 8)+
    theme(legend.position = "none",
          axis.line = element_line(size = 1.5),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5)) + xlab("PDI (scaled)") + ylab("Environmental volatility (\U1D717)")
  
  
  g_om_mean 
  
 
## similar results obtained when using models with parameters as DVs and including random effects per subject per session:
  if (FALSE) {
    model_name= "om2_hXpdiTXsess.rds"
    
    if (!file.exists(paste("Brms_models/", model_name,sep = ""))) {
      print(model_name)
      
      bay.model <- brm(th_h ~  pdi_t * session + (session|ID),
                       data = data_all,
                       prior =c(set_prior("normal(0,1)", class = "sd"),
                                set_prior("normal(0,1)", class = "b"),
                                set_prior("lkj(2)", class = "cor")),
                       warmup = 1000, iter =3000, chains =4,
                       control = list(adapt_delta = 0.95))
      saveRDS(bay.model,file =  paste("Brms_models/", model_name,sep = ""))
    }
   
    beh_model_th <- readRDS("Brms_models/om2_hXpdiTXsess.rds")
    beh_model_th %>% summary()
    post_th = beh_model_th %>% as_draws_df()
    post_th = post_th %>% mutate(sd_total = (sigma^2 +sd_ID__Intercept^2 + sd_ID__session2^2 + sd_ID__session3^2 + sd_ID__session4^2) %>% sqrt)
    (post_th$b_pdi_t) %>% sf(1)
    (post_th$b_pdi_t+ post_th$`b_pdi_t:session2`) %>% sf(1)
    (post_th$b_pdi_t+ post_th$`b_pdi_t:session3`) %>% sf(1)
    (post_th$b_pdi_t+ post_th$`b_pdi_t:session4`) %>% sf(1)
      
  }
  if (FALSE) {
    
    model_name= "oma_hXpdiTXsess.rds"
    
    if (!file.exists(paste("Brms_models/", model_name,sep = ""))) {
      print(model_name)
      
      bay.model <- brm(oma_h ~  pdi_t * session + (session|ID),
                       data = data_all,
                       prior =c(set_prior("normal(0,1)", class = "sd"),
                                set_prior("normal(0,1)", class = "b"),
                                set_prior("lkj(2)", class = "cor")),
                       warmup = 1000, iter =3000, chains =4,
                       control = list(adapt_delta = 0.95))
      saveRDS(bay.model,file =  paste("Brms_models/", model_name,sep = ""))
    }
    beh_model_oma <- readRDS("Brms_models/oma_hXpdiTXsess.rds")
    beh_model_oma %>% summary()
    post_oma = beh_model_oma %>% as_draws_df()
    post_oma = post_oma %>% mutate(sd_total = (sigma^2 +sd_ID__Intercept^2 + sd_ID__session2^2 + sd_ID__session3^2 + sd_ID__session4^2) %>% sqrt)
    (post_oma$b_pdi_t) %>% sf(1)
    (post_oma$b_pdi_t+ post_oma$`b_pdi_t:session2`) %>% sf(1)
    (post_oma$b_pdi_t+ post_oma$`b_pdi_t:session3`) %>% sf(1)
    (post_oma$b_pdi_t+ post_oma$`b_pdi_t:session4`) %>% sf(1)
    
  }
  if (FALSE) {
    model_name= "om_hXpdiTXsess.rds"
    
    if (!file.exists(paste("Brms_models/", model_name,sep = ""))) {
      print(model_name)
      
      bay.model <- brm(om_h ~  pdi_t * session + (session|ID),
                       data = data_all,
                       prior =c(set_prior("normal(0,1)", class = "sd"),
                                set_prior("normal(0,1)", class = "b"),
                                set_prior("lkj(2)", class = "cor")),
                       warmup = 1000, iter =3000, chains =4,
                       control = list(adapt_delta = 0.95))
      saveRDS(bay.model,file =  paste("Brms_models/", model_name,sep = ""))
    }
    
    beh_model_om <- readRDS("Brms_models/om_hXpdiTXsess.rds")
    beh_model_om %>% summary()
    post_om = beh_model_om %>% as_draws_df()
    post_om = post_om %>% mutate(sd_total = (sigma^2 +sd_ID__Intercept^2 + sd_ID__session2^2 + sd_ID__session3^2 + sd_ID__session4^2) %>% sqrt)
    (post_om$b_pdi_t) %>% sf(1)
    (post_om$b_pdi_t+ post_om$`b_pdi_t:session2`) %>% sf(1)
    (post_om$b_pdi_t+ post_om$`b_pdi_t:session3`) %>% sf(1)
    (post_om$b_pdi_t+ post_om$`b_pdi_t:session4`) %>% sf(1)
    
  }

  



## plot change points /trials with PPCs of linear models  ####
  data_trial$confidence_next_trial = data_trial$confidence %>% lead
  model_name= "conf_PE_pdi.rds"
  
  if (!file.exists(paste("Brms_models/", model_name,sep = ""))) {
    print(model_name)
    baym <- brm(confidence_next_trial ~ (PE_t_s)*pdi_t + (PE_t_s|ID),
                data = data_trial,
                prior =c(set_prior("cauchy(0,2)", class = "sd"),
                         set_prior("normal(0,1)", class = "b"),
                         set_prior("lkj(2)", class = "cor")),
                warmup = 1000, iter = 3000, chains =4)
    # control = list(adapt_delta = 0.95, max_treedepth = 12))
    saveRDS(baym,file =  paste("Brms_models/", model_name,sep = ""))
  }
  conf_model <- readRDS("Brms_models/conf_PE_pdi.rds")

  new_data_pdi <- fitted(conf_model , newdata = data_trial)
  
  data_trial$conf_next_trial_pp = new_data_pdi[,1]
  data_trial$conf_next_trial_lowci = new_data_pdi[,3]
  data_trial$conf_next_trial_highci = new_data_pdi[,4]
 
  

  
  data_trial_id <- data_trial   %>% filter(!is.na(pdi_total_cat), trials_in_block<10, trials_in_block>0, trials !=1)  %>% 
    mutate(trials_in_block =as.factor(trials_in_block)) %>% 
    group_by(ID, trials_in_block, pdi_total_cat) %>% 
    summarize(conf_id =  mean(confidence_next_trial,na.rm = TRUE),
              conf_pred_id =  mean(conf_next_trial_pp,na.rm = TRUE),
              conf_pred_lci =  mean(conf_next_trial_lowci,na.rm = TRUE),
              conf_pred_hci =  mean(conf_next_trial_highci,na.rm = TRUE)) %>%  
    group_by( trials_in_block, pdi_total_cat) %>% 
    summarize(N = n(),
              mean_conf_id =  median(conf_id,na.rm = TRUE),
              mean_conf_id_pp =  median(conf_pred_id,na.rm = TRUE),
              mean_conf_id_pp_lowci =  median(conf_pred_lci,na.rm = TRUE),
              mean_conf_id_pp_highci =  median(conf_pred_hci,na.rm = TRUE),
              se_conf_id =  sd(conf_id,na.rm = TRUE)/sqrt(N-1),
              se_conf_pred_id =  sd(conf_pred_id,na.rm = TRUE)/sqrt(N-1))
  
  
  
  g_conf_PDI <- data_trial_id %>% filter(!is.na(trials_in_block) ) %>% mutate(PDI = factor(pdi_total_cat, levels = c("HighPDI", "LowPDI"), labels = c("High", "Low"))) %>%
    ggplot(aes(x = trials_in_block, y = mean_conf_id, fill = PDI, colour = PDI)) + 
    geom_errorbar(aes(ymin =mean_conf_id - se_conf_id, ymax = mean_conf_id + se_conf_id ),  position = position_dodge(0.5), width = 0, size = 1.2)+
   
    geom_point(position = position_dodge(0.5), size = 3, shape= 21, colour = "black") +#, colour = "black") +
    geom_line(aes(x = trials_in_block, y = mean_conf_id_pp, group = PDI)) +
    geom_ribbon(aes(ymin = mean_conf_id_pp_lowci, ymax = mean_conf_id_pp_highci, fill = PDI, group = PDI), alpha = 0.1, colour = "black", size = 0.1 )+ # fill = "#E6E6E6",
    theme_Publication() +  scale_colour_manual(values =   pdi_colour )+
    scale_fill_manual(values =   pdi_colour)+
    theme(axis.title.y = element_text(size = 15),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(size = 15),
          axis.line = element_line(size = 1.5),
          legend.title = element_text("PDI"), 
          legend.background = element_rect(colour="black"),
          legend.position = c(0.8,0.2)) +
    # scale_x_discrete(labels = c("Low PDI", "High PDI")) +
    ylab("Confidence") + xlab("Trials after change point")
  g_conf_PDI
  
  ## stats for above
  post = posterior_samples(conf_model)
  (post$`b_PE_t_s:pdi_t`) %>% sf(1)
  # effect size:
  post$total_var = sqrt(post$sigma^2 + post$sd_ID__Intercept^2 + post$sd_ID__PE_t_s^2)
    post = post/post$total_var
  (post$`b_PE_t_s:pdi_t`) %>% sf(1)
  (post$b_pdi_t) %>% sf(1)
   
  
# plot figure 4 part 1  ####
    
g_fig4_part1 = plot_grid(g_pdi+ylab("") , g_om2_mean, g_conf_PDI, nrow = 1, align = "h")
g_fig4_part1

ggsave("g_fig4_part1.png", plot = g_fig4_part1, width = 13, height = 5, dpi = 300)
