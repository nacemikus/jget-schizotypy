data {
  int<lower=1> N;
  int<lower=1> T;
  // int<lower=0, upper=T> Tsubj[N,4];
  int<lower=-1, upper=100> outcome[N,4*T];  // 1: 1 - 4
  real<lower=-1, upper=100> response[N,4*T];  // 1: 1 - 4
  
  // real pdi[N];
  // int<lower=0, upper=1> ankk[N]; kakikivi
}
transformed data {
}
parameters {
  // Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[3] mu_p;
  
  vector<lower=0>[3] sigma_across_subjects; // one for each paramater 
  vector<lower=0>[3] sigma_within_subjects;
  // vector<lower=0>[3] sigma_ses;
  
  vector<lower=0>[1] mu_noise;
  vector<lower=0>[1] sigma_noise_across_subjects; // one for each paramater 
  vector<lower=0>[1] sigma_noise_within_subjects;
  // vector<lower=0>[1] sigma_noise_ses;
  
  vector[N] om_z;
  vector[N] th_z;
  vector[N] oma_z;
  vector[N] noise_z;
  
  matrix[3,N] delta_om_z; // N for each session after the first
  matrix[3,N] delta_th_z;
  matrix[3,N] delta_oma_z;
  matrix[3,N] delta_noise_z;
  
  vector[3] om_q; // one for each session after the first
  vector[3] th_q;
  vector[3] oma_q;
  vector[3] noise_q;
  
}
transformed parameters {
  // Transform subject-level raw parameters
  vector[N] om_s1;
  vector[N] th_s1;
  vector[N] oma_s1;
  vector[N] om_s2;
  vector[N] th_s2;
  vector[N] oma_s2;
  vector[N] om_s3;
  vector[N] th_s3;
  vector[N] oma_s3;
  vector[N] om_s4;
  vector[N] th_s4;
  vector[N] oma_s4;
  
  vector[N] noise_s1;
  vector[N] noise_s2;
  vector[N] noise_s3;
  vector[N] noise_s4;
  
  // matrix[5, N] r1;
  // 
  // r1 = (diag_pre_multiply(sigma,L_Omega) * z);
  // ordered[10] c[N];
  
  for (i in 1:N) {
    om_s1[i] =  mu_p[1] + sigma_across_subjects[1] *om_z[i]; // to draw the mean from N(-3,1)
    th_s1[i] = mu_p[2] + sigma_across_subjects[2] *th_z[i] - 6;
    oma_s1[i] =  mu_p[3] + sigma_across_subjects[3] *oma_z[i] - 6;
    noise_s1[i] =  mu_noise[1] + sigma_noise_across_subjects[1] *noise_z[i] +1;
    
    om_s2[i] =   om_s1[i] + delta_om_z[1,i]*sigma_within_subjects[1] + om_q[1];//*sigma_ses[1];
    th_s2[i] =  th_s1[i] + delta_th_z[1,i]*sigma_within_subjects[2] + th_q[1];//*sigma_ses[2];
    oma_s2[i] =  oma_s1[i] + delta_oma_z[1,i]*sigma_within_subjects[3] + oma_q[1];//*sigma_ses[3];
    noise_s2[i] =   noise_s1[i] + delta_noise_z[1,i]*sigma_noise_within_subjects[1] + noise_q[1];//*sigma_noise_ses[1];
    
    
    om_s3[i] =   om_s1[i]+ delta_om_z[2,i]*sigma_within_subjects[1] + om_q[2];//*sigma_ses[1];
    th_s3[i] =  th_s1[i] + delta_th_z[2,i]*sigma_within_subjects[2] + th_q[2];//*sigma_ses[2];
    oma_s3[i] =  oma_s1[i] + delta_oma_z[2,i]*sigma_within_subjects[3] + oma_q[2];//*sigma_ses[3];
    noise_s3[i] =  noise_s1[i] + delta_noise_z[2,i]*sigma_noise_within_subjects[1] + noise_q[2];//*sigma_noise_ses[1];
    
    om_s4[i] =   om_s1[i] + delta_om_z[3,i]*sigma_within_subjects[1] + om_q[3];//*sigma_ses[1];
    th_s4[i] =  th_s1[i] + delta_th_z[3,i]*sigma_within_subjects[2] + th_q[3];//*sigma_ses[2];
    oma_s4[i] =  oma_s1[i]+ delta_oma_z[3,i]*sigma_within_subjects[3] + oma_q[3];//*sigma_ses[3];
    noise_s4[i] =   noise_s1[i] + delta_noise_z[3,i]*sigma_noise_within_subjects[1] + noise_q[3];//*sigma_noise_ses[1];
    
    
  }
}
model {
  // Hyperparameters
  mu_p[1]  ~ normal(0, 1);//normal(0, 1);
  mu_p[2]  ~ normal(0, 1);//normal(0, 1);
  mu_p[3]  ~ normal(0, 1);//normal(0, 1);
  mu_noise[1]  ~ normal(0, 1);//normal(0, 1);
  
  
  sigma_across_subjects ~ normal(0, 1);//cauchy(0,5); // normal(0, 1);
  sigma_within_subjects ~ normal(0, 1);//cauchy(0,5); // normal(0, 1);
  // sigma_ses ~ normal(0, 1);//cauchy(0,5); // normal(0, 1);
  
  sigma_noise_across_subjects ~ normal(0, 1);//cauchy(0,5); // normal(0, 1);
  sigma_noise_within_subjects ~ normal(0, 1);//cauchy(0,5); // normal(0, 1);
  // sigma_noise_ses ~ normal(0, 1);//cauchy(0,5); // normal(0, 1);
  
  
  om_z~ normal(0, 1);
  th_z~ normal(0, 1);
  oma_z~ normal(0, 1);
  noise_z~ normal(0, 1);
  
  
  om_q~ normal(0, 1); // one for each session after the first
  th_q~ normal(0, 1);
  oma_q~ normal(0, 1);
  noise_q~ normal(0, 1);
  
  delta_om_z[1,] ~ normal(0, 1); // N for each session after the first
  delta_th_z[1,]~ normal(0, 1);
  delta_oma_z[1,]~ normal(0, 1);
  delta_noise_z[1,]~ normal(0, 1);
  
  delta_om_z[2,] ~ normal(0, 1); // N for each session after the first
  delta_th_z[2,]~ normal(0, 1);
  delta_oma_z[2,]~ normal(0, 1);
  delta_noise_z[2,]~ normal(0, 1);
  
  delta_om_z[3,] ~ normal(0, 1); // N for each session after the first
  delta_th_z[3,]~ normal(0, 1);
  delta_oma_z[3,]~ normal(0, 1);
  delta_noise_z[3,]~ normal(0, 1);
  
  
  // define variables
  // real mux[N,4*T];
  // real pix[N,4*T];
  // real mua[N,4*T];
  // real pia[N,4*T];
  // real mux2[N,4*T];
  // real pix2[N,4*T];
  // real mua2[N,4*T];
  // real pia2[N,4*T];
  // 
  // real muuhat[N,4*T];
  // real piuhat[N,4*T];
  // real muxhat[N,4*T];
  // real pixhat[N,4*T];
  // real muahat[N,4*T];
  // real piahat[N,4*T];
  // real daux[N,4*T];
  // real daua[N,4*T];
  // real wx[N,4*T];
  // real dax[N,4*T];
  // real wa[N,4*T];
  // real daa[N,4*T];
  
  for (i in 1:N)   {
    real mux[4*T,2];
    real pix[4*T,2];
    real mua[4*T,2];
    real pia[4*T,2];
    
    real muuhat[4*T];
    real piuhat[4*T];
    real muxhat[4*T,2];
    real pixhat[4*T,2];
    real muahat[4*T,2];
    real piahat[4*T,2];
    real daux[4*T];
    real daua[4*T];
    real wx[4*T];
    real dax[4*T];
    real wa[4*T];
    real daa[4*T];
    
    real session;
    real om_par;
    real oma_par;
    real th_par;
    real noise_par;
    
    real kax_par;
    real kaa_par;
    real tha_par;
    
    // placeholder variables 
    real vvx;
    real vva;
    real pixmhat;
    real piamhat;
    real wwx;
    real wwa;
    real rrx;
    real rra;
    real ddx;
    real dda;
    
    // real kau_good_par;
    
    for (t in 1:4*T)  {
      session = ceil(t/240.0);
      
      if (t%240 == 1 ) {
        // new session starts
        // print("trial ", t, " inital settings, for session", session);
        //initial variables
        mux[t,1] = outcome[i,t];
        mux[t,2] = 3.0;
        
        pix[t,1] = 1.0/3.0;
        pix[t,2] = 1.0/0.1;
        
        mua[t,1] = log(5.0^2);
        mua[t,2] = 1.0;
        
        pia[t,1] = 1./3;
        pia[t,2] = 1./0.1;
        
        if (session == 1) {
          om_par = om_s1[i];
          th_par = exp(th_s1[i]);
          oma_par = oma_s1[i];
          noise_par  = exp(noise_s1[i]);
          
          
        } else if (session == 2) {
          om_par = om_s2[i];
          th_par = exp(th_s2[i]);
          oma_par = oma_s2[i];
          noise_par  = exp(noise_s2[i]);
          
          
        } else if (session == 3) {
          om_par = om_s3[i];
          th_par = exp(th_s3[i]);
          oma_par = oma_s3[i];
          noise_par  = exp(noise_s3[i]);
        } else if (session == 4) {
          om_par = om_s4[i];
          th_par = exp(th_s4[i]);
          oma_par = oma_s4[i];
          noise_par  = exp(noise_s4[i]);
        }
        
        
        kax_par = 1;
        kaa_par = 1;
        tha_par = exp(-10);
        // kau_good_par = 1;
        
      } else {
        
        if (outcome[i,t] != -1 && response[i,t] != -1) {
          
          muuhat[t] = mux[t-1,1];
          
          // Precision of prediction
          piuhat[t] = 1/exp(mua[t-1,1]);
          
          // Mean prediction error
          daux[t] = outcome[i,t] - muuhat[t];
          
          // 1st level
          // ~~~~~~~~~
          // Predictions
          muxhat[t,1] = mux[t-1,1];
          muahat[t,1] = mua[t-1,1];
          
          // Precisions of predictions
          pixhat[t,1] = 1/(1/pix[t-1,1] +exp(kax_par *mux[t-1,2] +om_par));
          piahat[t,1] = 1/(1/pia[t-1,1] +exp(kaa_par *mua[t-1,2] +oma_par));
          
          if (pixhat[t,1]<0)
          reject("pixhat[t,1] <0 not allowed, trial ", t, " for subject ", i, " session ", session, " <0 here: ", pixhat[t,1],  ", ", pix[t-1,1]);
          
          // likelihood function
          
          
          response[i,t] ~ normal(muxhat[t,1], noise_par);
          
          // x-updates
          pix[t,1] = pixhat[t,1] +piuhat[t];
          mux[t,1] = muxhat[t,1] +piuhat[t]/pix[t,1] *daux[t];
          
          // Prediction error of input precision
          daua[t] = (1/pix[t,1] +(mux[t,1] -outcome[i,t])^2) *piuhat[t] -1;
          
          // alpha-updates
          pia[t,1] = piahat[t,1] +1.0/2 *(1 +daua[t]);
          mua[t,1] = muahat[t,1] +1.0/2 *1/pia[t,1] *daua[t];
          
          // Volatility prediction errors
          dax[t] = (1/pix[t,1] +(mux[t,1] -muxhat[t,1])^2) *pixhat[t,1] -1;
          daa[t] = (1/pia[t,1] +(mua[t,1] -muahat[t,1])^2) *piahat[t,1] -1;
          // Last level
          // ~~~~~~~~~~
          // Predictions
          muxhat[t,2] = mux[t-1,2];
          muahat[t,2] = mua[t-1,2];
          
          // Precision of prediction
          pixhat[t,2] = 1/(1/pix[t-1,2] +th_par);
          piahat[t,2] = 1/(1/pia[t-1,2] +tha_par);
          
          if (pixhat[t,2]<0)
          reject("pixhat[t,2] <0 not allowed, trial ", t, " for subject ", i, " session ", session, " <0 here: ", pixhat[t,2],  ", ", pix[t-1,2]);
          
          
          // Weighting factor
          wx[t] = exp(kax_par *mux[t-1,2] +om_par) *pixhat[t,1];
          wa[t] =exp(kaa_par *mua[t-1,2] +oma_par) *piahat[t,1];
          
          // Mean updates
          mux[t,2] = muxhat[t,2] +1.0/2 *1/pixhat[t,2] *kax_par *wx[t] *dax[t];
          mua[t,2] = muahat[t,2] +1.0/2 *1/piahat[t,2] *kaa_par *wa[t] *daa[t];
          
          
          // Ingredients of precision updates which depend on the mean
          // updates
          
          
          vvx = exp(kax_par *mux[t,2] +om_par);
          vva = exp(kaa_par *mua[t,2] +oma_par);
          pixmhat = 1/(1/pix[t-1,1] +vvx); 
          piamhat = 1/(1/pia[t-1,1] +vva); 
          wwx = vvx *pixmhat;
          wwa = vva *piamhat;
          rrx = (vvx -1/pix[t-1,1]) *pixmhat;
          rra = (vva -1/pia[t-1,1]) *piamhat;
          ddx = (1/pix[t,1] +(mux[t,1] -muxhat[t,1])^2) *pixmhat -1;
          dda = (1/pia[t,1] +(mua[t,1] -muahat[t,1])^2) *piamhat -1;
          
          // Precision update
          pix[t,2] = pixhat[t,2] +fmax(0.0, 1.0/2 *kax_par^2 *wwx*(wwx +rrx*ddx));
          pia[t,2] = piahat[t,2] +fmax(0.0, 1.0/2 *kaa_par^2 *wwa*(wwa +rra*dda));
          
          // Volatility prediction error
          // dax[t,2] = (1/pix[t,2] +(mux[t,2] -muxhat[t,2])^2) *pixhat[t,2] -1;
          // daa[t,2] = (1/pia[t,2] +(mua[t,2] -muahat[t,2])^2) *piahat[t,2] -1;
          // 
          if (is_nan(mux[t,2]) || is_nan(mux[t,1]) || is_nan(pix[t,1]) || is_nan(pix[t,2]) )
          reject("nans not allowed, trial ", t, " for subject ", i, " session ", session, " nans here: ", mux[t,1],
          ", ", mux[t,2], ", ", pix[t,1], ", ", pix[t,2], ", ",mua[t,1], ", ",mua[t,2], ", ", pixhat[t,2], ", ", om_par, ", ", oma_par, ", ", th_par);
          
          
          
        } else {
          mux[t,1] = mux[t-1,1];
          mux[t,2] =  mux[t-1,2];
          
          pix[t,1] = pix[t-1,1];
          pix[t,2] = pix[t-1,2];
          
          mua[t,1] =mua[t-1,1];
          mua[t,2] =  mua[t-1,2];
          
          pia[t,1] = pia[t-1,1];
          pia[t,2] =pia[t-1,2];
          
        }// end of !isnan(outcome(t)) block
      }// end of first session trial block
      
    }// end of t (trial) loop

  } // end of i (subject) loop
}

generated quantities {

//   // For group level parameters
//   // real mu_om_good;
//   // real mu_om_bad;
//   // real mu_beta;
//   // // real<lower=0,upper=1> mu_al;
//   // // real<lower=0>         mu_beta2;
//   // 
//   // real mu_mu0;
//   
//   
//   // For log likelihood calculation
//   real log_lik[N];
//   // int count_trial;
//   // For posterior predictive check
//   real y_pred[N,4*T];
//   
//   // For model diagnosis
//   real mux_mat[N,4*T];
//   real pix_mat[N,4*T];
//   real mua_mat[N,4*T];
//   real pia_mat[N,4*T];
//   real mux2_mat[N,4*T];
//   real pix2_mat[N,4*T];
//   real mua2_mat[N,4*T];
//   real pia2_mat[N,4*T];
//   
//   // real y_pred_step2[N,T];
//   
//   
//   // Set all posterior predictions to 0 (avoids NULL values)
//   for (i in 1:N) {
//     for (t in 1: 4*T) {
//       y_pred[i,t] = -1;
//       mux_mat[i,t] = -1;
//       pix_mat[i,t] = -1;
//       mua_mat[i,t] = -1;
//       pia_mat[i,t] = -1;
//       mux2_mat[i,t] = -1;
//       pix2_mat[i,t] = -1;
//       mua2_mat[i,t] = -1;
//       pia2_mat[i,t] = -1;
//       // y_pred_step2[i,t] = -1;
//     }
//   }
//   for (i in 1:N) {
//     log_lik[i]=0;
//   }
//   // Generate group level parameter values
//   // mu_om_good     = mu_p[1]-4;
//   // mu_om_bad     = mu_p[2]-4;
//   // mu_beta  =  mu_p[3] ;
//   // // mu_beta2  = exp( mu_p[4] );
//   // mu_mu0     =  mu_p[4];
//   
//   { // local section, this saves time and space
//   // count_trial= 1;
//   for (i in 1:N) {
//     
//     real mux[4*T,2];
//     real pix[4*T,2];
//     real mua[4*T,2];
//     real pia[4*T,2];
//     
//     real muuhat[4*T];
//     real piuhat[4*T];
//     real muxhat[4*T,2];
//     real pixhat[4*T,2];
//     real muahat[4*T,2];
//     real piahat[4*T,2];
//     real daux[4*T];
//     real daua[4*T];
//     real wx[4*T];
//     real dax[4*T];
//     real wa[4*T];
//     real daa[4*T];
//     
//     real session;
//     real om_par;
//     real oma_par;
//     real th_par;
//     real noise_par;
//     
//     real kax_par;
//     real kaa_par;
//     real tha_par;
//     
//     // placeholder variables 
//     real vvx;
//     real vva;
//     real pixmhat;
//     real piamhat;
//     real wwx;
//     real wwa;
//     real rrx;
//     real rra;
//     real ddx;
//     real dda;
//     
//     // real kau_good_par;
//     
//     for (t in 1:4*T)  {
//       if (outcome[i,t] != -1 && response[i,t] != -1) {
//         
//         session = ceil(t/240.0);
//         
//         if (t%240 == 1 ) {
//           // new session starts
//           
//           //initial variables
//           mux[t,1] =outcome[i,t];
//           mux[t,2] = 3.0;
//           
//           pix[t,1] = 1.0/3.0;
//           pix[t,2] = 1.0/0.1;
//           
//           mua[t,1] = log(5.0^2);
//           mua[t,2] = 1.0;
//           
//           pia[t,1] = 1./3;
//           pia[t,2] = 1./0.1;
//           
//           if (session == 1) {
//             om_par = om_s1[i];
//             th_par = exp(th_s1[i]);
//             oma_par = oma_s1[i];
//             noise_par = exp(noise_s1[i]);
//             
//           } else if (session == 2) {
//             om_par = om_s2[i];
//             th_par = exp(th_s2[i]);
//             oma_par = oma_s2[i];
//             noise_par = exp(noise_s2[i]);
//             
//             
//           } else if (session == 3) {
//             om_par = om_s3[i];
//             th_par = exp(th_s3[i]);
//             oma_par = oma_s3[i];
//             noise_par = exp(noise_s3[i]);
//           } else if (session == 4) {
//             om_par = om_s4[i];
//             th_par = exp(th_s4[i]);
//             oma_par = oma_s4[i];
//             noise_par = exp(noise_s4[i]);
//           }
//           
//           
//           kax_par = 1;
//           kaa_par = 1;
//           tha_par = exp(-10);
//           // kau_good_par = 1;
//           
//         } else {
//           
//           muuhat[t] = mux[t-1,1];
//           
//           // Precision of prediction
//           piuhat[t] = 1/exp(mua[t-1,1]);
//           
//           // Mean prediction error
//           daux[t] = outcome[i,t] - muuhat[t];
//           
//           // 1st level
//           // ~~~~~~~~~
//           // Predictions
//           muxhat[t,1] = mux[t-1,1];
//           muahat[t,1] = mua[t-1,1];
//           
//           // Precisions of predictions
//           pixhat[t,1] = 1/(1/pix[t-1,1] +exp(kax_par *mux[t-1,2] +om_par));
//           piahat[t,1] = 1/(1/pia[t-1,1] +exp(kaa_par *mua[t-1,2] +oma_par));
//           
//           
//           // likelihood function
//           log_lik[i] = log_lik[i] + normal_lpdf( response[i,t] | muxhat[t,1], noise_par );
//           y_pred[i,t] = normal_rng(muxhat[t,1], noise_par);
//           
//           // x-updates
//           pix[t,1] = pixhat[t,1] +piuhat[t];
//           mux[t,1] = muxhat[t,1] +piuhat[t]/pix[t,1] *daux[t];
//           
//           // Prediction error of input precision
//           daua[t] = (1/pix[t,1] +(mux[t,1] -outcome[i,t])^2) *piuhat[t] -1;
//           
//           // alpha-updates
//           pia[t,1] = piahat[t,1] +1.0/2 *(1 +daua[t]);
//           mua[t,1] = muahat[t,1] +1.0/2 *1/pia[t,1] *daua[t];
//           
//           // Volatility prediction errors
//           dax[t] = (1/pix[t,1] +(mux[t,1] -muxhat[t,1])^2) *pixhat[t,1] -1;
//           daa[t] = (1/pia[t,1] +(mua[t,1] -muahat[t,1])^2) *piahat[t,1] -1;
//           // Last level
//           // ~~~~~~~~~~
//           // Predictions
//           muxhat[t,2] = mux[t-1,2];
//           muahat[t,2] = mua[t-1,2];
//           
//           // Precision of prediction
//           pixhat[t,2] = 1/(1/pix[t-1,2] +th_par);
//           piahat[t,2] = 1/(1/pia[t-1,2] +tha_par);
//           
//           // Weighting factor
//           wx[t] = exp(kax_par *mux[t-1,2] +om_par) *pixhat[t,1];
//           wa[t] =exp(kaa_par *mua[t-1,2] +oma_par) *piahat[t,1];
//           
//           // Mean updates
//           mux[t,2] = muxhat[t,2] +1.0/2 *1/pixhat[t,2] *kax_par *wx[t] *dax[t];
//           mua[t,2] = muahat[t,2] +1.0/2 *1/piahat[t,2] *kaa_par *wa[t] *daa[t];
//           
//           // Ingredients of precision updates which depend on the mean
//           // updates
//           
//           
//           vvx = exp(kax_par *mux[t,2] +om_par);
//           vva = exp(kaa_par *mua[t,2] +oma_par);
//           pixmhat = 1/(1/pix[t-1,1] +vvx); 
//           piamhat = 1/(1/pia[t-1,1] +vva); 
//           wwx = vvx *pixmhat;
//           wwa = vva *piamhat;
//           rrx = (vvx -1/pix[t-1,1]) *pixmhat;
//           rra = (vva -1/pia[t-1,1]) *piamhat;
//           ddx = (1/pix[t,1] +(mux[t,1] -muxhat[t,1])^2) *pixmhat -1;
//           dda = (1/pia[t,1] +(mua[t,1] -muahat[t,1])^2) *piamhat -1;
//           
//           // Precision update
//           pix[t,2] = pixhat[t,2] +fmax(0.0, 1.0/2 *kax_par^2 *wwx*(wwx +rrx*ddx));
//           pia[t,2] = piahat[t,2] +fmax(0, 1.0/2 *kaa_par^2 *wwa*(wwa +rra*dda));
//           
//           // Volatility prediction error
//           // dax[t,2] = (1/pix[t,2] +(mux[t,2] -muxhat[t,2])^2) *pixhat[t,2] -1;
//           // daa[t,2] = (1/pia[t,2] +(mua[t,2] -muahat[t,2])^2) *piahat[t,2] -1;
//           // 
//           
//           // save for model diagnotics
//           
//           mux_mat[i,t] = mux[t,1];
//           pix_mat[i,t] = pix[t,1];
//           mua_mat[i,t] = mua[t,1];
//           pia_mat[i,t] = pia[t,1];
//           mux2_mat[i,t] = mux[t,2];
//           pix2_mat[i,t] =  pix[t,2];
//           mua2_mat[i,t] = mua[t,2];
//           pia2_mat[i,t] = pia[t,2];
//         }
//       } // end of !isnan(outcome(t)) block
//     }// end of t (trial) loop
//     
//     
//   } // end of i loop
//   
//   
//   
//   } // end local
}

