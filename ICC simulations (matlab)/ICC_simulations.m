%% Looking for outcome trajectories for the jumping gaussian estimation (predictive inference) task
% under which conditions are good test-retest scores of model parameters
% feasible?

% The aim is to :
% - simulate data from priors for two different reward trajectories (created by the function "create_rew_traj.m" and
% - look at correlations 

% - repeat the process for different trial numbers and levels of decison noise

% - repeat the proces for different models and outcome typs (e.g. binary vs
%   continuous)
%% Put things on path
get_path = mfilename('fullpath');
[GenFolder, name, ext] = fileparts(which(get_path));
% if any(ismember('P:\Projects\',path_list_cell))
%     rmpath(genpath('P:\Projects\'));
% end

ScriptsFolder = [GenFolder,'\hgf-toolbox']; % your path to tapas hgf
folder_save = 'Simulations_forTRT_analysis';

cd (GenFolder)
cd ..
mainFolder = pwd;

addpath(genpath(mainFolder));%'ReinfLearn')); % modelling code
addpath(genpath(ScriptsFolder));%


analysis_type = input('choose the analysis type: binary/cont_fixed_sd/cont/: ', 's');
% analysis_type = 'binary'; % {'binary', 'cont_fixed_sd', 'cont'}
name_analysis = input('choose the model type: rw/ehgf/ehgf_jget/nassar_delta: ', 's');
name_prefix = input('add a prefix: (eg. InferNoise)', 's');

n_trials_set = [16, 48, 120, 240,480];
n_trials_set = input('trial numbers? if empty :[16, 48, 120, 240,480]');



%% settings

%n_trials_set = [16,48,120,240, 480];

no_sim = 50;
trials_per_block = [3,30];%[49,51;49,51];% min max for each of the two volatility stages
%     sd_options = [5,15];
prob_set = [0.2,0.8,0.5];
p_vol = 0.1;
rew_interval = [10,90]; % [-6,6];
%{'rw', 'ehgf', 'nassar', 'ehgf_jget'}
% analysis_table = table({'rw', 'ehgf'})
%% binary
if strcmp(analysis_type, 'binary')
    noize_levels = [1,8,48];
    
    
    if strcmp(name_analysis, 'ehgf')
        trans_fu = @tapas_ehgf_binary_transp;
        c = tapas_ehgf_binary_config; %tapas_ehgf_binary_config;
    elseif strcmp(name_analysis, 'rw')
        trans_fu = @tapas_rw_binary_transp;
        c = tapas_rw_binary_config; %tapas_ehgf_binary_config;
    end
    
    for nt = 1: length(n_trials_set)
        %% set the analysis configs:
        
        
        %% run the simulation
        ntrials= n_trials_set(nt);% no of trials,
        
        
        mu_set = c.priormus;
        sa_set = c.priorsas;
        %     co = tapas_softmax_binary_config;
        pars_sim = NaN(no_sim,length(sa_set)+1,length(noize_levels));
        % noize_sim = NaN(no_sim,1);
        % pars_estuse= NaN(no_sim,16);
        % noize_estuse =  NaN(no_sim,1);
        pars_est1= NaN(no_sim,length(sa_set)+1,length(noize_levels));
        pars_est2= NaN(no_sim,length(sa_set)+1,length(noize_levels));
        u_traj = NaN(no_sim,ntrials);
        u2_traj = NaN(no_sim,ntrials);
        %
        % [u, p_traj] = create_rew_traj_binary(ntrials, prob_set, trials_per_block, p_vol);
        %     pars_sim(sim,:,n)=  tapas_eeehgf_jget_transp(levels_ph, pars_sim(sim,:,n));
        %
        
        %
        % mu_set(1) = 62;
        % plotGrahps = 0;
        %% simulate, estimate
        
        for n = 1 : length(noize_levels)
            pars_sim(:,end,n) = noize_levels(n);
            for sim = 1: no_sim
                pars_sim(sim,1:length(sa_set),n) = 0;
                
                sim_succesful = 0;
                while sim_succesful == 0
                    try
                        for pp = 1 :length(sa_set)
                            
                            
                            pars_sim(sim,pp,n) = randn(1,1)*sqrt(sa_set(pp)) +  mu_set(pp);
                            
                        end
                        
                        r_levels.c_prc.n_levels = 3;
                        pars_sim(sim,1:end-1,n)=  trans_fu(r_levels, pars_sim(sim,1:end-1,n));
                        %     noize_sim(sim) =  randi([1,100],1);%tapas_gaussian_obs_transp([], co.priormus + sqrt(co.priorsas)*randn(1));
                        %
                        u_traj(sim,:) = create_rew_traj_binary(ntrials, prob_set, trials_per_block, p_vol);
                        u2_traj(sim,:) = create_rew_traj_binary(ntrials, prob_set, trials_per_block, p_vol);
                        pars_sim_1 = pars_sim(sim,1:end-1,n);
                        pars_sim_2 = pars_sim(sim,1:end-1,n);
                        pars_sim_1(1) = u_traj(sim,1);
                        pars_sim_2(1) = u2_traj(sim,1);
                        
                        %                     u_traj(sim,1:ntrials/4) = create_rew_traj_binary(120, prob_set, trials_per_block, 0.01);
                        %                     u_traj(sim,ntrials/4+1 :ntrials/4*2) = create_rew_traj_binary(120, prob_set, trials_per_block, 0.1);
                        %                     u_traj(sim,ntrials/4*2+1:ntrials/4*3) = create_rew_traj_binary(120, prob_set, trials_per_block, 0.01);
                        %                     u_traj(sim,ntrials/4*3+1:ntrials/4*4) = create_rew_traj_binary(120, prob_set, trials_per_block, 0.1);
                        %                      est1= tapas_fitModel(sim1.y ,...
                        %                         u_traj(sim,:)',...
                        %                         ['tapas_',name_analysis,'_binary_config'],...% 'tapas_rw_binary',...
                        %                         'tapas_softmax_binary_config',...
                        %                         'tapas_quasinewton_optim_config');
                        %
                        %                     pars =  pars_sim(sim,1:end-1,n);
                        %                     pars(13:14) = [1,-4];
                        sim1 = tapas_simModel( u_traj(sim,:)',...
                            ['tapas_',name_analysis,'_binary'],...% 'tapas_rw_binary',...
                            pars_sim_1,...
                            'tapas_softmax_binary',...
                            noize_levels(n));
                        %                     tapas_ehgf_binary_plotTraj(sim1);
                        
                        
                        sim2 = tapas_simModel( u2_traj(sim,:)',...
                            ['tapas_',name_analysis,'_binary'],...% 'tapas_rw_binary',...
                            pars_sim_2,...
                            'tapas_softmax_binary',...
                            noize_levels(n));
                        if (any(isnan(sim1.traj.muhat(:,1))) || any(isnan(sim2.traj.muhat(:,1))))
                            sim_succesful = 0;
                        else
                            sim_succesful = 1;
                        end
                        
                        
                    catch
                        disp('sim caught')
                        
                        continue
                        
                    end
                end
                %            tapas_rw_binary_plotTraj(sim1)
                %         [est_sim, estuse_sim, parset]= estimate_models(sim1.y ,u);
                est_succesful = 0;
                while est_succesful== 0
                    try
                        est1= tapas_fitModel(sim1.y ,...
                            u_traj(sim,:)',...
                            ['tapas_',name_analysis,'_binary_config'],...% 'tapas_rw_binary',...
                            'tapas_softmax_binary_config',...
                            'tapas_quasinewton_optim_config');
                        
                        est2= tapas_fitModel(sim2.y ,...
                            u2_traj(sim,:)',...
                            ['tapas_',name_analysis,'_binary_config'],...% 'tapas_ehgf_binary_config',...
                            'tapas_softmax_binary_config',...
                            'tapas_quasinewton_optim_config');
                        est_succesful = 1;
                    catch
                        disp('est caught')
                        
                        continue
                    end
                end
                pars_est1(sim,:,n) =  [est1.p_prc.p, est1.p_obs.p];
                pars_est2(sim,:,n) = [est2.p_prc.p, est2.p_obs.p];
                
                %          tapas_ehgf_binary_plotTraj(est1)
            end
        end
        save([GenFolder,'/',folder_save,'/',name_analysis,'_binary_',num2str(ntrials),'trials.mat'], 'pars_est1','pars_est2','pars_sim', 'noize_levels', 'c', 'u_traj', 'u2_traj')
    end
    return
end


%% continuous fixed variance settings
% n_trials_set = [16,48,120, 240,480];%[16,48,120, 240,480];
if strcmp(analysis_type, 'cont_fixed_sd')
    
    noize_levels = [96,48,1];
    trials_per_block = [3,30];%[49,51;49,51];% min max for each of the two volatility stages
    sd_options = 5;%[5,15];
    %     prob_set = [0.2,0.8,0.5];
    p_vol = 0.1;
    sd_vol = 0;
    rew_interval = [10,90]; % [-6,6];
    
    if strcmp(name_analysis, 'ehgf_jget')
        trans_fu = @tapas_ehgf_jget_transp;
        c = tapas_ehgf_jget_config; %tapas_ehgf_binary_config;
        c.omasa(1) = 0; %fix omega noise to -3 and do not estimate it; 

        c = tapas_align_priors(c);
        
        
    elseif strcmp(name_analysis, 'nassar_delta')
        trans_fu = @tapas_nassar_delta_transp;
        c = tapas_nassar_delta_config; %tapas_ehgf_binary_config;
    end
    
    
    for nt = 1: length(n_trials_set)
        
        ntrials= n_trials_set(nt);% no of trials,
        
        
        
        
        mu_set = c.priormus;
        sa_set = c.priorsas;
     
        pars_sim = NaN(no_sim,length(sa_set)+1,length(noize_levels));
        pars_est1= NaN(no_sim,length(sa_set)+1,length(noize_levels));
        pars_est2= NaN(no_sim,length(sa_set)+1,length(noize_levels));
        u_traj = NaN(no_sim,ntrials);
        u2_traj = NaN(no_sim,ntrials);
        %
        % [u, p_traj] = create_rew_traj_binary(ntrials, prob_set, trials_per_block, p_vol);
        %     pars_sim(sim,:,n)=  tapas_eeehgf_jget_transp(levels_ph, pars_sim(sim,:,n));
        %
        
        %
        % mu_set(1) = 62;
        % plotGrahps = 0;
        %% simulate, estimate
        
        for n = 1 : length(noize_levels)
            pars_sim(:,end,n) = noize_levels(n);
            for sim = 1: no_sim
                pars_sim(sim,1:length(sa_set),n) = 0;
                
                sim_succesful = 0;
                while sim_succesful == 0
                    try
                        u_traj(sim,:) = create_rew_traj(ntrials, rew_interval, trials_per_block, sd_options, p_vol, sd_vol);
                        u2_traj(sim,:) = create_rew_traj(ntrials, rew_interval, trials_per_block, sd_options, p_vol, sd_vol);
                      
                        for pp = 1 :length(sa_set)
                            
                            
                            pars_sim(sim,pp,n) = randn(1,1)*sqrt(sa_set(pp)) +  mu_set(pp);
                            
                        end
%                       
                        r_levels.c_prc.n_levels = 2;
                        pars_sim(sim,1:end-1,n)=  trans_fu(r_levels, pars_sim(sim,1:end-1,n));
                       
                        
                        sim1 = tapas_simModel( u_traj(sim,:)',...
                            ['tapas_',name_analysis],...% 'tapas_rw_binary',...
                            pars_sim(sim,1:end-1,n),...
                            'tapas_gaussian_obs',...
                            noize_levels(n));
                        sim2 = tapas_simModel( u2_traj(sim,:)',...
                            ['tapas_',name_analysis],...% 'tapas_rw_binary',...
                            pars_sim(sim,1:end-1,n),...
                            'tapas_gaussian_obs',...
                            noize_levels(n));
                     
                        sim_succesful = 1;
                    
                        
                    catch
                        disp('sim caught')
                        
                        continue
                        
                    end
                end
                %            tapas_ehgf_jget_plotTraj(sim1)
                %         [est_sim, estuse_sim, parset]= estimate_models(sim1.y ,u);
                est_succesful = 0;
                while est_succesful== 0
                    try
                        est1= tapas_fitModel(sim1.y ,...
                            u_traj(sim,:)',...
                            ['tapas_',name_analysis,'_config'],...% 'tapas_rw_binary',...
                            'tapas_gaussian_obs_config',...
                            'tapas_quasinewton_optim_config');
                        
                        est2= tapas_fitModel(sim2.y ,...
                            u2_traj(sim,:)',...
                            ['tapas_',name_analysis,'_config'],...% 'tapas_rw_binary',...
                            'tapas_gaussian_obs_config',...
                            'tapas_quasinewton_optim_config');
                        est_succesful = 1;
                   
                    catch
                        disp('est caught')
                        
                        continue
                    end
                end
                pars_est1(sim,:,n) =  [est1.p_prc.p, est1.p_obs.p];
                pars_est2(sim,:,n) = [est2.p_prc.p, est2.p_obs.p];
                
               
            end
        end
        save([GenFolder,'/',folder_save,'/',name_analysis,'_cont_fixed_sd', num2str(sd_options),'_',num2str(ntrials),'trials',name_prefix,'.mat'], 'pars_est1','pars_est2','pars_sim', 'noize_levels', 'c', 'u_traj', 'u2_traj')
    end
    return
end
%% continuous + infer noise
if strcmp(analysis_type, 'cont')
    noize_levels = [96,48,1];
    trials_per_block = [3,30];%[49,51;49,51];% min max for each of the two volatility stages
  
    sd_options = [5,15];
    %     prob_set = [0.2,0.8,0.5];
    p_vol = 0.1;
    sd_vol = 0.4;
    rew_interval = [10,90]; % [-6,6];
    
    if strcmp(name_analysis, 'ehgf_jget')
        trans_fu = @tapas_ehgf_jget_transp;
        c = tapas_ehgf_jget_config; %tapas_ehgf_binary_config;
    elseif strcmp(name_analysis, 'nassar_delta')
        trans_fu = @tapas_nassar_delta_transp;
        c = tapas_nassar_delta_config; %tapas_ehgf_binary_config;
    end
    
   
    for nt = 1: length(n_trials_set)
        
        ntrials= n_trials_set(nt);% no of trials,
        
        
        mu_set = c.priormus;
        sa_set = c.priorsas;
 
        pars_sim = NaN(no_sim,length(sa_set)+1,length(noize_levels));
    
        pars_est1= NaN(no_sim,length(sa_set)+1,length(noize_levels));
        pars_est2= NaN(no_sim,length(sa_set)+1,length(noize_levels));
        u_traj = NaN(no_sim,ntrials);
        u2_traj = NaN(no_sim,ntrials);
     
        %% simulate, estimate
        
        for n = 1 : length(noize_levels)
            pars_sim(:,end,n) = noize_levels(n);
            for sim = 1: no_sim
                pars_sim(sim,1:length(sa_set),n) = 0;
                
                sim_succesful = 0;
                while sim_succesful == 0
                    try
                        
                        for pp = 2 :length(sa_set)
                            
                            
                            pars_sim(sim,pp,n) = randn(1,1)*sqrt(sa_set(pp)) +  mu_set(pp);
                            
                        end
                        
                        r_levels.c_prc.n_levels = 2;
                        pars_sim(sim,1:end-1,n)=  trans_fu(r_levels, pars_sim(sim,1:end-1,n));
                        %     noize_sim(sim) =  randi([1,100],1);%tapas_gaussian_obs_transp([], co.priormus + sqrt(co.priorsas)*randn(1));
                        %
                        u_traj(sim,:) = create_rew_traj(ntrials, rew_interval, trials_per_block, sd_options, p_vol, sd_vol);
                        u2_traj(sim,:) = create_rew_traj(ntrials, rew_interval, trials_per_block, sd_options, p_vol, sd_vol);
                        pars_sim_1 = pars_sim(sim,1:end-1,n);
                        pars_sim_2 = pars_sim(sim,1:end-1,n);
                        pars_sim_1(1) = u_traj(sim,1);
                        pars_sim_2(1) = u2_traj(sim,1);
                        %                     u_traj(sim,1:ntrials/4) = create_rew_traj_binary(120, prob_set, trials_per_block, 0.01);
                        %                     u_traj(sim,ntrials/4+1 :ntrials/4*2) = create_rew_traj_binary(120, prob_set, trials_per_block, 0.1);
                        %                     u_traj(sim,ntrials/4*2+1:ntrials/4*3) = create_rew_traj_binary(120, prob_set, trials_per_block, 0.01);
                        %                     u_traj(sim,ntrials/4*3+1:ntrials/4*4) = create_rew_traj_binary(120, prob_set, trials_per_block, 0.1);
                        %                      est1= tapas_fitModel(sim1.y ,...
                        %                         u_traj(sim,:)',...
                        %                         ['tapas_',name_analysis,'_binary_config'],...% 'tapas_rw_binary',...
                        %                         'tapas_softmax_binary_config',...
                        %                         'tapas_quasinewton_optim_config');
                        
                        
                        sim1 = tapas_simModel( u_traj(sim,:)',...
                            ['tapas_',name_analysis],...% 'tapas_rw_binary',...
                            pars_sim_1,...
                            'tapas_gaussian_obs',...
                            noize_levels(n));
                        
                        sim2 = tapas_simModel( u2_traj(sim,:)',...
                            ['tapas_',name_analysis],...% 'tapas_rw_binary',...
                            pars_sim_2,...
                            'tapas_gaussian_obs',...
                            noize_levels(n));
                        sim_succesful = 1;
                    catch
                        disp('sim caught')
                        
                        continue
                        
                    end
                end
                %            tapas_ehgf_jget_plotTraj(sim1)
                %         [est_sim, estuse_sim, parset]= estimate_models(sim1.y ,u);
                est_succesful = 0;
                while est_succesful== 0
                    try
                        est1= tapas_fitModel(sim1.y ,...
                            u_traj(sim,:)',...
                            ['tapas_',name_analysis,'_config'],...% 'tapas_rw_binary',...
                            'tapas_gaussian_obs_config',...
                            'tapas_quasinewton_optim_config');
                         %                     tapas_ehgf_jget_plotTraj(est1)
                        est2= tapas_fitModel(sim2.y ,...
                            u2_traj(sim,:)',...
                            ['tapas_',name_analysis,'_config'],...% 'tapas_ehgf_jget_config   _binary',...
                            'tapas_gaussian_obs_config',...
                            'tapas_quasinewton_optim_config');
                        est_succesful = 1;
                        %                     tapas_ehgf_jget_plotTraj(est2)
                    catch
                        disp('est caught')
                        
                        continue
                    end
                end
                pars_est1(sim,:,n) =  [est1.p_prc.p, est1.p_obs.p];
                pars_est2(sim,:,n) = [est2.p_prc.p, est2.p_obs.p];
%                 if (est1.p_obs.p > 20 | est2.p_obs.p >20)
%                     est1.p_obs.p
%                 end
%                 if (sum((est1.p_prc.p([13,14,15]) - est2.p_prc.p([13,14,15])).^2) < 5)
%                      tapas_ehgf_binary_plotTraj(est1)
%                 end%          tapas_ehgf_binary_plotTraj(est1)
            end
        end
        save([GenFolder,'/',folder_save,'/',name_analysis,'_cont_noise_',num2str(ntrials),'trials',name_prefix,'.mat'], 'pars_est1','pars_est2','pars_sim', 'noize_levels', 'c', 'u_traj', 'u2_traj')
    end
    return
end

%% plot the reliability results
% names = {'vhat_0', 'N_0', 'hazExp'};
% c = analysis_session_config_info
%%
remove_outliers = 1;
sd_outliers = 5;
% Model_set = {'ehgf_jget_cont_noise_480trials'};%{'cont'};
Model_set = {'ehgf_jget_cont_noise'};%{'cont'};
% cd (GenFolder)
% list_dir = dir([GenFolder, '\Make Trajectories\Simulations_forTRT_analysis\']);

for s = 1: length(Model_set)
    list_dir = dir([GenFolder, '\',folder_save,'\*',Model_set{s},'*']);
    
    for l = 1 : length(list_dir)
        load([GenFolder, '\',folder_save,'\', list_dir(l).name])
        
        par_index_set = find(c.priorsas ~= 0 & ~isnan(c.priorsas));
        % noize_levels = [1,50];
        figure;
        for n = 1 : length(noize_levels)
            
            count = 0;
            for pp = 1 :length(par_index_set)
                par_ind = par_index_set(pp);
                
                pars_est1_set = pars_est1(~isnan(pars_est1(:,par_ind)), par_ind,n);
                pars_est2_set = pars_est2(~isnan(pars_est2(:,par_ind)), par_ind,n);
                
                count = count +1 ;
                subplot(length(noize_levels),length(par_index_set),(n-1)*length(par_index_set)+count)
                % remove outliers
                outliers = zeros(length(pars_est1_set), 1);
                if remove_outliers == 1
                    outliers = or(abs(pars_est1_set- mean(pars_est1_set)) > sd_outliers*std(pars_est1_set ) ,...
                        abs(pars_est2_set- mean(pars_est2_set)) > sd_outliers*std(pars_est2_set ) );
                end
                %             
                %       pars_est_set  = pars_est_set(~isnan(pars_est_set));
%                 [r,p] = corr(pars_est1_set(~outliers),  pars_est2_set(~outliers));
%                 [r,p] = corr(pars_est1_set(~outliers),  pars_est2_set(~outliers));
%                 r = ICC( [pars_est1_set(~outliers),pars_est2_set(~outliers)])
                %                 figure;
                [r, LB, UB, F, df1, df2, p] = ICC([pars_est1_set(~outliers),pars_est2_set(~outliers)], 'C-k');
                plot(pars_est1_set(~outliers), pars_est2_set(~outliers), 'x');
                
                %             xlim([min([pars_sim_set; pars_est_set]),max([pars_sim_set; pars_est_set]) ])
                %             ylim([min([pars_sim_set; pars_est_set]),max([pars_sim_set; pars_est_set])] )
                %
                title({['corr = ', num2str(r)], ['pval = ', num2str(p)], ['par_i = ', num2str(par_ind)]});
                %     waitforbuttonpress
                %     close all
                
                
            end
            
            % subplot(2,2,6)
            % noize_sim_set  = noize_sim;
            % noize_est_set  = noize_est;
            % [r,p] = corr(noize_sim_set,  noize_est_set);
            %
            % %      figure;
            % plot(noize_sim_set, noize_est_set, 'x');
            %
            % xlim([min([noize_sim_set; noize_est_set]),max([noize_sim_set; noize_est_set]) ])
            % ylim([min([noize_sim_set; noize_est_set]),max([noize_sim_set; noize_est_set])] )
            %
            % title(['Noise, corr = ', num2str(r), ', pval = ', num2str(p)]);
%             ses_name = list_dir(l).name;
%             plot_name = [ses_name(1:end-4),', noise level ',num2str(noize_levels) ]
%             plot_name = strrep(plot_name, '_', ' ');
%             plot_name = strrep(plot_name, '   ', '  ')
%             plot_name = strrep(plot_name, '  ', ', ')
%             suptitle(plot_name);
            
            
        end
        waitforbuttonpress
    end
end

return
%