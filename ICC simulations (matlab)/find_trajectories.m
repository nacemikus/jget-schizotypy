%% find trajectories using 

function u_ibo = find_trajectories(ntrials, rew_interval, trials_per_block, sd_options, p_vol, sd_vol)



no_seq = 2;


u_ibo = NaN(ntrials, no_seq);

mds = ones(1,2)*NaN;
disp(' ')
disp('Mahalanobis distances:')
disp(mds)

%% find MD treshhold
MD_set = [];
if (exist(['MahalanobisD_dist_set_trials_',num2str(ntrials),'.mat']) == 0 )
    no_sims = 100;
else 
    no_sims = 1;

end
for sims = 1 : no_sims
    pars = [];
    precs = [];
    % Initialize priors
    priorpar = [];
    priorprec = [];
    % Initialize indices of optimized parameters
    opt_idx = [];
    for i = 1:2
        u_IBO_start = create_rew_traj(ntrials, rew_interval, trials_per_block, sd_options, p_vol, sd_vol);
        u_ibo(:,i)= u_IBO_start;
        %
        estuse = tapas_fitModel([],...
            u_IBO_start,...
            'tapas_ehgf_jget_config',...
            'tapas_bayes_optimal_config');
        opt_idx = estuse.c_prc.priorsas;
        opt_idx(isnan(opt_idx)) = 0;
        opt_idx = find(opt_idx);
        
        pars = [pars, estuse.p_prc.ptrans(opt_idx)'];
        precs = cat(3, precs, estuse.optim.H);
        priorpar = estuse.c_prc.priormus(opt_idx)';
        priorprec = diag(estuse.c_prc.priorsas(opt_idx));
    end
    
    mds(1) = md(pars(:,1), pars(:,2), precs(:,:,2));
    mds(2) = md(pars(:,2), pars(:,1), precs(:,:,1));
    MD_set =  [MD_set, mds];
end
if (exist(['MahalanobisD_dist_set_trials_',num2str(ntrials),'.mat']) == 0 )
    save(['MahalanobisD_dist_set_trials_',num2str(ntrials),'.mat'], 'MD_set');
end



%%
load(['MahalanobisD_dist_set_trials_',num2str(ntrials),'.mat'])

Md_th = quantile(MD_set, 0.05);%0.05;

while max(mds(1:2)) > Md_th %Chris?s settings were mean(mds) > 0.05 || max(mds) > 0.02
    % Loop through sequences trying to get their parameter values closer to the
    % mean of the values of the other sequences
    
    
    for i = 1:2
        %         count_reps(i) = count_reps(i)  + 1;
        % Remove this parameter vector from the parameter matrix
        parsred = pars;
        parsred(:,i) = [];
        % Remove this precision matrix from the precision array
        precsred = precs;
        precsred(:,:,i) = [];
        u_IBO_proposal = create_rew_traj(ntrials, rew_interval, trials_per_block, sd_options, p_vol, sd_vol);
        %
        estuse = tapas_fitModel([],...
            u_IBO_proposal,...
            'tapas_ehgf_jget_config',...
            'tapas_bayes_optimal_config');
        opt_idx = estuse.c_prc.priorsas;
        opt_idx(isnan(opt_idx)) = 0;
        opt_idx = find(opt_idx);
        
        
        newpar = estuse.p_prc.ptrans(opt_idx)';
        %             disp('newmd?');
        newmd = md(newpar, parsred, precsred);
        %             disp('newmd worked');
        
        
        oldmd = mds(i);
        %             mds(i) = oldmd;
        if newmd < oldmd
            %             fits_criteria
            %             if fits_criteria
            pars(:,i) = newpar;
            precs(:,:,i) = estuse.optim.H;
            u_ibo(:,i)= u_IBO_proposal;
            
            mds(1) = md(pars(:,1), pars(:,2), precs(:,:,2));
            mds(2) = md(pars(:,2), pars(:,1), precs(:,:,1));
            
            
            %             end
            
        end
        
        %         catch
        %
        %             disp(['failed for seq ', num2str(i)]);
        %         end
        
        
        %         if count_reps == stp
        %             count_reps
        %         end
    end
end





end