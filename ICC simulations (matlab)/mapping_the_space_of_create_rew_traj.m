
no_trajs= 8000;
AllVarsAll=NaN(no_trajs,9);

ntrials = 240;

no_sim = 50;
trials_per_block = [3,30];%[49,51;49,51];% min max for each of the two volatility stages
%     sd_options = [5,15];
prob_set = [0.2,0.8,0.5];
p_vol = 0.1;
rew_interval = [10,90]; % [-6,6];
sd_options = [5,15];
for i = 1: no_trajs
    if 100*i/no_trajs == round(100*i/no_trajs )
        disp([num2str(100*i/no_trajs), '%']);
    end
    [r,m,s] = create_rew_traj(ntrials, rew_interval, trials_per_block, sd_options);
    



AllVarsAll(i,:) = InspectRewTrajs(r,m,s);

% AllVars = [SDratio, avg_div_from_mean,...
%   avg_block_length, avg_block_length_15, avg_block_length_5,...
% length(reversals), min(reversals), max(reversals),mean(reversals)];
end
VarNames = {'SDratio', 'avg-div-from-mean',...
  'avg-block-length', 'avg-block-length-15', 'avg-block-length-5',...
'no-reversals', 'min-revs', 'max-revs', 'mean-revs'};
threshold = zeros(length(VarNames),2);
for i = 1: length(VarNames)
    
    figure;
    counts = hist(AllVarsAll(:,i),100);
    hist(AllVarsAll(:,i),100);
    hold on
    for s = 1:4
    plot([AllVars(s,i) AllVars(s,i)], [0, max(counts)]);
    end
    title(VarNames{i});
    threshold(i,1) = input('mean = ');
    threshold(i,2) = input('bound = ');
    close
end

save('thresholds.mat','threshold')

threshold = load('thresholds.mat');
%% given the thresholds how hard is it to find a trajectory??
findtarget =0;
while findtarget == 0
    [r,m,s] = create_rew_traj(ntrials, rew_interval, trials_per_block, sd_options);
    
    AllVarsAll_Target = InspectRewTrajs(r,m,s);
    
    if abs(AllVarsAll_Target(1) - threshold(1,1)) < 0.02;
        findtarget = 1;
        
    end
    
end
findproposal =0;
tic
while findproposal ==0;
       [r,m,s] = create_rew_traj(ntrials, rew_interval, trials_per_block, sd_options);
       
       AllVarsAll_Proposal = InspectRewTrajs(r,m,s);
       diff_vect = abs(AllVarsAll_Proposal - AllVarsAll_Target)';
       diff_vect(isnan(threshold(:,2))) = [];
       th = threshold(:,2);
       th(isnan(threshold(:,2))) = [];
       if  sum(diff_vect<th) == length(th)
           findproposal = 1;
       end
end
toc
       
