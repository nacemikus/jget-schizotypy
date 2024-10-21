function [rew_traj, mu_traj, sd_traj] = create_rew_traj(ntrials, rew_interval, trials_per_block, sd_options, p_vol, sd_vol)
if nargin < 5
    p_vol = 0.1;
end
if nargin < 6 
    sd_vol = 0.4;
end

mu_traj = zeros(ntrials,1);
sd_traj = zeros(ntrials,1);

% Blocks = randi(trials_per_block,[no_blocks-1,1]);

Blocks = geornd(p_vol, [100, 1])+ trials_per_block(1);
Blocks = Blocks(1:find(cumsum(Blocks)<= ntrials, 1,'last'));
if sum(Blocks) ~= ntrials
    
    Blocks = [Blocks; ntrials - sum(Blocks)];
    
end
no_blocks = length(Blocks);
% Blocks2 = randi(trials_per_block,[no_blocks,1]);


% 
% 
% while sum(Blocks) ~= ntrials
%     Blocks = randi(trials_per_block,[no_blocks,1]);
% end
% sd_set = ones(length(Blocks),1);
% sd_set(1:round(length(Blocks)/2)) = 0;
% sd_set = sd_set(randperm(length(sd_set)));
% sd_set(sd_set==0) = sd_options(1);
% sd_set(sd_set==1) = sd_options(2);

rew_traj = NaN(ntrials,1);

RevTrials = [0;cumsum(Blocks)];
mu_set = NaN(length(Blocks),1);
% sd_set = NaN(length(Blocks1),1);

mu_set = randi(rew_interval, [no_blocks,1]);
 sd_opt = randi(length(sd_options));
 sd =  sd_options(sd_opt);
for b = 1 : length(Blocks)
    mu = mu_set(b);
   
    
    rew_traj_block = round(randn([Blocks(b),1])*sd + mu);
    %rew_traj_block(rew_traj_block > 5) = 5 + (5 - rew_traj_block(rew_traj_block > 5));
    %rew_traj_block(rew_traj_block < -4) = -4 + (-4 - rew_traj_block(rew_traj_block < -4));
    rew_traj(RevTrials(b)+1: RevTrials(b+1)) = rew_traj_block ;
    mu_traj(RevTrials(b)+1: RevTrials(b+1)) = repmat(mu ,[Blocks(b),1]);
    sd_traj(RevTrials(b)+1: RevTrials(b+1)) = repmat(sd ,[Blocks(b),1]);
    
    if (rand(1) < sd_vol) && (length(sd_options) >1)
        if sd_opt == 1
            sd_opt= 2;
        else
            sd_opt =1;
        end
    
        sd =  sd_options(sd_opt);
    end
    
%     mu_new = mu;
%     while abs(mu_new - mu) < sd_set(b)
%         
%         mu_new = randi(rew_interval);
%     end
%     mu = mu_new;
%     sd = sd_set(randi([1,length(sd_set)]));
%   
%      sd = sd_set(randi([1,length(sd_set)]));
end

rew_traj(rew_traj<1) = 1;
rew_traj(rew_traj>100)= 100;

