function [u_traj, prob_traj] = create_rew_traj_binary(ntrials, prob_set, trials_per_block, p_vol)

if nargin < 4
    p_vol = 0.1;
end
% u_traj = zeros(ntrials,1);

% prob_set = [0.2,0.8,0.5];
% Blocks = randi(trials_per_block,[no_blocks-1,1]);

Blocks = geornd(p_vol, [100, 1])+ trials_per_block(1);
Blocks = Blocks(1:find(cumsum(Blocks)<= ntrials, 1,'last'));
if sum(Blocks) ~= ntrials
    
    Blocks = [Blocks; ntrials - sum(Blocks)];
    
end
no_blocks = length(Blocks);

u_traj = NaN(ntrials,1);
prob_traj =  NaN(ntrials,1);
RevTrials = [0;cumsum(Blocks)];

% p_set = randi([0,1], [no_blocks,1]);
p_set = randi(length(prob_set),  [no_blocks,1]) ;
p_set = prob_set(p_set);

for b = 1 : length(Blocks)
    prob = p_set(b);
   
    u_traj_block = binornd(1, prob*ones(Blocks(b),1));
    prob_traj_block = prob*ones(Blocks(b),1);
    %rew_traj_block(rew_traj_block > 5) = 5 + (5 - rew_traj_block(rew_traj_block > 5));
    %rew_traj_block(rew_traj_block < -4) = -4 + (-4 - rew_traj_block(rew_traj_block < -4));
    
    u_traj(RevTrials(b)+1: RevTrials(b+1)) = u_traj_block ;
      prob_traj(RevTrials(b)+1: RevTrials(b+1)) = prob_traj_block;
      
end

