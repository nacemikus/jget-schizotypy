function AllVars = InspectRewTrajs(rew,mu,sd)

% inspecting reward trajectories to compare them later:

jumps  = diff(mu);

block_length_temp = double(jumps == 0);
block_length = (diff([0;find(block_length_temp==0);240]));

avg_block_length = mean(block_length);
sd_all = sd(find(block_length_temp==0));
avg_block_length_15 = mean(block_length(sd_all==15));
avg_block_length_5 = mean(block_length(sd_all==5));


div_from_mean = abs(mu - rew);
avg_div_from_mean = mean(div_from_mean);
min_div_from_mean = min(div_from_mean);
max_div_from_mean= max(div_from_mean);

SDratio = sum(sd == 15)/240;

jumps(jumps==0) = [];
reversals = abs(jumps);
AllVars = [SDratio, avg_div_from_mean,...
  avg_block_length, avg_block_length_15, avg_block_length_5,...
length(reversals), min(reversals), max(reversals),mean(reversals)];


end