load /Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly2/atm/hist/lat2525/DJF/mcs_cluster
% percentage of large cluster over ocean/land (size defined by max lifetime cell number)
ilnd = 1; % land (1) or ocean (0)
sizethres = 2; % threshold for size
c_size = mcsnll(island==ilnd,:);
c_msize = max(c_size,[],2); % max size for each cluster
nclarge_ratio = sum(c_msize>sizethres)/numel(c_msize)
ncsmall_ratio = sum(c_msize<=sizethres)/numel(c_msize)




% average cluster grid cell coverage 
mean(c_size(c_size~=ilnd))
