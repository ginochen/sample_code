function [icr isr] = getRowIndices4Cluster(ir,mcs)
  % Purpose: get the mcs.lonlat row index of centroid (and surroundings) of the cluster
  % icl: 'cl'uster index
  % isr: 's'urrounding mcs.lonlat 'r'ow indices of ic-cluster
  % icr: 'c'entroid mcs.lonlat 'r'ow index of ic-cluster
  icr = find(ismember(mcs.lonlat, mcs.llcentroids(ir,:),'rows'));
  icl = mcs.clusterindex(icr); % cluster index of min distance centroid
  isr = find(mcs.clusterindex==icl);
