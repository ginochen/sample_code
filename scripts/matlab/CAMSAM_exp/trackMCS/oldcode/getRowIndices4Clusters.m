function [llcri_it, lonlatcri_it, lonlatsri_it] = getRowIndices4Clusters(t,nt,maxd)
  % Purpose: get the timestep-separated row indices of clusters
  % open two consecutive time maps and compare the centroid Euclidean distance (lon,lat)
  % If the two centroids are close enough between two timesteps (it,it+1), then save
  % the associated row index of mcs.lonlat in the two columns of
  % lonlatcri_it{it}, and all the surrounding mcs.lonlat row indices in lonlatsri_it.
  % Otherwise if no centroid is detected at the later timestep it+1, save 0
  % as the row index in the 2nd column of lonlatcri_it{it}
  % Variables:
  % llcri_it{it}(icluster,1:2): row indices of mcs.llcentroid for each cluster (2-tsteps)
  % lonlatcri_it{it}(icluster,1:2): 
  %   row indices of mcs.lonlat (not the actual longitude & latitude) for the
  %   two centroids (2-tsteps), the 2nd column could be zero
  %   (i.e., no row) if no nearest neighbor centroid is detected
  % lonlatsri_it{it}{icluster_rowset1,icluster_rowset2): 
  %   row indices of mcs.lonlat for the surrounding points of the two centroids
  %   (2-tsteps)
  for it = 1:nt-1
    eval(sprintf('tmp = load(''var_PC1_%d.mat'',''mcs'');',t(it)));
    mcs1 = tmp.mcs;
    eval(sprintf('tmp = load(''var_PC1_%d.mat'',''mcs'');',t(it+1)));
    mcs2 = tmp.mcs;
    for ic1=1:mcs1.ncentroids % llcentroid row index of 1
      [dist mind] = getCentroidMinDist(mcs1,mcs2,ic1); 
      [icr1 isr1] = getRowIndices4Cluster(ic1,mcs1); % centroid and surrounding mcs.lonlat row indices of 1 <--- needed to index varCAM0D (since mcs.lonlat rows are the same as varCAM0D rows)
      icr2=0; isr2=0; % default to not the same cluster for 1 and 2
      if (mind <= maxd) % less than 3deg ~= 300km 
        ic2 = find(dist==mind); % llcentroid row index of 2 (with min dist to 1 centroid)
        if (numel(ic2)>1); % if more than one mcs2-centroid that are close to mcs1-centroid
          ic2 = getCloserCluster(ic2,isr1,mcs1,mcs2);
        end
        if ic2 % with overlapped points, zero represents if cluster terminates at 2nd-step
          [icr2 isr2] = getRowIndices4Cluster(ic2,mcs2);
        end
      end
      llcri_it{it}(ic1,1:2) = [ic1,ic2]; % llcentroid row index associated to 1 and 2 for the same cluster
      lonlatcri_it{it}(ic1,1:2)=[icr1,icr2]; % lonlat centroid row index
      lonlatsri_it{ic1,it}= {isr1,isr2}; % row indices of lonlat for the surrounding points of the two centroids 
    end
  end
