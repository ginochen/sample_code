%%%% MCS / Cloud Cluster Centroid Tracking Algo %%%%%%%%%%%%%%
% Purpose: Obtain timeseries row indices (icr and isr) for the centroids of the
%          same clusters. Use the row indices to track a cluster's evolution
% Caution: Apply this code when varCAM0D-1D-2D variables are the same # of
%          lon-lat points as mcs.lonlat
%%%%%%%% PARAMETERS %%%%%%%%%%%%%
cdt = 10; % track cluster timestep gap
mint = 1; % minimum timesteps to qualify as mcs 3 steps = 3*5hrs = 15hrs 
maxd = 5; % max dist to remain qualified as the same cluster (between two timesteps)
% according to Corfidi 2003, max propagation speed is around 60 knot ~ 30 m/s ~ 108 km/hr
t = [11:cdt:401]; %:10:401
nt= numel(t);
if (nt <= mint); error(['specify more than ' num2str(mint) ' steps!']); end

[llcri_it, lonlatcri_it, lonlatsri_it] = getRowIndices4Clusters(t,nt,maxd);
[llcri, lonlatcri, lonlatsri]       = combineRowIndices4clusters(llcri_it,lonlatcri_it,lonlatsri_it,nt);

%%%%%%%%% Get Var with Row Indices %%%%%%%%%%%%%%%%
clt = sum((llcri~=0),2); % # of times in a cluster 
cli = find(clt>=mint)'; % index of clusters that are greater than mint times 
ntClust = clt(cli); % the clt that is actually associated to the saved varCAM0D_mcs
if isempty(cli); error(['no clusters are greater than steps = ' num2str(mint)]); end
ii=1;
[varCAM0D_mcs, varCRM0D_mcs, varCAM1D_mcs, varCRM1D_mcs, varCRM2D_mcs, lonlatmcs] = getClusterVar(cli,ntClust,lonlatsri,lonlatcri,llcri,t);

save('var_PC1_MCSt','varCAM0D_mcs','varCAM1D_mcs',...
     'varCRM0D_mcs','varCRM1D_mcs','varCRM2D_mcs',...
     't', 'maxd', 'mint','llcri','llcri_it',...
     'lonlatsri','lonlatsri_it','lonlatcri','lonlatcri_it','cli','clt','ntClust','lonlatmcs')


function [varCAM0D_mcs, varCRM0D_mcs, varCAM1D_mcs, varCRM1D_mcs, varCRM2D_mcs, lonlatmcs] = getClusterVar(cli,ntClust,lonlatsri,lonlatcri,llcri,t)
  ii=1;
  for ic = 1:numel(cli)
    iit = find(llcri(cli(ic),:)~=0); % llcri(icluster,it) saves the mcs.llcentroid row indices as non-zero values, so finding the non-zero columns of llcri is finding the timesteps assoc to a cluster 
    for it = 1:ntClust(ic)
      eval(sprintf('load(''var_PC1_%d.mat'',''varCAM0D'',''varCRM0D'',''varCAM1D'',''varCRM1D'',''varCRM2D'',''mcs'');',t(iit(it))));
      illr = lonlatsri{cli(ic),iit(it)}; % lonlat row indices in a cluster at a timestep
      lonlatmcs.centroid{ic}(:,it) = mcs.lonlat(lonlatcri(cli(ic),iit(it)),:); % save the centroid lon lat time-series for each cluster separately
      lonlatmcs.npts(ic,it) = numel(illr); % # of lon-lat points at a given tstep for a cluster
      for is = 1:lonlatmcs.npts(ic,it) % surroundings' row indices counter
        ir = illr(is); % lonlat row for surroundings' 
        % Analyze at every cluster/timestep 
        varCAM0D_mcs{ic}{it}(:,is)          = varCAM0D(ir,:);
        varCRM0D_mcs{ic}{it}(:,is)          = varCRM0D(ir,:);
        varCAM1D_mcs{ic}{it}(:,:,is)        = varCAM1D(:,ir,:);
        varCRM1D_mcs{ic}{it}(:,:,is)        = varCRM1D(:,ir,:);
        varCRM2D_mcs{ic}{it}(:,:,:,is)      = varCRM2D(:,:,ir,:);
        lonlatmcs.surrounding{ic}{it}(:,is) = mcs.lonlat(ir,:); % save the surrounding lon lat time-series for each cluster separately
        ii=ii+1;
      end
    end
  end
end



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
      if (mind <= maxd) % less than 3deg ~= 300km 
        ic2 = find(dist==mind); % llcentroid row index of 2 (with min dist to 1 centroid)
        if (numel(ic2)>1); % more than one mcs2-centroid that are close to mcs1-centroid
          ic2 = getCloserCluster(ic2,isr1,mcs1,mcs2); 
        end
        [icr2 isr2] = getRowIndices4Cluster(ic2,mcs2);
      else
        ic2 = 0; % zero represents if cluster terminates at 2nd-step
        icr2 = 0;
        isr2 = {};
      end
      llcri_it{it}(ic1,1:2) = [ic1,ic2]; % llcentroid row index associated to 1 and 2 for the same cluster
      lonlatcri_it{it}(ic1,1:2)=[icr1,icr2]; % lonlat centroid row index
      lonlatsri_it{ic1,it}= {isr1,isr2}; % row indices of lonlat for the surrounding points of the two centroids 
    end
  end
end

function ic2 = getCloserCluster(ic2,isr1,mcs1,mcs2); 
  % Purpose: get the closer cluster to mcs1 if there are two centroids in mcs2 that are the same min dist
  % Method: use 'pdist2' to eval the mean Euclidean dist of all surrounding points of mcs1 and mcs2
  for iii=1:numel(ic2)
    [icr2 isr2] = getRowIndices4Cluster(ic2(iii),mcs2);
    meanDist(iii) =mean(mean(pdist2(mcs1.lonlat(isr1,:),mcs2.lonlat(isr2,:))));
  end 
  ic2 = ic2(find(min(meanDist)));% eval which mean distance is closer
end

function [icr isr] = getRowIndices4Cluster(ir,mcs)
  % Purpose: get the mcs.lonlat row index of centroid (and surroundings) of the cluster
  % icl: 'cl'uster index
  % isr: 's'urrounding mcs.lonlat 'r'ow indices of ic-cluster
  % icr: 'c'entroid mcs.lonlat 'r'ow index of ic-cluster
  icr = find(ismember(mcs.lonlat, mcs.llcentroids(ir,:),'rows'));
  icl = mcs.clusterindex(icr); % cluster index of min distance centroid
  isr = find(mcs.clusterindex==icl);
end

function [dist mind] = getCentroidMinDist(mcs1,mcs2,ic1)
  % Purpose: get the minimum dist between mcs1 (ic1-th centroid) and mcs2 (all centroids)
  dist = sqrt((mcs1.llcentroids(ic1,1)-mcs2.llcentroids(:,1)).^2+...
         (mcs1.llcentroids(ic1,2)-mcs2.llcentroids(:,2)).^2);
  mind = min(dist);
end


function [llcri, lonlatcri, lonlatsri] = combineRowIndices4clusters(llcri_it,lonlatcri_it,lonlatsri_it,nt)
  % Purpose: combine the timestep-separated row indices into a single matrix to track a cluster
  % llcri(icluster,it): combined matrix of llcentroid row indices for each cluster (timeseries)
  % lonlatcri(icluster,it): combined matrix of mcs.lonlat centroid row indices for each cluster (timeseries)
  iic=1; % row counter for clusters
  for it = 1:nt-2 % only nt-1 llcri_it's, hence nt-1-1
    if it==1
      ind = 1:size(llcri_it{1},1); % use all row indices for the 41-tstep 
    else
      ind = setdiff([llcri_it{it}(:,1);0],llcri_it{it-1}(:,2))'; % starting a new cluster if the rows aren't used previously
    end
    for i = ind % index for all it=1 clusters
      if (i & llcri_it{it}(i,2)~=0)
        itt=1;
        llcri(iic,it:it+1) = llcri_it{it}(i,1:2);
        lonlatcri(iic,it:it+1) = lonlatcri_it{it}(i,1:2);
        lonlatsri{iic,it} = lonlatsri_it{i,it}{1};
        lonlatsri{iic,it+1} = lonlatsri_it{i,it}{2};
        while ( llcri_it{it+itt}(llcri(iic,it+itt),2)~=0 ) % if eq zero then go to the next i in ind
          it0 = it+itt; it1 = it+itt+1;
          llcri(iic,it1)     = llcri_it{it0}(llcri(iic,it0),2);
          lonlatcri(iic,it1) = lonlatcri_it{it0}(llcri(iic,it0),2);
          lonlatsri{iic,it1} = lonlatsri_it{llcri(iic,it0),it0}{2};
          itt=itt+1; 
          if (it+itt==nt); break; end
        end
        iic=iic+1;
      end
    end
  end
end
