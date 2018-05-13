function trackMCS(season)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % THIS PROGRAM (trackMCS.m) CALCULATES AND SAVES THE FOLLOWING VARIABLES
  % run after doing getMCSlonlat.m
  % mcslltC: 4 col (lon-lat-time-Clustind) matrix with samples in the row
  % mcsillt: 3 col (ilon-ilat-itime) matrix
  % mcslltClRowInd(ic,1:nmcs4Cl(ic)): mcsllt row indices (saved in col) associated to each cluster (saved in row)
  % mcsillt4Cl{ic}{iit}(samp,1:2): ilon-ilat samples for ic-cluster at iit-time
  % nt4Cl(ic): # of unique time for ic-cluster
  % nt4ClRowInd(i): row indices assoc with clusters with i-timesteps, 1 timestep = 3 hours
  % nCl4nt(i): total # of clusters with i-timesteps 
  % t4Cl(ic,1:nt4Cl(ic)): unique time indices for ic-cluster
  % illt4Cl{ic}(samp,1:3): ilon-ilat-itime samples for ic-cluster
  % nmcs4Cl(ic): # of mcs in ic-cluster
  % nCl: total # of clusters
  % t: a vector that saves time strings as '0001-02-01-10800' to load the nc files 
  % dt: time distance set to heirarchical clustering analysis cutoff length
  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load lonlat variables and add the time variable
%%%%% PARAMETER %%%%%%%%%%%%%%%%%%%%
  doecon = 1; % integer clustering, euclidean distance are eval in if loop and integer assigned to be actual distance
  ym = containers.Map({'MAM','JJA','SON','DJF'},{[1 3; 1 4; 1 5], [1 6; 1 7; 1 8], [1 9; 1 10; 1 11], [1 12; 2 1; 2 2]});
%  season = {'MAM','JJA','SON','DJF'};
  mfile = '~/scripts/matlab/CAMSAM_exp/trackMCS.m';
  display('make sure to use bigmem queue and allocate at least 70G of mem to calc 3 month of data points')
  run /nethome/gchen/scripts/matlab/startup.m
  spcase = 'F_2000_SPCAM_m2005_3hrly1';
  spHist.dir = ['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/' spcase '/atm/hist/'];
  fname = [spcase '_MCSlonlat'];
  mday = [31 28 31 30 31 30 31 31 30 31 30 31]; % # of days in a month
  % check if .mat file exists
  ii = 1;
  for iy = 1:1:2
    for im = 1:1:12 % start from 0001-03-01-0000 for MAM JJA SON DJF
      for id = 1:1:mday(im)
        for is = 0:10800:75600
          if ismember([iy im],ym(season),'rows') 
            time=num2str(sprintf('%04g-%02g-%02g-%05g',iy,im,id,is));
            if exist([spHist.dir fname '.' time '.mat'],'file')
              t{ii} = time;
              ii=ii+1;
            else
              error([time ' is missing, need it!! exiting ...'])
            end
          end
        end
      end
    end
  end
  if ( ~exist('t') )
    error('no file to track')
  end
  lon = ncread([spHist.dir spcase '.cam.h1.' t{1} '.nc'],'lon');
  lat = ncread([spHist.dir spcase '.cam.h1.' t{1} '.nc'],'lat');
  landfrac = squeeze(ncread([spHist.dir spcase '.cam.h1.' t{1} '.nc'],'LANDFRAC'));
  nt = numel(t);
  dx = diff(lon(1:2));
  dy = diff(lat(1:2));
  dt = sqrt(dx^2+dy^2);
  cutoff = dt + .01; %dy; don't get lower than dy, otherwise 2*dt away points may be artificially qualified as the same cluster
  rb = lon(end); % right boundary 358.25, this is for Re-Cluster step
  lb = lon(1); % left boundary 0
  %cutoff = sqrt(dx^2 + dt^2) + .01; %fig 4
  display(cutoff)
  % dlon = 1.25km, dlat = 0.94241, d = sqrt(1.25^2+0.94241^2) = 1.5654, with time dt =1, d = 1.85... (still less than 0.942*2~1.9)
  % To satisfy as the same 3hrly time cluster in space, there must be overlapped points between two time steps.
  % Here we relax to without requirement of overlapping, but simply one gridbox neighboring.
  % Therefore, the diagonal neighboring distance (1.565) is the furthest between two grid points in space
  % and furthest diagonal neighboring distance 1.85 in time
  % Hence, d = 2*1.565 ~ 3.0 or 2*0.94241 ~ 1.88 or 2*1.25 ~ 2.5 are the possible unqualified distance, hence if we define dt=1, and d=1.85
  % --------------------------------------
  %    |      |     |      |    |    |
  %  0.9424  1.25  1.85  1.88  2.5  3.0
  %  
  % where 1.85<1.88, hence, as long as cutoff not greater than 1.88, it's good 
  % to be considered the same time cluster
  %cutoff=2.8; % used for 'euclidean' 'single' distance
%%%%% Merge all timesteps' mcs.lonlat into one matrix %%%%%%%%%%%%%%%%
  mcsllt = []; mcsillt=[];
  mcst(1:nt) = dt*(1:nt); % time values of mcs at each time index
  display('start concatenating mcs') 
  for it = 1:nt
    eval(sprintf('load(''%s/%s.%s.mat'',''mcs'');',spHist.dir,fname,t{it}));
    mcsn(it) = size(mcs.lonlat,1); % # of mcs globally at each time 
    mcsllt  = cat(1,mcsllt, [[mcs.lonlat(:,1),mcs.lonlat(:,2)],  ones(size(mcs.lonlat,1),1)*mcst(it)]); % lon-lat-time of scattered cloud points saved into a 2D matrix mcsllt
    mcsillt = cat(1,mcsillt,[[mcs.ilonlat(:,1),mcs.ilonlat(:,2)],ones(size(mcs.lonlat,1),1)*it]); % ilon-ilat-itime of scattered cloud points saved into a 2D matrix mcsllt
  end
  display('end concatenating mcs')
  mcsN = size(mcsllt,1); % total number of LS MCS points in time and space
%%%%% Cluster MCSs in space and time %%%%%%%%%%%%
  display('start clustering MCSs') % try to break it down to clustering by months and then glue together
%  mcslltC = [];
%  tvec = {[1,sum(mcsn(1:35))],[sum(mcsn(1:25)),sum(mcsn(1:65))],[sum(mcsn(1:55)),sum(mcsn)]};
%  for im = 1:numel(tvec)
%    [C, ~] = myclusterdata(mcsllt(tvec{im}(1):tvec{im}(2),:),'euclidean','single','distance',cutoff,2); % C: cluster index array; Z: pairwise-linkage-distance matrix
  ny = size(mcsllt,1);
  ny = ny*(ny-1)/2;
  disp(sprintf('%16.f',ny));
  Y = zeros(1,ny,'single');
  [C, ~] = myclusterdata(mcsllt,Y,'euclidean','single','distance',cutoff,doecon); % C: cluster index array; Z: pairwise-linkage-distance matrix
%    maxC(im) = max(C); % save the maximum cluster index and add to the next month
  mcslltC = [mcsllt, C]; % cluster-labeled mcsllt
%  end
  display('finish clustering MCSs')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [mcslltC, C] = reCluster(mcsllt, mcst, C, cutoff, lb, rb, nt)
%   Re-Cluster MCS clusters with neighboring points at the periodic bdry
%   mcslltC: mcsllt with merged cluster index added to the 4-th column
%   C:       just the last column of mcslltC
  % find the same time-valued 360 rows and 
  for it = 1:nt
    tftmp = mcsillt(:,3)==it; % rows of the it-th time 
    itmp1 = find(tftmp & mcslltC(:,1)==rb); % find the right bdy rows of mcslltC
    mcslltC(itmp1,1)=rb-360; % translate rb to the closest value of lb, which is a neg value
    itmp2 = find(tftmp & mcslltC(:,1)==lb); % find the left  bdy rows of mcslltC
    % element-wise compare if they are neighbors
    for i = 1:numel(itmp1)
      for j = 1:numel(itmp2)
        if ( (pdist([mcslltC(itmp1(i),1:2);mcslltC(itmp2(j),1:2)]) < cutoff) & (diff(mcslltC([itmp1(i),itmp2(j)],4))) ) % if the two bdry points are neighbors, and different cluster
          Cs = sort([mcslltC(itmp1(i),4),mcslltC(itmp2(j),4)]); % sort the two cluster index (c-ind) by magnitude
          mcslltC(mcslltC(:,4)==Cs(2),4) = Cs(1); % assign the larger c-ind with smaller c-ind, essentially merging and reducing the number of clusters by one
        end
      end
    end
  end
  C = mcslltC(:,4); % new cluster index set after Re-Cluster
  Cu = unique(sort(C)); % unique cluster indices, e.g., unique(sort([1 4 3 4])) = [1 3 4]
  nCl = numel(Cu); % total clusters 
  island = zeros(nCl,1);  % 0 for 'o'cean
  for ic = 1:nCl% index of clusters
tic
    nmcs4Cl(ic) = sum(C==Cu(ic)); % # of mcs in a Cl
    mcslltClRowInd(ic,1:nmcs4Cl(ic)) = find(C==Cu(ic)); % mcsllt row indices associated to each cluster
    illt4Cl{ic}(:,1:3) = mcsillt(nonzeros(mcslltClRowInd(ic,:)),1:3); % the ilon-ilat-itime matrix for ic-th cluster
    indtmp = sub2ind(size(landfrac),illt4Cl{ic}(:,1),illt4Cl{ic}(:,2)); % lon-lat vector index for mcs
    if ( any(landfrac(indtmp) > 0 & landfrac(indtmp) < 1) )
      island(ic) = 2; % 2 for coastal 'l'and 
    end
    if ( all(landfrac(indtmp) == 1) ) % pure land swirling MCS
      island(ic) = 1; % 1 for 'l'and
    end
    tmp = unique(illt4Cl{ic}(:,3)); % find the unique time coord
    nt4Cl(ic) = numel(tmp);
    t4Cl(ic,1:nt4Cl(ic)) = tmp; 
    for iit=1:nt4Cl(ic) %loop over each time for a cluster, and find the centroid at each time
      mcsillt4Cl{ic}{iit}(:,1:3) = illt4Cl{ic}(find(illt4Cl{ic}(:,3)==t4Cl(ic,iit)),1:3); % mcsillt4Cl{ic}{iit} gives the ilon-ilat indices at a fixed time of a cluster
      mcsnll(ic,iit) = size(mcsillt4Cl{ic}{iit},1); % size (total # of points) for ic-cluster at iit-time
    end
    mcsnllave(ic) = mean(mcsnll(ic,1:nt4Cl(ic))); % time average size for ic-cluster  
toc
  end
  t4Cl = sparse(t4Cl);
  mcslltClRowInd = sparse(mcslltClRowInd); % save as sparse matrix 
  for i=1:max(nt4Cl); 
    nt4ClRowInd{i} = find(nt4Cl==i); % find all clusters with nt = i
    nCl4nt(i) = numel(nt4ClRowInd{i}); % # of clusters with nt = i 
  end
  ilatlim = [min(mcsillt(:,2)),max(mcsillt(:,2))];
  latlim = [min(mcslltC(:,2)), max(mcslltC(:,2))];
  
  save([spHist.dir season '/mcs_clusters.mat'],'mfile','mcslltC','mcsillt','mcsnll','mcsnllave','mcslltClRowInd','mcsillt4Cl','nt4Cl','nt4ClRowInd','nCl4nt','t4Cl','illt4Cl','ilatlim','island','nmcs4Cl','nCl','lon','lat','latlim','t','dt');
  % mcslltC: 4 col (lon-lat-time-Clustind) matrix with samples in the row
  % mcsillt: 3 col (ilon-ilat-itime) matrix
  % mcslltClRowInd(ic,1:nmcs4Cl(ic)): mcsllt row indices (saved in col) associated to each cluster (saved in row)
  % mcsillt4Cl{ic}{iit}(samp,1:2): ilon-ilat samples for ic-cluster at iit-time
  % nt4Cl(ic): # of unique time for ic-cluster
  % nt4ClRowInd(i): row indices assoc with clusters with i-timesteps, 1 timestep = 3 hours
  % nCl4nt(i): total # of clusters with i-timesteps 
  % t4Cl(ic,1:nt4Cl(ic)): unique time indices for ic-cluster
  % illt4Cl{ic}(samp,1:3): ilon-ilat-itime samples for ic-cluster
  % nmcs4Cl(ic): # of mcs in ic-cluster
  % nCl: total # of clusters
  % t: a vector that saves time strings as '0001-02-01-10800' to load the nc files 
  % dt: time distance set to heirarchical clustering analysis cutoff length
  display('run cluster_centroid.m after this code')
