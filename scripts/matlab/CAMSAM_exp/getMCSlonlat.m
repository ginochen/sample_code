function mcs = getMCSlonlat(fin,fout,season) 
  display('Specify fin as filename without the path, file path should be defined in this function manualy')
  tic
  %%%%%%%%%%%%%%%%%%%%%
  %      Paths        %
  %%%%%%%%%%%%%%%%%%%%%
%  time = fin(end-18:end-3); % get the time of nc file
%  compset     = 'm2005';
%  casei      = 'F_2000_SPCAM_m2005_3hrly_f09f09_branch' %'F_2000_SPCAM_m2005_3hrly1';
%  spArchive   ='/glade2/scratch2/ginochen/archive/matlab/'
%  spArchive   = ['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/'];
%  spHist.dir  = [spArchive casei '/atm/hist/'];
  %tmp         = dir([spHist.dir casei '.cam.h1.*-*-*-*.nc']); % list filenames in a dir and save it to a structure variable  
  %tmp = {tmp.name}; 
  %spHist.name = tmp{it};
%  spHist.name = [casei '.' time];
%  fout = [spHist.dir season '/MCSlonlat/' casei '_MCSlonlat.' time '.mat']; display(fout)
  disp(['input ncfile: ' fin])
  disp(['output mat: ' fout])
  if ~exist(fin,'file')%[spHist.dir season '/MCSlonlat/' casei '_MCSlonlat.' time '.mat'])
    error('input file non-exist, exiting...')
  end
  if exist(fout,'file')%[spHist.dir season '/MCSlonlat/' casei '_MCSlonlat.' time '.mat'])
    error('output file exists, exiting...')
  end
%  fin      = [spHist.dir spHist.name];
  display(fin)
  %%%%%%%%%%%%%%%%%%%%%%%%
  %    Load dimensions   %
  %%%%%%%%%%%%%%%%%%%%%%%%
  ncid = netcdf.open(fin,'NC_NOWRITE');
  dimid = netcdf.inqDimID(ncid,'crm_nx');
  [~, crm_nx] = netcdf.inqDim(ncid,dimid);
  p0   = ncread(fin,'P0'); % [Pa] in SI units not [hPa], can only use spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-03600, the other timesteps P0, PHIS are all zero
  phis = ncread(fin,'PHIS'); % [Pa] in SI units not [hPa], can only use spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-03600, the other timesteps P0, PHIS are all zero
  hyam = flipud(ncread(fin,'hyam'));
  hybm = flipud(ncread(fin,'hybm')); 
  hyai = flipud(ncread(fin,'hyai'));
  hybi = flipud(ncread(fin,'hybi'));
  lat  = ncread(fin,'lat');
  ilatlim = find(lat>=-60 & lat<=60)'; % tropical ilat
%  ilatlim = find(lat>=-25 & lat<=25)'; % tropical ilat
%  ilatlim = find(lat>=20 & lat<=50)' % north midlat ilat
%  ilatlim = find(lat>=0 & lat<=55)' % north midlat ilat
  lon  = ncread(fin,'lon');
  lev  = ncread(fin,'lev');  % level centers
  ilev = ncread(fin,'ilev'); % level interfaces, two interfaces surround one level center
                                  % index 1 (TOA pressure) to nlev (surface pressure) 
  qv_crm   = squeeze(ncread(fin,'CRM_QV')); % (lon, lat, nx, 1, nz) 
%  qc_crm   = squeeze(ncread(fin,'CRM_QC')); % (lon, lat, nx, 1, nz) 
  qi_crm   = squeeze(ncread(fin,'CRM_QI'));
%  w_crm    = squeeze(ncread(fin,'CRM_W'));
  prec_crm = ncread(fin,'CRM_PREC')*3.6e6;
  T_crm    = squeeze(ncread(fin,'CRM_T')); % (lon, lat, nx, 1, nz) 
  T_cam      = flipdim(ncread(fin,'T'),3); % temp in CAM [K]
  display('start MCS search')
  %%%%%%%%% Start MCS criteria search %%%%%%%%%%%%%%%%%%%%%
  illct = 1; % counter for the lon lat indices with zonal lons over ocn only
  minrc = 25;%11; % only if more than 4km*25 = 100km %12*4km=48km (approx pi*24^2 = 1.8e3km^2) of high cloud can be defined as HCS
  minmcs = 25;
  addix = 3; % additional index for the two boundaries
  ibd(1:2) = [addix+1,addix+crm_nx]; % the actual bdy index w.r.t locs 
  for ilon = 1:numel(lon)
    for ilat = ilatlim
      % decide if cloud-top temperature < 241K and if LS cloud fraction = 1
      %[~, pm, ~, ~]  = hybrid2p(p0, ps(ilon,ilat), hyai, hybi, hyam, hybm, nz);
      %%%%%% find high cloud zonal points and the associated cloud-top height %%%%%%%%%%
      leave = 0;
      ictop=NaN(1,crm_nx); tfcold=NaN(1,crm_nx); % predefine no cold
      for ix = 1:crm_nx
        tfcld = squeeze(qi_crm(ilon,ilat,ix,:)) > 1e-6;%| squeeze(qc_crm(ilon,ilat,ix,:)) > 1e-6; % true-false logic variable of cloud 1e-6kg/kg = 0.001g/kg according to Eun-Kyoung Seo and Guosheng Liu
% Yang and Houze 1995 seems favoring 0.5g/kg of nonprecipitating hydrometeor as cloud top
        if ( any(T_crm(ilon,ilat,ix,tfcld)<=260) ) % if high cloud
          ictop(ix) = max(find(tfcld)); % save cloud-top vertical index for ix point
          if (T_crm(ilon,ilat,ix,ictop(ix))<220) % high cloud with very cold top
%          if any(T_crm(ilon,ilat,ix,tfcld)<220) % high cloud with very cold top
             tfcold(ix) = 1;
          end
        end
      end
      ixhcl  = find(~isnan(ictop)); % zonal index of high cloud   % Ex. suppose ixhcl = {1 4 5 7 8 9 31 32} => {{1},{4,5},{7,8,9},{31,32}} => {{1,31,32},{4,5},{7,8,9}}
      nxhcl = numel(ixhcl); % total zonal index of high cloud
      %%%%%% group the zonal points into high cloud complexes (HCC) %%%%%%%%%%
%      ixpf = find(squeeze(prec_crm(ilon,ilat,:))>=1.0)'; % precipitation feature points
      ixpf = find(squeeze(prec_crm(ilon,ilat,:))>=0.1)'; % precipitation feature points
      if (nxhcl<minmcs | numel(ixpf)<minrc); continue; end 
      [ixhcc nxhcc] = continuousSet(ixhcl,nxhcl,1,crm_nx); % HCC set
      [ixhcc nxhcc] = mergeboundary(ixhcc,nxhcc,crm_nx); % HCC set
      for ihcc=1:numel(ixhcc) % go through HCC sets
        if nxhcc(ihcc)< minmcs; continue; end % subset with number of points less than minrc, skip
        Tctop = []; % cloud top temp
        for j=1:nxhcc(ihcc) % go through ith HCC ix indices
          Tctop(j) = T_crm(ilon,ilat,ixhcc{ihcc}(j),ictop(ixhcc{ihcc}(j)));
        end
%        if numel(findpeaks( -Tctop,1:nxhcc(ihcc),'MinPeakwidth',1,'MinPeakProminence',10 ))>1 % if multiple numbers of minimum cold tops
%          [pks locs] = findpeaks(Tctop,1:nxhcc(ihcc),'MinPeakwidth',1,'MinPeakProminence',10); % find the maximum points between cold tops
        if numel(findpeaks( -Tctop ))>1 % if multiple numbers of minimum cold tops
%          [pks locs] = findpeaks(Tctop,1:nxhcc(ihcc),'MinPeakwidth',2, 'MinPeakProminence',5); % find the maximum points between cold tops
          if nxhcc(ihcc)==crm_nx % if the entire 2D domain is being searched, add more points to the boundary points
            [~, locs] = findpeaks([Tctop(crm_nx-addix+1:crm_nx),Tctop,Tctop(1:addix)],1:crm_nx+2*addix,'MinPeakwidth',2, 'MinPeakProminence',10); % find the maximum points between cold tops, peak width intersects the half prominence
            if numel(locs(locs>=ibd(1) & locs<=ibd(2)))<=1 % 0 points or 1 point
              locs = [1,crm_nx]; % include both boundaries indices as two peaks so the peak temp to peak temp will represent an HCS with a trough temp in btw
            elseif any(ismember(locs,[ibd(1),ibd(2)])) % >=2 points and including bdy point(s)
              locs = [1, locs(locs>ibd(1) & locs<ibd(2))-addix, crm_nx];
            else % >=2 points and excluding bdy point(s)
              locs = locs(locs>ibd(1) & locs<ibd(2))-addix;
            end
          else
            [~, locs] = findpeaks(Tctop,1:nxhcc(ihcc),'MinPeakwidth',2, 'MinPeakProminence',10); 
            locs = [1,locs,nxhcc(ihcc)];
          end
          %plot(Tctop); hold on
          %plot(locs,Tctop(locs),'*');pause; hold off
          ixhcs=[]; nxhcs=[];
          for k=1:numel(locs)-1
            ixhcs{k} = ixhcc{ihcc}(locs(k):locs(k+1)); % subsetting the high cloud complex index according to the peak indices
            nxhcs(k) = numel(ixhcs{k});
          end
          if ~isempty(setdiff(ixhcc{ihcc},[ixhcs{:}])) % suppose nxhcc(ihcc)==crm_nx, and >2 points excluding bdy points, 
                                                       % then the rest of the points should be considered as another continuous set
            ixhcs{k+1} = setdiff(ixhcc{ihcc},[ixhcs{:}]);
            nxhcs(k+1) = numel(ixhcs{k+1});
%            [ixhcs{k+1} nxhcs(k+1)] = continuousSet(ixhcs{k+1},numel(ixhcs{k+1}),crm_nx);
          end
        else
          ixhcs = ixhcc(ihcc); % index set for ihcc-th hcs
          nxhcs = nxhcc(ihcc);
        end
        for k = 1:numel(ixhcs)

      %%%%%% find at least one MCS out of the HCSs, stop looping once it's found %%%%%%%%%%
%            if numel(ixpf)>=minrc
%              for i = 1:numel(ixhcc)% total # of HCCs
%                if nxhcs(k) < minrc; continue; end
                if nxhcs(k) < minmcs; continue; end
%                irain = intersect(ixhcc{i}, (find(squeeze(prec_crm(ilon,ilat,:))>0.3)')); % raining points 0.3 as drizzle
%                nrain = numel(irain); % total # of raining points
%                if nrain > minrc
%                  [irains, nrains] = continuousSet(irain,nrain,crm_nx); % raining points in the i-th HCS
%                  if max(nrains)>=minrc
                ixrc = intersect(ixhcs{k}, ixpf); % candidate points of RC before subsetting (precipitatino feature in an HCS)
                nrc = numel(ixrc); % total points
                if nrc<minrc; continue; end
                [ixrcs, nrcs] = continuousSet(ixrc,nrc,1,crm_nx); % RC subsets in the i-th HCS
                [ixrcs, nrcs] = mergeboundary(ixrcs,nrcs,crm_nx); % RC subsets in the i-th HCS
                nrc_max = max(nrcs); % the number of points with maximum RC in i-th HCS
                if nrc_max<minrc; continue; end
                inrc_max = find(nrcs==nrc_max);
                for j = 1:numel(inrc_max) % if more than one maximum RC with the same number of points, go through all of them
                  ixrc_max = ixrcs{inrc_max(j)}; % every HCS has at least one max RC zonal index set ixrc_max
                  if ( nrc_max/nrc>=0.7 & any(tfcold(ixrc_max)) & sum(prec_crm(ilon,ilat,ixrc_max)>6.0)/nrc_max>=0.1 ) % MCS conditions
                    minT = [];
                    for iirc = 1:numel(ixrc_max)
                      minT(iirc) = T_crm(ilon,ilat,ixrc_max(iirc),ictop(ixrc_max(iirc))); % minimum cloud top temp of the MCS
                    end
                    mcs.Tctopmin(illct) = min(minT);
                    mcs.lonlat(illct,1:2)  = [lon(ilon),lat(ilat)]; % remove all the overlapped lon-lat pairs
                    mcs.ilonlat(illct,1:2) = [ilon,ilat];
                    mcs.llx{illct} = ixhcs{k}; % the MCS zonal indices to use
                    mcs.ixrc_max{illct} = ixrc_max; % save the largest RC zonal indices in the MCS
                    mcs.ixrcs{illct} = ixrcs; % save all RC zonal index set in the MCS
                    illct = illct+1;
                    leave = 1;
                    break % once counted, leave this LS grid point and go to the next
                  end
                end    
%              end
%            end
%          end
          if leave; break; end
        end
        if leave; break; end
      end
    end %ilat
  end %ilon
  allvars = whos;
  memused = sum([allvars.bytes])/16/1e6; disp(['mem = ' num2str(memused)])
  display(['total detected MCS = ' num2str(illct-1)])
  display('end MCS search')
  display('saving mcs')
  save(fout,'mcs')
  display('done saving mcs')
  close all force
toc

