function [lonsquall,latsquall] = plot_mcs_xy_ave(vname,iic,ixshift,dosave)
% plot the grid points crossing the entire storm
%xshift=[-5:5];for i=1:numel(xshift);iic=19;dosave=0;[lonsq{i} latsq{i}] = plot_mcs_xy_ave({'spdtz'},iic,xshift(i),dosave);end
% ixshift: shifting the central storm point by a natural number
% iic=287; squall-line like cell orientation, not sure how to dissect the flow
% iic=19;vname={'LI','ps','du'};dosave=0;plot_mcs_xy_ave(vname,iic,dosave)
% iic=19;vname={'spdtz','camdtz','spdqz','camdqz','spd12z','camd12z'};dosave=0;plot_mcs_xy_ave(vname,iic,dosave)
% iic=19;vname={'nompdtz','mpdtz'}
% iic=18;vname={'LI','ps','du'};dosave=0;plot_mcs_xy_ave(vname,iic,dosave)
%dosave=1;for iv=[6];vn={'ps','camdt','nompdtz','mpdtz','du','LI','spdtz','camdtz','spdqz','camdqz','spd12z','camd12z'}; iic=14; ; for  itr=[1,3,5,7,21,23,38,40]; plot_mcs_xy_ave(vn{iv},iic,itr,dosave);pause(0.5);end;close all;end
%  dosave=1;iv=1;vn={'zmdt'}; iic=14; ; for  itr=[1,5,40]; plot_mcs_xy_ave(vn{iv},iic,itr,dosave);pause(0.5);end
% itr: range of time index average 
% 12,21: land
% 24: end over land
% 28: similar to backbuilding (BB) from Schumacher and Johnson 2005 <-- midlat not lat
% 19: similar to training line-Adjoining stratiform (TL/AS)
% 30: well organized event with shear in front
% 16: with stratiform downdraft
% 13 (it=12): stratiform downdraft rear inflow jet
% 14 (it=6): presquall low, stratiform rain-induced mesohigh (cooling), wake low (heating)
%    (it=14): coma!
% obs: 1) parallel shear (35,36,39,40,43) : 
%         - multiple cold pool (symmetric instability?)
%         - cold pool closer to shear
%      2) propagating away shear (36) : cold pool max overlaps precc max
%      3) propagating towards shear (9,10,13,14,16,37,38,41,47), cold pool max leads precc max
%      4) little shear (48)
% US lon lat:
% lon = 235 to 300 (-125 to -60) 
% lat = 25 to 50
% plot the x-z contour of variables for a cluster
doextrat=0;
if ~exist('ixshift') 
 ixshift=0;
end
%%%%%%%%%% PICK A LIFETIME INDEX SET TO PLOT %%%%%%%%%%%%%%%%%%%%%%%%
%itvec = [1,7,8,9,10,21,45,48]; % set lifetime index here, must be even number. e.g, growth(1,2),mature(3,4), decay(5,6)
%itvec = [7,21,40]; % set lifetime index here, must be 6. growth(1,2),mature(3,4), decay(5,6)
if any(ismember(vname,{'nompdtz','mpdtz'}))
  itvec = [1, 15, 40];
else
 % itvec = [1:2:40];
 % itvec = [1,3,4,5,6,7,10:16, 38, 40]; % [ 1 4 10 11 12 14 15 16 % 14 15 16 standard 3 branches
%  itvec = [ 1 3 4 5 6 10 11 12 16 38 40] % 14 15 16 standard 3 branches
%  itvec = [1 2 3 4 5 15 38 40]
%  itvec = [3]; % rotate pvecmcs by 20deg
  itvec = [1:7]; %rotate pvecmcs by 20deg
%  itvec = [1:3]; %[1:10]
% itvec=[16]
end
domcsave=0; % storm wind averaged over mcs cluster cells
if dosave
  docbar=0; %change docbar to 1 here 
else
  docbar=0;
end
if ~exist('dosave')
  dosave = 0;
end
if all(ismember(vname,{'spdtz','camdtz','spdqz','camdqz','spd12z','camd12z',...
                   'zmdtz','nompdtz','mpdtz','nompdqz','mpdqz','ucrmz','wcrmz'})) % vert crossections
  docross=1;
  docrosstype=2; % docrosstype=1 simple cross without multi-lat affixed to a
                 % lon, docrosstype=2 do multi-lat average
  icross =0;% crossections: (0) fix lat (1) top-left (2) bottom-left
else
  docross=0;
end
disp(['iic=' num2str(iic)])
load ~/comp
season = 'JJA';
%dlat = 'lat2050';
%dlat = 'lat2525';
dlat = 'lat6060';
xpts = 30;
dpts = [num2str(xpts) 'pts'];
dims = '2d';
casei = 'F_2000_SPCAM_m2005_3hrly2';
diro = ['/Users/g/archive/matlab/' casei '/figure/' dlat '/' season '/mcsmap/tave/'] % output
if strcmp(comp,'MAC')
  %cd /Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly1/atm/hist/lat2050/JJA/mcs_cluster_var/100pts
  diri = ['/Users/g/archive/matlab/' casei '/atm/hist/' dlat '/' season '/mcs_cluster_var/' dpts '/' dims]
  cd(diri)
  load ../../../mcs_cluster_parm.mat
  %cd /Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly1/atm/hist/lat2050/JJA/mcs_cluster_var/30pts
%  load ../../mcs_cluster_u_storm.mat
%  load ../../mcs_cluster_v_storm.mat
%  load ../../mcs_cluster_u_cam.mat
%  load ../../mcs_cluster_v_cam.mat
else
  diri = ['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/' casei '/atm/hist/' dlat '/' season '/mcs_cluster_var/' dpts '/' dims]
  cd(diri)
  parm = mcs_cluster_parm(['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/' casei '/atm/hist/' casei '.cam.h1.0002-06-01-00000.nc']);
%  error(['plot on MAC, not on' comp])
end
if doextrat
  eval(sprintf(['load mcs_cluster_xy_var_' sprintf('%02d',iic) '_extrat.mat']));
else
  eval(sprintf(['load mcs_cluster_xy_var_' sprintf('%02d',iic) '.mat']));
end
load ../../../mcs_clusters_1.mat 
%if ~exist('itr','var')
%  itr = floor(Nt(iic)/10):Nt(iic)-floor(Nt(iic)/10);
%end
scl=10;
usclz=0.01*scl; % 20deg/2000km=0.01, multiply by 10
wsclz=1*scl; 
%usclz = 1;
%wsclz=1;
uscl=0.1;
uscl2=2; % same as ncquiverref
wscl=100;
ip850=7;
ip200=18;
ip400=13;
ic=idxc(iic);
disp(['cluster # = ' num2str(ic)])
npanel = 1;
nprow = 1;
yl = 1:xpts;
xl = 1:xpts;
ext = 14;
%ixc = xpts/2-5:xpts/2+5;
%iyc = xpts/2-5:xpts/2+5;
icent = xpts/2;
ixc = icent-ext:icent+ext;
iyc = icent-ext:icent+ext;
ixcs = icent-3:icent+3;
iycs = icent-3:icent+3;
ixp  = icent-8:icent+8;
iyp  = icent-8:icent+8;
for itt = 1:numel(itvec)
itr = itvec(itt);
if itr>Nt(iic)
  error(['reset time index, ' num2str(itr) ' longer than' num2str(Nt(iic))])
end
if numel(itr)==1 % if itr is just one timestep, then center lon lat around this
  loncs = parm.lon(vo.ilons{ic}{itr}(ixcs));
  latcs = parm.lat(vo.ilat{ic}{itr}(iycs));
  lonc = parm.lon(vo.ilons{ic}{itr}(ixp));
  latc = parm.lat(vo.ilat{ic}{itr}(iyp));
else % if itr more than one timestep, then center at the first timestep
  loncs = parm.lon(vo.ilons{ic}{1}(ixcs));
  latcs = parm.lat(vo.ilat{ic}{1}(iycs));
  lonc = parm.lon(vo.ilons{ic}{1}(ixp));
  latc = parm.lat(vo.ilat{ic}{1}(iyp));
end
xp = 1:numel(ixp);
yp = 1:numel(iyp);
disp(['date = ' num2str(t{mcsillt4Cl{ic}{itr}(1,3)})])
disp(['halflife it' num2str(floor(Nt(iic)/2))])
disp(['max cluster size = ' num2str(max(mcsnll(ic,:)))])
nt = numel(itr); 
%  Nt = Nt(iic);
disp(['ntime=' num2str(itr)]);
disp(['nhour=' num2str(3*itr)]);
%for it=1:1:Nt
%  lonc = parm.lon(vo.ilons{ic}{it}(ixc));
%  latc = parm.lat(vo.ilat{ic}{it}(iyc));
%disp(['it= ' num2str(it)])
%disp(['time='  t{mcsillt4Cl{ic}{it}(1,3)}]); % time=20 the wake low appears with high entropy drawn from subsidence; two kinds of subsidence, rain (low entropy) & dry air (high entropy) 
disp(['starting lon lat=' num2str(parm.lon(mcsilltcentroids{ic}(1,1))) ', ' num2str(parm.lat(mcsilltcentroids{ic}(1,2)))]); % time=20 the wake low appears with high entropy drawn from subsidence; two kinds of subsidence, rain (low entropy) & dry air (high entropy) 
disp(['ending lon lat=' num2str(parm.lon(mcsilltcentroids{ic}(end,1))) ', ' num2str(parm.lat(mcsilltcentroids{ic}(end,2)))]); % time=20 the wake low appears with high entropy drawn from subsidence; two kinds of subsidence, rain (low entropy) & dry air (high entropy) 

figure('units','normalized','outerposition',[0 0 1 1]); 
%%%%%%%%%%%%%%%%%%%%%% CORR (spatial) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp(['du vs spdt3km ' num2str(corr(vo.dU03km_cam{ic}{it}(:),reshape(vo.spdt{ic}{it}(:,:,parm.i3km),[],1),'rows','complete')) ]) % spdt at 3km is similar to precc, so corr is close
%disp(['du vs precc ' num2str(corr(reshape(vo.dU03km_cam{ic}{it}(ixp,iyp),[],1),reshape(vo.precc{ic}{it}(ixp,iyp),[],1),'rows','complete')) ]) % mcs domain corr not obvious
%disp(['du vs CPI ' num2str(corr(vo.dU03km_cam{ic}{it}(:),vo.Bmax{ic}{it}(:),'rows','complete')) ])
%disp(['du vs pc ' num2str(corr(vo.dU03km_cam{ic}{it}(:),vo.precc{ic}{it}(:),'rows','complete')) ]); % 0.2 to 0.5 corr on the large-scale domain (3000km * 3000km)
%disp(['CPI vs pc ' num2str(corr(vo.Bmax{ic}{it}(:),vo.precc{ic}{it+1}(:),'rows','complete'))])
%disp(['pc vs CPI ' num2str(corr(vo.precc{ic}{it}(:),vo.Bmax{ic}{it+1}(:),'rows','complete'))])
%disp(['CPI vs pl ' num2str(corr(vo.Bmax{ic}{it}(:),vo.precl{ic}{it}(:),'rows','complete'))])
%disp(['pc vs LI ' num2str(corr(vo.precc{ic}{it}(:),vo.LIave{ic}{it}(:),'rows','complete'))]) % no corr
%%%%%%%%%%%%%%%%% 
  w = 0; u = 0; v = 0; thetaa=0; du=[]; pc=0; pl=0; us=zeros(xpts,xpts); 
  vs=zeros(xpts,xpts); LI=0; B=0; cth=0; pva4=0; umcs=0; vmcs=0; ws=0;
  ps=0; usf=0;vsf=0; theta=0; spdt=0; thetas=0; pva1=0; pva16=0; Ts=0; camdt=0;
  spdq=0; camdq=0; spdti=0; camdti=0;spdqi=0; camdqi=0; pt=0; w_precc1=0; w_precl1=0; us1=zeros(xpts,xpts); vs1=zeros(xpts,xpts);
  w_precc4=0; w_precl4=0;up=[];wz=[];thetaez=[]; thetaes=0;
  ucrmz=zeros(numel(ixp),parm.nx,parm.nzi,nt);wcrmz=ucrmz; ulow=0; vlow=0;
  spdtz=[];camdtz=[];spdqz=[];camdqz=[];spd12z=[];camd12z=[];pvaz=[];
  thetaaz=[];thetaz=[];qiz=[];qcz=[];qrz=[];mpdtz=[];nompdtz=[];

  for it = 1:nt %Nt(iic) % time average is taken in this loop automatically
    % obtain ilon ilat (according to vo.u_cam indices using 30 pts)
    if domcsave
      [v1,ip1] = ismember(mcsillt4Cl{ic}{itr(it)}(:,1),vo.ilons{ic}{itr(it)});
      [v2,ip2] = ismember(mcsillt4Cl{ic}{itr(it)}(:,2),vo.ilat{ic}{itr(it)});
      ip3 = mcsillt4Cl{ic}{itr(it)}(:,2)>=(max(mcsillt4Cl{ic}{itr(it)}(:,2))-10);
%      ip3 = mcsillt4Cl{ic}{itr(it)}(:,2)>=(max(mcsillt4Cl{ic}{itr(it)}(:,2))-1);
      utmp(it) = mean(mean(mean(vo.u_cam{ic}{itr(it)}(ip1(ip3),ip2(ip3),1:parm.i6km))));
      vtmp(it) = mean(mean(mean(vo.v_cam{ic}{itr(it)}(ip1(ip3),ip2(ip3),1:parm.i6km))));
    else
     % utmp(it) = mean(mean(mean(vo.u_cam{ic}{itr(it)}(ixcs,iycs,1:parm.i6km))));% mcs speed
     % vtmp(it) = mean(mean(mean(vo.v_cam{ic}{itr(it)}(ixcs,iycs,1:parm.i6km))));
%  propag vector/negative inflow jet from Corfidi paper
%  http://www.spc.noaa.gov/publications/corfidi/mccthes.htm
%        utt = vo.u_cam{ic}{itr(it)}(ixcs,iycs,1:3); 
%        vtt = vo.v_cam{ic}{itr(it)}(ixcs,iycs,1:3);
%        [~,iut] = max(abs(utt(:)));
%        [~,ivt] = max(abs(vtt(:)));
%        utmp(it) = -utt(iut); 
%        vtmp(it) = -vtt(ivt); 
        utmp(it) = mean(mean(mean(vo.u_cam{ic}{itr(it)}(ixcs,iycs,1:parm.i6km))));% mcs speed
        vtmp(it) = mean(mean(mean(vo.v_cam{ic}{itr(it)}(ixcs,iycs,1:parm.i6km))));
       % utmp(it) = -mean(mean(mean(vo.u_cam{ic}{itr(it)}(ixcs,iycs,1:2))));% mcs speed
       % vtmp(it) = -mean(mean(mean(vo.v_cam{ic}{itr(it)}(ixcs,iycs,1:2))));
%    utmp = vo.uclust{ic}{itr(it)};
%    vtmp = vo.vclust{ic}{itr(it)};
    end
    utmp1(it) = mean(mean(mean(vo.u_cam{ic}{itr(it)}(ixp,iyp,1:parm.i6km))));% synoptic low propagation speed
    vtmp1(it) = mean(mean(mean(vo.v_cam{ic}{itr(it)}(ixp,iyp,1:parm.i6km))));
    ulow = utmp1(it);
    vlow = vtmp1(it);
    umcs = utmp(it);
    vmcs = vtmp(it);
    for iiit=1:numel(vo.u{ic}) % calc lifecycle deep shear vector, for ic=19, shear is reduced over time, which is good for TC growth
      %dudeep(iiit) = mean(mean(vo.u{ic}{iiit}(ixp,iyp,ip200) - vo.u{ic}{iiit}(ixp,iyp,ip850),1),2);
      %dvdeep(iiit) = mean(mean(vo.v{ic}{iiit}(ixp,iyp,ip200) - vo.v{ic}{iiit}(ixp,iyp,ip850),1),2);
      %uall(:,iiit)= mean(mean(vo.u{ic}{iiit}(ixp,iyp,:),1),2);
      %vall(:,iiit)= mean(mean(vo.v{ic}{iiit}(ixp,iyp,:),1),2);
      dudeep(iiit) = mean(mean(vo.u{ic}{iiit}(:,:,ip200) - vo.u{ic}{iiit}(:,:,ip850),1),2);
      dvdeep(iiit) = mean(mean(vo.v{ic}{iiit}(:,:,ip200) - vo.v{ic}{iiit}(:,:,ip850),1),2);
      [ix200,iy200] = find(vo.pva{ic}{iiit}(:,:,7)==max(max(vo.pva{ic}{iiit}(:,:,7))));
      [ix850,iy850] = find(vo.pva{ic}{iiit}(:,:,1)==max(max(vo.pva{ic}{iiit}(:,:,1))));
      dx(iiit) = parm.lon(vo.ilons{ic}{iiit}(ix200))-parm.lon(vo.ilons{ic}{iiit}(ix850));
      dy(iiit) = parm.lat(vo.ilat{ic}{iiit}(iy200))-parm.lat(vo.ilat{ic}{iiit}(iy850));
      %quiver(vo.u{ic}{iiit}(:,:,ip850)',vo.v{ic}{iiit}(:,:,ip850)');pause
      pvecshear(:,iiit)=[dudeep(iiit),dvdeep(iiit)]/norm([dudeep(iiit),dvdeep(iiit)]); % proj storm unit vector, neg for inflow
      pvectilt(:,iiit) = [dx(iiit),dy(iiit)]/norm([dx(iiit),dy(iiit)]);
      Degrees(iiit) = atan2d(norm(cross([dx(iiit),dy(iiit),1],[pvecshear(:,iiit)',1])),dot([dx(iiit),dy(iiit)],pvecshear(:,iiit)));
    end
    %quiver(pvectilt(1,:),pvectilt(2,:)); hold on
    %quiver(pvecshear(1,:),pvecshear(2,:)); pause
    %plot(Degrees); hold on
    %contour(abs(uall),30)
    %plot(dudeep,'b'); hold on
    %plot(dvdeep,'y'); 
    %plot(sqrt(dudeep.^2+dvdeep.^2),'r');pause % time series of TC deep shear
    pvecmcs=[dudeep(itr(it)),dvdeep(itr(it))]/norm([dudeep(itr(it)),dvdeep(itr(it))]); % proj storm unit vector, neg for inflow
    pvecmcs = rotate_coord(20,pvecmcs)'; % rotate the pvecmcs  cyclonic (+) or anticyclonic (-)
    %pvecmcs =  [-1,1]/norm([-1 1]);
    %pvecmcs=[umcs,vmcs]/norm([umcs,vmcs]); % proj storm unit vector, neg for inflow
    pveclow=[ulow,vlow]/norm([ulow,vlow]);
    if docross
      %for i = 1:numel(ixp)
%   ppvec=[0.5;-0.5]; %projection unit vector, diag of [1 -1]
      %pvec=[0.7985;-0.6020]; %projection unit vector, dlon ~ 139130m, dlat ~ 104910m
      if docrosstype == 1 
        nj = numel(ixp);
      elseif docrosstype == 2
        [ilonlat,xsquall,xlonsquall,ylatsquall,lonsquall,latsquall]=findsquall(lonc,latc,pvecmcs',ixshift);pause(0.5);%[umcs,vmcs]); 
        nj = size(ilonlat,1);
%        ilontmp = unique(ilonlat(:,1)); % meridionally average
%        nj = numel(ilontmp); % unique # of lon, meridionally average the variables  
      end
%      for i=1:nj  
      for i=1:nj
        if docrosstype==1
          switch icross
          case 0
            iloncs(i) = ixcs(i);
            ilatcs(i) = xpts/2; % fix lat crossection
          case 1
            iloncs(i) = ixcs(i);
            ilatcs(i) = iycs(end-i+1); % diag crossection from top left
          case 2
            iloncs(i) = ixcs(i);
            ilatcs(i) = iycs(i); % diag crossection from bottom left
          end
          ix = iloncs(i);
          iy = ilatcs(i);
        elseif docrosstype==2 % cut through mcs propag dir
          ix = ixp(ilonlat(i,1));
          iy = iyp(ilonlat(i,2));
%          if ilonlat(i+1,1)==ilonlat(i,1)
%          ix = ixp(ilontmp(i));
%          iy = iyp(ilonlat(ilonlat(:,1)==ilontmp(i),2)); % multiple lat indices at a fixed lon
          iloncs(i) = ix;
          ilatcs(i) = iy;
%          ilatcs(i) = iy(1);
        end
        % mean is taken over meridional index for the multi-lat cross section
 %       ucrmz(i,:,:,it) = mean(vo.u_crmi{ic}{itr(it)}(ix,iy,:,:),2); 
 %       wcrmz(i,:,:,it) = mean(vo.w_crmi{ic}{itr(it)}(ix,iy,:,:),2);
        thetaez(i,:,it) = mean(vo.thetaei{ic}{itr(it)}(ix,iy,:),2);
        camdtz(i,:,it)  = mean(vo.camdti{ic}{itr(it)}(ix,iy,:),2);
        spdtz(i,:,it)   = mean(vo.spdti{ic}{itr(it)}(ix,iy,:),2);
        camdqz(i,:,it)  = mean(vo.camdqi{ic}{itr(it)}(ix,iy,:)*(-2260*86400),2);
        spdqz(i,:,it)   = mean(vo.spdqi{ic}{itr(it)}(ix,iy,:)*(-2260*86400),2);
        mpdtz(i,:,it)   = mean(vo.mpdt{ic}{itr(it)}(ix,iy,:),2);
        nompdtz(i,:,it) = mean(vo.zmdt{ic}{itr(it)}(ix,iy,:)+...
                          vo.evaptzm{ic}{itr(it)}(ix,iy,:)+...
                          vo.cmfdt{ic}{itr(it)}(ix,iy,:)+...
                          vo.macpdt{ic}{itr(it)}(ix,iy,:),2);
        mpdqz(i,:,it)   = mean(vo.mpdq{ic}{itr(it)}(ix,iy,:)+...
                          vo.macpdq{ic}{itr(it)}(ix,iy,:),2);
        nompdqz(i,:,it) = mean(vo.zmdq{ic}{itr(it)}(ix,iy,:)+...
                        vo.evapqzm{ic}{itr(it)}(ix,iy,:)+...
                        vo.cmfdq{ic}{itr(it)}(ix,iy,:),2);%+...
                    %    vo.macpdq{ic}{itr(it)}(ix,iy,:),2);
%         zmdtz(i,:,it)= mean(vo.zmdq{ic}{itr(it)}(ix,iy,:)*(-2260)+...
%                        vo.evapqzm{ic}{itr(it)}(ix,iy,:)*(-2260)+...
%                        vo.cmfdq{ic}{itr(it)}(ix,iy,:)*(-2260)+...
%                        vo.macpdq{ic}{itr(it)}(ix,iy,:)*(-2260),2);%+...
        thetaaz(i,:,it) = mean(vo.thetaa{ic}{itr(it)}(ix,iy,:),2);
        thetaz(i,:,it)  = mean(vo.theta{ic}{itr(it)}(ix,iy,:),2);
        pvaz(i,1:19,it) = mean(vo.pva{ic}{itr(it)}(ix,iy,:),2);
        pvaz(i,20,it)   = 0; % 20th level missing due to vert grad calc
        uz(i,:,it)      = mean(vo.u_cam{ic}{itr(it)}(ix,iy,:),2)-umcs;
        vz(i,:,it)      = mean(vo.v_cam{ic}{itr(it)}(ix,iy,:),2)-vmcs;  
        wz(i,:,it)      = mean(vo.w{ic}{itr(it)}(ix,iy,:),2);
        qiz(i,:,it)     = mean(vo.qii{ic}{itr(it)}(ix,iy,:),2);
        qcz(i,:,it)     = mean(vo.qci{ic}{itr(it)}(ix,iy,:),2);
        qrz(i,:,it)     = mean(vo.qri{ic}{itr(it)}(ix,iy,:),2);
        qvz(i,:,it)     = mean(vo.qvi{ic}{itr(it)}(ix,iy,:),2);
        up(i,:,it)      = -[uz(i,:,it)',vz(i,:,it)']*pvecmcs'; % project the u,v vector onto the unit vector
      end
    else
%    ps = vo.ps{ic}{itr(it)}(ixp,iyp)/nt-mean(mean(vo.ps{ic}{itr(it)}(ixp,iyp)))/nt + ps; 
%    ps = vo.ps{ic}{itr(it)}(ixp,iyp)/nt-parm.p0/nt + ps; 
      ps = vo.ps{ic}{itr(it)}(ixp,iyp)/nt + ps; 
%    spdt = vo.spdt{ic}{itr(it)}(ixp,iyp,4)/nt + spdt; 
%    camdt = vo.camdt{ic}{itr(it)}(ixp,iyp,4)/nt + camdt; 
%    spdq = vo.spdq{ic}{itr(it)}(ixp,iyp,4)/nt + spdq; 
%    camdq = vo.camdq{ic}{itr(it)}(ixp,iyp,4)/nt + camdq; 
      spdti = vo.spdti{ic}{itr(it)}(ixp,iyp,4)/nt + spdti; 
      camdti = vo.camdti{ic}{itr(it)}(ixp,iyp,4)/nt + camdti; 
      spdqi = vo.spdqi{ic}{itr(it)}(ixp,iyp,4)*(-2260*86400)/nt + spdqi; % 9 (7km 400hPa) 10 (8.5km 350hPa) is stronger than 6 (5km 550hPa local dip), 4 (3km 700hPa)
      camdqi = vo.camdqi{ic}{itr(it)}(ixp,iyp,4)*(-2260*86400)/nt + camdqi; 
      vo.precc{ic}{itr(it)}(isnan(vo.precc{ic}{itr(it)}(:)))=0;
      vo.precl{ic}{itr(it)}(isnan(vo.precl{ic}{itr(it)}(:)))=0;
      pc = vo.precc{ic}{itr(it)}(ixp,iyp)/nt + pc;
      pl = vo.precl{ic}{itr(it)}(ixp,iyp)/nt + pl;
      pt = vo.prect{ic}{itr(it)}(ixp,iyp)*3.6e6/nt + pt;
      cth = vo.cthick{ic}{itr(it)}(ixp,iyp)/nt + cth;
      B = vo.Bave{ic}{itr(it)}(ixp,iyp)/nt + B;
      pva16 = vo.pva{ic}{itr(it)}(ixp,iyp,16)/nt + pva16; % 250 hPa at 12km
      pva4 = vo.pva{ic}{itr(it)}(ixp,iyp,parm.i3km)/nt + pva4; % at 3km
      pva1 = vo.pva{ic}{itr(it)}(ixp,iyp,1)/nt + pva1;
      LI = vo.LIave{ic}{itr(it)}(ixp,iyp)/nt + LI;
%      ws = vo.w{ic}{itr(it)}(ixp,iyp,1)/nt + ws; % surface w
      ws = vo.ws{ic}{itr(it)}(ixp,iyp)/nt + ws; % surface w
      w = vo.w{ic}{itr(it)}(ixp,iyp,4)/nt + w; % midlev w
      for ilo=1:numel(ixp)
        for ila=1:numel(iyp)
          w_theta(ilo,ila,:,it) = interp1(squeeze(vo.theta{ic}{itr(it)}(ixp(ilo),iyp(ila),:)), squeeze(vo.w{ic}{itr(it)}(ixp(ilo),iyp(ila),:)), [299:1:350]','linear'); 
        end; 
      end
%      vo.w_precc{ic}{itr(it)}(isnan(vo.w_precc{ic}{itr(it)}(:)))=0;
%      vo.w_precl{ic}{itr(it)}(isnan(vo.w_precl{ic}{itr(it)}(:)))=0;
%      w_precc1 = vo.w_precc{ic}{itr(it)}(ixp,iyp,2)/nt + w_precc1; % midlev w
%      w_precl1 = vo.w_precl{ic}{itr(it)}(ixp,iyp,2)/nt + w_precl1; % midlev w
%      w_precc4 = vo.w_precc{ic}{itr(it)}(ixp,iyp,4)/nt + w_precc4; % midlev w
%      w_precl4 = vo.w_precl{ic}{itr(it)}(ixp,iyp,4)/nt + w_precl4; % midlev w
%      thetaa = vo.thetaa{ic}{itr(it)}(ixp,iyp,1)/nt + thetaa;
%      theta = vo.theta{ic}{itr(it)}(ixp,iyp,1)/nt + theta;
%      thetas = vo.thetas{ic}{itr(it)}(ixp,iyp)/nt + thetas;
      thetaes = vo.thetae{ic}{itr(it)}(ixp,iyp,1)/nt + thetaes;
%      Ts = vo.Ts{ic}{itr(it)}(ixp,iyp)/nt + Ts;
      u = vo.u_cam{ic}{itr(it)}(ixp,iyp,parm.i3km)/nt-utmp(it)/nt + u; %storm relative flow
      v = vo.v_cam{ic}{itr(it)}(ixp,iyp,parm.i3km)/nt-vtmp(it)/nt + v;
      usf = vo.u_cam{ic}{itr(it)}(ixp,iyp,1)/nt + usf; %storm relative flow
      vsf = vo.v_cam{ic}{itr(it)}(ixp,iyp,1)/nt + vsf;
%      usf = vo.u_cam{ic}{itr(it)}(ixp,iyp,1)/nt-utmp(it)/nt + usf; %storm relative flow
%      vsf = vo.v_cam{ic}{itr(it)}(ixp,iyp,1)/nt-vtmp(it)/nt + vsf;
%      du = vo.dU03km_cam{ic}{itr(it)}(ixp,iyp)/nt + du;
      du(:,:,it) = (diff(vo.u_cam{ic}{itr(it)}(ixp,iyp,[1,parm.i3km]),1,3));
%      dv(:,:,it) = (diff(vo.v_cam{ic}{itr(it)}(ixp,iyp,[1,parm.i3km]),1,3));
%      dU(:,:,it) = sqrt(du(:,:,it)^2+dv(:,:,it)^2);
    end
  end
  if docross
    thetaeaz  = thetaez - mean(mean(thetaez,1),3);
  end
%  ucrmz = reshape(ucrmz,size(ilonlat,1)*parm.nx,parm.nzi,nt);
%  wcrmz = reshape(wcrmz,size(ilonlat,1)*parm.nx,parm.nzi,nt);
  if docross
    if numel(itr)==1 % if itr is just one timestep, then center lon lat around this
      lonc = parm.lon(vo.ilons{ic}{itr}(iloncs));
      latc = parm.lat(vo.ilat{ic}{itr}(ilatcs));
    else % if itr more than one timestep, then center at the first timestep
      lonc = parm.lon(vo.ilons{ic}{1}(iloncs));
      latc = parm.lat(vo.ilat{ic}{1}(ilatcs));
    end 
  end
  du = mean(du,3);
%  dU = mean(dU,3);
  us(xpts/2,xpts/2) = pvecmcs(1);%umcs;%+us(xpts/2,xpts/2);
  vs(xpts/2,xpts/2) = pvecmcs(2);%vmcs;%+vs(xpts/2,xpts/2);
  us1(xpts/2,xpts/2) = pveclow(1);%umcs;%+us(xpts/2,xpts/2);
  vs1(xpts/2,xpts/2) = pveclow(2);%vmcs;%+vs(xpts/2,xpts/2);
  disp(['mcs speed = ' num2str(sqrt(umcs^2+vmcs^2))])
  disp(['mcs [u,v] = ' num2str(umcs) ' ' num2str(vmcs)])
  disp(['low speed = ' num2str(sqrt(ulow^2+vlow^2))])
  disp(['low [u,v] = ' num2str(ulow) ' ' num2str(vlow)])
%  us(xpts/2,xpts/2) = mean(utmp);%+us(xpts/2,xpts/2);
%  vs(xpts/2,xpts/2) = mean(vtmp);%+vs(xpts/2,xpts/2);
%  us(xpts/2,xpts/2) = median([vo.uclust{ic}{itr(1:end)}]);%+us(xpts/2,xpts/2);
%  vs(xpts/2,xpts/2) = median([vo.vclust{ic}{itr(1:end)}]);%+vs(xpts/2,xpts/2);
%  contourf(ws',250,'linestyle','none'); colorbar; 
%  caxis([-0.01 0.01]); colormap(hh(ip),cmap(24)); hold on; 
%  maxw=round(max(ws(:)),2)-0.02;
cs=3; % linewidth of all the contours
%%%%%%%%% plot isosurface of pva in 3D %%%%%%%%%%%%%%%%%%%%%
lonc1 = parm.lon(vo.ilons{ic}{itr}(ixp));
latc1 = parm.lat(vo.ilat{ic}{itr}(iyp));
%{
figure;
p = patch(isosurface(lonc1-360,latc1,parm.zint/1000,permute(vo.thetaei{ic}{itr}(ixp,iyp,:)-mean(mean(mean(vo.thetaei{ic}{itr}(ixp,iyp,:)))),[2,1,3]),1)); hold on
%isonormals(lonc1-360,latc1,1:28,permute(vo.thetae{ic}{itr(it)}(ixp,iyp,:)-350,[2,1,3]),p)
set(p,'FaceColor','red');
%  contourf(lonc1-360,latc1,ps'/100,[998:.5:1013],'linestyle','none'); hold on
%  contour(lonc1-360,latc1,ps'/100,[998:.5:1013],'color',[0 0.5 1]); hold on
%p1= patch(isosurface(lonc1-360,latc1,parm.zint,permute(vo.thetaei{ic}{itr(it)}(ixp,iyp,:)-350,[2,1,3]),-8)); hold on
%set(p1,'FaceColor','yellow');
p2=patch(isosurface(lonc1-360,latc1,parm.zint/1000,permute(vo.qci{ic}{itr}(ixp,iyp,:),[2,1,3]),0.05e-3)); hold on
set(p2,'FaceColor','blue');
%isosurface(lonc1-360,latc1,parm.zint,permute(vo.qci{ic}{itr(it)}(ixp,iyp,:),[2,1,3]),0.03e-3); hold on
p3=patch(isosurface(lonc1-360,latc1,parm.zint/1000,permute(vo.qii{ic}{itr}(ixp,iyp,:),[2,1,3]),0.07e-3))
set(p3,'FaceColor','blue');
stormarrow_plan(lonc1-360,latc1,us,vs,uscl2,ixp,iyp,[0 0 0])
crossline(-pvecmcs,lonc1-360,latc1);
pause
%}
%{
figure;
p4 = patch(isosurface(lonc1-360,latc1,parm.zint(1:end-1)/1000,permute(vo.pva{ic}{itr}(ixp,iyp,:),[2,1,3]),0.5e-6)); hold on
set(p4,'FaceColor','red')
stormarrow_plan(lonc1-360,latc1,us,vs,uscl2,ixp,iyp,[0 0 0])
crossline(-pvecmcs,lonc1-360,latc1);
pause
%}
%  contourf(lonc1-360,latc1,thetaes'-350,50); pause
%  contour(lonc1-360,latc1,vo.thetae{ic}{itr(it)}(ixp,iyp,2)'-350,50); pause
%%%%%%%% plot 2D slices of horiz velocity field with height %%%%%%%%%%%%%%%%%
%{
for i=1:numel(parm.zint)
  figure;
%  [x1,y1]=meshgrid(parm.lon(vo.ilons{ic}{itr}(ixp)),parm.lat(vo.ilat{ic}{itr}(iyp)));
%  quiver(x1,y1,vo.u_cam{ic}{itr(it)}(ixp,iyp,i)',vo.v_cam{ic}{itr(it)}(ixp,iyp,i)'); hold on
%  contour(lonc1-360,latc1,vo.pv{ic}{itr(it)}(ixp,iyp,i)',[.5e-6 .5e-6]); 
%  contourf(lonc1-360,latc1,vo.w{ic}{itr(it)}(ixp,iyp,i)'); caxis([-0.02 0.2])
  %dp = vo.pm{ic}{itr(it)}(ixp,iyp,i)-vo.ps{ic}{itr(it)}(ixp,iyp); plot pressure on hybridsurf?
  %phi = -dp./vo.rho{ic}{itr(it)}(ixp,iyp,i)/10;
  z1d = vo.z1d{ic}{itr}(ixp,iyp,i);
  contourf(lonc1-360,latc1,z1d',80);colormap(cmap(1))
 % contourf(lonc1-360,latc1,vo.pm{ic}{itr(it)}(ixp,iyp,i)',80);colormap(cmap(1))
%  caxis([-15 15])
%  crossline(-pvecmcs,lonc1-360,latc1);
end
pause
%}
%%%%%%% plot 3D velocity field %%%%%%%%%%%%%%%%%%%%
%{
[x1,y1,z1]=meshgrid(parm.lon(vo.ilons{ic}{itr}(ixp))-360,parm.lat(vo.ilat{ic}{itr}(iyp)),parm.zint/1000);
ipos = vo.w{ic}{itr(it)}>0;
ineg = vo.w{ic}{itr(it)}<0;
wpos = vo.w{ic}{itr(it)};
wneg = vo.w{ic}{itr(it)};
wpos(ineg)=NaN;
wneg(ipos)=NaN;
quiver3(x1,y1,z1,permute(vo.u_cam{ic}{itr(it)}(ixp,iyp,:),[2,1,3])*0.03,permute(vo.v_cam{ic}{itr(it)}(ixp,iyp,:),[2,1,3])*0.05,permute(wpos(ixp,iyp,:),[2,1,3])*30,'autoscale','off','color','g'); % plot only positive w quivers
%quiver3(x1,y1,z1,permute(vo.u_cam{ic}{itr(it)}(ixp,iyp,:),[2,1,3])*0.03,permute(vo.v_cam{ic}{itr(it)}(ixp,iyp,:),[2,1,3])*0.05,permute(wneg(ixp,iyp,:),[2,1,3])*30,'autoscale','off','color',[0 0 1]); % plot only negative w quivers
pause
%}
%%%%%%% plot 3D streamline %%%%%%%%%%%%%%%%%%%%
%{
figure
%[x1,y1,z1]=meshgrid(parm.lon(vo.ilons{ic}{itr}(ixp)),parm.lat(vo.ilat{ic}{itr}(iyp)),parm.zint);
[x1,y1,z1]=meshgrid([1:17],[1:17],[1:20]);
streamline(x1,y1,z1,permute(vo.u_cam{ic}{itr(it)}(ixp,iyp,:),[2,1,3])*0.1,permute(vo.v_cam{ic}{itr(it)}(ixp,iyp,:),[2,1,3])*0.1,permute(vo.w{ic}{itr(it)}(ixp,iyp,:),[2,1,3])*10,x1(:,:,1),y1(:,:,1),1*ones(1,numel(z1(:,:,1)))); pause
%}
for ivv = 1:numel(vname)
%pause
clf
ip=1;
hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
switch vname{ivv}
%%%%%%%%%%%%%%%% SHEAR %%%%%%%%%%%%%%%%
case 'du'
  contourf(lonc,latc,du',[-50:2:50],'linestyle','none'); 
  caxis([-20 20]); cm = colormap(hh(ip),cmap(24,0,20)); hold on; 
  contour(lonc,latc,du',[-50:2:0],'b:'); 
  contour(lonc,latc,du',[0:2:50],'-','color',[0.81 0.53 0.53]); 
  cl = [0,0,0];
%  contour(lonc,latc,w_theta(:,:,29)',[-0.5:0.002:0],':','color',cl,'linewidth',cs); 
%  contour(lonc,latc,w_theta(:,:,29)',[0:0.002:0.5],'-','color',cl,'linewidth',cs); 
  contour(lonc,latc,ws',[-0.02:0.002:0],':','color',cl,'linewidth',cs); 
  contour(lonc,latc,ws',[0:0.002:0.02],'-','color',cl,'linewidth',cs); 
  %contour(lonc,latc,w_precl1',[-1:0.1:0],':','color',cl,'linewidth',cs); 
  %contour(lonc,latc,w_precl1',[0:0.1:1],'-','color',cl,'linewidth',cs); 
%  contour(lonc,latc,ws',[0 0],'-.','color',cl,'linewidth',cs); 
  contour(lonc,latc,pc',[1:3:50],'r','linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
  contour(lonc,latc,B',[40:10:200],'c','linewidth',cs);
  sublabel_plan(itt,{'c','d'},'\DeltaU')
  set_mcs_domain(parm,mcsillt4Cl{ic}{itr(it)},lonc,latc,[1 1 1])
  stormarrow_plan(lonc,latc,us,vs,uscl2,ixp,iyp,[0 0 0])
  %stormarrow_plan(lonc,latc,us1,vs1,uscl2,ixp,iyp,[0.5 0.5 0.5])
  crossline(-pvecmcs,lonc,latc);
  cbar(docbar,cm,[-20 20],[-20:5:20],vname{ivv},'[m s^{-1}]',iic)
  ip=ip+1;
%%%%%%%%%%%%%% uv %%%%%%%%%%%%%%%%% 
case 'uv'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip); hold on
  ih = 1;
  %pcolor(lonc,latc,vo.prect{ic}{itr(it)}(ixp,iyp)'*3.6e6); hold on
  contour(lonc,latc,vo.prect{ic}{itr(it)}(ixp,iyp)'*3.6e6,[1:1:30],'r','linewidth',cs); hold on
%  contour(lonc,latc,pl',[1:3:30],'y','linewidth',cs); hold on
%  contour(lonc,latc,pc',[1:3:50],'r','linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
  %contour(lonc,latc,B',[40:10:200],'c','linewidth',cs);
%  contourf(lonc,latc,vo.z1d{ic}{itr}(ixp,iyp,ih)',20);colormap(cmap(1))
%  cm = colormap(hh(ip),cmap(1,0,16)); hold on; 
%  quiver(lonc,latc,vo.u_cam{ic}{itr(it)}(ixp,iyp,ih)'*0.1,vo.v_cam{ic}{itr(it)}(ixp,iyp,ih)'*0.1,'autoscale','off')
  %sublabel_plan(itt,{'e','f'},'(u,v)')
  set_mcs_domain(parm,mcsillt4Cl{ic}{itr(it)},lonc,latc,[0 0 0])
  stormarrow_plan(lonc,latc,us,vs,uscl2,ixp,iyp,[0 0 0])
  %stormarrow_plan(lonc,latc,us1,vs1,uscl2,ixp,iyp,[0.5 0.5 0.5])
  crossline(-pvecmcs,lonc,latc);
  ip=ip+1;
%%%%%%%%%%%%%% LI %%%%%%%%%%%%%%%%% 
case 'LI'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  contourf(lonc,latc,LI',linspace(-10,0,16),'linestyle','none');  caxis([-8 0]); 
  cm = colormap(hh(ip),cmap(15,1,16)); hold on; 
  contourf(lonc,latc,LI',linspace(-10,0,16),'linestyle',':','color',[1 0.2 0]);  
  contour(lonc,latc,w',[0 0],'k-.','linewidth',cs); 
  contour(lonc,latc,w',[0:0.02:0.2],'k','linewidth',cs); 
  contour(lonc,latc,w',[-0.2:0.02:0],'k:','linewidth',cs); 
  cl4 = [.25 .75 .9];
  cl1 = [0 221 0]/255;
  contour(lonc,latc,pva4',[-5e-6:5e-7:0],':','linecolor',cl4,'linewidth',cs); % 1PVU = 10-6 m-2 s-1 K kg-1
  contour(lonc,latc,pva4',[0:5e-7:5e-6],'-','linecolor',cl4,'linewidth',cs); 
%  contour(lonc,latc,pva16',[5e-7:5e-7:5e-6],'-','linecolor',cl1,'linewidth',cs); 
%  contour(lonc,latc,pva16',[0 0],'-.','linecolor',cl1,'linewidth',cs); 
%  contour(lonc,latc,pva16',[-5e-6:5e-7:-5e-7],':','linecolor',cl1,'linewidth',cs); 
  ncquiverref(lonc,latc,usf',vsf','m/s',10,'true',[0 0 0]);
%  quiver(lonc,latc,vo.u_cam{ic}{itr(it)}(ixp,iyp,ih)'*0.1,vo.v_cam{ic}{itr(it)}(ixp,iyp,ih)'*0.1,'autoscale','off')
  sublabel_plan(itt,{'e','f'},'LI')
  set_mcs_domain(parm,mcsillt4Cl{ic}{itr(it)},lonc,latc,[1 1 1])
  stormarrow_plan(lonc,latc,us,vs,uscl2,ixp,iyp,[0 0 0])
  %stormarrow_plan(lonc,latc,us1,vs1,uscl2,ixp,iyp,[0.5 0.5 0.5])
  crossline(-pvecmcs,lonc,latc);
  cbar(docbar,cm,[-8 0],[-8:0],vname{ivv},'[K]',iic)
  ip=ip+1;

%%%%%%%%%%%%%%%%% PS %%%%%%%%%%%%%%%%%%%% 
case 'ps'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  psm = mean(ps(:))/100;
%  ps1=1004; ps2=1013;
  ps1=psm-200; ps2=psm+200;
  contourf(lonc,latc,ps'/100,[psm-200:psm+200],'linestyle','none'); hold on
  contour(lonc,latc,ps'/100,[psm-200:psm+200],'color',[0 0.5 1]); hold on
%  contourf(lonc,latc,ps'/100,[998:.5:1013],'linestyle','none'); hold on
%  contour(lonc,latc,ps'/100,[998:.5:1013],'color',[0 0.5 1]); hold on
  % presquall low (front), mesohigh (stratiform rain-induced), wake low (high-entropy subsidence)
  caxis([998 ps2]); cm = colormap(hh(ip),cmap(27,1,20)); hold on; 
  clev = [-2:3:15];
  tmp = cmap(17,0,numel(clev));
  thetaes(ceil(size(thetaes,1)/2),ceil(size(thetaes,2)/2))
  for ii=1:numel(clev)
    contour(lonc,latc,thetaes'-350,[clev(ii) clev(ii)],'linecolor',tmp(ii,:),'linewidth',cs);
%    contour(lonc,latc,thetas',[clev(ii) clev(ii)],'linecolor',tmp(ii,:),'linewidth',cs);
%    contour(lonc,latc,Ts',[clev(ii) clev(ii)],'linecolor',tmp(ii,:),'linewidth',cs);
  end
  contour(lonc,latc,pl',[1:3:30],'y','linewidth',cs); hold on
  set_mcs_domain(parm,mcsillt4Cl{ic}{itr(it)},lonc,latc,[1 1 1])
  stormarrow_plan(lonc,latc,us,vs,uscl2,ixp,iyp,[0 0 0])
  %stormarrow_plan(lonc,latc,us1,vs1,uscl2,ixp,iyp,[0.5 0.5 0.5])
  crossline(-pvecmcs,lonc,latc);
  sublabel_plan(itt,{'a','b'},'p_s')
  cbar(docbar,cm,[ps1 ps2],[ps1:3:ps2],vname{ivv},'[hPa]',iic)
  ip=ip+1;
%%%%%%%%%%%%%%%%% SPDT %%%%%%%%%%%%%%%%%%%% 
case 'spdt'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  contourf(lonc,latc,spdti',[-100:5:100],'linestyle','none'); hold on
  contour(lonc,latc,spdti',[-100:5:0],'b:'); hold on
  contour(lonc,latc,spdti',[0:5:100],'linestyle','-','color',[1,.5,0]); hold on
  contour(lonc,latc,pt',[1:2:16],'m','linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
  caxis([-100 100]); cm = colormap(hh(ip),cmap(10,0,20)); hold on; 
  set_mcs_domain(parm,mcsillt4Cl{ic}{itr(it)},lonc,latc,[0 0 0])
  stormarrow_plan(lonc,latc,us,vs,uscl2,ixp,iyp,[0 0 0])
  crossline(-pvecmcs,lonc,latc);
  cbar(docbar,cm,[-100 100],[-100:20:100],vname{ivv},'[K/day]',iic)
  ip=ip+1;
%%%%%%%%%%%%%%%%% CAMDT %%%%%%%%%%%%%%%%%%%% 
case 'camdt'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  contourf(lonc,latc,camdti',[-100:5:100],'linestyle','none'); hold on
  contour(lonc,latc,camdti',[-100:5:0],'b:'); hold on
  contour(lonc,latc,camdti',[0:5:100],'linestyle','-','color',[1,.5,0]); hold on
  caxis([-100 100]); cm = colormap(hh(ip),cmap(10,0,20)); hold on; 
  set_mcs_domain(parm,mcsillt4Cl{ic}{itr(it)},lonc,latc,[0 0 0])
  stormarrow_plan(lonc,latc,us,vs,uscl2,ixp,iyp,[0 0 0])
  crossline(-pvecmcs,lonc,latc);
  cbar(docbar,cm,[-100 100],[-100:20:100],vname{ivv},'[K/day]',iic)
  ip=ip+1;
%%%%%%%%%%%%%%%%% SPDQ %%%%%%%%%%%%%%%%%%%% 
case 'spdq'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  contourf(lonc,latc,spdqi',[-100:5:100],'linestyle','none'); hold on
  caxis([-100 100]); cm = colormap(hh(ip),cmap(10,0,20)); hold on; 
  contour(lonc,latc,spdqi',[-100:5:0],'b:'); hold on
  contour(lonc,latc,spdqi',[0:5:100],'linestyle','-','color',[1,.5,0]); hold on
  contour(lonc,latc,pt',[1:2:16],'m','linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
  set_mcs_domain(parm,mcsillt4Cl{ic}{itr(it)},lonc,latc,[0 0 0])
  stormarrow_plan(lonc,latc,us,vs,uscl2,ixp,iyp,[0 0 0])
  crossline(-pvecmcs,lonc,latc);
  cbar(docbar,cm,[-100 100],[-100:20:100],vname{ivv},'[K/day]',iic)
  ip=ip+1;
%%%%%%%%%%%%%%%%% CAMDQ %%%%%%%%%%%%%%%%%%%% 
case 'camdq'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  contourf(lonc,latc,camdqi',[-100:5:100],'linestyle','none'); hold on
  caxis([-100 100]); cm = colormap(hh(ip),cmap(10,0,20)); hold on; 
  contour(lonc,latc,camdqi',[-100:5:0],'b:'); hold on
  contour(lonc,latc,camdqi',[0:5:100],'linestyle','-','color',[1,.5,0]); hold on
%  contour(lonc,latc,camdqi',[0 0],'k-.'); hold on
  set_mcs_domain(parm,mcsillt4Cl{ic}{itr(it)},lonc,latc,[0 0 0])
  stormarrow_plan(lonc,latc,us,vs,uscl2,ixp,iyp,[0 0 0])
  crossline(-pvecmcs,lonc,latc);
  cbar(docbar,cm,[-100 100],[-100:20:100],vname{ivv},'[K/day]',iic)
  ip=ip+1;
case 'spdtz'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  [x y] = domeshz(xsquall,parm);
  clev = [-100:5:100];
  contourf(x,y,mean(spdtz,3)',clev,'linestyle','none'); hold on
  contour(x,y, mean(spdtz,3)',clev(clev>=0),'-','color',[.8 .5 .5]); hold on
  contour(x,y, mean(spdtz,3)',clev(clev<0),':','color',[.3 .6 .6]); hold on
  caxis([clev(1) clev(end)]); cm = colormap(hh(ip),cmap(10,0,numel(clev)-1)); hold on; 
%%{
  contour(x,y,mean(pvaz,3)',[.5e-6:1e-6:5e-6],'color',[0,0,.4],'linestyle','-','linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
  contour(x,y,mean(pvaz,3)',[-5e-6:1e-6:-.5e-6],'linestyle',':','color',[0 0 .4],'linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
%  contour(x,y,mean(thetaaz,3)','color',[.9,.0,.0],'linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
%  contour(x,y,mean(thetaaz,3)',[0.5:1:3],'color','m','linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
%  contour(x,y,mean(thetaaz,3)',[-3:1:-0.5],'m--','linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
  contour(x,y,mean(thetaeaz,3)',[1:3:20],'m','linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
  contour(x,y,mean(thetaeaz,3)',[-20:3:-1],'m--','linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
  %contour(x,y,mean(qvz,3)'*1e3,[0.01:0.05:10],'color',[.7,.7,.7],'linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
  quiver_vert(x,y,up,usclz,wz,wsclz); hold on % up is the projected velocity onto the zonal vertical plane
%%}
%  streamline(x,y,mean(up,3)'*usclz,mean(wz,3)'*wsclz,x,y)
%  plot_env_wind(x,y,up,usclz)
  annot([.15,.8,.1,.1],'(a) Q1c');
  stormarrow_vertcross()
  set_mcs_domain_vert(parm,xsquall)
  cbar(docbar,cm,[-100 100],[-100:20:100],vname{ivv},'[K/day]',iic)
  ip=ip+1;
case 'camdtz'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  [x y] = domeshz(xsquall,parm);
  clev = [-100:5:100];
  contourf(x,y,mean(camdtz,3)',clev,'linestyle','none'); hold on
  contour(x,y,mean(camdtz,3)',clev(clev>=0),'-','color',[.8 .5 .5]); hold on
  contour(x,y,mean(camdtz,3)',clev(clev<0),':','color',[.3 .6 .6]); hold on
  caxis([clev(1) clev(end)]); cm = colormap(hh(ip),cmap(10,0,numel(clev)-1)); hold on; 
%  quiver_vert(x,y,up,usclz,wz,wsclz); hold on
  annot([.15,.8,.1,.1],'(b) Q1c');
  stormarrow_vertcross()
  set_mcs_domain_vert(parm,xsquall)
  cbar(docbar,cm,[-100 100],[-100:20:100],vname{ivv},'[K/day]',iic)
  ip=ip+1;
case 'spdqz'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  [x y] = domeshz(xsquall,parm);
  clev = [-100:5:100];
  contourf(x,y,mean(spdqz,3)',clev,'linestyle','none'); hold on
  contour(x,y,mean(spdqz,3)',clev(clev>=0),'-','color',[.8 .5 .5]); hold on
  contour(x,y,mean(spdqz,3)',clev(clev<0),':','color',[.3 .6 .6]); hold on
%{
  contour(x,y,mean(qiz,3)'*1e3,[0.01:0.1:10],'color',[.7,.7,.7],'linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
  contour(x,y,mean(qcz,3)'*1e3,[0.01:0.1:10],'color',[.5,.5,.5],'linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
  contour(x,y,mean(qrz,3)'*1e3,[0.01:0.1:10],'k','linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
%}
%  contour(x,y,mean(spdtz-spdqz,3)',clev(clev>=0),'-','color',[1,0,1],'linewidth',cs); hold on
%  contour(x,y,mean(spdtz-spdqz,3)',clev(clev<0),':','color',[1,0,1],'linewidth',cs); hold on
  caxis([clev(1) clev(end)]); cm = colormap(hh(ip),cmap(10,0,numel(clev)-1)); hold on; 
%  contour(x,y,mean(qiz,3)'*1e3,[0.01:0.05:10],'color',[.7,.7,.7],'linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
%  contour(x,y,mean(qcz,3)'*1e3,[0.01:0.05:10],'color',[.5,.5,.5],'linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
%  contour(x,y,mean(qrz,3)'*1e3,[0.01:0.05:10],'k','linewidth',cs);  % 1[mm/hr] = 0.1[cm/hr] = 2.4[cm/day] => 5[mm/hr] = 12[cm/day]
%  quiver_vert(x,y,up,usclz,wz,wsclz); hold on
  annot([.15,.8,.1,.1],'(c) Q2')
%  stormarrow_vertcross()
  set_mcs_domain_vert(parm,xsquall)
  ip=ip+1;
case 'camdqz'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  [x y] = domeshz(xsquall,parm);
  clev = [-100:5:100];
  contourf(x,y,mean(camdqz,3)',clev,'linestyle','none'); hold on
  contour(x,y,mean(camdqz,3)',clev(clev>=0),'-','color',[.8 .5 .5]); hold on
  contour(x,y,mean(camdqz,3)',clev(clev<0),':','color',[.3 .6 .6]); hold on
  caxis([clev(1) clev(end)]); cm = colormap(hh(ip),cmap(10,0,numel(clev)-1)); hold on; 
  quiver_vert(x,y,up,usclz,wz,wsclz); hold on
  annot([.15,.8,.1,.1],'(d) Q2');
  stormarrow_vertcross()
  set_mcs_domain_vert(parm,xsquall)
  ip=ip+1;
case 'spd12z'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  [x y] = domeshz(xsquall,parm);
  clev = [-100:5:100];
  contourf(x,y,mean(spdtz-spdqz,3)',clev,'linestyle','none'); hold on
  contour(x,y,mean(spdtz-spdqz,3)',clev(clev>=0),'-','color',[.8 .5 .5]); hold on
  contour(x,y,mean(spdtz-spdqz,3)',clev(clev<0),':','color',[.3 .6 .6]); hold on
  caxis([clev(1) clev(end)]); cm = colormap(hh(ip),cmap(10,0,numel(clev)-1)); hold on; 
  quiver_vert(x,y,up,usclz,wz,wsclz); hold on
  annot([.18,.8,.1,.1],'(e) Q1c-Q2');
  stormarrow_vertcross()
  set_mcs_domain_vert(parm,xsquall)
  ip=ip+1;
case 'camd12z'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  [x y] = domeshz(xsquall,parm);
  clev = [-100:5:100];
  contourf(x,y,mean(camdtz-camdqz,3)',clev,'linestyle','none'); hold on
  contour(x,y,mean(camdtz-camdqz,3)',clev(clev>=0),'-','color',[.8 .5 .5]); hold on
  contour(x,y,mean(camdtz-camdqz,3)',clev(clev<0),':','color',[.3 .6 .6]); hold on
  caxis([clev(1) clev(end)]); cm = colormap(hh(ip),cmap(10,0,numel(clev)-1)); hold on; 
  %quiver_vert(x,y,up,usclz,wz,wsclz)
  annot([.18,.8,.1,.1],'(f) Q1c-Q2');
  stormarrow_vertcross()
  set_mcs_domain_vert(parm,xsquall)
  ip=ip+1;
case 'nompdtz'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  [x y] = domeshz(xsquall,parm);
  clev = [-100:5:100];
  contourf(x,y,mean(nompdtz*86400,3)',clev,'linestyle','none'); hold on
  contour(x,y,mean(nompdtz*86400,3)',clev(clev>=0),'-','color',[.8 .5 .5]); hold on
  contour(x,y,mean(nompdtz*86400,3)',clev(clev<0),':','color',[.3 .6 .6]); hold on
  caxis([clev(1) clev(end)]); cm = colormap(hh(ip),cmap(10,0,numel(clev)-1)); hold on; 
  %quiver_vert(x,y,up,usclz,wz,wsclz)
  sublabel_mpdt(itt,itr,{'b','d','f'})
%  stormarrow_vertcross()
  set_mcs_domain_vert(parm,xsquall)
  cbar(docbar,cm,[-100 100],[-100:20:100],vname{ivv},'[K/day]',iic)
  ip=ip+1;
case 'mpdtz'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  [x y] = domeshz(xsquall,parm);
  clev = [-100:5:100];
  contourf(x,y,mean(mpdtz*86400,3)',clev,'linestyle','none'); hold on
  contour(x,y,mean(mpdtz*86400,3)',clev(clev>=0),'-','color',[.8 .5 .5]); hold on
  contour(x,y,mean(mpdtz*86400,3)',clev(clev<0),':','color',[.3 .6 .6]); hold on
  caxis([clev(1) clev(end)]); cm = colormap(hh(ip),cmap(10,0,numel(clev)-1)); hold on; 
  %quiver_vert(x,y,up,usclz,wz,wsclz)
  sublabel_mpdt(itt,itr,{'a','c','e'})
  stormarrow_vertcross()
  set_mcs_domain_vert(parm,xsquall)
  cbar(docbar,cm,[-100 100],[-100:20:100],vname{ivv},'[K/day]',iic)
  ip=ip+1;

case 'nompdqz'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  [x y] = domeshz(xsquall,parm);
  clev = [-100:5:100];
  contourf(x,y,mean(nompdqz*(-2260*86400),3)',clev,'linestyle','none'); hold on
  contour(x,y,mean(nompdqz*(-2260*86400),3)',clev(clev>=0),'-','color',[.8 .5 .5]); hold on
  contour(x,y,mean(nompdqz*(-2260*86400),3)',clev(clev<0),':','color',[.3 .6 .6]); hold on
  caxis([clev(1) clev(end)]); cm = colormap(hh(ip),cmap(10,0,numel(clev)-1)); hold on; 
  %quiver_vert(x,y,up,usclz,wz,wsclz)
  sublabel_mpdt(itt,itr,{'b','d','f'})
  stormarrow_vertcross()
  set_mcs_domain_vert(parm,xsquall)
  cbar(docbar,cm,[-100 100],[-100:20:100],vname{ivv},'[K/day]',iic)
  ip=ip+1;
case 'mpdqz'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  [x y] = domeshz(xsquall,parm);
  clev = [-100:5:100];
  contourf(x,y,mean(mpdqz*(-2260*86400),3)',clev,'linestyle','none'); hold on
  contour(x,y,mean(mpdqz*(-2260*86400),3)',clev(clev>=0),'-','color',[.8 .5 .5]); hold on
  contour(x,y,mean(mpdqz*(-2260*86400),3)',clev(clev<0),':','color',[.3 .6 .6]); hold on
  caxis([clev(1) clev(end)]); cm = colormap(hh(ip),cmap(10,0,numel(clev)-1)); hold on; 
  %quiver_vert(x,y,up,usclz,wz,wsclz)
  sublabel_mpdt(itt,itr,{'a','c','e'})
  stormarrow_vertcross()
  set_mcs_domain_vert(parm,xsquall)
  cbar(docbar,cm,[-100 100],[-100:20:100],vname{ivv},'[K/day]',iic)
  ip=ip+1;
case 'ucrmz'
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  [x y] = domeshz(xsquall,parm);
usclz=0.4;wsclz=2;
  contourf(squeeze(mean(ucrmz,2))');colorbar
%  contourf(squeeze(mean(wcrmz,3))');colorbar
%  quiver_vert(x,y,up,usclz,wz,wsclz)
  if itr==1; str='(a) 0h 6.13N';pos=[.15,.8,.2,.1]; 
  elseif itr==5;  str='(c) 12h 8.01N'; pos=[.15,.8,.2,.1];
  elseif itr==23;  str='(c) 63h 14.61N'; pos=[.15,.8,.2,.1];
  elseif itr==40; str='(e) 117h 19.32N';pos=[.20,.8,.2,.1];end
  annot(pos,str);
%  stormarrow()
%  set_mcs_domain_vert(parm,xsquall)
%  cbar(docbar,cm,[-100 100],[-100:20:100],vname,'[K/day]',iic)
  ip=ip+1;

%  contourf(x,y,camdtz',[-100:10:100],'linestyle','none'); hold on
%  streamline(x,y,uu'*uscl,ww',ones(size(y))*30,y)
%  streamline(x,y,squeeze(vo.u_cam{ic}{itr(it)}(ixp,ilatcs,:))'*uscl,squeeze(vo.w{ic}{itr(it)}(ixp,ilatcs,:))');% fix lat crossec
%{
%%%%%%%%%%%%%%%%% 
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  contourf(thetaa',250,'linestyle','none'); colorbar; caxis([-3 3]); colormap(hh(ip),cmap(3)); hold on; 
  axis square
  hold off;
  ip=ip+1;
%%%%%%%%%%%%%%%%% 
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  contourf(du',250,'linestyle','none'); colorbar; caxis([-15 15]); colormap(hh(ip),cmap(3)); hold on; 
  axis square
  hold off;
  ip=ip+1;
%%%%%%%%%%%%%%%%% 
  hh(ip) = subplot(nprow,ceil(npanel/nprow),ip);
  contourf(B',250,'linestyle','none'); colorbar; caxis([0 100]); colormap(hh(ip),flipud(cmap(2))); hold on; 
  axis square
  hold off;
  ip=ip+1;
%%%%%%%%%%%%%%%% COMPUTE CLUSTER MEAN FLOW %%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
    ncc=30
    for iic=1:ncc 
      eval(sprintf(['load mcs_cluster_xy_var_' sprintf('%02d',iic) '.mat']))
      ic=idxc(iic);
      for i = 1:numel(ixp)
        switch cross
        case 0
          iloncs(i) = ixp(i);
          ilatcs(i) = xpts/2; % fix lat crossection
        case 1
          iloncs(i) = ixp(i);
          ilatcs(i) = iyp(end-i+1); % diag crossection from top left
        case 2
          iloncs(i) = ixp(i);
          ilatcs(i) = iyp(i); % diag crossection from bottom left
        end
        iit = ceil(nt4Cl(ic)/2);
%        umcs = mean(mean(mean(vo.u_cam{ic}{iit}(ixcs,iycs,1:parm.i6km))));
%        vmcs = mean(mean(mean(vo.v_cam{ic}{iit}(ixcs,iycs,1:parm.i6km))));
        umcs = -max(max(max(vo.u_cam{ic}{iit}(ixcs,iycs,1:3))));
        vmcs = -max(max(max(vo.v_cam{ic}{iit}(ixcs,iycs,1:3))));
        pvec=-[umcs/(umcs^2+vmcs^2);vmcs/(umcs^2+vmcs^2)]; % proj storm unit vector
        uz = squeeze(vo.u_cam{ic}{iit}(iloncs(i),ilatcs(i),:))-umcs;
        vz = squeeze(vo.v_cam{ic}{iit}(iloncs(i),ilatcs(i),:))-vmcs;  
        up(i,:) = [uz,vz]*pvec/ncc + up(i,:)'; % project the u,v vector onto the unit vector
        wz(i,:) = squeeze(vo.w{ic}{iit}(iloncs(i),ilatcs(i),:))/ncc + wz(i,:)';
      end
    end
  quiver(up',wz','color',[.3 .6 .6],'linewidth',3,'Autoscale','off'); 
'done'
  pause
%}
%}
%%%%%%%%%%

otherwise
end
savef(hh(ip-1),dosave,diro,iic,itr,vname{ivv})
end
end

%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sublabel_mpdt(iit,itr,labstr)
  str = ['(' labstr{iit} ') ' num2str(3*(itr-1)) 'h'];
  if itr<=34 % less than or equal to double digit
    pos=[.10,.8,.2,.1];
  else
    pos=[.11,.8,.2,.1];
  end
  annot(pos,str);

function sublabel_plan(itt,labstr,labstr2)
  if mod(itt,2) % odd # (a) (c) (e) 
    str = ['(' labstr{1} ') ' labstr2];
  else % even # (b) (d) (f)
    str = ['(' labstr{2} ') ' labstr2];
  end
  pos=[.25,.8,.2,.1];
  annot(pos,str);
function sublabel_vert(iit,itr,labstr)
  if mod(itt,2) % odd # (a) (c) (e) 
    str = ['(' labstr{1} ') ' num2str(3*(itr-1)) 'h'];
  else % even # (b) (d) (f)
    str = ['(' labstr{2} ') ' num2str(3*(itr-1)) 'h'];
  end
  if iit<=34 % less than or equal to double digit
    pos=[.15,.8,.2,.1];
  else
    pos=[.20,.8,.2,.1];
  end
  annot(pos,str);
function annot(pos,text)
  annotation('textbox',pos,'String',text,'FitBoxToText','on','Margin',0,'fontsize',44,'fontweight','bold','linestyle','none','backgroundColor','w','HorizontalAlignment','center')

function quiver_vert(x,y,up,usclz,wz,wsclz) 
  quiver(x,y,mean(up,3)'*usclz,mean(wz,3)'*wsclz,'color',[.3 .6 .6],'linewidth',3,'Autoscale','off'); hold on

function [x y]=domeshz(xsquall,parm)
  [x y] = meshgrid(xsquall/100,parm.zint/1000);

function savef(h,dosave,diro,iic,itr,vname)
if dosave
  if numel(itr)==1
    saveas(h,[diro '/mcs_c' sprintf('%02d',iic) '_t' sprintf('%02d',itr) '_' vname '.png'])
    crop([diro '/mcs_c' sprintf('%02d',iic) '_t' sprintf('%02d',itr) '_' vname '.png'])
  else
    saveas(h,[diro '/mcs_c' sprintf('%02d',iic) '_tave' vname  '.png'])
    crop([diro '/mcs_c' sprintf('%02d',iic) '_tave' vname  '.png'])
  end
  disp('finished saving')
end

function set_mcs_domain_vert(parm,xsq)
  ylim([-.5 parm.zint(end)/1000]);
  xtcks = round([xsq(1):200:xsq(end)]/100); % label every 200km, divide by 100 to make quiver arrow show up 
  set(gca,'XTick',xtcks);
  for i=1:numel(xtcks)
    if mod(i,2)
      xticks{i} = num2str((xtcks(i)-xtcks(ceil(numel(xtcks)/2)))*100);
    else 
      xticks{i} = ''; 
    end
  end
  set(gca,'XTicklabel',xticks);
  xlim([xsq(1),xsq(end)]/100);
%{
  xlim([lonc(1) lonc(end)])
  xtcks = round(lonc(1)):4:round(lonc(end));
  set(gca,'xtick',xtcks(1:numel(xtcks)-1));
  if xtcks(1)>180
    xtcks = 360-xtcks;
    for i = 1:numel(xtcks)
      xtc{i} = [num2str(xtcks(i)) 'W'];
    end
  else
    for i = 1:numel(xtcks)
      xtc{i} = [num2str(xtcks(i)) 'E'];
    end
  end  
  set(gca,'XTicklabel',xtc)
  disp(['lat centered at ' num2str(latc(1))])
%}
  %for i=1:numel(lonc)
  %  str{i} = ['(' sprintf('%0.2f',round(latc(i),2)) 'N,' sprintf('%0.2f',360-lonc(i)) 'W)'];
  %end
  %set(gca,'XTick',lonc(4:5:numel(lonc)-2))
  %set(gca,'XTicklabel',str(4:5:numel(lonc)-2))
  ylabel('[km]');
  set(gca,'YTick',1:3:18); %parm.zint/1000)
%  set(gca,'YTicklabel',2:2:18)
%  set(gca,'TickDir','out')
  set(gca,'fontsize',44);

function set_mcs_domain(parm,mcsillt4Cl,lonc,latc,colr)
  scatter(parm.lon(mcsillt4Cl(:,1)),parm.lat(mcsillt4Cl(:,2)),2700,colr,'s','linewidth',2); 
  coast_centered(0)
  xlim([lonc(1) lonc(end)]);
  ylim([latc(1) latc(end)]);
  %xticks((round(lonc(1)):4:round(lonc(end))))
  xtcks = round(lonc(1)):4:round(lonc(end));
  set(gca,'xtick',xtcks(1:numel(xtcks)-1));
  if xtcks(1)>180
    xtcks = 360-xtcks;
    for i = 1:numel(xtcks)
      xtc{i} = [num2str(xtcks(i)) 'W'];
    end
  else
    for i = 1:numel(xtcks)
      xtc{i} = [num2str(xtcks(i)) 'E'];
    end
  end 
  set(gca,'XTicklabel',xtc)
%  xTick=get(gca,'xtick'); 
%  newXTick=linspace(min(xTick),max(xTick),numel(ticks)); 
%  set(gca,'xtick',newXTick(1:numel(ticks)-1));
  ytcks = round(latc(1)):4:round(latc(end));
  for i=1:numel(ytcks)
    ytc{i} = [num2str(ytcks(i)) 'N'];
  end
  set(gca,'YTicklabel',ytc) 
  set(gca,'ytick',ytcks(1:numel(ytcks)));
  set(gca,'fontsize',44)
  axis square

function cbar(docbar,cm,climit,clev,vname,units,iic)
if docbar
  figure('units','normalized','outerposition',[0 0 1 1]); 
  axis off; 
  horz = 1; % make a horizontal colorbar, otherwise vertical as normal
  if horz
    hc = colorbar('southoutside'); 
  else
    hc = colorbar('eastoutside'); 
  end
  set(hc,'ytick',clev)
  colormap(cm); 
  caxis(climit);
 % colormap(cmap(ncmap(ivv(vname)),flipm(ivv(vname)),10)); 
 % caxis(caxislim{ivv(vname)}); 
  if horz %horizontal cbar
    hc.Position = [0.05 0.2170 0.7750 0.0181]; % 1&2 coord, 3&4 width&height
    hc.Label.String = units; 
    hc.Label.HorizontalAlignment = 'right';
    hc.Label.Position = [climit(2)+0.2*diff(climit) 0 0];
  else % vertical cbar
    hc.Position = [0.9 0.15 0.0181 0.7750]; % 1&2 coord, 3&4 width&height
  %  ylabel(hc,units) % label on the side
    set(get(hc,'title'),'string',units,'horizontalAlignment', 'left'); % label on the top
  end
  set(gca,'Fontsize',44); 
  saveas(gcf,['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/mcs_c' num2str(iic) '_' vname '.png']); 
  crop(['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/mcs_c' num2str(iic) '_' vname '.png']); 
end

function stormarrow_plan(lonc,latc,us,vs,uscl2,ixp,iyp,color,lw)
if ~exist('lw')
  lw=6;
end
  quiver(lonc,latc,us(ixp,iyp)'*uscl2,vs(ixp,iyp)'*uscl2,'color',color,'linewidth',lw,'autoscale','off');

function stormarrow_vertcross
  xstart = 0.77;
  xend = xstart-0.05;
  ystart = 0.125; 
  yend=ystart;
%  annotation('textarrow',[xstart,xend],[ystart,yend],'String','shear_{200-850}','fontsize',30,'linewidth',5,'headwidth',25,'HeadStyle','plain')


function plot_env_wind(x,y,up,usclz)
  tmp = zeros(size(x)); % vert-zonal matrix
  tmp(:,end) = mean(mean(up,3),1); % time and zonal average, and put at right
  quiver(x,y,tmp*usclz,zeros(size(x)),'k','linewidth',5,'Autoscale','off','clipping','off'); hold on
%  xlim([x(1),x(end)])

function lat = crossline(pvec,lon,lat)
  latc=lat(ceil(numel(lat)/2));lonc=lon(ceil(numel(lat)/2));
  for i=1:numel(lon)
    a = (lon(i)-lonc)/pvec(1);
    lati(i) = latc + pvec(2)*a;
  end
  plot(lon,lati,'k--','linewidth',2); hold on
%  plot(lati,lon,'k*'); hold on

