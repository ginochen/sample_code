% time-corr map of e and dtcond
%load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/e_eabs_p14p30-ave.mat','e_pm','eabs_pm','lat','lon')
load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/e_eabs_p14p21m.mat','eabs_p14p21m','lat','lon');
load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/e_eabs_p14p30m.mat','e_p14p21m','lat','lon');
idv.cmfmcdzm=6; varname{idv.cmfmcdzm} = 'CMFMCDZM'; cmfmcdzm_pi29pi30 = loadvar(varname{idv.cmfmcdzm},camCase,1,1,[1,1,29,1],[Inf,Inf,2,Inf]); % just read level from 850 to 985 skip 1000
idv.cmfmc=7; varname{idv.cmfmc} = 'CMFMC'; cmfmc_pi29pi30 = loadvar(varname{idv.cmfmc},camCase,1,1,[1,1,29,1],[Inf,Inf,2,Inf]); % just read level from 850 to 985 skip 1000

ilat = find(lat>=-35 & lat<=35);
ilon = 1:numel(lon);
icase =30;
vname = {'e_dtcond','eabs_dtcond','e_spdt','eabs_spdt',...
         'e_preccdzm','eabs_preccdzm', 'e_cmfmcdzm_pi24pi30m','eabs_cmfmcdzm_pi24pi30m',...
         'e_cmfmc_pi24pi30m','eabs_cmfmc_pi24pi30m','eabs_p24p30m_cmfmc_pi24pi30m','eabs_p24p30m_cmfmcdzm_pi24pi30m',...
         'eabs_p14p30m_preccsh','eabs_p24p30m_relhum_p24p30m','eabs_p14p24m_relhum_p14p24m','eabs_p14p24m_cmfmc_pi24pi30m',...
         'eabs_p14p24m_cmfmcdzm_pi24pi30m','eabs_p14p24m_cmfmc_pi14pi24m','eabs_p14p24m_cmfmcdzm_pi14pi24m','eabs_p24p30m_preccsh',...
         'eabs_p24p30m_preccdzm','eabs_p14p24m_relhum_p24p30m','eabs_p14p24m_preccdzm','eabs_p14p24m_preccsh',...
         'eabs_p14p30m_preccdzm','eabs_p14p21m_preccdzm','eabs_p14p21m_preccsh','eabs_p21p30m_preccdzm',...
         'eabs_p21p30m_preccsh','eabs_p14p21m_cmfmc_pi14pi21m','eabs_p14p21m_cmfmc_pi21pi30m','eabs_p21p30m_cmfmc_pi14pi21m',...
         'eabs_p21p30m_cmfmc_pi21pi30m','eabs_p14p21m_cmfmcdzm_pi14pi21m','eabs_p14p21m_cmfmcdzm_pi21pi30m','eabs_p21p30m_cmfmcdzm_pi14pi21m','eabs_p21p30m_cmfmcdzm_pi21pi30m',...
         'eabs_p14p21m_cmfmcdzm_pi29','eabs_p14p21m_cmfmc_pi29','eabs_p14p21m_cmfmc_pi30','eabs_p14p30m_cmfmc_pi29','eabs_p14p30m_cmfmcdzm_pi29','e_p14p30m_cmfmcdzm_pi29',...
         'e_p14p30m_cmfmc_pi29','e_p14p30m_pzm'...
}; %33 counts

%
%evalc(['load(''/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/e_' casename{sdcase} '.mat'',''' casename{sdcase} ''',''lat'',''lon'')']);
switch icase
case 45 % Pdeep vs dQ_200-1000
  [rhomap.e_p14p30m_preccdzm,pmap.e_p14p30m_preccdzm] = corrcoef_time_map(e_p14p30m,var{idv.pzm},ilat,ilon,1);
   figure; hold on
   title('Time-Correlation dQ_{200-1000hPa} v.s. P_{Deep}')
case 44 % cloud-base mass flux corr dQ_200-1000
  [rhomap.e_p14p30m_cmfmc_pi29,pmap.e_p14p30m_cmfmc_pi29] = corrcoef_time_map(e_pm,squeeze(cmfmc_pi29pi30(:,:,1,:)),ilat,ilon,1); % 30th layer of ilev = 985
   figure; hold on
   title('Time-Correlation dQ_{200-1000hPa} v.s. MF_{Shallow 967hPa}')
case 43 % cloud-base mass flux corr dQ_200-1000
  [rhomap.e_p14p30m_cmfmcdzm_pi29,pmap.e_p14p30m_cmfmcdzm_pi29] = corrcoef_time_map(e_pm,squeeze(cmfmcdzm_pi29pi30(:,:,1,:)),ilat,ilon,1); % 30th layer of ilev = 985
   figure; hold on
   title('Time-Correlation dQ_{200-1000hPa} v.s. MF_{Deep 967hPa}')
case 42 % cloud-base mass flux corr |dQ|_200-1000
  [rhomap.eabs_p14p30m_cmfmcdzm_pi29,pmap.eabs_p14p30m_cmfmcdzm_pi29] = corrcoef_time_map(eabs_p14p30m,squeeze(cmfmcdzm_pi29pi30(:,:,1,:)),ilat,ilon,1); % 30th layer of ilev = 985
   figure; hold on
   title('Time-Correlation |dQ|_{200-1000hPa} v.s. MF_{Deep 967hPa}')
case 41 % cloud-base mass flux corr |dQ|_200-1000
  [rhomap.eabs_p14p30m_cmfmc_pi29,pmap.eabs_p14p30m_cmfmc_pi29] = corrcoef_time_map(eabs_p14p30m,squeeze(cmfmc_pi29pi30(:,:,1,:)),ilat,ilon,1); % 30th layer of ilev = 985
   figure; hold on
   title('Time-Correlation |dQ|_{200-1000hPa} v.s. MF_{Shallow 967hPa}')
case 40 % cloud-base mass flux corr |dQ|
  [rhomap.eabs_p14p21m_cmfmc_pi30,pmap.eabs_p14p21m_cmfmc_pi30] = corrcoef_time_map(eabs_p14p21m,squeeze(cmfmc_pi29pi30(:,:,2,:)),ilat,ilon,1); % 30th layer of ilev = 985
   figure; hold on
   title('Time-Correlation |dQ|_{200-700hPa} v.s. MF_{Shallow 985hPa}')
case 39 % cloud-base mass flux corr |dQ|
  [rhomap.eabs_p14p21m_cmfmc_pi29,pmap.eabs_p14p21m_cmfmc_pi29] = corrcoef_time_map(eabs_p14p21m,squeeze(cmfmc_pi29pi30(:,:,1,:)),ilat,ilon,1); % 29th layer of ilev = 967
   figure; hold on
   title('Time-Correlation |dQ|_{200-700hPa} v.s. MF_{Shallow 967hPa}')
case 38 % cloud-base mass flux corr |dQ|
  [rhomap.eabs_p14p21m_cmfmcdzm_pi29,pmap.eabs_p14p21m_cmfmcdzm_pi29] = corrcoef_time_map(eabs_p14p21m,squeeze(cmfmcdzm_pi29pi30(:,:,1,:)),ilat,ilon,1); % 29th layer of ilev = 967
   figure; hold on
   title('Time-Correlation |dQ|_{200-700hPa} v.s. MF_{Deep 967hPa}')
case 37
  [rhomap.eabs_p21p30m_cmfmcdzm_pi21pi30m,pmap.eabs_p21p30m_cmfmcdzm_pi21pi30m] = corrcoef_time_map(eabs_p21p30m,cmfmcdzm_pi21pi30m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{700-1000mb} v.s. MF_{Deep 700-1000mb}')
case 36
  [rhomap.eabs_p21p30m_cmfmcdzm_pi14pi21m,pmap.eabs_p21p30m_cmfmcdzm_pi14pi21m] = corrcoef_time_map(eabs_p21p30m,cmfmcdzm_pi14pi21m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{700-1000mb} v.s. MF_{Deep 200-700mb}')
case 35
  [rhomap.eabs_p14p21m_cmfmcdzm_pi21pi30m,pmap.eabs_p14p21m_cmfmcdzm_pi21pi30m] = corrcoef_time_map(eabs_p14p21m,cmfmcdzm_pi21pi30m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{200-700mb} v.s. MF_{Deep 700-1000mb}')
case 34
  [rhomap.eabs_p14p21m_cmfmcdzm_pi14pi21m,pmap.eabs_p14p21m_cmfmcdzm_pi14pi21m] = corrcoef_time_map(eabs_p14p21m,cmfmcdzm_pi14pi21m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{200-700mb} v.s. MF_{Deep 200-700mb}')
case 33
  [rhomap.eabs_p21p30m_cmfmc_pi21pi30m,pmap.eabs_p21p30m_cmfmc_pi21pi30m] = corrcoef_time_map(eabs_p21p30m,cmfmc_pi21pi30m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{700-1000mb} v.s. MF_{Shallow 700-1000mb}')
case 32
  [rhomap.eabs_p21p30m_cmfmc_pi14pi21m,pmap.eabs_p21p30m_cmfmc_pi14pi21m] = corrcoef_time_map(eabs_p21p30m,cmfmc_pi14pi21m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{700-1000mb} v.s. MF_{Shallow 200-700mb}')
case 31
  [rhomap.eabs_p14p21m_cmfmc_pi21pi30m,pmap.eabs_p14p21m_cmfmc_pi21pi30m] = corrcoef_time_map(eabs_p14p21m,cmfmc_pi21pi30m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{200-700mb} v.s. MF_{Shallow 700-1000mb}')
case 30
  [rhomap.eabs_p14p21m_cmfmc_pi14pi21m,pmap.eabs_p14p21m_cmfmc_pi14pi21m] = corrcoef_time_map(eabs_p14p21m,cmfmc_pi14pi21m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{200-700mb} v.s. MF_{Shallow 200-700mb}')
case 29 
  [rhomap.eabs_p21p30m_preccsh,pmap.eabs_p21p30m_preccsh] = corrcoef_time_map(eabs_p21p30m,var{idv.psh},ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |dQ|_{700-1000mb} v.s. P_{Shallow}')
case 28 
  [rhomap.eabs_p21p30m_preccdzm,pmap.eabs_p21p30m_preccdzm] = corrcoef_time_map(eabs_p21p30m,var{idv.pzm},ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{700-1000mb} v.s. P_{Deep}')
case 27 
  [rhomap.eabs_p14p21m_preccsh,pmap.eabs_p14p21m_preccsh] = corrcoef_time_map(eabs_p14p21m,var{idv.psh},ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{200-700mb} v.s. P_{Shallow}')
case 26
   [rhomap.eabs_p14p21m_preccdzm,pmap.eabs_p14p21m_preccdzm] = corrcoef_time_map(eabs_p14p21m,var{idv.pzm},ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{200-700mb} v.s. P_{Deep}')
case 25
   [rhomap.eabs_p14p30m_preccdzm,pmap.eabs_p14p30m_preccdzm] = corrcoef_time_map(eabs_p14p30m,var{idv.pzm},ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |dQ|_{200-1000mb} v.s. P_{Deep}')
case 24
   [rhomap.eabs_p14p24m_preccsh,pmap.eabs_p14p24m_preccsh] = corrcoef_time_map(eabs_p14p24m,var{idv.psh},ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{200-850mb} v.s. P_{Shallow}')
case 23
   [rhomap.eabs_p14p24m_preccdzm,pmap.eabs_p14p24m_preccdzm] = corrcoef_time_map(eabs_p14p24m,var{idv.pzm},ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{200-850mb} v.s. P_{Deep}')
case 22
   [rhomap.eabs_p14p24m_relhum_p24p30m,pmap.eabs_p14p24m_relhum_p24p30m] = corrcoef_time_map(eabs_p14p24m,relhum_p24p30m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of |Q_{SPCAM} - Q_{CAM}|_{200-850mb} v.s. RH_{850-1000mb}')

case 21
   [rhomap.eabs_p24p30m_preccdzm,pmap.eabs_p24p30m_preccdzm] = corrcoef_time_map(eabs_p24p30m,var{idv.pzm},ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{850-1000mb} v.s. P_{Deep}')
case 20
   [rhomap.eabs_p24p30m_preccsh,pmap.eabs_p24p30m_preccsh] = corrcoef_time_map(eabs_p24p30m,var{idv.psh},ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{850-1000mb} v.s. P_{Shallow}')
case 19
   [rhomap.eabs_p14p24m_cmfmcdzm_pi14pi24m,pmap.eabs_p14p24m_cmfmcdzm_pi14pi24m] = corrcoef_time_map(eabs_p14p24m,cmfmcdzm_pi14pi24m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{200-850mb} v.s. MF_{Deep 200-850mb}')
case 18
   [rhomap.eabs_p14p24m_cmfmc_pi14pi24m,pmap.eabs_p14p24m_cmfmc_pi14pi24m] = corrcoef_time_map(eabs_p14p24m,cmfmc_pi14pi24m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{200-850mb} v.s. MF_{Shallow 200-850mb}')
case 16
   [rhomap.eabs_p14p24m_cmfmc_pi24pi30m,pmap.eabs_p14p24m_cmfmc_pi24pi30m] = corrcoef_time_map(eabs_p14p24m,cmfmc_pi24pi30m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{200-850mb} v.s. MF_{Shallow 850-1000mb}')
case 17
   [rhomap.eabs_p14p24m_cmfmcdzm_pi24pi30m,pmap.eabs_p14p24m_cmfmcdzm_pi24pi30m] = corrcoef_time_map(eabs_p14p24m,cmfmcdzm_pi24pi30m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation |Q_{SPCAM} - Q_{CAM}|_{200-850mb} v.s. MF_{Deep 850-1000mb}')
case 13
   [rhomap.eabs_p14p30m_preccsh,pmap.eabs_p14p30m_preccsh] = corrcoef_time_map(eabs_p14p30m,var{idv.psh},ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of |dQ|_{200-1000mb} v.s. P_{Shallow}')
case 14
   [rhomap.eabs_p24p30m_relhum_p24p30m,pmap.eabs_p24p30m_relhum_p24p30m] = corrcoef_time_map(eabs_p24p30m,relhum_p24p30m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of |Q_{SPCAM} - Q_{CAM}|_{850-1000mb} v.s. RH_{850-1000mb}')
case 15
   [rhomap.eabs_p14p24m_relhum_p14p24m,pmap.eabs_p14p24m_relhum_p14p24m] = corrcoef_time_map(eabs_p14p24m,relhum_p14p24m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of |Q_{SPCAM} - Q_{CAM}|_{200-850mb} v.s. RH_{200-850mb}')
case 1
%   [rhomap.e_dtcond,pmap.e_dtcond] = corrcoef_time_map(e_pm,dtcond_pm,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of pressure-weighted-averaged Q_{SPCAM} - Q_{CAM} v.s. Q_{CAM}')
case 2
%   [rhomap.eabs_dtcond,pmap.eabs_dtcond] = corrcoef_time_map(eabs_pm,dtcond_pm,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of pressure-weighted-averaged |Q_{SPCAM} - Q_{CAM}| v.s. Q_{CAM}')
case 3
%   [rhomap.e_spdt,pmap.e_spdt] = corrcoef_time_map(e_pm,spdt_pm,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of pressure-weighted-averaged Q_{SPCAM} - Q_{CAM} v.s. Q_{SPCAM}')
case 4
%   [rhomap.eabs_spdt,pmap.eabs_spdt] = corrcoef_time_map(eabs_pm,spdt_pm,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of pressure-weighted-averaged |Q_{SPCAM} - Q_{CAM}| v.s. Q_{SPCAM}')
case 5
%   [rhomap.e_preccdzm,pmap.e_preccdzm] = corrcoef_time_map(e_pm,var{idv.pzm},ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of pressure-weighted-averaged Q_{SPCAM} - Q_{CAM} v.s. P_{Deep}')
case 6
%   [rhomap.eabs_preccdzm,pmap.eabs_preccdzm] = corrcoef_time_map(eabs_pm,var{idv.pzm},ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of pressure-weighted-averaged |Q_{SPCAM} - Q_{CAM}| v.s. P_{Deep}')
case 7
   [rhomap.e_cmfmcdzm_pi24pi30m, pmap.e_cmfmcdzm_pi24pi30m] = corrcoef_time_map(e_pm,cmfmcdzm_pi24pi30m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of Q_{SPCAM} - Q_{CAM} v.s. MF_{Deep 850-1000mb}')
case 8
   [rhomap.eabs_cmfmcdzm_pi24pi30m, pmap.eabs_cmfmcdzm_pi24pi30m] = corrcoef_time_map(eabs_pm,cmfmcdzm_pi24pi30m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of |Q_{SPCAM} - Q_{CAM}| v.s. MF_{Deep 850-1000mb}')
case 9
   [rhomap.e_cmfmc_pi24pi30m, pmap.e_cmfmc_pi24pi30m] = corrcoef_time_map(e_pm,cmfmc_pi24pi30m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of Q_{SPCAM} - Q_{CAM} v.s. MF_{Shallow 850-1000mb}')
case 10
   [rhomap.eabs_cmfmc_pi24pi30m, pmap.eabs_cmfmc_pi24pi30m] = corrcoef_time_map(eabs_pm,cmfmc_pi24pi30m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of |Q_{SPCAM} - Q_{CAM}| v.s. MF_{Shallow 850-1000mb}')
case 11
   [rhomap.eabs_p24p30m_cmfmc_pi24pi30m, pmap.eabs_p24p30m_cmfmc_pi24pi30m] = corrcoef_time_map(eabs_p24p30m,cmfmc_pi24pi30m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of |Q_{SPCAM} - Q_{CAM}|_{850-1000mb} v.s. MF_{Shallow 850-1000mb}')
case 12
   [rhomap.eabs_p24p30m_cmfmcdzm_pi24pi30m, pmap.eabs_p24p30m_cmfmcdzm_pi24pi30m] = corrcoef_time_map(eabs_p24p30m,cmfmcdzm_pi24pi30m,ilat,ilon,1);
   figure; hold on
   title('Time-Correlation of |Q_{SPCAM} - Q_{CAM}|_{850-1000mb} v.s. MF_{Deep 850-1000mb}')
end
%
%
plotmap(lat,lon,20,1); hold on;
% find insignificant index in rhomap and set to NaN
rhotmp = eval(['rhomap.' vname{icase}]);
insig = find(eval(['pmap.' vname{icase}])>0.05);
rhotmp(insig)=NaN;
m_contourf(lon(ilon),lat(ilat),rhotmp',10,'linestyle','none');
caxis([-1 1]);
m_contour(lon(ilon),lat(ilat),mean(eabs_p14p21m(ilon,ilat,:),3)',[1 1],'r','linewidth',2)
set(gcf,'position',[-5 419 1397 383]);
set(gcf,'color',[1 1 1])
set(gca,'fontsize',20)
%m_contour(lon(ilon),lat(ilat),mean(eabs_p14p24m(ilon,ilat,:),3)',[3 3],'r','linewidth',2)
%m_contour(lon(ilon),lat(ilat),mean(eabs_p24p30m(ilon,ilat,:),3)',[7 7],'r','linewidth',2)
%m_contour(lon(ilon),lat(ilat),mean(eabs_p21p30m(ilon,ilat,:),3)',[5 5],'r','linewidth',2)
%m_contour(lon(ilon),lat(ilat),mean(eabs_p14p30m(ilon,ilat,:),3)',[3 3],'r','linewidth',2)
% savefig('/projects/rsmas/kirtman/gchen/archive/matlab/figure/rho_map/E/p-ave/rhomap_e_preccdzm_t1t754_p1p30-ave.fig')
% save('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/p-ave/rhomap_p_map_e_eabs_var_p1p30-ave.mat','rhomap','pmap')
% save('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/p-ave/rhomap_p_map_e_eabs_cmfmcdzm_cmfmc_relhum_p-ave.mat','rhomap','pmap')

