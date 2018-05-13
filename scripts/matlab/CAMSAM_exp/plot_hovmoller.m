%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   plot p-wgt e eabs latitudinal hovmoller 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it=1:numel(itimes)
%   etmp(:,:,it)  = pwgtave(var{31}(:,:,:,itimes(it)),iilev,dlev);
%   eatmp(:,:,it) = pwgtave(var{32}(:,:,:,itimes(it)),iilev,dlev);
%   vtmp1(:,:,it)  = pwgtave(var{46}(:,:,:,itimes(it)),iilev,dlev);
   vtmp2(:,:,it)  = pwgtave(var{47}(:,:,:,itimes(it)),iilev,dlev);
end
stdmap.etmp=std(etmp,0,3);
stdmap.eatmp=std(eatmp,0,3);
e_lonave = squeeze(mean(etmp,1));
eabs_lonave = squeeze(mean(eatmp,1));
% 200 400 700 900
subplot(1,2,1)
contourf(squeeze(e_lonave),30,'linestyle','none')
colormap(cmap(1)); caxis([-2,2])
subplot(1,2,2)
contourf(squeeze(eabs_lonave(:,:)),30,'linestyle','none')
colormap(cmap(1)); caxis([-2,2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot p-wgt time-zonal-ave e vs var
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot(lat,mean(mean(etmp,1),3),'b'); hold on
plot(lat,mean(mean(eatmp,1),3),'g'); hold on
plot(lat,mean(mean(var{10},1),3),'r')
plot(lat,mean(mean(var{11},1),3),'m')
%plot(lat,nanmean(nanmean(var37,1),3),'k') % CLOUD
plot(lat,mean(mean(var{18},1),3),'k')
legend('EABS', 'PRECCDZM', 'PRECSH','SOLIN')
xlim([-90 90])
xlabel('latitude')
title('Conditioned on nonzero PRECCDZM')
savefig('e_eabs_preccdzm_xpt-ave.fig')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot p-wgt time-zonal-ave e vs var (conditioned on nonzero)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ix iy iz it] = findnzero(var{10},lat,lon,1,itimes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e eabs RELHUM latitudinal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zetmplat=zscore(nanmean(nanmean(etmp,1),3));
zeatmplat=zscore(nanmean(nanmean(eatmp,1),3));
zvtmplat =zscore(nanmean(nanmean(vtmp,1),3));
plot(lat,zetmplat,'b');hold on
plot(lat,zeatmplat,'g')
plot(lat,zvtmplat,'r')
legend('E', 'EABS', 'RELHUM')
xlim([-90 90])
xlabel('latitude')
savefig('e_eabs_relhum_xpt-ave.fig')
title('Conditioned on RELHUM nonzero')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e eabs preccdzm precsh latitudinal 
% (conditioned on preccdzm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it=1:numel(itimes)
   etmp(:,:,it)  = pwgtave(var{31}(:,:,:,itimes(it)),iilev,dlev);
   eatmp(:,:,it) = pwgtave(var{32}(:,:,:,itimes(it)),iilev,dlev);
%   var3(:,:,it)  = pwgtave(var{3}(:,:,:,itimes(it)),iilev,dlev); % RELHUM
%   var37(:,:,it)  = pwgtave(var{37}(:,:,:,itimes(it)),ilev,dlev);
end
ilat = [1:192];
var10 = nzero2nan(var{10},ilon,ilat,1,itimes);
var11 = nzero2nan(var{11},ilon,ilat,1,itimes);
% use the NaN in var10 to make etmp NaN at the same indices
var30     = var3 + var10; var30 = var30 - var10;
var80     = var8 + var10; var80 = var80 - var10; 
etmp0     = etmp + var10; etmp0 = etmp0 - var10;
etmp0lat  = nanmean(nanmean((etmp0),1),3);
eatmp0    = eatmp + var10; eatmp0 = eatmp0 - var10;
eatmp0lat = nanmean(nanmean((eatmp0 - var10),1),3);
var110    = var11 + var10; var110 = var110 - var10;
var110lat = nanmean(nanmean((var110 - var10),1),3);
var180    = var18 + var10; var180 = var180 - var10; % SOLIN
var180lat = nanmean(nanmean((var180 - var10),1),3);
%var370 = var37 + var10;
%var370 = nanmean(nanmean((var370 - var10),1),3);
%plot(lat,nanzscore(etmp0'),'b');hold on
doznorm=1;
if (doznorm)
   plot(lat,nanzscore(eatmp0'),'g');hold on
   plot(lat,nanzscore(nanmean(nanmean(var10,1),3)'),'r');hold on
   plot(lat,nanzscore(var110'),'m');
   plot(lat,nanzscore(var180'),'k');
%   plot(lat,nanzscore(var370'),'k');
%plot(lat,nanzscore(var180'),'k') ; % RELHUM doesn't change lat-shape when
else
   plot(lat,(etmp0'),'b');hold on
   plot(lat,(eatmp0'),'g');hold on
   plot(lat,(nanmean(nanmean(var10,1),3)'),'r');hold on
   plot(lat,(var110'),'m');hold off
end
legend('EABS', 'PRECCDZM', 'PRECSH', 'CLOUDFRAC')
xlim([-90 90])
xlabel('latitude')
title('Conditioned on PRECCDZM nonzero')
savefig('eabs_preccdzm_precsh_xpt-ave_znorm.fig')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% e eabs preccdzm precsh latitudinal 
% (conditioned on precsh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var10 = nzero2nan(var{10},ilon,ilat,1,itimes);
var11 = nzero2nan(var{11},ilon,ilat,1,itimes);
etmp1 = etmp + var11;
etmp1 = etmp1 - var11;
eatmp1 = eatmp + var11;
eatmp1 = eatmp1 - var11;
var111 = var11 + var11;
var111 = var111 - var11;
plot(lat,nanmean(nanmean(etmp1,1),3),'b');hold on
plot(lat,nanmean(nanmean(eatmp1,1),3),'g')
plot(lat,nanmean(nanmean(var10,1),3),'r')
plot(lat,nanmean(nanmean(var111,1),3),'m')
legend('E', 'EABS', 'PRECCDZM', 'PRECSH')
xlim([-90 90])
xlabel('latitude')
title('Conditioned on PRECSH nonzero')
%e_lonave = squeeze(mean(var{31},1));
%eabs_lonave = squeeze(mean(var{32},1));
iilev = [14 18 22 26];
for i=1:4
   subplot(2,4,i)
   contourf(squeeze(e_lonave(:,iilev(i),:)),30,'linestyle','none')
   colormap(cmap(1)); caxis([-10,10])
   subplot(2,4,i+4)
   contourf(squeeze(eabs_lonave(:,iilev(i),:)),30,'linestyle','none')
end
