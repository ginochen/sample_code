% 1. mean vertical wind shear (200 '400' 500 vs 850 '1000'), 400hPa is
% close to 6km which is where the storm shear below is estimated
% 2. SST DTCOND
% 3. time lag corr
% 4. LTS (lower atmospheric stability), theta_700 - theta_0 
% 5. zenith angle
%
%
% 3D-spatial-cross-correlation coefficient matrix time series
%
[R3,P3,ivvec] = corrcoef_3d(var,ivarsall,ilat,itimes);
%
% 2D-spatial-cross-correlation coefficient matrix time series
% (includes pwt'd average of 3D variables)
[R2,P2,v1]    = corrcoef_2d(var,ivarsall,ilat,ilon,ilev,dlev,itimes);
%
%for iv=1:numel(ivv);   vname{iv} = varname{ivv(iv)}; end
%
%
%   v(:,:,31) = var{31}(:,ilat,end,itime);
%   v(:,:,17) = var{18}(:,ilat,19,itime);
%   v(:,:,18) = var{18}(:,ilat,20,itime);
%   v(:,:,19) = var{18}(:,ilat,21,itime);
%   v(:,:,21) = ((var{21}(:,ilat,iilev(1),itime)-var{21}(:,ilat,iilev(3),itime))/(lev(iilev(3))-lev(iilev(1)))) + v(:,:,21); 
%   v(:,:,22) = ((var{22}(:,ilat,iilev(1),itime)-var{22}(:,ilat,iilev(3),itime))/(lev(iilev(3))-lev(iilev(1)))) + v(:,:,22); 
%   for j=[iilev(2),iilev(3)-1];
%% u mean vertical shear
%      v(:,:,21) = ((var{21}(:,ilat,j,itime)-var{21}(:,ilat,j+1,itime))/(lev(j)-lev(j+1))) + v(:,:,21); 
%% v mean vertical shear
%      v(:,:,22) = ((var{22}(:,ilat,j,itime)-var{22}(:,ilat,j+1,itime))/(lev(j)-lev(j+1))) + v(:,:,22); 
%      if (j==numel(lev)-1)
%         v(:,:,21) = v(:,:,21)/j; % average
%         v(:,:,22) = v(:,:,22)/j; 
%      end   
%   end
% mean vetical shear magnitude
%   v(:,:,23) = sqrt((v(:,:,21)).^2 + (v(:,:,22)).^2);    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  R table of 2D spatial correlation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          'DTCOND'  'OMEGAT'  'RELHUM'  'SPDT'    'e'        'e_abs'  'PRECC'  'PRECCDZM' 'PRECSH'  'SOLIN'
%'DTCOND'   1.0000   -0.4369    0.3473    0.4515   -0.0105    0.3123    0.7217    0.6381    0.4361    0.0636
%'OMEGAT'  -0.4369    1.0000   -0.5168   -0.7006   -0.5583   -0.6134   -0.5263   -0.5470   -0.1018   -0.0539
%'RELHUM'   0.3473   -0.5168    1.0000    0.4033    0.2725    0.4568    0.5444    0.5722    0.0904    0.0361
%'SPDT'     0.4515   -0.7006    0.4033    1.0000    0.8855    0.8102    0.5417    0.5568    0.1209    0.0171
%'e'       -0.0105   -0.5583    0.2725    0.8855    1.0000    0.7473    0.2347    0.2953   -0.0906   -0.0143
%'e_abs'    0.3123   -0.6134    0.4568    0.8102    0.7473    1.0000    0.5955    0.5932    0.1850    0.1021
%'PRECC'    0.7217   -0.5263    0.5444    0.5417    0.2347    0.5955    1.0000    0.9422    0.4489    0.0735
%'PRECCDZM' 0.6381   -0.5470    0.5722    0.5568    0.2953    0.5932    0.9422    1.0000    0.1269    0.0658
%'PRECSH'   0.4361   -0.1018    0.0904    0.1209   -0.0906    0.1850    0.4489    0.1269    1.0000    0.0422
%'SOLIN'    0.0636   -0.0539    0.0361    0.0171   -0.0143    0.1021    0.0735    0.0658    0.0422    1.0000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  R table of 3D spatial correlation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           'DTCOND'  'OMEGAT'  'RELHUM'  'SPDT'    'e'       'e_abs'
%'DTCOND'    1.0000   -0.1378    0.0068    0.1917   -0.4667   -0.0437
%'OMEGAT'   -0.1378    1.0000   -0.1795   -0.2753   -0.1561   -0.2068
%'RELHUM'    0.0068   -0.1795    1.0000    0.1072    0.0913    0.2189
%'SPDT'      0.1917   -0.2753    0.1072    1.0000    0.7732   -0.2247
%'e'        -0.4667   -0.1561    0.0913    0.7732    1.0000   -0.1824
%'e_abs '   -0.0437   -0.2068    0.2189   -0.2247   -0.1824    1.0000

%
%
% normalized dtcond
%   norm_dtcond(:,1)=(v1(:,1))/std(v1(:,1));
%   norm_dtcond(:,2)=(v1(:,2))/std(v1(:,2));
%   norm_dtcond(:,1) = (v1(:,1) - mean(v1(:,1)))/std(v1(:,1));
%   norm_dtcond(:,2) = (v1(:,2) - mean(v1(:,2)))/std(v1(:,2));
%   norm_dtcond(:,3) = mv2(v(:,:,3)); % e_abs
%  norm_dtcond(:,4) = abs(norm_dtcond(:,1)-norm_dtcond(:,2)); % e_norm
%
%
%ivars=[1,2,3,4,17,18,19,21,22,23,28];
ivars=[2,4,10,11,15,26,29,31]; %2,3: cam - crm, PRECC_CAM, PRECCDZM_CAM, PRECSH_CAM, OMEGAT_CAM, SOLIN_CAM, CMFDT_CAM, |crm-cam|
nlead=0; % lead 2 give the highest corr between PRECCDZM, and e (.3648) or e_abs (.5939)
itt=1;
for itime = itimes
   [R(:,:,itt), P(:,:,itt)] = corrcoef(horzcat([v1(:,ivars,itime), v1(:,3,itime)])); % make sure itimes is 1:700 no skipping indices
%   R(:,:,itt) = corrcoef(horzcat([v1(:,ivars,itime+nlead),v1(:,10,itime)]));
%   R(:,:,itt) = corrcoef(horzcat([v3(:,ivars,itime)]));
   itt=itt+1;
end
Rm=mean(R,3)
Rm(10,:)
%
itt=1;
for itime = itimes
% search indices
%   inzero0 = find(v1(:,4,itime)~=0); % find non-zero indices for PREC_CAM
%   inzero1 = find(v1(:,10,itime)~=0); % find non-zero indices for PREC_CAM
%   inzero2 = find(v1(:,11,itime)~=0); % find non-zero indices for PREC_CAM
%   innan{1} = find(isnan(v1(:,6,itime))==0); % find non-nan indices for CAPE_CAM
%   innan{2} = find(isnan(v1(:,7,itime))==0); % find non-nan indices for CAPE_SPCAM
%   innan{3} = find(isnan(v1(:,8,itime))==0); % find non-nan indices for CIN_CAM
%   innan{4} = find(isnan(v1(:,9,itime))==0); % find non-nan indices for CIN_SPCAM
%   innan{5} = find(isnan(v1(:,16,itime))==0);
%   innan{6} = find(isnan(v1(:,17,itime))==0);
%   ithres{1} = find(v1(:,3,itime)>1);  % e_abs (with threshold)
%   ithres{2} = find(v1(:,3,itime)>5);  % e_abs (with threshold)
%   ithres{3} = find(v1(:,3,itime)>10); % e_abs (with threshold)
%    ithres{1} = find(v1(:,27,itime)>200); % solar insolation
%
% v3 is a counter map to mark where the threshold is reached
%   v2=zeros(nlon*nlat,3);
%   v2(ithres{1},1)=1;
%   v2(ithres{2},2)=1;
%   v2(ithres{3},3)=1;
%   v3=reshape(v2,nlon,nlat,3);
%
% correlation values (average over all 700 times, cutoff_lat = 65,
% conditioned on PRECC_inzero doesn't lead to better corr for all corr)
%   [rho(itt,1),  pval(itt,1)]  =   corr(norm_dtcond(:,1),norm_dtcond(:,2)); % relative dtcond (CAM vs SPCAM) <---- .55 
%   [rho(itt,2),  pval(itt,2)]  =   corr(norm_dtcond(:,1),norm_dtcond(:,3)); % norm_dtcond vs e_norm_dtcond (CAM vs e_SPCAM-CAM) <----- .43
%   [rho(itt,1),  pval(itt,1)]  = corr(v1(:,1,itime),v1(:,2,itime)); % norm_dtcond vs e_norm_dtcond (CAM vs e_SPCAM-CAM) <----- .43
%   [rho(itt,2),  pval(itt,2)]  = corr(v1(:,1,itime),abs(v1(:,3,itime))); % norm_dtcond vs e_norm_dtcond (CAM vs e_SPCAM-CAM) <----- .43
%   [rho(itt,3),  pval(itt,3)]  = corr(norm_dtcond(:,2),norm_dtcond(:,3)); % norm_dtcond vs e_norm_dtcond (SPCAM vs e_SPCAM-CAM) <----- .72
%   [rho(itt,4),  pval(itt,4)]  = corr(v1(:,4,itime),v1(:,3,itime)); % PREC vs e_norm_dtcond (CAM vs e_SPCAM-CAM) <--- .54(e_norm) .57 (e_abs)
%   [rho(itt,5),  pval(itt,5)]  = corr(v1(:,5,itime),v1(:,3,itime)); % PREC vs e_norm_dtcond (CAM vs e_SPCAM-CAM) <--- .54
%   [rho(itt,6),  pval(itt,6)]  = corr(v1(innan{1},6),v1(innan{1},3,itime)); % CAPE vs e_norm_dtcond (CAM   vs e_SPCAM-CAM) <--- .28
%   [rho(itt,7),  pval(itt,7)]  = corr(v1(innan{2},7),v1(innan{2},3,itime)); % CAPE vs e_norm_dtcond (SPCAM vs e_SPCAM-CAM) <---- .23
%   [rho(itt,8),  pval(itt,8)]  = corr(v1(innan{3},8,itime),v1(innan{3},3,itime)); % CIN vs e_norm_dtcond (CAM   vs e_SPCAM-CAM) <--- no corr
%   [rho(itt,9),  pval(itt,9)]  = corr(v1(innan{4},9,itime),v1(innan{4},3,itime)); % CIN vs e_norm_dtcond (SPCAM vs e_SPCAM-CAM) <--- no corr
%   [rho(itt,10), pval(itt,10)] = corr(v1(ithres{1},4),norm_dtcond(ithres{1},3));
%   [rho(itt,11), pval(itt,11)] = corr(v1(inzero,4),norm_dtcond(inzero,3));
%   [rho(itt,8),  pval(itt,8)]  = corr(v1(:,10,itime),v1(:,1,itime)); % PRECCDZM vs DTCOND (CAM vs CAM) <--- .47
%   [rho(itt,8),  pval(itt,8)]  = corr(v1(:,10,itime),v1(:,3,itime)); % PRECCDZM vs e (CAM vs |SPCAM-CAM|) <--- .30
%   [rho(itt,8),  pval(itt,8)]  = corr(v1(:,10,itime),v1(:,3,itime+2)); % PRECCDZM vs e (CAM vs SPCAM-CAM) <--- .3648
%   [rho(itt,8),  pval(itt,8)]  = corr(v1(:,10,itime),v1(:,3,itime)); % PRECCDZM vs e_abs (CAM vs |SPCAM-CAM|) <--- .56
%   [rho(itt,8),  pval(itt,8)]  = corr(v1(:,10,itime),v1(:,3,itime+2)); % PRECCDZM vs e_abs (CAM vs |SPCAM-CAM|) <--- .5939
%ePRECCDZM e_abs
%e_abs sign is wrong predicting PRECCDZM
%all physics process going into DTCOND (list)
%PRECL
%DTCOND_SPCAM PRECCDZM_CAM
%horiz conv at the bdry layer (Quo)
%f(SOLIN)=[SHEAR DTCOND] related to eDTCOND
%   [rho(itt,9),  pval(itt,9)]  = corr(v1(:,11,itime),v1(:,3,itime)); % PRECSH vs e_abs (CAM vs |SPCAM-CAM|) <---- .21
%   [rho(itt,10), pval(itt,10)] = corr(v1(inzero1,10,itime),v1(inzero1,3,itime)); % PRECCDZM vs e_abs (CAM vs |SPCAM-CAM|) <--- .51 (worse than no condition)
%   [rho(itt,11), pval(itt,11)] = corr(v1(inzero2,11,itime),v1(inzero2,3,itime)); % PRECSH vs e_abs (CAM vs |SPCAM-CAM|) <---- .32 (better than no condition)
%   [rho(itt,12), pval(itt,12)] = corr(v1(:,12,itime),v1(:,3,itime)); % PRECL vs e_abs (CAM vs |SPCAM-CAM|) <---- .16
%   [rho(itt,13), pval(itt,13)] = corr(v1(:,13,itime),v1(:,3,itime)); % PRECT vs e_abs (CAM vs |SPCAM-CAM|) <---- .48
%   [rho(itt,14), pval(itt,14)] = corr(v1(:,14,itime),v1(:,3,itime));%-v1(:,1,itime)); % OMEGA vs e_abs <----- -.44 (time series of rho overlaps with OMEGAT)
%   [rho(itt,15), pval(itt,15)] = corr(v1(:,15,itime),v1(:,3,itime));%-v1(:,1,itime)); % OMEGAT vs e_abs <---- -.45
%   [rho(itt,16), pval(itt,16)] = corr(v1(innan{5},16,itime),v1(innan{5},3,itime)); % TMQ vs e_abs
%   [rho(itt,17), pval(itt,17)] = corr(v1(innan{6},17,itime),v1(innan{6},3,itime)); % ATMEINT vs e_abs
%   [rho(itt,18), pval(itt,18)] = corr(v1(:,18,itime),v1(:,3,itime)); % RELHUM vs e_abs
%   [rho(itt,19), pval(itt,19)] = corr(v1(:,21,itime),v1(:,3,itime)); % RELHUM vs e_abs <------ .22 (higher than any level above (20,19))
%    [rho(itt,20),pval(itt,20)] = corr(v1(:,27,itime),v1(:,3,itime)); % SOLIN vs e <---- .33 (don't do thresholding, it's much smaller conditioned on larger (>200) SOLIN)
%    [rho(itt,21),pval(itt,21)] = corr(v1(:,28,itime),v1(:,3,itime)); % SRFRAD vs e <---- .27 
%    SOLIN 
    itt=itt+1;
end
%%%%%%%%%%%%%%%%%%%%%%%
% 2-D maps of variables
%%%%%%%%%%%%%%%%%%%%%%%
   % subplots of vertical time mean maps
   figure;
   contour_z_subplot(mean(var{31},4),[1:3:15],2,cmap(1),[-2 2])
   figure;
   contour_z_subplot(mean(var{31},4),[18:30],2,cmap(1),[-4 4])
   figure;
   contourf(mean(mean(var{31}(:,:,14:22,:),3),4)',100,'linestyle','none'); colormap(cmap(1)); caxis([-4 4]); colorbar
   % plevel 232 to 763 according to corrcoef of e and eabs with
   % RELHUM, PRECCDZM, OMEGAT doing better there (rho ~ .5)
   contourf(v(:,:,4)',100,'linestyle','none'); hold on;  contour(v3(:,:,1)','k.'); hold off;
   subplot(3,3,1)
   contourf(v(:,:,1)',100,'linestyle','none'); caxis([-10 60]); colorbar; title('P-weighted DTCOND (CAM)')
   subplot(3,3,2)
   contourf(v(:,:,4)',100,'linestyle','none'); colorbar; title('Convective Precipitation Rate (mm/day) (CAM)');
   subplot(3,3,3)
   contourf(v(:,:,6)',100,'linestyle','none'); colorbar; title('CAPE (CAM)');
   subplot(3,3,4)
   contourf(v(:,:,2)',100,'linestyle','none'); caxis([-10 60]); colorbar; title('P-weighted DTCOND (SPCAM)');
   subplot(3,3,5)
   contourf(v(:,:,5)',100,'linestyle','none');  colorbar; title('Convective Precipitation Rate (mm/day) (SPCAM)');
   subplot(3,3,6)
   contourf(v(:,:,7)',100,'linestyle','none'); colorbar; title('CAPE (SPCAM)');
   subplot(3,3,7)
   contourf(v(:,:,3)',100,'linestyle','none'); caxis([-10 60]); colorbar; title('P-weighted eDTCOND (SPCAM-CAM)');
%   subplot(3,3,6)
%   hist(v1(:,3),100); title('P-weighted eDTCOND')
   subplot(3,3,8)
   contourf(v(:,:,8)',100,'linestyle','none'); colorbar; title('CIN (CAM)');
   subplot(3,3,9)
   contourf(v(:,:,9)',100,'linestyle','none'); colorbar; title('CIN (SPCAM)');
%
% density contour over scatter plot of two variables
   figure;
   subplot(1,3,1)
   nbins = [100,100];
   [bin_counts bin_centers] = hist3([v1(:,1),v1(:,4)],nbins); % DTCOND_CAM vs PRECC_CAM
   var_axis(:,1)=linspace(min(v1(:,1)),max(v1(:,1)),nbins(1));
   var_axis(:,2)=linspace(min(v1(:,4)),max(v1(:,4)),nbins(2));
   contourf(var_axis(:,1),var_axis(:,2),log(1+bin_counts)',100,'linestyle','none'); colorbar% log scale of bin_counts, since the density of too focused locally
   xlim([-10 50]); ylim([0 240]);
   xlabel('Pressure weighted average DTCOND (K/day)'); ylabel('Convective Precipitation Rate (mm/day)'); title('CAM')
   subplot(1,3,2)
   [bin_counts bin_centers] = hist3([v1(:,2),v1(:,5)],nbins); % DTCOND_SPCAM vs PRECC_SPCAM
   var_axis(:,3)=linspace(min(v1(:,2)),max(v1(:,2)),nbins(1));
   var_axis(:,4)=linspace(min(v1(:,5)),max(v1(:,5)),nbins(2));
   contourf(var_axis(:,3),var_axis(:,4),log(1+bin_counts)',100,'linestyle','none'); colorbar% log scale of bin_counts, since the density of too focused locally
   xlim([-10 50]); ylim([0 240]);
   xlabel('Pressure weighted average DTCOND (K/day)'); ylabel('Convective Precipitation Rate (mm/day)'); title('SPCAM')
   subplot(1,3,3)
   [bin_counts bin_centers] = hist3([v1(:,3),v1(:,4)],nbins); % eDTCOND vs PRECC_CAM
   var_axis(:,5)=linspace(min(v1(:,3)),max(v1(:,3)),nbins(1));
   var_axis(:,6)=linspace(min(v1(:,4)),max(v1(:,4)),nbins(2));
   contourf(var_axis(:,5),var_axis(:,6),log(1+bin_counts)',100,'linestyle','none'); colorbar% log scale of bin_counts, since the density of too focused locally
%   xlim([-10 100]); ylim([0 240]);
   xlabel('Pressure weighted average eDTCOND (K/day)'); ylabel('Convective Precipitation Rate (mm/day)'); title('SPCAM - CAM')
%
% density contour over scatter plot with PRECC_CAM non-zero indices
   figure;
   subplot(2,2,1)
   nbins = [100,100];
   [bin_counts bin_centers] = hist3([v1(inzero,1),v1(inzero,4)],nbins); % DTCOND_CAM vs PRECC_CAM (conditioned on nonzero index of PRECC_CAM)
   var_axis(:,1)=linspace(min(v1(inzero,1)),max(v1(inzero,1)),nbins(1));
   var_axis(:,2)=linspace(min(v1(inzero,4)),max(v1(inzero,4)),nbins(2));
   contourf(var_axis(:,1),var_axis(:,2),log(1+bin_counts)',100,'linestyle','none'); colorbar% log scale of bin_counts, since the density of too focused locally
   xlim([-10 50]); ylim([0 240]);
   xlabel('P-weighted DTCOND (K/day)'); ylabel('Convective Precipitation Rate (mm/day)'); title('CAM (conditioned)')
   subplot(2,2,2)
   [bin_counts bin_centers] = hist3([v1(:,2),v1(:,5)],nbins); % DTCOND_SPCAM vs PRECC_SPCAM
   var_axis(:,3)=linspace(min(v1(:,2)),max(v1(:,2)),nbins(1));
   var_axis(:,4)=linspace(min(v1(:,5)),max(v1(:,5)),nbins(2));
   contourf(var_axis(:,3),var_axis(:,4),log(1+bin_counts)',100,'linestyle','none'); colorbar% log scale of bin_counts, since the density of too focused locally
   xlim([-10 50]); ylim([0 240]);
   xlabel('P-weighted DTCOND (K/day)'); ylabel('Convective Precipitation Rate (mm/day)'); title('SPCAM')
   subplot(2,2,3)
   [bin_counts bin_centers] = hist3([v1(inzero,3),v1(inzero,4)],nbins); % eDTCOND vs PRECC_CAM (conditioned)
   var_axis(:,5)=linspace(min(v1(inzero,3)),max(v1(inzero,3)),nbins(1));
   var_axis(:,6)=linspace(min(v1(inzero,4)),max(v1(inzero,4)),nbins(2));
   contourf(var_axis(:,5),var_axis(:,6),log(1+bin_counts)',100,'linestyle','none'); colorbar% log scale of bin_counts, since the density of too focused locally
%   xlim([-10 100]); ylim([0 240]);
   xlabel('P-weighted eDTCOND (K/day)'); ylabel('Convective Precipitation Rate (mm/day)'); title('SPCAM - CAM (conditioned)')
   subplot(2,2,4)
   nbins = [100,100];
   [bin_counts bin_centers] = hist3([v1(inzero,1),v1(inzero,2)],nbins); % DTCOND_CAM vs PRECC_CAM
   var_axis(:,7)=linspace(min(v1(inzero,1)),max(v1(inzero,1)),nbins(1));
   var_axis(:,8)=linspace(min(v1(inzero,2)),max(v1(inzero,2)),nbins(2));
   contourf(var_axis(:,7),var_axis(:,8),log(1+bin_counts)',100,'linestyle','none'); colorbar% log scale of bin_counts, since the density of too focused locally
%   xlim([-10 50]); ylim([0 240]);
   xlabel('P-weighted DTCOND (CAM)'); ylabel('P-weighted DTCOND (SPCAM)'); title('CAM vs SPCAM (conditioned)')
%
% scatter plot of two dtcond conditioned on PRECC_CAM nonzero indices
   figure;  
   subplot(2,1,1); histogram(v1(:,1),dtcondbin,'Normalization','probability','EdgeColor','none'); title('P-weighted DTCOND (CAM)'); xlim([-20 20]); ylabel('Probability')
   subplot(2,1,2); histogram(v1(:,2),dtcondbin,'Normalization','probability','EdgeColor','none'); title('P-weighted DTCOND (SPCAM)'); xlim([-20 20]); ylabel('Probability')
   xlabel('P-weighted eDTCOND (K/day)'); ylabel('Convective Precipitation Rate (mm/day)'); title('CAM - SPCAM (conditioned)')
   scatter(v1(inzero,1),v1(inzero,2)); xlabel('P-weighted DTCOND (CAM)'); ylabel('P-weighted DTCOND (SPCAM)'); title('CAM')
   scatter(mv2(v(:,:,3)),mv2(v(:,:,4))); xlabel('Pressure weighted average DTCOND'); ylabel('PREC'); title('CAM')
%   subplot(2,1,2)
%   scatter(mv2(v(:,:,2)),mv2(v(:,:,5))); xlabel('Pressure weighted average DTCOND'); ylabel('PREC'); title('SPCAM')
%
% plot time series of e_abs e_norm e_dtcond
end
