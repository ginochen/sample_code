% merge all variable time and spatial index into one dimension
fname = {'div', 'vort', 'u', 'w', 'T', 'qT', 'qc'};
casename = '15km_lead8-24';
%casename = '15km';
ii(1:4)=1;ii2(1:4)=0;idx.u = 1; idx.w = 2; idx.T = 3; idx.qT = 4; idx.qc = 5; idx.rho = 6;
docamsclvar=0;
docrmsclvar=1;
%lagind = [48:-1:-48]; % 12 hours lead
lagind=[-8:-1:-24];
tind = [11:10:491];
for it = tind
   %load(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km/var_PC1_' num2str(it) '.mat'],'var','div','vort','n2crm','FFTke','parm'); 
   %load(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km/var_PC1_' num2str(it) '.mat'],'vort','div'); 
   %load(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/' casename '/var_PC1_' num2str(it) '.mat'],'var','div','vort','n2crm','FFTke','qr'); 
   load(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/' casename '/var_PC1_' num2str(it) '.mat'],'qr'); 
   load(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km/var_PC1_' num2str(it) '.mat'],'parm'); % used by 15km_lead 
   for ib=1:4
      ii2(ib) = ii(ib)+numel(parm.iCAMll{ib})-1; 
      if (docamsclvar)
         iilon{ib} = parm.iillt_zm_PC{ib}(1,parm.iCAMll{ib});
         iilat{ib} = parm.iillt_zm_PC{ib}(2,parm.iCAMll{ib});
         for i = 1:numel(iilon{ib})
            dtcond{ib}(1:30,:,ii(ib)+i-1) = varCAM{idv.dt}(iilon{ib}(i),iilat{ib}(i),1:30,it-lagind-1); precc{ib}(:,ii(ib)+i-1)  = varCAM{idv.precc}(iilon{ib}(i),iilat{ib}(i),it-lagind-1);
            sppflx{ib}(:,ii(ib)+i-1) = varCAM{idv.sppflx}(iilon{ib}(i),iilat{ib}(i),end,it-lagind-1);   spdt{ib}(1:30,:,ii(ib)+i-1) = varCAM{idv.spdt}(iilon{ib}(i),iilat{ib}(i),1:30,it-lagind-1);
         end
      end
      if (docrmsclvar)
         %u{ib}(:,:,:,ii(ib):ii2(ib)) = var{ib}(:,:,:,:,idx.u);   w{ib}(:,:,:,ii(ib):ii2(ib)) = var{ib}(:,:,:,:,idx.w);   T{ib}(:,:,:,ii(ib):ii2(ib)) = var{ib}(:,:,:,:,idx.T);
         %qT{ib}(:,:,:,ii(ib):ii2(ib)) = var{ib}(:,:,:,:,idx.qT); qc{ib}(:,:,:,ii(ib):ii2(ib)) = var{ib}(:,:,:,:,idx.qc); rho{ib}(:,:,:,ii(ib):ii2(ib)) = var{ib}(:,:,:,:,idx.rho);
         %vortTime{ib}(:,:,:,ii(ib):ii2(ib)) = vort{ib};          divTime{ib}(:,:,:,ii(ib):ii2(ib)) = div{ib};            FFTkeTime{ib}(:,:,:,ii(ib):ii2(ib)) = FFTke{ib};
         %iCAMll{ib}(ii(ib):ii2(ib)) = parm.iCAMll{ib};           n2{ib}(:,:,:,ii(ib):ii2(ib)) = n2crm{ib};               
         qrTime{ib}(:,:,ii(ib):ii2(ib)) = qr{ib};
      end
      ii(ib) = ii2(ib)+1;
   end
end
vort = vortTime; div = divTime; FFTke = FFTkeTime; qr = qrTime; 
% for all time 11:10:491
save('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/camsclvar_alltime_PC1.mat','dtcond','spdt','sppflx','precc','lagind')
save(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/' casename '/div_all_PC1.mat'],'div')
save(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/' casename '/vort_all_PC1.mat'],'vort')
save(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/' casename '/u_all_PC1.mat'],'u')
save(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/' casename '/w_all_PC1.mat'],'w')
save(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/' casename '/T_all_PC1.mat'],'T')
save(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/' casename '/qT_all_PC1.mat'],'qT')
save(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/' casename '/qc_all_PC1.mat'],'qc')
save(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/' casename '/qr_all_PC1.mat'],'qr')
%
% load and concatenate the lag-lead dimension of the variable 
iv=2; fname{iv}
eval(sprintf('%s1 = load(''/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km/%s_all_PC1.mat'',''%s'')',fname{iv}, fname{iv}, fname{iv}));
eval(sprintf('%s2 = load(''/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km_lead8-24/%s_all_PC1.mat'',''%s'')',fname{iv}, fname{iv}, fname{iv}));
eval(sprintf('for ib=1:4; %s{ib} = cat(3,%s1.%s{ib},%s2.%s{ib}); end',fname{iv}, fname{iv}, fname{iv}, fname{iv}, fname{iv}));
%s2 =load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km_lead8-24/%s_all_PC1.mat','div'); 
vort0 = load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km/vort_all_PC1.mat','vort'); 
u0 = load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km/u_all_PC1.mat','u'); 
w0 = load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km/w_all_PC1.mat','w'); 
T0 = load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km/T_all_PC1.mat','T'); 
qT0 = load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km/qT_all_PC1.mat','qT'); 
qc0 = load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km/qc_all_PC1.mat','qc'); 
for ib = 1:4
%   divtmp{ib} = cat(3,div0.div{ib},div{ib});
%   vorttmp{ib} = cat(3,vort0.vort{ib},vort{ib});
%   utmp{ib} = cat(3,u0.u{ib},u{ib});
%   wtmp{ib} = cat(3,w0.w{ib},w{ib});
%   Ttmp{ib} = cat(3,T0.T{ib},T{ib});
%   qTtmp{ib} = cat(3,qT0.qT{ib},qT{ib});
%   qctmp{ib} = cat(3,qc0.qc{ib},qc{ib});
end
save('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km_lag7-lead24/div_all_PC1.mat')
save('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km_lag7-lead24/vort_all_PC1.mat')

% concat threshold indices in time
nthr = 4;
thres_bigPC = [20:10:60]; 
ii(1:4,1:nthr)=1; % counter for iCAMll_thres
jj(1:4)=1; % counter for iCAMll
for it = 11:10:491
   load(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km/var_PC1_' num2str(it) '.mat'],'parm'); 
   for ib=1:4
      for ithr = 1:nthr
         ii2(ib,ithr) = ii(ib,ithr)+numel(nonzeros(parm.iCAMll_thres{ib}(:,ithr)))-1;
         iCAMll_thres{ib}(ii(ib,ithr):ii2(ib,ithr),ithr) = nonzeros(parm.iCAMll_thres{ib}(:,ithr)) + jj(ib) - 1 ;
         ii(ib,ithr) = ii2(ib,ithr)+1;
      end
      jj(ib) = length(parm.iCAMll{ib})+jj(ib);
   end
end
save('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km/thres_index.mat','iCAMll_thres','thres_bigPC','iCAMll','-append')

%load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km/var_allTime_PC1.mat','iCAMll_thres','u','w','vortTime','FFTkeTime','divTime')
%
% plot zone
run ~/scripts/matlab/CAMSAM_exp/main_stoch
spd = 86400;
load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/param.mat','lagind','lev','ilon','ilat')
x = [1:nx/2] - 1; % -1 to start at zero
t = -[nlag:-1:-nlag]*0.5; % 0.5 hrs
idv.dt=1;      varname{idv.dt}     = 'DTCOND'; varCAM{idv.dt}     = loadvar(varname{idv.dt},camCase,spd,1); %(K/day) <----- only moist parameterized heating
idv.spdt=30;   varname{idv.spdt}   = 'SPDT';   varCAM{idv.spdt}   = loadvar(varname{idv.spdt},spcamCase,spd,0);% (K/day) T tendency due to CRM (without QRL+QRS radiative heating)
idv.precc=4;   varname{idv.precc}  = 'PRECC';  varCAM{idv.precc}  = loadvar(varname{idv.precc},camCase,spd*1000,1);% (mm/day) 
idv.sppflx=46; varname{idv.sppflx} = 'SPPFLX'; varCAM{idv.sppflx} = loadvar(varname{idv.sppflx},spcamCase,spd*1000,0); % Precip flux for CRM
%
%
figure
ithr=0; dokestd=0; for ib=1:4; subplot(4,1,ib); if ithr==0; idx_thr = ':'; else; idx_thr = nonzeros(iCAMll_thres{ib}(:,ithr)); end; yyaxis left; contourf(squeeze(median(FFTke{ib}(2:10,1,:,   idx_thr),4)),30,'linestyle','none'); colorbar; caxis([0 0.02]); yyaxis right; plot(squeeze(median(mean(mean(u{ib}(:,:,:,idx_thr).^2+w{ib}(:,:,:,idx_thr).^2,1),2),4))); hold off; disp(size(vort{ib}(:,:,:,idx_thr),4)); end
%
%
figure
ithr=0; for ib=1:4; subplot(4,1,ib); if ithr==0; idx_thr = ':'; else; idx_thr = nonzeros(iCAMll_thres{ib}(:,ithr)); end;  plot(squeeze(median(sum(sum(vort{ib}(:,:,:,idx_thr).^2,1),2),4))); hold off; disp(size(vort{ib}(:,:,:,idx_thr),4)); end
%
% plot circulation scalar (threshold (20:10:100) of PC1)
figure; 
ithr=0; for ib=1:4; subplot(4,1,ib); if ithr==0; idx_thr = ':'; else; idx_thr = nonzeros(iCAMll_thres{ib}(:,ithr)); end; plot(squeeze(median(abs(sum(sum(vort{ib}(:,:,:,idx_thr),1),2)),4)));    hold on; plot(squeeze(median(abs(sum(sum(vort{ib}(:,:,:,idx_thr),1),2)),4))); hold off; disp(size(vort{ib}(:,:,:,idx_thr),4)); end
%
%

for ill = iCAMll{ib}%1:2=lon,lat 
   iilon{ib} = parm.iillt_zm_PC{ib}(1,iCAMll{ib});
   iilat{ib} = parm.iillt_zm_PC{ib}(2,iCAMll{ib});
end
sub2ind(size(varCAM{idv.dt}),iilon{ib},iilat{ib},1:30,it
figure; yyaxis left; contourf(t,lev,(squeeze(varCAM{idv.dt}(ilon,ilat,1:30,it-lagind-1))),'linestyle','none'); set(gca,'YDir','Reverse'); yyaxis right;  plot(t,squeeze(varCAM{idv.precc}(ilon,ilat,it-lagind-1)),'w','linewidth',2.5); title('CAM heating & rain rate')

figure; yyaxis left; contourf(t,lev,(squeeze(varCAM{idv.spdt}(ilon,ilat,1:30,it-lagind-1))),'linestyle','none'); set(gca,'YDir','Reverse'); yyaxis right;  plot(t,squeeze(varCAM{idv.sppflx}(ilon,ilat,end,it-lagind-1)),'w','linewidth',2.5); title('SPCAM heating & rain rate')

% plot large-scale kinetic energy growth
figure; 
ib=4;plot(squeeze(median(FFTke{ib}(2,1,:,:),4)))

% plot large-scale rain-fall mxing ratio growth <--- do fft on crm precip <---- no! rainfall is very discontinuous, flat spectrum from fft
% so plot its max variance mean svd mode <--- result shows 25 points decaying slowly to 22 points, so rainfall still exists
for ib=1:4; for i=1:size(qr{ib},2); for j=1:size(qr{ib},3);  qrr{ib}(i,j) = numel(nonzeros(qr{ib}(:,i,j))); end; end; [svec, sval, pc]=svdVar(qrr{ib}); ivec= 1; plot(svec(:,ivec)*median(pc(:,ivec))); pause; end
% plot median lag-lead time series 
for ib=1:4; for i=1:size(qr{ib},2); for j=1:size(qr{ib},3);  qrr{ib}(i,j) = numel(nonzeros(qr{ib}(:,i,j))); end; end; plot(median(qrr{ib},2)); pause; end
% plot actual rainfall lag-lead time series
for ib=1:4; for i=1:size(qr{ib},2); for j=1:size(qr{ib},3);  qrr{ib}(i,j) = mean(qr{ib}(:,i,j)); end; end; plot(median(qrr{ib},2)); pause; end

% plot large-scale (>100km^2) rainfall pattern 

% plot circulation (sum of vort) <---- see if it's steady after 3.5
% hrs


% Calculate large-scale MCS index and see if those time-space indices
% are on PC1 indices
