% define "it" before running this script
doloadvar = 0
if (doloadvar)
   run ~/scripts/matlab/startup.m
   load(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/15km/var_PC1_' num2str(it) '.mat'],'var','div','vort','n2crm','FFTke','parm');
   load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/param.mat'); 
   ib = 1;
   basin = {'IND'    'WPAC'    'EPAC'    'ATL'    'ALL'};
   Case.dir = '/projects/rsmas/kirtman/gchen/cesm_spcam/archive/spcam_actual_m2005_f09f09_branch/atm/rest/';
   load('param')
   %nxy = length(iCAMll); % # of lat-lon indeices at it-timestep
   km2m=1000;
   xwavelength = nx*dx./(1:nx/2)/km2m; % a_0 + b_0 + \sum ( a_n * cos(2pi/(L/n)*x) + b_n * sin(2pi/(L/n)*x) ), where L=nx*dx, n=1:nx/2
                           % 128.0000   64.0000   42.6667   32.0000 25.6000   21.3333   18.2857   16.0000
                           % 14.2222    12.8000   11.6364   10.6667 9.8462    9.1429    8.5333    8.0000 km
   zwavelength = parm.nzi*parm.dzi./(1:nzi/2)/km2m; % the entire height is 28km, first wave covers entire 28km with one wavelength, so zwavelength=nzi*dzi
   freq = 1./xwavelength; % [km^-1]
   x = [1:nx/2] - 1; % -1 to start at zero
   t = -[nlag:-1:-nlag]*0.5; % 0.5 hrs
   nt = length(t);
   %z = 0:nzi/2-1;
   izcoef = 1;
   ixcoef = 1;
   ixcoefs = ixcoef+1:1:ixcoef+16;
   spd = 86400; %seconds per day
   %%%%%%%%%%%%%%%%%%%%%%%%
   %   Load 2D CAM Var    %
   %%%%%%%%%%%%%%%%%%%%%%%%
   %for itlag = 1:nt
   run ~/scripts/matlab/CAMSAM_exp/main_stoch % includes camCase and spcamCase
   % cam   .rh0. var starts from 0001-01-14-03600 to 0001-01-29-64800
   % cam   .r.   var starts from 0001-01-14-03600 to 0001-01-29-64800
   % spcam .rh0. var starts from 0001-01-14-03600 to 0001-01-29-63000
   % spcam .r.   var starts from 0001-01-14-01800 to 0001-01-29-63000 <--- CRM_U CRM_W 
   % loadvar is to load .rh0. files, thus when comparing with spcam .r.
   % files, a timestep need to be removed when plotting
   idv.dt=1;  varname{idv.dt} = 'DTCOND';   varCAM{idv.dt}   = loadvar(varname{idv.dt},camCase,spd,1); %(K/day) <----- only moist parameterized heating
   idv.spdt=30; varname{idv.spdt} = 'SPDT'; varCAM{idv.spdt} = loadvar(varname{idv.spdt},spcamCase,spd,0);% (K/day) T tendency due to CRM (without QRL+QRS radiative heating)
   idv.precc=4;     varname{idv.precc} = 'PRECC';       varCAM{idv.precc}   = loadvar(varname{idv.precc},camCase,spd*1000,1);% (mm/day) 
   idv.sppflx=46; varname{idv.sppflx} = 'SPPFLX';   varCAM{idv.sppflx} = loadvar(varname{idv.sppflx},spcamCase,spd*1000,0); % Precip flux for CRM
end
%for itlag = 1:nt
   %prec(:,:,itlag)  = ncread([Caserh0.dir Caserh0.name{it-lagind(itlag)}],'PRECC',[1 1],[inf inf]);
   %dtcond(1:nlon,1:nlat,1:30,itlag)  = ncread([Caserh0.dir Caserh0.name{it-lagind(itlag)}],'DTCOND',[1 1 1],[inf inf inf]);
   %dtcond(1:nlon,1:nlat,1:30,itlag)  = ncread([Caserh0.dir Caserh0.name{it+itlag}],'DTCOND',[1 1 1],[inf inf inf]);
   %spdt(1:nlon,1:nlat,1:30,itlag)  = ncread([Caserh0.dir Caserh0.name{it-lagind(itlag)}],'SPDT',[1 1 1],[inf inf inf]);
   %sppflx(1:nlon,1:nlat,1:30,itlag)  = ncread([Caserh0.dir [Caserh0.name { num2str(it+itlag) }],'SPPFLX',[1 1 1],[inf inf inf]);
   %crm_qr(1:nlon,1:nlat,:,itlag)  = ncread([Case.dir Case.name{it-lagind(itlag)}],'CRM_QR',[1 1 1],[inf inf inf]); % from spcam .r. data 
   %preccdzm(1:nlon,1:nlat,itlag)  = ncread([Caserh0.dir Caserh0.name{it-lagind(itlag)}],'PRECCDZM',[1 1],[inf inf]);
   %qrl(1:nlon,1:nlat,1:30,itlag)  = ncread([Caserh0.dir Caserh0.name{it-lagind(itlag)}],'QRL',[1 1 1],[inf inf inf]);
%end
%crm_qr = reshape(crm_qr,nlon,nlat,nx,nz,2*nlag+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the var on only the PC1 lon-lat indices   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for ib = 1:4;
%   for i = 1:length(iillt_zm_PC{ib}(1,:)); 
%      dtcond_PC1{ib}(1:30,1:nt,i) = squeeze(dtcond(iillt_zm_PC{ib}(1,i),iillt_zm_PC{ib}(2,i),1:30,1:nt)); 
%   end
%end  
%prec  = ncread([Case.dir Case.name{it-lagind(iit)}],'PREC',[1 1 1],[inf inf inf]);
%
idx.u = 1; idx.w = 2; idx.T = 3; idx.qT = 4; idx.qc = 5; idx.rho = 6;
%
% plot enstrophy (rotational tendency) what does increasing enstrophy
% mean? more smaller scale rotation?
figure
ithr=0; for ib=1:4; subplot(4,1,ib); if ithr==0; idx_thr = ':'; else; idx_thr = nonzeros(parm.iCAMll_thres{ib}(:,ithr)); end; plot(squeeze(mean(sum(sum(vort{ib}(:,:,:,idx_thr).^2,1),2),4))); hold on; plot(squeeze(median(sum(sum(vort{ib}(:,:,:,idx_thr).^2,1),2),4))); hold off; disp(size(vort{ib}(:,:,:,idx_thr),4)); end
%
% plot circulation scalar (threshold (20:10:100) of PC1)
ithr=0; for ib=1:4; subplot(4,1,ib); if ithr==0; idx_thr = ':'; else; idx_thr = nonzeros(parm.iCAMll_thres{ib}(:,ithr)); end; plot(squeeze(mean(abs(sum(sum(vort{ib}(:,:,:,idx_thr),1),2)),4))); hold on; plot(squeeze(median(abs(sum(sum(vort{ib}(:,:,:,idx_thr),1),2)),4))); hold off; disp(size(vort{ib}(:,:,:,idx_thr),4)); end
%
% plot ke-spectrum (2:smallest coeff) and total ke (threshold (10:10:100) of PC1  
% median shows better results than mean
% variance between different case went down 
ithr=1; dokestd=0; for ib=1:4; subplot(4,1,ib); if ithr==0; idx_thr = ':'; else; idx_thr = nonzeros(parm.iCAMll_thres{ib}(:,ithr)); end; yyaxis left; contourf(squeeze(median(FFTke{ib}(2:10,1,:,idx_thr),4)),30,'linestyle','none'); colorbar; caxis([0 0.02]); yyaxis right; if (dokestd) plot(squeeze(std(mean(mean(var{ib}(:,:,:,idx_thr,idx.u).^2+var{ib}(:,:,:,idx_thr,idx.w).^2,1),2),0,4))); hold on; end; plot(squeeze(median(mean(mean(var{ib}(:,:,:,idx_thr,idx.u).^2+var{ib}(:,:,:,idx_thr,idx.w).^2,1),2),4))); hold off; disp(size(vort{ib}(:,:,:,idx_thr),4)); end
%
% plot ke only
ithr=0; for ib=1:4; subplot(4,1,ib); if ithr==0; idx_thr = ':'; else; idx_thr = nonzeros(parm.iCAMll_thres{ib}(:,ithr)); end; plot(squeeze(median(mean(mean(var{ib}(:,:,:,idx_thr,idx.u).^2+var{ib}(:,:,:,idx_thr,idx.w).^2,1),2),4))); hold on; plot(squeeze(median(mean(mean(var{ib}(:,:,:,idx_thr,idx.u).^2+var{ib}(:,:,:,idx_thr,idx.w).^2,1),2),4))); hold off; disp(size(vort{ib}(:,:,:,idx_thr),4)); end
%
doplot=1
if (doplot)
   nv = 8;
   for ib=1:4
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % circulation over 2D CRM slice %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subplot(4,1,ib);
      plot(squeeze(mean(abs(sum(sum((vort{ib}(:,:,:,:)),1),2)),4)));
      iCAMll{ib} = find(parm.iillt_zm_PC{ib}(3,:)==it);
      ii=1; % counter for saved PC1 indices
      for ill = iCAMll{ib}%1:2=lon,lat
         iiv=1; % counter for subplots 
         ilon = parm.iillt_zm_PC{ib}(1,ill); lon(ilon)
         ilat = parm.iillt_zm_PC{ib}(2,ill); lat(ilat)
   %   display(['lon: ' num2str(lon(ill(1))) ' lat: ' num2str(lat(ill(2))) ' PC1 mag: ' pc_spz{ib}(:,1)] )
   %   contourf(x,t,squeeze(mean(mean(FFT.ke{ib}(2:nx/2+1,izcoef,1:nlag+1,:),4),2))');
      %contourf(x,t,squeeze(mean(mean(FFT.ke{ib}(2:nx/2+1,izcoef,:,:),4),2))');
         for ixx = ixcoefs 
            fftke_std(ixx) = std(squeeze(FFTke{ib}(ixx,izcoef,1:nt,ii)));
         end
         for itlag = 1:nt
            fftke(ixcoefs,itlag) = squeeze((FFTke{ib}(ixcoefs,izcoef,itlag,ii)-mean(FFTke{ib}(ixcoefs,izcoef,1:nt,ii),3)))./fftke_std(ixcoefs)';
         end
         %for ixx = ixcoefs
         %   plot(t,fftke(ixx,:),'ko','markersize',(ixcoefs(end)-ixx+1)/2); hold on; % ylim([0 .7])
         %end; ylim([-3 3]); hold off
         subplot(nv,1,iiv);iiv=iiv+1;
         contourf(t,ixcoefs,fftke(ixcoefs,:),'linestyle','none'); % ylim([0 .7])
         set(gca,'YDir','Reverse')
   %      lh = legend(num2str(ixcoefs(:))); 
   %      set(lh,'location','northeastoutside');
   %      xlabel('time from large-PC1 event (hrs)')
         title(['KE normalized spectrum' 'k_x:' num2str(ixcoef) 'k_z' num2str(izcoef)])
         %
         subplot(nv,1,iiv);iiv=iiv+1;
         yyaxis left
         contourf(t,ixcoefs,squeeze(FFTke{ib}(ixcoefs,izcoef,:,ii)),'linestyle','none'); 
         set(gca,'YTick',ixcoefs)
         set(gca,'YTickLabel',round(xwavelength,1),'FontSize',6)
         %set(ax,'FontSize',6)  
         set(gca,'YDir','Reverse');
         yyaxis right
         plot(t,squeeze(FFTke{ib}(1,izcoef,:,ii)),'w')
         title(['KE spectrum'  'k_z' num2str(izcoef)])
         %
   %      subplot(nv,1,iiv); iiv=iiv+1;
   %      plot(t,squeeze(mean(mean(FFTke{ib}(ixcoef+1,izcoef,1:nt,ii),2),4))); % ylim([0 .7])
   %      title(['KE spec' 'k_x:' num2str(ixcoef+1) 'k_z' num2str(izcoef)])
   %   set(gca,'XTickLabel',t)
   %   xlabel('Wavelength (km)')
   %   set(gca,'XTick',x);
   %   set(gca,'XTickLabel',round(xwavelength,1))
         %
   %      subplot(nv,1,iiv);iiv=iiv+1;
      % spcam large-scale prec
   %      plot(t,squeeze(mean(mean(spprec{ib}(end,1:nt,ii),1),3))'-squeeze(prec(ilon,ilat,1:nt))); % surface spprec equals precc
   %      title('rain rate')
         %
         subplot(nv,1,iiv);iiv=iiv+1;
         yyaxis left
         contourf(t,lev,(squeeze(varCAM{idv.dt}(ilon,ilat,1:30,it-lagind-1))),'linestyle','none'); %caxis([-7 7]);
         set(gca,'YDir','Reverse')
         yyaxis right
         plot(t,squeeze(varCAM{idv.precc}(ilon,ilat,it-lagind-1)),'w','linewidth',2.5); %ylim([0 150]); % -1 is to compare .rh0. with spcam .r.
         title('CAM heating & rain rate')
         %
         subplot(nv,1,iiv);iiv=iiv+1;
         yyaxis left
         contourf(t,lev,(squeeze(varCAM{idv.spdt}(ilon,ilat,1:30,it-lagind-1))),'linestyle','none');% caxis([-7 7]); 
         set(gca,'YDir','Reverse')
         yyaxis right
         plot(t,squeeze(varCAM{idv.sppflx}(ilon,ilat,end,it-lagind-1)),'w','linewidth',2.5); %ylim([0 150]); % -1 is to compare .rh0. with spcam .r.
         title('SPCAM heating & rain rate')
         %
         subplot(nv,1,iiv); iiv=iiv+1;
         contourf(t,parm.zint/1000,squeeze(mean(vort{ib}(:,:,:,ii),1)),'linestyle','none'); ylim([0 10])
         %plot(t,abs(squeeze(mean(u_p{ib}(1:2,1:nt,ii),1))),'g'); xlim([-3.5 3.5]);% hold on   % 976:820
         title('SPCAM vort')
         %
         subplot(nv,1,iiv); iiv=iiv+1;
         contourf(t,parm.zint/1000,squeeze(mean(div{ib}(:,:,:,ii),1)),'linestyle','none'); hold on; ylim([0 10])
         quiver(t,parm.zint/1000,squeeze(mean(var{ib}(:,:,:,ii,idx.u),1)),squeeze(mean(var{ib}(:,:,:,ii,idx.w),1)),0.3,'w'); xlim([-3.5 3.5]); hold off
         title('SPCAM div')
         %
         subplot(nv,1,iiv); iiv=iiv+1;
         %contourf(1:32,parm.zint/1000,uw{ib}(:,:,1,nlag+1,ii)','linestyle','none')
         contourf(1:32,parm.zint/1000,u_tlag_eof{ib}(:,:,ii)','linestyle','none'); ylim([0 10])
         title('SPCAM u_eof1 field over nlag*2+1 time')
         %
         subplot(nv,1,iiv); iiv=iiv+1;
         plot(linspace(0,10,parm.nzi),mean(umap,1))
         title('SPCAM u_eof1 mean vert profile')
   %      subplot(nv,1,iiv); iiv=iiv+1;
   %      contourf(t,parm.zint/1000,squeeze(mean(uw{ib}(:,:,1,:,ii),1)),'linestyle','none');
   %      title('SPCAM u')
   %      plot(t,abs(squeeze(mean(dRhoU{ib}(1,1:nt,ii),1)))); xlim([-3.5 3.5]);    % 976:820
   %      contourf(t,crmlev(1:nz-1),(squeeze((dRhoU{ib}(1:nz-1,1:nt,ii)))),'linestyle','none'); xlim([-3.5 3.5]);    % 976:820
         %
   %      subplot(nv,1,iiv); iiv=iiv+1;
   %      plot(t,abs(squeeze(mean(dRhoW{ib}(1,1:nt,ii),1)))); xlim([-3.5 3.5]);   % 976:820
   %      contourf(t,crmlev(1:nz-1),(squeeze(-(dRhoW{ib}(1:nz-1,1:nt,ii)))),'linestyle','none'); xlim([-3.5 3.5]);    % 976:820
   %      plot(t,abs(squeeze(mean(u_p{ib}(9:18,1:nt,ii),1))),'r'); xlim([-3.5 3.5]); hold off; % 763:197
   %      subplot(nv,1,iiv); iiv=iiv+1;
   %      contour(t,crmlev(1:nz-1),abs(squeeze((u_p{ib}(1:nz-1,1:nt,ii))))); %ylim([0.0 0.02])
   %      title('mean vert wind shear u_p')
         %
   %      subplot(nv,1,iiv);iiv=iiv+1;
   %      contour(t,crmlev(2:nz),squeeze(mean(N2{ib}(1:nx,1:nz-1,1:nt,ii),1)));
   %      title('Buoyancy frequency')
         %
   %      subplot(nv,1,iiv);iiv=iiv+1;
   %      contour(t,crmlev(2:nz),(squeeze(squeeze(mean(N2{ib}(1:nx,1:nz-1,1:nt,ii),1))./u_p{ib}(1:nz-1,1:nt,ii).^2)),[0.1,0.15,0.25]);% caxis([0 .05])
   %      title('Richardson #') % bigger (smaller) means more (less) stable env less (more) shear highly (lowly) stratified
         %subplot(5,1,5)
         %contourf(t,crmlev,squeeze(mean(crm_qr(ilon,ilat,1:nx,1:28,1:nt),3)),'linestyle','none'); 
      % cam large-scale prec
         ii=ii+1;
         pause
      end
   end
end
