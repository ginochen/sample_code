function h = plot_waterphase_all(var,parm)
% plot basinwise spatially averaged organized cloud water over lead time
% plot TS,LFLX in rh0 to show that the difference of strength is due
% to SST
% plot EOF1 vert profile for each basin at the side
% 
% parm.lagind*30min
%
%load fname
%
run ~/scripts/matlab/startup
%
qthres=[0.01 0.01]; 
nlons=size(var{1},1); 
nt=size(var{1},4); 
ii=':'; % spatial indices for average
%
[X Z]=meshgrid([1:parm.nx],[1:parm.nzi]);
%X = [1:parm.nx];
%Z = [1:parm.nzi];
% vertical wind profile
startz = linspace(1,parm.nzi,20); 
startx = ones(size(startz))*1;
for ib=1:4;
   tic 
   figure; [ha ~] = tight_subplot(nt,nlons,[0.01,0.0],[.05 .05],[.05 .05]); % (Nh,Nw,gap,marg_h,marg_w)
   iii=1;
   for iit=1:nt;
      for iil = 1:nlons
         axes(ha(iii));
         contourf(X,Z,squeeze(var{ib}(iil,:,:,iit,4)-var{ib}(iil,:,:,iit,5))'*1000,qthres,'linestyle','none','facecolor',[.7 .7 .7]); hold on
         contourf(X,Z,squeeze(var{ib}(iil,:,:,iit,6))'*1000,qthres,'linestyle','none','facecolor',[.5 .5 .5]); 
         h=streamslice(X,Z,squeeze(var{ib}(iil,:,:,iit,1))',squeeze(var{ib}(iil,:,:,iit,2))');
         set(h,'color',[.6 .6 .6]);
         %streamline(stream2(X,Z,squeeze(mean(var{ib}(iil,:,:,iit,ii,1),5))',squeeze(mean(var{ib}(iil,:,:,iit,ii,2),5))',X,Z)); % last two arg are startx startz
%         u_zave = zave(squeeze(mean(mean(var{ib}(iil,:,:,iit,ii,1),2),5)),parm.nx,parm.nzi);
%         w_zave = sparse(parm.nx,parm.nzi);
%         quiver(X,Z,squeeze(mean(var{ib}(iil,:,:,iit,ii,1),5))',squeeze(mean(var{ib}(iil,:,:,iit,ii,2),5))',0.3); hold off
%          quiver(X,Z,u_zave',w_zave',3); hold off; % this is for MATLAB2015
         % quiver(X,Z,u_zave',w_zave',3,'HeadStyle','plain'); hold off; pause % this is for MATLAB2016b
         iii=iii+1;
   end; pause;end; 
   set(ha(1:end),'XTickLabel','')
   set(ha(setdiff([1:nlons*nt],[1:nlons:nlons*nt])),'YTickLabel',''); % remove all Ytick except the left side subplots
%   print(gcf,['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/18km_lons11_lead0-10-lead48/organizeQcQiQr_basin' num2str(ib) '_t' num2str(it)],'-djpeg99')
%   savefig(gcf,['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/18km_lons11_lead0-10-lead48/organizeQcQiQr_basin' num2str(ib) '_t' num2str(it) '.fig'])
   toc
end


function  var = zave(varin,nx,nz)
   ix = ceil(nx/2);
   var = sparse(nx,nz);
   var(ix,:) = varin;
