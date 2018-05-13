% merge all variable time and spatial index into one dimension
vname = {'div', 'vort', 'u', 'w', 'T', 'qT', 'qc','qr'};
archive = '/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/';
casename = '18km_lons11_lead0-10-lead48';%'15km_lead8-24'; %casename = '15km';
vind = [1,2,4,5,6]; % 'u' 'w' 'qT' 'qv' 'qr'
ii(1:4)=1;ii2(1:4)=0;
docamsclvar=0;
docrmsclvar=1;
tind = [11:10:491];
for it = tind
   load([ archive casename '/var_PC1_' num2str(it) '.mat'],'var','parm'); 
   for ib=1:4
      ii2(ib) = ii(ib)+size(var{ib},5)-1; 
      if (docamsclvar)
         iilon{ib} = parm.iillt_zm_PC{ib}(1,parm.iCAMll{ib});
         iilat{ib} = parm.iillt_zm_PC{ib}(2,parm.iCAMll{ib});
         for i = 1:numel(iilon{ib})
            dtcond{ib}(1:30,:,ii(ib)+i-1) = varCAM{idv.dt}(iilon{ib}(i),iilat{ib}(i),1:30,it-parm.lagind-1); precc{ib}(:,ii(ib)+i-1)  = varCAM{idv.precc}(iilon{ib}(i),iilat{ib}(i),it-parm.lagind-1);
            sppflx{ib}(:,ii(ib)+i-1) = varCAM{idv.sppflx}(iilon{ib}(i),iilat{ib}(i),end,it-parm.lagind-1);   spdt{ib}(1:30,:,ii(ib)+i-1) = varCAM{idv.spdt}(iilon{ib}(i),iilat{ib}(i),1:30,it-parm.lagind-1);
         end
      end
      if (docrmsclvar)
         for iv = vind 
            evalc(sprintf('%s{ib}(:,1:parm.nx,1:parm.nzi,1:numel(parm.lagind),ii(ib):ii2(ib)) = var{ib}(:,:,:,:,:,iv)',parm.var{iv}));
         end
      end
      ii(ib) = ii2(ib)+1;
   end
end

% sum over all basin space-time indices
for iv = vind 
   for ib=1:4
      eval(sprintf('%s_sum{ib} = sum(%s{ib},5)', parm.var{iv}, parm.var{iv}));
      nSum(ib) = size(u{ib},5);
   end
   eval( sprintf(' save(''%s/%s/%s_sum_PC1.mat'',''%s_sum'',''nSum'',''-v7.3'') ',archive, casename, parm.var{iv}, parm.var{iv}) )
end
% save all basin space-time indices for a variable, FILES are HUGE
for iv = vind  
   eval( sprintf(' save(''%s/%s/%s_all_PC1.mat'',''%s'',''-v7.3'') ',archive, casename, parm.var{iv}, parm.var{iv}) )
end
