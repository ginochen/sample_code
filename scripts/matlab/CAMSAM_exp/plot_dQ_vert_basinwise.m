%idv.e=31; varname{idv.e} = 'E';        var{idv.e} = loadvar(varname{idv.e},camCase,1,0);% nt=754, b/c var{30} - var{1}(:,:,:,1:end-1); 
%idv.ea=32; varname{idv.ea} = 'EABS';     var{idv.ea} = abs(var{idv.e});
%idv.dt=1;  varname{idv.dt} = 'DTCOND';   var{idv.dt} = loadvar(varname{idv.dt},camCase,spd,1);% (K/day) <---- only moist parametrized heating
%idv.spdt=30; varname{idv.spdt} = 'SPDT';     var{idv.spdt} = loadvar(varname{idv.spdt},spcamCase,spd,0);% T tendency due to CRM (without QRL+QRS radiative heating)
%idv.pzm=10;      varname{idv.pzm} = 'PRECCDZM';      var{idv.pzm}     = loadvar(varname{idv.pzm},camCase,spd*1000,1);% (mm/day) 
%idv.psh=11;      varname{idv.psh} = 'PRECSH';        var{idv.psh}     = loadvar(varname{idv.psh},camCase,spd*1000,1);% (mm/day) 

basin{1}='IND';  % indian ocean
cutoff_lat{1} = [[-20, 0],   [0,   20] ]; % lon=[40 to 120], lat=[-20 to 0], etc
cutoff_lon{1} = [[40,  120], [40,  100]];
basin{2}='WPAC'; % western pacific
cutoff_lat{2} = [[-20, 0],   [0,   20] ];
cutoff_lon{2} = [[120, 180], [100, 180]];
basin{3}='EPAC'; % eastern pacific
cutoff_lat{3} = [[-20, 10],  [10,  20] ];
cutoff_lon{3} = [[180, 300], [180, 271]];
basin{4}='ATL';  % atlantic
cutoff_lat{4} = [[-20, 10],  [-20, 10],  [10,  20]];
cutoff_lon{4} = [[  0, 20],  [290, 360], [271, 350]];  
basin{5}='ALL';
%
landfrac = loadvar('LANDFRAC',camCase,1,0);
%load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/index/basin_index.mat','iiocn');
vartmp = var{idv.pzm};

for ib = 1:numel(basin)-1
   iii=1; % counter for a basin-wise index (lon lat time)
   for i = 1:numel(cutoff_lat{ib})/2 % a basin may be separated into two parts or more
      iilon = find(lon>=cutoff_lon{ib}((i-1)*2+1) & lon<=cutoff_lon{ib}((i-1)*2+2))';
      iilat = find(lat>=cutoff_lat{ib}((i-1)*2+1) & lat<=cutoff_lat{ib}((i-1)*2+2))';
      for it=1:numel(itimes)
         for ii = iilon
            for jj = iilat
               if (landfrac(ii,jj) == 0 & vartmp(ii,jj,it)>0) % if ocean index, if precip is deep
                  etmp{ib}(:,iii) = var{idv.ea}(ii,jj,:,it);
                  iillt{ib}(:,iii) = [ii,jj,it]; % index of lon lat time when the condition is satisfied
                  iii=iii+1;
               end
            end
         end
      end
   end
end
%
eabs_p1p30_basin=etmp;
save('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/eabs_vert/cond/eabs_p1p30_basin_preccdzm_gt0.mat',...
      'eabs_p1p30_basin','cutoff_lat','cutoff_lon','basin','lat','lon','lev','itimes','iillt')

