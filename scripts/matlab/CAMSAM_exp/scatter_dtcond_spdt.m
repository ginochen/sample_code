idv.dt=1;    varname{idv.dt} = 'DTCOND';   var{idv.dt} = loadvar(varname{idv.dt},camCase,spd,1);% (K/day) <---- only moist parametrized heating
idv.spdt=30; varname{idv.spdt} = 'SPDT';   var{idv.spdt} = loadvar(varname{idv.spdt},spcamCase,spd,0);% T te
%idv.qrs=6;   varname{idv.qrs}  = 'QRS';    var{idv.qrs} = loadvar(varname{idv.qrs},spcamCase,spd,0); % (mm/day)
%idv.qrl=7;   varname{idv.qrl}  = 'QRL';    var{idv.qrl} = loadvar(varname{idv.qrl},spcamCase,spd,0); %(mm/day) <--- is this SPQRL?
itimes = [1:size(var{idv.spdt},4)];
for it=1:numel(itimes)
   dtcond_pm(:,:,it) = pwgtave(var{idv.dt},iilev,dlev);
   spdt_pm(:,:,it)   = pwgtave(var{idv.spdt},iilev,dlev);
%   qrs(:,:,it)    = pwgtave(var{idv.qrs},iilev,dlev);
%   qrl(:,:,it)    = pwgtave(var{idv.qrl},iilev,dlev);
end
dtcond_tpm= mean(dtcond_pm,3);
spdt_tpm  = mean(spdt_pm,3);
%qrsm   = mean(qrs,3);
%qrlm   = mean(qrl,3);
for ib=1:numel(basin)
   subplot(2,2,1)
   hist(spdt_tpm(iiocn{ib}))
   xlim([-5 35]);
   subplot(2,2,4)
   hist(dtcond_tpm(iiocn{ib}))
   xlim([-5 35]);
   subplot(2,2,2)
   scatter(dtcond_tpm(iiocn{ib}),spdt_tpm(iiocn{ib}));
   xlim([-5 35]);
   ylim([-5 35]);
   pause
end
save('dt_spdt_p_t_ave.mat','dtcond_pm','spdt_pm','dtcond_tpm','spdt_tpm','iilev','itimes')
