function plot_mcs_time_hist()
% plot of the histogram of MCS time-of-day conditioned on land/ocean,
% size, and duration.  This is to see if there are preference for diurnal cycle to peak or bottom
% So, x-axis is the time-of-day, y-axis is the total # of MCSs 
  
  load mcs_clusters
  time = [];
  for ilnd = 0:1
    ii = 1;
    mcsci = find([mcsnllave'>=1 & island==ilnd & nt4Cl'>=1])';
%    mcsci = find([mcsnllave'>=1 & island==ilnd])'; % mcs cluster index set
    nmcss = numel(mcsci)
    for ic = 1:nmcss 
%      for iit = 1:nt4Cl(mcsci(ic))
%        time{ilnd+1}{ic}{iit} = str2num(t{mcsillt4Cl{mcsci(ic)}{iit}(3)}(end-4:end))/3600
%      [~,iit] = max(mcsnll(ic,:));
iit=1;
        time{ilnd+1}(ii) = str2num(t{mcsillt4Cl{mcsci(ic)}{iit}(3)}(end-4:end))/3600;
        ii=ii+1;
%      end
    end
    tvec = 0:3:21
    for it = 1:numel(tvec) 
      nt(it) = sum(time{ilnd+1}==tvec(it));
    end
    plot(tvec,nt/sum(nt),'o-'); hold on; 
pause
  end

% this is to see if the 
  load mcs_clusters
  time1 = [];
  time2 = [];
  for ilnd = 0:1
    mcsci = find([mcsnllave'>=1 & island==ilnd & nt4Cl'>=1])';
    nmcss = numel(mcsci);
    for ic = 1:nmcss 
%      for iit = 1:nt4Cl(mcsci(ic))
        time1{ilnd+1}(ic) = str2num(t{mcsillt4Cl{mcsci(ic)}{1}(3)}(end-4:end))/3600;
        time2{ilnd+1}(ic) = str2num(t{mcsillt4Cl{mcsci(ic)}{end}(3)}(end-4:end))/3600;
%      end
    end
  end
%{
    figure
    binranges = -3:3:24;
    bincounts = histc(time{ilnd+1},binranges);
    bar(binranges,bincounts,'histc')
%}
  end
  %%%%%%%%%%% plot and see if large clusters (>6 cells) are over land in the afternoon and over ocean in the early morning %%%%
  hist(time(mcsnllave>=6 & island==0),30)
  island % use island to find the ic-cluster that are land to make the hist of time
