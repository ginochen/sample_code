function [s n] =  mergeboundary(s,n,nx)
  % Purpose: used with continuousSet() to merge the two subsets containing boundary points for periodicity
  %%%%%% merge the two boundary systems into one due to the periodic domain %%%%%%%%%%%%
  if ( numel(s)>1 & any(s{1}==1) & any(s{end}==nx) ) % if 1 and 32 is in the set, then concat those two into one set
    s{1} = ([s{end},s{1}]); 
    s = s(1:end-1); % delete the last cell
    n(1) = n(1)+n(end);
    n = n(1:end-1);
  end
