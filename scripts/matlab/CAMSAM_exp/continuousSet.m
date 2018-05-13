function [s n] = continuousSet(S,N,doboundary,nx)
%function [s n] = continuousSet(S,N,doboundary,nx)
% Purpose: group a set of number into neighboring number subsets and regroup the boundary subsets together
% Condition: total set is 1:nx, S could be just a subset of this set (e.g., (30,31,32,1,2,3))
% S(:): a set of numbers that are ordered and needs to be binned into neighboring subset of numbers
% s{:}: subsets neighboring numbers of S 
% N: numel(S)
% n(:): numel(s{j})
% nx: 32 for crm
  i = 1; j= 1;
  if doboundary
    dif = [1 nx-1]; % difference between two numbers in S, neighbors if dif = 1, 
                    % if periodic boundary is considered neighbors then the difference is nx-1
  else
    dif = 1; 
  end
  while (i <= N)
    i1 = i;
    if (i~=N) % if the last index hasn't reached
      while (ismember(abs(S(i+1)-S(i)),dif)) % 1 is the difference for neighbor points, if periodic, then 1-nx is also neighbor 
        i=i+1; 
        if (i==N) 
          break; 
        end % if the last index is reached
      end
    end
    s{j} = S(i1:i); 
    n(j) = numel(s{j}); % numel of each subset
    j=j+1;
    i=i+1; 
  end

