function [svec, variance, pc, pc_mean, varmean]=svdVar(varin,doecon,dosignflip,dopca,ignoreNaNrandvar,discardNaNsamp);
%function [svec, variance, pc, pc_mean, varmean]=svdVar(varin,doecon,dosignflip,dopca,ignoreNaNrandvar,discardNaNsamp);
% varin(nsample,nrandomvar)
% variance: an array of variance of each mode 
% svec: svec = V,  a matrix of singular vector with NaN (column vectors), saves
% only the top ilim vectors that has 90% of total variance
% pc: pc=U*S,  a matrix of principle component time series column vectors.
%    Average of pc(:,im) will be near zero since its the anomaly magnitude
% pc_mean: X_mean projected onto svec, used for plotting CAFD purpose
%    X_mean = sum(pc_mean(im)*svec(:,im))
% Why center at zero fro PCA? (From Wikipedia https://en.wikipedia.org/wiki/Principal_component_analysis): 
% Mean subtraction (a.k.a. "mean centering") is necessary for performing PCA to
% ensure that the first principal component describes the direction of maximum
% variance. If mean subtraction is not performed, the first principal component
% might instead correspond more or less to the mean of the data. A mean of zero
% is needed for finding a basis that minimizes the mean square error of the
% approximation of the data.
% 
% For PCA:
% Suppose X(z,t) is the time anomaly, then time mean value of pc for any mode must be zero. 
% Proof: 
% X(z,t) = \sum(pc_k(t)*svec_k(z)) => ( project X onto svec_k using orthonormality of svec_k )
% pc_k(t) = \int(X(z,t) svec_k(z))dz => ( take time mean operator 1/T*\int[ ]dt )
% 1/T * \int[pc_k(t)]dt = 1/T * \int[\int(X(z,t) * svec_k(z))dz]dt => ( only X is a function of t )
%                       = \int(\int[X(z,t)]dt * svec_k(z))dz = 0 ( since X is centered at zero in time ) 
% 
% Since pc_k(t) = \int(X(z,t) * svec_k(z))dz  = \int((X_original(z,t) - X_mean(z))*svec_k(z))dz =>
% pc_k(t) = \int((X_original(z,t)*svec_k(z))dz - \int(X_mean(z)*svec_k(z))dz =>
% pc_k(t) = pc_k_original(t) - pc_k_mean =>
% X_original(z,t) = \sum(pc_k(t) + pc_k_mean)*svec_k(z) = 
% Therefore we can plot the decomposed X_original, for instance, mode-1
% timeseries (pc_1(t)+pc_1_mean)*svec_k(z) centered around pc_1_mean*svec_k(z)if
%%%%%% set default arg %%%%%%%%%%%%%%%%%%%
disp(['matrix(randvar,samp) size = ' num2str(size(varin))]);
ids = 1;
idv = 2; 
if size(varin,ids)<size(varin,idv)
  display('warning, samples are less than the r.v. dimension, may want to skip pca for this data...')
end
if ~exist('dopca')
  dopca = 1; % set default to do PCA
end
if ~exist('doecon')
  doecon = 0; % set default to no economical svd, econ will truncate the modes to the shortest dimension, so if nsample is shorter it'll just do the modes limited to nsamp instead of nrandomvar
              % The econ modes will be the same as non econ modes, saving the bigger modes and truncated to lesser modes
end
disp(['doecon = ' num2str(doecon)])
if ~exist('ignoreNaNrandvar')
  ignoreNaNrandvar = 0; % default to 0, suppose one random variable is always NaN, then find it out and don't do svd on that dimension, but put the NaNs back in the end
end
if ~exist('discardNaNsamp')
  discardNaNsamp = 1;
end
if ~exist('dosignflip')
  dosignflip=0;
end
disp(['dosignflip = ' num2str(dosignflip)])

if (ignoreNaNrandvar)
  display('Discard the NaN variables dimension, put the NaNs back in the end')
  innan= find(~isnan(varin(1,:))); % select a sample, and check which randvars are nan
  inan = find(isnan(varin(1,:)));
  varin=varin(:,innan); 
else
  display('Use all variable dimensions, specify the non-NaN variables indices before calling svdVar()')
  innan = 1:size(varin,idv);
  inan = [];
end
if (discardNaNsamp)
  inansamp=find(isnan(sum(varin,idv)));
  innansamp=find(~isnan(sum(varin,idv)));
  display(['Discard sample indices ' num2str(inansamp') ', which contain at least a NaN in one variable dimension'])
  varin = varin(innansamp,:); % find the samples that doesn't contain NaN in any of the variable dimension, 
     % this may leed to all samples being discarded if one dimension is always NaN, use carefully
end
if (dopca)
  varmean = mean(varin,ids); % X_mean
  varin = bsxfun(@minus,varin,varmean);
else
  display('If doing PCA, make sure to set dopca to 1 to ensure sample mean to be subtracted from data and watch out for the eof vector sign mismatch')
end
%
disp('start svd')
if (doecon)
  [U S V] = svd(varin,'econ'); % A = U*S*V', singular vectors of U & V are the columns, and U (V) columns is a time (spacial) series 
  % varin(nsample,nrandvar)
else
  [U S V] = svd(varin); % A = U*S*V', singular vectors of U & V are the columns, and U (V) columns is a time (spacial) series 
end
disp('done svd')
if size(S,2)>1
  variance=diag(S).^2; % eigenvalue of the covariance matrix, same as third output arg (latent) of pca()
else
  variance = 0; % one sample doesn't have variance
end
%fracvar=eigencov/sum(eigencov); % fractional variance
%
pc = U*S; % same as the second output argument (score) of pca()

%{
if dosignflip
  disp('start sign flip')
% flip signs
  [sgns,loads] = sign_flip({pc,V},varin);
  pc = loads{1};
  svec = loads{2};
  disp('end sign flip')
end
%}
svec(innan,:)= V;
svec(inan,:) = NaN; % insert the NaNs back 
% calc pc_mean (X_mean projected onto svec(:,im))
vmean(1,innan) = varmean;
vmean(1,inan) = NaN;
inan
for im = 1:size(svec,2) % modes
  pc_mean(im) = vmean*svec(:,im); % amplitude of mean projected onto eof modes
end
% xsamp(1:nz,im,it) = (pc_mean(im) + pc(it,im))*svec(1:nz,im) is the vertical profile 
% sampled around xmean(1:nz,im) = pc_mean(im)*svec(1:nz,im)
%






% var = U*S*V' = pc*svec'; 
% V'V = I; % demonstration of unitary
% pc(:,:) = var*V = U*S; % the PC time-series vector pc(itime,imode)
% variance on the i-th mode V(:,i) is the S(i,i)
% TODO:
% plot the ith time-series vector of pc(:) associated to each V(:,i)
%
%
%contourf(reshape(Vtmp,dim(1),dim(2)),100,'linestyle','none')

