function [R,P,ivvec] = corrcoef_3d(var,ivarsall,ilat,itimes)
itt=1;
for it = itimes
   ivv = 1;
   for iv = ivarsall
      if (ndims(var{iv})==4)
        v1(:,ivv,itt) = mv3(var{iv}(:,ilat,:,it));
        ivvec(ivv)=iv; % output the variables that are eval
        ivv=ivv+1;
      end
   end
   [R(:,:,itt) P(:,:,itt)] = corrcoef(v1(:,:,itt));
   itt=itt+1;
end
