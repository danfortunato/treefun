function vals = coeffs2checkvals(coeffs,x,y)
%
%

persistent Evalx Evaly
[p,~,nd] = size(coeffs);
ncheckpts = numel(x); % hopefully not too large...

Evalx = ones(ncheckpts, p);
Evaly = ones(ncheckpts, p);
Evalx(:,2) = x(:);
Evaly(:,2) = y(:);
for k=3:p
  Evalx(:,k) = 2*x(:).*Evalx(:,k-1)-Evalx(:,k-2);
  Evaly(:,k) = 2*y(:).*Evaly(:,k-1)-Evaly(:,k-2);
end
% coeffs to value map
% vals = zeros(nd,ncheckpts);
% for k=1:ncheckpts
%   tmp1 = permute(tensorprod(Evalx(k,:),coeffs,2,1),[2 1 3]);
%   vals(:,k) = squeeze(permute(tensorprod(Evaly(k,:),tmp1,2,1),[2 1 3]));
% end
vals = zeros(nd,ncheckpts);
for k=1:ncheckpts
  tmp1 = permute(tensorprod(Evaly(k,:),coeffs,2,1),[2 1 3]);
  vals(:,k) = squeeze(permute(tensorprod(Evalx(k,:),tmp1,2,1),[2 1 3]));
end

end