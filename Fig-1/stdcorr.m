function y = nancorr(x,dim,n)

% N is number of degrees of freedom
% Only n-1 DOF for N variables

chival = sqrt(2./(n-1)).*gamma(n/2)./gamma((n-1)/2);

y = bsxfun(@rdivide,squeeze(nanstd(x,[],dim)),chival);
y(abs(imag(y)) > 0) = nan; 
y = real(y); 




