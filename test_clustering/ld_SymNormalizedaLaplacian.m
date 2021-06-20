function [Lsym,u,e] = ld_SymNormalizedaLaplacian(A)

A = A.*(A>0);   A(find(isnan(A)))=0;   A(find(isinf(A)))=0;
degs = sum(A,2)+eps;   % eps to avoid degree =0
L = diag(degs)-A; % L = D-A;

invdegs = (degs).^(-0.5);
Dn = diag(invdegs);
Lsym = Dn * L * Dn;    %L = D^(-1/2)*(D-W)*D^(-1/2);
[u,e] = eig(Lsym);
u= bsxfun(@rdivide, u, sqrt(sum(u.^2, 2)));

if ~issorted(diag(e))
    [e,index] = sort(diag(e));
    e = diag(e);
    u = u(:,index);
    
end


end                 
                    
  