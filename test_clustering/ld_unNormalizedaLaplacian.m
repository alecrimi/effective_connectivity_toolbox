function [Lap,u,e] = ld_unNormalizedaLaplacian(A)

A = A.*(A>0);   A(find(isnan(A)))=0;   A(find(isinf(A)))=0;
degs = sum(A,2)+eps;   % eps to avoid degree =0
D = diag(degs);
Lap = D-A;
[u,e] = eig(Lap);

if ~issorted(diag(e))
    [e,index] = sort(diag(e));
    e = diag(e);
    u = u(:,index);
    
end

end





  