function [L,U,E] = ld_ComputeLaplacianGraph(W,caseLap)

L = zeros(size(W));
U = zeros(size(W)); 
E = zeros(size(W));

for i = 1:size(W,3);
    subject = squeeze(W(:,:,i));
    caseLap
    switch(caseLap);
        case 'Normalized'
            [l,u,e] = ld_SymNormalizedaLaplacian(subject);
        case 'UnNormalized'
            [l,u,e] = ld_unNormalizedaLaplacian(subject);
    end
    
    L(:,:,i) = l;
    U(:,:,i) = u;
    E(:,:,i) =e;
end
    

end