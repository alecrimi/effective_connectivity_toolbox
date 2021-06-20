function [Vsort,Lambda_tilde,lambdaSort,lambdaID]= ld_reorderJointEigenspace(V,L)

% Output:
% Vsort: joint Eigenspace Sorted
% Lambda: Apporiximate Eigenvalues



Lambda_tilde = zeros(size(L));
for i = 1:size(L,3)
    l = squeeze(L(:,:,i));
     lambda_tilde = diag(V'*l*V); % Approximate Eigenvalues l
     Lambda_tilde(:,:,i) = diag(lambda_tilde);
end

% Compute Average Approximate eigenvalues
lambdaAvg = mean(Lambda_tilde,3);
% Sort Eigenvalues ( from smallest to the biggest one)
[lambdaSort,lambdaID] = sort(diag(lambdaAvg),'ascend');
% Sort Joint Eigenspace
Vsort = V(:,lambdaID);
end