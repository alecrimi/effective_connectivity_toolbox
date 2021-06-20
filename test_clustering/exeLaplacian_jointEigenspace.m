%%
caseLap = 'Normalized';
%caseLap = 'unNormalized'
% W set of connectivity matrices
[L,U,E] = ld_ComputeLaplacianGraph(W,'caseLap');
%% JointDiag 
Q = [];
for i = 1:size(L,3)
Q = [Q L(:,:,i)];
end
thr = 10^-8;%% We used 10^-8. Higher thr makes the algorithm faster 
% but the joint eigenspace is less close to the original eispaces of laplacians
[ V ]= jointD(Q,thr); % V is the Joint Eigenspace between two or more laplacians

%%  Reorder V for spectral clustering
[Vsort,Lambda_tilde,lambdaSort,lambdaID]= ld_reorderJointEigenspace(V,L);
plot(lambdaSort,'.b')

%%
K = 3; % Choose the number of cluster (through Spectral Gap of the average
% approximate eigenvalues)
[clusterID,C] = kmeans(Vsort(:,1:K),K, 'Maxiter', 200, 'Replicates', 500 ,'display','iter');
    