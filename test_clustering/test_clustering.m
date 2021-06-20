% Luca Dodero's code IIT 

load('matrices.mat');
 
M_odd = matrices(1:2:end,:,:);
M_even = matrices(2:2:end,:,:);
Mtemp = cat(1,M_odd,M_even);
M_odd = Mtemp(:,1:2:end,:);
M_even = Mtemp(:,2:2:end,:);
matrices = cat(2,M_odd,M_even);
 
caseLap = 'Normalized';
%caseLap = 'unNormalized'
% W set of connectivity matrices
[L,U,E] = ld_ComputeLaplacianGraph(matrices,caseLap);
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
%plot(lambdaSort,'.b')

%%
K = 4; % Choose the number of cluster (through Spectral Gap of the average
% approximate eigenvalues)
[clusterID,C] = kmeans(Vsort(:,1:K),K, 'Maxiter', 200, 'Replicates', 500 ,'display','iter');

%Load atlas
folder_atlas = '';%'/home/alex/Desktop/ICBM/MNI_0579/BOLD_MOSAIC_64_resting_st./2006-09-14_14_11_40.0/S59582/';
atlas_file = load_untouch_nii([folder_atlas 'atlas_reg.nii.gz']);
atlas = atlas_file.img;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save the clustered brain
clustered_brain = zeros(size(atlas));
brain_par = 90;
[r, c, d, t] = size(clustered_brain); 

%For each area of the atlas
for b = 1 : brain_par 
  
        %Find indices
         for all_d = 1 : d
             [indices_x,indices_y] = find( atlas(:,:,all_d) == b ); 
             
             %Read index and set on a copy the value
             for kk = 1: length(indices_x)
                      clustered_brain( indices_x(kk),indices_y(kk),all_d)= clusterID(b);
             end
         end
end

%Save clustered altas
atlas_file.img = clustered_brain;
save_untouch_nii(atlas_file, 'clustered_original.nii');
