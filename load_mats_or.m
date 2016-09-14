%function [ M, matrices] = load_mats_or

matrices = zeros(90,90,20);

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/1427581/experiments/MR_1427581_3b35/MR_1427581_3b35/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,1) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/1523112/experiments/MR_1523112_43e6/MR_1523112_43e6/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,2) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/1616468/experiments/MR_1616468_54cf/MR_1616468_54cf/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,3) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/2021454/experiments/MR_2021454_4ade/MR_2021454_4ade/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,4) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/2328270/experiments/MR_2328270_6ebf/MR_2328270_6ebf/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,5) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/2505567/experiments/MR_2505567_3f9e/MR_2505567_3f9e/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,6) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/2541244/experiments/MR_2541244_5224/MR_2541244_5224/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,7) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/2784584/experiments/MR_2784584_1a4d/MR_2784584_1a4d/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,8) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/2861923/experiments/MR_2861923_4f24/MR_2861923_4f24/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,9) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/3374719/experiments/MR_3374719_104f/MR_3374719_104f/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,10) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/3466763/experiments/MR_3466763_386d/MR_3466763_386d/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,11) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/3611212/experiments/MR_3611212_36f0/MR_3611212_36f0/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,12) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/4097119/experiments/MR_4097119_886b/MR_4097119_886b/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,13) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/4143704/experiments/MR_4143704_4f12/MR_4143704_4f12/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,14) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/4290056/experiments/MR_4290056_39b5/MR_4290056_39b5/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,15) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/4323037/experiments/MR_4323037_4dd1/MR_4323037_4dd1/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,16) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/6539040/experiments/MR_6539040_7cc6/MR_6539040_7cc6/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,17) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/6692978/experiments/MR_6692978_5233/MR_6692978_5233/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,18) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/6913939/experiments/MR_6913939_6047/MR_6913939_6047/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,19) = M;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/9421819/experiments/MR_9421819_1328/MR_9421819_1328/SCANS/9/NIFTI/';
M = csvread([connectome_folder 'foo.csv']);  
matrices(:,:,20) = M;

save('matrices.mat');

%{
%M = mean(matrices,3);

matrices_bin = matrices > 0;
%tot = sum(matrices_bin,3);
%M = tot.*(tot>9);
M = sum(matrices_bin,3);

M = M > 9;
%}
 

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
