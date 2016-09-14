function [M, etem] = structured_G_causality(name_session, name_fold, norm_opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% These scripts aim at computing the effective connectivity by using functional and structural MRI data. 
% The main function to compute the autoregressive model as in A.Crimi et al. "Effective Brain Connectivity Through a Constrained Autoregressive Model" MICCAI 2016
% 
% The path is assumed to be of the data from the NKI 1000Connectome project 
%
% Input: 
% name_session: the folder with the subjects for the same session (or folder with all subjects)
% name_fold: main folder for single subject e.g.'MR_9421819_1328';
% norm_opt: = 0 means no normalization, = 1 means binarize, = 2 means divide by the maximum
%
% Output:
% M: the resulting effective matrix
% etem: the reconstruction error over the iterations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of brain areas from the atlas
brain_par = 90; % 90 because we are used the main part 
% Max number of iterations for the gradient descent
maxiter = 5000;  
% Initial error of reconstruction
error_des = 1e10; %Set it initially very high
% Acceptable error threshold for the reconstruction error
th = 0.016; 
% Learning rate of the gradient descent
eta = 0.00005;%0.0000000005; %0.000005; %nozeroing 0.00005; 0.005;%
% Initial frames of BOLD sequence to be skiped 
skip_bold = 5;
% Mimimum number of fibers to be kept
min_tract_num = 2;

% Add path for libraries, e.g. loading NIFTI files
addpath('libs');


%%%%%%%%%%%%%% FUNCTIONAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take part of the name, it is needed for subfolders
pre = name_fold(4:10);
folder_bold = [name_session '/' pre '/experiments/' name_fold '/' name_fold '/SCANS/4/NIFTI/' ];
bold_file = load_untouch_nii([folder_bold 'BOLDrestingCAP_pp_final.nii']);   
bold_sequence = bold_file.img; 
[r, c, d, t] = size(bold_sequence); 
averaged_bold_seq =  zeros(t,brain_par);

% Load atlas registered to the BOLD file
atlas_file = load_untouch_nii([folder_bold 'atlas_reg.nii.gz']);
atlas = atlas_file.img;

% Averaging Bold signal
% For all frames
for frame = 1 :  t
    
    bold_frame = bold_sequence(:,:,:,frame);
    % For all areas
    for b = 1 : brain_par 
         temp = zeros(2,1);
         % For all axial plane
         for all_d = 1 : d
             %Find indices          
             [indices_x,indices_y] = find( atlas(:,:,all_d) == b ); 
            
             for kk = 1: length(indices_x)
                    temp(end+1) =  bold_frame(indices_x(kk),indices_y(kk),all_d);
             end
         end
         % Remove empty voxels
         temp(temp==0) = [];      
         % Save indices voxels with mean for the specific brain area
         averaged_bold_seq(frame,b) =   mean( temp );
    end
end

% Remove mean from the already averaged bold sequences
mean_sig = mean( averaged_bold_seq(:,:),2);
for b = 1 : brain_par
    averaged_bold_seq(:,b) =  averaged_bold_seq(:,b) - mean_sig ; 
end


%%%%%%%%%%%%%%%% STRUCTURAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load structural matrix
connectome_folder = [name_session '/' pre '/experiments/' name_fold '/' name_fold '/SCANS/9/NIFTI/' ];
M =  csvread([connectome_folder 'foo.csv']);   %load_mats_or; %

% Filter out all connection with less than n fiber reaching
M = M.*( M > ( min_tract_num - 1 ) );

% Either normalize or binarize the matrix
switch norm_opt
    case 0 
        disp('No normalization of the initial matrix');
    case 1
        disp('Initial matrix binarized');
        M = M >0;
    case 2
        disp('Initial matrix normalized according to its largest value');
        max_val = max(max(M));
        M = M/max_val;
    otherwise
        disp('No normalization of the initial matrix')

end

% Construct a matrix with only zeros and ones to be used to reinforce the zero connection 
noconn_M = ~~M; 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradient descent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 0; % Counter of the current iteration
etem = zeros(1,1); % Reconstruction error along iterations

%Until criteria are met
while  (  error_des > th  )  && ( iter < maxiter )
        error_des = 0; 
        gradient = zeros(brain_par,brain_par);
        
        %Error to be used as a criterion
        for jj = skip_bold : t -1
             error_des  = error_des +   norm(  M * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
        end
        error_des = 0.5 * error_des;
        
        etem(end+1) =  error_des; % TODO: THIS IS NOT CLEAN, IT HAS TO BE DONE BETTER!
            
        % Gradient computation 
        for jj = skip_bold : t -1
            gradient = gradient + ( ( M*averaged_bold_seq(jj,:)'*averaged_bold_seq(jj,:))  - (averaged_bold_seq(jj+1,:)'*averaged_bold_seq(jj,:)) ); 
            % gradient = gradient + ( (averaged_bold_seq(jj,:)' - M*averaged_bold_seq(jj+1,:)')   * averaged_bold_seq(jj,:)  )       ;
             %gradient = gradient + ( (averaged_bold_seq(jj+1,:)' - M*averaged_bold_seq(jj,:)')   * averaged_bold_seq(jj,:)  )       ;
             % gradient = ( ( M*averaged_bold_seq(jj,:)' - averaged_bold_seq(jj+1,:)')   * averaged_bold_seq(jj,:)  )       ;
        end  
        %   imagesc(gradient); pause % 
        %gradient
        M =  M - eta * gradient;
        
        %Reinforce where there was no connection at the beginning
        M = M .* noconn_M;
        %Remove negative values
        M = M .*( M>0);
        iter = iter + 1 ;
end 

% Remove negative values at the end of the process
%M = M .*( M>0);
etem(etem==0) = [];    
%figure; imagesc(M)
%figure; plot(etem)
save([connectome_folder  'eff_conn0']);
save([ name_fold  '_seq'],'averaged_bold_seq'); 