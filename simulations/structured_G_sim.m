function [M, etem] = structured_G_sim(init, ts)
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
norm_opt = 0;

%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of brain areas from the atlas
[t, brain_par] = size(ts);
% Max number of iterations for the gradient descent
maxiter = 1000;  
% Initial error of reconstruction
error_des = 1e10; %Set it initially very high
% Acceptable error threshold for the reconstruction error
th = 0.016; 
% Learning rate of the gradient descent
eta = 0.0000000005;%0.00005;%0.0000000005; %0.000005; %nozeroing 0.00005; 0.005;%
% Initial frames of BOLD sequence to be skiped 
skip_bold = 5;
% Mimimum number of fibers to be kept
min_tract_num = 2;

% Add path for libraries, e.g. loading NIFTI files
addpath('libs');


%%%%%%%%%%%%%% FUNCTIONAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
averaged_bold_seq = ts;

% Remove mean from the already averaged bold sequences
mean_sig = mean( averaged_bold_seq(:,:),2);
size(mean_sig)
for b = 1 : brain_par
    averaged_bold_seq(:,b) =  averaged_bold_seq(:,b) - mean_sig ; 
end


%%%%%%%%%%%%%%%% STRUCTURAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load structural matrix
M = init;

% Filter out all connection with less than n fiber reaching
%M = M.*( M > ( min_tract_num - 1 ) );

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
            % gradient = gradient + ( (averaged_bold_seq(jj+1,:)' - M*averaged_bold_seq(jj,:)')   * averaged_bold_seq(jj,:)  )       ;
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
 
