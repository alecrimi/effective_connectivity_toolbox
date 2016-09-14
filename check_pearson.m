function m = check_pearson (bold_seq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function compute the functional connectivity given an fMRI series,
% and plot it to check visually if it is meaningful.
%
% Input:
%       bold_seq = matrix with each row an fMRI sequence for a brain area
% Output:
%       m = functional connectivity matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

RtoZ=0;

[m,~] = corr(bold_seq);
m(abs(m)<.01)=0; % set small values to zero
m(isnan(m))=0; % remove any NaNs

%if Abs, CM=abs(CM); end
m=triu(m,1)+triu(m,1)';

if RtoZ==1;
       m = atanh(m);
end

imagesc(m);