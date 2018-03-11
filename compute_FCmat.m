RtoZ=0;
folder =pwd;
%data = dir(fullfile(folder,'*.mat'));

%for c = 1:size(data,1); 
c=1;
     yrs= load('MR_1427581_3b35_seq.mat');
     [m,~] = corr(yrs.averaged_bold_seq);
     m(abs(m)<.01)=0; % set small values to zero
     m(isnan(m))=0; % remove any NaNs
     %if Abs, CM=abs(CM); end
     m=triu(m,1)+triu(m,1)';

     if RtoZ==1;
         m = atanh(m);
     end
     fMRI(:,:,c) = m;
     clear m
%end
