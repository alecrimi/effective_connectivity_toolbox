function [granger, granger_F, granger_pvalue, granger_num, granger_den, granger_num_dof, granger_den_dof, granger_dir, granger_inst]=etc_granger(v,ar_order,varargin)
% etc_granger   use AR model for granger causality test on timeseries
%
% [granger, granger_F, granger_pvalue]=etc_granger(v,ar_order,[option, option_value,...]);
%   v: [n,m] timeseries of n-timepoint and m-nodes
%   n-timepoint and k-observations
%
%   ar_order: the order of the AR model
%   option:
%       'ts_name': names for time series, in string cells
%       'g_threshold': the threshold of the granger causlity to be
%       considered significant.
%       'g_threshold_p': the p-value threshold of the granger causlity to be
%       considered significant.
%   granger: [m,m] granger causality matrix. each entry is the variance ratio: (res_orig^2/res_additional_timeseries^2)
%
% fhlin@aug 8 2005
%

flag_display=1;
ts_name={};
g_threshold=[];
g_threshold_p=[];
flag_pairwise=1;
flag_robust=0;

for i=1:length(varargin)/2
    option=varargin{i*2-1};
    option_value=varargin{i*2};

    switch option
        case 'flag_display'
            flag_display=option_value;
        case 'ts_name'
            ts_name=option_value;
        case 'g_threshold'
            g_threshold=option_value;
        case 'g_threshold_p'
            g_threshold_p=option_value;
        case 'flag_pairwise'
            flag_pairwise=option_value;
        case 'flag_robust'
            flag_robust=option_value;
    end;
end;

if(~iscell(v)) %univariate version
    for to_idx=1:size(v,2)
        for from_idx=1:size(v,2)
            if(to_idx==from_idx)
                granger_num(to_idx,from_idx)=nan;
                granger_den(to_idx,from_idx)=nan;
                granger(to_idx,from_idx)=0;
                granger_F(to_idx,from_idx)=1.0;
                granger_pvalue(to_idx,from_idx)=0.5;
            else

                [w_1,A_1,C_1,SBC_1,FPE_1,th_1]=arfit(v(:,to_idx),ar_order,ar_order);
                [w_12,A_12,C_12,SBC_12,FPE_12,th_12]=arfit([v(:,to_idx), v(:,from_idx)],ar_order,ar_order);
                
                m     = 1;
                p     = size(A_1,2)/m;
                [siglev_1,res_1]=arres(w_1,A_1,v(:,to_idx),p+1);
                m     = 2;
                p     = size(A_12,2)/m;
                [siglev_12,res_12]=arres(w_12,A_12,[v(:,to_idx), v(:,from_idx)],p+1);


                T=size(v,1);
                if(flag_robust)
                    s=sort(res_1(:,1));
                    s=s(ceil(length(s)*0.25):floor(length(s)*0.95));
                    granger_num(to_idx,from_idx)=s(:)'*s(:);
                    
                    s=sort(res_12(:,1));
                    s=s(ceil(length(s)*0.25):floor(length(s)*0.95));
                    granger_den(to_idx,from_idx)=s(:)'*s(:);
                else
                    %m1=res_1(:,1:size(v{to_idx},2))'*res_1(:,1:size(v{to_idx},2))
                    %m2=res_12(:,1:size(v{to_idx},2))'*res_12(:,1:size(v{to_idx},2))
                    granger_num(to_idx,from_idx)=(res_1(:,1)'*res_1(:,1));
                    granger_den(to_idx,from_idx)=(res_12(:,1)'*res_12(:,1));
                end;
                granger_num_dof=(T-2*ar_order);
                granger_den_dof=(T-ar_order);
                granger(to_idx,from_idx)=log((res_1(:,1)'*res_1(:,1))/(res_12(:,1)'*res_12(:,1))*(T-2*ar_order)/(T-ar_order));
                granger_dir(to_idx,from_idx)=log((res_1(:,1)'*res_1(:,1))/(res_12(:,1)'*res_12(:,1)));
                granger_inst(to_idx,from_idx)=log(abs(C_12(1,1)).*abs(C_12(2,2))./det(C_12));
                granger_F(to_idx,from_idx)=(T-ar_order)./(ar_order).*(exp(granger(to_idx,from_idx))-(T-2*ar_order)/(T-ar_order));
                granger_pvalue(to_idx,from_idx)=1-fcdf(granger_F(to_idx,from_idx), ar_order, T-2*ar_order);
                
%close;
%plot(res_1(:,1),'r'); hold on;
%plot(res_12(:,1),'b'); hold on;

%granger
%save epi_vismotor_roi_050714 res_1 res_12 v granger_F granger  granger_pvalue;
%save ini_vismotor_roi_052314.mat v res_1 res_12 granger_F granger  granger_pvalue;
%show_epi_vismotor_roi_granger_050714
%
%keyboard;
            end;
        end;
    end;
else
     for to_idx=1:size(v,2)
        for from_idx=1:size(v,2)
            if(to_idx==from_idx)
                granger_num(to_idx,from_idx)=nan;
                granger_den(to_idx,from_idx)=nan;
                granger(to_idx,from_idx)=0;
                granger_F(to_idx,from_idx)=1.0;
                granger_pvalue(to_idx,from_idx)=0.5;
            else
                [w_1,A_1,C_1,SBC_1,FPE_1,th_1]=arfit(v{to_idx},ar_order,ar_order);
                [w_12,A_12,C_12,SBC_12,FPE_12,th_12]=arfit([v{to_idx}, v{from_idx}],ar_order,ar_order);
                
                m     = 1;
                p     = size(A_1,2)/m;
                [siglev_1,res_1]=arres(w_1,A_1,v{to_idx},p+1);
                m     = 2;
                p     = size(A_12,2)/m;
                [siglev_12,res_12]=arres(w_12,A_12,[v{to_idx} v{from_idx}],p+1);

                T=size(v{to_idx},1);
                
                
                if(flag_robust)
                    [s,idx]=sort(max(abs(res_1(:,1:size(v{to_idx},2))),[],2));
                    rr=res_1(:,1:size(v{to_idx},2));
                    rr=rr(idx(1:floor(length(idx)*0.95)),:);
                    granger_num(to_idx,from_idx)=det(rr'*rr);
                    
                    [s,idx]=sort(max(abs(res_12(:,1:size(v{to_idx},2))),[],2));
                    rr=res_12(:,1:size(v{to_idx},2));
                    rr=rr(idx(1:floor(length(idx)*1)),:);
                    granger_den(to_idx,from_idx)=det(rr'*rr);
                else
                    %m1=res_1(:,1:size(v{to_idx},2))'*res_1(:,1:size(v{to_idx},2))
                    %m2=res_12(:,1:size(v{to_idx},2))'*res_12(:,1:size(v{to_idx},2))
                    %det(m1)
                    %det(m2)
                    granger_num(to_idx,from_idx)=det(res_1(:,1:size(v{to_idx},2))'*res_1(:,1:size(v{to_idx},2)));
                    granger_den(to_idx,from_idx)=det(res_12(:,1:size(v{to_idx},2))'*res_12(:,1:size(v{to_idx},2)));
                end;
                %size(v{to_idx},2)
                %ar_order
                granger_num_dof=(T-2*ar_order*size(v{to_idx},2));
                granger_den_dof=(T-ar_order*size(v{to_idx},2));
                granger(to_idx,from_idx)=log(granger_num(to_idx,from_idx)/granger_den(to_idx,from_idx)*(T-2*ar_order*size(v{to_idx},2))/(T-ar_order*size(v{to_idx},2)));
                %granger(to_idx,from_idx)=log(granger_num(to_idx,from_idx)/granger_den(to_idx,from_idx)*(T-ar_order)/(T-ar_order))
                granger_dir(to_idx,from_idx)=log(granger_num(to_idx,from_idx)/granger_den(to_idx,from_idx));
                granger_inst(to_idx,from_idx)=log(det(C_12(1:size(v{to_idx},2),1:size(v{to_idx},2))).*det(C_12(size(v{to_idx},2)+1:end,size(v{to_idx},2)+1:end))./det(C_12));
                granger_F(to_idx,from_idx)=(T-ar_order)./(ar_order).*(exp(granger(to_idx,from_idx))-(T-2*ar_order)/(T-ar_order));
                granger_pvalue(to_idx,from_idx)=1-fcdf(granger_F(to_idx,from_idx), ar_order, T-2*ar_order);

            end;
        end;
    end;

    if(flag_display) fprintf('\n'); end;
end;

if(flag_display)
    if(isempty(ts_name))
        for i=1:size(granger,1)
            ts_name{i}=sprintf('node%2d',i);
        end;
    end;

    fprintf('\n\nSummarizing static Granger Causality...\n');
    fprintf('list from the most significant connection...\n\n');

    GR=granger;
    GR_p=granger_pvalue;

    if(isempty(g_threshold_p))
        g_threshold_p=0.01;
    end;

    [s_granger,s_idx]=sort(GR(:));
    for i=length(s_idx):-1:1
        [r(i),c(i)]=ind2sub(size(GR),s_idx(i));
        if(s_granger(i)>g_threshold_p)
            fprintf('<<%d>> from [%s] ---> to [%s] [(%d,%d)]: %3.3f (p-value=%3.3f)\n',length(s_idx)-i+1,ts_name{c(i)},ts_name{r(i)},r(i),c(i),s_granger(i),GR_p(s_idx(i)));
        end;
    end;
end;
