function [granger, granger_F, granger_pvalue, granger_num, granger_den, granger_num_dof, granger_den_dof, granger_dir]= etc_granger_constrained(v,ar_order,A,varargin)

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

flag_display=1;
ts_name={};
g_threshold=[];
g_threshold_p=[];
flag_pairwise=1;
flag_robust=0;

% Time points of the series
T=size(v,1);  
 

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


%Consider two models, 1 and 2, where model 1 is 'nested' within model 2. Model 1 is the restricted model, and model 2 is the unrestricted one. That is, model 1 has p1 parameters, and model 2 has p2  parameters, where p1 < p2, and for any choice of parameters in model 1, the same regression curve can be achieved by some choice of the parameters of model 2.

% The assumption is that model1 is the univariate autoregression for 1 time-series (model1), and model2 is the multivariate regression for the time-series of 2 ROIs (model2).

if(~iscell(v)) 
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
 
            %    m     = 1;
            %    p     = size(A_1,2)/m;
            %    [siglev_1,res_1]=arres(w_1,A_1,v(:,to_idx),p+1);
            %    m     = 2;
            %    p     = size(A_12,2)/m;
            %    [siglev_12,res_12]=arres(w_12,A_12,[v(:,to_idx), v(:,from_idx)],p+1);
            
            % Compute residuals for the known coefficient matrix 
            res_1 = zeros(T-1,1);
            for jj = 1 :  T -1
                res_1(jj)  =     norm( A(to_idx,to_idx) * v(jj,to_idx)- v(jj+1,to_idx) , 2 ) ;
            end   
 
            %A_12 = [ A(to_idx,to_idx) A(to_idx,to_idx); A(to_idx,from_idx) A(to_idx,from_idx)   ];
            %    res_12 = arres(0,A_12,[v(:,to_idx), v(:,from_idx)],2);

            res_12 = zeros(T-1,1);
            for jj = 1 :  T -1
                res_12(jj)  =   norm( A(to_idx,to_idx) * v(jj,to_idx) +  A(to_idx,from_idx) * v(jj,from_idx) - v(jj+1,to_idx) , 2 ) ;
            end  
            
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



                granger_dir(to_idx,from_idx)=log((res_1(:,1)'*res_1(:,1))/(res_12(:,1)'*res_12(:,1)));
                granger_inst(to_idx,from_idx)=log(abs(C_12(1,1)).*abs(C_12(2,2))./det(C_12));

                % F-statistics
                granger(to_idx,from_idx)=log( (res_1(:,1)'*res_1(:,1))/(res_12(:,1)'*res_12(:,1))*(T-2*ar_order)/(T-ar_order) );
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
        g_threshold_p=0.05;
    end;

    [s_granger,s_idx]=sort(GR(:));
    for i=length(s_idx):-1:1
        [r(i),c(i)]=ind2sub(size(GR),s_idx(i));
        if(s_granger(i)>g_threshold_p)
            fprintf('<<%d>> from [%s] ---> to [%s] [(%d,%d)]: %3.3f (p-value=%3.3f)\n',length(s_idx)-i+1,ts_name{c(i)},ts_name{r(i)},r(i),c(i),s_granger(i),GR_p(s_idx(i)));
        end;
    end;
end;
