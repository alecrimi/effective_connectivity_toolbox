%% test of GC and eGC analysis on simulated data
% This implements the simulation of Schiatti L, Nollo G, Rossato G, and Faes L 2015 Extended Granger causality: a new tool to identify the structure of physiological networks. Physiological Measurement

clear all; close all; clc;
M=3; %number of simulated time series
N=300; %length of simulated time series
pmax=2; % maximum lag in the simulation
insteff='n'; % 'n' for strictly causal, 'y' for extended with instantaneous effects
numsimu=10; %number of simulations

alpha=0.01; % statistical significance for F-test on GC and eGC

% Parameters of eGC analysis
numboot=1000;   % number of bootstrap samples
LRmethod=1; % which pairwise measure to be chosen [see Hyvarinen and Smith 2013]

%% eMVAR and MVAR parameters
r1=0.8; f1=0.1; % AR process AR for Y2
r2=0.9; f2=0.3; % AR process AR for Y1

a1=2*r1*cos(2*pi*f1);
a2=-r1^2;
b1=2*r2*cos(2*pi*f2);
b2=-r2^2;
c1=0.5;
Sw=eye(M);
%sets instantaneous effects
if insteff=='y'
    b0=0.8; 
else
    b0=0;
end
% matrix of instantaneous effects
B0(1,:)=[0 b0 0];
B0(2,:)=[0 0 0];
B0(3,:)=[b0 0 0];
% matrices of lagged effects
Bk=NaN*ones(M,M,pmax);
% effects at lag 1
Bk(1,:,1)=[b1 0 0];
Bk(2,:,1)=[0 a1 0];
Bk(3,:,1)=[c1 0 0];
% effects at lag 2
Bk(1,:,2)=[b2 0 0];
Bk(2,:,2)=[0.5 a2 c1];
Bk(3,:,2)=[0 0 0];
% concatenation in Bm matrix
Bm=[];
for kk=1:pmax
    Bm=[Bm Bk(:,:,kk)];
end

M_or = ones(3,3) - diag([1,1,1]); 


% corresponding strictly causal parameters (Eq. 13 Phys Meas paper)
[Am, Su]=egc_diag_coeff_rev(Bm,B0,Sw);

% stability check
p=pmax;
E=eye(M*p);AA=[Am;E(1:end-M,:)];
lambda=eig(AA);lambdamax=max(abs(lambda));

%% theoretical values of GC and eGC
numlags=p; %number of lags for which the atocovariance of the VAR process is computed
GCtrue=NaN*ones(M,M); eGCtrue=GCtrue;
for jj=1:M
    for ii=1:M
        if ii~=jj
            ret1 = egc_CElin_analytic(Am,Su,jj,ii,numlags);
            GCtrue(jj,ii)=2*(ret1.Hy_yz-ret1.Hy_yzx);
            
            tmp=B0(jj,[ii setdiff(1:M,[ii jj])]);
            flag0=tmp>0;
            ret2 = egc_CElin_analytic_0(Am,Su,jj,ii,numlags,flag0);
            eGCtrue(jj,ii)=2*(ret2.Hy_yz-ret2.Hy_yzx);
        end
    end
end

%% estimation from several realizations

for i=1:numsimu
    disp(['Simulation ' int2str(i) ' of ' int2str(numsimu)]);

    % generation of simulated time series
    U=egc_InstModelfilter(N,Sw,'ExtendedNonGauss',B0); % non gaussian innovations with covariance Su
    Y=egc_MVARfilter(Am,U); % realization 

    %estimation of GC
   % [GC1,p_val1,Am1,Su1,SigmaR1,Res]=egc_gcMVAR(Y,p);
    GC1= simple_structured_g_causality(Y,M_or);
%Add our Constrained GC

%Add stochastic DCM
    
    % estimation of extended GC
  %  [EGC1,pval_e1,Bm1,B01,Sw1,SigmaR2,ZeroLag1,LRpruned,LR]=egc_gceMVAR(Y,p,Res,LRmethod,numboot);
     
    GC(:,:,i)=GC1;
  %  EGC(:,:,i)=EGC1;
  %  p_GC(:,:,i)=p_val1;
  %  p_EGC(:,:,i)=pval_e1;
   % ZeroLag(:,:,i)=ZeroLag1;
end
%{
%% counts nonzero causality values and tests differences btw GC and eGC
ind=0;
for i=1:M
    for j=1:M
    ind=ind+1;
    p_GC_sum(i,j)=sum(p_GC(i,j,:)<alpha);
    p_EGC_sum(i,j)=sum(p_EGC(i,j,:)<alpha);
    end
end
%}

GCtrue(isnan(GCtrue))=0;

dif = zeros(1,numsimu);

for kk = 1 : numsimu
    dif(kk) = norm(GCtrue - GC(:,:,kk) );
end
 

