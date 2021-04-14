	close all; clear all;

n=1000;
I=0.5; %strong influence
A=[[-0.9 0; I -0.9],[-0.3 0; I*0.5 -0.3]];
S=[1 0; 0 1];
D=10; %additional time points to reach "steady state"

%calling AR model to simulate the time series;
v=arsim(zeros(1,2),A,S,n+2000+D);
v=v(1:n,:);

v(:,1)=fmri_scale(v(:,1),1,-1);
v(:,2)=fmri_scale(v(:,2),1,-1);

%simulate noise time series;
vn=randn(size(v));
vn(:,1)=fmri_scale(vn(:,1),1,-1);
vn(:,2)=fmri_scale(vn(:,2),1,-1);

%generate dynamics
envelop1=1./(1+exp(1.*([1:n]-600)./10)); 
envelop2=1-1./(1+exp(1.*([1:n]-300)./10)); 
envelop=(envelop1.*envelop2)';
hold on; plot(envelop,'r');

%generate dynamic models
v(:,1)=envelop.*v(:,1)+(1-envelop).*vn(:,1);
v(:,2)=envelop.*v(:,2)+(1-envelop).*vn(:,2);

[AIC, AIC1, BIC]=etc_ar_sure(v,1, [100 300 800],[1:6]);

return;
subplot(211); hold on; plot(v(:,1),'r'); %time series "X";
legend({'X'});
subplot(212); hold on; plot(v(:,2),'b'); %time series "Y";
legend({'Y'});
xlabel('time (sample)');

%calculate the granger
[granger, granger_F, granger_p]=etc_granger(v,1);
