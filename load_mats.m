

matrices = zeros(90,90,20);

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/1427581/experiments/MR_1427581_3b35/MR_1427581_3b35/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,1) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,1) = 0.5 * error_des;%etem(end);
 
%%etem_list(1,:) = etem;
bold_seqs(:,:,1) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/1523112/experiments/MR_1523112_43e6/MR_1523112_43e6/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,2) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,2) =  0.5 * error_des;%etem(end);
%%etem_list(2,:) = etem;
bold_seqs(:,:,2) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/1616468/experiments/MR_1616468_54cf/MR_1616468_54cf/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,3) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,3) =  0.5 * error_des;%etem(end);
%%etem_list(3,:) = etem;
bold_seqs(:,:,3) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/2021454/experiments/MR_2021454_4ade/MR_2021454_4ade/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,4) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,4) =  0.5 * error_des;%etem(end);
%%etem_list(4,:) = etem;
bold_seqs(:,:,4) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/2328270/experiments/MR_2328270_6ebf/MR_2328270_6ebf/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,5) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,5) =  0.5 * error_des;%etem(end);
%%etem_list(5,:) = etem;
bold_seqs(:,:,5) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/2505567/experiments/MR_2505567_3f9e/MR_2505567_3f9e/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,6) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,6) =  0.5 * error_des;%etem(end);
%%etem_list(6,:) = etem;
bold_seqs(:,:,6) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/2541244/experiments/MR_2541244_5224/MR_2541244_5224/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,7) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,7) =  0.5 * error_des;%etem(end);
%etem_list(7,:)= etem;
bold_seqs(:,:,7) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/2784584/experiments/MR_2784584_1a4d/MR_2784584_1a4d/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,8) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,8) =  0.5 * error_des;%etem(end);
%etem_list(8,:) = etem;
bold_seqs(:,:,8) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/2861923/experiments/MR_2861923_4f24/MR_2861923_4f24/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,9) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,9) =  0.5 * error_des;%etem(end);
%etem_list(9,:) = etem;
bold_seqs(:,:,9) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/3374719/experiments/MR_3374719_104f/MR_3374719_104f/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,10) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,10) = 0.5 * error_des;%etem(end);
%etem_list(10,:) = etem;
bold_seqs(:,:,10) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/3466763/experiments/MR_3466763_386d/MR_3466763_386d/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,11) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,11) =  0.5 * error_des;%etem(end);
%etem_list(11,:) = etem;
bold_seqs(:,:,11) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/3611212/experiments/MR_3611212_36f0/MR_3611212_36f0/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,12) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,12) =  0.5 * error_des;%etem(end);
%etem_list(12,:) = etem;
bold_seqs(:,:,12) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/4097119/experiments/MR_4097119_886b/MR_4097119_886b/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,13) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,13) =  0.5 * error_des;%etem(end);
%etem_list(13,:) = etem;
bold_seqs(:,:,13) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/4143704/experiments/MR_4143704_4f12/MR_4143704_4f12/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,14) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,14) =  0.5 * error_des;%etem(end);
%etem_list(14,:) = etem;
bold_seqs(:,:,14) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/4290056/experiments/MR_4290056_39b5/MR_4290056_39b5/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,15) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,15) =  0.5 * error_des;%etem(end);
%etem_list(15,:) = etem;
bold_seqs(:,:,15) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/4323037/experiments/MR_4323037_4dd1/MR_4323037_4dd1/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,16) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,16) =  0.5 * error_des;%etem(end);
%etem_list(16,:) = etem;
bold_seqs(:,:,16) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/6539040/experiments/MR_6539040_7cc6/MR_6539040_7cc6/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,17) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,17) =  0.5 * error_des;%etem(end);
%etem_list(17,:) = etem;
bold_seqs(:,:,17) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/6692978/experiments/MR_6692978_5233/MR_6692978_5233/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,18) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,18) =  0.5 * error_des;%etem(end);
%etem_list(18,:) = etem;
bold_seqs(:,:,18) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/6913939/experiments/MR_6913939_6047/MR_6913939_6047/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,19) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,19) =  0.5 * error_des;%etem(end);
%etem_list(19,:) = etem;
bold_seqs(:,:,19) = averaged_bold_seq;

connectome_folder='/media/acrimi/TOSHIBA EXT1/1000Connectome/NKI.105.001.NIFTI/9421819/experiments/MR_9421819_1328/MR_9421819_1328/SCANS/9/NIFTI/';
load([connectome_folder 'eff_conn0.mat']);
matrices(:,:,20) = M;
error_des= 0;
 for jj = 5 :  t -1
       error_des  = error_des +   norm(  M  * averaged_bold_seq(jj,:)'- averaged_bold_seq(jj+1,:)' , 2 ) ;
 end
el(:,20) = 0.5 * error_des;%etem(end);
%etem_list(20,:) = etem;
bold_seqs(:,:,20) = averaged_bold_seq;
save('bold_seqs.mat','bold_seqs');
save('after_autoregresive.mat','matrices');
%{
matrices_bin = matrices>0;
tot_mat = sum(matrices_bin,3);
%figure; imagesc(tot_mat);
figure; imagesc(mean(matrices,3))



plot(el)


matrices(1:20,:,20)= matrices(:,1:20,20)';
matrices(1:20,:,19)= matrices(:,1:20,19)';
matrices(1:20,:,18)= matrices(:,1:20,18)';
matrices(1:20,:,17)= matrices(:,1:20,17)';
matrices(1:20,:,16)= matrices(:,1:20,16)';
matrices(1:20,:,15)= matrices(:,1:20,15)';
matrices(1:20,:,14)= matrices(:,1:20,14)';
matrices(1:20,:,13)= matrices(:,1:20,13)';
matrices(1:20,:,12)= matrices(:,1:20,12)';
matrices(1:20,:,11)= matrices(:,1:20,11)';
matrices(1:20,:,10)= matrices(:,1:20,10)';
matrices(1:20,:,9)= matrices(:,1:20,9)';
matrices(1:20,:,8)= matrices(:,1:20,8)';
matrices(1:20,:,7)= matrices(:,1:20,7)';
matrices(1:20,:,6)= matrices(:,1:20,6)';
matrices(1:20,:,5)= matrices(:,1:20,5)';
matrices(1:20,:,4)= matrices(:,1:20,4)';
matrices(1:20,:,3)= matrices(:,1:20,3)';
matrices(1:20,:,2)= matrices(:,1:20,2)';
matrices(1:20,:,1)= matrices(:,1:20,1)';
save('after_autoregresive.mat','matrices');
 

y = prctile(mean_M(:),90)
mask = mean_M>y;


[phat, pci] = gamfit(vec);
x = gaminv(0.1, phat(1), phat(2));
y= x

 
hold on
plot(el,'go','MarkerSize',8)
plot(el_grang-50,'bx','MarkerSize',8)
plot(el_dti,'rs','MarkerSize',8)
xlabel('Subjects') % x-axis label
ylabel('Reconstruction Error') % y-axis label
box on
 %}
