%Script to test the Autoregressive model on a list of subjects

name_session = '../NKI.105.001.NIFTI';

list(1,:) = 'MR_1427581_3b35';
list(2,:) = 'MR_1523112_43e6';
list(3,:) = 'MR_1616468_54cf';
list(4,:) = 'MR_2021454_4ade';
list(5,:) = 'MR_2328270_6ebf';
list(6,:) = 'MR_2505567_3f9e';
list(7,:) = 'MR_2541244_5224';
list(8,:) = 'MR_2784584_1a4d';
list(9,:) = 'MR_2861923_4f24';
list(10,:) = 'MR_3374719_104f';   
list(11,:) = 'MR_3466763_386d';   
list(12,:) = 'MR_3611212_36f0';   
list(13,:) = 'MR_4097119_886b';    
list(14,:) = 'MR_4143704_4f12';    
list(15,:) = 'MR_4290056_39b5';     
list(16,:) = 'MR_4323037_4dd1';      
list(17,:) = 'MR_6539040_7cc6';       
list(18,:) = 'MR_6692978_5233';       
list(19,:) = 'MR_6913939_6047';           
list(20,:) = 'MR_9421819_1328';        

[r c] = size(list); 


 
for ll = 1 : r
    ll
    structured_G_causality( name_session, list(ll,:) , 2 );
end
