 #!/bin/bash  

## Set high pass and low pass cutoffs for temporal filtering
hp=0.005
lp=0.1


for d in */ ; do
    echo "-------------------------------------------------"
    echo "entering subfolder" "$d"
    cd $d
       cd experiments
       sub1=*/
       echo $sub1
          cd $sub1
              cd $sub1
                 cd SCANS

                 ################## Pre-processing BOLD  ###############
                 echo " ################################################"
                 echo "Pre-processing fMRI"
                 cd 4    
                 cd NIFTI   
                 rest="BOLDrestingCAP"
 
                 #echo "Deobliquing"
                 3drefit -deoblique ${rest}.nii.gz       

                 echo "Motion correcting"
                 mcflirt -in ${rest}.nii.gz -o ${rest}_mc.nii.gz

                 echo "Skull stripping"
                 bet ${rest}_mc.nii.gz  brain.nii.gz -m -R -f 0.4
                 fslmaths ${rest}_mc.nii.gz -mas brain_mask.nii.gz ${rest}_bet.nii.gz 

                 echo "Band-pass filtering"
                 3dFourier -lowpass ${lp} -highpass ${hp} -retrend -prefix ${rest}_filt.nii.gz ${rest}_bet.nii.gz #era gms

                 echo "Removing linear and quadratic trends"
                 3dTstat -mean -prefix ${rest}_filt_mean.nii.gz ${rest}_filt.nii.gz
                 3dDetrend -polort 2 -prefix ${rest}_dt.nii.gz ${rest}_filt.nii.gz #era 
                 3dcalc -a ${rest}_filt_mean.nii.gz -b ${rest}_dt.nii.gz -expr 'a+b' -prefix ${rest}_pp.nii.gz
 
                 echo "Generating mask of preprocessed data"
                 fslmaths ${rest}_pp.nii.gz  -Tmean -bin ${rest}_pp_mask.nii.gz -odt char  


                 echo "Nuissance removal"
                 sh ../../../../../../../nuissance.sh                                 

                 echo "Framewise displacement"
                 #fsl_motion_outliers -i ${rest}_pp_final.nii -o ${rest}_fd_final.nii --fd
                
                 
                 #  generate a registration map using the refence volume of the atlas
                 flirt -ref ${rest}_pp_final.nii  -in  ../../../../../../../../MNI152lin_T1_1mm_brain.nii.gz -omat my_mat.mat

                 # register the atlas using the generated matrix
                 flirt -ref ${rest}_pp_final.nii  -in  ../../../../../../../../aal.nii -applyxfm -init my_mat.mat -out atlas_reg.nii.gz -interp nearestneighbour


                 cd ..
                 cd ..  
                 #############################################
 
                 ################## Pre-processing DTI  ###############
                 echo " ################################################"
                 echo "Pre-processing DTI"
                 cd 9
                 cd NIFTI   
                 rest="DTICAP"

                 echo "Eddy current correction"
                 eddy_correct ${rest}.nii.gz  ${rest}_eddycorrected.nii.gz   -interp trilinear

                 echo "Deobliquing ${subject}"
                 3drefit -deoblique ${rest}.nii.gz

                 echo "Reorienting ${subject}"
                 3dresample -orient RPI -inset ${rest}.nii.gz -prefix ${rest}_ro.nii.gz

                 echo "Skull stripping ${subject}"
                 #bet  ${rest}_eddycorrected.nii.gz    brain.nii.gz -m -R
                 bet  ${rest}_ro.nii.gz    brain.nii.gz -m -R -f 0.4
                 fslmaths ${rest}_eddycorrected.nii.gz   -mas brain_mask.nii.gz ${rest}_bet.nii.gz

                 echo "generate a registration map using the refence volume of the atlas"
                 flirt -ref ${rest}_bet.nii.gz -in  ../../../../../../../../MNI152lin_T1_1mm_brain.nii.gz -omat my_mat.mat

                 echo "register the atlas using the generated matrix"
                 flirt -ref ${rest}_bet.nii.gz -in  ../../../../../../../../aal.nii -applyxfm -init my_mat.mat -out atlas_reg.nii.gz -interp nearestneighbour

 
 
                 echo "Generate Connectome matrix"
                 cp ../../../../../../../../*.py .
                 python  visualisation.py 
                 cd ..
                 cd ..
 
                 cd ..
              cd ..
          cd ..
      
       cd ..
    cd ..
done
