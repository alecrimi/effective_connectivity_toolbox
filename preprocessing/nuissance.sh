#!/usr/bin/env bash

nuisance_dir='nuissance'
segment_dir='segment'
anat_dir='anat'
#standard='../../../../templates/MNI152_T1_3mm_brain.nii.gz'
rest="BOLDrestingCAP"
nuisance_template='../../../../../../../../templates/nuisance.fsf'
      
## 1. make nuisance directory
rm -rf ${nuisance_dir}
mkdir -p ${nuisance_dir}

## 2. Generate Motion correction timeseries values
echo "Motion correction data"
rm ${rest}_mc.nii.gz 
3dTstat -mean -prefix ${rest}_ro_mean.nii.gz ${rest}_pp.nii.gz 
3dvolreg -Fourier -twopass -base ${rest}_ro_mean.nii.gz -zpad 4 -prefix ${rest}_mc.nii.gz -1Dfile ${rest}_mc.1D ${rest}_pp.nii.gz

## 3. Seperate motion parameters into seperate files
echo "Splitting up motion parameters"
awk '{print $1}' ${rest}_mc.1D > ${nuisance_dir}/mc1.1D
awk '{print $2}' ${rest}_mc.1D > ${nuisance_dir}/mc2.1D
awk '{print $3}' ${rest}_mc.1D > ${nuisance_dir}/mc3.1D
awk '{print $4}' ${rest}_mc.1D > ${nuisance_dir}/mc4.1D
awk '{print $5}' ${rest}_mc.1D > ${nuisance_dir}/mc5.1D
awk '{print $6}' ${rest}_mc.1D > ${nuisance_dir}/mc6.1D
 
: '
## 4. Global regression
echo "Extracting global signal"
#fslmaths ${rest}_pp.nii.gz -Tmin -bin ${rest}_pp_mask.nii.gz -odt char
3dcopy  ${rest}_pp_mask.nii.gz  ${segment_dir}/global_mask.nii.gz #
3dmaskave -mask ${segment_dir}/global_mask.nii.gz -quiet ${rest}_pp.nii.gz > ${nuisance_dir}/global.1D
'
 
################################################Generate regressors and perform nuissance regression #############
#Clean previous data
rm -rf ${segment_dir}
mkdir ${segment_dir}

echo "Register T1 to data and segment"
bet  ../../11/DEFACED_NIFTI/defaced_MPRAGE.nii.gz    ${rest}_T1_bet.nii.gz -m -R
flirt -ref ${rest}_pp.nii.gz -in ${rest}_T1_bet.nii.gz -o ${segment_dir}/${rest}_T1_reg.nii.gz
fast -g --nopve -t 1  ${segment_dir}/${rest}_T1_reg.nii.gz 

## 19. csf
echo "Extracting signal from csf "
3dmaskave -mask ${segment_dir}/${rest}_T1_reg_seg_0.nii.gz -quiet ${rest}_pp.nii.gz > ${nuisance_dir}/csf.1D

## 20. wm
echo "Extracting signal from white matter  "
3dmaskave -mask ${segment_dir}/${rest}_T1_reg_seg_2.nii.gz -quiet ${rest}_pp.nii.gz > ${nuisance_dir}/wm.1D

## 21. Generate mat file (for use later)
## create fsf file
echo "Modifying model file"
sed -e s:nuisance_dir:"${nuisance_dir}":g <${nuisance_template} >${nuisance_dir}/temp1
sed -e s:nuisance_model_outputdir:"${nuisance_dir}/residuals.feat":g <${nuisance_dir}/temp1 >${nuisance_dir}/temp2
sed -e s:nuisance_model_TR:"${TR}":g <${nuisance_dir}/temp2 >${nuisance_dir}/temp3
sed -e s:nuisance_model_numTRs:"${n_vols}":g <${nuisance_dir}/temp3 >${nuisance_dir}/temp4
sed -e s:nuisance_model_input_data:"${func_dir}/${rest}_pp.nii.gz":g <${nuisance_dir}/temp4 >${nuisance_dir}/nuisance.fsf 

## Framewise displacement from Matlab
cp ../../../../../../../../*.m .
matlab -nodesktop -nosplash - nojvm -r "run bramila_framewiseDisplacement.m" 

rm ${nuisance_dir}/temp*
3dDeconvolve \
     -input ${rest}_pp.nii.gz -polort 0 \
     -mask ${rest}_pp_mask.nii.gz \
     -num_stimts 9 \
     -stim_file 1 ${nuisance_dir}/mc1.1D -stim_label 1 mot0    -stim_base 1 \
     -stim_file 2 ${nuisance_dir}/mc2.1D -stim_label 2 mot1    -stim_base 2 \
     -stim_file 3 ${nuisance_dir}/mc3.1D -stim_label 3 mot2    -stim_base 3 \
     -stim_file 4 ${nuisance_dir}/mc4.1D -stim_label 4 mot3    -stim_base 4 \
     -stim_file 5 ${nuisance_dir}/mc5.1D -stim_label 5 mot4    -stim_base 5 \
     -stim_file 6 ${nuisance_dir}/mc6.1D -stim_label 6 mot5    -stim_base 6 \
     -stim_file 7 ${nuisance_dir}/wm.1D -stim_label 7 WM -stim_base 7 \
     -stim_file 8 ${nuisance_dir}/csf.1D -stim_label 8 CSF -stim_base 8 \
     -stim_file 9 FD.1D -stim_label 9 FD -stim_base 9 \
     -x1D ${nuisance_dir}/nuisance.fsf \
     -x1D_stop
#     -stim_file 9 ${nuisance_dir}/global.1D -stim_label 9 global -stim_base 9 \
#     -stim_file 10 FD.1D -stim_label 10 FD -stim_base 10 \

 
echo "Running 3dREMlfit to get residuals"
3dREMLfit -input ${rest}_pp.nii.gz \
    -matrix ${nuisance_dir}/nuisance.fsf.xmat.1D \
    -mask ${rest}_pp_mask.nii.gz \
    -Rerrts ${nuisance_dir}/res4d.$cor.nii.gz \
    -GOFORIT

echo "Demeaning residuals "
3dTstat -mean -prefix ${nuisance_dir}/res4d_mean.$cor.nii.gz ${nuisance_dir}/res4d.$cor.nii.gz
3dcalc -a ${nuisance_dir}/res4d.$cor.nii.gz -b ${nuisance_dir}/res4d_mean.$cor.nii.gz -expr '(a-b)' -prefix ${rest}_pp_final.nii    #?


