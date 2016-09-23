%function [fwd,rms]=bramila_framewiseDisplacement(cfg)
function [fwd,rms]=bramila_framewiseDisplacement()
% BRAMILA_FRAMEWISEDISPLACEMENT - Computes the framewise displacement
% metric as described in 
% Power et al. (2012) doi:10.1016/j.neuroimage.2011.10.018 and also 
% Power et al. (2014) doi:10.1016/j.neuroimage.2013.08.048
%   - Usage:
%   fwd = bramila_framewiseDisplacement(cfg)
%   - Input:
%   cfg is a struct with following parameters
%       cfg.motionparam = the 6 time series of motion parameters (time in 1st dimension)
%       cfg.prepro_suite = 'fsl-fs', 'spm' (default fsl-fs, fs = freesurfer)
%       cfg.radius = radius of sphere in mm to convert degrees to motion,
%       default = 50 as in Power et al 2014
%   - Output:
%       fwd = framewise displacement timeseries
%   - Notes:
%   Need to check that spm is indeed different, see end of Yan 2013 10.1016/j.neuroimage.2013.03.004 
%
%	Last edit: EG 2010-01-10
    
     
    fprintf('Computing framewise displacement...');
    cfg = 0;
    temp_cfg=[];
    ts  = importdata('BOLDrestingCAP_mc.1D'); %    ts = load(cfg.motionparam);
    temp_cfg.vol=double(ts);
    temp_cfg.write=0;
    temp_cfg.detrend_type='linear-demean';
    ts=bramila_detrend(temp_cfg);   % demean and detrend as specified in Power et al 2014
    
    if(size(ts,2)~=6)
        error(['The motion time series must have 6 motion parameters in 6 columns; the size of the input given is ' num2str(size(ts))])
    end
    prepro_suite='fsl-fs'; % 1 is FSL, 2 is SPM
    if(isfield(cfg,'prepro_suite'))
        prepro_suite = cfg.prepro_suite;
    end
    
    radius=50; % default radius
    if(isfield(cfg,'radius'))
        radius = cfg.radius;
    end    
    
    if(strcmp(prepro_suite,'fsl-fs'))
        % convert radians into motion in mm
		% in FSL the first 3 columns are rotations, in spm its viceversa
        temp=ts(:,1:3);
        temp=radius*temp;
        ts(:,1:3)=temp;
    else % SPM way
        % convert degrees into motion in mm;
        temp=ts(:,4:6);
        %temp=(2*radius*pi/360)*temp; % UPDATE: it seems they are in radians afterall
        temp=radius*temp;
        ts(:,4:6)=temp;
    end
    
    dts=diff(ts);
    dts=[
        zeros(1,size(dts,2)); 
        dts
        ];  % first element is a zero, as per Power et al 2014
    fwd=sum(abs(dts),2);
    rms=sqrt(mean(ts.^2,1));    % root mean square for each column
    
    fileID = fopen('FD.1D','w');
    fprintf(fileID,'%f\n', fwd);
    fclose(fileID);
    fprintf(' done\n');
    exit
