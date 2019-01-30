function MRSGetTissueCompFromMask_UI()

% MRSGetTissueFromMask_UI is designed to calculate the amount of CSF and 
% grey and white matter based on a MRS voxel mask (nifti) and the
% segmentation nifti files for GM, WM and CSF.
% 
% INPUT (via User Interfaces)
%      - nifti file of MRS voxle mask
%      - nifti file of grey matter mask of segmentation
%      - nifti file of white matter mask of segmentation
%      - nifti file of CSF mask of segmentation
%      - output folder
% 
% OUTPUT
%   - Text file containing the values:
%       GM_per: percentage of gray matter in the MRS voxel
%       WM_per: percentage of white matter in the MRS voxel
%       CSF_per: percentage of CSF in the MRS voxel
%       GM_per: Fraction of gray matter in the tissue of the MRS voxel
%       WM_per: Fraction of white matter in the tissue of the MRS voxel
% 
% LICENSE
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License (GPLv3) as published
% by the Free Software Foundation;
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY.
%
% AUTHOR
% Michael Lindner
% University of Reading, 2019
% School of Psychology and Clinical Language Sciences
% Center for Integrative Neuroscience and Neurodynamics
% https://www.reading.ac.uk/cinn/cinn-home.aspx


[MASKFile,MASKPath]=uigetfile('*.nii','Select MRS mask image');
[GMFile,GMPath]=uigetfile('*.nii','Select GM image',MASKPath);
[WMFile,WMPath]=uigetfile('*.nii','Select WM image',MASKPath);
[CSFFile,CSFPath]=uigetfile('*.nii','Select CSF image',MASKPath);
[outpath] = uigetdir(MASKPath, 'Select output folder');


mask_nii=load_untouch_nii([MASKPath MASKFile]);
GM_nii = load_untouch_nii([GMPath GMFile]);
WM_nii = load_untouch_nii([WMPath WMFile]);
CSF_nii = load_untouch_nii([CSFPath CSFFile]);



fprintf('\nMask segmentation data');
gm_mat = GM_nii.img .* mask_nii.img;
wm_mat = WM_nii.img .* mask_nii.img;
csf_mat = CSF_nii.img .* mask_nii.img;

fprintf('\nCalculate tissue composition');
gm = sum(gm_mat(:));
wm = sum(wm_mat(:));
csf = sum(csf_mat(:));

gm_per = gm / (gm+wm+csf);
wm_per = wm / (gm+wm+csf);
csf_per = csf / (gm+wm+csf);

gm_fra = gm / (gm+wm);
wm_fra = wm /(gm+wm);

fprintf(['\n\n', MASKFile]);
fprintf('\nTissue composition of the voxel :');
fprintf('\nGM_per = %1.2f',gm_per);
fprintf('\nWM_per = %1.2f',wm_per);
fprintf('\nCSF_per = %1.2f',csf_per);
fprintf('\nGM_fra = %1.2f',gm_fra);
fprintf('\nWM_fra = %1.2f\n',wm_fra);

txt_tiss_filename=[outpath,filesep, MASKFile(1:end-4), '_tiss_comp.gm'];
fileID=fopen(txt_tiss_filename,'wt');
fprintf(fileID,'Tissue composition of voxel: ' );
fprintf(fileID, MASKFile );
fprintf(fileID,'\n');
fprintf(fileID,'GM_per = %1.2f\n',gm_per);
fprintf(fileID,'WM_per = %1.2f\n',wm_per);
fprintf(fileID,'CSF_per = %1.2f\n',csf_per);
fprintf(fileID,'GM_fra = %1.2f\n',gm_fra);
fprintf(fileID,'WM_fra = %1.2f\n',wm_fra);

fclose(fileID);

GMper = gm_per;
WMper = wm_per;
CSFper = csf_per;
GMfrac = gm_fra;
WMfrac = wm_fra;
        

mat_out_filename = [outpath,filesep, MASKFile(1:end-4) '_tiss_comp.mat'];
save(mat_out_filename, 'GMper', 'WMper', 'CSFper',...
    'GMfrac', 'WMfrac' )


