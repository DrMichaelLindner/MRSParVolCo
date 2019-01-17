function MRSGetTissueComp(batch)

% MRSGetTissueComp is designed to calculate the amount of CSF and grey and 
% white matter based on a structural nifti files and SIEMENS dicom or rda 
% MRS file.
%
% INPUT
%   Matlab structure of the following format:
%       batch(x).T1file = nifti file of structrual scan after brain 
%                           extraction (bet)
%       batch(x).GMfile = nifti file of grey matter mask of segmentation
%       batch(x).WMfile = nifti file of white matter mask of segmentation
%       batch(x).CSFfile = nifti file of CSF mask of segmentation
%       batch(x).MRSfiles{y} = SIEMENS rda or dicom MRS file
%       batch(x).outpath = output folder
% 
%   MRSParVolCo_1_get_tissue can loop over multiple subjects (x) and 
%   multiple MRS files per subject (y).
%
% OUTPUT
%   - A Mask of the MRS sinlge voxel 
%   - Text file containing the center coordinates of the MRS voxel (in
%     single subject space)
%   - Text file containing the values:
%       GM_per: percentage of gray matter in the MRS voxel
%       WM_per: percentage of white matter in the MRS voxel
%       CSF_per: percentage of CSF in the MRS voxel
%       GM_frac: Fraction of gray matter in the tissue of the MRS voxel
%       WM_frac: Fraction of white matter in the tissue of the MRS voxel
%   - Matlab file containing the vectors GM_per, WM_per, CSF_per, GM_frac
%       WM_frac(over subjects). They can be uses as input for 
%       MRSParVolCo_2_correct.
%
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
%


outcount = 0;
GMper = [];
WMper = [];
CSFper = [];
GMfrac = [];
GMfrac = [];

all_tiss_filename=[batch(1).outpath,filesep, 'all_tiss_comp.txt'];
fileIDall=fopen(all_tiss_filename,'wt');       
fprintf(fileIDall,'GMper\tWMper\tCSFper\tGMfrac\tWMfrac\n');

for ii=1:length(batch)
    
    % get files
    t1File = batch(ii).T1file;
    GMFile = batch(ii).GMfile;
    WMFile = batch(ii).WMfile;
    CSFFile = batch(ii).CSFfile;
    MRSFiles = batch(ii).MRSfiles;
    outpath = batch(ii).outpath;
    
    fprintf('\nGet info from T1 file');
    [n1,n2,n3,p1,p2,p3,ps_r,ps_c,r,c] = get_t1_infos(t1File);
    
    fprintf('\nLoad segmentation data');
    GM_nii = load_untouch_nii(GMFile);
    WM_nii = load_untouch_nii(WMFile);
    CSF_nii = load_untouch_nii(CSFFile);
    
    
    for fff=1:length(MRSFiles)
        
        MRSFile = MRSFiles{fff};
        [filepath,MRSfilename,ext] = fileparts(MRSFile);
        
        fprintf('\nGet info from MRS ');
        if strcmp(MRSFile(end-2:end),'dcm')
            fprintf('dicom file');
            a = get_dicom_info(MRSFile);
        elseif strcmp(MRSFile(end-2:end),'rda')
            fprintf('rda file');
            a = get_rda_info(MRSFile);
        end
        
        fprintf('\n\nMRS coordinates of voxel:');
        fprintf('\nSag: %2.4f',-a.pdSag);
        fprintf('\nCor: %2.4f',-a.pdCor);
        fprintf('\nTra: %2.4f\n',a.pdTra);
        
        fprintf('\n\nWrite text file');
        txtfilename=[outpath,filesep, MRSfilename, '_mrs_coord.txt'];
        fileID=fopen(txtfilename,'wt');
        fprintf(fileID,'MRS coordinates of voxel: ');
        fprintf(fileID, [filepath, filesep, MRSfilename, filesep, ext] );
        fprintf(fileID,'\n');
        fprintf(fileID,'%2.4f\n',-a.pdSag);
        fprintf(fileID,'%2.4f\n',-a.pdCor);
        fprintf(fileID,'%2.4f\n',a.pdTra);
        fclose(fileID);
        
        fprintf('\nGet coordinates of the VOI corners.');
        oct = create_VOI_nifti_siemens(a);
        
        fprintf('\nCalculate the coordinates of the intersections between VOI and the image planes.');
        co = intersection_nifti_siemens(n3,p1,p2,p3,oct);
        
        fprintf('\nCalculate the transition to the closest coordinates.');
        G = transition_nifti(p1,n1,n2,n3,ps_r,ps_c,r,c,co);
        
        fprintf('\nGenerate Mask');
        mask = create_filter_nifti(G,n1,n2,n3);
        
        fprintf('\nSave Mask');
        nii=load_untouch_nii(t1File);
        nii.img=mask;
        
        savename = [outpath, filesep 'Mask_',MRSfilename,ext];
        save_untouch_nii(nii,savename);
        
        fprintf('\nMask segmentation data');
        gm_mat = GM_nii.img .* mask;
        wm_mat = WM_nii.img .* mask;
        csf_mat = CSF_nii.img .* mask;
        
        fprintf('\nCalculate tissue composition');
        gm = sum(gm_mat(:));
        wm = sum(wm_mat(:));
        csf = sum(csf_mat(:));
        
        gm_per = gm / (gm+wm+csf);
        wm_per = wm / (gm+wm+csf);
        csf_per = csf / (gm+wm+csf);
        
        gm_fra = gm / (gm+wm);
        wm_fra = wm /(gm+wm);
        
        fprintf(['\n\n', MRSFile]);
        fprintf('\nTissue composition of the voxel :');
        fprintf('\nGM_per = %1.2f',gm_per);
        fprintf('\nWM_per = %1.2f',wm_per);
        fprintf('\nCSF_per = %1.2f',csf_per);
        fprintf('\nGM_fra = %1.2f',gm_fra);
        fprintf('\nWM_fra = %1.2f\n',wm_fra);
        
        txt_tiss_filename=[outpath,filesep, MRSfilename, '_tiss_comp.gm'];
        fileID=fopen(txt_tiss_filename,'wt');
        fprintf(fileID,'Tissue composition of voxel: ' );
        fprintf(fileID, MRSFile );
        fprintf(fileID,'\n');
        fprintf(fileID,'GM_per = %1.2f\n',gm_per);
        fprintf(fileID,'WM_per = %1.2f\n',wm_per);
        fprintf(fileID,'CSF_per = %1.2f\n',csf_per);
        fprintf(fileID,'GM_fra = %1.2f\n',gm_fra);
        fprintf(fileID,'WM_fra = %1.2f\n',wm_fra);
        
        fprintf(fileIDall,'%1.2f\t',gm_per);
        fprintf(fileIDall,'%1.2f\t',wm_per);
        fprintf(fileIDall,'%1.2f\t',csf_per);
        fprintf(fileIDall,'%1.2f\t',gm_fra);
        fprintf(fileIDall,'%1.2f\n',wm_fra);
        fclose(fileID);
        
        outcount = outcount + 1;
        order{outcount,1} = MRSFiles{fff};
        GMper(outcount,1) = gm_per;
        WMper(outcount,1) = wm_per;
        CSFper(outcount,1) = csf_per;
        GMfrac(outcount,1) = gm_fra;
        WMfrac(outcount,1) = wm_fra;
        
    end
end

fclose(fileIDall);

mat_out_filename = [outpath,filesep, 'all_tiss_comp.mat'];
save(mat_out_filename, 'order', 'GMper', 'WMper', 'CSFper',...
    'GMfrac', 'WMfrac' )


% -------------------------------------------------------------------------
% -------------------------- Nested functions -----------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% get T1 info
% -------------------------------------------------------------------------

    function [n1,n2,n3,p1,p2,p3,ps_r,ps_c,r,c] = get_t1_infos(t1File)
        
        nii=load_untouch_nii(t1File);
        
        n3=size(nii.img,3);
        n1=size(nii.img,1);
        n2=size(nii.img,2);
        ps_r=nii.hdr.dime.pixdim(2);
        ps_c=nii.hdr.dime.pixdim(3);
        mat=[nii.hdr.hist.srow_x;nii.hdr.hist.srow_y;nii.hdr.hist.srow_z;[0 0 0 1]];
        
        [R,C,P]  = ndgrid(1:n1,1:n2,1:n3);
        RCP      = [R(:)';C(:)';P(:)'];
        clear R C P
        RCP(4,:) = 1;
        XYZ      = mat(1:3,:)*RCP;
        
        l=length(XYZ);
        count=1:l;
        rcount=reshape(count,n1,n2,n3);
        
        p1=[];
        p2=[];
        p3=[];
        
        for i=1:n3
            p1=[p1 XYZ(:,rcount(1,1,i))]; %#ok<*AGROW>
            p2=[p2 XYZ(:,rcount(n1,1,i))];
            p3=[p3 XYZ(:,rcount(1,n2,i))];
        end
        
        r=(p1(:,1)-p2(:,1))/norm(p1(:,1)-p2(:,1));
        c=(p1(:,1)-p3(:,1))/norm(p1(:,1)-p3(:,1));
    end

% -------------------------------------------------------------------------
% get dicom info
% -------------------------------------------------------------------------

    function [a] = get_dicom_info(filename)
        
        dcm = dicominfo(filename);
        
        a.TR = dcm.SharedFunctionalGroupsSequence.Item_1.MRTimingAndRelatedParametersSequence.Item_1.RepetitionTime;
        a.TE = dcm.SharedFunctionalGroupsSequence.Item_1.MREchoSequence.Item_1.EffectiveEchoTime;
        a.SeriesDescription = dcm.SeriesDescription;
        a.PatientName = dcm.PatientID;
        a.radian = 0;
        a.VOIReadoutFOV = dcm.SharedFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.SliceThickness;
        a.VOIThickness = dcm.SharedFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing(1);
        a.VOIPhaseFOV = dcm.SharedFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing(2);
        a.pdSag = dcm.VolumeLocalizationSequence.Item_1.MidSlabPosition(1);
        a.pdCor = dcm.VolumeLocalizationSequence.Item_1.MidSlabPosition(2);
        a.pdTra = dcm.VolumeLocalizationSequence.Item_1.MidSlabPosition(3);
        a.RowVector = [dcm.SharedFunctionalGroupsSequence.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient(1),...
            dcm.SharedFunctionalGroupsSequence.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient(2),...
            dcm.SharedFunctionalGroupsSequence.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient(3)];
        a.ColumnVector = [dcm.SharedFunctionalGroupsSequence.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient(4),...
            dcm.SharedFunctionalGroupsSequence.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient(5),...
            dcm.SharedFunctionalGroupsSequence.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient(6)];
        a.ndSag = dcm.VolumeLocalizationSequence.Item_1.SlabOrientation(1);
        a.ndCor = dcm.VolumeLocalizationSequence.Item_1.SlabOrientation(2);
        a.ndTra = dcm.VolumeLocalizationSequence.Item_1.SlabOrientation(3);
    end

% -------------------------------------------------------------------------
% get rda info
% -------------------------------------------------------------------------

    function a = get_rda_info(File)
        fid = fopen(File);
        
        head_start_text = '>>> Begin of header <<<';
        head_end_text   = '>>> End of header <<<';
        
        tline = fgets(fid);
        
        while (isempty(strfind(tline , head_end_text))) %#ok<*STREMP>
            
            tline = fgets(fid);
            
            if ( isempty(strfind (tline , head_start_text)) +...
                    isempty(strfind (tline , head_end_text )) == 2)
                
                occurence_of_colon = findstr(':',tline); %#ok<*FSTR>
                variable = tline(1:occurence_of_colon-1) ;
                value    = tline(occurence_of_colon+1 : length(tline)) ;
                
                switch variable
                    case {  'PatientName' , 'SeriesDescription' }
                        eval(['a.' , variable , ' = value;']);
                    case { 'TR' , 'TE' }
                        eval(['a.' , variable , ' = str2num(value);']);
                    case {'VOIPositionSag' }
                        a.pdSag = str2num(value); %#ok<*ST2NM>
                    case {'VOIPositionCor' }
                        a.pdCor = str2num(value);
                    case {'VOIPositionTra' }
                        a.pdTra = str2num(value);
                    case {'VOINormalSag' }
                        a.ndSag = str2num(value);
                    case {'VOINormalCor' }
                        a.ndCor = str2num(value);
                    case {'VOINormalTra' }
                        a.ndTra = str2num(value);
                    case {'VOIRotationInPlane' }
                        a.radian = str2num(value);
                    case {'VOIThickness' }
                        a.VOIThickness = str2num(value);
                    case {'VOIPhaseFOV' }
                        a.VOIPhaseFOV = str2num(value);
                    case {'VOIReadoutFOV' }
                        a.VOIReadoutFOV = str2num(value);
                    case {'RowVector[0]' }
                        a.RowVector(1) = str2num(value);
                    case {'RowVector[1]' }
                        a.RowVector(2) = str2num(value);
                    case {'RowVector[2]' }
                        a.RowVector(3) = str2num(value);
                    case {'ColumnVector[0]' }
                        a.ColumnVector(1) = str2num(value);
                    case {'ColumnVector[1]' }
                        a.ColumnVector(2) = str2num(value);
                    case {'ColumnVector[2]' }
                        a.ColumnVector(3) = str2num(value);
                end
%             else
            end
        end
    end

% -------------------------------------------------------------------------
% return the coordinates of the VOI corners
% -------------------------------------------------------------------------

    function [oct] = create_VOI_nifti_siemens(a)
        
        ap=a.VOIReadoutFOV;
        lr=a.VOIThickness;
        cc=a.VOIPhaseFOV;
        aoc=a.pdCor;
        loc=a.pdSag;
        coc=a.pdTra;
        row=a.RowVector;
        col=a.ColumnVector;
        normal=[a.ndSag a.ndCor a.ndTra];
        
        % The matrix "oct" is created and consists of the coordinates of the VOI-
        % corners as defined around the origin of the coordinate system.
        oct_1 = (lr/2)*(row) + (ap/2)*(col) + (cc/2)*(normal);
        oct_2 = (lr/2)*(-row) + (ap/2)*(col) + (cc/2)*(normal);
        oct_3 = (lr/2)*(-row) + (ap/2)*(-col) + (cc/2)*(normal);
        oct_4 = (lr/2)*(row) + (ap/2)*(-col) + (cc/2)*(normal);
        oct_5 = (lr/2)*(row) + (ap/2)*(col) + (cc/2)*(-normal);
        oct_6 = (lr/2)*(-row) + (ap/2)*(col) + (cc/2)*(-normal);
        oct_7 = (lr/2)*(-row) + (ap/2)*(-col) + (cc/2)*(-normal);
        oct_8 = (lr/2)*(row) + (ap/2)*(-col) + (cc/2)*(-normal);
        
        oct1=[oct_1 1; oct_2 1; oct_3 1; oct_4 1; oct_5 1; oct_6 1; oct_7 1; oct_8 1];
        
        At = [1 0 0 loc;0 1 0 aoc;0 0 1 coc;0 0 0 1];
        f=At*oct1';
        
        e = f';
        oct{1} = [-e(1,1) -e(1,2) e(1,3)];
        oct{2} = [-e(2,1) -e(2,2) e(2,3)];
        oct{3} = [-e(3,1) -e(3,2) e(3,3)];
        oct{4} = [-e(4,1) -e(4,2) e(4,3)];
        oct{5} = [-e(5,1) -e(5,2) e(5,3)];
        oct{6} = [-e(6,1) -e(6,2) e(6,3)];
        oct{7} = [-e(7,1) -e(7,2) e(7,3)];
        oct{8} = [-e(8,1) -e(8,2) e(8,3)];
        
    end

% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------

    function [co] = intersection_nifti_siemens(b,p1,p2,p3,oct)
        
        co_1 = zeros(b,3);
        for k = 1:b
            
            X_1 = [p1(1,k) p1(2,k) p1(3,k)];
            X_2 = [p2(1,k) p2(2,k) p2(3,k)];
            X_3 = [p3(1,k) p3(2,k) p3(3,k)];
            X_4 = oct{1};
            X_5 = oct{2};
            
            A = [1 1 1 1;X_1(1) X_2(1) X_3(1) X_4(1);...
                X_1(2) X_2(2) X_3(2) X_4(2);X_1(3) X_2(3) X_3(3) X_4(3)];
            B = [1 1 1 0;X_1(1) X_2(1) X_3(1) X_5(1)-X_4(1);...
                X_1(2) X_2(2) X_3(2) X_5(2)-X_4(2);...
                X_1(3) X_2(3) X_3(3) X_5(3)-X_4(3)];
            t = - det(A)/det(B);
            
            if t>=0 && t<=1
                co_1(k,:) = [X_4(1)+(X_5(1)-X_4(1))*t X_4(2)+(X_5(2)-X_4(2))...
                    *t X_4(3)+(X_5(3)-X_4(3))*t];
            end
        end
        
        co_2 = zeros(b,3);
        for k = 1:b
            X_1 = [p1(1,k) p1(2,k) p1(3,k)];
            X_2 = [p2(1,k) p2(2,k) p2(3,k)];
            X_3 = [p3(1,k) p3(2,k) p3(3,k)];
            X_4 = oct{1};
            X_5 = oct{4};
            
            A = [1 1 1 1;X_1(1) X_2(1) X_3(1) X_4(1);...
                X_1(2) X_2(2) X_3(2) X_4(2);X_1(3) X_2(3) X_3(3) X_4(3)];
            B = [1 1 1 0;X_1(1) X_2(1) X_3(1) X_5(1)-X_4(1);...
                X_1(2) X_2(2) X_3(2) X_5(2)-X_4(2);...
                X_1(3) X_2(3) X_3(3) X_5(3)-X_4(3)];
            t = - det(A)/det(B);
            
            if t>=0 && t<=1
                co_2(k,:) = [X_4(1)+(X_5(1)-X_4(1))*t X_4(2)+(X_5(2)-X_4(2))...
                    *t X_4(3)+(X_5(3)-X_4(3))*t];
            end
        end
        
        co_3 = zeros(b,3);
        for k = 1:b
            X_1 = [p1(1,k) p1(2,k) p1(3,k)];
            X_2 = [p2(1,k) p2(2,k) p2(3,k)];
            X_3 = [p3(1,k) p3(2,k) p3(3,k)];
            X_4 = oct{1};
            X_5 = oct{5};
            
            A = [1 1 1 1;X_1(1) X_2(1) X_3(1) X_4(1);...
                X_1(2) X_2(2) X_3(2) X_4(2);X_1(3) X_2(3) X_3(3) X_4(3)];
            B = [1 1 1 0;X_1(1) X_2(1) X_3(1) X_5(1)-X_4(1);...
                X_1(2) X_2(2) X_3(2) X_5(2)-X_4(2);...
                X_1(3) X_2(3) X_3(3) X_5(3)-X_4(3)];
            t = - det(A)/det(B);
            
            if t>=0 && t<=1
                co_3(k,:) = [X_4(1)+(X_5(1)-X_4(1))*t X_4(2)+(X_5(2)-X_4(2))...
                    *t X_4(3)+(X_5(3)-X_4(3))*t];
            end
        end
        
        co_4 = zeros(b,3);
        for k = 1:b
            X_1 = [p1(1,k) p1(2,k) p1(3,k)];
            X_2 = [p2(1,k) p2(2,k) p2(3,k)];
            X_3 = [p3(1,k) p3(2,k) p3(3,k)];
            X_4 = oct{2};
            X_5 = oct{3};
            
            A = [1 1 1 1;X_1(1) X_2(1) X_3(1) X_4(1);...
                X_1(2) X_2(2) X_3(2) X_4(2);X_1(3) X_2(3) X_3(3) X_4(3)];
            B = [1 1 1 0;X_1(1) X_2(1) X_3(1) X_5(1)-X_4(1);...
                X_1(2) X_2(2) X_3(2) X_5(2)-X_4(2);...
                X_1(3) X_2(3) X_3(3) X_5(3)-X_4(3)];
            t = - det(A)/det(B);
            
            if t>=0 && t<=1
                co_4(k,:) = [X_4(1)+(X_5(1)-X_4(1))*t X_4(2)+(X_5(2)-X_4(2))...
                    *t X_4(3)+(X_5(3)-X_4(3))*t];
            end
        end
        
        co_5 = zeros(b,3);
        for k = 1:b
            X_1 = [p1(1,k) p1(2,k) p1(3,k)];
            X_2 = [p2(1,k) p2(2,k) p2(3,k)];
            X_3 = [p3(1,k) p3(2,k) p3(3,k)];
            X_4 = oct{2};
            X_5 = oct{6};
            
            A = [1 1 1 1;X_1(1) X_2(1) X_3(1) X_4(1);...
                X_1(2) X_2(2) X_3(2) X_4(2);X_1(3) X_2(3) X_3(3) X_4(3)];
            B = [1 1 1 0;X_1(1) X_2(1) X_3(1) X_5(1)-X_4(1);...
                X_1(2) X_2(2) X_3(2) X_5(2)-X_4(2);...
                X_1(3) X_2(3) X_3(3) X_5(3)-X_4(3)];
            t = - det(A)/det(B);
            
            if t>=0 && t<=1
                co_5(k,:) = [X_4(1)+(X_5(1)-X_4(1))*t X_4(2)+(X_5(2)-X_4(2))...
                    *t X_4(3)+(X_5(3)-X_4(3))*t];
            end
        end
        
        co_6 = zeros(b,3);
        for k = 1:b
            X_1 = [p1(1,k) p1(2,k) p1(3,k)];
            X_2 = [p2(1,k) p2(2,k) p2(3,k)];
            X_3 = [p3(1,k) p3(2,k) p3(3,k)];
            X_4 = oct{5};
            X_5 = oct{6};
            
            A = [1 1 1 1;X_1(1) X_2(1) X_3(1) X_4(1);...
                X_1(2) X_2(2) X_3(2) X_4(2);X_1(3) X_2(3) X_3(3) X_4(3)];
            B = [1 1 1 0;X_1(1) X_2(1) X_3(1) X_5(1)-X_4(1);...
                X_1(2) X_2(2) X_3(2) X_5(2)-X_4(2);...
                X_1(3) X_2(3) X_3(3) X_5(3)-X_4(3)];
            t = - det(A)/det(B);
            
            if t>=0 && t<=1
                co_6(k,:) = [X_4(1)+(X_5(1)-X_4(1))*t X_4(2)+(X_5(2)-X_4(2))...
                    *t X_4(3)+(X_5(3)-X_4(3))*t];
            end
        end
        
        co_7 = zeros(b,3);
        for k = 1:b
            X_1 = [p1(1,k) p1(2,k) p1(3,k)];
            X_2 = [p2(1,k) p2(2,k) p2(3,k)];
            X_3 = [p3(1,k) p3(2,k) p3(3,k)];
            X_4 = oct{5};
            X_5 = oct{8};
            
            A = [1 1 1 1;X_1(1) X_2(1) X_3(1) X_4(1);...
                X_1(2) X_2(2) X_3(2) X_4(2);X_1(3) X_2(3) X_3(3) X_4(3)];
            B = [1 1 1 0;X_1(1) X_2(1) X_3(1) X_5(1)-X_4(1);...
                X_1(2) X_2(2) X_3(2) X_5(2)-X_4(2);...
                X_1(3) X_2(3) X_3(3) X_5(3)-X_4(3)];
            t = - det(A)/det(B);
            
            if t>=0 && t<=1
                co_7(k,:) = [X_4(1)+(X_5(1)-X_4(1))*t X_4(2)+(X_5(2)-X_4(2))...
                    *t X_4(3)+(X_5(3)-X_4(3))*t];
            end
        end
        
        co_8 = zeros(b,3);
        for k = 1:b
            X_1 = [p1(1,k) p1(2,k) p1(3,k)];
            X_2 = [p2(1,k) p2(2,k) p2(3,k)];
            X_3 = [p3(1,k) p3(2,k) p3(3,k)];
            X_4 = oct{3};
            X_5 = oct{4};
            
            A = [1 1 1 1;X_1(1) X_2(1) X_3(1) X_4(1);...
                X_1(2) X_2(2) X_3(2) X_4(2);X_1(3) X_2(3) X_3(3) X_4(3)];
            B = [1 1 1 0;X_1(1) X_2(1) X_3(1) X_5(1)-X_4(1);...
                X_1(2) X_2(2) X_3(2) X_5(2)-X_4(2);...
                X_1(3) X_2(3) X_3(3) X_5(3)-X_4(3)];
            t = - det(A)/det(B);
            
            if t>=0 && t<=1
                co_8(k,:) = [X_4(1)+(X_5(1)-X_4(1))*t X_4(2)+(X_5(2)-X_4(2))...
                    *t X_4(3)+(X_5(3)-X_4(3))*t];
            end
        end
        
        co_9 = zeros(b,3);
        for k = 1:b
            X_1 = [p1(1,k) p1(2,k) p1(3,k)];
            X_2 = [p2(1,k) p2(2,k) p2(3,k)];
            X_3 = [p3(1,k) p3(2,k) p3(3,k)];
            X_4 = oct{3};
            X_5 = oct{7};
            
            A = [1 1 1 1;X_1(1) X_2(1) X_3(1) X_4(1);...
                X_1(2) X_2(2) X_3(2) X_4(2);X_1(3) X_2(3) X_3(3) X_4(3)];
            B = [1 1 1 0;X_1(1) X_2(1) X_3(1) X_5(1)-X_4(1);...
                X_1(2) X_2(2) X_3(2) X_5(2)-X_4(2);...
                X_1(3) X_2(3) X_3(3) X_5(3)-X_4(3)];
            t = - det(A)/det(B);
            
            if t>=0 && t<=1
                co_9(k,:) = [X_4(1)+(X_5(1)-X_4(1))*t X_4(2)+(X_5(2)-X_4(2))...
                    *t X_4(3)+(X_5(3)-X_4(3))*t];
            end
        end
        
        co_10 = zeros(b,3);
        for k = 1:b
            X_1 = [p1(1,k) p1(2,k) p1(3,k)];
            X_2 = [p2(1,k) p2(2,k) p2(3,k)];
            X_3 = [p3(1,k) p3(2,k) p3(3,k)];
            X_4 = oct{4};
            X_5 = oct{8};
            
            A = [1 1 1 1;X_1(1) X_2(1) X_3(1) X_4(1);...
                X_1(2) X_2(2) X_3(2) X_4(2);X_1(3) X_2(3) X_3(3) X_4(3)];
            B = [1 1 1 0;X_1(1) X_2(1) X_3(1) X_5(1)-X_4(1);...
                X_1(2) X_2(2) X_3(2) X_5(2)-X_4(2);...
                X_1(3) X_2(3) X_3(3) X_5(3)-X_4(3)];
            t = - det(A)/det(B);
            
            if t>=0 && t<=1
                co_10(k,:) = [X_4(1)+(X_5(1)-X_4(1))*t X_4(2)+(X_5(2)-X_4(2))...
                    *t X_4(3)+(X_5(3)-X_4(3))*t];
            end
        end
        
        co_11 = zeros(b,3);
        for k = 1:b
            X_1 = [p1(1,k) p1(2,k) p1(3,k)];
            X_2 = [p2(1,k) p2(2,k) p2(3,k)];
            X_3 = [p3(1,k) p3(2,k) p3(3,k)];
            X_4 = oct{7};
            X_5 = oct{8};
            
            A = [1 1 1 1;X_1(1) X_2(1) X_3(1) X_4(1);...
                X_1(2) X_2(2) X_3(2) X_4(2);X_1(3) X_2(3) X_3(3) X_4(3)];
            B = [1 1 1 0;X_1(1) X_2(1) X_3(1) X_5(1)-X_4(1);...
                X_1(2) X_2(2) X_3(2) X_5(2)-X_4(2);...
                X_1(3) X_2(3) X_3(3) X_5(3)-X_4(3)];
            t = - det(A)/det(B);
            
            if t>=0 && t<=1
                co_11(k,:) = [X_4(1)+(X_5(1)-X_4(1))*t X_4(2)+(X_5(2)-X_4(2))...
                    *t X_4(3)+(X_5(3)-X_4(3))*t];
            end
        end
        
        co_12 = zeros(b,3);
        for k = 1:b
            X_1 = [p1(1,k) p1(2,k) p1(3,k)];
            X_2 = [p2(1,k) p2(2,k) p2(3,k)];
            X_3 = [p3(1,k) p3(2,k) p3(3,k)];
            X_4 = oct{6};
            X_5 = oct{7};
            
            A = [1 1 1 1;X_1(1) X_2(1) X_3(1) X_4(1);...
                X_1(2) X_2(2) X_3(2) X_4(2);X_1(3) X_2(3) X_3(3) X_4(3)];
            B = [1 1 1 0;X_1(1) X_2(1) X_3(1) X_5(1)-X_4(1);...
                X_1(2) X_2(2) X_3(2) X_5(2)-X_4(2);...
                X_1(3) X_2(3) X_3(3) X_5(3)-X_4(3)];
            t = - det(A)/det(B);
            
            if t>=0 && t<=1
                co_12(k,:) = [X_4(1)+(X_5(1)-X_4(1))*t X_4(2)+(X_5(2)-X_4(2))...
                    *t X_4(3)+(X_5(3)-X_4(3))*t];
            end
        end
        
        % All twelve matrices are stored in one big matrix called "co" before
        % returned to position.m.
        co = [co_1 co_2 co_3 co_4 co_5 co_6 co_7 co_8 co_9 co_10 co_11 co_12];
        
    end

% -------------------------------------------------------------------------
%   create mask matrix
% -------------------------------------------------------------------------

    function G = transition_nifti(p1,n1,n2,n3,ps_r,ps_c,r,c,co)
        
        P_1 = [p1(1,1) p1(2,1) p1(3,1)];
        P_2 = [p1(1,2) p1(2,2) p1(3,2)];
        P = P_2-P_1;
        
        ss = norm(P);
        
        d = cross(c,r);
        
        A = [ss*d(1) ps_r*r(1) ps_c*c(1);ss*d(2) ps_r*r(2) ps_c*c(2);...
            ss*d(3) ps_r*r(3) ps_c*c(3)];
        
        G = ones(n1,n2,n3);
        
        for s = 1:n3
            for t = 1:3:34
                koord = co(s,t:t+2);
                if norm(koord,2) > 0
                    b = [koord(1)+ss*(c(2)*r(3)-c(3)*r(2))+ps_r*r(1)+ps_c...
                        *c(1)-p1(1,1);koord(2)+ss*(c(3)*r(1)-c(1)*r(3))+ps_r...
                        *r(2)+ps_c*c(2)-p1(2,1);koord(3)+ss*(c(1)*r(2)-c(2)...
                        *r(1))+ps_r*r(3)+ps_c*c(3)-p1(3,1)];
                    %                     losn = inv(A)*b;
                    losn=A\b;
                    k = round(losn(1));
                    i = round(losn(2));
                    j = round(losn(3));
                    
                    if k < 0
                        k = -k;
                    end
                    if j < 0
                        j = -j;
                    end
                    if i < 0
                        i = -i;
                    end
                    G(i,j,k) = 65000;
                end
            end
        end
    end

% -------------------------------------------------------------------------
%   create mask
% -------------------------------------------------------------------------

    function mask = create_filter_nifti(G,n1,n2,n3)
        
        n1 = double(n1);
        n2 = double(n2);
        
        for k = 1:n3
            %             disp(max(max(G(:,:,k))))
            if max(max(G(:,:,k))) == 1
                G(:,:,k) = 0;
            else
                [i, j] = find(G(:,:,k)>1);
                
                if length(i)>2
                    
                    K = convhull(i,j);
                    if length(K) == 4
                        i = [i(K(1));i(K(2));i(K(3))];
                        j = [j(K(1));j(K(2));j(K(3))];
                    elseif length(K) == 5
                        i = [i(K(1));i(K(2));i(K(3));i(K(4))];
                        j = [j(K(1));j(K(2));j(K(3));j(K(4))];
                    elseif length(K) == 6
                        i = [i(K(1));i(K(2));i(K(3));i(K(4));i(K(5))];
                        j = [j(K(1));j(K(2));j(K(3));j(K(4));j(K(5))];
                    elseif length(K) == 7
                        i = [i(K(1));i(K(2));i(K(3));i(K(4));i(K(5));i(K(6))];
                        j = [j(K(1));j(K(2));j(K(3));j(K(4));j(K(5));j(K(6))];
                    elseif length(K) == 8
                        i = [i(K(1));i(K(2));i(K(3));i(K(4));i(K(5));...
                            i(K(6));i(K(7))];
                        j = [j(K(1));j(K(2));j(K(3));j(K(4));j(K(5));...
                            j(K(6));j(K(7))];
                    elseif length(K) == 9
                        i = [i(K(1));i(K(2));i(K(3));i(K(4));i(K(5));i(K(6));...
                            i(K(7));i(K(8));];
                        j = [j(K(1));j(K(2));j(K(3));j(K(4));j(K(5));j(K(6));...
                            j(K(7));j(K(8));];
                    end
                else
                    G(:,:,k) = 0;
                end
                
                bw = poly2mask(j,i,n1,n2);
                G(:,:,k) = bw;
                
                
            end
        end
        
        mask = G;
    end

end
