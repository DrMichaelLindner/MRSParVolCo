function MRSParVolCo_UI

% MRSParVolCo_UI corrects the given MRS metabolite amplitudes depending on 
% the given tissue composition values of the voxel.
%
% INPUT (via GUI) 
% 	- text or excel files including values of 
%		the metabolites NAA, CRE, CHO ,GLU, INS and/or GABA
%       Text files can contain comma separated columns for different 
%       amplitudes of the same metabolite (e.g. from different estimation 
%       approaches)
%   - text file of tissue composition values (e.g. from MRSGetTissueComp) 
%       with columns in the order %GM, %WM and %CSF
% !!!!!!!!! number of rows in metabolite file and tissue composition value 
% file need to be the same and the rows need to correspond (e.g. 
% subjects s: metabolite value in row s as well as tissue composition value
% is row s)
%	- TE
%	- TR
%	- outputprefix
%
% OUTPUT
%   .txt textfile with input amplitudes and corrected amplitudes
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

% get files
[NAAFile,NAAPath]=uigetfile({'*.txt';'*.xls';'*.xlsx'},'Select NAA matrix file');
[CREFile,CREPath]=uigetfile({'*.txt';'*.xls';'*.xlsx'},'Select CRE matrix file',NAAPath);
[CHOFile,CHOPath]=uigetfile({'*.txt';'*.xls';'*.xlsx'},'Select CHO matrix file',NAAPath);
[GLUFile,GLUPath]=uigetfile({'*.txt';'*.xls';'*.xlsx'},'Select GLU matrix file',NAAPath);
[INSFile,INSPath]=uigetfile({'*.txt';'*.xls';'*.xlsx'},'Select INS matrix file',NAAPath);
[GABAFile,GABAPath]=uigetfile({'*.txt';'*.xls';'*.xlsx'},'Select GABA matrix file',NAAPath);

[GMFile,GMPath]=uigetfile('*.txt','Tissue composition file',NAAPath);

prompt = {'TR:','TE:','Output filename prefix:'};
dlg_title = 'Input parameter';
num_lines = 1;
def = {'1500','68','Metabolites_corrected'};
answer = inputdlg(prompt,dlg_title,num_lines,def);

outputpath = uigetdir(GMPath,'Select output folder');

savename = answer{3};

TR=str2double(answer{1});
TE=str2double(answer{2});


%load metabolite files
count = 0;
if ischar(NAAFile)
    [Pathstr,name,ext] = fileparts(NAAFile); %#ok<*ASGLU>
    if strcmp(ext,'.txt')
        NAA=importdata([NAAPath,NAAFile]);
    elseif strcmp(ext,'.xls')
        NAA=xlsread([NAAPath,NAAFile]);
    end
    count = count + 1;
    lenvec(count)=size(NAA,1);
    savename = [savename, '_NAA'];
end

if ischar(CREFile)
    [Pathstr,name,ext] = fileparts(CREFile);
    if strcmp(ext,'.txt')
        CRE=importdata([CREPath,CREFile]);
    elseif strcmp(ext,'.xls')
        CRE=xlsread([CREPath,CREFile]);
    end
    count = count + 1;
    lenvec(count)=size(CRE,1);
    savename = [savename, '_CRE'];
end

if ischar(CHOFile)
    [Pathstr,name,ext] = fileparts(CHOFile);
    if strcmp(ext,'.txt')
        CHO=importdata([CHOPath,CHOFile]);
    elseif strcmp(ext,'.xls')
        CHO=xlsread([CHOPath,CHOFile]);
    end
    count = count + 1;
    lenvec(count)=size(CHO,1);
    savename = [savename, '_CHO'];
end

if ischar(GLUFile)
    [Pathstr,name,ext] = fileparts(GLUFile);
    if strcmp(ext,'.txt')
        GLU=importdata([GLUPath,GLUFile]);
    elseif strcmp(ext,'.xls')
        GLU=xlsread([GLUPath,GLUFile]);
    end
    count = count + 1;
    lenvec(count)=size(GLU,1);
    savename = [savename, '_GLU'];
end

if ischar(INSFile)
    [Pathstr,name,ext] = fileparts(INSFile);
    if strcmp(ext,'.txt')
        INS=importdata([INSPath,INSFile]);
    elseif strcmp(ext,'.xls')
        INS=xlsread([INSPath,INSFile]);
    end
    count = count + 1;
    lenvec(count)=size(INS,1);
    savename = [savename, '_INS'];
end


if ischar(GABAFile)
    [Pathstr,name,ext] = fileparts(GABAFile);
    if strcmp(ext,'.txt')
        GABA=importdata([GABAPath,GABAFile]);
    elseif strcmp(ext,'.xls')
        GABA=xlsread([GABAPath,GABAFile]);
    end
    count = count + 1;
    lenvec(count)=size(GABA,1);
    savename = [savename, '_GABA'];
end

if mean(lenvec)~=min(lenvec)
    error('ERROR: vectors have different lengths')
end



% load GM file

GWC=importdata([GMPath,GMFile]);

GM=GWC(:,1);
WM=GWC(:,2);
CSF=GWC(:,3);


if mean(lenvec)<length(GM)
    error('ERROR: more tissue files selected than values in vectors')
elseif mean(lenvec)>length(GM)
    error('ERROR: less tissue files selected than values in vectors')
end


% Tissue and metabolite relaxation time in 3T:

% Tissue relaxation times in 3T:
GM_T1 = 1500;
GM_T2 = 63;
WM_T1 = 1000;
WM_T2 = 50;
CSF_T1 = 4000;
CSF_T2 = 200;

% Metabolite relaxation times in 3T
GM_T1_NAA = 1403;
GM_T2_NAA = 272;
GM_T1_CHO = 1182;
GM_T2_CHO = 217;
GM_T1_CRE = 1320;
GM_T2_CRE = 146;
GM_T1_GLU = 1220;
GM_T2_GLU = 185;
GM_T1_INS = 1100;
GM_T2_INS = 200;
GM_T1_GABA = 1310;
GM_T2_GABA = 88;


% Water visibility correction
WVC_GM = 0.779407794;
WVC_WM = 0.645846458;
WVC_CSF = 0.97;


R_GM = exp(-TE/GM_T2)*(1-exp(-TR/GM_T1));
R_WM = exp(-TE/WM_T2)*(1-exp(-TR/WM_T1));
R_CSF = exp(-TE/CSF_T2)*(1-exp(-TR/CSF_T1));

for gg = 1:length(GM)
    ATT_H2O(gg) = (GM(gg)*WVC_GM*R_GM + WM(gg)*WVC_WM*R_WM + CSF(gg)*WVC_CSF*R_CSF) / (GM(gg)*WVC_GM + WM(gg)*WVC_WM + CSF(gg)*WVC_CSF); %#ok<*AGROW>
    XX(gg) = 1 - CSF(gg)*WVC_CSF/(GM(gg)*WVC_GM + WM(gg)*WVC_WM + CSF(gg)*WVC_CSF);
end

% create output text file
cd(outputpath)
txt_tiss_filename=[savename, '.txt'];
fileID=fopen(txt_tiss_filename,'wt');



fprintf(fileID,'Input values for:\n' );
fprintf('Input values for:\n' );
if ischar(NAAFile)
    fprintf(fileID,'NAA = \n');
    fprintf('NAA = \n');
    for gg = 1:length(GM)
        for ii = 1:size(NAA,2)
            fprintf(fileID,'%1.4f\t',NAA(gg,ii));
            fprintf('%1.4f\t',NAA(gg,ii));
        end
        fprintf(fileID,'\n');
        fprintf('\n');
    end
    Metabolites_orig.NAA = NAA;
    fprintf('\n');
end
if ischar(CHOFile)
%     CHO = str2num(CHOc{1});
    fprintf(fileID,'CHO = \n');
    fprintf('CHO = \n');
    for gg = 1:length(GM)
        for ii = 1:size(CHO,2)
            fprintf(fileID,'%1.4f\t',CHO(gg,ii));
            fprintf('%1.4f\t',CHO(gg,ii));
        end
        fprintf(fileID,'\n');
        fprintf('\n');
    end
    Metabolites_orig.CHO = CHO;
    fprintf('\n');
end
if ischar(CREFile)
    fprintf(fileID,'CRE = \n');
    fprintf('CRE = \n');
    for gg = 1:length(GM)
        for ii = 1:size(CRE,2)
            fprintf(fileID,'%1.4f\t',CRE(gg,ii));
            fprintf('%1.4f\t',CRE(gg,ii));
        end
        fprintf(fileID,'\n');
        fprintf('\n');
    end
    Metabolites_orig.CRE = CRE;
    fprintf('\n');
end
if ischar(GLUFile)
    fprintf(fileID,'GLU = \n');
    fprintf('GLU = \n');
    for gg = 1:length(GM)
        for ii = 1:size(GLU,2)
            fprintf(fileID,'%1.4f\t',GLU(gg,ii));
            fprintf('%1.4f\t',GLU(gg,ii));
        end
        fprintf(fileID,'\n');
        fprintf('\n');
    end
    Metabolites_orig.GLU = GLU;
    fprintf('\n');
end
if ischar(INSFile)
    fprintf(fileID,'INS = \n');
    fprintf('INS = \n');
    for gg = 1:length(GM)
        for ii = 1:size(INS,2)
            fprintf(fileID,'%1.4f\t',INS(gg,ii));
            fprintf('%1.4f\t',INS(gg,ii));
        end
        fprintf(fileID,'\n');
        fprintf('\n');
    end
    Metabolites_orig.INS = INS; %#ok<*STRNU>
    fprintf('\n');
end

fprintf(fileID,'\n\nVoxel tissue values (in percent):\n' );
fprintf(fileID,'GM\tWMr\tCSFr\n');
for gg = 1:length(GM)
    fprintf(fileID,'%1.2f\t',GM(gg));
    fprintf(fileID,'%1.2f\t',WM(gg));
    fprintf(fileID,'%1.2f\n',CSF(gg));
end

fprintf('\n\nVoxel tissue values (in percent):\n' );
fprintf('GM\tWM\tCSF\n');
for gg = 1:length(GM)
    fprintf('%1.2f\t',GM(gg));
    fprintf('%1.2f\t',WM(gg));
    fprintf('%1.2f\n',CSF(gg));
end

fprintf(fileID,'\n\nCorrected values for:\n\n' );
fprintf('\n\nCorrected values for:\n\n' );
if ischar(NAAFile)
    fprintf(fileID,'NAA\n');
    fprintf('NAA\n');
    R_NAA = exp(-TE/GM_T2_NAA)*(1-exp(-TR/GM_T1_NAA));
    for gg = 1:length(GM)
        ATT_NAA(gg) = R_NAA*XX(gg);
        for ii = 1:size(NAA,2)
            COR_NAA(gg,ii) = NAA(gg,ii) * ATT_H2O(gg)/ATT_NAA(gg);
            fprintf(fileID,'%1.4f\t',COR_NAA(gg,ii));
            fprintf('%1.4f\t',COR_NAA(gg,ii));
        end
        fprintf(fileID,'\n');
        fprintf('\n');
    end
    Metabolites_corrected.NAA = COR_NAA;
    fprintf('\n');
end

if ischar(CHOFile)
    fprintf(fileID,'CHO\n');
    fprintf('CHO\n');
    R_CHO = exp(-TE/GM_T2_CHO)*(1-exp(-TR/GM_T1_CHO));
    for gg = 1:length(GM)
        ATT_CHO(gg) = R_CHO*XX(gg);
        for ii = 1:size(CHO,2)
            COR_CHO(gg,ii) = CHO(gg,ii) * ATT_H2O(gg)/ATT_CHO(gg);
            fprintf(fileID,'%1.4f\t',COR_CHO(gg,ii));
            fprintf('%1.4f\t',COR_CHO(gg,ii));
        end
        fprintf(fileID,'\n');
        fprintf('\n');
    end
    Metabolites_corrected.CHO = COR_CHO;
    fprintf('\n');
end

if ischar(CREFile)
    fprintf(fileID,'CRE\n');
    fprintf('CRE\n');
    R_CRE = exp(-TE/GM_T2_CRE)*(1-exp(-TR/GM_T1_CRE));
    for gg = 1:length(GM)
        ATT_CRE(gg) = R_CRE*XX(gg);
        for ii = 1:size(CRE,2)
            COR_CRE(gg,ii) = CRE(gg,ii) * ATT_H2O(gg)/ATT_CRE(gg);
            fprintf(fileID,'%1.4f\t',COR_CRE(gg,ii));
            fprintf('%1.4f\t',COR_CRE(gg,ii));
        end
        fprintf(fileID,'\n');
        fprintf('\n');
    end
    Metabolites_corrected.CRE = COR_CRE;
    fprintf('\n');
end

if ischar(GLUFile)
    fprintf(fileID,'GLU\n');
    fprintf('GLU\n');
    R_GLU = exp(-TE/GM_T2_GLU)*(1-exp(-TR/GM_T1_GLU));
    for gg = 1:length(GM)
        ATT_GLU(gg) = R_GLU*XX(gg);
        for ii = 1:size(GLU,2)
            COR_GLU(gg,ii) = GLU(gg,ii) * ATT_H2O(gg)/ATT_GLU(gg);
            fprintf(fileID,'%1.4f\t',COR_GLU(gg,ii));
            fprintf('%1.4f\t',COR_GLU(gg,ii));
        end
        fprintf(fileID,'\n');
        fprintf('\n');
    end
    Metabolites_corrected.GLU = COR_GLU;
    fprintf('\n');
end

if ischar(INSFile)
    fprintf(fileID,'INS\n');
    fprintf('INS\n');
    R_INS = exp(-TE/GM_T2_INS)*(1-exp(-TR/GM_T1_INS));
    for gg = 1:length(GM)
        ATT_INS(gg) = R_INS*XX(gg);
        for ii = 1:size(INS,2)
            COR_INS(gg,ii) = INS(gg,ii) * ATT_H2O(gg)/ATT_INS(gg);
            fprintf(fileID,'%1.4f\t',COR_INS(gg,ii));
            fprintf('%1.4f\t',COR_INS(gg,ii));
        end
        fprintf(fileID,'\n');
        fprintf('\n');
    end
    Metabolites_corrected.INS = COR_INS;
end


if ischar(GABAFile)
    fprintf(fileID,'GABA\n');
    fprintf('GABA\n');
    R_GABA = exp(-TE/GM_T2_GABA)*(1-exp(-TR/GM_T1_GABA));
    for gg = 1:length(GM)
        ATT_GABA(gg) = R_GABA*XX(gg);
        for ii = 1:size(GABA,2)
            COR_GABA(gg,ii) = GABA(gg,ii) * ATT_H2O(gg)/ATT_GABA(gg);
            fprintf(fileID,'%1.4f\t',COR_GABA(gg,ii));
            fprintf('%1.4f\t',COR_GABA(gg,ii));
        end
        fprintf(fileID,'\n');
        fprintf('\n');
    end
    Metabolites_corrected.INS = COR_GABA;
end


fclose(fileID);

cd(outputpath)
savename = [savename, '.mat'];
save(savename, 'Metabolites_orig','Metabolites_corrected')





