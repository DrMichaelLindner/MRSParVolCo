# MRSParVolCo - MRS Partial Volume Correction


MRSParVolCo is designed to correct estimated metabolite amplitudes of single voxel MRS scans acquired at SIEMENS 3T MR scanners. 
MRSParVolCo corrects the estimated/quantified amplitudes of the following brain metabolites depending on the different tissue concentrations.
- N-Acetylaspartic acid (NAA)
- Creatine (CRE)
- Choline (CHO)
- Myo-inositol (INS)
- Glutamate (GLU)
- Gamma-Aminobutyric acid (GABA)


MRSParVolCo uses the following values for tissue and metabolite relaxation time in 3T found in the literature:

- Tissue relaxation times in 3T:
  - GM_T1 = 1500;
  - GM_T2 = 63;
  - WM_T1 = 1000;
  - WM_T2 = 50;
  - CSF_T1 = 4000;
  - CSF_T2 = 200;

- Metabolite relaxation times in 3T:
  - GM_T1_NAA = 1403;
  - GM_T2_NAA = 272;
  - GM_T1_CHO = 1182;
  - GM_T2_CHO = 217;
  - GM_T1_CRE = 1320;
  - GM_T2_CRE = 146;
  - GM_T1_GLU = 1220;
  - GM_T2_GLU = 185;
  - GM_T1_INS = 1100;
  - GM_T2_INS = 200;
  - GM_T1_GABA = 1310;
  - GM_T2_GABA = 88;

  - Water visibility correction:
  - WVC_GM = 0.779407794;
  - WVC_WM = 0.645846458;
  - WVC_CSF = 0.97;


These values can be found and changed in MRSParVolCo in lines 77ff and in MRSParVolCo_UI in line 157ff.


## *Install*  
Copy the MRSParVolCo folder in a folder of your choice on your system and add the directory to your matlab path.


## *Dependencies*  
- For the usage of MRSParVolCo you need to have Matlab (www.mathworks.com) installed. MRSParVolCo is developed on MatLab2017b but it should work with other Matlab distributions later than 2015b as well.
- For loading nifti files into Matlab MRSParVolCo needs the "Tools for NIfTI and ANALYZE image" toolbox (https://uk.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) in the Matlab path.


## *Usage*  
Two version of MRSParVolCo are available: a scriptable version (MRSParVolCo) and a version using user interfaces for inputs (MRSParVolCo_UI).

The estimated/quantified amplitudes of the metabolites (e.g. from LCModal, jMRUI, Tarquin, etc) are corrected based on the tissue concentration in the acquired single MRS voxel. 



### MRSParVolCo:  
A scriptable version using input vectors for the metabolite amplitude values and percentages of tissue concentration .

#### USAGE
MRSParVolCo(GM, WM, CSF, TR, TE, ouputfolder, NAA, CRE, CHO ,GLU, INS, GABA)

#### INPUT
- vector (size nx1) of amplitudes of the metabolites NAA, CRE, CHO ,GLU, INS and/or GABA. Use [] for the metabolites you do not want to correct.  For using m different amplitudes of the same metabolite (e.g. from different estimation approaches) as input you can use a matrix (size nxm)
- all tissue composition percent vectors for GM, WM and CSF (e.g. from MRSGetTissueComp) (size nx1)
- TE
- TR
- outputfolder


!!!!!!!!! metabolite vectors and tissue composition vectors needs to have the same length n and the rows need to correspond (e.g. subjects s: metabolite value in row s as well as tissue composition value is row s)

### MRSParVolCo_UI:  
A version using simple input dialogues for manually selecting text files.

#### USAGE
MRSParVolCo_UI

#### INPUT 
The input can be selected via User Interfaces
- text or excel files including columns of values of the metabolites NAA, CRE, CHO ,GLU, INS and/or GABA. Text files can contain comma separated columns for different amplitudes of the same metabolite (e.g. from different estimation approaches)
- text file of tissue composition values (e.g. from MRSGetTissueComp)
- TE
- TR
- outputfolder


!!!!!!!!! number of rows in metabolite file and tissue composition value file need to be the same and the rows need to correspond (e.g. subjects s: metabolite value in row s as well as tissue composition value is row s)


### OUTPUT of MRSParVolCo and MRSParVolCo_UI
Both functions produce text files as output containing input amplitudes and corrected amplitudes of each given metabolite.


### Getting tissue composition for the input of MRSParVolCo:
Tissue composition represents the percentage of gray and white matter and cerebrospinal fluid in the acquired single voxel MRS. For getting the tissue composition  you can use the function MRSGetTissueComp before using MRSParVolCo.
The tissue composition is calculated based on a structural nifti files and a SIEMENS dicom or rda MRS file. (tested on MRS data from  SIEMENS single voxel PRESS and MEGA-PRESS sequences from SIEMENS TRIO and PRISMAfit)


Again a scriptable as well as a UI version are available:

### MRSGetTissueComp:
A scriptable version for looping over multiple subjects and MRS files.

#### USAGE
MRSGetTissueComp(batch)

#### INPUT
Matlab structure of the following format:
- batch(x).t1file = filename with path of nifti file of structrual scan after brain extraction (bet)
- batch(x).GMfile = filename with path of nifti file of grey matter mask of segmentation
- batch(x).WMfile = filename with path of nifti file of white matter mask of segmentation
- batch(x).CSFfile = filename with path of nifti file of CSF mask of segmentation
- batch(x).MRSfiles{y} = SIEMENS rda or dicom MRS file
- batch(x).outpath = output folder
MRSGetTissueComp can loop over multiple subjects (x) and multiple MRS files per subject (y).


### MRSGetTissueComp_UI:
A version using simple input dialogues for a single file.

#### USAGE
MRSGetTissueComp_UI

#### INPUT
The input can be selected via User Interfaces
- nifti file of structrual scan after brain extraction (bet)
- nifti file of grey matter mask of segmentation
- nifti file of white matter mask of segmentation
- nifti file of CSF mask of segmentation
- SIEMENS rda or dicom MRS file
- output folder


### OUTPUT of MRSGetTissueComp and MRSGetTissueComp_UI
Both functions produce the following outputs:

- A Mask of the MRS single voxel 
- Text file containing the center coordinates of the MRS voxel (in single subject space)
- Text file containing the values:
  - GM_per: percentage of gray matter in the MRS voxel
  - WM_per: percentage of white matter in the MRS voxel
  - CSF_per: percentage of CSF in the MRS voxel
  - GM_frac: Fraction of gray matter in the tissue of the MRS voxel
  - WM_frac: Fraction of white matter in the tissue of the MRS voxel
- Matlab file containing the following vectors (They can be uses as input for MRSParVolCo_2_correct):
  - GMper: percentage of gray matter 
  - WMper: percentage of white matter
  - CSFper: percentage of cerebrospinal fluid 
  - GMfrac: fraction of grey matter
  - WMfrac: fraction of white matter 
  

## *Licence*  
MRSParVolCo is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License (GPLv3) as published
by the Free Software Foundation;
  
  
## *Author*
Michael Lindner  
University of Reading, 2018  
School of Psychology and Clinical Language Sciences  
Centre for Integrative Neuroscience and Neurodynamics
