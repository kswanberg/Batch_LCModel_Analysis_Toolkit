function extract_fid_from_INSPECTOR_basis( ) 
% Batch LCModel toolkit: Extract individual FIDs from INSPECTOR basis
% Did a colleague give you a .MAT INSPECTOR basis without providing the
% original simulation files? How rude. 
% 
% Dependencies: None
% 
% Recipe: 
% 1. Put name of basis file to be decomposed into line 17 (must be in MATLAB path) 
% 2. GRATTIS now you have a bunch of MAT files 
% 
% K. Swanberg, 2020-22

load('YOUR_FAVORITE_INSPECTOR_BASIS.mat'); % Put name of basis file here 

% Define fundamental aspects of basis functions 
larmor_freq = lcmBasis.sf*1000000;
spectral_width = lcmBasis.sw_h; 
num_basis_functions = length(lcmBasis.data); 

% Walk through basis set to extract all signals one at a time
for i=1:num_basis_functions 
    
% Extract name of basis function 
basis_function_name = lcmBasis.data{i}{1};

% Extract basis function FID 
basis_function_fid = lcmBasis.data{i}{4};

% Insert final spectrum into structure readable by INSPECTOR
exptDat.fid = basis_function_fid;
exptDat.sf = larmor_freq / 1000000; %MHz
exptDat.sw_h = spectral_width; 
exptDat.nspecC = length(lcmBasis.data{1,1}{1,4}); % Extracts FID length from first FID in basis

% Save spectra to individual MAT files 
extracted_fid = sprintf('%s.mat', basis_function_name); % Specify file name for spectrum to be saved 
save(extracted_fid, 'exptDat');
end
end 