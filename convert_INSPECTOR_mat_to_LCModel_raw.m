function convert_INSPECTOR_mat_to_LCModel_raw()
% 
% Converts a folder of INSPECTOR-format spectral .mat files to .raw files
% readable by LCModel 
%
% Inputs: Function run prompts user input of folder containing only 
% INSPECTOR-format .mat files to be converted to LCModel-readable .raw
% files.
% 
% Outputs: LCModel-readable .raw files containing same time-domain free induction decay (FID)
% data as original INSPECTOR-format .mat files 
% 
% Author: Kelley Swanberg (Columbia University, 2023) 
% 
% Acknowledgements: create_fmt() function adapted from materials 
% provided by MathWorks Support team available at 
% https://www.mathworks.com/matlabcentral/answers/464993-how-can-i-write-a-matrix-to-a-text-file-with-special-delimiters
% 
% Locate correct folder; must contain only INSPECTOR-readable .mat files to convert 
current_folder = uigetdir(); 

% Find .mat files to convert
list_of_cases_struct = dir(current_folder);
list_of_cases_struct_clean = list_of_cases_struct(~ismember({list_of_cases_struct.name},{'.','..'}));
list_of_cases = {list_of_cases_struct_clean.name}'; 

% Determine appropriate size of table to build
num_cases = length(list_of_cases);

for i = 1:num_cases
% Convert directory of INSPECTOR-format MAT files to LCModel RAW files 
mat_filename = list_of_cases{i}; 
data_mat = load(mat_filename);

created_date =  datetime('today'); 
raw_ID = strrep(mat_filename, '.mat', ''); 
raw_volume = '1.00000E+00'
raw_tramp = '1.00000E+00'

% Convert complex FID data into real-valued matrix of shape readable by LCModel 
raw_fid = data_mat.exptDat.fid; 
length_fid_complex = data_mat.exptDat.nspecC; 
length_fid_values = 2 * length_fid_complex; 
num_rows_raw = length_fid_values/8; 

raw_fid_real_imag_sep = [real(raw_fid) -imag(raw_fid)]; % Negative imaginary component because otherwise they will end up upside-down in LCModel 
raw_fid_real_imag_sep_transpose = raw_fid_real_imag_sep'; 
raw_fid_real_imag_sep_vector = raw_fid_real_imag_sep_transpose(:); 
raw_fid_matrix = reshape(raw_fid_real_imag_sep_vector, 8, num_rows_raw); 

%Create raw file 
raw_filename = strrep(mat_filename, '.mat', '.raw'); 

% Write header information to raw file 
fileID = fopen(raw_filename,'w');
fprintf(fileID,' Created: %s\n', created_date);
fprintf(fileID,' $NMID\n');
fprintf(fileID,' ID=''%s''\n', raw_ID);
fprintf(fileID,' FMTDAT=''(8E13.5)''\n');
fprintf(fileID,' VOLUME=   %s\n', raw_volume);
fprintf(fileID,' TRAMP=   %s\n', raw_tramp);
fprintf(fileID,' $END\n');
fprintf(fileID,' ');

% Write FID data to raw file; adapted from MathWorks support team 
% https://www.mathworks.com/matlabcentral/answers/464993-how-can-i-write-a-matrix-to-a-text-file-with-special-delimiters 
precision = '%+5.5E'; % desired precision for values in A
delimiter = ' ';
line_terminator = '\n ';
format = [create_fmt(precision, delimiter, size(raw_fid_matrix)') line_terminator];
fprintf(fileID, format, raw_fid_matrix);

fclose(fileID);
close all; 

end
end

% Taken from MathWorks support team 
% https://www.mathworks.com/matlabcentral/answers/464993-how-can-i-write-a-matrix-to-a-text-file-with-special-delimiters
function s = create_fmt(prec, dlm, n_fmt)
s = prec;
for i = 1:2*(n_fmt-1)
    if mod(i, 2) == 1
        s = [s dlm];
    else
        s = [s prec];
    end
end
end