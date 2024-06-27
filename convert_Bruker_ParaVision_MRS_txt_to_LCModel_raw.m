function convert_Bruker_ParaVision_MRS_txt_to_LCModel_raw()
%
% Converts text files exported from Bruker ParaVision
% 360.3.3 spectral processing tool into RAW-format files suitable for analysis
% directly in LCModel 
% 
%
% Inputs: Upon function run the user will be prompted to select a root 
% directory containing .TXT files from the Bruker ParaVision spectral processing
% tool marked for conversion to .RAW files
% 
% Note that this function rewinds a frequency shift present in the
% .TXT-based FIDs that may not generalize to all versions of ParaVision. 
% 
% Outputs: LCModel-readable .RAW files for use in spectral quantification 
% 
% Author: Kelley Swanberg (Columbia University, 2023) 
%
% Syntax for defining precision of FID data output borrowed from % https://www.mathworks.com/matlabcentral/answers/
% 464993-how-can-i-write-a-matrix-to-a-text-file-with-special-delimiters
% Written for MATLAB 2018a.
%

% % Locate correct folder; must contain only TXT files to convert 
current_folder = uigetdir(); 

% Find files on which to perform conversion to RAW
list_of_cases_struct = dir(current_folder);
list_of_cases_struct_clean = list_of_cases_struct(~ismember({list_of_cases_struct.name},{'.','..'}));
list_of_cases = {list_of_cases_struct_clean.name}'; 

% Determine appropriate size of case table to build
num_cases = length(list_of_cases);

for i = 1:num_cases
% Convert directory of Bruker-exported TXT files to LCModel RAW files 
txt_filename = list_of_cases{i}; 
fileOpenID = fopen(txt_filename,'r');
formatSpec = '%c'
data_txt = fscanf(fileOpenID,formatSpec); 
data_text_str = splitlines(convertCharsToStrings(data_txt)); 

created_date =  datetime('today'); 
raw_ID = strrep(txt_filename, '.mat', ''); 
raw_volume = '1.00000E+00'
raw_tramp = '1.00000E+00'

% Convert complex FID data into real-valued matrix of shape readable by LCModel
fid_string = data_text_str(11:end); 
fid_num = str2double(fid_string); 
fid_num = fid_num(~isnan(fid_num));
raw_fid = fid_num; 
length_fid_values = length(raw_fid); 
num_rows_raw = length_fid_values/8; 

% Convert spectrum to FID 
raw_fid_complex = flipud(complex(raw_fid(1:2:end), raw_fid(2:2:end))); 
raw_fid_ifft_unshifted = -flipud(ifft(ifftshift(raw_fid_complex))); 

% Manipulate FID for frequency-domain shift by 2 ppm (seems to be needed 
% for preprocessed text files under current parameters; note that this may 
% change for different versions of ParaVision) 
f0_ppm = -2; %ppm; probably specific to acquisition software processing pipeline 
B0 = 300.3; %MHz; obviously 7 T
f0_Hz = f0_ppm * B0; 
dwell_time = 0.003; % s; dependent on receive bandwidth of acquisition protocol 
time_vector_length = length_fid_values/2;
time_vector = [linspace(0, time_vector_length*dwell_time, time_vector_length)]'; 

raw_fid_ifft = exp(j*2*pi*f0_Hz*time_vector).*raw_fid_ifft_unshifted; 

raw_fid_real_imag_sep = [real(raw_fid_ifft) imag(raw_fid_ifft)]; 
raw_fid_real_imag_sep_transpose = raw_fid_real_imag_sep'; 
raw_fid_real_imag_sep_vector = raw_fid_real_imag_sep_transpose(:); 
raw_fid_matrix = reshape(raw_fid_real_imag_sep_vector, 8, num_rows_raw); 

%Create raw file 
raw_filename = strrep(txt_filename, '.txt', '.raw'); 

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