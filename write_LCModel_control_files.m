function write_LCModel_control_files()
% 
% Duplicates an LCModel CONTROL file named template.control over a directory
% of LCModel-readable .raw data files for easy batch analysis of multiple
% LCModel .raw datasets 
%
% Inputs: 
% 
%       1) LCModel CONTROL file, dictating control parameters to be duplicated over the whole analysis, 
%       entitled 'template.control', inside folder named 'CONTROL_Files'. Note
%       that .raw-file-specific intputs for CONTROL file must read as follows: 
%
%           FILRAW='/home/CASEHERE.raw'
%
%       where 'home' can be replaced with the directory structure leading to the
%       location of the .raw files to be analyzed. The corresponding outputs
%       are similarly structured: 
% 
%           ltable=7, filtab='/home/CASEHERE.table'
%           lprint=6, filpri='/home/CASEHERE.dump'
%           lcsv=11, filcsv='/home//CASEHERE.csv'
%           FILPS='/home/CASEHERE.ps'
% 
%       where, again, 'home' can be replaced with the actual directory structure
%       leading to the desired output files from the LCModel analysis.
% 
%       2) LCModel .raw files for which CONTROL files are to be written, inside
%       folder named 'Data_RAW'
% 
% Function run prompts user to select root folder containing both CONTROL_FIles and Data_RAW directories 
% 
% Outputs: LCModel-readable CONTROL files containing same parameters as
% template.control; one for each .raw file in Data_RAW folder 
% 
% Author: Kelley Swanberg (Columbia University, 2023) 
% 
% 
current_folder_root = uigetdir(); 
current_folder = sprintf('%s\\Data_RAW\', current_folder_root); 
export_folder = sprintf('%s\\CONTROL_Files\', current_folder_root); 

% Find .raw files for which to write CONTROL files
list_of_cases_struct = dir(current_folder);
list_of_cases_struct_clean = list_of_cases_struct(~ismember({list_of_cases_struct.name},{'.','..'}));
list_of_cases = {list_of_cases_struct_clean.name}'; 

% Now we have the list of cases in memory as a vertical vector 
num_cases = length(list_of_cases);  

for i = 1:num_cases
    
    case_name = extractBefore(list_of_cases{i},'.raw');
    
    fid  = fopen('template.control','r');
    f=fread(fid,'*char')';
    fclose(fid);
    f = strrep(f,'CASEHERE', case_name);

    filename = sprintf('%s\\%s_export.control', export_folder, case_name); 
    fid  = fopen(filename,'w');
    fprintf(fid,'%s',f);
    fclose(fid);
    
end

end