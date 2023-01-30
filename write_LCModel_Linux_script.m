function write_LCModel_Linux_script()
% 
% Composes a .sh shell script for use on Linux systems to batch-schedule
% LCModel runs over a directory of CONTROL files 
%
% Inputs: LCModel CONTROL files to include in Linux shell script 
% 
% Function run prompts user to select root folder containing directory of
% CONTROL_Files for which to write shell script 
% 
% Outputs: .sh script running LCModel once for each CONTROL file in the input directory 
% Note: To run this script, a user must first set execute permissions via
% 'chmod -x [filename including path].sh'
% 
% Author: Kelley Swanberg (Columbia University, 2023) 
%
% Prompt user to select root directory containing folder of CONTROL files
% for which to write .sh script 
current_folder_root = uigetdir(); 
current_folder = sprintf('%s\\CONTROL_Files\', current_folder_root); 

% Find CONTROL files for which to write script 
list_of_cases_struct = dir(current_folder);
list_of_cases_struct_clean = list_of_cases_struct(~ismember({list_of_cases_struct.name},{'.','..'}));
list_of_cases = {list_of_cases_struct_clean.name}'; 

% Now we have the list of cases in memory as a vertical vector 
num_cases = length(list_of_cases);  

% Set up .sh script 
sh_script_filename = 'LCModel_Control_File_Run.sh'; 
fileID = fopen(sh_script_filename,'w');
fprintf(fileID,'#/usr/bin/bash\n');

% Write LCModel calls for each CONTROL file in input directory 
for i = 1:num_cases
    fprintf(fileID,'$HOME/.lcmodel/bin/lcmodel <%s;\n',list_of_cases{i});    
end

fclose(fileID);

end