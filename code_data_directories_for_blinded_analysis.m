function code_data_directories_for_blinded_analysis()
% 
% Enables self-blinded data analysis by coding the root names of a 
% directory of data subdirectories, copying all
% subdirectory contents into coded folders, and outputting a CSV key.
%
% Inputs: Function run prompts user input of folder containing only
% subdirectories of data to be coded 
%   Example /Input_Directory/
%   /Input_Directory/Non-coded_data_ID_01/Data_Contents
%   /Input_Directory/Non-coded_data_ID_02/Data_Contents
%   /Input_Directory/Non-coded_data_ID_03/Data_Contents
%   ...
%   /Input_Directory/Non-coded_data_ID_n/Data_Contents
% 
% Outputs: 
% 
% 1. Output directory "Data_Coded" containing each data subdirectory with
% coded name plus original contents (note: Contents are not themselves coded) 
% Example: /Data_Coded/
%   /Data_Coded/71936/Data_Contents
%   /Data_Coded/96381/Data_Contents
%   /Data_Coded/93810/Data_Contents
%   ...
%   /Data_Coded/02947/Data_Contents
% 
% 2. Data coding key matching 'Non-coded_data_ID_n' with each code: 'Coded_Data_Key.csv'
% 
% Author: Kelley Swanberg (Lunds universitet, 2024) 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Locate correct folder; must contain only subdirectories to be coded
    data_to_code_path = uigetdir(); 
    
    data_to_code = dir(data_to_code_path); 
    data_to_code_clean = data_to_code(~ismember({data_to_code.name},{'.','..'}));
    data_to_code = {data_to_code_clean.name}'; 
    
    data_codes = data_to_code; 
    
    data_coded_path = 'Data_Coded'; 
    mkdir(data_coded_path); 
    random_num = zeros(1,5);
  
%% Create randomized data codes and folders; copy data contents to new coded directories
    for ii=1:length(data_codes)
        for jj=1:5
            random_num(1,jj) = randi(9); 
        end
        random_num_string = strrep(strcat(num2str(random_num(1,:))), ' ', ''); 
        data_codes{ii} = random_num_string; 
    
        oldpath = sprintf('%s//%s', data_to_code_path, data_to_code{ii}); 
        newpath = sprintf('%s//%s', data_coded_path, data_codes{ii}); 
        copyfile(oldpath, newpath); 
    end
    
 %% Output code key 
    data_coded = cell2table([data_to_code, data_codes]); 
    tablename = 'Coded_Data_Key.csv'; 
    writetable(data_coded, tablename); 

end

