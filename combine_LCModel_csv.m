function combine_LCModel_csv()
% 
% Combines a folder of LCModel CSV outputs with the same output structure
% (i.e. consistent number and order of output metabolites) into a master
% document for easier downstream data handling 
%
% Inputs: Function run prompts user input of folder containing only CSVs to
% be collated 
% 
% Outputs: Single CSV file 'Collated_LCModel_Results.csv'
% 
% Author: Kelley Swanberg (Columbia University, 2023) 
% 
% Locate correct folder; must contain only CSV files to combine
current_folder = uigetdir(); 

% Find CSV files to combine 
list_of_cases_struct = dir(current_folder);
list_of_cases_struct_clean = list_of_cases_struct(~ismember({list_of_cases_struct.name},{'.','..'}));
list_of_cases = {list_of_cases_struct_clean.name}'; 

% Determine appropriate size of table to build
num_cases = length(list_of_cases);

combined_table = []; 
row_names = {};

for i = 1:num_cases
% Read and combine CSVs from input folder into table for later writing
csv_filename = list_of_cases{i};

case_ID = strrep(csv_filename, '.csv', ''); 
csv_mat = readtable(csv_filename);
row_names{i} = case_ID; 
combined_table = [combined_table; csv_mat]; 
end

% Add row names 
combined_table.Properties.RowNames = row_names; 

% Write final table to combined CSV
CSV_name = 'Collated_LCModel_Results.csv'
writetable(combined_table, CSV_name, 'WriteRowNames', true); 

end
