function calculate_LCModel_resid_normality()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Calculates and reports residual normality from directory of COORD outputs 
% 
%
% Inputs: LCModel COORD files to use in residual calculation and normality assessment 
% 
% Function run prompts user to select root folder containing directory of
% COORD files from which to calculate residuals and normality thereof
% 
% Outputs: 'Residual_normality_outputs.csv' wherein data are reported in
% five columns: 
% 
% Column 0: Row names based on coord file name
% Column 1: Fit residual K-S test p-value < default alpha? 1=yes; 0=no 
% Column 2: Fit residual K-S test p-value 
% Column 3: Fit residual K-S test statistic
% Column 4: Fit residual K-S test critical value
% Column 5: Fit residual kurtosis value 
% 
% Will optionally plot residual and save figure as png (see plot inputs) 
% 
% Author: Kelley Swanberg (Lunds universitet, 2024) 
% swanberg@post.harvard.edu 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Input values
    plot_residual = 0; 
    
    %% Prompt user to select root directory containing folder of COORD files
    % for which to write .sh script 
    current_folder = uigetdir(); 
    
    %% Find COORD files for which to write script 
    list_of_cases_struct = dir(current_folder);
    list_of_cases_struct_clean = list_of_cases_struct(~ismember({list_of_cases_struct.name},{'.','..'}));
    list_of_cases = {list_of_cases_struct_clean.name}'; 
    % https://se.mathworks.com/matlabcentral/answers/431023-list-all-and-only-files-with-no-extension#answer_348046
    list_of_cases_notcoord = contains(list_of_cases, '.');
    list_of_cases(list_of_cases_notcoord) = []; 
    
    % Now we have the list of cases in memory as a vertical vector 
    num_coords = length(list_of_cases);  
    
    %% Prepare table for outputs 
    combined_array = zeros(num_coords, 5); 
    combined_table = array2table(combined_array); 
    row_names = {num_coords};
    
    %% For each COORD file in directory calculate residual and test normality 
    for ii = 1:num_coords
    
    % Find and load data coords
    case_name = list_of_cases{ii}; 
    
    fid=fopen(case_name,'r');
    f=fread(fid,'*char')';
    fclose(fid);
    
    f_lines = splitlines(f); 
    
    % Find first line of data
    start_index = find(contains(f_lines,'points on ppm-axis'));
    end_index = find(contains(f_lines,'phased data points follow'));
    num_data_lines = end_index - start_index - 1; 
    
    % Find and prepare fit ppm axis 
    ppm_axis_line_first = start_index+1; 
    ppm_axis_line_last = end_index-1; 
    ppm_axis_lines = f_lines(ppm_axis_line_first:ppm_axis_line_last); 
    
    ppm_axis_vector = coord_lines_to_vector(ppm_axis_lines); 
    
    % Find and prepare data coords 
    data_coords_line_first = end_index+1; 
    data_coords_line_last = data_coords_line_first + num_data_lines-1; 
    data_coords_lines = f_lines(data_coords_line_first:data_coords_line_last);
    
    data_coords_vector = coord_lines_to_vector(data_coords_lines); 
    
    % Find and load fit coords 
    fit_coords_line_first = data_coords_line_last+2; 
    fit_coords_line_last = fit_coords_line_first + num_data_lines-1; 
    fit_coords_lines = f_lines(fit_coords_line_first:fit_coords_line_last);
    
    fit_coords_vector = coord_lines_to_vector(fit_coords_lines); 
    
    % Calculate residual coords 
    residual_coords_vector = data_coords_vector - fit_coords_vector; 
    
    % Plot residual and export figure as sanity check  
    if plot_residual==1
        figure(1); 
        plot(ppm_axis_vector, residual_coords_vector); 
        set(gca, 'XDir','reverse'); 
        filename = sprintf('%s_resid.png', case_name); 
        saveas(gcf, filename); 
    end
    
    % Calculate normality of residual using one-sample Kolmogorov-Smirnov test
    [h,p,ksstat,cv] = kstest(residual_coords_vector); 
    k = kurtosis(residual_coords_vector); 
    
    row_names{ii} = case_name; 
    combined_table{ii, 1} = h; 
    combined_table{ii, 2} = p; 
    combined_table{ii, 3} = ksstat; 
    combined_table{ii, 4} = cv; 
    combined_table{ii, 5} = k; 
    
    end 
    
    %% Add row names to normality table 
    combined_table.Properties.RowNames = row_names; 
    combined_table.Properties.VariableNames = {'Res K-S stat < alpha?','Res K-S stat p-value', 'Res K-S stat', 'Res K-S stat critical value', 'Res kurtosis'}; 

    CSV_name = 'Residual_normality_outputs.csv'; 
    
    %% Output CSV of residual normality results 
    writetable(combined_table, CSV_name, 'WriteRowNames', true); 

end 

function [coord_vector] = coord_lines_to_vector(coord_lines)

    coord_lines_split = split(coord_lines(1:end-1)); 
    coord_lines_split_1d = reshape(coord_lines_split.',1,[]);

    % Remove empty cells: https://se.mathworks.com/matlabcentral/answers/
    % 209-how-do-i-remove-the-empty-cells-from-a-vector-of-cells#answer_253
    coord_lines_split_1d_clean = coord_lines_split_1d(~cellfun('isempty',coord_lines_split_1d))'; 

    coord_lines_split_last = split(coord_lines(end)); 
    coord_lines_split_last_clean = coord_lines_split_last(~cellfun('isempty',coord_lines_split_last)); 

    coord_vector =str2double([coord_lines_split_1d_clean; coord_lines_split_last_clean]); 

end