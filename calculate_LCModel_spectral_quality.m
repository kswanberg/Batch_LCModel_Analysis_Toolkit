function calculate_LCModel_spectral_quality()
% Calculates and reports spectral and fit quality values from directory of COORD outputs 
% 
% Inputs: LCModel COORD files to use in spectral and fit quality
% asssessments 
% 
% plot_residual - 0 or 1 to plot residual  
% plot_baseline - 0 or 1 to plot baseline model; 
% plot_fwhm - 0 or 1 to plot FWHM parameters; 
% ppm_end - Lowest ppm value in 'signal' range (e.g. 2.8 for tCr) for SNR/FWHM 
% ppm_start -  Highest ppm value in 'signal range' (e.g. 3.15 for tCr) for SNR/FWHM 
% ppm_base_end - Lowest ppm value to define baseline range (e.g. 0)
% ppm_base_start - Lowest ppm value to define baseline range (e.g. 0.5)
% sf - Larmor frequency in MHz (e.g. 399.5 for 9.4 T)
% 
% Function run prompts user to select root folder containing directory of
% COORD files from which to calculate residuals and normality thereof
% 
% Outputs: 'Spectral_and_fit_quality_outputs.csv' wherein data are reported in
% five columns: 
% 
% Column 0: Row names based on coord file name
% Column 1: Fit residual K-S test p-value < default alpha? 1=yes; 0=no 
% Column 2: Fit residual K-S test p-value 
% Column 3: Fit residual K-S test statistic
% Column 4: Fit residual K-S test critical value
% Column 5: Fit residual kurtosis value 
% Column 6: Fit residual standard deviation
% Column 7: Noise (frequency domain) standard deviation
% Column 8: Fit quality number (Near et al., 2022 doi: 10.1002/nbm.4257) 
% Column 9: Baseline-corrected signal maximum 
% Column 10: SNR
% Column 11: Signal FWHM (ppm) 
% Column 12: Signal FWHM (Hz) 
% 
% Will optionally plot residuals, baseline model for signal correction in SNR calculation, 
% and FWHM calculation parameters and save figure as png (see plot inputs) 
% 
% Author: Kelley Swanberg (Lunds universitet, 2024)
% 
% swanberg@post.harvard.edu 
% 
% Written and tested in MATLAB 2023b

    %% Input values
    plot_residual = 1; 
    plot_baseline = 1; 
    plot_fwhm = 1; 
    ppm_end = 2.8; 
    ppm_start = 3.15;
    ppm_base_end = 0; 
    ppm_base_start = 0.5; 
    sf = 399.515512048197; 
    
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
    combined_array = zeros(num_coords, 12); 
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

    % Find and load baseline coords 
    baseline_coords_line_first = find(contains(f_lines,'background values follow'))+1;
    baseline_coords_line_last = baseline_coords_line_first + num_data_lines-1; 
    baseline_coords_lines = f_lines(baseline_coords_line_first:baseline_coords_line_last);
        
    baseline_coords_vector = coord_lines_to_vector(baseline_coords_lines); 
    
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
    
    % Calculate amplitude of signal in ppm range 
    % https://se.mathworks.com/matlabcentral/answers/
    % 152301-find-closest-value-in-array#answer_336210
    [minValue_start, ppm_start_index] = min(abs(ppm_axis_vector - ppm_start)); 
    [minValue_end, ppm_end_index] = min(abs(ppm_axis_vector - ppm_end)); 
    signal_ppm_selection = data_coords_vector(ppm_start_index:ppm_end_index); 
    
    % Calculate baseline of signal 
    [base_minValue_start, ppm_base_start_index] = min(abs(ppm_axis_vector - ppm_base_start)); 
    [base_minValue_end, ppm_base_end_index] = min(abs(ppm_axis_vector - ppm_base_end)); 
    baseline_ppm_selection = data_coords_vector(ppm_base_start_index:ppm_base_end_index); 
    ppm_selection_x = ppm_axis_vector(ppm_base_start_index:ppm_base_end_index); 
    p = polyfit(ppm_selection_x,baseline_ppm_selection,0);
    p_baseline_model = polyval(p,ppm_axis_vector); 

    % Calculate noise from baseline area 
    noise_baseline_corr = data_coords_vector(ppm_base_start_index:ppm_base_end_index) - p_baseline_model(ppm_base_start_index:ppm_base_end_index); 
    noise_std = std(noise_baseline_corr); 
 
    % Calculate standard deviation of residual 
    resid_std = std(residual_coords_vector); 
    
    % Calculate normality (one-sample Kolmogorov-Smirnov test and kurtosis
    % assessment) of residual 
    [res_h,res_p,res_ksstat,res_cv] = kstest(residual_coords_vector); 
    res_k = kurtosis(residual_coords_vector); 
    [k,res_autocorr] = ac2rc(xcorr(residual_coords_vector)); 

    % Calculate fit quality number (residual / noise std) 
    fqn = resid_std / noise_std; 

    % Calculate signal max and FWHM using interpolated function for high
    % precision 
    number_of_elements_new = 5000; 
    ppm_selected = ppm_axis_vector(ppm_start_index:ppm_end_index);
    ppm_selected_interpolated = linspace(ppm_start, ppm_end, 5000); 
    signal_baseline_corr = signal_ppm_selection - p_baseline_model(ppm_start_index:ppm_end_index); 
    number_of_elements_old = length(signal_baseline_corr); 
    interpolated_domain = linspace(1, number_of_elements_old, number_of_elements_new); 
    signal_baseline_corr_interpolated = interp1(signal_baseline_corr, interpolated_domain, 'linear');
    [signal_min_baseline_corr_nosub, signal_min_ppm] = min(signal_baseline_corr_interpolated); 
    [signal_max_baseline_corr_nosub, signal_max_ppm] = max(signal_baseline_corr_interpolated); 
    signal_max_baseline_corr = signal_max_baseline_corr_nosub - signal_min_baseline_corr_nosub; 
    signal_half_max_baseline_corr = signal_max_baseline_corr / 2 + signal_min_baseline_corr_nosub; 
    [signal_half_max_ppm_bigger, signal_half_max_ppm_bigger_index] = min(abs(signal_baseline_corr_interpolated(signal_max_ppm:end) - signal_half_max_baseline_corr)); 
    [signal_half_max_ppm_smaller, signal_half_max_ppm_smaller_index] = min(abs(signal_baseline_corr_interpolated(1:signal_max_ppm) - signal_half_max_baseline_corr));
    signal_fwhm_ppm = ppm_selected_interpolated(signal_half_max_ppm_smaller_index) - ppm_selected_interpolated(signal_half_max_ppm_bigger_index + signal_max_ppm - 1); 
    signal_fwhm_Hz = signal_fwhm_ppm * sf; 

    % Calculate SNR 
    signal_noise_ratio = signal_max_baseline_corr / noise_std; 

    row_names{ii} = case_name; 
    combined_table{ii, 1} = res_h; 
    combined_table{ii, 2} = res_p; 
    combined_table{ii, 3} = res_ksstat; 
    combined_table{ii, 4} = res_cv; 
    combined_table{ii, 5} = res_k; 
    combined_table{ii, 6} = resid_std;  
    combined_table{ii, 7} = noise_std;  
    combined_table{ii, 8} = fqn;  
    combined_table{ii, 9} = signal_max_baseline_corr; 
    combined_table{ii, 10} = signal_noise_ratio; 
    combined_table{ii, 11} = signal_fwhm_ppm;
    combined_table{ii, 12} = signal_fwhm_Hz; 

    % Plot residual and export figure as sanity check  
    if plot_residual==1
        figure(1); 
        plot(ppm_axis_vector, residual_coords_vector); 
        set(gca, 'XDir','reverse'); 
        xlim('tight'); 
        filename = sprintf('%s_resid.png', case_name); 
        saveas(gcf, filename); 
    end

    % Plot baseline and export figure as sanity check  
    if plot_baseline==1
        figure(2); 
        plot(ppm_axis_vector, data_coords_vector, 'Color', 'k');
        hold on; 
        plot(ppm_axis_vector, p_baseline_model, 'Color', 'r')
        xline(ppm_base_start, 'Color', 'r'); 
        xline(ppm_base_end, 'Color', 'r'); 
        hold off; 
        set(gca, 'XDir','reverse'); 
        xlim('tight'); 
        filename = sprintf('%s_baseline_model.png', case_name); 
        saveas(gcf, filename); 
    end

    if plot_fwhm==1
        figure(3); 
        plot(ppm_selected_interpolated, signal_baseline_corr_interpolated, 'Color', 'k');
        hold on; 
        yline(signal_half_max_baseline_corr, 'Color', 'c'); 
        xline(ppm_selected_interpolated(signal_half_max_ppm_smaller_index), 'Color', 'c'); 
        xline(ppm_selected_interpolated(signal_half_max_ppm_bigger_index + signal_max_ppm - 1), 'Color', 'c'); 
        hold off; 
        set(gca, 'XDir','reverse'); 
        xlim('tight'); 
        filename = sprintf('%s_fwhm_base_corrected.png', case_name); 
        saveas(gcf, filename); 
    end
   
    end 
    
    %% Add row names to output table 
    combined_table.Properties.RowNames = row_names; 
    combined_table.Properties.VariableNames = {'Res K-S stat < alpha?','Res K-S stat p-value', 'Res K-S stat', ...
    'Res K-S stat critical value', 'Res kurtosis', 'Res std', 'Noise std', 'FQN', 'Signal max', ...
    'SNR', 'Signal FWHM (ppm)', 'Signal FWHM (Hz)'}; 

    CSV_name = 'Spectral_and_fit_quality_outputs.csv'; 
    
    %% Output CSV of quality assessment results 
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