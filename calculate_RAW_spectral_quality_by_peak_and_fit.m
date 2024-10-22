% calculate_RAW_spectral_quality
function calculate_RAW_spectral_quality_by_peak_and_fit() 
% Calculates and reports spectral and fit quality values from directory of RAW outputs 
% 
% Inputs: LCModel RAW files to use in spectral quality asssessments 
% 
% plot_baseline - 0 or 1 to plot baseline model; 
% plot_fwhm - 0 or 1 to plot FWHM parameters; 
% ppm_end - Lowest ppm value in 'signal' range (e.g. 2.8 for tCr) for SNR/FWHM 
% ppm_start -  Highest ppm value in 'signal range' (e.g. 3.15 for tCr) for SNR/FWHM 
% ppm_base_end - Lowest ppm value to define noise baseline range (e.g. 0)
% ppm_base_start - Lowest ppm value to define noise baseline range (e.g. 0.5)
% ppm_subtract_end - Lowest ppm value to define signal baseline range (e.g. 2.75 for tCr with macromolecules)
% ppm_subtract_start - Lowest ppm value to define signal baseline range (e.g. 2.85 for tCr with macromolecules)
% sf - Larmor frequency in MHz (e.g. 399.5 for 9.4 T)
% fid_cut - Where to cut the data for SNR and FWHM assessment (should
% reflect how data were modeled in spectral analysis step) 
% fid_zf - How much to ZF data for SNR and FWHM assessment 
% fid_lb - How much to Lorentzian line-broaden in Hz the spectra before FWHM but
% not SNR assessment (added broadening is subtracted from final FWHM value; useful for
% very low-SNR data for which noise may otherwise confound meaningful peak-picking) 
% 
% Function run prompts user to select root folder containing directory of
% LCModel RAW files from which to calculate SNR and FWHM
% 
% Outputs: 'RAW_Spectral_quality_outputs.csv' wherein input and output data are reported in
% twenty-eight columns: 
% 
% Column 1: Row names based on coord file name
% Column 2: sf input (Larmor frequency in MHz (e.g. 399.5 for 9.4 T))
% Column 3: sw input (spectral width in Hz (e.g. 5000 Hz))
% Column 4: ppm_calib input (expected ppm value for water (e.g. 4.7 ppm)) 
% Column 5: ppm_signal input (expected ppm value for fitted signal (e.g.
% 3.03 ppm for tCr)
% Column 6: fid_cut input (number of non-zero-filled complex points considered in FID)
% Column 7: fid_zf input (final number of zero-filled complex points in FID)
% Column 8: fid_lb input (Lorentzian line-broadening in Hz applied to
% Column 9: ppm_end input (Lowest ppm value in 'signal' range (e.g. 2.8 for
% tCr) for SNR/FWHM) 
% Column 10: ppm_start input (Highest ppm value in 'signal range' (e.g. 3.15
% for tCr) for SNR/FWHM)
% Column 11: ppm_base_end - input (Lowest ppm value to define noise baseline range
% (e.g. 0))
% Column 12: ppm_base_start input (Lowest ppm value to define noise baseline range
% (e.g. 0.5))
% Column 13: ppm_subtract_end input (Lowest ppm value to define signal baseline range
% (e.g. 2.75))
% Column 14: ppm_subtract_start input (Lowest ppm value to define signal baseline range
% (e.g. 2.85))
% Column 15: fit_peak_initial_FWHM input (starting value in Hz for signal peak fit) 
% Column 16: lambda_2 input (peak fit regularization lambda on frequency shift term difference from initial condition) 
% Column 17: lambda_3 input (peak fit regularization lambda on Lorentzian broadening term difference from initial condition) 
% Column 18: Baseline-corrected signal maximum (for peak-based FWHM but not SNR assessment) 
% Column 19: Noise standard deviation 
% Column 20: Peak-based SNR 
% Column 21: Peak-based FWHM range ppm end as calculated from interpolated spectrum for maximum precision  
% Column 22: Peak-based FWHM range ppm start as calculated from interpolated spectrum for maximum precision 
% Column 23: Peak-based signal FWHM (ppm) as calculated from line-broadened and interpolated spectrum for maximum precision 
% Column 24: Peak-based signal FWHM (Hz) as calculated from line-broadened and interpolated spectrum for maximum precision
% Column 25: Fit-based signal amplitude (a.u.) as calculated from fitted Lorentzian signal peak 
% Column 26: Fit-based SNR as calculated from fitted Lorentzian signal peak
% Column 27: Fit-based signal FWHM (Hz) as calculated from line-broadened and interpolated spectrum for maximum precision 
% 
% Will optionally plot baseline model for signal correction in SNR and FWHM calculation as well as 
% FWHM calculation parameters and save figures as png (see plot inputs) 
% 
% Author: Kelley Swanberg (Lunds universitet, 2024)
% 
% swanberg@post.harvard.edu 
% 
% Written and tested in MATLAB 2023b

%%************************************************************************
%% INPUTS
%%************************************************************************

    %% FILE HANDLING INPUTS
    output_dir = 'RAW_spectral_quality_outputs'; 
    
    %% PLOTTING INPUTS 
    plot_fit = 1; 
    plot_baseline = 1; 
    plot_fwhm = 1; 

    %% DATA INPUTS 
    sf = 399.515512048197; % MHz
    sw = 5000; % Hz 
    ppm_calib = 4.7; % ppm 
    ppm_signal = 3.03; 

    %% DATA PREPROCESSING INPUTS 
    fid_cut = 1472; 
    fid_zf = 16384; 
    fid_lb = 5; % Hz 

    %% PEAK-PICKING INPUTS
    ppm_end = 2.85; % Signal peak
    ppm_start = 3.1; % Signal peak
    ppm_base_end = 0; % Noise range
    ppm_base_start = 0.5; % Noise range
    ppm_subtract_end = 2.8; % Signal baseline for data-based (not fit-based) SNR and FWHM measurements
    ppm_subtract_start = 2.85; % Signal baseline for data-based (not fit-based) SNR and FWHM measurements

    %% FITTING INPUTS
    fit_peak_initial_FWHM = 12; % Starting FWHM for fitted peak in Hz
    lambda_2 = 1; % Regularization lambda on frequency shift term difference from initial condition 
    lambda_3 = 1; % Regularization lambda on line-broadening term difference from initial condition 
    
    %% DERIVED INPUTS (do not change)
    dwell_time = 1/sw; % s
    bw_ppm = sw / sf; 
   
%% ************************************************************************
%% FILE HANDLING 
%% ************************************************************************

    % Prompt user to select root directory containing folder of RAW files
    current_folder = uigetdir(); 
    
    % Find RAW files for which to write script 
    list_of_cases_struct = dir(current_folder);
    list_of_cases_struct_clean = list_of_cases_struct(~ismember({list_of_cases_struct.name},{'.','..'}));
    list_of_cases = {list_of_cases_struct_clean.name}'; 
    %list_of_cases = contains(list_of_cases, '.raw');
    
    % Now we have the list of cases in memory as a vertical vector 
    num_raw = length(list_of_cases);  
    
    % Prepare table for outputs 
    combined_array = zeros(num_raw, 26); 
    combined_table = array2table(combined_array); 
    row_names = {num_raw};
    
    mkdir(output_dir); 
 
%% ************************************************************************
%% DATA- AND FIT-BASED SPECTRAL QUALITY ESTIMATION 
%% ************************************************************************

    % Prepare ppm axis
    ppm_axis_vector_st = ppm_calib + bw_ppm/2; 
    ppm_axis_vector_end = ppm_calib - bw_ppm/2; 
    ppm_axis_vector = linspace(ppm_axis_vector_st, ppm_axis_vector_end, fid_zf)'; 

    % Prepare time-domain vector for use in fitting
    time_domain_vector = linspace(0, (fid_zf-1)*dwell_time, fid_zf)'; 

    % For each COORD file in directory calculate FWHM and SNR 
    for ii = 1:num_raw
    
    
        %% ************************************************************************
        %% RAW FILE CONVERSION TO COMPLEX FID 
        %% ************************************************************************

        % Find and load data RAW file
        case_name = list_of_cases{ii}; 
        fid=fopen(case_name,'r');
        f=fread(fid,'*char')';
        fclose(fid);
        f_lines = splitlines(f); 
        
        % Find first line of data
        start_index = find(contains(f_lines,'$END'));
        end_index = length(f_lines);
        num_data_lines = end_index - start_index - 1; 
        
        data_coords_line_first = start_index+1; 
        data_coords_line_last = data_coords_line_first + num_data_lines-1; 
        data_coords_lines = f_lines(data_coords_line_first:data_coords_line_last);
            
        % Convert lines in RAW file to comprehensible complex FID
        data_coords_vector_1d = coord_lines_to_vector(data_coords_lines); 
        data_coords_real = data_coords_vector_1d(1:2:end); 
        data_coords_imag = data_coords_vector_1d(2:2:end); 
        data_coords_vector_unprocessed = complex(data_coords_real, data_coords_imag);

        %% ************************************************************************
        %% FID PREPROCESSING
        %% ************************************************************************

        % Preprocess FID according to inputs for more reliable peak-picking
        % and fitting 
        data_length_unprocessed = length(data_coords_vector_unprocessed); 
        padsize = fid_zf - fid_cut; 
        fid_vector_unbroadened = padarray(data_coords_vector_unprocessed(1:fid_cut), padsize, 'post'); 
        fid_vector = exp(-time_domain_vector*pi*fid_lb) .* fid_vector_unbroadened;  
        
        % Find and prepare unprocessed ppm axis 
        ppm_axis_vector_unprocessed = linspace(ppm_axis_vector_st, ppm_axis_vector_end, data_length_unprocessed); 
        
        % Prepare unbroadened and broadened spectra 
        data_spec_unbroadened = fftshift(fft(fid_vector_unbroadened)); 
        data_spec_real_unbroadened = real(data_spec_unbroadened); 
        data_spec = fftshift(fft(fid_vector)); 
        data_spec_real = real(data_spec); 
        
        %% ************************************************************************
        %% RAW PEAK-BASED SNR and FWHM ESTIMATION
        %% ************************************************************************

        % Calculate amplitude of signal in ppm range 
        % https://se.mathworks.com/matlabcentral/answers/
        % 152301-find-closest-value-in-array#answer_336210
        [minValue_start, ppm_start_index] = min(abs(ppm_axis_vector - ppm_start)); 
        [minValue_end, ppm_end_index] = min(abs(ppm_axis_vector - ppm_end)); 
        signal_ppm_selection = data_spec_real(ppm_start_index:ppm_end_index);
        signal_ppm_selection_unbroadened = data_spec_real_unbroadened(ppm_start_index:ppm_end_index); % Do not line-broaden for this
        
        % Calculate baseline of signal (use broadened spectrum)
        [base_minValue_start, ppm_base_start_index] = min(abs(ppm_axis_vector - ppm_base_start)); 
        [base_minValue_end, ppm_base_end_index] = min(abs(ppm_axis_vector - ppm_base_end)); 
        baseline_ppm_selection = data_spec_real(ppm_base_start_index:ppm_base_end_index); 
        ppm_selection_x = ppm_axis_vector(ppm_base_start_index:ppm_base_end_index); 
        p = polyfit(ppm_selection_x,baseline_ppm_selection,0);
        p_baseline_model = polyval(p,ppm_axis_vector); 
        
        % Calculate noise from baseline area (use unbroadened spectrum)
        noise_baseline_corr = data_spec_real_unbroadened(ppm_base_start_index:ppm_base_end_index) - p_baseline_model(ppm_base_start_index:ppm_base_end_index); 
        noise_std = std(noise_baseline_corr);
    
        % Calculate subtraction baseline for FWHM and signal amplitude 
        [subtract_minValue_start, ppm_subtract_start_index] = min(abs(ppm_axis_vector - ppm_subtract_start)); 
        [subtract_minValue_end, ppm_subtract_end_index] = min(abs(ppm_axis_vector - ppm_subtract_end)); 
        subtract_ppm_selection = data_spec_real(ppm_subtract_start_index:ppm_subtract_end_index); 
        subtract_ppm_selection_x = ppm_axis_vector(ppm_subtract_start_index:ppm_subtract_end_index); 
        subtract_p = polyfit(subtract_ppm_selection_x,subtract_ppm_selection,0);
        subtract_p_baseline_model = polyval(subtract_p,ppm_axis_vector); 
        subtract_amp_min = mean(subtract_p_baseline_model(ppm_subtract_start_index:ppm_subtract_end_index)); 
        subtract_amp_min_bc = mean(subtract_amp_min - p_baseline_model(ppm_start_index:ppm_end_index)); 
    
         % Calculate signal max and SNR from unbroadened spectrum for SNR calculation 
        signal_baseline_corr_unbroadened = signal_ppm_selection_unbroadened - p_baseline_model(ppm_start_index:ppm_end_index);
        [signal_min_baseline_corr_unbroadened, signal_min_unbroadened_ppm] = min(signal_baseline_corr_unbroadened); 
        [signal_max_baseline_corr_unbroadened, signal_max_unbroadened_ppm] = max(signal_baseline_corr_unbroadened);  
        signal_amp_baseline_corr_unbroadened = signal_max_baseline_corr_unbroadened - subtract_amp_min_bc; 
        signal_noise_ratio = signal_amp_baseline_corr_unbroadened / noise_std; 

        % Calculate signal max and FWHM from preprocessed spectrum using interpolated function for high precision 
        number_of_elements_new = 50000; 
        ppm_selected = ppm_axis_vector(ppm_start_index:ppm_end_index);
        ppm_selected_interpolated = linspace(ppm_start, ppm_end, number_of_elements_new); 
        signal_baseline_corr = signal_ppm_selection - p_baseline_model(ppm_start_index:ppm_end_index); 
        number_of_elements_old = length(signal_baseline_corr); 
        interpolated_domain = linspace(1, number_of_elements_old, number_of_elements_new); 
        signal_baseline_corr_interpolated = interp1(signal_baseline_corr, interpolated_domain, 'linear');
        %[signal_min_baseline_corr_nosub, signal_min_ppm] = min(signal_baseline_corr_interpolated); 
        [signal_max_baseline_corr_nosub, signal_max_ppm] = max(signal_baseline_corr_interpolated); 
        signal_max_baseline_corr = signal_max_baseline_corr_nosub - subtract_amp_min_bc; 
        signal_half_max_baseline_corr = signal_max_baseline_corr / 2 + subtract_amp_min_bc; 
        [signal_half_max_ppm_bigger, signal_half_max_ppm_bigger_index] = min(abs(signal_baseline_corr_interpolated(signal_max_ppm:end) - signal_half_max_baseline_corr)); 
        [signal_half_max_ppm_smaller, signal_half_max_ppm_smaller_index] = min(abs(signal_baseline_corr_interpolated(1:signal_max_ppm) - signal_half_max_baseline_corr));
        signal_fwhm_ppm = ppm_selected_interpolated(signal_half_max_ppm_smaller_index) - ppm_selected_interpolated(signal_half_max_ppm_bigger_index + signal_max_ppm - 1) - fid_lb/sf; 
        signal_fwhm_Hz = signal_fwhm_ppm * sf; 

        %% ************************************************************************
        %% LORENTZIAN FIT-BASED SNR and FWHM ESTIMATION
        %% ************************************************************************

        % Calculate signal max and FWHM by fitting Lorentzian signal to preprocessed spectrum within input ppm range   
        % Fitting parameters (Levenberg-Marquardt with regularization on soft constraints for similarity to
        % LCModel) 
        f = @(x, y)LorentzReg(x, y, noise_std, lambda_2, lambda_3, ppm_start_index, ppm_end_index, sf); 
        opts = optimoptions(@lsqcurvefit); % reuse the options
        opts.Algorithm = 'levenberg-marquardt'; 
        opts.FunctionTolerance = 1e-10; 
        opts.OptimalityTolerance = 1e-10; 
        opts.MaxIterations=1000; 

        % Initial conditions 
        x0_1 = max(data_spec_real)/sqrt(fid_zf); % Signal amplitude in time domain 
        x0_2 = (ppm_signal - ppm_calib); % Signal ppm shift from water in frequency domain 
        x0_3 = fit_peak_initial_FWHM + fid_lb; % Signal line width in Hz in frequency domain 
        x0 = [x0_1, x0_2, x0_3]; 
    
        % Set up x and y for fit, including real and imaginary components
        % plus regularization terms 
        xdata = time_domain_vector; 
        ydata_real = real(data_spec);
        ydata_imag = imag(data_spec); 
        ydata_reg = noise_std*[lambda_2*x0_2; lambda_3*x0_3]; 
        ydata = [ydata_real(ppm_start_index:ppm_end_index); ydata_imag(ppm_start_index:ppm_end_index); ydata_reg]; 
    
        % Sweet fitting action 
        [vestimated,resnorm] = lsqcurvefit(f,x0,xdata,ydata, [], [], opts); 
        
        % Reconstruct fitted peak (in time domain) from fit outputs 
        fitted_peak = vestimated(1) * exp(-2*pi*sf*vestimated(2)*j*xdata) .* exp(-xdata*pi*(vestimated(3)- fid_lb)); 
        fitted_peak_spec_real = real(fftshift(fft(fitted_peak))); 
        
        % Signal and FWHM outputs by fitting 
        int_fitted = sqrt(fid_zf)*vestimated(1);
        signal_fitted = max(fitted_peak_spec_real) - min(fitted_peak_spec_real);
        SNR_fitted = signal_fitted/noise_std; 
        fwhm_fitted_Hz = vestimated(3) - fid_lb;

        %% ************************************************************************
        %% RESULTS PLOTTING 
        %% ************************************************************************

        % Plot fit result as sanity check 
        if plot_fit==1
            figure(1); 
            plot(ppm_axis_vector, data_spec_real_unbroadened, 'Color', 'k', 'DisplayName','Spectrum (cut and zero-filled)');
            hold on; 
            plot(ppm_axis_vector, fitted_peak_spec_real, 'Color', 'r', 'DisplayName','Fitted peak'); 
            hold off; 
            set(gca, 'XDir','reverse'); 
            xlim('tight'); 
            xlabel('Chemical shift (ppm)');
            ylabel('Signal (a.u.)'); 
            title(['Fitted peak SNR: ',num2str(SNR_fitted, '%.1f'), '; FWHM: ' num2str(fwhm_fitted_Hz, '%.1f'), ' Hz']); 
            legend('Location','southoutside'); 
            filename = sprintf('%s//%s_fitted_peak.png', output_dir, case_name); 
            saveas(gcf, filename); 
        end

        % Plot baseline and export figure as sanity check  
        if plot_baseline==1
            figure(2); 
            plot(ppm_axis_vector, data_spec_real, 'Color', 'k', 'DisplayName','Spectrum (cut, zero-filled, line-broadened)');
            hold on; 
            plot(ppm_axis_vector, p_baseline_model, 'Color', 'r', 'DisplayName', 'Noise baseline'); 
            xline(ppm_base_start, 'Color', 'r', 'LineStyle', '--', 'DisplayName','Noise baseline limits'); 
            xline(ppm_base_end, 'Color', 'r', 'LineStyle', '--', 'DisplayName','Noise baseline limits'); 
            yline(subtract_amp_min, 'Color', 'b', 'DisplayName','Peak baseline'); 
            hold off; 
            set(gca, 'XDir','reverse'); 
            xlim('tight'); 
            xlabel('Chemical shift (ppm)'); 
            ylabel('Signal (a.u.)'); 
            title('Preprocessed spectrum and noise baseline'); 
            legend('Location','southoutside'); 
            filename = sprintf('%s//%s_baseline_model.png', output_dir, case_name); 
            saveas(gcf, filename); 
        end
        
        % Plot FWHM ppm range and baseline as sanity check 
        if plot_fwhm==1
            figure(3); 
            plot(ppm_selected_interpolated, signal_baseline_corr_interpolated, 'Color', 'k', 'DisplayName', 'Peak (cut, zero-filled, line-broadened, noise baseline subtracted)');
            hold on; 
            yline(signal_half_max_baseline_corr, 'Color', 'b', 'DisplayName', 'Peak half maximum'); 
            xline(ppm_selected_interpolated(signal_half_max_ppm_smaller_index), 'Color', 'b', 'DisplayName', 'Peak half maximum ppm value'); 
            xline(ppm_selected_interpolated(signal_half_max_ppm_bigger_index + signal_max_ppm - 1), 'Color', 'b', 'DisplayName', 'Peak half maximum ppm value'); 
            yline(subtract_amp_min_bc, 'Color', 'r', 'DisplayName', 'Peak baseline'); 
            hold off; 
            set(gca, 'XDir','reverse'); 
            xlim('tight'); 
            xlabel('Chemical shift (ppm)'); 
            ylabel('Signal (a.u.)'); 
            title(['Peak SNR: ',num2str(signal_noise_ratio, '%.1f'), '; FWHM: ' num2str(signal_fwhm_Hz, '%.1f'), ' Hz']); 
            legend('Location','southoutside'); 
            filename = sprintf('%s//%s_fwhm_base_corrected.png', output_dir, case_name); 
            saveas(gcf, filename); 
        end
        

        %% ************************************************************************
        %% OUTPUT TABLE PREPARATION 
        %% ************************************************************************

        % Output spectral quality estimates into table for writing as CSV 
        row_names{ii} = case_name; 
        combined_table{ii, 1} = sf; 
        combined_table{ii, 2} = sw; 
        combined_table{ii, 3} = ppm_calib;
        combined_table{ii, 4} = ppm_signal;
        combined_table{ii, 5} = fid_cut;  
        combined_table{ii, 6} = fid_zf;  
        combined_table{ii, 7} = fid_lb;
        combined_table{ii, 8} = ppm_end; 
        combined_table{ii, 9} = ppm_start; 
        combined_table{ii, 10} = ppm_base_end; 
        combined_table{ii, 11} = ppm_base_start; 
        combined_table{ii, 12} = ppm_subtract_end; 
        combined_table{ii, 13} = ppm_subtract_start; 
        combined_table{ii, 14} = fit_peak_initial_FWHM; 
        combined_table{ii, 15} = lambda_2; 
        combined_table{ii, 16} = lambda_3; 
        combined_table{ii, 17} = signal_amp_baseline_corr_unbroadened;
        combined_table{ii, 18} = noise_std;  
        combined_table{ii, 19} = signal_noise_ratio; 
        combined_table{ii, 20} = ppm_selected_interpolated(signal_half_max_ppm_bigger_index + signal_max_ppm - 1); 
        combined_table{ii, 21} = ppm_selected_interpolated(signal_half_max_ppm_smaller_index);
        combined_table{ii, 22} = signal_fwhm_ppm;
        combined_table{ii, 23} = signal_fwhm_Hz; 
        combined_table{ii, 24} = signal_fitted; 
        combined_table{ii, 25} = SNR_fitted; 
        combined_table{ii, 26} = fwhm_fitted_Hz;
    end
      
    %% Add row names to output table 
    combined_table.Properties.RowNames = row_names; 
    combined_table.Properties.VariableNames = {'Larmor frequency (MHz)', 'Spectral width (Hz)', 'Water chemical shift (ppm)', 'Signal chemical shift (ppm)', ...
    'Complex FID pts analyzed', 'Zero-filled length', 'Lorentzian line broadening applied (Hz)', ...
    'Signal end (ppm)', 'Signal start (ppm)', 'Noise baseline end (ppm)', 'Noise baseline start (ppm)', 'Signal calculation baseline end (ppm)', 'Signal calculation baseline start (ppm)', ...
    'Peak fit initial FWHM (Hz)', 'Peak fit frequency term lambda', 'Peak fit Lorentzian broadening term lambda', ... 
    'Signal peak amplitude', 'Noise std', 'SNR', 'Signal peak FWHM end (ppm)', 'Signal peak FWHM start (ppm)', ...
    'Signal peak FWHM (ppm)', 'Signal peak FWHM (Hz)', 'Fitted signal amplitude', 'Fitted SNR', 'Fitted FWHM (Hz)'}; 
    
    CSV_name = sprintf('%s//RAW_Spectral_and_fit_quality_outputs.csv', output_dir); 
    
    %% Output CSV of quality assessment results 
    writetable(combined_table, CSV_name, 'WriteRowNames', true); 
    
    close all; 

end

function y = LorentzReg(v, xdata, noise_std, lambda_2, lambda_3, ppmstart, ppmend, sf)

   output = fftshift(fft(v(1) * exp(-2*pi*sf*v(2)*j*xdata) .* exp(-xdata*pi*v(3))));
   y_real = real(output);
   y_imag = imag(output);
   y_reg = noise_std*[lambda_2*v(2); lambda_3*v(3)]; 
   y = [y_real(ppmstart:ppmend); y_imag(ppmstart:ppmend); y_reg]; 
   
end

% Convert LCModel RAW lines to readable complex vector 
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


