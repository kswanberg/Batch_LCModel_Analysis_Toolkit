function add_LB(input_folder, LB_factor, spectral_width_for_raw)
% Add a constant Lorentzian line broadening factor to a directory of FIDs
% 
% Loads a directory of time-domain FIDs in LCModel .raw and/or INSPECTOR
% .mat format, adds the line broadening factor specified in Hz, and outputs
% the result in the same file format as the input. 
%
% Dependencies: None; support functions read_LCModel_raw() and create_fmt()
% included here
% 
% Inputs: 
% 
%  1. input_folder: Full path including directory of .mat and/or .raw FIDs
%  to preprocess with time-domain exponential decay / frequency-domain Lorentzian line broadening 
% 
%  2. LB_factor: Amount to broaden input FIDs in Hz
%
%  3. spectral_width_for_raw: Non-optional input because necessary for .raw
%  files. Spectral bandwidth in Hz (1/dwell time) necessary to compute
%  appropriate time domain for exponential decay / Lorentzian line
%  broadening application.
% 
% Outputs: 
% 
% Line-broadened fids in same format as inputs (.raw and/or
% INSPECTOR .mat) in a directory named after the input directory, as
% follows: '[input_folder]_with_LB_[LB_factor]_Hz'. 
% 
% Example usage: add_LB('F:\BasisSets\SummedSpins_for_MARSSinput\', 1, 5000)
% 
% Author: Kelley Swanberg (Lund University, 2024) 
% 
% Written and tested in MATLAB 2023b

    %% Load input directory and prepare output directory 

    % Load all files in folder 
    input_folder_split = strsplit(input_folder, '\'); 
    input_folder_target = input_folder_split{end-1}; 
    input_folder_root_split = strsplit(input_folder, input_folder_target); 
    input_folder_root = input_folder_root_split{1}; 
    
    % Create output folder
    output_folder = sprintf('%s\\%s_with_LB_%d_Hz\\', input_folder_root, input_folder_target, LB_factor); 
    mkdir(output_folder);
    
    % Determine if .mat or .raw 
    input_files_pathinfo = [dir(fullfile(input_folder, '*.mat')), dir(fullfile(input_folder, '*.raw'))];
    input_files_num = length(input_files_pathinfo); 
    
    %% Loop through input FIDs, add exponential decay corresponding with specified Lorentzian broadening in Hz, and output result as same format 
    for ii=1:input_files_num
        
        %% Load FID to broaden 
        % Specify input file
        input_file = input_files_pathinfo(ii).name
        input_path = input_files_pathinfo(ii).folder;
        input_file_path = sprintf('%s\\%s', input_path, input_file);
        [~,file_ID,ext] = fileparts(input_file_path);
        
        % If mat just load it 
        if ext == '.mat'
            fid_complex_struct = load(input_file_path); 
            fid_complex = fid_complex_struct.exptDat.fid; 
            spectral_width = fid_complex_struct.exptDat.sw_h; 
        % If raw use my loading function 
        else
            fid_complex = read_LCModel_raw(input_file_path); 
            spectral_width = spectral_width_for_raw; 
        end
    
        %% Broaden FID in time domain 
        % Figure out number of points 
        num_complex_points = length(fid_complex); 
    
        % Add Lorentzian line-broadening
        acq_time = (1/spectral_width)*num_complex_points; 
        time_domain = linspace(0, acq_time, num_complex_points); 
        fid_complex_broadened = exp(-pi * LB_factor * time_domain)' .* fid_complex; 
        
        %% Export broadened result 
        % Set up output file
        output_file = sprintf('%s\\%s', output_folder, input_file); 
    
        % If input was mat write mat file 
        if ext == '.mat'
            exptDat = fid_complex_struct.exptDat; 
            exptDat.fid = fid_complex_broadened; 
            save(output_file, 'exptDat'); 
    
        % If input was raw write raw file
        else
            created_date =  datetime('today'); 
            raw_volume = '1.00000E+00'; 
            raw_tramp = '1.00000E+00'; 
            % Write header information to raw file 
            fileID = fopen(output_file,'w');
            fprintf(fileID,' Created: %s\n', created_date);
            fprintf(fileID,' $NMID\n');
            fprintf(fileID,' ID=''%s''\n', file_ID);
            fprintf(fileID,' FMTDAT=''(8E13.5)''\n');
            fprintf(fileID,' VOLUME=   %s\n', raw_volume);
            fprintf(fileID,' TRAMP=   %s\n', raw_tramp);
            fprintf(fileID,' $END\n');
            fprintf(fileID,' ');
            
            % Prepare data vector for export as RAW
            raw_fid = fid_complex_broadened; 
            length_fid_complex = num_complex_points; 
            length_fid_values = 2 * length_fid_complex; 
            num_rows_raw = length_fid_values/8; 
    
            raw_fid_real_imag_sep = [real(raw_fid) +imag(raw_fid)]; % Negative imaginary component because otherwise they will end up upside-down in LCModel 
            raw_fid_real_imag_sep_transpose = raw_fid_real_imag_sep'; 
            raw_fid_real_imag_sep_vector = raw_fid_real_imag_sep_transpose(:); 
            raw_fid_matrix = reshape(raw_fid_real_imag_sep_vector, 8, num_rows_raw); 
    
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

    %% Support function: Read LCModel RAW file 
    % Also included separately in https://github.com/kswanberg/Batch_LCModel_Analysis_Toolkit/
    function fid_complex = read_LCModel_raw(filepath)
    
        file = filepath;
        
        rawfile = fileread(file);
        rawfile_split = strsplit(rawfile, '$END'); 
        
        file_fid = str2num(rawfile_split{1,2}); 
        file_fid_size = size(file_fid)
        file_fid_length = file_fid_size(1) * file_fid_size(2) / 2; % Assumes complex vector 
        file_fid_reshaped = reshape(file_fid', 2, [])'; 
        
        fid_complex = complex(file_fid_reshaped(:,1), file_fid_reshaped(:,2)); 
    
    end 

    %% Support function: Write FID values to specified precision for LCModel RAW format 
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

end