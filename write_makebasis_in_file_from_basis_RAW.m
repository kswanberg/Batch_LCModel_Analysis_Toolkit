function write_makebasis_in_file_from_basis_RAW()
% 
% Converts a folder of LCModel-format spectral .raw files to a single
% makebasis.in file for use in LCModel fits 
%
% Inputs: Function run prompts user input of folder containing only 
% LCModel-format .raw files to be converted to a LCModel-readable makebasis.in
% files according to additional inputs provided by user in 'General Inputs'
% and 'Raw File Inputs' (assumed constant for all bases in set) (see
% comments for both input sets) 
% 
% Outputs: LCModel-readable makebasis.in file written according to RAW
% files in chosen directory and additional function inputs as described
% above
% 
% Author: Kelley Swanberg (Lunds universitet, 2024) 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% General Inputs
    filename_ID = 'Varian_14p1T_13C_ISIS_10ms'; % Will form root of basis file names 
    sequence_name = 'ISIS'; % Pulse sequence 
    echo_time = '10.'; % ms
    fwhm_basis = '.049'; % ppm 
    sf = '150.732'; % MHz
    dwell_time = '0.0000248'; % s
    num_complex_pts = '8192'; % Number of complex points in basis FIDs
    basis_output_dir = '/home/nmr/Documents/Basis/' % In LCModel console 
    autosc_tf = 'false'; % Auto-scale basis functions in LCModel: true; do not auto-scale: false 
    autoph_tf = 'false'; % Auto-phase basis functions in LCModel: true; do not auto-scale: false 
    basis_id = 'Varian 14.1 T 13C ISIS (TE 10 ms) basis prepared in MATLAB'; % Will appear in LCModel makebasis.in output files
    basis_input_dir = '/home/nmr/Documents/Basis_Raw/Varian_14p1T_13C_ISIS_10ms/'; % Where the RAW bases live in the LCModel console 
    
    %% Raw File Inputs 
    zero_order_ph = '0.'; % Inherent zero-order phase of raw basis; assumes constant for all 
    first_order_ph = '0.'; % Inherent first-order phase of raw basis; assumes constant for all 
    scaling = '1.'; % Inherent scaling of raw basis; assumes constant for all 
    noshift_tf = 'TRUE'; % Should LCModel frequency-shift this basis? TRUE or FALSE; assumes constant for all 
    ppm_app_higher = '0.'; % Upper ppm limit for basis reference peak; assumes constant for all 
    ppm_app_lower = '-.4'; % Lower ppm limit for basis reference peak; assumes constant for all 
    
    % Find .raw files to convert into makebasis.in
    current_folder = uigetdir(); 
    list_of_cases_struct = dir(current_folder);
    list_of_cases_struct_clean = list_of_cases_struct(~ismember({list_of_cases_struct.name},{'.','..'}));
    list_of_cases = {list_of_cases_struct_clean.name}'; 
    
    % Determine appropriate size of table to build
    num_raw = length(list_of_cases);
    
    % Write header information to makebasis
    makebasisin_filename = sprintf('%s_makebasis.in', filename_ID); 
    fileID = fopen(makebasisin_filename,'w');
    fprintf(fileID,'$seqpar\n');
    fprintf(fileID,'seq=''%s''\n', sequence_name);
    fprintf(fileID,'echot=%s\n', echo_time);
    fprintf(fileID,'fwhmba=%s\n', fwhm_basis);
    fprintf(fileID,'$end\n');
    
    fprintf(fileID,'\n');
    
    fprintf(fileID,'$nmall\n');
    fprintf(fileID,'hzpppm=%s\n', sf);
    fprintf(fileID,'deltat=%s\n', dwell_time);
    fprintf(fileID,'nunfil=%s\n', num_complex_pts);
    fprintf(fileID, 'filbas=''%s%s.basis''\n', basis_output_dir, filename_ID); 
    fprintf(fileID, 'filps=''%s%s.ps''\n', basis_output_dir, filename_ID); 
    fprintf(fileID, 'autosc=.%s.\n', autosc_tf); 
    fprintf(fileID, 'autoph=.%s.\n', autoph_tf); 
    fprintf(fileID, 'idbasi=''%s''\n', basis_id); 
    fprintf(fileID,'$end\n');
    
    fprintf(fileID,'\n');
    
    % Now add metabolite-specific information to makebasis.in (note: 'Raw
    % File Inputs' assumed constant across metabolites here) 
    for ii = 1:num_raw
        raw_filename = list_of_cases{ii}; 
        raw_filename_fullfile = sprintf('%s/%s', current_folder, raw_filename); 
        [path, raw_name, ext] = fileparts(raw_filename_fullfile); 
        fprintf(fileID,'$nmeach\n');
        fprintf(fileID,'filraw=''%s%s''\n', basis_input_dir, raw_filename);
        fprintf(fileID,'metabo=''%s''\n', raw_name);
        fprintf(fileID,'degzer=%s\n', zero_order_ph);
        fprintf(fileID,'degppm=%s\n', first_order_ph);
        fprintf(fileID,'conc=%s\n', scaling);
        fprintf(fileID,'NOSHIF=.%s\n', noshift_tf);
        fprintf(fileID,'ppmapp=%s, %s\n', ppm_app_higher, ppm_app_lower);
        fprintf(fileID,'$end\n');
        fprintf(fileID,'\n');
    end
    
    % Finish the makebasis.in writing process 
    fclose(fileID);
    close all; 

end
