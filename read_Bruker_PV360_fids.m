%% Convert Bruker PV 360 FIDs to structs readable by MATLAB and INSPECTOR

% Read Bruker FID file 
fid_to_open = 'fid_proc.64'; % File to convert 
larmor_freq = 3.995162426879109e+02; % MHz; method file $PVM_FrqWork
rx_bw = 8333.33333; % Hz; method file $PVM_SpecSWH
grpdly_shift = 76; % acqus file $GRPDLY whole number
grpdly_phase_value = .0833333 % acqus file $GRPDLY decimal
zero_order_phase = 20.3; % Degrees; specific to your data
ppm_calib = 2.7; % ppm 
plot_spectrum = 1; 

% Open reference file. We know this is a little-endian representation of
% 64-bit floating-point values 
fid = fopen(fid_to_open, 'rb', 'ieee-le');   % little-endian

% Read the whole file as doubles 
raw = fread(fid, Inf, 'double=>double');
fclose(fid);

% Reshape into a complex FID
raw = reshape(raw, 2, []).';
fid_complex = complex(raw(:,1), raw(:,2));
length_fid_complex = length(fid_complex); 

% Prepare ppm axis for spectral visualization
dwell_time = 1/rx_bw; 
bw_ppm = rx_bw / larmor_freq; 
ppm_axis_vector_st = ppm_calib + bw_ppm/2; 
ppm_axis_vector_end = ppm_calib - bw_ppm/2; 
ppm_axis_vector = linspace(ppm_axis_vector_st, ppm_axis_vector_end, length_fid_complex)'; 

% Prepare time-domain vector for use in fitting
ppm_axis_vector_unprocessed = linspace(ppm_axis_vector_st, ppm_axis_vector_end, length_fid_complex); 

% Correct for Bruker's group delay
grpdly_phase = grpdly_phase_value*360*pi/180;
length_fid_complex_phase_corrected = length_fid_complex - grpdly_shift; 
grpdly_phase_vector = linspace(0, grpdly_phase, length_fid_complex_phase_corrected); 
fid_complex_shifted = fid_complex(grpdly_shift+1:end); 
fid_complex_phase_corrected = ifft(ifftshift(fftshift(fft(fid_complex_shifted)) .* exp(j*grpdly_phase_vector'))); 
fid_complex_phase_corrected_zf = padarray(fid_complex_phase_corrected, grpdly_shift, 'post');  
spec_to_plot = real(exp(j*zero_order_phase)*fftshift(fft(fid_complex_phase_corrected_zf))); 

% Sanity check 
if plot_spectrum==1
       figure(1); 
       plot(ppm_axis_vector, flipud(spec_to_plot), 'Color', 'k', 'DisplayName','Spectrum');
       set(gca, 'XDir','reverse'); 
       xlim('tight'); 
       xlabel('Chemical shift (ppm)');
       ylabel('Signal (a.u.)');  
       legend('Location','southoutside'); 
end

% Put data into struct for reading by MATLAB and INSPECTOR
exptDat.sf = larmor_freq; 
exptDat.sw_h = rx_bw; 
exptDat.nspecC = length_fid_complex; 
exptDat.fid = fid_complex_phase_corrected_zf; 

% Save end result 
save('fid.mat', 'exptDat'); 