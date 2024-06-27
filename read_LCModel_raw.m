function fid_complex = read_LCModel_raw(filepath)
% 
% Read LCModel RAW files into MATLAB for your processing enjoyment 
% 
% Inputs: 
% 
%  1. filepath: Full file path to LCModel RAW file to be loaded
% 
% Outputs: 
% 
%  1. fid_complex: n x 1 complex vector for FID inside LCModel RAW
% 
% Example usage: read_LCModelRAW('F:\BasisSets\SummedSpins_for_MARSSinput\Hom.raw')
% 
% Author: Kelley Swanberg (Lund University, 2024) 
% 
% Written and tested in MATLAB 2023b

    % Load the string data into a variable and extract the good part 
    rawfile = fileread(filepath);
    rawfile_split = strsplit(rawfile, '$END'); 
    
    % Convert FID into array and turn into a single complex vector 
    file_fid = str2num(rawfile_split{1,2}); 
    file_fid_size = size(file_fid)
    file_fid_length = file_fid_size(1) * file_fid_size(2) / 2; % Assumes complex vector 
    file_fid_reshaped = reshape(file_fid', 2, [])'; 
    fid_complex = complex(file_fid_reshaped(:,1), file_fid_reshaped(:,2)); 
    
    % Uncomment to view spectrum 
    %plot(real(fftshift(fft(fid_complex))))

end 