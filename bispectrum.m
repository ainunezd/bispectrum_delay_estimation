%Bispectra Fast Fourier transform on each segment.
function [B_xxx, B_yyy, B_xyx, frequencies] = bispectrum(K_records, M_samples, X, Y)
%X and Y are the signals to be compare in the cross-bispectrum.
%x_segmented and y_segmented must have the same size.
%M is the lenght of the segment (segleng).
%K is the number of segments to be analyze. 
%It is set to have a 50% overlap (segshift).

total_segments = K_records - 1;
segshift = M_samples/2;
frequencies = M_samples;

B_xxx = zeros(frequencies, frequencies);
B_yyy = zeros(frequencies, frequencies);
B_xyx = zeros(frequencies, frequencies);

%Create a matrix for the last term of the bispectrum.
mask = hankel([1:frequencies],[frequencies,1:frequencies-1] ); 

for ii=1:total_segments

    n_segment_X = X((ii-1) * segshift +1 : (ii+1) * segshift);
    n_segment_Y = Y((ii-1) * segshift +1 : (ii+1) * segshift);
    
    X_i_lamb = fft(n_segment_X - mean(n_segment_X), frequencies) / M_samples; % The time series must have mean = 0-
    Y_i_lamb = fft(n_segment_Y - mean(n_segment_Y), frequencies) / M_samples;
    CX_i_lamb = conj(X_i_lamb);
    CY_i_lamb = conj(Y_i_lamb);
    

    B_xxx = B_xxx + (X_i_lamb * X_i_lamb.') .* CX_i_lamb(mask);
    B_yyy = B_yyy + (Y_i_lamb * Y_i_lamb.') .* CY_i_lamb(mask);
    B_xyx = B_xyx + (X_i_lamb * Y_i_lamb.') .* CX_i_lamb(mask);
    
end

%Divide by total number of segments to obtain the mean.
B_xxx = B_xxx / total_segments;
B_yyy = B_yyy / total_segments;
B_xyx = B_xyx / total_segments;

end