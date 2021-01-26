%Parametric method
function [tau_values, A] = method_parametric_2(K_records, M_samples, X, Y, P_value, rho_values)
%x_segmented and y_segmented must have the same size
%M is the lenght of the segment (segleng)
%X and Y are the signals to be compare in the cross-bispectrum.
%P_value is the maximum number of samples to take in to the account for the
%delay window. It creates an array from -P  to P in steps of 1
%rho_values is a list of values for rho. For example [-1, 1, 0]

total_segments = K_records - 1;
segshift = M_samples/2;
%frequencies = M_samples;
tau_values = -P_value: 1: P_value;
t = length(tau_values);
p = length(rho_values);

mask = toeplitz(-P_value:P_value);
mask2 = -P_value : P_value;
%k = P_value + 1 : M_samples - P_value;
R_xxx = zeros(t, t, p);
R_xyx = zeros(t, p);

for rr = 1 : p
    rho = rho_values(rr);
    for ii=1:total_segments
    
        n_segment_X = X((ii-1) * segshift +1 : (ii+1) * segshift);
        n_segment_Y = Y((ii-1) * segshift +1 : (ii+1) * segshift);
        rxxx = zeros(t, t);
        rxyx = zeros(t, 1);
        
        for kk = P_value + 1 : M_samples - P_value            
            rxxx = rxxx + n_segment_X(kk) * n_segment_X(kk + mask) * n_segment_X(kk + rho);
            %[size(n_segment_X(kk) * n_segment_Y(kk + mask2) * n_segment_X(kk + rho))]
            rxyx = rxyx + n_segment_X(kk) * n_segment_Y(kk + mask2) * n_segment_X(kk + rho);
           
        end
      
        
        R_xxx(:,:,rr) = R_xxx(:,:,rr) + rxxx / M_samples;
        R_xyx(:,rr) = R_xyx(:,rr) + rxyx / M_samples;
      
    end
    R_xxx(:,:, rr) = R_xxx(:,:, rr) / total_segments;
    R_xyx(:, rr) = R_xyx(:, rr) / total_segments;
    
end

R_xxx = cat(1, R_xxx(:,:, 1), R_xxx(:,:, 2), R_xxx(:,:, 3));
R_xyx = cat(1, R_xyx(:, 1), R_xyx(:, 2), R_xyx(:, 3));

A = inv(R_xxx.' * R_xxx) * R_xxx.' * R_xyx;
A = fftshift(A);



end
