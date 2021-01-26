function [tau, T] = method_4(B_xxx,B_xyx, B_yyy, maxfreqbins)
%Method 4 from Nikias(1988) paper.
%B_xxx is the cross-bispectrum from x, x, x
%B_xyx is the cross-bispectrum from x, y, x
%maxfreqbins is the total number of samples in the fourier transform array.
%It is equivalent to M_samples from the bispectrum function.

phi = angle(B_xyx)- 0.5 * (angle(B_xxx) + angle(B_yyy));
I = (abs(B_xyx) .* exp(j*phi)) ./ (sqrt(abs(B_xxx)) .* sqrt(abs(B_yyy)));

tau = -maxfreqbins/2 : maxfreqbins/2 - 1;

T = sum(I, 1); %Sum over lambda2 or omega2
T = abs(fftshift(ifft(T)));

end
