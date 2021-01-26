function [tau, T] = method_3(B_xxx,B_xyx, maxfreqbins)
%Method 3 from Nikias(1988) paper.
%B_xxx is the cross-bispectrum from x, x, x
%B_xyx is the cross-bispectrum from x, y, x
%maxfreqbins is the total number of samples in the fourier transform array.
%It is equivalent to M_samples from the bispectrum function.

I = B_xyx ./ B_xxx;

tau = -maxfreqbins/2 : maxfreqbins/2 - 1;

T = sum(I, 1); %Sum over lambda2 or omega2
T = abs(fftshift(ifft(T)));

end