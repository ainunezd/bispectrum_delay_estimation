function[X,Y] = create_signals_2(N, delay, signal, snr)
%Create signals X and Y from a common one but delayed with respect to each
%other
%N: Number of samples
%delay: Number of samples to be delayed
%Common signal
%snr: Siganl to noise ratio for the additive uncorrelated Gaussian noise

S = transpose(signal(1:(N+delay)));
S = S/std(S);
S = S - mean(S);
[~, argmax] = max(S);
S(argmax) = 2 * S(argmax);

time = 0:1:N;
W1 = wgn(length(time)-1, 1, 0);
W2 = wgn(length(time)-1, 1, 0);

Y = snr * S(1:N) + (1-snr) * W1;
X = snr * (S(delay+1:(N+delay))) + (1-snr) * W2;

end