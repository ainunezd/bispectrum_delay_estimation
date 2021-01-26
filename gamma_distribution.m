function [ISI_matrix] = gamma_distribution(N, fs)
%Create ISI_matrix signal from a gamma distribution
%N is the number of samples for the vector
%fs is the sampling rate asssuming to be 100 from firing rate and the
%ISI_matrix are the random interspike intervals 

gamma = 2;
rand_numbers = rand(N*10, 1);
%Values in ms for interspike interval
ISIS = (-log(rand_numbers) / fs) * 1000; 
spikes_times = cumsum(ISIS);
indices = gamma-1 : gamma : gamma * N ; 
spikes = spikes_times(indices);
ISI_matrix = transpose(diff(spikes));


end
