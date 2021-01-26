%% My own implementation for several runs or experiments and changing the methods
%For gaussian noise
for snr = 0.9
    [snr]
    N = 96000; %In samples
    %snr = 0.4;
    delay = 16;
    %[signal] = gamma_distribution(1000000, 100);
    signal = spikes_signal(896000:1024000)';
    records = 375;
    samples = 256;
    tau_value = 30;
    runs = 1;
    T_array = zeros(runs, samples);
    T_param = zeros(runs, tau_value *2 +1);

    for r = 1:runs
        disp(r)
        [X,Y] = create_signals_2(N, delay, signal, snr);
        [Bxxx, Byyy, Bxyx, frequencies] = bispectrum(records, samples, X, Y);
        [taus, T_array(r, :)] = method_3(Bxxx,Bxyx, frequencies); 
        %[taus, T_array(r, :)] = method_4(Bxxx,Bxyx,Byyy, frequencies); 
        %[taus_param, T_param(r, :)] = method_parametric(records, samples, X, Y, tau_value, [-1,0,1]);
        

    end


    fig = figure();
    plot(taus, mean(abs(T_array), 1))
    title(sprintf('Method 3. Delay %.0f. EEG with SNR = %.1f', delay, snr))
    %plot(taus_param, mean(abs(T_param), 1))
    %title(sprintf('Method parametric. Delay %.0f. SNR = %.1f.', delay, snr))
    xlim([-tau_value tau_value]);
    grid on
    xlabel('Tau (samples)')
    saveas(fig, sprintf('/home/ana/Documents/Lab rotation/Matlab_figure_ifft/EEG_method_3_%.0f_%.0f_delay%.0f_trials%.0f_snr%.2f.fig', records, samples, delay, runs, snr));
    %clc;clear;
end

%% My own implementation
%For correlated gaussian noise noise
for snr = 0.9:-0.1:0.1
    [snr]
    N = 18000; %In samples
    %snr = 0.4;
    delay = 16;
    [signal] = gamma_distribution(1000000, 100);
    records = 128;
    samples = 128;
    tau_value = 30;
    runs = 1;
    T_array = zeros(runs, samples);

    for r = 1:runs
        disp(r)
        [X,Y] = create_signals_correlated_noise(N, delay, signal, snr);
        %[Bxxx, Byyy, Bxyx, frequencies] = bispectrum(records, samples, X, Y);
        %[taus, T_array(r, :)] = method_3(Bxxx,Bxyx, frequencies); 
        %[taus, T_array(r, :)] = method_4(Bxxx,Bxyx,Byyy, frequencies); 
        [taus_param, T_param(r, :)] = method_parametric(records, samples, X, Y, tau_value, [-1,0,1]);


    end


    fig = figure();
    %plot(taus, mean(abs(T_array), 1))
    %title(sprintf('Method 3. Correlated noise. Delay %.0f. %.0f x %.0f. SNR = %.1f', delay, records, samples, snr))
    plot(taus_param, mean(abs(T_param), 1))
    title(sprintf('Method parametric. Correlated noise. Delay %.0f. SNR = %.1f.', delay, snr))
    xlim([-tau_value tau_value]);
    grid on
    xlabel('Tau (samples)')
    saveas(fig, sprintf('/home/ana/Documents/Lab rotation/Matlab_figure_ifft/method_param_correlated_noise_%.0f_%.0f_delay%.0f_trials%.0f_snr%.2f.fig', records, samples, delay, runs, snr));
    %clc;clear;
end