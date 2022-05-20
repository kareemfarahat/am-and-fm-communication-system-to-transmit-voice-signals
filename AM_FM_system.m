
%% Read in the file
[m, Fs] = audioread('test_message.mp3');
T = 1 / Fs;             
L = length(m);          
t = (0:L-1) * T;        
%% Create object for playing audio
pl = audioplayer(m,Fs);  % original signal
%pl.play;
%% Plot audio 
N = size(m, 1);
figure;  stem(t, m);
title('message Time-domain');  xlabel('time(seconds)');
%% Plot the spectrum
dm = Fs / N;
w = (-(N/2):(N/2)-1)*dm;
y = fft(m) / N;         % For normalizing
y2 = fftshift(y);
figure;  plot(w, abs(y2));
title(' Amplitude Spectrum');  xlabel('Frequency(Hz)');
%% filtering at 4K Hz (Low pass filter)
lowpass(m,4000,Fs);
filter_m = lowpass(m,4000,Fs);
%% the average error between the original signal and the band limited signal.
err = immse(m,filter_m);
disp("the error =")
disp (err)
%% AM modulation of the recorded signal 
AM_modulate_signal = ammod(m,Fs,1000000);
%% Spectrum of modulated signal
spectrumAM = fft(AM_modulate_signal);
lengthOfSignal = length(m);
normalizedSpectrumAm= abs(spectrumAM/lengthOfSignal);
frequencySpectrumAm= normalizedSpectrumAm(1:lengthOfSignal/2+1);
frequencySpectrumAm(2:end-1) = 2*frequencySpectrumAm(2:end-1);
figure;  plot(frequencySpectrumAm);
title(' Spectrum of AM signal');  xlabel('Frequency(Hz)');
%% Envelope Detection based on Hilbert Transform and then FFT
envelope = abs (hilbert (AM_modulate_signal));
figure;  plot(envelope);
title(' Spectrum of AM signal');
%% the average error between the original signal and the transmitted message.
err = immse(m,envelope);
disp("the error =")
disp (err)
%% signal with noise
Noise_signal = awgn(m,10);
%% AM modulation of the recorded signal + Noise 
AM_modulate_noise_signal = ammod(Noise_signal,Fs,1000000);
%% Spectrum of modulated signal
spectrumAM = fft(AM_modulate_noise_signal);
lengthOfSignal = length(m);
normalizedSpectrumAm= abs(spectrumAM/lengthOfSignal);
frequencySpectrumAm= normalizedSpectrumAm(1:lengthOfSignal/2+1);
frequencySpectrumAm(2:end-1) = 2*frequencySpectrumAm(2:end-1);
figure;  plot(frequencySpectrumAm);
title(' Spectrum of AM noise signal');  xlabel('Frequency(Hz)');
%% Envelope Detection based on Hilbert Transform and then FFT for noise signal
envelope = abs (hilbert (AM_modulate_noise_signal));
figure;  plot(envelope);
title('AM noise signal');
%% the average error between the original signal and the transmitted message.
err = immse(m,envelope);
disp("the error =")
disp (err)
%% FM modulatiom
%% Read in the file
[m, Fs] = audioread('test_message.mp3');
T = 1 / Fs;             
L = length(m);          
t = (0:L-1) * T;        
%% Create object for playing audio
pl = audioplayer(m,Fs);  % original signal
%pl.play;
%% Plot audio 
N = size(m, 1);
figure;  stem(t, m);
title('message Time-domain');  xlabel('time(seconds)');
%% Plot the spectrum
dm = Fs / N;
w = (-(N/2):(N/2)-1)*dm;
y = fft(m) / N;         % For normalizing
y2 = fftshift(y);
figure;  plot(w, abs(y2));
title(' Amplitude Spectrum');  xlabel('Frequency(Hz)');
%% filtering at 4K Hz (Low pass filter)
lowpass(m,4000,Fs);
filter_m = lowpass(m,4000,Fs);
%% the average error between the original signal and the band limited signal.
err = immse(m,filter_m);
disp("the error =")
disp (err)
%% FM modulation of the recorded signal 
FM_modulate_signal = fmmod(m,Fs,1000000,2);
%% Spectrum of modulated signal
spectrumFM = fft(FM_modulate_signal);
lengthOfSignal = length(m);
normalizedSpectrumFm= abs(spectrumFM/lengthOfSignal);
frequencySpectrumFm= normalizedSpectrumFm(1:lengthOfSignal/2+1);
frequencySpectrumFm(2:end-1) = 2*frequencySpectrumFm(2:end-1);
figure;  plot(frequencySpectrumFm);
title(' Spectrum of FM signal');  xlabel('Frequency(Hz)');
%% Envelope Detection based on Hilbert Transform and then FFT
envelope = abs (hilbert (FM_modulate_signal));
figure;  plot(envelope);
title(' Spectrum of FM signal');
%% with noise
Noise_signal = awgn(m,10);
%% AM modulation of the recorded signal + Noise 
FM_modulate_noise_signal = fmmod(Noise_signal,Fs,1000000,2);
%% Spectrum of modulated signal
spectrumFM = fft(AM_modulate_noise_signal);
lengthOfSignal = length(m);
normalizedSpectrumFm= abs(spectrumFM/lengthOfSignal);
frequencySpectrumFm= normalizedSpectrumFm(1:lengthOfSignal/2+1);
frequencySpectrumFm(2:end-1) = 2*frequencySpectrumFm(2:end-1);
figure;  plot(frequencySpectrumFm);
title(' Spectrum of FM noise signal');  xlabel('Frequency(Hz)');
%% Envelope Detection based on Hilbert Transform and then FFT for noise signal
envelope = abs (hilbert (AM_modulate_noise_signal));
figure;  plot(envelope);
title('FM noise signal');

