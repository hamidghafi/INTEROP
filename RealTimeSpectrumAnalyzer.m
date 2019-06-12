%% load Measurement data
% Y.mat matlab variable consists of measured IQ data variables and
% properties. XDelta is the sampling period and Y is the vector of measured
% IQ data.
addpath('D:\HAMIAD_FILES\Measurement')
load Y.mat;
fSample = 1/XDelta;                                                        % sampling frequency of measured data
% Truncate the signal to reduce memory usage.
processTime = 200e-3;                                                      % time period (length) of processing
signalLength_index = fix(processTime*fSample);                             % number of indexes corresponding to processing time
IQSignal = Y(1:signalLength_index);
clear Y

%% Real Time Spectrum Analyzer Parameter Setting 
acqTime_list = [40e-6 100e-6 200e-6 500e-6];                               % acquisition time
fftLength = 1024;                                                          % number of samples used in fft processing 
numOfFreqPoints = 32;
freq = linspace(0,fSample,numOfFreqPoints);

%% Equivalent Signal Reconstruction
% Reconstruct signal for different time acquisition values
reconsSignal = zeros(length(acqTime_list),signalLength_index);             

for i = 1:length(acqTime_list)
    
    acqTime = acqTime_list(i);
    % Call spectrum function to derive spectrum of the measured signal
    [exactAcqTime, eqSpectrum] = spectrogram(IQSignal,fSample,acqTime,fftLength,'spectrumtype', 'maximum');
    eqSpectrumReduced = [freq.' fftLength/numOfFreqPoints*spline(eqSpectrum(:,1).',eqSpectrum(:,2:end).',freq).'];
    % Reconstruct equivalent signal
    y = eqsignal(eqSpectrumReduced,exactAcqTime);
    calFactor = 1;
    y = calFactor*y;
    reconsSignal(i,1:length(y)) = y;
end

%% Resalmpling Measured and Modeled Signal (after filtering)
intSignals = [IQSignal(1:signalLength_index).';reconsSignal];
% Define new smaller sampling frequency. This reduces the rest simulation
% time.
fSampleNew = 10e6;                                                         % new sampling frequency 
[n,d] = rat(fSampleNew/(1/XDelta), 1e-6);                                  % upsampling and downsampling factors.
resampleSignal = resample(double(intSignals).',n,d).';                     % change the interference signal sampling rate to new defined sampling rate

%% BER calculation
gainListdB = -100:1:-58;
BER = calculateBER(resampleSignal,fSampleNew,gainListdB);

%% Generate Plots 
% BER plot
figure; 
semilogy(gainListdB, max(BER(:, 1),eps), 'linewidth', 1)
hold on;
semilogy(gainListdB, max(BER(:,2), eps), '*', 'linewidth', 1)
hold on;
semilogy(gainListdB, max(BER(:,3),eps), '-.', 'linewidth',1)
hold on;
semilogy(gainListdB, max(BER(:,4),eps), '--', 'linewidth',1)
hold on;
semilogy(gainListdB, max(BER(:,5),eps), '-o', 'linewidth',1)
ylim([0.0001 0.1])
xlabel('Gain (dB)', 'Interpreter', 'latex')
ylabel('BER', 'Interpreter', 'latex')
legendCell = {'Measurement', 'Acqtime = 40 $\mu$s', ...
    'Acqtime = 100 $\mu$s', 'Acqtime = 200 $\mu$s', 'Acqtime = 500 $\mu$s'};
legend(legendCell, 'Location', 'southwest', 'Interpreter', 'latex')
grid on

