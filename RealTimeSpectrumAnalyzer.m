%% load Measurement data
% Y.mat matlab variable consists of measured IQ data variables and
% properties. XDelta is the sampling period and Y is the vector of measured
% IQ data.
addpath('D:\Messung Mediamarkt Graz\Recordings')
addpath('P:\PRJ EFRE Interreg InterOP\HK\HK_Codes\Measurement Analysis\AUXFUNCTIONS')
addpath('P:\PRJ EFRE Interreg InterOP\HK\HK_Codes\BLE_Functions')
tic;
fprintf('loading data...\n');
load  20170704T182310755.mat;
fprintf('data loaded. loading time %f s\n',toc);
tic
fprintf('calculating equivalent signal ...\n');
fSample = 1/XDelta;                                                        % sampling frequency of measured data
% Truncate the signal to reduce memory usage.
processTime = 200e-3;                                                      % time period (length) of processing
signalLength = fix(processTime*fSample);                             % number of indexes corresponding to processing time
% Pick a fraction of measured signal
nSlot = 6;                                                                % nth slot of the measured signal
IQSignal = Y((nSlot-1)*signalLength+1:(nSlot)*signalLength);
clear Y

%% Real Time Spectrum Analyzer Parameter Setting 
acqTime_list = [10e-6 40e-6 100e-6];                               % acquisition time
fftLength = 1024;                                                          % number of samples used in fft processing 
numOfFreqPoints = 1024;
freq = linspace(0,fSample,numOfFreqPoints);

%% Equivalent Signal Reconstruction
% Reconstruct signal for different time acquisition values
reconsSignal = zeros(length(acqTime_list),signalLength);             

for i = 1:length(acqTime_list)
    
    acqTime = acqTime_list(i);
    % Call spectrum function to derive spectrum of the measured signal
    [exactAcqTime, eqSpectrum] = spectrogram(IQSignal,fSample,acqTime,fftLength,'spectrumtype', 'maximum');
    eqSpectrumReduced = [freq.' fftLength/numOfFreqPoints*abs(spline(eqSpectrum(:,1).',eqSpectrum(:,2:end).',freq)).'];
    % Reconstruct equivalent signal
    y = generate_eq_signal(eqSpectrumReduced,exactAcqTime);
    calFactor = 1;
    y = calFactor*y;
    reconsSignal(i,1:length(y)) = y;
end

%% Resalmpling Measured and Modeled Signal (after filtering)
intSignals = [IQSignal.';reconsSignal];
% Define new smaller sampling frequency. This reduces the rest simulation
% time.
fSampleNew = 10e6;                                                         % new sampling frequency 
[n,d] = rat(fSampleNew/(1/XDelta), 1e-6);                                  % upsampling and downsampling factors.
resampleSignal = resample(double(intSignals).',n,d).';                     % change the interference signal sampling rate to new defined sampling rate
fprintf('equivalent signal was calculated. run time %f s\n',toc);
%% BER calculation
gainListdB = -100:1:-58;
BER = compareBER(resampleSignal,fSampleNew,gainListdB);

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

