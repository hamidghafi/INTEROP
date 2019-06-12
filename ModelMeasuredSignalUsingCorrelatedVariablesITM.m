%   In this code, firstly, measured IQ data is used to find the interference
%   signals. Then, inverse transform method is employed to find the
%   statistics of interference signals. This statistics (CDF function of randome variable)
%   are used to model inteference. The performance of the model is assessed using BER calculation.
%% load IQ data
% Y.mat matlab variable consists of measured IQ data variables and
% properties. XDelta is the sampling period and Y is the vector of measured
% IQ data.
addpath(genpath(cd));
load  Y.mat;
fSample = 1/XDelta;
% Truncate the signal to reduce memory usage.
processTime = 600e-3;                                                      % time period (length) of processing
signalLength = fix(processTime*fSample);                                   % number of indexes corresponding to processing time
% Pick a fraction of measured signal
nSlot = 6;                                                                % nth slot of the measured signal
IQSignal = Y((nSlot-1)*signalLength+1:(nSlot)*signalLength);
clear Y

%% Resampling and Filtering
fSampleNew = 10e6;                                                         % new sampling frequency
% Do upsampling and downsampling 
[n,d] = rat(fSampleNew /fSample,1e-6);
IQSignalRes = resample(double(IQSignal).',n,d);
% Design low pass filter 
N = 96;                                                                    % FIR filter order
Bw = 2e6;                                                                  % low pass filter bandwidth
eqNum = fir1(N,Bw/fSampleNew);                                             % eqnum = vec of coeffs
% Apply filtering               
IQSignalFilt = conv(IQSignalRes, eqNum/sum(eqNum),'same');                 % filter signals out of BLE band.
% Find noise floor of filtered signal 
noisePow = var(IQSignalFilt(1:10000));

%% Find the Interferences 
threshold   = sqrt(2*fSampleNew/Bw*noisePow);                              % sets a threshold to detect interference signals (3dB above the noise floor)
% Find the starting index of first interference signal
startIndex = find(abs(IQSignalFilt)>threshold,1);                          
intParams = int_detect(IQSignalFilt,fSampleNew,threshold);

%% Inverse Transform Method
% Use inverse transform method to find an equivalent signal for
% interference.
% Calculate CDF of random variables. pulseWidth, ampLevel, tau are supposed
% to be random variables.

%% Genrate Equivalent Signal
eqSignal = generate_CITM_eq_signal(intParams,fSampleNew,noisePow,processTime);
eqSignal = (sqrt(fSampleNew/Bw))*eqSignal; 
eqSignalFilt = conv(eqSignal, eqNum/sum(eqNum),'same');                    % filter signals out of BLE band.
intSignals = [IQSignalFilt;eqSignalFilt];

%% Generate Plot
gainListdB = -110:1:-58;
BER = compareBER(intSignals,fSampleNew,gainListdB);
figure; 
semilogy(gainListdB, max(BER(:, 1),eps), 'linewidth', 1)
hold on;
semilogy(gainListdB, max(BER(:,2), eps), '*', 'linewidth', 1)
ylim([0.0001 0.1])
xlabel('Gain (dB)', 'Interpreter', 'latex')
ylabel('BER', 'Interpreter', 'latex')
legendCell = {'Measurement', 'Correlated Variables Inverse Transform Method'};
legend(legendCell, 'Location', 'southwest', 'Interpreter', 'latex')
grid on

