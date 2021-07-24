%All plots are commented Except last one.
%Some lines of codes were commented due to being unnecessary or certain
%conditions and questions.

%I

%The signals

[y,fs] = audioread('Short_BBCArabic2.wav'); [y2,fs2] = audioread('Short_QuranPalestine.wav');
%1st message is BBC 'x1'
%Here I split the 2 cols then add them into one
X1 = y(:,1); X2 = y(:,2); x1 = X1+X2;
%2nd Message
Z1 = y2(:,1); Z2 = y2(:,2); x2 = Z1+Z2;
%Getting the length of both
N1 = length(x1);
N2 = length(x2);
%To know which message will be padded with zero *Smaller one*
if(N1 > N2)
    Max = N1;
    Min = N2;
    zero1 = zeros(Max-Min,1);
    zero1 = zero1';
    x2 = x2';
    x2 = [x2 zero1];
    x2 = x2';
else
    Max = N2;
    Min = N1;
    zero1 = zeros(Max-Min,1);
    zero1 = zero1';
    x1 = x1';
    x1 = [x1 zero1];
    x1 = x1';
end
%interping at first since I know it will need interp anyway
x1 = interp(x1,10);
x2 = interp(x2,10);
N = length(x1);
X = fft(x1);
Z = fft(x2);
%to adjust the x axis
k = -N/2:N/2-1;
% subplot(2,1,1); plot(k*fs/N, fftshift(abs(X))); title('1st Msg');
% subplot(2,1,2); plot(k*fs/N, fftshift(abs(Z))); title('2nd Msg');

%II

%The AM modulator

fn1 = 100*10^3;
fn2 = 150*10^3;
ts = 1/fs;
n = 1:N;
%Since interp is used, like the project slides
%Fs is multiplied by 10, while ts is divided by 10
carrier1 = cos(2*pi*fn1*n*ts/10)';
carrier2 = cos(2*pi*fn2*n*ts/10)';
%modulated message 1,2
modmsg1 = x1.*carrier1;
modmsg2 = x2.*carrier2;
modmsg1 = fft(modmsg1);
modmsg2 = fft(modmsg2);
%Multiplixing both
modmsg = modmsg1 + modmsg2;
%The spectrum of the output of the transmitter
% subplot(3,1,1); plot(k*10*fs/N, fftshift(abs(modmsg1))); title('1st Modulated Msg'); xlabel('Frequency (Hz)'); ylabel('Amplitude');
% subplot(3,1,2); plot(k*10*fs/N, fftshift(abs(modmsg2))); title('2nd Modulated Msg'); xlabel('Frequency (Hz)'); ylabel('Amplitude');
% subplot(3,1,3); plot(k*10*fs/N, fftshift(abs(modmsg))); title('Both Mod Msgs'); xlabel('Frequency (Hz)'); ylabel('Amplitude');

%III

%The wireless channel
%nothing

user = input("Please enter 1 or 2 for the desired channel : ");

if ( user == 1)
    
%IV
%The RF
%Start of the band passfilter

 F_pass1 = fn1-10000; F_stop1 = fn1-11000;
 F_pass2 = fn1+10000; F_stop2 = fn1+11000;
 A_stop2 = 60; A_pass = 1; A_stop1 = 60;
 % d = fdesign.bandpass;
 BandPassSpecObj = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
  		F_stop1, F_pass1, F_pass2, F_stop2, A_stop1, A_pass, ...
  		A_stop2, 10*fs);
 BandPassFilt = design(BandPassSpecObj, 'equiripple');
 %The next 3 line of codes are not so useful
 designmethods(BandPassSpecObj,'SystemObject',true);
 BandPassFiltSysObj = design(BandPassSpecObj,...
 'equiripple','SystemObject',true);
%  fvtool(BandPassFilt); %plotting the Band Pass Filter
% End
% Filter function takes the modulated msg in time domain only
 modmsgtime = ifft(modmsg);
% Y is just the signal name that's all
Y = filter(BandPassFilt, modmsgtime);
% This is the new modulated msg
Y = fft(Y);
% the output of the RF filter (before the mixer)
% subplot(1,1,1); plot(k*10*fs/N, fftshift(abs(Y))); title('Filtered Msg at WN'); xlabel('Frequency (Hz)'); ylabel('Amplitude');

%V

%The Oscillator

%Just like Lec5 

wif = 25*10^3;
wosc1 = fn1 + wif;

% wosc2 = fn2 + wif;
%Testing for shifting
% wosc1 = wosc1 + 1000;
% wosc1 = wosc1 + 100;

carrierLO1 = cos(2*pi*wosc1*n*ts/10)';

% carrierLO2 = cos(2*pi*wosc2*n*ts/10)';
%to test without RF, replace Y with modmsg
% Ytimedomain = ifft(modmsg);

Ytimedomain = ifft(Y);
Ynew = Ytimedomain.*carrierLO1;
Ynew = fft(Ynew);
%The output of the mixer
% subplot(1,1,1); plot(k*10*fs/N, fftshift(abs(Ynew))); title('Modulated Msg at Osc'); xlabel('Frequency (Hz)'); ylabel('Amplitude');

%VI

%The IF stage

%Start of the band passfilter, same as previous one

F_pass3 = wif-10000; F_stop3 = wif-15000;
F_pass4 = wif+10000; F_stop4 = wif+15000;
A_stop4 = 60; A_pass2 = 1; A_stop3 = 60;

% d1 = fdesign.bandpass;

BandPassSpecObj1 = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
		F_stop3, F_pass3, F_pass4, F_stop4, A_stop3, A_pass2, ...
		A_stop4, 10*fs);
BandPassFilter = design(BandPassSpecObj1, 'equiripple');
designmethods(BandPassSpecObj1,'SystemObject',true);
BandPassFiltSysObj1 = design(BandPassSpecObj1,...
'equiripple','SystemObject',true);
% fvtool(BandPassFilter);
%End 
modmsgIF = ifft(Ynew);
YIF = filter(BandPassFilter,modmsgIF);
%This is the new modulated msg
YIF = fft(YIF);
%Output of the IF filter
% subplot(1,1,1); plot(k*10*fs/N, fftshift(abs(YIF))); title('Filtered Msg at WIF'); xlabel('Frequency (Hz)'); ylabel('Amplitude');


%VII

%The Baseband detection

%Creating a carrier in the IF 
carrierIF = cos(2*pi*wif*n*ts/10)';
% subplot(2,1,1); plot(k*10*fs/N, fftshift(abs(fft(carrierIF))));
% subplot(2,1,1); plot(k*10*fs/N, fftshift(abs(YIF)))
YNEWTD = ifft(YIF);
YNEW = YNEWTD.*carrierIF;
YNEW = fft(YNEW);
%Output of the mixer (before the LPF)
% subplot(1,1,1); plot(k*10*fs/N, fftshift(abs(YNEW))); title('Mixer O/p'); xlabel('Frequency (Hz)'); ylabel('Amplitude');
%Start of the lowpass filter
FP = 10000; FS = 11000;
AP = 1; AS = 60;
designSpecs = fdesign.lowpass('Fp,Fst,Ap,Ast', FP, FS, AP, AS, fs*10);
lpFilter = design(designSpecs,'equiripple');
% fvtool(lpFilter)

%Last Y Time Domain
LASTYTD = ifft(YNEW);
LastY = filter(lpFilter,LASTYTD);
LastY = fft(LastY);
%Output of the LPF
subplot(1,1,1); plot(k*10*fs/N, fftshift(abs(LastY))); title('Final Msg - BBC'); xlabel('Frequency (Hz)'); ylabel('Amplitude');
%For testing audio
LastY2 = ifft(LastY);
%by the same concept of interp, decimate is the inverse
LastY2 = decimate(LastY2,10);
sound(LastY2,fs);
% pause(17);
%For original audio
% sound(y,fs);


elseif(user == 2)
%IV

%The RF
%Start of the band passfilter 
    
    
F_pass1 = fn2-10000; F_stop1 = fn2-11000;
F_pass2 = fn2+10000; F_stop2 = fn2+11000;
A_stop2 = 60; A_pass = 1; A_stop1 = 60;
% d = fdesign.bandpass;
BandPassSpecObj = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
 		F_stop1, F_pass1, F_pass2, F_stop2, A_stop1, A_pass, ...
 		A_stop2, 10*fs);
BandPassFilt = design(BandPassSpecObj, 'equiripple');
%The next 3 line of codes are not so useful
designmethods(BandPassSpecObj,'SystemObject',true);
BandPassFiltSysObj = design(BandPassSpecObj,...
'equiripple','SystemObject',true);
% fvtool(BandPassFilt); %plotting the Band Pass Filter
%End
%Filter function takes the modulated msg in time domain only
modmsgtime = ifft(modmsg);
%Y is just the signal name that's all
Y = filter(BandPassFilt, modmsgtime);
%This is the new modulated msg
Y = fft(Y);

% subplot(1,1,1); plot(k*10*fs/N, fftshift(abs(Y))); title('Filtered Msg at WN'); xlabel('Frequency (Hz)'); ylabel('Amplitude');


%V

%The Oscillator

%Just like Lec5 

wif = 25*10^3;
wosc1 = fn1 + wif;

% wosc2 = fn2 + wif;
%Testing for shifting
% wosc1 = wosc1 + 1000;
% wosc1 = wosc1 + 100;

carrierLO1 = cos(2*pi*wosc1*n*ts/10)';

% carrierLO2 = cos(2*pi*wosc2*n*ts/10)';
%to test without RF
% Ytimedomain = ifft(modmsg);

Ytimedomain = ifft(Y);
Ytimedomain2 = Ytimedomain;
Ynew = Ytimedomain.*carrierLO1;
Ynew = fft(Ynew);
% subplot(1,1,1); plot(k*10*fs/N, fftshift(abs(Ynew))); title('Modulated Msg at Osc'); xlabel('Frequency (Hz)'); ylabel('Amplitude');


%VI

%The IF stage

%Start of the band passfilter, same as previous one

F_pass3 = wif-10000; F_stop3 = wif-15000;
F_pass4 = wif+10000; F_stop4 = wif+15000;
A_stop4 = 60; A_pass2 = 1; A_stop3 = 60;
% d1 = fdesign.bandpass;
BandPassSpecObj1 = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
		F_stop3, F_pass3, F_pass4, F_stop4, A_stop3, A_pass2, ...
		A_stop4, 10*fs);
BandPassFilter = design(BandPassSpecObj1, 'equiripple');
designmethods(BandPassSpecObj1,'SystemObject',true);
BandPassFiltSysObj1 = design(BandPassSpecObj1,...
'equiripple','SystemObject',true);
% fvtool(BandPassFilter);
%End 
modmsgIF = ifft(Ynew);
YIF = filter(BandPassFilter,modmsgIF);
%This is the new modulated msg

YIF = fft(YIF);
% subplot(1,1,1); plot(k*10*fs/N, fftshift(abs(YIF))); title('Filtered Msg at WIF'); xlabel('Frequency (Hz)'); ylabel('Amplitude');


%VII

%The Baseband detection

%Creating a carrier in the IF 
carrierIF = cos(2*pi*wif*n*ts/10)';
% subplot(2,1,1); plot(k*10*fs/N, fftshift(abs(fft(carrierIF))));
% subplot(2,1,1); plot(k*10*fs/N, fftshift(abs(YIF)))
YNEWTD = ifft(YIF);
YNEW = YNEWTD.*carrierIF;
YNEW = fft(YNEW);
% subplot(1,1,1); plot(k*10*fs/N, fftshift(abs(YNEW))); title('Mixer O/p'); xlabel('Frequency (Hz)'); ylabel('Amplitude');
%Start of the lowpass filter
FP = 10000; FS = 11000;
AP = 1; AS = 60;
designSpecs = fdesign.lowpass('Fp,Fst,Ap,Ast', FP, FS, AP, AS, fs*10);
lpFilter = design(designSpecs,'equiripple');
% fvtool(lpFilter)
%Last Y Time Domain
LASTYTD = ifft(YNEW);
LastY = filter(lpFilter,LASTYTD);
LastY = fft(LastY);
subplot(1,1,1); plot(k*10*fs/N, fftshift(abs(LastY))); title('Final Msg - Short Quraan'); xlabel('Frequency (Hz)'); ylabel('Amplitude');
%For testing audio
LastY2 = ifft(LastY);
%by the same concept of interp, decimate is the inverse
LastY2 = decimate(LastY2,10);
sound(LastY2,fs);
% pause(17);
%For original audio
% sound(y2,fs);
end


%All done
