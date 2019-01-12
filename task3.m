close all
clear variables
clc

fs = 44100;
N = 4096;

G_db = [-15, -5, 5, 15];
fc = 3000;

%unit pulse
x_n1=[1, zeros(1,N-1)];

figure('Name','Shelving Filter');
for gain = 1:length(G_db)
    y_n1 = shelvingfilter(x_n1, G_db(gain), fc,1);

    mag_Y1 = abs(fft(y_n1,N));
    fvec_y1 = fs*(0:length(mag_Y1)-1)/length(mag_Y1);                    
    Y_log1 = 20*log10(mag_Y1);                                   
    semilogx(fvec_y1,Y_log1);
    hold on;
end
legend(['gain = ' num2str(G_db(1))],['gain = ' num2str(G_db(2))],['gain = ' num2str(G_db(3))],['gain = ' num2str(G_db(4))]);
xlabel('frequency in Hz');
ylabel('frequency response');
title('Different Gains to Shelving Filter');
hold off;

%% Peak Filter
clear variables
fs = 44100;
fc =1000;
N = 4096;
Q= 4;
G_db = [2, 4, 6,8];
%unit pulse
x_n1=[1, zeros(1,N-1)];

figure('Name','Peak Filter');
for gain = 1:length(G_db)
    y_n1 = peakfilter(x_n1, G_db(gain), fc,Q);

    mag_Y1 = abs(fft(y_n1,N));
    fvec_y1 = fs*(0:length(mag_Y1)-1)/length(mag_Y1);                    
    Y_log1 = 20*log10(mag_Y1);                                   
    semilogx(fvec_y1,Y_log1);
    hold on;
end
legend(['gain = ' num2str(G_db(1))],['gain = ' num2str(G_db(2))],['gain = ' num2str(G_db(3))],['gain = ' num2str(G_db(4))]);
xlabel('frequency in Hz');
ylabel('frequency response in dB');
title('Different Gains to Peak Filter');
hold off;

%% Parametric Filter
clc

[x_signal, fs] =audioread('./audio/drumloop.wav');
gainLP = 6; %Change value to boost or dampen
gainHP = 6; %Change value to boost or dampen
gainP1 = 6; %Change value to boost or dampen
gainP2 = 6; %Change value to boost or dampen
fc = 5000;

LP = shelvingfilter(x_signal,gainLP,fc,1);
HP=shelvingfilter(LP,gainHP,fc,-1);
Peak_1= peakfilter(HP, gainP1, fc , Q);
Peak_2= peakfilter(Peak_1, gainP2, fc , Q);

sound(Peak_2,fs);

%b)
[x_popsong, fs] =audioread('./audio/popsong.wav');
n = 0:length(x_popsong)-1;
f = 440;
x_noise = 0.5*sin(n*2*pi*f/fs);
x_noisy = x_popsong + x_noise';
%to remove it apply peak filter with fc 440
fc = 440;
gain = -10;
recover = peakfilter(x_noisy, gain, fc , Q);
%sound(recover,fs);

%c)
gain =6;
fc = 4000;
Q = 6;

x_vocals = audioread('./audio/vocals.wav');
[y_peak2, y_peak1, y_high, y_low] = parametricEqualiser(x_vocals, gain, fc , Q);
%y_low = shelvingfilter(x_vocals,gain,fc,1);
%y_high = shelvingfilter(y_low,gain,fc,-1);
%y_peak1 = peakfilter(y_high,gain,fc,Q);
%y_peak2 = peakfilter(y_peak1,gain,fc,Q);

%sound(x_vocals,fs);
%sound(y_low,fs);
%sound(y_high,fs);
%sound(y_peak1,fs);
%sound(y_peak2,fs);


%% Functions
function y_n = shelvingfilter(x_signal, gain, fc, HPF_or_LPH) %High = -1 Low = 1 %what is the use of fc?!?!
    fs = 44100;
    
    %linear Gain
    H0 = 10^(gain/20)-1;
    V0 = H0 + 1;
    wc = 2*pi*fc/fs;
    if HPF_or_LPH == -1 %High Pass Filter
        if gain >=0
            a = (1-tan(wc/2)) / (tan(wc/2)+1);   %boost
        else
            a = (1-V0*tan(wc/2)) / (1+V0*tan(wc/2));   %cut
        end
    else %LOW Pass Filter
        if gain >=0
            a = (1-tan(wc/2)) / (tan(wc/2)+1);   %boost
        else
            a = (V0-tan(wc/2)) / (V0+tan(wc/2));   %cut
        end
    end
    x_hist=0;
    y_hist=0;
    y_ap1 = zeros(1,length(x_signal));
    y_n = zeros(1,length(x_signal));
    for i=1:length(x_signal)
       y_ap1(i)=-a*x_signal(i) + x_hist + a*y_hist;
       %history values
       x_hist=x_signal(i);
       y_hist=y_ap1(i);
       %Full diagram
       y_n(i) = 0.5*H0*(x_signal(i)+HPF_or_LPH*y_ap1(i))+x_signal(i);
    end 
end

function y_n = peakfilter(x_signal, gain, fc , Q)
    fs = 44100;
    %linear Gain
    H0 = 10^(gain/20)-1;
    V0 = H0 + 1;
    wc = 2*pi*fc/fs;
    wb = wc/Q;
    b = cos(wb);
    if gain >=0
        a = (1-tan(wb/2)) / (tan(wb/2)+1);   %boost
    else
        a = (V0-tan(wb/2)) / (V0+tan(wb/2));   %cut
    end
    c = b*(1+a); 
    x_hist=[0, 0];
    y_hist=[0, 0];
    y_ap2 = zeros(1,length(x_signal));
    y_n = zeros(1,length(x_signal));
    for i=1:length(x_signal)
%       y_ap2(i)=-a*x_signal(i) + x_hist(2) + a*y_hist(2);
       y_ap2(i)= a*x_signal(i)-c*x_hist(1) +x_hist(2) + c*y_hist(1) - a*y_hist(2);
       %history values
       x_hist(2)=x_hist(1);
       y_hist(2)=y_hist(1);       
       x_hist(1)=x_signal(i);
       y_hist(1)=y_ap2(i);
       %Full diagram
       y_n(i) = 0.5*H0*(x_signal(i)-y_ap2(i))+x_signal(i);
    end 
end

function [P2, P1, LP, HP] = parametricEqualiser(x_signal, gain, fc , Q)
    LP = shelvingfilter(x_signal,gain,fc,1);
    HP=shelvingfilter(LP,gain,fc,-1);
    P1= peakfilter(HP, gain, fc , Q);
    P2= peakfilter(P1, gain, fc , Q);
end


