close all;
clear

%% loading

load data6.mat
Fs= 1000

%% plotting of ECG, PCG & Spectrogram

%% patient_1

fig1=figure
subplot(3,1,1)
plot(data(1).t(1000:3000),data(1).ECG(1000:3000))
xlabel('patient1 ECG')

subplot(3,1,2)
plot(data(1).t(1000:3000),data(1).PCG(1000:3000))
xlabel('patient1 PCG')

subplot(3,1,3)
spectrogram(data(1).PCG(1000:3000),50,45,[],'MinThreshold',-10,'yaxis')
xlabel('patient1 Spectrogram')


%% patient_2
fig2=figure
subplot(3,1,1)
plot(data(2).t(1500:4500),data(2).ECG(1500:4500))
xlabel('patient2 ECG')

subplot(3,1,2)
plot(data(2).t(1500:4500),data(2).PCG(1500:4500))
xlabel('patient2 PCG')

subplot(3,1,3)
spectrogram(data(2).PCG(1500:4500),50,45,[],'MinThreshold',-10,'yaxis')
xlabel('patient2 Spectrogram')


%% patient_3
fig3=figure
subplot(3,1,1)
plot(data(3).t(1000:2400),data(3).ECG(1000:2400))
xlabel('patient3 ECG')

subplot(3,1,2)
plot(data(3).t(1000:2400),data(3).PCG(1000:2400))
xlabel('patient3 PCG')

subplot(3,1,3)
spectrogram(data(3).PCG(1000:2400),50,45,[],'MinThreshold',-10,'yaxis')
xlabel('patient3 Spectrogram')


%% patient_4
fig4=figure
subplot(3,1,1)
plot(data(4).t(1000:1800),data(4).ECG(1000:1800))
xlabel('patient4 ECG')

subplot(3,1,2)
plot(data(4).t(1000:1800),data(4).PCG(1000:1800))
xlabel('patient4 PCG')

subplot(3,1,3)
spectrogram(data(4).PCG(1000:1800),50,45,[],'MinThreshold',-10,'yaxis')
xlabel('patient4 Spectrogram')


%% patient_5

fig5=figure
subplot(3,1,1)
plot(data(5).t(1000:3000),data(5).ECG(1000:3000))
xlabel('patient5 ECG')

subplot(3,1,2)
plot(data(5).t(1000:3000),data(5).PCG(1000:3000))
xlabel('patient5 PCG')

subplot(3,1,3)
spectrogram(data(5).PCG(1000:3000),50,45,[],'MinThreshold',-10,'yaxis')
xlabel('patient5 Spectrogram')


%% Part_4_ patient_1
 
ECG_new_P1=data(1).ECG-data(1).ECG(1);
 ECG_new_P1_res=resample(ECG_new_P1,200,1000);
% Low pass filter

b=1/32*[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a=[1 -2 1]
Low_pass_filter=filter(b,a,ECG_new_P1_res);

% Bandpass Filter

b_h=1/32*[-1 zeros(1,15) 32 -32 zeros(1,14) 1];  
a_h=[1 -1]
band_pass_filter=filter(b_h,a_h,Low_pass_filter); 

% Derivative filter

b_drv=1/8*[1 2 0 -2 -1]; 
a_drv=[1]
derivative_filter=filter(b_drv,a_drv,band_pass_filter);

% squaring
square_1=(derivative_filter.^2);

% integrating
a_intg=[1]
b_intg=1/30*ones(1,30);
Integrated=filter(b_intg,a_intg,square_1);
% detecting QRS
[peakOnsets_1, peakOffsets_1] = detectPeaks(Integrated);

peakOnsets_1_Original= 5*peakOnsets_1;


%% part_4_patient_2

ECG_new_P2=data(2).ECG-data(2).ECG(1);
ECG_new_P2_res=resample(ECG_new_P2,200,1000);
% Low pass filter
b=1/32*[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a=[1 -2 1]
Low_pass_filter_2=filter(b,a,ECG_new_P2_res);
% Bandpass Filter
b_h=1/32*[-1 zeros(1,15) 32 -32 zeros(1,14) 1];  
a_h=[1 -1]
band_pass_filter_2=filter(b_h,a_h,Low_pass_filter_2); 
% Derivative filter
b_drv=1/8*[1 2 0 -2 -1]; 
a_drv=[1]
derivative_filter_2=filter(b_drv,a_drv,band_pass_filter_2);
% squaring
square_2=(derivative_filter_2.^2);
% integrating
a_intg=[1]
b_intg=1/30*ones(1,30);
Integrated_2=filter(b_intg,a_intg,square_2);
% detecting QRS
[peakOnsets_2, peakOffsets_2] = detectPeaks(Integrated_2);

Original_ECG_2= 5*peakOnsets_2;
%% Part_4_patient_3

ECG_new_P3=data(3).ECG-data(3).ECG(1);
ECG_new_P3_res=resample(ECG_new_P3,200,1000);
% Low pass filter
b=1/32*[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a=[1 -2 1]
Low_pass_filter_3=filter(b,a,ECG_new_P3_res);
% Bandpass Filter
b_h=1/32*[-1 zeros(1,15) 32 -32 zeros(1,14) 1];  
a_h=[1 -1]
band_pass_filter_3=filter(b_h,a_h,Low_pass_filter_3); 
% Derivative filter
b_drv=1/8*[1 2 0 -2 -1]; 
a_drv=[1]
derivative_filter_3=filter(b_drv,a_drv,band_pass_filter_3);
% squaring
square_3=(derivative_filter_3.^2);
% integrating
a_intg=[1]
b_intg=1/30*ones(1,30);
Integrated_3=filter(b_intg,a_intg,square_3);
% detecting QRS
[peakOnsets_3, peakOffsets_3] = detectPeaks(Integrated_3);

Original_ECG_3= 5*peakOnsets_3;

%% part_4_patient_4

ECG_new_P4=data(4).ECG-data(4).ECG(1);
ECG_new_P4_res=resample(ECG_new_P4,200,1000);
% Low pass filter
b=1/32*[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a=[1 -2 1]
Low_pass_filter_4=filter(b,a,ECG_new_P4_res);
% Bandpass Filter
b_h=1/32*[-1 zeros(1,15) 32 -32 zeros(1,14) 1];  
a_h=[1 -1]
band_pass_filter_4=filter(b_h,a_h,Low_pass_filter_4); 
% Derivative filter
b_drv=1/8*[1 2 0 -2 -1]; 
a_drv=[1]
derivative_filter_4=filter(b_drv,a_drv,band_pass_filter_4);
% squaring
square_4=(derivative_filter_4.^2);
% integrating
a_intg=[1]
b_intg=1/30*ones(1,30);
Integrated_4=filter(b_intg,a_intg,square_4);
% detecting QRS
[peakOnsets_4, peakOffsets_4] = detectPeaks(Integrated_4);

Original_ECG_4= 5*peakOnsets_4;

%% part_4_patient_5

ECG_new_P5=data(5).ECG-data(5).ECG(1);
ECG_new_P5_res=resample(ECG_new_P4,200,1000);
% Low pass filter
b=1/32*[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a=[1 -2 1]
Low_pass_filter_5=filter(b,a,ECG_new_P5_res);
% Bandpass Filter
b_h=1/32*[-1 zeros(1,15) 32 -32 zeros(1,14) 1];  
a_h=[1 -1]
band_pass_filter_5=filter(b_h,a_h,Low_pass_filter_5); 
% Derivative filter
b_drv=1/8*[1 2 0 -2 -1]; 
a_drv=[1]
derivative_filter_5=filter(b_drv,a_drv,band_pass_filter_5);
% squaring
square_5=(derivative_filter_5.^2);
% integrating
a_intg=[1]
b_intg=1/30*ones(1,30);
Integrated_5=filter(b_intg,a_intg,square_5);
% detecting QRS
[peakOnsets_5, peakOffsets_5] = detectPeaks(Integrated_5);

Original_ECG_5= 5*peakOnsets_5;


%% patient_1

fig6=figure
subplot(2,1,1)
plot(data(1).t,data(1).ECG)

subplot(2,1,2)
plot(data(1).t,data(1).PCG)

subplot(2,1,1)
hold on
plot(data(1).t(peakOnsets_1_Original),data(1).ECG(peakOnsets_1_Original),'r*')
hold on
plot(data(1).t(peakOnsets_1_Original+300),data(1).ECG(peakOnsets_1_Original+300),'r*')

subplot(2,1,2)
hold on
plot(data(1).t(peakOnsets_1_Original),data(1).PCG(peakOnsets_1_Original),'r*')
hold on
plot(data(1).t(peakOnsets_1_Original+300),data(1).PCG(peakOnsets_1_Original+300),'r*')




%% patient_2
fig7=figure
subplot(2,1,1)
plot(data(2).t,data(2).ECG)

subplot(2,1,2)
plot(data(2).t,data(2).PCG)

subplot(2,1,1)
hold on
plot(data(2).t(Original_ECG_2),data(2).ECG(Original_ECG_2),'r*')
hold on
plot(data(2).t(Original_ECG_2+300),data(2).ECG(Original_ECG_2+300),'r*')

subplot(2,1,2)
hold on
plot(data(2).t(Original_ECG_2),data(2).PCG(Original_ECG_2),'r*')
hold on
plot(data(2).t(Original_ECG_2+300),data(2).PCG(Original_ECG_2+300),'r*')

%% patient_3
fig8=figure
subplot(2,1,1)
plot(data(3).t,data(3).ECG)

subplot(2,1,2)
plot(data(3).t,data(3).PCG)

subplot(2,1,1)
hold on
plot(data(3).t(Original_ECG_3),data(3).ECG(Original_ECG_3),'r*')
hold on
plot(data(3).t(Original_ECG_3+300),data(3).ECG(Original_ECG_3+300),'r*')

subplot(2,1,2)
hold on
plot(data(3).t(Original_ECG_3),data(3).PCG(Original_ECG_3),'r*')
hold on
plot(data(3).t(Original_ECG_3+300),data(3).PCG(Original_ECG_3+300),'r*')
%% patient_4

fig9=figure
subplot(2,1,1)
plot(data(4).t,data(4).ECG)

subplot(2,1,2)
plot(data(4).t,data(4).PCG)

subplot(2,1,1)
hold on
plot(data(4).t(Original_ECG_4),data(4).ECG(Original_ECG_4),'r*')
hold on
plot(data(4).t(Original_ECG_4+300),data(4).ECG(Original_ECG_4+300),'r*')

subplot(2,1,2)
hold on
plot(data(4).t(Original_ECG_4),data(4).PCG(Original_ECG_4),'r*')
hold on
plot(data(4).t(Original_ECG_4+300),data(4).PCG(Original_ECG_4+300),'r*')

%% patient_5
fig10=figure
subplot(2,1,1)
plot(data(5).t,data(5).ECG)

subplot(2,1,2)
plot(data(5).t,data(5).PCG)

subplot(2,1,1)
hold on
plot(data(5).t(Original_ECG_5),data(5).ECG(Original_ECG_5),'r*')
hold on
plot(data(5).t(Original_ECG_5+300),data(5).ECG(Original_ECG_5+300),'r*')

subplot(2,1,2)
hold on
plot(data(5).t(Original_ECG_5),data(5).PCG(Original_ECG_5),'r*')
hold on
plot(data(5).t(Original_ECG_5+300),data(5).PCG(Original_ECG_5+300),'r*')


%% PDG



segments1=zeros(22,50);
for i=1:length(peakOnsets_1_Original)
 segments1(i,:)=data(1).PCG(peakOnsets_1_Original(i)+300:peakOnsets_1_Original(i)+349);

end


segments2=zeros(16,50);
for i=1:length(Original_ECG_2)
 segments2(i,:)=data(2).PCG(Original_ECG_2(i)+300:Original_ECG_2(i)+349);

end

segments3=zeros(27,50);
for i=1:length(Original_ECG_3)
 segments3(i,:)=data(3).PCG(Original_ECG_3(i)+300:Original_ECG_3(i)+349);

end

segments4=zeros(39,50);
for i=1:length(Original_ECG_4)
 segments4(i,:)=data(4).PCG(Original_ECG_4(i)+300:Original_ECG_4(i)+349);

end

segments5=zeros(30,50);
for i=1:length(Original_ECG_5)
 segments5(i,:)=data(5).PCG(Original_ECG_5(i)+300:Original_ECG_5(i)+349);

end
%% PSD part 


[Pxx1,f1]=pwelch(segments1',[],[],[],Fs);
Average1=mean(Pxx1,2);
freq1 = meanfreq(f1,Fs) 

[Pxx2,f2]=pwelch(segments2',[],[],[],Fs);
Average2=mean(Pxx2,2);
freq2 = meanfreq(f2,Fs) 


[Pxx3,f3]=pwelch(segments3',[],[],[],Fs);
Average3=mean(Pxx3,2);
freq3 = meanfreq(f3,Fs) 


[Pxx4,f4]=pwelch(segments4',[],[],[],Fs);
Average4=mean(Pxx4,2);
freq4 = meanfreq(f4,Fs) 

[Pxx5,f5]=pwelch(segments5',[],[],[],Fs);
Average5=mean(Pxx5,2);
freq5 = meanfreq(f5,Fs) 





%% plot 
%first patiant 
fig11=figure
subplot(5,3,1)
plot(segments1(1,:))

subplot(5,3,2)
plot(Pxx1(:,1))

subplot(5,3,3)
plot(Average1)

%second patiant
figure(fig11)
subplot(5,3,4)
plot(segments2(1,:))

subplot(5,3,5)
plot(Pxx2(:,1))

subplot(5,3,6)
plot(Average2)


%third paitiant 
figure(fig11)
subplot(5,3,7)
plot(segments3(1,:))

subplot(5,3,8)
plot(Pxx3(:,1))

subplot(5,3,9)
plot(Average3)

%fourth patiant
figure(fig11)
subplot(5,3,10)
plot(segments4(1,:))

subplot(5,3,11)
plot(Pxx4(:,1))

subplot(5,3,12)
plot(Average4)

%fifth patiant
figure(fig11)
subplot(5,3,13)
plot(segments5(1,:))

subplot(5,3,14)
plot(Pxx5(:,1))

subplot(5,3,15)
plot(Average5)
