%Phase-shifting Digital Holographic Microscopy using PCA
%PCA-based method was proposed by Vargas et al. in 2011 [https://opg.optica.org/ol/abstract.cfm?uri=ol-36-8-1326]
%Code created by ADoblas, February 13, 2023

clc;
clear all;
close all;
% ____________________Parameters_____________________________
N=1024;               %Number of pixels
wavelength=0.532;   %Source's wavelenght in microns
k=2*pi/wavelength;        %Wavenumber in microns
DeltaX=4.65;     %Pitch in X direction in microns
DeltaY=4.65;     %Pitch in Y direction in microns

[X Y]=meshgrid(-N/2+1:N/2,-N/2+1:N/2); %Dicrete transversal coordinates

%% _____________Simulated Object______________

A=phantom(N);%phantom
A=imresize(A,[N N]);
A=A/max(A(:));
%figure;colormap gray;imagesc(A); axis image

Obj=exp(1i.*A); %create a phase phantom

% 
% figure,
% subplot(121);imshow(abs(Obj),[]); axis image;
% subplot(122);imshow(angle(Obj),[]); axis image;

FTObj=FT(Obj); %computation of the Fourier Transform
%figure;colormap gray;imagesc(log(abs(FTObj)));

%True phase distribution
Phiobj=angle(Obj);
Phiobj=(Phiobj-Phiobj(1,1))/max(max(Phiobj)); %normalization by subtractign the minimum value and dividing by the maximum one

%% ____________________ Reference Wave __________________________

%Generation of N-th number of reference waves with arbitrary phase steps

nophaseshift = 3; %number of total phase steps 

for k =1:1%100000   %repeat process to evaluate the robustness
    k
alpha = randi([0 360], 1,nophaseshift)*pi/180;
    %alpha = 2*pi*rand(1,2); %phase steps in rad between 0 and 2pi
phi = alpha;

for i = 1:length(phi)
    Cref=round(N/2)+1;Rref=round(N/2)+1;
    Cmax = Cref; Rmax = Rref; %completely in-line DHM system since there is no angle between the two waves
    ThetaXM=asin((Cref-Cmax)*wavelength/(N*DeltaX));
    ThetaYM=asin((Rref-Rmax)*wavelength/(N*DeltaX));
    R(:,:,i)=exp(1i*k*(sin(ThetaXM)*X*DeltaX+sin(ThetaYM)*Y*DeltaX)+1i*phi(i)); %different reference waves with different phi values
end


% Generation of the interference between the reference wave and the object waves
for i = 1:length(phi)   
    H(:,:,i)=abs(squeeze(R(:,:,i))+Obj).^2; %Hologram = In number of interferograms
end

% figure;colormap gray;imagesc(squeeze(H(:,:,1)));axis image
% figure;colormap gray;imagesc(squeeze(H(:,:,2)));axis image
% figure;colormap gray;imagesc(squeeze(H(:,:,3)));axis image


FTH1=FT(H(:,:,1)); %Fourier Transform of the first hologram to check that the system is in-line
% figure;colormap gray;imagesc(log(abs(FTH1)));axis image

%% _ Implement the PCA method
%We have modified the orignal pcaDemod.m code developed by Vargas et al.

[pw,Mod,U1,U2,V]=pcaDemod_ad(H);
% figure;colormap gray; imagesc(pw);axis image


%The estimated phase from pca (pw in code) is compared with the true
%values of the phase (Phiobj defined in line 36). Note that Phiobj is
%normalized so you need to normalize pw.
pw=(pw-pw(1,1))/max(max(pw)); %normalization

% figure;
% subplot(121);colormap gray; imagesc(Phiobj);axis image
% subplot(122);colormap gray; imagesc(pw);axis image

%For comparison purposes, estimate the MSE (mean square error) using the
%MATLAB built-in function immse(Phiobj,pw)

err(k) = immse(Phiobj,pw);
corr(k) = corr2(Phiobj,pw);


end

