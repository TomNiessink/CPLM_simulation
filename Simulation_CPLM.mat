%% CPLM simulator V2
% This script can calculate the effect of various optical elements,
% including crystals, on polarized light. All simulations are based
% on a mathematical model called 'Jones Calculus'. 
% Written by Tom Niessink, Univeristy of Twente
% E-mail: t.niessink@utwente.nl
% 3 April 2024
clear all; close all; clc;
%% Definition of optical elements
% Sources:  [1] Introduction to optics by Pedrotti^3 
%           [2] https://www.researchgate.net/publication/235963739_Obtainment_of_the_polarizing_and_retardation_parameters_of_a_non-depolarizing_optical_system_from_the_polar_decomposition_of_its_Mueller_matrix

% These are all Jones matrices.
H = [1 0; 0 0];         % Horizontal Polarizer [1]
V = [0 0; 0 1];         % Vertical Polarizer [1]
RHCP = .5*[1 1i;-1i 1];     % Right Hand Circular Polarizer [1], not used
LHCP = .5*[1 -1i;1i 1];     % Left Hand Circular Polarizer [1], not used
QWPH = exp(1i*pi/4)*[1 0;0 -1i];    % Quarter wave plate horizontal [1], not used
QWPV = exp(-1i*pi/4)*[1 0;0 1i];     % Quarter wave plate vertical [1], not used
LP = @(theta) [cos(theta)^2 cos(theta)*sin(theta);cos(theta)*sin(theta) sin(theta)^2]; % Linear polarizer (angular)
BFO = @(theta,eta) exp(-1*1i*eta/2)*[cos(theta)^2+exp(1i*eta)*sin(theta)^2 ...
    (1-exp(1i*eta))*cos(theta)*sin(theta) ; (1-exp(1i*eta))*cos(theta)*sin(theta) ...
    sin(theta)^2+exp(1i*eta)*cos(theta)^2];     % Arbirtrary birefringencent material as phase retarder (Linear) [2]


%% load Led_spectrum.csv
% LED data retrieved from EMPIR 15SIB07 PhotoLED, LED s57 was used.
load("LED_spectrum.csv")
lambda = (LED_spectrum(:,1))/1E9; % Wavelengths array
intensities = (LED_spectrum(:,2)); % Intensities per wavelength
intensities = intensities/max(intensities); % These are normalized for easy access

% There is a two-fold nested loop to calculate the image. By downsampling
% the LED spectrum we increase the the loop time with the downsampling
% factor
lambda = downsample(lambda,4);
intensities = downsample(intensities,4);


%% Define crystal
BF = 0.015; % Birefringence (delta n); -0.30 for MSU and 0.0155 for CPP
r = 1E-6; % Radiance of crystal (m)
l = 5E-6; % Length of crystal (m)
angle = pi/4; % Angle of crystal in radians with respect to the polarizer

Gridsize_px=100; % Gridsize in pixels
Gridsize_m = 10E-6; % Size of pixels in grid

[d,ang] = rod(r,l,angle,Gridsize_m,Gridsize_px); % Calls a function 'rod' which is in a separate file

%% Calculate image
I_out_Hg = zeros(size(lambda,1),100,100); % Prelocate data to save time


Sz_Hg = size(d); % The size of the image

for i=1:size(lambda,1) % Loop over all wavelengths in the wavelength table
I_in(:) = intensities(i); % Retrieve intensity of certain wavelength
H_dphi = (2*pi/lambda(i))*(BF)*d; % Calculate the retardance for the full image
dphiqwl = (2*pi/lambda(i))*550E-9; % Calculate the retardance for the additional quarter lambda plate

for j = 1:Sz_Hg(1) % Loop over one axis
    for k = 1:Sz_Hg(2) % Loop over other axis
        dphi=H_dphi(j,k); % Take the retardance of the current pixel
        %*BFO(pi/4,dphiqwl)
        % Calculate output
        I_out_Hg(i,j,k) = norm(I_in*H*BFO(ang(j,k),dphi)*LP(pi/2));
    end
end

end


%% color analysis
image = zeros(100,100,3); % Prelocate data space

for j = 1:Sz_Hg(1)
    for k = 1:Sz_Hg(2)
    Spectrum = I_out_Hg(:,j,k);
    image(j,k,1) = sum(Spectrum(find(lambda == 612E-9):find(lambda == 800E-9))); % Red bin
    image(j,k,2) = sum(Spectrum(find(lambda == 512E-9):find(lambda == 600E-9))); % Green bin
    image(j,k,3) = sum(Spectrum(find(lambda == 380E-9):find(lambda == 500E-9))); % Blue bin
    end
end



% The gain set with R, G, B is determined such that there was a 
% 'warm' white brightfield image
image(:,:,1)=image(:,:,1);
image(:,:,2)=image(:,:,2)*1.1;
image(:,:,3)=image(:,:,3)*4.2;
image = image/15;

figure()
imagesc(image)
