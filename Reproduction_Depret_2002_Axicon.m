%Simulate the GB propagation through a converging lens and conical lens --
%reproduce the Depret 2002 Fig 3 & Fig 5
%Completed on 01/04/2021


% Clear the memory, close all figures and clear the screen
clear; close; clc;

%Input Coordinate System
x_range = 20000 * 10 ^ (-6);
y_range = 20000 * 10 ^ (-6);
N = 1000;
period = x_range/N;
lambda = 852 * 10 ^ (-9);

nPixel = (-N/2:(N/2-1)); % create an array of the sample number, with zero in the centre
x = nPixel/N*x_range; % calculate the spatial coordinates
y = nPixel/N*y_range; % calculate the spatial coordinates
[X,Y] = meshgrid(x,y);

%Fourier transformed coordinate system
dx = x_range/N;      %spatial period x-direction
dy = y_range/N;      %spatial period y-direction

fx=-1/(2*dx):1/x_range:1/(2*dx)-1/x_range;  % calculate the FT coordinates(Scaled)
fy=-1/(2*dy):1/y_range:1/(2*dy)-1/y_range;  % calculate the FT coordinates(Scaled)
FT_period = 1/x_range;
FT_full_range = N/x_range;

[FX,FY] = meshgrid(fx,fy);

%Spatial Fourier transformed coordinate system--Observ Plane
z_propafteroptics = 115 * 10 ^ (-3);
x_prime = fx.*lambda.*z_propafteroptics;    % calculate the FT spatial coordinates(Scaled)
y_prime = fy.*lambda.*z_propafteroptics;    % calculate the FT spatial coordinates(Scaled)
[X_prime,Y_prime] = meshgrid(x_prime,y_prime);


%Gaussian Input Beam Parameters
k = 2*pi/lambda;
A = 1;          %amplitude of GB
w0 = 645 * 10 ^ (-6);  %Beam waist following paper
zR  = pi * w0 ^ 2/lambda;   %Rayleigh range

z = 0;         %distance(in m) from the beam waist at the position of the lens = z_bef
t = 0;

rho = sqrt(X.^2 + Y.^2);      %radial distance squared note the .^ rather than ^


if z == 0
    Rz = Inf; %at beam waist plane wave
else   
    Rz = z + zR^2/z;%at z =! 0
end

[E_GB, I_GB, power] = getGB(A,z,lambda,w0,zR, t, rho, Rz);%Obtain GB profile and its total power


%{
%%%%%%Trial TEST%%%%%
%circular aperture
ySlit = zeros(N);
ySlit(sqrt(X.^2 + Y.^2) < (80 * 10^(-6))) = 1;

%rectangular aperture
ySlit2 = zeros(N);
ySlit2(abs(X) < (80 * 10^(-6)) & abs(Y) < (200 * 10^(-6))) = 1;
%}

%plot GB input
figure(1);
pcolor(X, Y, I_GB);shading interp;colorbar;axis equal;
%pcolor(X, Y, ySlit);shading interp;colorbar;axis equal;
title('Input GB Intensity Profile');
xlabel('x(m)')
ylabel('y(m)')
zlabel('Amplitude')

%%%%%%Propagation through converging lens + Axicon%%%%%%
r_axi = 5 * 10 ^ (-3);            %radius of axicon lens = 5mm
alpha = 2 * pi / 180;             %2 degrees
f = 100 * 10 ^ (-3);               %converging lens focal length
R_exp = 3.5 * 10 ^(-3);           %curvature of the axicon lens -- for practical axicon
e0 = 0 * 10 ^ (-3);               %plane axicon thickness
n = 1.51;                         %refractive index of the axicon lens and converging lens

wz = w0 * sqrt( 1 + z^2 / zR^2);    %beam radius

%%%%%Comparing with the geometric approximation -- Edmund optics%%%%%
ring_thick_geo = wz * sqrt(1 - n^2 * (sin(alpha))^2)/((cos(alpha)) * (n*(sin(alpha))^2 + (cos(alpha)) * sqrt(1 - n^2 * (sin(alpha))^2)));

ring_outer_radius_geo = z_propafteroptics * (sin(alpha) * (n* cos(alpha) - sqrt(1 - n^2*(sin(alpha))^2)))/(n*(sin(alpha))^2 + cos(alpha) * sqrt(1 - n^2 * (sin(alpha))^2));

ring_inner_radius_geo = ring_outer_radius_geo - ring_thick_geo;


%Trial test
%U4 = propTF2dim_Isabel(ySlit,X,Y,lambda,n,r_axi,alpha,z_propafteroptics,X_prime,Y_prime);
%I_U4 = abs(U4.^2);

%{
%Ideal case nooptics FFT
U4_Fraunhofer = propTF2dim_nooptics_Isabel(E_GB,X,Y,k,n,r_axi,alpha,f, R_exp,z_propafteroptics,X_prime,Y_prime);
I_Fraunhofer = abs(U4_Fraunhofer.^2);
%}


%{
%Ideal case Axicon doublet FFT
U4_Fraunhofer = propTF2dim_ideal_Isabel(E_GB,X,Y,k,n,r_axi,alpha,f, R_exp,z_propafteroptics,X_prime,Y_prime);
I_Fraunhofer = abs(U4_Fraunhofer.^2);
Mall_ideal = max(I_Fraunhofer,[],'all');
I_Fraunhofer = I_Fraunhofer/Mall_ideal;
%}



%Practical case Axicon doublet FFT
U4_Fraunhofer_prac = propTF2dim_practical_Isabel(E_GB,X,Y,k,n,r_axi,alpha,f, R_exp,e0,z_propafteroptics,X_prime,Y_prime);
I_Fraunhofer_prac = abs(U4_Fraunhofer_prac.^2);


%Manual Integration(prac)
surf_element = (1/x_range*lambda.*z_propafteroptics) * (1/y_range*lambda.*z_propafteroptics);
FFT_power_manual = sum(I_Fraunhofer_prac .* surf_element,'all');

%trapz integration(prac)
FFT_power=trapz(x_prime,trapz(y_prime,I_Fraunhofer_prac,2));

%%%%rescale the output beam intensity so that the total power of the output
%%%%beam = input GB power%%%%
I_Fraunhofer_prac_scl = I_Fraunhofer_prac .* power/FFT_power_manual;
FFT_power_manual_scl = sum(I_Fraunhofer_prac_scl .* surf_element,'all');  %should = GB power


%plot axicon diffraction output
figure(3);
%pcolor(X_prime,Y_prime, I_U4);shading interp;colorbar;axis equal;
%pcolor(X_prime * 1000,Y_prime * 1000, I_Fraunhofer);shading interp;colorbar;axis equal;
pcolor(X_prime * 1000,Y_prime * 1000, I_Fraunhofer_prac_scl);shading interp;colorbar;axis equal;
drawellipse('Center',[X_prime(501,501) * 1000, Y_prime(501,501) * 1000],'SemiAxes',[ring_outer_radius_geo * 1000, ring_outer_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");

drawellipse('Center',[X_prime(501,501) * 1000, Y_prime(501,501) * 1000],'SemiAxes',[ring_inner_radius_geo * 1000, ring_inner_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");
title('Fraunhofer Intensity after PRAC Axicon doublet at waist, zprop = 115mm, f = 100mm,Rexp = 3.5mm');
xlabel('xprime(mm)')
ylabel('yprime(mm)')
zlabel('Amplitude')

%plot axicon diffraction output 1D
figure(4);
plot(x_prime * 1000,I_Fraunhofer_prac_scl(:,501),'LineWidth',1); 
xline(ring_outer_radius_geo * 1000, '-r','LineWidth',1)
xline(ring_inner_radius_geo * 1000, '-g','LineWidth',1)
%plot(x_prime,I_U4(:,1000)); 
title('1D Fraunhofer Intensity after PRAC Axicon doublet at waist, zprop = 115mm, f = 100mm,Rexp = 3.5mm');
xlabel('xprime(mm)')
ylabel('Amplitude')



%{
%%%%%%Variation of Central Intensity(Ring centre)%%%%%%%%
zrange = 0:0.1:2*zR;
central_intensity = zeros(1,length(zrange));
for i = 1:length(zrange)
    if zrange(i) == 0
        Rz = Inf; %at beam waist plane wave
    else   
        Rz = zrange(i) + zR^2/zrange(i);%at z =! 0
    end
    [E_GB, I_GB] = getGB(A,zrange(i),lambda,w0,zR, t, rho, Rz);
    
    U4_Fraunhofer_prac = propTF2dim_practical_Isabel(E_GB,X,Y,k,n,r_axi,alpha,f, R_exp,e0,z_propafteroptics,X_prime,Y_prime);
    I_Fraunhofer_prac = abs(U4_Fraunhofer_prac.^2);
    Mall_ideal = max(I_Fraunhofer_prac,[],'all');
    I_Fraunhofer_prac = I_Fraunhofer_prac/Mall_ideal;
    central_intensity(i) = I_Fraunhofer_prac(501,501);
end

%}

%{
%%%%%%Plot z_beforelens vs Central Intensity%%%%%%
figure(5);

plot(zrange,central_intensity,'b','LineWidth',1); 

title('Ring centre intensity, zprop = 115mm, f = 100mm, zbeforelens = 0:0.1:2zR');
xlabel('zbreforelens(m)')
ylabel('Central normalized intensity')
%}

function[E_GB, I_GB, power] = getGB(A,z,lambda,w0,zR, t, rho, Rz)
%%%return GB field, intensity & total input power, power will be compared
%%%with the rescaled output FFT beam
    %A: GB constant
    %z: z_bef location of the beam at first optics relative to the beam
    %waist
    %w0: beam waist
    %zR: Rayleigh range
    %t: time
    %rho: polor coordinate
    %Rz: beam curvature
    c = 3 * 10 ^8;
    k = 2*pi/lambda;            %wavevector
    angularfreq = 2 * pi * c / lambda;
    wz = w0 * sqrt( 1 + z^2 / zR^2);    %beam radius
    phi = atan(z/zR);
    E_GB = A/sqrt(1 + z^2/zR^2) .* exp(1i * (k*z - angularfreq*t)) .* exp( - rho.^2./wz^2) .* exp(1i * (k * rho.^2./(2*Rz))) .* exp(-1i * phi);
    I_GB = abs(E_GB.^2);
    
    Integration_intensity = @(rho,theta) rho.*abs((A/sqrt(1 + z^2/zR^2) .* exp(1i * (k*z - angularfreq*t)) .* exp( - rho.^2./wz^2) .* exp(1i * (k * rho.^2./(2*Rz))) .* exp(-1i * phi)).^2);
    power = integral2(Integration_intensity,0,Inf,0,2*pi);
end

function[U4_Fraunhofer]=propTF2dim_nooptics_Isabel(u1,X,Y,k,n,r_axi,alpha,f, R_exp, z_propafteroptics,X_prime,Y_prime)
% this funtion calculates de propagation of the field u over a distance z
% to the plane observation plane, the window size is the same for source
% and obesrvation plane
% Fraunhofer Diffraction Integral

% propagation - transfer function approach
% using fresnel approximation
% assumes equal side lengths for x and y
% uniform sampling
% u1 - source plane field
% x - input beam coordinate
% y - input beam coordinate
% lambda - wavelength (monochromatic approximation)
% n - refractive index of optics
% r_axi - radius of axicon
% alpha - conical lens base angle
% z_propafteroptics - propagation distance
% u2 - observation plane field

H = 1;
U2=fftshift(H.*u1);                            %multiply U2: field after optics H
U3=fftshift(fft2(U2));             %inverse fft and center observ field

U4_Fraunhofer = -1i.*k/(2 * pi).*exp(1i .* k .* z_propafteroptics)./z_propafteroptics * exp((1i.*k./(2*z_propafteroptics)).*(X_prime.^2 + Y_prime.^2)).*U3;
end

function[U4_Fraunhofer]=propTF2dim_ideal_Isabel(u1,X,Y,k,n,r_axi,alpha,f, R_exp, z_propafteroptics,X_prime,Y_prime)
% this funtion calculates de propagation of the field u over a distance z
% to the plane observation plane, the window size is the same for source
% and obesrvation plane
% Fraunhofer Diffraction Integral

% propagation - transfer function approach
% using fresnel approximation
% assumes equal side lengths for x and y
% uniform sampling
% u1 - source plane field
% x - input beam coordinate
% y - input beam coordinate
% lambda - wavelength (monochromatic approximation)
% n - refractive index of optics
% r_axi - radius of axicon
% alpha - conical lens base angle
% z_propafteroptics - propagation distance
% u2 - observation plane field

H1 = exp( -1i* k* 0.5* (X.^2+Y.^2)./f);                     %transfer function lens
H2 = exp(-1i*k*(n-1)*(sqrt(X.^2+Y.^2)).* tan(alpha));      %transfer function conical lens on paper

H = H1.*H2;

U2=fftshift(H.*u1);                            %multiply U2: field after optics H
U3=fftshift(fft2(U2));             %inverse fft and center observ field

U4_Fraunhofer = -1i.*k/(2 * pi).*exp(1i .* k .* z_propafteroptics)./z_propafteroptics * exp((1i.*k./(2*z_propafteroptics)).*(X_prime.^2 + Y_prime.^2)).*U3;
end


function[U4_Fraunhofer]=propTF2dim_practical_Isabel(u1,X,Y,k,n,r_axi,alpha,f, R_exp,e0, z_propafteroptics,X_prime,Y_prime)
% this funtion calculates de propagation of the field u over a distance z
% to the plane observation plane, the window size is the same for source
% and obesrvation plane
% Fraunhofer Diffraction Integral

% propagation - transfer function approach
% using fresnel approximation
% assumes equal side lengths for x and y
% uniform sampling
% u1 - source plane field
% X - input beam coordinate
% Y - input beam coordinate
% k - wavevector (monochromatic approximation)
% n - refractive index of optics
% r_axi - radius of axicon
% alpha - conical lens base angle
% f - focal length
% R_exp - curvature of the pratical axicon
% e0 - plane axicon thickness
% z_propafteroptics - propagation distance
% X_prime - observation plane coordinate
% Y_prime - observation plane coordinate

factor4 = exp(0.5.*1i .* k ./z_propafteroptics.*(X.^2 + Y.^2));

H1 = exp( -1i* k* 0.5* (X.^2+Y.^2)./f);                     %transfer function lens
H2 = exp(-1i*k*(n-1)*(sqrt(X.^2+Y.^2)).* tan(alpha));      %worktransfer function conical lens on paper
thickness_function = e0 - R_exp .* (tan(alpha))^2 .* sqrt(1 + (X.^2 + Y.^2)./(R_exp .* tan(alpha))^2);
H3 = exp(1i*k*(n-1).*thickness_function);  %work

H = H3;
u1 = H.*u1;
U2=fftshift(factor4.*u1);                            %multiply U2: field after optics H
U3=fftshift(fft2(U2));             %inverse fft and center observ field

U4_Fraunhofer = -1i.*k/(2 * pi).*exp(1i .* k .* z_propafteroptics)./z_propafteroptics * exp((1i.*k./(2*z_propafteroptics)).*(X_prime.^2 + Y_prime.^2)).*U3;
end


