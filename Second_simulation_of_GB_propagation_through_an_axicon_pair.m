% Clear the memory, close all figures and clear the screen
clear; close; clc;

%Input Coordinate System
x_range = 10000 * 10 ^ (-6);
y_range = 10000 * 10 ^ (-6);
N = 4000;
period = x_range/N;
lambda = 828 * 10 ^ (-9);

nPixel = (-N/2:(N/2-1)); % create an array of the sample number, with zero in the centre
x = nPixel/N*x_range; % calculate the spatial coordinates
y = nPixel/N*y_range; % calculate the spatial coordinates
[X,Y] = meshgrid(x,y);


%Second Axicon Coordinate System
dx5 = period;      %spatial period x-direction
dy5 = period;       %spatial period y-direction

fx5=-1/(2*dx5):1/x_range:1/(2*dx5)-1/x_range;  % calculate the FT coordinates(Scaled)
fy5=-1/(2*dy5):1/y_range:1/(2*dy5)-1/y_range;  % calculate the FT coordinates(Scaled)
FT5_period = 1/y_range;
FT5_full_range = N/y_range;


[FX5,FY5] = meshgrid(fx5,fy5);

z4 = 120 * 10 ^(-3);           % separation between two axicons
x_prime5 = fx5.*lambda.*z4;    % calculate the FT spatial coordinates(Scaled)
y_prime5 = fy5.*lambda.*z4;    % calculate the FT spatial coordinates(Scaled)
[X_prime5,Y_prime5] = meshgrid(x_prime5,y_prime5);
xprime5_period = FT5_period.*lambda.*z4; 
xprime5_full_range = FT5_full_range.*lambda.*z4; 


%OBSERV plane Coordinate System
dx6 = xprime5_period;      %spatial period x-direction
dy6 = xprime5_period;       %spatial period y-direction

fx6=-1/(2*dx6):1/xprime5_full_range:1/(2*dx6)-1/xprime5_full_range;  % calculate the FT coordinates(Scaled)
fy6=-1/(2*dy6):1/xprime5_full_range:1/(2*dy6)-1/xprime5_full_range;  % calculate the FT coordinates(Scaled)
FT6_period = 1/xprime5_full_range;
FT6_full_range = N/xprime5_full_range;


[FX6,FY6] = meshgrid(fx6,fy6);

z5 = 0.2;                      % position of the OBSERV plane
x_prime6 = fx6.*lambda.*z5;    % calculate the FT spatial coordinates(Scaled)
y_prime6 = fy6.*lambda.*z5;    % calculate the FT spatial coordinates(Scaled)
[X_prime6,Y_prime6] = meshgrid(x_prime6,y_prime6);
xprime6_period = FT6_period.*lambda.*z5; 
xprime6_full_range = FT6_full_range.*lambda.*z5; 

%Gaussian Input Beam
k = 2*pi/lambda;
A = 1;
w0 = 2700/2 * 10 ^ (-6);
zR  = pi * w0 ^ 2/lambda;   %Rayleigh range

z = 0;         %distance(in m) from the beam waist at the position of the lens
t = 0;

rho = sqrt(X.^2 + Y.^2);      %radial distance squared note the .^ rather than ^

if z == 0
    Rz = Inf; %at beam waist plane wave
else   
    Rz = z + zR^2/z;%at z =! 0
end

[E_GB, I_GB, power] = getGB(A,z,lambda,w0,zR, t, rho, Rz);%Obtain GB profile and its total power



%plot GB input
figure(1);
pcolor(X, Y, I_GB);shading interp;colorbar;axis equal;
title('Input GB Intensity Profile, w0 = 2.7mm/2, zbef = 0m');
xlabel('x(m)')
ylabel('y(m)')
zlabel('Amplitude')

%%%%%%Propagation through Axicon%%%%%%
r_axi = 25.4 /2* 10 ^ (-3);       %radius of axicon lens = 25.4mm
alpha = 5 * pi / 180;             %5 degrees
R_exp = 1.0 * 10 ^(-3);           %curvature of the axicon lens
e0 = 6.1 * 10 ^ (-3);             %plane axicon thickness
n = 1.45;                         %refractive index of the axicon lens and converging lens

wz = w0 * sqrt( 1 + z^2 / zR^2);    %beam radius

%%%%%Comparing with the geometric approximation -- Edmund optics%%%%%
ring_thick_geo = wz * sqrt(1 - n^2 * (sin(alpha))^2)/((cos(alpha)) * (n*(sin(alpha))^2 + (cos(alpha)) * sqrt(1 - n^2 * (sin(alpha))^2)));

ring_outer_radius_geo = z4 * (sin(alpha) * (n* cos(alpha) - sqrt(1 - n^2*(sin(alpha))^2)))/(n*(sin(alpha))^2 + cos(alpha) * sqrt(1 - n^2 * (sin(alpha))^2));

ring_inner_radius_geo = ring_outer_radius_geo - ring_thick_geo;


%Practical case Axicon doublet FFT
[U15,U19] = propTF2dim_DoublePracAxicon(E_GB,X,Y,X_prime5,Y_prime5,X_prime6,Y_prime6,k,n,r_axi,alpha,R_exp,e0,z4,z5);
I15= abs(U15.^2);   %Intensity profile at the entrance of the second axicon
I19= abs(U19.^2);   %Intensity profile at the OBSERV plane


%Manual Integration(prac)
surf_element = xprime6_period * xprime6_period;
FFT_power_manual = sum(I19 .* surf_element,'all');  

I19 = I19 .* power/FFT_power_manual; %Rescaled intensity profile at OBSERV plane so the total power = GB power
FFT_power_manual_scl = sum(I19 .* surf_element,'all');
Power_check = sum(I19 .* surf_element,'all'); %should = GB power

%Fresnel Intensity Before Second PRAC Axicon
figure(5);
pcolor(X_prime5 * 1000,Y_prime5 * 1000, I15);shading interp;colorbar;axis equal;
title('Fresnel Intensity Before Second PRAC Axicon, zbef = 0,z4 = 120mm,z5 = 0.2m,Rexp = 0.6mm');
xlabel('xprime5(mm)')
ylabel('yprime5(mm)')
zlabel('Amplitude')

%Fresnel Intensity at OBSERV plane
figure(6);
pcolor(X_prime6 * 1000,Y_prime6 * 1000, I19);shading interp;colorbar;axis equal;
drawellipse('Center',[X_prime6(2001,2001) * 1000, Y_prime6(2001,2001) * 1000],'SemiAxes',[ring_outer_radius_geo * 1000, ring_outer_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");
drawellipse('Center',[X_prime6(2001,2001) * 1000, Y_prime6(2001,2001) * 1000],'SemiAxes',[ring_inner_radius_geo * 1000, ring_inner_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");
title('Fresnel Intensity at OBSERV plane, zbef = 0, z4 = 120mm,z5 = 0.2m,Rexp = 1.0mm');
xlabel('xprime6(mm)')
ylabel('yprime6(mm)')
zlabel('Amplitude')
caxis([0 0.25])


%1D Fresnel Intensity bef Second Axicon
figure(10);
plot(x_prime6 * 1000,I15(:,2001),'LineWidth',1); 
xline(ring_outer_radius_geo * 1000, '-r','LineWidth',1)
xline(ring_inner_radius_geo * 1000, '-g','LineWidth',1)
title('1D Fresnel Intensity bef Second Axicon, IDEAL Axicon, zbef = 0, z4 = 120mm,z5 = 0.2m');
xlabel('xprime6(mm)')
ylabel('Amplitude')


%1D Fresnel Intensity at OBSERV plane
figure(7);
plot(x_prime6 * 1000,I19(:,2001),'LineWidth',1); 
xline(ring_outer_radius_geo * 1000, '-r','LineWidth',1)
xline(ring_inner_radius_geo * 1000, '-g','LineWidth',1)
title('1D Fresnel Intensity at OBSERV plane,IDEAL Axicon, zbef = 0, z4 = 120mm,z5 = 0.2m');
xlabel('xprime6(mm)')
ylabel('Amplitude')


%%%%%%plot experimental data with simulation data%%%%%%%
load('Iexpdata','xaxis2','Icrop1vert');
load('Iexpdata2','Icrop950');
load('Iideal.mat','Iideal');
figure(8);
xaxis = linspace(-4.60 * 10 ^ (-3),4.50 * 10 ^ (-3),2010);  %Reshape of exp data because the exp data ring size cannot be accurately measured
plot(xaxis*1000, 0.25/0.4 * Icrop1vert,'m','LineWidth',1)    %Intensity of the exp data reduced to match the simulated data
hold on
plot(x_prime6 * 1000,I19(:,2001),'b','LineWidth',1); 
hold on
plot(x_prime6 * 1000,0.25/0.12 * Iideal,'k--','LineWidth',1); 
xline(ring_outer_radius_geo * 1000, '-r','LineWidth',1)
xline(ring_inner_radius_geo * 1000, '-g','LineWidth',1)
legend('Exp Seed Input Vertical','Simulation 2nd axicon decenter = 0.03axicon,CORRECTED axicon, zbef = 0, z4 = 120mm,z5 = 0.2m, Rexp = 1.0mm','Simulation no decenter,IDEAL axicon, zbef = 0, z4 = 120mm,z5 = 0.2m,')
title('Comparison between Exp and Sim with central maxima');
xlabel('xprime6(mm)')
ylabel('Amplitude')
xlim([-9.00,9.00])
ylim([0, 0.28])

%%%%%Ring power Content%%%%%
I19_crop = I19;
I19_crop(sqrt(X_prime6.^2 + Y_prime6.^2) < 0.3*ring_outer_radius_geo) = 0;


%Cropped Intensity at OBSERV plane
figure(9);
pcolor(X_prime6 * 1000,Y_prime6 * 1000, I19_crop);shading interp;colorbar;axis equal;
drawellipse('Center',[X_prime6(2001,2001) * 1000, Y_prime6(2001,2001) * 1000],'SemiAxes',[ring_outer_radius_geo * 1000, ring_outer_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");

drawellipse('Center',[X_prime6(2001,2001) * 1000, Y_prime6(2001,2001) * 1000],'SemiAxes',[0.3*ring_outer_radius_geo * 1000, 0.3*ring_outer_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");
title('Fresnel Intensity at OBSERV plane, ring power = 99.9%,zbef = 0, z4 = 120mm,z5 = 0.2m,Rexp = 0.1mm');
xlabel('xprime6(mm)')
ylabel('yprime6(mm)')
zlabel('Amplitude')

%Obtain ring power%
surf_element = xprime6_period * xprime6_period;
FFT_power_manual_crop = sum(I19_crop .* surf_element,'all');
ring_power = FFT_power_manual_crop/Power_check;

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
    E_GB = A./sqrt(1 + z.^2/zR.^2) .* exp(1i .* (k.*z - angularfreq.*t)) .* exp( - rho.^2./wz.^2) .* exp(1i .* (k .* rho.^2./(2.*Rz))) .* exp(-1i .* phi);
    I_GB = abs(E_GB.^2);
    
    Integration_intensity = @(rho,theta) rho.*abs((A/sqrt(1 + z^2/zR^2) .* exp(1i * (k*z - angularfreq*t)) .* exp( - rho.^2./wz^2) .* exp(1i * (k * rho.^2./(2*Rz))) .* exp(-1i * phi)).^2);
    power = integral2(Integration_intensity,0,Inf,0,2*pi);
end




function[U15,U19]=propTF2dim_DoublePracAxicon(u1,X_prime4,Y_prime4,X_prime5,Y_prime5,X_prime6,Y_prime6,k,n,r_axi,alpha,R_exp,e0,z4,z5)
% this funtion calculates de propagation of the field u over a distance z
% to the plane observation plane, the window size is the same for source
% and obesrvation plane
% Fraunhofer Diffraction Integral

% propagation - transfer function approach
% using fresnel approximation
% assumes equal side lengths for x and y
% uniform sampling
% u1 - source plane field
% X_prime4 - input beam coordinate
% y_prime4 - input beam coordinate
% X_prime5 - beam coordinate at second axicon
% y_prime5 - beam coordinate at second axicon
% X_prime6 - beam coordinate at OBSERV plane
% y_prime6 - beam coordinate at OBSERV plane
% k - wavevector (monochromatic approximation)
% n - refractive index of optics
% r_axi - radius of axicon
% alpha - conical lens base angle
% f - focal length
% R_exp - curvature of the pratical axicon
% e0 - plane axicon thickness
% z4 - separation between axicons
% z5 - position of OBSERV plane
 

%%%Transmission_through_axicon_first
%decenter = 0.03 * r_axi;
H2 = exp(-1i*k*(n-1)*(sqrt(X_prime4.^2 + Y_prime4.^2)).* tan(alpha));      %worktransfer function conical lens on paper
%thickness_function1 = e0 - R_exp .* (tan(alpha)).^2 .* sqrt(1 + ((X_prime4+decenter).^2 + Y_prime4.^2)./(R_exp .* tan(alpha)).^2);
thickness_function1 = e0 - R_exp .* (tan(alpha)).^2 .* sqrt(1 + ((X_prime4).^2 + Y_prime4.^2)./(R_exp .* tan(alpha)).^2);
H3 = exp(1i.*k.*(n-1).*thickness_function1);                                     %transfer function of a single axicon
U12 = u1.*H3;

%%%Prop_to_second_axicon

factor4 = exp(0.5.*1i .* k ./z4.*(X_prime4.^2 + Y_prime4.^2));
U13 = fftshift(U12.*factor4);                                                         %field after first axicon
U14 = fftshift(fft2(U13));                                                      %FFT of U2
U15 =  -1i.*k./(2 .* pi).*exp(1i .* k .* z4)./z4 * exp((1i.*k./(2.*z4)).*(X_prime5.^2 + Y_prime5.^2)).*U14;


%%%Transmission_through_axicon_second
decenter = 0.03 * r_axi;
H4 = exp(-1i*k*(n-1)*(sqrt((X_prime5).^2 + Y_prime5.^2)).* tan(alpha));      %worktransfer function conical lens on paper
%thickness_function2 = e0 - R_exp .* (tan(alpha)).^2 .* sqrt(1 + ((X_prime5 + decenter).^2 + Y_prime5.^2)./(R_exp .* tan(alpha)).^2);
thickness_function2 = e0 - R_exp .* (tan(alpha)).^2 .* sqrt(1 + ((X_prime5).^2 + Y_prime5.^2)./(R_exp .* tan(alpha)).^2);
H5 = exp(1i.*k.*(n-1).*thickness_function2);                                     %transfer function of a single axicon
U16 = U15.*H5; 

%%%Prop_to_observ
factor5 = exp(0.5.*1i .* k ./z5.*(X_prime5.^2 + Y_prime5.^2));
U17=fftshift(U16.*factor5);                            %multiply U2: field after optics H
U18=fftshift(fft2(U17));             %inverse fft and center observ field

U19 = -1i.*k./(2 .* pi).*exp(1i .* k .* z5)./z5 .* exp((1i.*k./(2.*z5)).*(X_prime6.^2 + Y_prime6.^2)).*U18;
end
