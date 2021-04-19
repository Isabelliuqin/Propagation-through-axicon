%Simulate the GB propagation through an axicon pair(2 conical lens) --
%Modification of Matlab File: Reproduction_Depret_2002_Axicon.m to include
%another conical lens and remove the converging lens
%Completed on 01/04/2021


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

%Prime1 Coordinate System
dx1 = period;      %spatial period x-direction
dy1 = period;      %spatial period y-direction

fx1=-1/(2*dx1):1/x_range:1/(2*dx1)-1/x_range;  % calculate the FT coordinates(Scaled)
fy1=-1/(2*dy1):1/y_range:1/(2*dy1)-1/y_range;  % calculate the FT coordinates(Scaled)
FT_period = 1/x_range;
FT_full_range = N/x_range;

[FX1,FY1] = meshgrid(fx1,fy1);


z1 = 15 * 10 ^ (-3);
x_prime1 = fx1.*lambda.*z1;    % calculate the FT spatial coordinates(Scaled)
y_prime1 = fy1.*lambda.*z1;    % calculate the FT spatial coordinates(Scaled)
[X_prime1,Y_prime1] = meshgrid(x_prime1,y_prime1);
xprime1_period = FT_period.*lambda.*z1; 
xprime1_full_range = FT_full_range.*lambda.*z1; 

%Prime2 Coordinate System
dx2 = xprime1_period;      %spatial period x-direction
dy2 = xprime1_period;      %spatial period y-direction

fx2=-1/(2*dx2):1/xprime1_full_range:1/(2*dx2)-1/xprime1_full_range;  % calculate the FT coordinates(Scaled)
fy2=-1/(2*dy2):1/xprime1_full_range:1/(2*dy2)-1/xprime1_full_range;  % calculate the FT coordinates(Scaled)
FT2_period = 1/xprime1_full_range;
FT2_full_range = N/xprime1_full_range;


[FX2,FY2] = meshgrid(fx2,fy2);

z2 = 0.6;
x_prime2 = fx2.*lambda.*z2;    % calculate the FT spatial coordinates(Scaled)
y_prime2 = fy2.*lambda.*z2;    % calculate the FT spatial coordinates(Scaled)
[X_prime2,Y_prime2] = meshgrid(x_prime2,y_prime2);
xprime2_period = FT2_period.*lambda.*z2; 
xprime2_full_range = FT2_full_range.*lambda.*z2; 

%Prime3 Coordinate System
x_prime3 = x_prime2*2;    % calculate the FT spatial coordinates(Scaled)
y_prime3 = y_prime2*2;    % calculate the FT spatial coordinates(Scaled)
[X_prime3,Y_prime3] = meshgrid(x_prime3,y_prime3);
xprime3_period = xprime2_period*2; 
xprime3_full_range = xprime2_full_range*2; 

%Prime4 Coordinate System
dx4 = xprime3_period;      %spatial period x-direction
dy4 = xprime3_period;       %spatial period y-direction

fx4=-1/(2*dx4):1/xprime3_full_range:1/(2*dx4)-1/xprime3_full_range;  % calculate the FT coordinates(Scaled)
fy4=-1/(2*dy4):1/xprime3_full_range:1/(2*dy4)-1/xprime3_full_range;  % calculate the FT coordinates(Scaled)
FT4_period = 1/xprime3_full_range;
FT4_full_range = N/xprime3_full_range;


[FX4,FY4] = meshgrid(fx4,fy4);

z3 = 0.2;
x_prime4 = fx4.*lambda.*z3;    % calculate the FT spatial coordinates(Scaled)
y_prime4 = fy4.*lambda.*z3;    % calculate the FT spatial coordinates(Scaled)
[X_prime4,Y_prime4] = meshgrid(x_prime4,y_prime4);
xprime4_period = FT4_period.*lambda.*z3; 
xprime4_full_range = FT4_full_range.*lambda.*z3; 

%Prime5 Coordinate System
dx5 = xprime4_period;      %spatial period x-direction
dy5 = xprime4_period;       %spatial period y-direction

fx5=-1/(2*dx5):1/xprime4_full_range:1/(2*dx5)-1/xprime4_full_range;  % calculate the FT coordinates(Scaled)
fy5=-1/(2*dy5):1/xprime4_full_range:1/(2*dy5)-1/xprime4_full_range;  % calculate the FT coordinates(Scaled)
FT5_period = 1/xprime4_full_range;
FT5_full_range = N/xprime4_full_range;


[FX5,FY5] = meshgrid(fx5,fy5);

z4 = 120 * 10 ^(-3);
x_prime5 = fx5.*lambda.*z4;    % calculate the FT spatial coordinates(Scaled)
y_prime5 = fy5.*lambda.*z4;    % calculate the FT spatial coordinates(Scaled)
[X_prime5,Y_prime5] = meshgrid(x_prime5,y_prime5);
xprime5_period = FT5_period.*lambda.*z4; 
xprime5_full_range = FT5_full_range.*lambda.*z4; 


%Prime6 Coordinate System
dx6 = xprime5_period;      %spatial period x-direction
dy6 = xprime5_period;       %spatial period y-direction

fx6=-1/(2*dx6):1/xprime5_full_range:1/(2*dx6)-1/xprime5_full_range;  % calculate the FT coordinates(Scaled)
fy6=-1/(2*dy6):1/xprime5_full_range:1/(2*dy6)-1/xprime5_full_range;  % calculate the FT coordinates(Scaled)
FT6_period = 1/xprime5_full_range;
FT6_full_range = N/xprime5_full_range;


[FX6,FY6] = meshgrid(fx6,fy6);

z5 = 0.2;
x_prime6 = fx6.*lambda.*z5;    % calculate the FT spatial coordinates(Scaled)
y_prime6 = fy6.*lambda.*z5;    % calculate the FT spatial coordinates(Scaled)
[X_prime6,Y_prime6] = meshgrid(x_prime6,y_prime6);
xprime6_period = FT6_period.*lambda.*z5; 
xprime6_full_range = FT6_full_range.*lambda.*z5; 

%Gaussian Input Beam
k = 2*pi/lambda;
A = 1;
w0 = 5 * 10 ^ (-6);
zR  = pi * w0 ^ 2/lambda;   %Rayleigh range

z = 0;         %distance(in m) from the beam waist at the position of the lens
t = 0;

rho = sqrt(X.^2 + Y.^2);      %radial distance squared note the .^ rather than ^


if z == 0
    Rz = Inf; %at beam waist plane wave
else   
    Rz = z + zR^2/z;%at z =! 0
end

[E_GB, I_GB, power] = getGB(A,z,lambda,w0,zR, t, rho, Rz);

%Integration_intensity = @(rho,theta) rho.*abs((A/sqrt(1 + z^2/zR^2) .* exp(1i * (k*z - angularfreq*t)) .* exp( - rho.^2./wz^2) .* exp(1i * (k * rho.^2./(2*Rz))) .* exp(-1i * phi)).^2);
%q = integral2(Integration_intensity,0,Inf,0,2*pi);


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
%surf(X,Y,I1);
title('Input GB Intensity Profile');
xlabel('x(m)')
ylabel('y(m)')
zlabel('Amplitude')

%%%%%%Propagation through converging lens + Axicon%%%%%%
r_axi = 25.4 /2* 10 ^ (-3);            %radius of axicon lens = 5mm
alpha = 2 * pi / 180;             %2 degrees
f = 15 * 10 ^ (-3);               %converging lens focal length
R_exp = 0.1 * 10 ^(-3);           %curvature of the axicon lens
e0 = 5 * 10 ^ (-3);               %plane axicon thickness
n = 1.51;                         %refractive index of the axicon lens and converging lens

wz = w0 * sqrt( 1 + z^2 / zR^2);    %beam radius

%{
ring_thick_geo = wz * sqrt(1 - n^2 * (sin(alpha))^2)/((cos(alpha)) * (n*(sin(alpha))^2 + (cos(alpha)) * sqrt(1 - n^2 * (sin(alpha))^2)));

ring_outer_radius_geo = z_propafteroptics * (sin(alpha) * (n* cos(alpha) - sqrt(1 - n^2*(sin(alpha))^2)))/(n*(sin(alpha))^2 + cos(alpha) * sqrt(1 - n^2 * (sin(alpha))^2));

ring_inner_radius_geo = ring_outer_radius_geo - ring_thick_geo;
%}

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
[U4,U5, U8,U11,U15,U19] = propTF2dim_DoublePracAxicon(E_GB,X,Y,X_prime1,Y_prime1,X_prime2,Y_prime2,X_prime3,Y_prime3,X_prime4,Y_prime4,X_prime5,Y_prime5,X_prime6,Y_prime6,k,n,f,r_axi,alpha,R_exp,e0,z1,z2,z3,z4,z5);
I4= abs(U4.^2);
I5 = abs(U5.^2);
I8= abs(U8.^2);
I11= abs(U11.^2);
I15= abs(U15.^2);
I19= abs(U19.^2);
%Mall_ideal = max(I_Fraunhofer_prac2,[],'all');
%I_Fraunhofer_prac = I_Fraunhofer_prac/Mall_ideal;


%{
%Manual Integration(prac)
surf_element = (1/x_range*lambda.*z_propafteroptics) * (1/y_range*lambda.*z_propafteroptics);
FFT_power_manual = sum(I_Fraunhofer_prac .* surf_element,'all');

%trapz integration(prac)
FFT_power=trapz(x_prime,trapz(y_prime,I_Fraunhofer_prac,2));


I_Fraunhofer_prac_scl = I_Fraunhofer_prac .* power/FFT_power_manual;
FFT_power_manual_scl = sum(I_Fraunhofer_prac_scl .* surf_element,'all');

%}
%{
%%%%%%%trial tests%%%%%%
FFT_COS = fftshift(fft2(fftshift(U4_Fraunhofer)));

FFT_slit = fftshift(fft2(fftshift(ySlit)));
%}

figure(2);
%pcolor(X_prime,Y_prime, I_U4);shading interp;colorbar;axis equal;
%pcolor(X_prime * 1000,Y_prime * 1000, I_Fraunhofer);shading interp;colorbar;axis equal;
pcolor(X_prime1 .* 1000,Y_prime1 .* 1000, I4);shading interp;colorbar;axis equal;
%circles(x,y,ring_outer_radius_geo * 1000,'color','black')
%drawellipse('Center',[X_prime2(501,501) * 1000, Y_prime2(501,501) * 1000],'SemiAxes',[ring_outer_radius_geo * 1000, ring_outer_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");

%drawellipse('Center',[X_prime2(501,501) * 1000, Y_prime2(501,501) * 1000],'SemiAxes',[ring_inner_radius_geo * 1000, ring_inner_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");
title('Fraunhofer Intensity after PRAC Axicon doublet at waist, zprop = 115mm, f = 100mm,Rexp = 3.5mm');
xlabel('xprime(mm)')
ylabel('yprime(mm)')
zlabel('Amplitude')

figure(8);
%pcolor(X_prime,Y_prime, I_U4);shading interp;colorbar;axis equal;
%pcolor(X_prime * 1000,Y_prime * 1000, I_Fraunhofer);shading interp;colorbar;axis equal;
pcolor(X_prime1 .* 1000,Y_prime1 .* 1000, I5);shading interp;colorbar;axis equal;
%circles(x,y,ring_outer_radius_geo * 1000,'color','black')
%drawellipse('Center',[X_prime2(501,501) * 1000, Y_prime2(501,501) * 1000],'SemiAxes',[ring_outer_radius_geo * 1000, ring_outer_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");

%drawellipse('Center',[X_prime2(501,501) * 1000, Y_prime2(501,501) * 1000],'SemiAxes',[ring_inner_radius_geo * 1000, ring_inner_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");
title('Fraunhofer Intensity after PRAC Axicon doublet at waist, zprop = 115mm, f = 100mm,Rexp = 3.5mm');
xlabel('xprime(mm)')
ylabel('yprime(mm)')
zlabel('Amplitude')



%plot axicon diffraction output
figure(3);
%pcolor(X_prime,Y_prime, I_U4);shading interp;colorbar;axis equal;
%pcolor(X_prime * 1000,Y_prime * 1000, I_Fraunhofer);shading interp;colorbar;axis equal;
pcolor(X_prime2 * 1000,Y_prime2 * 1000, I8);shading interp;colorbar;axis equal;
%circles(x,y,ring_outer_radius_geo * 1000,'color','black')
%drawellipse('Center',[X_prime2(501,501) * 1000, Y_prime2(501,501) * 1000],'SemiAxes',[ring_outer_radius_geo * 1000, ring_outer_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");

%drawellipse('Center',[X_prime2(501,501) * 1000, Y_prime2(501,501) * 1000],'SemiAxes',[ring_inner_radius_geo * 1000, ring_inner_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");
title('Fraunhofer Intensity after PRAC Axicon doublet at waist, zprop = 115mm, f = 100mm,Rexp = 3.5mm');
xlabel('xprime(mm)')
ylabel('yprime(mm)')
zlabel('Amplitude')


figure(9);
%pcolor(X_prime,Y_prime, I_U4);shading interp;colorbar;axis equal;
%pcolor(X_prime * 1000,Y_prime * 1000, I_Fraunhofer);shading interp;colorbar;axis equal;
pcolor(X_prime3 .* 1000,Y_prime3 .* 1000, I8);shading interp;colorbar;axis equal;
%circles(x,y,ring_outer_radius_geo * 1000,'color','black')
%drawellipse('Center',[X_prime2(501,501) * 1000, Y_prime2(501,501) * 1000],'SemiAxes',[ring_outer_radius_geo * 1000, ring_outer_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");

%drawellipse('Center',[X_prime2(501,501) * 1000, Y_prime2(501,501) * 1000],'SemiAxes',[ring_inner_radius_geo * 1000, ring_inner_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");
title('Fraunhofer Intensity after PRAC Axicon doublet at waist, zprop = 115mm, f = 100mm,Rexp = 3.5mm');
xlabel('xprime(mm)')
ylabel('yprime(mm)')
zlabel('Amplitude')


%plot axicon diffraction output
figure(4);
%pcolor(X_prime,Y_prime, I_U4);shading interp;colorbar;axis equal;
%pcolor(X_prime * 1000,Y_prime * 1000, I_Fraunhofer);shading interp;colorbar;axis equal;
pcolor(X_prime4 .* 1000,Y_prime4 .* 1000, I11);shading interp;colorbar;axis equal;
%circles(x,y,ring_outer_radius_geo * 1000,'color','black')
%drawellipse('Center',[X_prime2(501,501) * 1000, Y_prime2(501,501) * 1000],'SemiAxes',[ring_outer_radius_geo * 1000, ring_outer_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");

%drawellipse('Center',[X_prime2(501,501) * 1000, Y_prime2(501,501) * 1000],'SemiAxes',[ring_inner_radius_geo * 1000, ring_inner_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");
title('Fraunhofer Intensity after PRAC Axicon doublet at waist, zprop = 115mm, f = 100mm,Rexp = 3.5mm');
xlabel('xprime(mm)')
ylabel('yprime(mm)')
zlabel('Amplitude')


%plot axicon diffraction output
figure(5);
%pcolor(X_prime,Y_prime, I_U4);shading interp;colorbar;axis equal;
%pcolor(X_prime * 1000,Y_prime * 1000, I_Fraunhofer);shading interp;colorbar;axis equal;
pcolor(X_prime5 * 1000,Y_prime5 * 1000, I15);shading interp;colorbar;axis equal;
%circles(x,y,ring_outer_radius_geo * 1000,'color','black')
%drawellipse('Center',[X_prime2(501,501) * 1000, Y_prime2(501,501) * 1000],'SemiAxes',[ring_outer_radius_geo * 1000, ring_outer_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");

%drawellipse('Center',[X_prime2(501,501) * 1000, Y_prime2(501,501) * 1000],'SemiAxes',[ring_inner_radius_geo * 1000, ring_inner_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");
title('Fraunhofer Intensity after PRAC Axicon doublet at waist, zprop = 115mm, f = 100mm,Rexp = 3.5mm');
xlabel('xprime(mm)')
ylabel('yprime(mm)')
zlabel('Amplitude')

%plot axicon diffraction output
figure(6);
%pcolor(X_prime,Y_prime, I_U4);shading interp;colorbar;axis equal;
%pcolor(X_prime * 1000,Y_prime * 1000, I_Fraunhofer);shading interp;colorbar;axis equal;
pcolor(X_prime6 * 1000,Y_prime6 * 1000, I19);shading interp;colorbar;axis equal;
%circles(x,y,ring_outer_radius_geo * 1000,'color','black')
%drawellipse('Center',[X_prime2(501,501) * 1000, Y_prime2(501,501) * 1000],'SemiAxes',[ring_outer_radius_geo * 1000, ring_outer_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");

%drawellipse('Center',[X_prime2(501,501) * 1000, Y_prime2(501,501) * 1000],'SemiAxes',[ring_inner_radius_geo * 1000, ring_inner_radius_geo * 1000],'StripeColor','w','HandleVisibility',"off","InteractionsAllowed","none");
title('Fraunhofer Intensity after PRAC Axicon doublet at waist, zprop = 115mm, f = 100mm,Rexp = 3.5mm');
xlabel('xprime(mm)')
ylabel('yprime(mm)')
zlabel('Amplitude')






%plot axicon diffraction output 1D
figure(7);
%plot(x_prime * 1000,I_Fraunhofer(:,501)); 
plot(x_prime6 * 1000,I19(:,2001),'LineWidth',1); 
%xline(ring_outer_radius_geo * 1000, '-r','LineWidth',1)
%xline(ring_inner_radius_geo * 1000, '-g','LineWidth',1)
%plot(x_prime,I_U4(:,1000)); 
title('1D Fraunhofer Intensity after PRAC Axicon doublet at waist, zprop = 115mm, f = 100mm,Rexp = 3.5mm');
xlabel('xprime(mm)')
ylabel('Amplitude')





%%%%%%RING LOCATION%%%%%%%%%
%{
R_exprange = 0.1:0.1:5;
firstminloc = zeros(1,length(R_exprange));
firstminamp = zeros(1,length(R_exprange));

for i = 1:length(R_exprange)

    [E_GB, I_GB] = getGB(A,z,lambda,w0,zR, t, rho, Rz);
    
    U4_Fraunhofer_prac = propTF2dim_practical_Isabel(E_GB,X,Y,k,n,r_axi,alpha,f, R_exprange(i),e0,z_propafteroptics,X_prime,Y_prime);
    I_Fraunhofer_prac = abs(U4_Fraunhofer_prac.^2);
    Mall_ideal = max(I_Fraunhofer_prac,[],'all');
    I_Fraunhofer_prac = I_Fraunhofer_prac/Mall_ideal;
    
    linedata = I_Fraunhofer_prac(:,501);
    TF = islocalmin(linedata);
    minamp = linedata(TF);
    firstminloc = minloc(1);
    firstminamp = minamp(1);
end

linedata = I_Fraunhofer_prac(:,501);
TF = islocalmin(linedata);
minloc = x_prime(TF) * 1000;
minamp = linedata(TF);
firstminloc = minloc(1);
firstminamp = minamp(1);
%plot(x_prime * 1000,linedata,x_prime(TF) * 1000,linedata(TF),'r*')
plot(x_prime * 1000,linedata,firstminloc,firstminamp,'r*')
%}


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



function[U4,U5,U8,U11,U15,U19]=propTF2dim_DoublePracAxicon(u1,X,Y,X_prime1,Y_prime1,X_prime2,Y_prime2,X_prime3,Y_prime3,X_prime4,Y_prime4,X_prime5,Y_prime5,X_prime6,Y_prime6,k,n,f,r_axi,alpha,R_exp,e0, z1,z2,z3,z4,z5)
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


%%%Prop_to_lens
factor1 = exp(0.5.*1i .* k ./z1.*(X.^2+Y.^2));
U2 = fftshift(u1);
U3 = fftshift(fft2(U2.*factor1)); 
U4 =  -1i.*k./(2 .* pi).*exp(1i .* k .* z1)./z1 .* exp((1i.*k./(2.*z1)).*(X_prime1.^2 + Y_prime1.^2)).*U3;

%%%Transmission through lens
H1 = exp( -1i.* k.* 0.5.* (X_prime1.^2+Y_prime1.^2)./f);                     %transfer function lens
U5 = U4.*H1;

%%%Prop_to_expander
factor2 = exp(0.5.*1i .* k ./z2.*(X_prime1.^2+Y_prime1.^2));
U6 = fftshift(U5.*factor2);                                                         %field after first axicon
U7 = fftshift(fft2(U6));                                                      %FFT of U2
U8 =  -1i.*k./(2 .* pi).*exp(1i .* k .* z2)./z2 .* exp((1i.*k./(2.*z2)).*(X_prime2.^2 + Y_prime2.^2)).*U7;

%%%Prop_to_first_axicon
factor3 = exp(0.5.*1i .* k ./z3.*(X_prime3.^2+Y_prime3.^2));
U9 = fftshift(U8.*factor3);                                                         %field after first axicon
U10 = fftshift(fft2(U9));                                                      %FFT of U2
U11 =  -1i.*k./(2 .* pi).*exp(1i .* k .* z3)./z3 .* exp((1i.*k./(2.*z3)).*(X_prime4.^2 + Y_prime4.^2)).*U10;

%%%Transmission_through_axicon_first
H2=exp(-1i*k*(n-1)*(r_axi - sqrt(X_prime4.^2 + Y_prime4.^2))* tan(alpha));
thickness_function1 = e0 - R_exp .* (tan(alpha)).^2 .* sqrt(1 + (X_prime4.^2 + Y_prime4.^2)./(R_exp .* tan(alpha)).^2);
H3 = exp(i.*k.*(n-1).*thickness_function1);                                     %transfer function of a single axicon
U12 = U11.*H3;

%%%Prop_to_second_axicon
factor4 = exp(0.5.*1i .* k ./z4.*(X_prime4.^2 + Y_prime4.^2));
U13 = fftshift(U12.*factor4);                                                         %field after first axicon
U14 = fftshift(fft2(U13));                                                      %FFT of U2
U15 =  -1i.*k./(2 .* pi).*exp(1i .* k .* z4)./z4 * exp((1i.*k./(2.*z4)).*(X_prime5.^2 + Y_prime5.^2)).*U14;


%%%Transmission_through_axicon_second
H4=exp(-1i*k*(n-1)*(r_axi - sqrt(X_prime5.^2 + Y_prime5.^2))* tan(alpha));
thickness_function2 = e0 - R_exp .* (tan(alpha)).^2 .* sqrt(1 + (X_prime5.^2 + Y_prime5.^2)./(R_exp .* tan(alpha)).^2);
H5 = exp(i.*k.*(n-1).*thickness_function2);                                     %transfer function of a single axicon
U16 = U15.*H5; 

%%%Prop_to_observ
factor5 = exp(0.5.*1i .* k ./z5.*(X_prime5.^2 + Y_prime5.^2));
U17=fftshift(U16.*factor5);                            %multiply U2: field after optics H
U18=fftshift(fft2(U17));             %inverse fft and center observ field

U19 = -1i.*k./(2 .* pi).*exp(1i .* k .* z5)./z5 .* exp((1i.*k./(2.*z5)).*(X_prime6.^2 + Y_prime6.^2)).*U18;
end


