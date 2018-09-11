clc;
clear all;
close all;

%Simulation Box Measure
offset = 2;
ind = 0;
file = fopen('Ringkern.txt','r');
Rk_Simulation = textscan(file, '%s %s %s %s %s');
fclose(file);




for param_num = 2:4
    for length_num = 1:numel(Rk_Simulation{:,param_num})
        A = Rk_Simulation{:,param_num}{length_num};
        ind = strfind(A, '"');
        A(ind) = [];
        ind = strfind(A, ',');
        A(ind) = [];
        Num(length_num,param_num-1) = str2double(A);
    end    
end

%Reverse engineering of Jens' measured data

f = Num(:,1).*1e6;
w = 2.*pi.*f;
mus = Num(:,2);
muss = Num(:,3);
l = 0.025;
ra = 0.248;
ri = 0.126;

Lrk = mus.*l.*log(ra./ri)./2./pi;
Rrk = muss.*w.*Lrk./mus;
R_test = 3.5069e4;

Zrk = Rrk + 1i.*Lrk.*w;

L_test = 5.4e-7;
C_test = 11.2e-12; %10.74e-12

Zc_test = -1i./(C_test.*w);
Zl_test = 1i.*w.*L_test;
par = @(x,y)((x.*y)./(x+y));

Zges = par(Zrk + Zl_test,Zc_test);


%comparison with the own data
Box_measure4_1104 = load('Mess4_110418.dat');
f_measure = Box_measure4_1104(:,1);
Z_measure = Box_measure4_1104(:,2) + 1i.*Box_measure4_1104(:,3);
w_meas = 2.*pi.*f_measure;

%compensation of the testbox
%Zrk_meas = (Z_measure.*1i.*w_meas.*L_test.*(1./R_test + 1i.*w_meas.*C_test)-1i.*w_meas.*L_test)./...
%    (1 - Z_measure.*(1./R_test + 1i.*w_meas.*C_test));
Zrk_meas = (Z_measure.*(1 - w_meas.^2.*L_test.*C_test) - 1i.*w_meas.*L_test)./...
            (1 - 1i.*w_meas.*Z_measure.*C_test);
Zrk_meas_zaehler = Z_measure.*(w_meas.^2.*R_test.*L_test.*C_test -1i.*w_meas.*L_test - R_test)...
    + 1i.*w_meas.*R_test.*L_test;
Zrk_meas_nenner = Z_measure + 1i.*w_meas.*R_test.*C_test.*Z_measure - R_test;        
        
Lrk_meas = imag(Zrk_meas)./w_meas;
Rrk_meas = real(Zrk_meas);

%figure
%plot(f_measure,real(Zrk_meas),f_measure,imag(Zrk_meas),f_measure, abs(Zrk_meas), f, abs(Zges));
%legend('real','imag');

figure
plot(f,imag(Zrk))
title('Imag Zrk Data Denys')

figure
title('Comparison Zrk');
plot(f_measure,abs(Zrk_meas),f, abs(Zrk));
legend('measurement + testbox compensation','Reverse engineered data from Denys');
ylim([0,7e2])

figure
plot(f_measure(1:1e3),abs(Zrk_meas(1:1e3)),f_measure(1e3:end), abs(Zrk_meas_zaehler(1e3:end)./Zrk_meas_nenner(1e3:end)));
legend('measurement + testbox compensation', 'another formula');
ylim([0,7e2]);

mus_meas = Lrk_meas.*2.*pi./(l.*log(ra./ri).*4.*pi.*1e-7);
muss_meas = Rrk_meas.*mus_meas./(w_meas.*Lrk_meas);

figure
plot(f, mus, f, muss, f_measure,mus_meas,f_measure, muss_meas);
%xlim([0,3e6]);
legend('mus','muss','mus_meas','muss_meas');

figure
plot(f_measure,mus_meas,f_measure, muss_meas);


%Save material
offset = 2;
ind = 0;
A = [f_measure.*1e-6, mus_meas,muss_meas]';
file = fopen('RK_ME990.txt','w');
fprintf(file, '%f %f %f \n',A);
fclose(file);


