% Script to load the results into Octave
clear all
RHO_all = load("../RUN/rho.out");
Ux1_all = load("../RUN/ux1.out");
Ux2_all = load("../RUN/ux2.out");
P_all = load("../RUN/P.out");
E_all = load("../RUN/E.out");
A_all = load("../RUN/A.out");
ux_all = load("../RUN/ux.out");
uz_all = load("../RUN/uz.out");
vort_all = load("../RUN/vort.out");
X1 = load("../RUN/x1.out");
X2 = load("../RUN/x2.out");
X = load("../RUN/x.out");
Z = load("../RUN/z.out");
%

% #######################################################
% Read some of the parameters from init.params
nx1 = 50
nx2 = 50
outf = 200
% #######################################################
%
nt = length(RHO_all)/(nx1*nx2)
for i=1:nx1
for j=1:nx2
x1(i,j) = X1(nx2*(i-1)+j);
x2(i,j) = X2(nx2*(i-1)+j);
x(i,j) = X(nx2*(i-1)+j);
z(i,j) = Z(nx2*(i-1)+j);
endfor
endfor
for k=1:nt
for i=1:nx1
for j=1:nx2
rho(i,j,k) = RHO_all(nx1*nx2*(k-1) + nx2*(i-1)+j);
ux1(i,j,k) = Ux1_all(nx1*nx2*(k-1) + nx2*(i-1)+j);
ux2(i,j,k) = Ux2_all(nx1*nx2*(k-1) + nx2*(i-1)+j);
p(i,j,k) = P_all(nx1*nx2*(k-1) + nx2*(i-1)+j);
E(i,j,k) = E_all(nx1*nx2*(k-1) + nx2*(i-1)+j);
A(i,j,k) = A_all(nx1*nx2*(k-1) + nx2*(i-1)+j);
ux(i,j,k) = ux_all(nx1*nx2*(k-1) + nx2*(i-1)+j);
uz(i,j,k) = uz_all(nx1*nx2*(k-1) + nx2*(i-1)+j);
vort(i,j,k) = vort_all(nx1*nx2*(k-1) + nx2*(i-1)+j);
endfor
endfor
endfor
little_e = E./rho .- 0.5.*(ux1.**2.+ux2.**2);
umag=(ux1.**2 .+ux2.**2).**0.5;
%PovPinfi = (1/Pinf).*Pi;
%
T_all = load("../RUN/dt.out");
T = T_all(:,1);
dt = T_all(:,2);
R1_all = load('../RUN/interface.out');
R1 = R1_all(:,2);
dRdt(1) = 0.0
for i=2:length(R1)
dRdt(i) = (R1(i)-R1(i-1))/dt(i);
end
%
