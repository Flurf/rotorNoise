clc
clear all
%##################################################
%% propeller data definition
Rt = 8.5344;
R0 = 0.8;
b = 1;  %number of blades, 1 is correct just in hovering, total thrust has to be multiplied for 4
rn = 80;
chord = 0.4166616;
rootpitch = 13.908;
twist = 8;
airfoil = 'NACA0012';


rotor = Propeller(Rt,R0,b,rn,chord,rootpitch,twist,airfoil);
% %##################################################
%% flight conditions data definition

Texp = 53000; %[N]
RPM = 212;
Vf = 0;       %forward flight
deltapsi = 360;     %360 is correct for hovering
calc = TTcalculator(rotor,RPM,Vf,deltapsi,Texp);
%%##################################################
%% acoustic solvers definitions
c=340 ;                     % speed of sound, [m/s]

n = 10;                     % harmonics, total number
b=4 ;                       % no of blades
omega=RPM*(2*pi/60) ;       % rotational speed, [rad/s]
R = 228.6;                  %observer distance, [m]
sigma= -42;                  
teta= -63.43;

deltapsi = 1;
psi = linspace(0,360,360/deltapsi);                %hovering
l = zeros(length(rotor.r),length(psi));
for i = 1:length(psi)
        l(:,i) = calc.L(1,:,1,1);
end
beta = zeros(length(rotor.r),length(psi));
for j = 1:length(psi)
    beta(:,j) = rotor.beta;
end
[p_rms, db] = schlegel(rotor.r,chord,psi,beta,l,R,sigma,teta,n,b,omega,c);

% l = zeros(length(rotor.r),length(calc.psi));   %forward flight
% for i = 1:length(calc.psi)
%         l(:,i) = calc.L(1,:,i,1);
% end
% for j = 1:length(calc.psi)
%     beta(:,j) = rotor.beta;
% end
% [p_rms, db] = schlegel(rotor.r,chord,calc.psi,beta,l,R,sigma,teta,n,b,omega,c);

% custom loads, from schlegel experimental reference
psi2pa = 6894.76;
inch2m = 0.0254;
dr = 0.01;   
[ls, beta_harm, r, Lphi60] =load_2_order(dr, deltapsi);
r = r.*Rt;      %radius in meters


beta_custom = zeros(length(r),length(psi));
for i = 1:length(r)
    for j = 1:length(psi)
        beta_custom(i,j) = beta_harm(j) - twist*((r(i) - R0)/(Rt - R0));
    end
end
ls = ls.*(psi2pa * inch2m);   %from [lb/inch] to [Pa * m]
%[p_rms_custom, db_custom] = schlegel(r,chord,psi,beta_custom,ls,R,sigma,teta,n,b,omega,c);
%##################################################
rsch = [0.093,0.25,0.4,0.55,0.75,0.85,0.90,0.95,1].*Rt;
% figure (2)
% plot(r,ls(:,60),'blue','LineWidth', 1.5)
% hold on
% for i = 1:length(psi)
%     plot(r,l(:,i),'Color','black')
%     line(rsch(i),Lphi60(i).*psi2pa*inch2m,'LineWidth', 1.5)
%     hold on
% end

legend('experimental data')
% plot([0.093,0.25,0.4,0.55,0.75,0.85,0.90,0.95,1].*Rt,Lphi60(i).*psi2pa*inch2m,'o-','LineWidth','1.5')
ylabel('Load [Pa * m]','FontSize', 16)
xlabel('radius [m]','FontSize', 16)
grid on
lgd= legend;
lgd.FontSize = 16;
xlim([0 9])
ylim([-100 8300])
ax = gca; 
ax.FontSize = 16; 


