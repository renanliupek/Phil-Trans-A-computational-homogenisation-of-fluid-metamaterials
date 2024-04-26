% DESCRIPTION
% Renan Liupekevicius Carnielli TU/e
% start in 20-12-2023

% Compute frequency-dependent material properties.



%% ADAPTING NOTATION TO PAPER 2 (ROYAL TRANSACTIONS A)

% local resonance material properties
omega_lr = 2*pi*freqs;
c        =  heL;
gamma    =  hel;
d        =  hvL;
delta    =  hvl;

% d        =  1.0e-03 * [ 0.0000*e1 - 0.0000*e2
%                        -0.1966*e1 - 0.1966*e2
%                        -0.1966*e1 + 0.1966*e2
%                        -0.0000*e1 + 0.0000*e2
%                         0.0000*e1 + 0.0000*e2];
% delta    =  V*d;

% longwave material properties

Ce       = heC;
De       = heD;

CU       = hvC;
DU       = hvD;
HU       = hvH;



%% EFFECTIVE BULK MODULUS

f     = (0:.1:1084)'; %column

omega = 2*pi*f;

KM       =   -1./( Ce + ...
            -omega.^2./( -omega.^2+omega_lr(1)^2 )  * c(1)*gamma(1) ...
            -omega.^2./( -omega.^2+omega_lr(2)^2 )  * c(2)*gamma(2) ...
            -omega.^2./( -omega.^2+omega_lr(3)^2 )  * c(3)*gamma(3) ...
            -omega.^2./( -omega.^2+omega_lr(4)^2 )  * c(4)*gamma(4) ...
            -omega.^2./( -omega.^2+omega_lr(5)^2 )  * c(5)*gamma(5) );

SM       =    ( De*ones(length(omega),1) + ...
            -omega.^2./( -omega.^2+omega_lr(1)^2 )  * c(1)*delta(1) ...
            -omega.^2./( -omega.^2+omega_lr(2)^2 )  * c(2)*delta(2) ...
            -omega.^2./( -omega.^2+omega_lr(3)^2 )  * c(3)*delta(3) ...
            -omega.^2./( -omega.^2+omega_lr(4)^2 )  * c(4)*delta(4) ...
            -omega.^2./( -omega.^2+omega_lr(5)^2 )  * c(5)*delta(5) );

SM_T     =   ( CU*ones(length(omega),1) + ...
            -omega.^2./( -omega.^2+omega_lr(1)^2 )  * d(1)*gamma(1) ...
            -omega.^2./( -omega.^2+omega_lr(2)^2 )  * d(2)*gamma(2) ...
            -omega.^2./( -omega.^2+omega_lr(3)^2 )  * d(3)*gamma(3) ...
            -omega.^2./( -omega.^2+omega_lr(4)^2 )  * d(4)*gamma(4) ...
            -omega.^2./( -omega.^2+omega_lr(5)^2 )  * d(5)*gamma(5) );

RHOM_inv =    -(-omega.^2*DU + HU*ones(length(omega),1) +...
                -omega.^4./( -omega.^2+omega_lr(1)^2 )  * d(1)*delta(1) ...
                -omega.^4./( -omega.^2+omega_lr(2)^2 )  * d(2)*delta(2) ...
                -omega.^4./( -omega.^2+omega_lr(3)^2 )  * d(3)*delta(3) ...
                -omega.^4./( -omega.^2+omega_lr(4)^2 )  * d(4)*delta(4) ...
                -omega.^4./( -omega.^2+omega_lr(5)^2 )  * d(5)*delta(5) );

% get component 1 and 11
SM1    = dot(SM,e1);
SM2    = dot(SM,e2);
SM2_T    = dot(SM_T,e2);
SM1_T    = dot(SM_T,e1);
SM2    = dot(SM,e2);
RHOM11 = 1./dot(e1,RHOM_inv,e1);
RHOM22 = 1./dot(e2,RHOM_inv,e2);


%% Plot normalized

% figure(1)
% plot(f,KM/KM(1))
% ylim([-10 10])
% grid on
% title('effective macroscopic bulk modulus')
% 
% figure(2)
% plot(f,RHOM11/RHOM11(1))
% % hold on
% % plot(f,RHOM22/RHOM22(1))
% ylim([-100 100])
% hold off
% grid on
% % legend 
% title('effective macroscopic density')
% 
% % note diff between x and y because De is different in x and y (zero
% % machine)
% figure(3)
% plot(f,SM1/SM1(1))
% % hold on
% % plot(f,SM2/SM2(1))
% ylim([-10 10])
% hold off
% grid on
% % legend 
% title('effective macroscopic Willis coupling')
% 
% % figure(4)
% % plot(f,SM1_T/SM1_T(1))
% % hold on
% % plot(f,SM2_T/SM2_T(1))
% % ylim([-10 10])
% % hold off
% % grid on
% % legend 
% % title('effective macroscopic Willis coupling(transpose)')


%% logscale plot (to neglect SM plot in the paper)

% figure(1)
% plot(f,log10(abs(KM)))
% % ylim([-10 10])
% grid on
% title('effective macroscopic bulk modulus')
% 
% figure(2)
% plot(f,log10(abs(RHOM11)))
% % hold on
% % plot(f,RHOM22/RHOM22(1))
% % ylim([-100 100])
% % hold off
% grid on
% % legend 
% title('effective macroscopic density')
% %
% 
% % note diff between x and y because De is different in x and y (zero
% % machine)
% figure(3)
% plot(f,log10(abs(SM1)))
% % hold on
% % plot(f,SM2/SM2(1))
% % ylim([-10 10])
% hold off
% grid on
% % legend 
% title('effective macroscopic Willis coupling')


%% Plot for paper



% figure(1)
% plot(KM(1)*ones(length(f),1),f,'black--', 'LineWidth', 2)
% hold on
% plot(KM,f,'red', 'LineWidth', 2)
% hold off
% xlim([-10*KM(1) 10*KM(1)])
% ylim([0 730])
% grid on
% title('effective macroscopic bulk modulus')
% 
% figure(2)
% plot(RHOM11(1)*ones(length(f),1),f,'black--', 'LineWidth', 2)
% hold on
% plot(RHOM11,f,'red', 'LineWidth', 2)
% hold off
% % xlim([-1*RHOM11(1) 2*RHOM11(1)])
% ylim([0 730])
% hold off
% grid on
% % legend 
% title('effective macroscopic density')

% figure(3)
% % plot(f,SM1/SM1(1))
% plot(SM1(1)*ones(length(f),1),f,'black--', 'LineWidth', 2)
% hold on
% plot(SM2(1)*ones(length(f),1),f,'black--', 'LineWidth', 2)
% plot(SM1,f,'red', 'LineWidth', 2)
% plot(SM2,f,'blue', 'LineWidth', 2)
% hold off
% xlim([10*SM1(1) -10*SM1(1)])
% ylim([0 730])
% hold off
% grid on
% % legend 
% title('effective macroscopic willis coupling coefficient')


% % figure(3)
% % plot(f,SM1/SM1(1))
% 
% absS_LW = sqrt(SM1(1)^2+SM2(1)^2);
% 
% absS = sqrt(SM1.^2+SM2.^2);
% hold on
% plot(log10(absS_LW)*ones(length(f),1),f,'black--', 'LineWidth', 2)
% plot(log10(absS),f,'blue', 'LineWidth', 2)
% hold off
% % xlim([10*absS_LW -10*absS_LW])
% ylim([0 740])
% hold off
% grid on
% % legend 
% title('effective macroscopic willis coupling coefficient')


%% x is freq axis

% figure(3)
% plot(f,SM1/SM1(1))

absS_LW = sqrt(SM1(1)^2+SM2(1)^2);

absS = sqrt(SM1.^2+SM2.^2);
hold on
plot(f,log10(absS_LW)*ones(length(f),1),'black--', 'LineWidth', 2)
plot(f,log10(absS),'blue', 'LineWidth', 2)
hold off
% xlim([10*absS_LW -10*absS_LW])
xlim([0 740])
hold off 
grid on
% legend 
title('effective macroscopic willis coupling coefficient')
