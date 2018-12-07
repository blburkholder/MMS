% load('MMS_matfiles/cs_Bpar.mat')
% load('MMS_matfiles/cs_Bperp.mat')
% load('MMS_matfiles/cs_walen1.mat')
% load('MMS_matfiles/cs_walen2.mat')
% load('MMS_matfiles/cs_walen3.mat')
% load('MMS_matfiles/walegood1.mat')
% load('MMS_matfiles/walegood2.mat')
% load('MMS_matfiles/walegood3.mat')
% load('MMS_matfiles/cs_eigrats.mat')
% load('MMS_matfiles/cs_vx.mat')
% load('MMS_matfiles/cs_jumps_tot')
% load('MMS_matfiles/v_jumps_totx')
% load('MMS_matfiles/va_jumps_totx')
% load('MMS_matfiles/v_jumps_toty')
% load('MMS_matfiles/va_jumps_toty')
% load('MMS_matfiles/v_jumps_totz')
% load('MMS_matfiles/va_jumps_totz')
% load('MMS_matfiles/cs_jumps1')
% load('MMS_matfiles/cs_jumps2')
% load('MMS_matfiles/cs_jumps3')
% load('MMS_matfiles/cs_shear_tot')
% load('MMS_matfiles/cs_shear1')
% load('MMS_matfiles/cs_shear2')
% load('MMS_matfiles/cs_shear3')
% load('MMS_matfiles/cs_shflow_stab')
% %load('ion_inlength')
% load('MMS_matfiles/ion_beta')
% load('MMS_matfiles/tags')
load('cs_Bpar.mat')
load('cs_Bperp.mat')
load('cs_walen1.mat')
load('cs_walen2.mat')
load('cs_walen3.mat')
load('walegood1.mat')
load('walegood2.mat')
load('walegood3.mat')
load('cs_eigrats.mat')
load('cs_vx.mat')
load('cs_jumps_tot')
load('v_jumps_totx')
load('va_jumps_totx')
load('v_jumps_toty')
load('va_jumps_toty')
load('v_jumps_totz')
load('va_jumps_totz')
load('cs_jumps1')
load('cs_jumps2')
load('cs_jumps3')
load('cs_shear_tot')
load('cs_shear1')
load('cs_shear2')
load('cs_shear3')
load('cs_shflow_stab')
%load('ion_inlength')
load('ion_beta')
load('tags')
mu = 4*pi*10^(-7);

v_jumps_tot = sqrt(v_jumps_totx.^2+v_jumps_toty.^2+v_jumps_totz.^2);

cs_walen3(cs_walen3 ~=0) = (cs_walen3(cs_walen3 ~= 0)+1)/4;
cs_walen3 = -cs_vx.*cs_walen3;
cs_walen2(cs_walen2 ~=0) = (cs_walen2(cs_walen2 ~= 0)+1)/4;
cs_walen2 = -cs_vx.*cs_walen2;
cs_walen1(cs_walen1 ~=0) = (cs_walen1(cs_walen1 ~= 0)+1)/4;
cs_walen1 = -cs_vx.*cs_walen1;
cs_Bmag = sqrt(cs_Bpar.^2+cs_Bperp.^2);

% figure; hold on
% scatter(cs_walen3(walegood3 >= 0.5),ion_inlength(walegood3 >= 0.5)/1000,'b.');
% ylabel('c/\omega_{pi} km');
% title('Solar Wind Current Sheets');
% xlabel('length 0.9 < \Delta v / \Delta v_A < 1.1 (km)');

figure; hold on
scatter(sqrt(1.67/1.75)*va_jumps_totx,v_jumps_totx,'.')
scatter(sqrt(1.67/1.75)*va_jumps_toty,v_jumps_toty,'.')
scatter(sqrt(1.67/1.75)*va_jumps_totz,v_jumps_totz,'.')
xlabel('\Delta v_{Ai}')
ylabel('\Delta v_i')
legend('x','y','z')
title('Solar Wind Current Sheets')
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';

%slope
B = [ones(3*length(va_jumps_totx),1),[sqrt(1.67/1.75)*va_jumps_totx;sqrt(1.67/1.75)*va_jumps_toty;sqrt(1.67/1.75)*va_jumps_totz]]\...
    [v_jumps_totx;v_jumps_toty;v_jumps_totz]

% figure; hold on
% scatter(cs_Bpar./cs_Bmag,cs_shflow_stab/1000,'b.');
% ylabel('\Delta v_{anti-||}/v_A km/s');
% title('Solar Wind Current Sheets');
% xlabel('B_n/|B|');

figure
hold on;
scatter(cs_walen3(walegood3>=0.5),cs_jumps3(walegood3>=0.5),[100],...
    cs_Bpar(walegood3>=0.5)./cs_Bmag(walegood3>=0.5),'r.');
scatter(cs_walen3(cs_eigrats>10 & walegood3>=0.5 & cs_Bpar./cs_Bmag < 0.2),cs_jumps3(cs_eigrats>10 & walegood3>=0.5 & cs_Bpar./cs_Bmag < 0.2),[100],...
    cs_Bpar(cs_eigrats>10 & walegood3>=0.5 & cs_Bpar./cs_Bmag < 0.2)./cs_Bmag(cs_eigrats>10 & walegood3>=0.5 & cs_Bpar./cs_Bmag < 0.2),'b.');
scatter(cs_walen3(cs_eigrats>10 & walegood3>=0.5 & cs_Bpar./cs_Bmag < 0.2 & ion_beta < 10),cs_jumps3(cs_eigrats>10 & walegood3>=0.5 & cs_Bpar./cs_Bmag < 0.2 & ion_beta < 10),[100],...
    cs_Bpar(cs_eigrats>10 & walegood3>=0.5 & cs_Bpar./cs_Bmag < 0.2 & ion_beta < 10)./cs_Bmag(cs_eigrats>10 & walegood3>=0.5 & cs_Bpar./cs_Bmag < 0.2 & ion_beta < 10),'g.');
ylabel('\Delta B_{walen} (nT)');
ax = gca;
title('Solar Wind Current Sheets');
ax.YScale = 'log';
xlabel('length 0.9 < \Delta v / \Delta v_A < 1.1 (km)');
legend('R^2 > 0.5','R^2 > 0.5 & B_n/|B| < 0.2 & \lambda_2/\lambda_1 > 10','R^2 > 0.5 & B_n/|B| < 0.2 & \lambda_2/\lambda_1 > 10 & \beta < 10')
hold off;
%axis([0 800 0.7 30]);
% 
% figure
% [slop1,inter1] = logfit(cs_walen3,cs_jumps3,'logy');
% figure
% [slop2,inter2] = logfit([cs_walen3 cs_walen2 cs_walen1],[cs_jumps3 cs_jumps2 cs_jumps1],'logy');

figure
hold on;
scatter(cs_walen1(walegood1>=0.5),cs_jumps1(walegood1>=0.5),[100],...
    cs_Bpar(walegood1>=0.5)./cs_Bmag(walegood1>=0.5),'r.');
scatter(cs_walen1(cs_eigrats>10 & walegood1>=0.5 & cs_Bpar./cs_Bmag < 0.2),cs_jumps1(cs_eigrats>10 & walegood1>=0.5 & cs_Bpar./cs_Bmag < 0.2),[100],...
    cs_Bpar(cs_eigrats>10 & walegood1>=0.5 & cs_Bpar./cs_Bmag < 0.2)./cs_Bmag(cs_eigrats>10 & walegood1>=0.5 & cs_Bpar./cs_Bmag < 0.2),'b.');
scatter(cs_walen1(cs_eigrats>10 & walegood1>=0.5 & cs_Bpar./cs_Bmag < 0.2 & ion_beta < 10),cs_jumps1(cs_eigrats>10 & walegood1>=0.5 & cs_Bpar./cs_Bmag < 0.2 & ion_beta < 10),[100],...
    cs_Bpar(cs_eigrats>10 & walegood1>=0.5 & cs_Bpar./cs_Bmag < 0.2 & ion_beta < 10)./cs_Bmag(cs_eigrats>10 & walegood1>=0.5 & cs_Bpar./cs_Bmag < 0.2 & ion_beta < 10),'g.');
ylabel('\Delta B_{walen} (nT)');
ax = gca;
title('Solar Wind Current Sheets');
ax.YScale = 'log';
xlabel('length 0.75 < \Delta v / \Delta v_A < 1.25 (km)');
legend('R^2 > 0.5','R^2 > 0.5 & B_n/|B| < 0.2 & \lambda_2/\lambda_1 > 10','R^2 > 0.5 & B_n/|B| < 0.2 & \lambda_2/\lambda_1 > 10 & \beta < 10')
hold off;
%axis([0 800 0.7 30]);

figure
hold on
s1 = scatter(cs_walen1(walegood1 >= 0.5),cs_jumps1(walegood1 >= 0.5),'r.');
s2 = scatter(cs_walen2(walegood2 >= 0.5),cs_jumps2(walegood2 >= 0.5),'g.');
s3 = scatter(cs_walen3(walegood3 >= 0.5),cs_jumps3(walegood3 >= 0.5),'b.');
%plot([0 max(cs_walen1)],(10^inter1)*(10^slop1).^([0 max(cs_walen1)]),'k--')
%plot([0 max(cs_walen1)],(10^inter2)*(10^slop2).^([0 max(cs_walen1)]),'k')
ylabel('\Delta B_{walen} (nT)');
ax = gca;
title('Solar Wind Current Sheets');
ax.YScale = 'log';
xlabel('length 1-\epsilon < \Delta v / \Delta v_A < 1+\epsilon (km)');
hold off;
%axis([0 1600 0.01 11])
legend('\epsilon = 0.25','\epsilon = 0.2','\epsilon = 0.1','Location','southeast');

% figure
% hold on
% scatter(cs_Bperp(cs_eigrats > 10)./cs_Bmag(cs_eigrats > 10),cs_Bpar(cs_eigrats > 10),'b.')
% ylabel('B_n (nT)')
% ax = gca;
% title('Solar Wind Current Sheets')
% ax.YScale = 'log';
% xlabel('B_t/|B|')
% %axis([0.3 1 .01 30])

figure; hold on;
scatter(cs_jumps1(walegood1 >= 0.5),cs_jumps_tot(walegood1 >= 0.5),'r.')
scatter(cs_jumps2(walegood2 >= 0.5),cs_jumps_tot(walegood2 >= 0.5),'g.')
scatter(cs_jumps3(walegood3 >= 0.5),cs_jumps_tot(walegood3 >= 0.5),'b.')
ylabel('\Delta B_{cs} (nT)')
title('Solar Wind Current Sheets')
xlabel('\Delta B_{walen} (nT)')
ax = gca;
ax.YScale = 'log';
ax.XScale = 'log';
legend('\epsilon = 0.25','\epsilon = 0.2','\epsilon = 0.1','Location','Southeast')

figure; hold on;
scatter(v_jumps_tot,cs_jumps_tot,'b.')
ylabel('\Delta B_{cs} (nT)')
title('Solar Wind Current Sheets')
xlabel('\Delta v (km)')
ax = gca;
ax.YScale = 'log';
ax.XScale = 'log';

% figure
% scatter(cs_shear_tot,cs_jumps_tot,'b.')
% ylabel('\Delta B_{cs} (nT)')
% title('Solar Wind Current Sheets')
% xlabel('cos\theta_{cs}')
% ax = gca;
% ax.YScale = 'log';
