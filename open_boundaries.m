% load('MMS_matfiles_ion_cadence/walegood.mat')
% load('MMS_matfiles_ion_cadence/cs_vx.mat')
% load('MMS_matfiles_ion_cadence/cs_vy.mat')
% load('MMS_matfiles_ion_cadence/cs_vz.mat')
% load('MMS_matfiles_ion_cadence/cs_jumps_tot')
% load('MMS_matfiles_ion_cadence/v_jumps_totx')
% load('MMS_matfiles_ion_cadence/va_jumps_totx')
% load('MMS_matfiles_ion_cadence/v_jumps_toty')
% load('MMS_matfiles_ion_cadence/va_jumps_toty')
% load('MMS_matfiles_ion_cadence/v_jumps_totz')
% load('MMS_matfiles_ion_cadence/va_jumps_totz')
% load('MMS_matfiles_ion_cadence/v_jumps1')
% load('MMS_matfiles_ion_cadence/v_jumps2')
% load('MMS_matfiles_ion_cadence/va_jumps1')
% load('MMS_matfiles_ion_cadence/va_jumps2')
% load('MMS_matfiles_ion_cadence/tags')
% load('MMS_matfiles_ion_cadence/cs_loc')
% load('MMS_matfiles_ion_cadence/cs_Bpar')    
% load('MMS_matfiles_ion_cadence/cs_Bperp')    
% load('MMS_matfiles_ion_cadence/cs_eigrats')    
% load('MMS_matfiles_ion_cadence/cs_n1_nt')
% load('MMS_matfiles_ion_cadence/cs_cond')    
% load('MMS_matfiles_ion_cadence/cs_rho_enh')
% load('MMS_matfiles_ion_cadence/cs_temp_enh') 
% load('MMS_matfiles_ion_cadence/cs_norm_vx')
% load('MMS_matfiles_ion_cadence/cs_wal_len')
% load('MMS_matfiles_ion_cadence/cs_beta')

load('walegood.mat')
load('cs_vx.mat')
load('cs_vy.mat')
load('cs_vz.mat')
load('cs_jumps_tot')
load('v_jumps_totx')
load('va_jumps_totx')
load('v_jumps_toty')
load('va_jumps_toty')
load('v_jumps_totz')
load('va_jumps_totz')
load('v_jumps1')
load('v_jumps2')
load('va_jumps1')
load('va_jumps2')
load('tags')
load('cs_loc')
load('cs_Bpar')    
load('cs_Bperp')    
load('cs_eigrats')    
load('cs_n1_nt')
load('cs_cond')    
load('cs_rho_enh')
load('cs_temp_enh')  
load('cs_norm_vx')
load('cs_wal_len')
load('cs_beta')

mu = 4*pi*10^(-7);
cs = 363;

ion_length = (3e8*1./sqrt(min(cs_rho_enh,[],2)*(1.6e-19)^2/(1.67e-27*8.9e-12)))/1000;

% va_jumps_totx = (10^-12)*va_jumps_totx./sqrt(mu*(1.67*10^(-27)));
% va_jumps_toty = (10^-12)*va_jumps_toty./sqrt(mu*(1.67*10^(-27)));
% va_jumps_totz = (10^-12)*va_jumps_totz./sqrt(mu*(1.67*10^(-27)));

v_jumps_tot = sqrt(v_jumps_totx.^2+v_jumps_toty.^2+v_jumps_totz.^2);
cs_v = sqrt(cs_vx.^2+cs_vy.^2+cs_vz.^2);
cs_Bmag = sqrt(cs_Bpar.^2+cs_Bperp.^2);

% figure; hold on
% cooond = abs(cs_n1_nt)>0.7 & cs_eigrats > 5;
% histogram((cs_inds(cs_v<450 & cooond,1) - cs_inds(cs_v<450 & cooond,2))/128.*...
%     cs_vx(cs_v<450 & cooond).*sin(cs_norm_vx(cs_v<450 & cooond))./ion_length(cs_v<450 & cooond),16,'BinLimits',[0 80])
% histogram((cs_inds(cs_v>450 & cooond,1) - cs_inds(cs_v>450 & cooond,2))/128.*...
%     cs_vx(cs_v>450 & cooond).*sin(cs_norm_vx(cs_v>450 & cooond))./ion_length(cs_v>450 & cooond),16,'BinLimits',[0 80])
% legend('v<450','v>450')
% xlabel('\lambda_i')
% title('current sheet width')
% 
% figure; hold on
% histogram((cs_inds(cs_v<450,1) - cs_inds(cs_v<450,2))/128.*...
%     cs_vx(cs_v<450)./ion_length(cs_v<450),16)%,'BinLimits',[0 80])
% histogram((cs_inds(cs_v>450,1) - cs_inds(cs_v>450,2))/128.*...
%     cs_vx(cs_v>450)./ion_length(cs_v>450),16)%,'BinLimits',[0 80])
% legend('v<450','v>450')
% xlabel('\lambda_i')
% title('current sheet width no angle')

figure
histogram(cs_jumps_tot)
title('current sheets')
xlabel('$(\Delta B_x^2+\Delta B_y^2+\Delta B_z^2)^{1/2}$','interpreter','latex')

figure
histogram(cs_beta(cs_beta < 100))
title('plasma Beta')
xlabel('$\beta$','interpreter','latex')

figure
scatter(cs_loc(:,1)/6371,cs_loc(:,2)/6371,100,cs_loc(:,3)/6371,'.'); hold on;
colormap(jet)
h = colorbar;
ylabel(h, 'z (R_E)')
scatter(0,0,[1000],'b.')
xlabel('x (R_E)')
ylabel('y (R_E)')
theta = -pi/1.5:pi/32:pi/1.5;
r = (2./(1+cos(theta))).^0.6;
plot(7*r.*cos(theta),7*r.*sin(theta),'k--')
plot(10*r.*cos(theta),10*r.*sin(theta),'k--')
plot(13*r.*cos(theta),13*r.*sin(theta),'k--')
title('solar wind current sheets')
annotation('textbox',[.2 .5 .25 .25],'EdgeColor','w','String','r_0 = 7','FitBoxToText','on');
annotation('textbox',[.15 .45 .37 .37],'EdgeColor','w','String','r_0 = 10','FitBoxToText','on');
annotation('textbox',[.12 .42 .48 .48],'EdgeColor','w','String','r_0 = 13','FitBoxToText','on');

[ord_vajump_tot,sortorder] = sort([va_jumps_totx,va_jumps_toty,va_jumps_totz],2);
unord_vjump_tot = [v_jumps_totx,v_jumps_toty,v_jumps_totz];
ord_vajump = zeros(cs,3); ord_vjump = ord_vajump; or_jump_tot = ord_vajump;
for i = 1:cs
    ord_vajump(i,:) = va_jumps2(i,sortorder(i,:));
    ord_vjump(i,:) = v_jumps2(i,sortorder(i,:));
    ord_vjump_tot(i,:) = unord_vjump_tot(i,sortorder(i,:));
end

figure; hold on;
diff_S = (10^4*cs_temp_enh(:,1)./cs_rho_enh(:,1).^(2/3)-10^4*cs_temp_enh(:,2)./cs_rho_enh(:,2).^(2/3));
diff_temp = 2*(cs_temp_enh(:,1)-cs_temp_enh(:,2))./(cs_temp_enh(:,1)+cs_temp_enh(:,2));
diff_rho = 2*(cs_rho_enh(:,1)-cs_rho_enh(:,2))./(cs_rho_enh(:,1)+cs_rho_enh(:,2));
scatter(diff_temp(abs(diff_temp) < 0.5 & abs(diff_rho) < 0.5 & diff_S > 0),...
    diff_rho(abs(diff_temp) < 0.5 & abs(diff_rho) < 0.5 & diff_S > 0),100,'r.')
scatter(diff_temp(abs(diff_temp) < 0.5 & abs(diff_rho) < 0.5 & diff_S < 0),...
    diff_rho(abs(diff_temp) < 0.5 & abs(diff_rho) < 0.5 & diff_S < 0),100,'b.')
xlabel('(\Delta T)/T_0')
ylabel('(\Delta \rho)/\rho_0')
plot([0 0],[-0.5 0.5],'k')
plot([-0.5 0.5],[0 0],'k')
legend('increase','decrease')
title('\Delta T/n^{\gamma-1} all current sheets')

figure; hold on;
cooond3 = ((diff_temp > 0 & diff_rho > 0) | (diff_temp < 0 & diff_rho < 0))  & abs(cs_n1_nt)>0.7 & cs_eigrats > 5;
histogram(2*cs_Bpar(cooond3)./(cs_jumps_tot(cooond3)),'BinWidth',0.2)
title('B$_n/(0.5 \Delta \mathbf{B)}$ correlated $\Delta \rho$ and $\Delta T$','interpreter','latex')

figure; hold on;
cooond = ord_vajump(:,3)>0 & ord_vjump(:,3) > 0;
scatter(ord_vajump(cooond,3),ord_vjump(cooond,3),100,'r.')
%scatter(ord_vajump(cooond2,2),ord_vjump(cooond2,2),100,'g.')
%scatter(ord_vajump(cooond2,1),ord_vjump(cooond2,1),'b.')
scatter(ord_vajump_tot(cooond,3),ord_vjump_tot(cooond,3),'k.')
slopex2 = ord_vajump(cooond,3)\ord_vjump(cooond,3);
fz = slopex2.*ord_vjump(cooond,3);
[gz ,~] = rsquare(ord_vajump(cooond,3),fz);
%plot([0 80],[0 slopex2*80],'b')
%plot([0 70],[0 70],'k')
plot([0 100],[0 100],'k')
plot([0 100],[0 slopex2*100],'r')
daspect([1 1 1])
annotation('textbox',[.3 .5 .4 .4],'EdgeColor','w','String',['y =',num2str(slopex2(1),2),'x'],'FitBoxToText','on','Color','r');
annotation('textbox',[.3 .5 .3 .3],'EdgeColor','w','String','y=x','FitBoxToText','on');
xlabel('\Delta v_A (km/s)')
ylabel('\Delta v (km/s)')
%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
title('Walen Relation')
axis([0 50 0 50])
legend('sub-layer','whole layer','Location','Southeast')

figure; hold on;
diff_S = (10^4*cs_temp_enh(cooond,1)./cs_rho_enh(cooond,1).^(2/3)-10^4*cs_temp_enh(cooond,2)./cs_rho_enh(cooond,2).^(2/3));
diff_temp = 2*(cs_temp_enh(cooond,1)-cs_temp_enh(cooond,2))./(cs_temp_enh(cooond,1)+cs_temp_enh(cooond,2));
diff_rho = 2*(cs_rho_enh(cooond,1)-cs_rho_enh(cooond,2))./(cs_rho_enh(cooond,1)+cs_rho_enh(cooond,2));
scatter(diff_temp(abs(diff_temp) < 0.5 & abs(diff_rho) < 0.5 & diff_S > 0),...
    diff_rho(abs(diff_temp) < 0.5 & abs(diff_rho) < 0.5 & diff_S > 0),100,'r.')
scatter(diff_temp(abs(diff_temp) < 0.5 & abs(diff_rho) < 0.5 & diff_S < 0),...
    diff_rho(abs(diff_temp) < 0.5 & abs(diff_rho) < 0.5 & diff_S < 0),100,'b.')
xlabel('(\Delta T)/T_0')
ylabel('(\Delta \rho)/\rho_0')
plot([0 0],[-0.5 0.5],'k')
plot([-0.5 0.5],[0 0],'k')
legend('increase','decrease')
title('\Delta T/n^{\gamma-1} alfvenic current sheets')

figure; hold on;
cooond2 = ord_vajump(:,3)>0 & ord_vjump(:,3) > 0 & abs(cs_n1_nt)>0.7 & cs_eigrats > 5;
histogram(2*cs_Bpar(cooond2)./(cs_jumps_tot(cooond2)),'BinWidth',0.2)
title('B$_n/(0.5 \Delta \mathbf{B)}$ Alfvenic layers','interpreter','latex')

figure; hold on
histogram(-cs_wal_len(cooond2)/6.67.*cs_vx(cooond2).*...
    sin(cs_norm_vx(cooond2))./ion_length(cooond2))%,16,'BinLimits',[0 80])
xlabel('\lambda_i')
title('good walen relation interval')

% figure; hold on;
% diff_S = (10^4*cs_temp_enh(cooond2,1)./cs_rho_enh(cooond2,1).^(2/3)-10^4*cs_temp_enh(cooond2,2)./cs_rho_enh(cooond2,2).^(2/3));
% diff_temp = 2*(cs_temp_enh(cooond2,1)-cs_temp_enh(cooond2,2))./(cs_temp_enh(cooond2,1)+cs_temp_enh(cooond2,2));
% diff_rho = 2*(cs_rho_enh(cooond2,1)-cs_rho_enh(cooond2,2))./(cs_rho_enh(cooond2,1)+cs_rho_enh(cooond2,2));
% bn = 2*cs_Bpar(cooond2)./cs_jumps_tot(cooond2);
% scatter(diff_temp(bn <= 0.4 & diff_S > 0),diff_rho(bn <= 0.4 & diff_S > 0),'r.')
% scatter(diff_temp(bn <= 0.4 & diff_S < 0),diff_rho(bn <= 0.4 & diff_S < 0),'b.')
% xlabel('(\Delta T)/T_0')
% ylabel('(\Delta \rho)/\rho_0')
% plot([0 0],[-0.5 0.5],'k')
% plot([-0.5 0.5],[0 0],'k')
% title('\Delta T/n^{\gamma-1} alfvenic layers B_n < 0.4')

