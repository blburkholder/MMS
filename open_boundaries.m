% load('walen_analysis/cs_vx.mat')
% load('walen_analysis/cs_vy.mat')
% load('walen_analysis/cs_vz.mat')
% load('walen_analysis/cs_jumps_tot')
% load('walen_analysis/v_jumps_totx')
% load('walen_analysis/va_jumps_totx')
% load('walen_analysis/v_jumps_toty')
% load('walen_analysis/va_jumps_toty')
% load('walen_analysis/v_jumps_totz')
% load('walen_analysis/va_jumps_totz')
% load('walen_analysis/v_jumps2')
% load('walen_analysis/va_jumps2')
% load('walen_analysis/tags')
% load('walen_analysis/cs_loc')
% load('walen_analysis/cs_Bpar')    
% load('walen_analysis/cs_Bperp')    
% load('walen_analysis/cs_eigrats')    
% load('walen_analysis/cs_n1_nt')   
% load('walen_analysis/cs_rho_enh')
% %load('walen_analysis/cs_temp_enh') 
% load('walen_analysis/cs_p_enh') 
% load('walen_analysis/cs_norm_vx')
% load('walen_analysis/cs_wal_len')
% load('walen_analysis/cs_beta')
% load('walen_analysis/press_balls_std')
% load('walen_analysis/cs_Bmag')
% load('walen_analysos/cs_shear')

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
load('v_jumps2')
load('va_jumps2')
load('tags')
load('cs_loc')
load('cs_Bpar')    
load('cs_Bperp')    
load('cs_eigrats')    
load('cs_n1_nt')   
load('cs_rho_enh')
load('cs_p_enh') 
load('cs_norm_vx')
load('cs_wal_len')
load('cs_beta')
load('press_balls_std')
load('cs_Bmag')
load('cs_shear')

mu = 4*pi*10^(-7);
cs = 381;

v_jumps_tot = sqrt(v_jumps_totx.^2+v_jumps_toty.^2+v_jumps_totz.^2);
cs_v = sqrt(cs_vx.^2+cs_vy.^2+cs_vz.^2);

press_balls_left = (1e-18)*cs_Bmag(:,1)/(2*mu) + (1.602e-19)*cs_rho_enh(:,1).*cs_temp_enh(:,1);
press_balls_right = (1e-18)*cs_Bmag(:,2)/(2*mu) + (1.602e-19)*cs_rho_enh(:,2).*cs_temp_enh(:,2);
%press_balls_left = (1e-18)*cs_Bmag(:,1)/(2*mu) + (1e-9)*cs_p_enh(:,1);
%press_balls_right = (1e-18)*cs_Bmag(:,2)/(2*mu) + (1e-9)*cs_p_enh(:,2);

lpu = (1e-18)*cs_Bmag(:,1)/(2*mu) + (1.602e-19)*cs_rho_enh(:,1).*(cs_temp_enh(:,1)+press_balls_std(:,1));
lpd = (1e-18)*cs_Bmag(:,1)/(2*mu) + (1.602e-19)*cs_rho_enh(:,1).*(cs_temp_enh(:,1)-press_balls_std(:,1));
rpu = (1e-18)*cs_Bmag(:,2)/(2*mu) + (1.602e-19)*cs_rho_enh(:,2).*(cs_temp_enh(:,2)+press_balls_std(:,2));
rpd = (1e-18)*cs_Bmag(:,2)/(2*mu) + (1.602e-19)*cs_rho_enh(:,2).*(cs_temp_enh(:,2)-press_balls_std(:,2));
% lpu = (1e-18)*cs_Bmag(:,1)/(2*mu) + (1e-9)*(cs_p_enh(:,1)+press_balls_std(:,1));
% lpd = (1e-18)*cs_Bmag(:,1)/(2*mu) + (1e-9)*(cs_p_enh(:,1)-press_balls_std(:,1));
% rpu = (1e-18)*cs_Bmag(:,2)/(2*mu) + (1e-9)*(cs_p_enh(:,2)+press_balls_std(:,2));
% rpd = (1e-18)*cs_Bmag(:,2)/(2*mu) + (1e-9)*(cs_p_enh(:,2)-press_balls_std(:,2));

press_balls_left_errup = abs(press_balls_left - lpu);
press_balls_left_errdown = abs(press_balls_left - lpd);
press_balls_left_err = max([press_balls_left_errup,press_balls_left_errdown]');

press_balls_right_errup = abs(press_balls_right-rpu);
press_balls_right_errdown = abs(press_balls_right - rpd);
press_balls_right_err = max([press_balls_right_errup,press_balls_right_errdown]');

figure; hold on
scatter(press_balls_left,press_balls_right,'k.')
errorbar(press_balls_left,press_balls_right,...
    press_balls_right_err,'.')
herrorbar(press_balls_left,press_balls_right,...
    press_balls_left_err,'.')
daspect([1 1 1])
plot([0 1.4e-9],[0 1.05*1.4e-9],'k')
plot([0 1.4e-9],[0 0.95*1.4e-9],'k')
xlabel('(p+B^2/2\mu)_{1}')
ylabel('(p+B^2/2\mu)_{2}')
title('pressure balance')

consistent = zeros(cs,1);
for i = 1:cs
    if press_balls_left(i)/press_balls_right(i) < 1.05 && press_balls_left(i)/press_balls_right(i) > 0.95
        consistent(i) = 1;
    else
        %see if box defined by errorbars intersects y = 0.95x or y = 1.05x
        top = [[press_balls_left(i) - press_balls_left_err(i) press_balls_left(i) + press_balls_left_err(i)],...
            [press_balls_right(i) + press_balls_right_err(i) press_balls_right(i) + press_balls_right_err(i)]];
        bottom = [[press_balls_left(i) - press_balls_left_err(i) press_balls_left(i) + press_balls_left_err(i)],...
            [press_balls_right(i) - press_balls_right_err(i) press_balls_right(i) - press_balls_right_err(i)]];
        left = [[press_balls_left(i) - press_balls_left_err(i) press_balls_left(i) - press_balls_left_err(i)],...
            [press_balls_right(i) - press_balls_right_err(i) press_balls_right(i) + press_balls_right_err(i)]];
        right = [[press_balls_left(i) + press_balls_left_err(i) press_balls_left(i) + press_balls_left_err(i)],...
            [press_balls_right(i) - press_balls_right_err(i) press_balls_right(i) + press_balls_right_err(i)]];

        [x1,y1] = polyxpoly([top(1),top(2)],[top(3),top(4)],[0 1.4e-9],[0 1.05*1.4e-9]); 
        [x2,y2] = polyxpoly([top(1),top(2)],[top(3),top(4)],[0  1.4e-9],[0 0.95*1.4e-9]);
        if isempty(x1) x1 = 0; end;         if isempty(x2) x2 = 0; end

        [x3,y3] = polyxpoly([bottom(1),bottom(2)],[bottom(3),bottom(4)],[0 1.4e-9],[0 1.05*1.4e-9]); 
        [x4,y4] = polyxpoly([bottom(1),bottom(2)],[bottom(3),bottom(4)],[0  1.4e-9],[0 0.95*1.4e-9]);
        if isempty(x3) x3 = 0; end;         if isempty(x4) x4 = 0; end

        [x5,y5] = polyxpoly([left(1),left(2)],[left(3),left(4)],[0 1.4e-9],[0 1.05*1.4e-9]); 
        [x6,y6] = polyxpoly([left(1),left(2)],[left(3),left(4)],[0  1.4e-9],[0 0.95*1.4e-9]);
        if isempty(x5) x5 = 0; end;         if isempty(x6) x6 = 0; end

        [x7,yi7] = polyxpoly([right(1),right(2)],[right(3),right(4)],[0 1.4e-9],[0 1.05*1.4e-9]); 
        [x8,y8] = polyxpoly([right(1),right(2)],[right(3),right(4)],[0  1.4e-9],[0 0.95*1.4e-9]);
        if isempty(x7) x7 = 0; end;         if isempty(x8) x8 = 0; end

        if sum(x1+x2+x3+x4+x5+x6+x7+x8) > 0
            consistent(i) = 1;
        end
    end       
end

% 
% figure
% scatter(cs_v,press_balls_left./press_balls_right,'k.');

pleft = cs_temp_enh(:,1).*cs_rho_enh(:,1);
pright = cs_temp_enh(:,2).*cs_rho_enh(:,2);
bstronger = cs_Bmag(:,1) > cs_Bmag(:,2);

pweaker = pleft < pright;
fixable = pweaker == bstronger;

% set(gca,'XScale','log')
% 
fix_fact = (1/(2*mu*(1.602e-19)))*((1e-18)*cs_Bmag(fixable,2) - (1e-18)*cs_Bmag(fixable,1))./...
    (cs_rho_enh(fixable,1).*cs_temp_enh(fixable,1)-cs_rho_enh(fixable,2).*cs_temp_enh(fixable,2))-1;
% 
% figure
% scatter(fix_fact,cs_rho_enh(fixable,3).*cs_temp_enh(fixable,3),'k.')
% figure
% scatter(fix_fact,cs_rho_enh(fixable,3),'k.')
% figure
% scatter(fix_fact,cs_vx(fixable),'k.')
% figure
% scatter(fix_fact,(cs_Bmag(fixable,1)+cs_Bmag(fixable,2))/2,'k.')

% press_balls_left(fixable) = (1e-18)*cs_Bmag(fixable,1)/(2*mu) + (1.602e-19)*cs_rho_enh(fixable,1).*(cs_temp_enh(fixable,1).*(1+fix_fact));
% press_balls_right(fixable) = (1e-18)*cs_Bmag(fixable,2)/(2*mu) + (1.602e-19)*cs_rho_enh(fixable,2).*(cs_temp_enh(fixable,2).*(1+fix_fact));
%  figure
%  scatter(press_balls_left(cs_Bmag(:,1) > cs_Bmag(:,2)),press_balls_right(cs_Bmag(:,1) > cs_Bmag(:,2)),'k.')
%  scatter(press_balls_right(cs_Bmag(:,1) < cs_Bmag(:,2)),press_balls_left(cs_Bmag(:,1) < cs_Bmag(:,2)),'k.')
%  daspect([1 1 1])
%  xlabel('(p+B^2/2\mu)_{stronger}')
%  ylabel('(p+B^2/2\mu)_{weaker}')
%  title('fixed pressure balance')
pressure_balance = press_balls_left./press_balls_right;
%cs_temp_enh(fixable,1) = cs_temp_enh(fixable,1).*(1+fix_fact);
%cs_temp_enh(fixable,2) = cs_temp_enh(fixable,2).*(1+fix_fact);

figure
scatter(max(cs_Bmag'),min(cs_Bmag'),'k.')
daspect([1 1 1])
xlabel('B^2_{stronger} (nT^2)')
ylabel('B^2_{weaker} (nT^2)')
title('B^2')

orang = [246/255 103/255 51/255];
purp = [82/255 45/255 128/255];

ion_length = (3e8*1./sqrt(min(cs_rho_enh,[],2)*(1.6e-19)^2/(1.67e-27*8.9e-12)))/1000;

% va_jumps_totx = (10^-12)*va_jumps_totx./sqrt(mu*(1.67*10^(-27)));
% va_jumps_toty = (10^-12)*va_jumps_toty./sqrt(mu*(1.67*10^(-27)));
% va_jumps_totz = (10^-12)*va_jumps_totz./sqrt(mu*(1.67*10^(-27)));

% figure; hold on
% cooond = abs(cs_n1_nt)>0.7 & cs_eigrats > 5;
% histogram((cs_inds(cs_v<450 & cooond,1) - cs_inds(cs_v<450 & cooond,2))/128.*...
%     cs_vx(cs_v<450 & cooond).*sin(cs_norm_vx(cs_v<450 & cooond))./ion_length(cs_v<450 & cooond),16,'BinLimits',[0 80])
% histogram((cs_inds(cs_v>450 & cooond,1) - cs_inds(cs_v>450 & cooond,2))/128.*...
%     cs_vx(cs_v>450 & cooond).*sin(cs_norm_vx(cs_v>450 & cooond))./ion_length(cs_v>450 & cooond),16,'BinLimits',[0 80])
% legend('v<450','v>450')
% xlabel('\lambda_i')
% title('current sheet width')

figure
histogram(180*acos(cs_shear)/pi)
title('magnetic shear')
xlabel('$\theta$','interpreter','latex')

bigbeta = cs_beta(:,3);
figure
histogram(bigbeta(bigbeta < 40),'BinWidth',2)
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
diff_rho = 2*(cs_rho_enh(:,1)-cs_rho_enh(:,2))./(cs_rho_enh(:,1)+cs_rho_enh(:,2));
diff_p = 2*(cs_temp_enh(:,1).*cs_rho_enh(:,1)-cs_temp_enh(:,2).*cs_rho_enh(:,2))./...
    (cs_temp_enh(:,1).*cs_rho_enh(:,1)+cs_temp_enh(:,2).*cs_rho_enh(:,2));
cooond4 = ((diff_rho > 0 & diff_p > 0 & cs_Bmag(:,2) > cs_Bmag(:,1)) | (diff_rho < 0 & diff_p < 0  & cs_Bmag(:,1) > cs_Bmag(:,2)));
coooond4 = pressure_balance < 1.05 & pressure_balance > 0.95;
%coooond4 = ones(size(cooond4));
scatter(diff_p(~cooond4 & coooond4 & diff_rho > 0),diff_rho(~cooond4 & coooond4 & diff_rho > 0),100,'r.')
scatter(diff_p(cooond4 & coooond4 &diff_rho > 0),diff_rho(cooond4 & coooond4 & diff_rho > 0),100,'b.')
scatter(-diff_p(~cooond4 & coooond4 & diff_rho < 0),-diff_rho(~cooond4 & coooond4 & diff_rho < 0),100,'r.')
scatter(-diff_p(cooond4 & coooond4 &diff_rho < 0),-diff_rho(cooond4 & coooond4 & diff_rho < 0),100,'b.')
scatter(diff_p(cooond4 & coooond4 &diff_rho > 0 & abs(diff_p) > 0.05 & abs(diff_rho) > 0.05),diff_rho(cooond4 & coooond4 & diff_rho > 0 & abs(diff_p) > 0.05 & abs(diff_rho) > 0.05),100,'g.')
scatter(-diff_p(cooond4 & coooond4 &diff_rho < 0 & abs(diff_p) > 0.05 & abs(diff_rho) > 0.05),-diff_rho(cooond4 & coooond4 & diff_rho < 0 & abs(diff_p) > 0.05 & abs(diff_rho) > 0.05),100,'g.')
xlabel('\Delta p/p_0')
ylabel('(n_{larger} - n_{smaller})/n_0')
plot([0.05 0.05],[0.0001 1],'k')
plot([-0.6 0.6],[0.05 0.05],'k')
axis([-0.6 0.6 0 1])
set(gca,'YScale','log')
legend('\Delta |B| > 0','\Delta |B| < 0')
title('changes in average quantities')

% figure; hold on;
% diff_rho = 2*(cs_rho_enh(:,1)-cs_rho_enh(:,2))./(cs_rho_enh(:,1)+cs_rho_enh(:,2));
% diff_p = 2*(cs_p_enh(:,1)-cs_p_enh(:,2))./(cs_p_enh(:,1)+cs_p_enh(:,2));
% cooond4 = ((diff_rho > 0 & diff_p > 0 & cs_Bmag(:,2) > cs_Bmag(:,1)) | (diff_rho < 0 & diff_p < 0  & cs_Bmag(:,1) > cs_Bmag(:,2)));
% coooond4 = pressure_balance < 1.05 & pressure_balance > 0.95;
% scatter(diff_p(~cooond4 & coooond4 & diff_rho > 0),diff_rho(~cooond4 & coooond4 & diff_rho > 0),100,'r.')
% scatter(diff_p(cooond4 & coooond4 &diff_rho > 0),diff_rho(cooond4 & coooond4 & diff_rho > 0),100,'b.')
% scatter(-diff_p(~cooond4 & coooond4 & diff_rho < 0),-diff_rho(~cooond4 & coooond4 & diff_rho < 0),100,'r.')
% scatter(-diff_p(cooond4 & coooond4 &diff_rho < 0),-diff_rho(cooond4 & coooond4 & diff_rho < 0),100,'b.')
% xlabel('\Delta p/p_0')
% ylabel('(n_{larger} - n_{smaller})/n_0')
% plot([0 0],[-1 1],'k')
% axis([-0.4 0.4 0 1])
% set(gca,'YScale','log')
% legend('\Delta |B| > 0','\Delta |B| < 0')
% title('changes in average quantities')

figure; hold on;
cooond = ord_vajump(:,3)>0 & ord_vjump(:,3) > 0;
sum(cooond)
'number of walen relationers'
scatter(ord_vajump(cooond,3),ord_vjump(cooond,3),100,'b.')
scatter(ord_vajump_tot(cooond,3),ord_vjump_tot(cooond,3),100,'r.')
slopex2 = ord_vajump(cooond,3)\ord_vjump(cooond,3);
fz = slopex2.*ord_vjump(cooond,3);
[gz ,~] = rsquare(ord_vajump(cooond,3),fz);
plot([0 100],[0 100],'Color','k')
plot([0 100],[0 slopex2*100],'Color','b')
daspect([1 1 1])
annotation('textbox',[.3 .5 .4 .4],'EdgeColor','w','String',['y =',num2str(slopex2(1),2),'x'],'FitBoxToText','on','Color','b');
annotation('textbox',[.3 .5 .3 .3],'EdgeColor','w','String','y=x','FitBoxToText','on','Color','k');
xlabel('\Delta v_A (km/s)')
ylabel('\Delta v (km/s)')
%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
title('Walen Relation')
axis([0 50 0 50])
legend('sub-layer','whole layer','Location','Southeast')

figure; hold on;
cooond2 = ((abs(cs_n1_nt)>0.7 & cs_eigrats > 5) | cs_eigrats > 10);
'number walen relationers and good boundary normal'
sum(cooond & cooond2)
blarger = max(cs_Bmag');
histogram(cs_Bpar(cooond & cooond2)./(blarger(cooond & cooond2)'),'BinWidth',0.03)
title('B$_n/|\mathbf{B}|$ Alfvenic layers','interpreter','latex')

figure; hold on;
histogram(cs_Bpar(cooond & cooond2 & cooond4 & coooond4)./(blarger(cooond & cooond2 & cooond4 & coooond4)'),'BinWidth',0.03)
title('B$_n/| \mathbf{B}|$ All reconnection signatures','interpreter','latex')

figure; hold on
histogram(-cs_wal_len(cooond & cooond2)/6.67.*cs_vx(cooond & cooond2).*...
    sin(cs_norm_vx(cooond & cooond2))./ion_length(cooond & cooond2))%,16,'BinLimits',[0 80])
xlabel('\lambda_i')
title('good walen relation interval')

figure; hold on
histogram(-cs_wal_len(cooond & cooond2)/6.67.*cs_vx(cooond & cooond2).*...
    sin(cs_norm_vx(cooond & cooond2))./(sqrt(cs_beta(cooond & cooond2,3)).*ion_length(cooond & cooond2)))%,16,'BinLimits',[0 80])
xlabel('r_{gi}')
title('good walen relation interval')

figure; hold on
 diff_S = 2*(cs_temp_enh(:,1)./cs_rho_enh(:,1).^(2/3)-cs_temp_enh(:,2)./cs_rho_enh(:,2).^(2/3))./...
     (cs_temp_enh(:,1)./cs_rho_enh(:,1).^(2/3)+cs_temp_enh(:,2)./cs_rho_enh(:,2).^(2/3));
%  diff_S = 2*((1e-9)*cs_p_enh(:,1)./cs_rho_enh(:,1).^(5/3)-(1e-9)*cs_p_enh(:,2)./cs_rho_enh(:,2).^(5/3))./...
%      ((1e-9)*cs_p_enh(:,1)./cs_rho_enh(:,1).^(5/3)+(1e-9)*cs_p_enh(:,2)./cs_rho_enh(:,2).^(5/3));
coooond4 = pressure_balance < 1.05 & pressure_balance > 0.95;
scatter(cs_beta(coooond4 & cs_beta(:,3)<30,3),abs(diff_S(coooond4 & cs_beta(:,3)<30)),[],'b.')
%xlim([0 30])
%ylim([0 0.5])
title('Entropy Changes')

% for i = 1:12
%     cond = cs_beta(:,3) > (i-1)*2.5 & cs_beta(:,3) < i*2.5;
%     smean(i) = std(abs(diff_S(cond & coooond4)));
% end

for i = 1:12
    cond = cs_beta(:,3) > (i-1)*2.5 & cs_beta(:,3) < i*2.5;
    smean(i) = std(abs(diff_S(cond & coooond4)));
end

f = fit([1.25:2.5:28.75]',smean','exp1');
plot(f,1.25:2.5:28.75,smean)
ylabel('\Delta S/S_0')
xlabel('\beta')





