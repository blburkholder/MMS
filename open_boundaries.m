load('walen_analysis/cs_vx.mat')
load('walen_analysis/cs_vy.mat')
load('walen_analysis/cs_vz.mat')
load('walen_analysis/cs_jumps_tot')
load('walen_analysis/v_jumps_totx')
load('walen_analysis/va_jumps_totx')
load('walen_analysis/v_jumps_toty')
load('walen_analysis/va_jumps_toty')
load('walen_analysis/v_jumps_totz')
load('walen_analysis/va_jumps_totz')
load('walen_analysis/v_jumps2')
load('walen_analysis/va_jumps2')
load('walen_analysis/tags')
load('walen_analysis/cs_loc')
load('walen_analysis/cs_Bpar')    
load('walen_analysis/cs_Bperp')    
load('walen_analysis/cs_eigrats')    
load('walen_analysis/cs_n1_nt')   
load('walen_analysis/cs_rho_enh')
load('walen_analysis/cs_temp_enh') 
load('walen_analysis/cs_p_enh') 
load('walen_analysis/cs_norm_vx')
load('walen_analysis/cs_wal_len')
load('walen_analysis/press_balls_std')
load('walen_analysis/cs_Bmag')
load('walen_analysis/cs_shear')
load('walen_analysis/cs_len')

mu = 4*pi*10^(-7);
cs = length(cs_shear);

v_jumps_tot = sqrt(v_jumps_totx.^2+v_jumps_toty.^2+v_jumps_totz.^2);
cs_v = sqrt(cs_vx.^2+cs_vy.^2+cs_vz.^2);

press_balls_left = (1e-18)*cs_Bmag(:,1)/(2*mu) + (1e-9)*cs_p_enh(:,1);
press_balls_right = (1e-18)*cs_Bmag(:,2)/(2*mu) + (1e-9)*cs_p_enh(:,2);
pressure_balance = press_balls_left./press_balls_right;

lpu = (1e-18)*cs_Bmag(:,1)/(2*mu) + (1e-9)*(cs_p_enh(:,1)+press_balls_std(:,1));
lpd = (1e-18)*cs_Bmag(:,1)/(2*mu) + (1e-9)*(cs_p_enh(:,1)-press_balls_std(:,1));
rpu = (1e-18)*cs_Bmag(:,2)/(2*mu) + (1e-9)*(cs_p_enh(:,2)+press_balls_std(:,2));
rpd = (1e-18)*cs_Bmag(:,2)/(2*mu) + (1e-9)*(cs_p_enh(:,2)-press_balls_std(:,2));

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

pleft = cs_p_enh(:,1);
pright = cs_p_enh(:,2);
bstronger = cs_Bmag(:,1) > cs_Bmag(:,2);
pweaker = pleft < pright;
fixable = pweaker == bstronger  & ~(pressure_balance > 0.95 & pressure_balance < 1.05);

fix_fact = (1/(2*mu*(1e-9)))*((1e-18)*cs_Bmag(fixable,2) - (1e-18)*cs_Bmag(fixable,1))./...
     (cs_p_enh(fixable,1)-cs_p_enh(fixable,2))-1;

press_balls_left(fixable) = (1e-18)*cs_Bmag(fixable,1)/(2*mu) + (1e-9)*cs_p_enh(fixable,1).*(1+fix_fact);
press_balls_right(fixable) = (1e-18)*cs_Bmag(fixable,2)/(2*mu) + (1e-9)*cs_p_enh(fixable,2).*(1+fix_fact);
fixed_pressure_balance = press_balls_left./press_balls_right;

figure
scatter(press_balls_left(cs_Bmag(:,1) > cs_Bmag(:,2)),press_balls_right(cs_Bmag(:,1) > cs_Bmag(:,2)),'k.')
scatter(press_balls_right(cs_Bmag(:,1) < cs_Bmag(:,2)),press_balls_left(cs_Bmag(:,1) < cs_Bmag(:,2)),'k.')
daspect([1 1 1])
xlabel('(p+B^2/2\mu)_{stronger}')
ylabel('(p+B^2/2\mu)_{weaker}')
title('fixed pressure balance')

cs_p_enh(fixable,1) = cs_p_enh(fixable,1).*(1+fix_fact);
cs_p_enh(fixable,2) = cs_p_enh(fixable,2).*(1+fix_fact);

 fixed_cs_beta(:,1) = (1e-9)*2*mu*cs_p_enh(:,1)./((1e-18)*cs_Bmag(:,1));
 fixed_cs_beta(:,2) = (1e-9)*2*mu*cs_p_enh(:,2)./((1e-18)*cs_Bmag(:,2));
 beta_fmean = mean(fixed_cs_beta,2);

figure
scatter(max(cs_Bmag'),min(cs_Bmag'),'k.')
daspect([1 1 1])
xlabel('B^2_{stronger} (nT^2)')
ylabel('B^2_{weaker} (nT^2)')
title('B^2')

orang = [246/255 103/255 51/255];
purp = [82/255 45/255 128/255];

ion_length = (3e8*1./sqrt(mean(cs_rho_enh,2)*(1.6e-19)^2/(1.67e-27*8.9e-12)))/1000;

figure
histogram(180*acos(cs_shear)/pi)
title('magnetic shear')
xlabel('$\theta$','interpreter','latex')

figure
histogram(beta_fmean(beta_fmean < 30),'BinWidth',2)
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
diff_B = 2*(sqrt(cs_Bmag(:,1))-sqrt(cs_Bmag(:,2)))./(sqrt(cs_Bmag(:,1))+sqrt(cs_Bmag(:,2)));
outflow_cond = ((diff_rho > 0 & diff_B < 0) | (diff_rho < 0 & diff_B > 0));
fpbal_cond = fixed_pressure_balance < 1.05 & fixed_pressure_balance > 0.95;
coooond5 = pressure_balance < 1.05 & pressure_balance > 0.95;
scatter(diff_B(~outflow_cond & fpbal_cond & diff_rho > 0),diff_rho(~outflow_cond & fpbal_cond & diff_rho > 0),100,'r.')
scatter(diff_B(outflow_cond & fpbal_cond & diff_rho > 0),diff_rho(outflow_cond & fpbal_cond & diff_rho > 0),100,'g.')
scatter(diff_B(outflow_cond & coooond5 & diff_rho > 0),diff_rho(outflow_cond & coooond5 & diff_rho > 0),100,'b.')
scatter(-diff_B(~outflow_cond & fpbal_cond & diff_rho < 0),-diff_rho(~outflow_cond & fpbal_cond & diff_rho < 0),100,'r.')
scatter(-diff_B(outflow_cond & fpbal_cond &diff_rho < 0),-diff_rho(outflow_cond & fpbal_cond & diff_rho < 0),100,'g.')
scatter(-diff_B(outflow_cond & coooond5 &diff_rho < 0),-diff_rho(outflow_cond & coooond5 & diff_rho < 0),100,'b.')
xlabel('\Delta B/B_0')
ylabel('(n_{larger} - n_{smaller})/n_0')
axis([-1 1 0 0.25])
title('changes in average quantities')
legend('inconsistent','p correction - consistent','consistent')

% figure; hold on;
% scatter(diff_B(fpbal_cond & diff_rho > 0 & walen_cond),diff_rho(fpbal_cond & diff_rho > 0 & walen_cond),[],...
%     log10(cs_Bpar(fpbal_cond & diff_rho > 0 & walen_cond)./blarger(fpbal_cond & diff_rho > 0 & walen_cond)),'o','fillled')
% scatter(diff_B(fpbal_cond & diff_rho > 0 & ~walen_cond),diff_rho(fpbal_cond & diff_rho > 0 & ~walen_cond),[],...
%     log10(cs_Bpar(fpbal_cond & diff_rho > 0 & ~walen_cond)./blarger(fpbal_cond & diff_rho > 0 & ~walen_cond)),'p','fillled')
% scatter(-diff_B(fpbal_cond & diff_rho < 0 & walen_cond),-diff_rho(fpbal_cond & diff_rho < 0 & walen_cond),[],...
%     log10(cs_Bpar(fpbal_cond & diff_rho < 0 & walen_cond)./blarger(fpbal_cond & diff_rho < 0 & walen_cond)),'o','filled')
% scatter(-diff_B(fpbal_cond & diff_rho < 0 & ~walen_cond),-diff_rho(fpbal_cond & diff_rho < 0 & ~walen_cond),[],...
%     log10(cs_Bpar(fpbal_cond & diff_rho < 0 & ~walen_cond)./blarger(fpbal_cond & diff_rho < 0 & ~walen_cond)),'p','filled')
% axis([-1 1 0 0.25])


figure; hold on;
walen_cond = ord_vajump(:,3)>0 & ord_vjump(:,3) > 0;
sum(walen_cond)
'number of walen relationers'
scatter(ord_vajump(walen_cond,3),ord_vjump(walen_cond,3),100,'b.')
scatter(ord_vajump_tot(walen_cond,3),ord_vjump_tot(walen_cond,3),100,'r.')
slopex2 = ord_vajump(walen_cond,3)\ord_vjump(walen_cond,3);
fz = slopex2.*ord_vjump(walen_cond,3);
[gz ,~] = rsquare(ord_vajump(walen_cond,3),fz);
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
dir_cond = ((abs(cs_n1_nt)>0.7 & cs_eigrats > 5) | cs_eigrats > 10);
'number walen relationers and good boundary normal'
sum(walen_cond& dir_cond)
blarger = max(cs_Bmag,[],2);
histogram(cs_Bpar(walen_cond& dir_cond)./(blarger(walen_cond& dir_cond)),'BinWidth',0.03)
title('B$_n/|\mathbf{B}|$ Alfvenic layers','interpreter','latex')

figure; hold on;
histogram(cs_Bpar(walen_cond & dir_cond & outflow_cond & fpbal_cond)./(blarger(walen_cond& dir_cond & outflow_cond & fpbal_cond)),'BinWidth',0.03)
title('B$_n/| \mathbf{B}|$ All reconnection signatures','interpreter','latex')


figure; hold on
histogram(-cs_wal_len(walen_cond& dir_cond)/6.67.*cs_vx(walen_cond& dir_cond).*...
    sin(cs_norm_vx(walen_cond& dir_cond))./ion_length(walen_cond& dir_cond))%,16,'BinLimits',[0 80])
xlabel('\lambda_i')
title('good walen relation interval')

figure; hold on
histogram(-cs_wal_len(walen_cond& dir_cond)/6.67.*cs_vx(walen_cond& dir_cond).*...
    sin(cs_norm_vx(walen_cond& dir_cond))./(sqrt(beta_fmean(walen_cond& dir_cond)).*ion_length(walen_cond& dir_cond)))%,16,'BinLimits',[0 80])
xlabel('r_{gi}')
title('good walen relation interval')

figure; hold on
diff_S = 2*((1e-9)*cs_p_enh(:,1)./cs_rho_enh(:,1).^(5/3)-(1e-9)*cs_p_enh(:,2)./cs_rho_enh(:,2).^(5/3))./...
      ((1e-9)*cs_p_enh(:,1)./cs_rho_enh(:,1).^(5/3)+(1e-9)*cs_p_enh(:,2)./cs_rho_enh(:,2).^(5/3));
scatter(beta_fmean(fpbal_cond & beta_fmean<30),abs(diff_S(fpbal_cond & beta_fmean<30)),[],'b.')
title('Entropy Changes')
beta_range = 2:2:20;
smean = zeros(length(beta_range),1);
stmean =  zeros(length(beta_range),1);
for i = 1:length(beta_range)
    cond = beta_fmean > beta_range(1)*(i-1) & beta_fmean < beta_range(1)*i;
    if sum(cond) > 1
        smean(i) = mean(abs(diff_S(cond & fpbal_cond & walen_cond & outflow_cond)));
        stmean(i) = std(abs(diff_S(cond & fpbal_cond & walen_cond & outflow_cond)));
    else
        smean(i) = smean(i-1);
        stmean(i) = stmean(i-1);
    end
end 
gloop = find(isnan(smean));
smean(gloop) = (smean(gloop-1)+smean(gloop+1))/2;
f = fit(beta_range',smean,'exp1');
errorbar(beta_range,f.a*exp(beta_range*f.b),stmean,'r')
ylabel('\Delta S/S_0')
xlabel('\beta')
ylim([0 0.4])
xlim([0 20])

figure; hold on
scatter(beta_fmean(fpbal_cond & beta_fmean<30),...
    abs(diff_rho(fpbal_cond & beta_fmean<30)),[],'b.')
title('Density Changes')
for i = 1:length(beta_range)
    cond = beta_fmean > beta_range(1)*(i-1) & beta_fmean < beta_range(1)*i;
    if sum(cond) > 1
        smean(i) = mean(abs(diff_rho(cond & fpbal_cond)));
        stmean(i) = std(abs(diff_rho(cond & fpbal_cond)));
    else
        smean(i) = smean(i-1);
        stmean(i) = stmean(i-1);
    end
end
smean(gloop) = (smean(gloop-1)+smean(gloop+1))/2;
f = fit(beta_range',smean,'exp1');
errorbar(beta_range,f.a*exp(beta_range*f.b),stmean,'r')
ylabel('\Delta \rho/\rho_0')
xlabel('\beta')
ylim([0 0.4])
xlim([0 20])

figure; hold on
scatter(cs_Bpar(dir_cond & walen_cond)./blarger(dir_cond & walen_cond),...
    abs(diff_rho(dir_cond & walen_cond)),'b.')
scatter(cs_Bpar(dir_cond & ~walen_cond)./blarger(dir_cond & ~walen_cond),...
    abs(diff_rho(dir_cond & ~walen_cond)),'r.')
bn_range = 0.01:.01:0.1;
smean = zeros(length(bn_range),1);
ssmean = zeros(length(bn_range),1);
stmean =  zeros(length(bn_range),1);
sstmean = zeros(length(bn_range),1);
for i = 1:length(bn_range)
    cond = cs_Bpar./blarger > bn_range(1)*(i-1) & cs_Bpar./blarger < bn_range(1)*i;
    if sum(cond & walen_cond) > 1
        smean(i) = mean(abs(diff_rho(cond & walen_cond)));
        stmean(i) = std(abs(diff_rho(cond & walen_cond)));
    else
        smean(i) = smean(i-1);
        stmean(i) = stmean(i-1);
    end
    if sum(cond & ~walen_cond) > 1
        ssmean(i) = mean(abs(diff_rho(cond & ~walen_cond)));
        sstmean(i) = std(abs(diff_rho(cond & ~walen_cond)));
    else
        ssmean(i) = ssmean(i-1);
        sstmean(i) = sstmean(i-1);
    end
end
f = fit(bn_range',smean,'exp1');
plot([0,bn_range],f.a*exp([0,bn_range]*f.b),'b')
errorbar(bn_range,f.a*exp(bn_range*f.b),stmean.^1.5,'b')
ff = fit(bn_range',ssmean,'exp1');
plot([0,bn_range],ff.a*exp([0,bn_range]*ff.b),'r')
errorbar(bn_range,ff.a*exp(bn_range*ff.b),sstmean.^1.5,'r')
xlabel('B_n/|B|')
ylabel('\Delta \rho/\rho_0')
legend('good Walen relation','no Walen relation')
axis([0 0.11 0 0.4])

figure; hold on
scatter(cs_Bpar(dir_cond & walen_cond & fpbal_cond)./blarger(dir_cond & walen_cond & fpbal_cond),...
    abs(diff_S(dir_cond & walen_cond & fpbal_cond)),'b.')
scatter(cs_Bpar(dir_cond & ~walen_cond & fpbal_cond)./blarger(dir_cond & ~walen_cond & fpbal_cond),...
    abs(diff_S(dir_cond & ~walen_cond & fpbal_cond)),'r.')
bn_range = 0.01:.01:0.1;
smean = zeros(length(bn_range),1);
ssmean = zeros(length(bn_range),1);
stmean =  zeros(length(bn_range),1);
sstmean = zeros(length(bn_range),1);
for i = 1:length(bn_range)
    cond = cs_Bpar./blarger > bn_range(1)*(i-1) & cs_Bpar./blarger < bn_range(1)*i & fpbal_cond;
    if sum(cond & walen_cond) > 1
        smean(i) = mean(abs(diff_S(cond & walen_cond)));
        stmean(i) = std(abs(diff_S(cond & walen_cond)));
    else
        smean(i) = smean(i-1);
        stmean(i) = stmean(i-1);
    end
    if sum(cond & ~walen_cond) > 1
        ssmean(i) = mean(abs(diff_S(cond & ~walen_cond)));
        sstmean(i) = std(abs(diff_S(cond & ~walen_cond)));
    else
        ssmean(i) = ssmean(i-1);
        sstmean(i) = sstmean(i-1);
    end
end
f = fit(bn_range',smean,'exp1');
plot([0,bn_range],f.a*exp([0,bn_range]*f.b),'b')
errorbar(bn_range,f.a*exp(bn_range*f.b),stmean.^2,'b')
ff = fit(bn_range',ssmean,'exp1');
plot([0,bn_range],ff.a*exp([0,bn_range]*ff.b),'r')
errorbar(bn_range,ff.a*exp(bn_range*ff.b),sstmean.^2,'r')
xlabel('B_n/|B|')
ylabel('\Delta S/S_0')
legend('good Walen relation','no Walen relation')
axis([0 0.11 0 0.4])

figure; hold on
scatter(cs_Bpar(dir_cond & walen_cond)./blarger(dir_cond & walen_cond),...
    abs(diff_B(dir_cond & walen_cond)),'b.')
scatter(cs_Bpar(dir_cond & ~walen_cond)./blarger(dir_cond & ~walen_cond),...
    abs(diff_B(dir_cond & ~walen_cond)),'r.')
bn_range = 0.01:.01:0.1;
smean = zeros(length(bn_range),1);
ssmean = zeros(length(bn_range),1);
stmean =  zeros(length(bn_range),1);
sstmean = zeros(length(bn_range),1);
% for i = 1:length(bn_range)
%     cond = cs_Bpar./blarger > bn_range(1)*(i-1) & cs_Bpar./blarger < bn_range(1)*i;
%     if sum(cond & walen_cond) > 1
%         smean(i) = mean(abs(diff_B(cond & walen_cond)));
%         stmean(i) = std(abs(diff_B(cond & walen_cond)));
%     else
%         smean(i) = smean(i-1);
%         stmean(i) = stmean(i-1);
%     end
%     if sum(cond & ~walen_cond) > 1
%         ssmean(i) = mean(abs(diff_B(cond & ~walen_cond)));
%         sstmean(i) = std(abs(diff_B(cond & ~walen_cond)));
%     else
%         ssmean(i) = ssmean(i-1);
%         sstmean(i) = sstmean(i-1);
%     end
% end
% f = fit(bn_range',smean,'exp1');
% plot([0,bn_range],f.a*exp([0,bn_range]*f.b),'b')
% errorbar(bn_range,f.a*exp(bn_range*f.b),stmean.^2,'b')
% ff = fit(bn_range',ssmean,'exp1');
% plot([0,bn_range],ff.a*exp([0,bn_range]*ff.b),'r')
% errorbar(bn_range,ff.a*exp(bn_range*ff.b),sstmean.^2,'r')
xlabel('B_n/|B|')
ylabel('\Delta B/B_0')
legend('good Walen relation','no Walen relation')
plot([0 0.11],[0.1 0.1],'k')
plot([0.02 0.02],[0.1 0.6],'k')
axis([0 0.11 0 0.6])

