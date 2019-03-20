load('current_sheets.mat')
dist_to_xline = zeros(21,1);
ex = 0;
for jjj = [10,13,20,25,32,39,42,45,49,129,132,135,178,202,216,219,220,235,270,272,284]
    ex = ex+1
    BB = current_sheets(jjj).mag_data;
    vv = current_sheets(jjj).v_data;
    rho = current_sheets(jjj).rho_data;
    temp = current_sheets(jjj).temp_data;
    ti = current_sheets(jjj).ti;
    te = current_sheets(jjj).te;
    jt = current_sheets(jjj).jt;

    bx_int1 = BB(:,1,1); by_int1 = BB(:,1,2); bz_int1 = BB(:,1,3);
    vx = vv(:,1,1); vy = vv(:,1,2); vz = vv(:,1,3);
    ti = ti(:,1);
    te = te(:,1);

    %figure('Visible','off','rend','painters','pos',[10 10 1000 800])
    figure('rend','painters','pos',[10 10 1000 800])
    subplot(4,1,1); hold on;
    %plot_all_mms_B(t,bx_int,by_int,bz_int,sc)
    plot(jt,bx_int1); plot(jt,by_int1); plot(jt,bz_int1);
    lg = legend('B_x','B_y','B_z','Location','northeastoutside');
    ylabel('nT')
    subplot(4,1,2); hold on
    plot(ti(ti~=0),vx(ti ~= 0)-nanmean(vx(ti ~= 0)),'r')
    plot(ti(ti~=0),vy(ti ~= 0)-nanmean(vy(ti ~= 0)),'g','HandleVisibility','off')
    plot(ti(ti~=0),vz(ti ~= 0)-nanmean(vz(ti ~= 0)),'b','HandleVisibility','off')
    lg = legend('v_x','v_y','v_z','Location','northeastoutside');
    ylabel('km/s')
    subplot(4,1,3); hold on
    plot(te(te~=0),rho(te~=0,1),'k')
    legend('\rho^-','Location','northeastoutside')
    ylabel('m^{-3}')
    subplot(4,1,4)
    plot(ti(ti~=0),temp(ti~=0,1),'k')
    legend('T','Location','northeastoutside')
    ylabel('eV')
    x = ginput(2);
    t1 = find(abs(jt-x(1,1)) == min(abs(x(1,1)-jt)));
    t2 = find(abs(jt-x(2,1)) == min(abs(x(2,1)-jt)));

    %distance to x-line is 10 times outflow length for petschek
    dist_to_xline(ex) = 10*mean(vx)*(t1 - t2)/6.67;
    close all
end

figure
histogram(dist_to_xline/150000000,'BinWidth',0.01)