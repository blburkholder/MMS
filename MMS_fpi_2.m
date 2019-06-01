tic
%evaluates walen relation at cadence of ion FPI measurements
mu = 4*pi*10^(-7);
%this file contains the start index and end index for current sheet
%intervals within the mms data files
load('current_sheets.mat')
%there are 381 current sheets within the 303 files
cs_vx = zeros(381,1);
cs_vy = zeros(381,1);
cs_vz = zeros(381,1);
v_jumps_totx = zeros(381,1);
v_jumps_toty = zeros(381,1);
v_jumps_totz = zeros(381,1);
va_jumps_totx = zeros(381,1);
va_jumps_toty = zeros(381,1);
va_jumps_totz = zeros(381,1);
cs_jumps_tot = zeros(381,1);
v_jumps2 = zeros(381,3);
va_jumps2 = zeros(381,3);
cs_loc = zeros(381,3);
tags = zeros(381,1);
cs_eigrats = zeros(381,1);
cs_Bperp = zeros(381,1);
cs_Bpar = zeros(381,1);
cs_Bmag = zeros(381,2);
cs_n1_nt = zeros(381,1);
cs_rho_enh = zeros(381,3);
cs_temp_enh = zeros(381,3);
cs_p_enh = zeros(381,3);
cs_norm_vx = zeros(381,1);
cs_wal_len = zeros(381,1);
press_balls_std_T = zeros(381,2);
press_balls_std = zeros(381,2);
cs_shear = zeros(381,1);
cs_len = zeros(381,1);
cs = 0;

%201 & 225 missing
for jjj = [1:200,202:224,226:303]
%for jjj = 137
    jjj
    BB = current_sheets(jjj).mag_data;
    vv = current_sheets(jjj).v_data;
    pospos = current_sheets(jjj).pos;
    rho_comp = current_sheets(jjj).rho_data;
    temppar_comp = current_sheets(jjj).temp_data;
    etemppar_comp = current_sheets(jjj).etemp_data;
    tempperp_comp = current_sheets(jjj).tempperp_data;
    etempperp_comp = current_sheets(jjj).etempperp_data;
    p_comp = current_sheets(jjj).p_data;
    ep_comp = current_sheets(jjj).ep_data;
    ti = current_sheets(jjj).ti;
    te = current_sheets(jjj).te;
    jt = current_sheets(jjj).jt;
    fpi_sc = current_sheets(jjj).fpi_sc;
    sc = current_sheets(jjj).sc;
    missing = current_sheets(jjj).m;
    current_density = current_sheets(jjj).curdens;

    %interpolating FPI data to same time stamps
     if fpi_sc > 1

        bx_int1 = BB(:,:,1); by_int1 = BB(:,:,2); bz_int1 = BB(:,:,3);
        vx_comp = vv(:,:,1); vy_comp = vv(:,:,2); vz_comp = vv(:,:,3);
        x_int1 = pospos(:,:,1); y_int1 = pospos(:,:,2); z_int1 = pospos(:,:,3);

        [vx_int1,vy_int1,vz_int1,it] =...
            interp_time_v(vx_comp,vy_comp,vz_comp,ti,fpi_sc);
        [rho_int1,~,~,et] =...
            interp_time_v(rho_comp,rho_comp,rho_comp,te,fpi_sc);
         [temppar_int1,tempperp_int1,~,~] =...
             interp_time_v(temppar_comp,tempperp_comp,temppar_comp,ti,fpi_sc);
         [etemppar_int1,etempperp_int1,~,~] =...
             interp_time_v(etemppar_comp,etempperp_comp,etemppar_comp,te,fpi_sc);
        [p_int1,~,~,~] =...
            interp_time_v(p_comp,p_comp,p_comp,ti,fpi_sc);
        [ep_int1,~,~,~] =...
            interp_time_v(ep_comp,ep_comp,ep_comp,te,fpi_sc);

        ti1 = find(it > max(jt(1),et(1)),1);        ti2 = find(it < min(jt(end),et(end)));
        tii2 = ti2(end);    t = it(ti1:tii2);

        bx_int = zeros(length(t),sc);
        by_int = zeros(length(t),sc);
        bz_int = zeros(length(t),sc);
        x_int = zeros(length(t),sc);
        y_int = zeros(length(t),sc);
        z_int = zeros(length(t),sc);
        vx_int = zeros(length(t),fpi_sc);
        vy_int = zeros(length(t),fpi_sc);
        vz_int = zeros(length(t),fpi_sc);
        rho_int = zeros(length(t),fpi_sc);
        temppar_int = zeros(length(t),fpi_sc);
        tempperp_int = zeros(length(t),fpi_sc);
        etemppar_int = zeros(length(t),fpi_sc);
        etempperp_int = zeros(length(t),fpi_sc);
        p_int = zeros(length(t),fpi_sc);
        ep_int = zeros(length(t),fpi_sc);

        %interpolating mag and position data to same time stamps
        for h = 1:sc
            bx_int(:,h) = interp1(jt,smooth(bx_int1(:,h),30),t);
            by_int(:,h) = interp1(jt,smooth(by_int1(:,h),30),t);
            bz_int(:,h) = interp1(jt,smooth(bz_int1(:,h),30),t);
            x_int(:,h) = interp1(jt,x_int1(:,h),t);
            y_int(:,h) = interp1(jt,y_int1(:,h),t);
            z_int(:,h) = interp1(jt,z_int1(:,h),t);
        end
        for h = 1:fpi_sc
            vx_int(:,h) = interp1(it,vx_int1(:,h),t);
            vy_int(:,h) = interp1(it,vy_int1(:,h),t);
            vz_int(:,h) = interp1(it,vz_int1(:,h),t);
            rho_int(:,h) = interp1(et,rho_int1(:,h),t);
            temppar_int(:,h) = interp1(it,temppar_int1(:,h),t);
            tempperp_int(:,h) = interp1(it,tempperp_int1(:,h),t);
            etemppar_int(:,h) = interp1(et,etemppar_int1(:,h),t);
            etempperp_int(:,h) = interp1(et,etempperp_int1(:,h),t);
            p_int(:,h) = interp1(it,p_int1(:,h),t);
            ep_int(:,h) = interp1(et,ep_int1(:,h),t);
        end

        %interpolating measurements from all available spacecraft to
        %center of formation
        if sc == 4 && fpi_sc == 4
            [px,py,pz] = nsimplex_4vec(bx_int,by_int,bz_int,x_int,y_int,z_int);
            [vx,vy,vz] = nsimplex_4vec(vx_int,vy_int,vz_int,x_int,y_int,z_int);
            rho = nsimplex_4(rho_int,x_int,y_int,z_int);
            temppar = nsimplex_4(temppar_int,x_int,y_int,z_int);
            tempperp = nsimplex_4(tempperp_int,x_int,y_int,z_int);
            etemppar = nsimplex_4(etemppar_int,x_int,y_int,z_int);
            etempperp = nsimplex_4(etempperp_int,x_int,y_int,z_int);
            p = nsimplex_4(p_int,x_int,y_int,z_int);
            ep = nsimplex_4(ep_int,x_int,y_int,z_int);


            [vax,vay,vaz] = nsimplex_4vec([bx_int(:,1)./sqrt(rho_int(:,1)),bx_int(:,2)./sqrt(rho_int(:,2)),bx_int(:,3)./sqrt(rho_int(:,3)),bx_int(:,4)./sqrt(rho_int(:,4))],...
                [by_int(:,1)./sqrt(rho_int(:,1)),by_int(:,2)./sqrt(rho_int(:,2)),by_int(:,3)./sqrt(rho_int(:,3)),by_int(:,4)./sqrt(rho_int(:,4))],...
                [bz_int(:,1)./sqrt(rho_int(:,1)),bz_int(:,2)./sqrt(rho_int(:,2)),bz_int(:,3)./sqrt(rho_int(:,3)),bz_int(:,4)./sqrt(rho_int(:,4))],x_int,y_int,z_int);
        elseif sc == 3 && fpi_sc == 3
            [px,py,pz] = nsimplex_3vec(bx_int,by_int,bz_int,x_int,y_int,z_int);  
            [vx,vy,vz] = nsimplex_3vec(vx_int,vy_int,vz_int,x_int,y_int,z_int);
            rho = nsimplex_3(rho_int,x_int,y_int,z_int);
            temppar = nsimplex_3(temppar_int,x_int,y_int,z_int);
            tempperp = nsimplex_3(tempperp_int,x_int,y_int,z_int);
            etemppar = nsimplex_3(etemppar_int,x_int,y_int,z_int);
            etempperp = nsimplex_3(etempperp_int,x_int,y_int,z_int);
            p = nsimplex_3(p_int,x_int,y_int,z_int);
            ep = nsimplex_3(ep_int,x_int,y_int,z_int);

            [vax,vay,vaz] = nsimplex_3vec([bx_int(:,1)./sqrt(rho_int(:,1)),bx_int(:,2)./sqrt(rho_int(:,2)),bx_int(:,3)./sqrt(rho_int(:,3))],...
                [by_int(:,1)./sqrt(rho_int(:,1)),by_int(:,2)./sqrt(rho_int(:,2)),by_int(:,3)./sqrt(rho_int(:,3))],...
                [bz_int(:,1)./sqrt(rho_int(:,1)),bz_int(:,2)./sqrt(rho_int(:,2)),bz_int(:,3)./sqrt(rho_int(:,3))],x_int,y_int,z_int);

        elseif sc == 4 && fpi_sc == 3
            [px,py,pz] = nsimplex_4vec(bx_int,by_int,bz_int,x_int,y_int,z_int);
            [vx,vy,vz] = nsimplex_3vec(vx_int,vy_int,vz_int,x_int,y_int,z_int);
            rho = nsimplex_3(rho_int,x_int,y_int,z_int);
            temppar = nsimplex_3(temppar_int,x_int,y_int,z_int);
            tempperp = nsimplex_3(tempperp_int,x_int,y_int,z_int);
            etemppar = nsimplex_3(etemppar_int,x_int,y_int,z_int);
            etempperp = nsimplex_3(etempperp_int,x_int,y_int,z_int);
            p = nsimplex_3(p_int,x_int,y_int,z_int);
            ep = nsimplex_3(ep_int,x_int,y_int,z_int);

            fpi = find([1 2 3 4] ~= missing);
            [vax,vay,vaz] = nsimplex_3vec([bx_int(:,fpi(1))./sqrt(rho_int(:,1)),bx_int(:,fpi(2))./sqrt(rho_int(:,2)),bx_int(:,fpi(3))./sqrt(rho_int(:,3))],...
                [by_int(:,fpi(1))./sqrt(rho_int(:,1)),by_int(:,fpi(2))./sqrt(rho_int(:,2)),by_int(:,fpi(3))./sqrt(rho_int(:,3))],...
                [bz_int(:,fpi(1))./sqrt(rho_int(:,1)),bz_int(:,fpi(2))./sqrt(rho_int(:,2)),bz_int(:,fpi(3))./sqrt(rho_int(:,3))],x_int,y_int,z_int);
        elseif sc == 4 && fpi_sc == 2
            [px,py,pz] = nsimplex_4vec(bx_int,by_int,bz_int,x_int,y_int,z_int);
            vx = nanmean(vx_int,2);
            vy = nanmean(vy_int,2);
            vz = nanmean(vz_int,2);
            rho = nanmean(rho_int,2);
            temppar = nanmean(temppar_int,2);
            tempperp = nanmean(tempperp_int,2);
            etemppar = nanmean(etemppar_int,2);
            etempperp = nanmean(etempperp_int,2);
            p = nanmean(p_int,2);
            ep = nanmean(ep_int,2);

            fpi = find([1 2 3 4] ~= missing);
            vax = nanmean([bx_int(:,fpi(1))./sqrt(rho_int(:,1)),bx_int(:,fpi(2))./sqrt(rho_int(:,2))],2);
            vay = nanmean([by_int(:,fpi(1))./sqrt(rho_int(:,1)),by_int(:,fpi(2))./sqrt(rho_int(:,2))],2);
            vaz = nanmean([bz_int(:,fpi(1))./sqrt(rho_int(:,1)),bz_int(:,fpi(2))./sqrt(rho_int(:,2))],2);
        end
        temp = ((temppar+etemppar) + 2*(tempperp+etempperp))/3;
        p = (p + ep);

        bsqr = px.^2+py.^2+pz.^2;
        % only do this when picking out current sheets
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %current_sheets(jjj).intervals = find_cs(jt,bx_int1(:,1),by_int1(:,1),bz_int1(:,1));
        %save('current_sheets','current_sheets')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        intervals = current_sheets(jjj).intervals;
        [csc,~] = size(intervals);

        time_norms = current_sheets(jjj).timing_normals;
        %time_norms = zeros(csc,3); 
        press_balance_intervals = current_sheets(jjj).press_bal_ints; 
        %press_balance_intervals = zeros(csc,4);

%%%%%%%%%% just a test to see what happens with 1 spacecraft%%%%%%%%
%         fpi = find([1 2 3 4] ~= missing);
%         px = bx_int(:,fpi(1)); py = by_int(:,fpi(1)); pz = bz_int(:,fpi(1));
%         vx = vx_int(:,1); vy = vy_int(:,1); vz = vz_int(:,1);
%         rho = rho_int(:,1); temp = temp_int(:,1);
%         vax = bx_int(:,fpi(1))./sqrt(rho_int(:,fpi(1)));
%         vay = by_int(:,fpi(1))./sqrt(rho_int(:,fpi(1)));
%         vaz = bz_int(:,fpi(1))./sqrt(rho_int(:,fpi(1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %calculate de-hoffman teller velocity and subtract off
%         [V_DHT] = get_DHT_frame(px,py,pz,vx,vy,vz);
%         vx = vx - V_DHT(1);
%         vy = vy - V_DHT(2);
%         vz = vz - V_DHT(3);
% % 
%          figure('Visible','off','rend','painters','pos',[10 10 1000 800])
% %         %figure('rend','painters','pos',[10 10 1000 800])
%          subplot(4,1,1); hold on;
%          plot(t,sqrt(mean(bx_int.^2+by_int.^2+bz_int.^2,2)),'k')
%          plot_all_mms_B(t,bx_int,by_int,bz_int,sc)
%          lg = legend('|B|','B_x','B_y','B_z');
%          plot_help([datestr(todatenum(cdfepoch(t(1)))),' to ',...
%              datestr(todatenum(cdfepoch(t(end))))],t(1),t(end))
%          ylabel('nT')
%          subplot(4,1,2); hold on
%          plot(t,vx-nanmean(vx),'r')
%          plot(t,vy-nanmean(vy),'g','HandleVisibility','off')
%          plot(t,vz-nanmean(vz),'b','HandleVisibility','off')
%          plot(t,(1/1000)*(10^-9)*(vax-mean(vax))/sqrt(mu*(1.67*10^(-27))),'r--')
%          plot(t,(1/1000)*(10^-9)*(vay-mean(vay))/sqrt(mu*(1.67*10^(-27))),'g--','HandleVisibility','off')
%          plot(t,(1/1000)*(10^-9)*(vaz-mean(vaz))/sqrt(mu*(1.67*10^(-27))),'b--','HandleVisibility','off')
%          lg = legend('v_x','v_{Ax}');
%          ylabel('km/s')
%          plot_help('v and v_A',t(1),t(end))
%          subplot(4,1,3); hold on
%          plot(t,rho,'k')
%          plot_help('electron number density',t(1),t(end))
%          ylabel('m^{-3}')
%          subplot(4,1,4); hold on
%          %plot(t,(1e-18)*bsqr/(2*mu),'k')
%          %plot(t,(1.602e-19)*rho.*temp,'k')
%          plot(t,1e-9*p,'k')
%          %plot(t,(1e-18)*bsqr/(2*mu)+(1.602e-19)*rho.*temp,'k')
%          plot_help_nice('thermal pressure',t(1),t(end))
%          %legend('B^2/2\mu_0','p','p+B^2/2\mu_0')
%          ylabel('Pa')
%          subplot(5,1,5); hold on
%          if sc == 4
%             jx_int = interp1(jt,smooth(current_density(:,1),30),t);
%             jy_int = interp1(jt,smooth(current_density(:,2),30),t);
%             jz_int = interp1(jt,smooth(current_density(:,3),30),t);
%             jpar = sum([jx_int,jy_int,jz_int].*[px,py,pz],2)./sqrt(px.^2+py.^2+pz.^2);
%             jperp = sqrt(sum(cross([jx_int,jy_int,jz_int],[px,py,pz]).^2,2)./(px.^2+py.^2+pz.^2));
%             plot(t,jpar)
%             plot(t,jperp)
%          end
%          lg = legend('j_{||}','j_\perp');
         %plot_help_nice('current density magnetic field aligned coordinates',t(1),t(end))
            
%          subplot(4,1,1) 
%           fill([jti(w_ind111) jti(w_ind222) jti(w_ind222) jti(w_ind111)],[-5 -5  10 10],'c')
%           set(gca,'children',flipud(get(gca,'children')))
        %do analysis on each current sheet in this mms file
        for j = 1:csc
            cs = cs+1;
            tags(cs) = jjj;

            %get start and end indices for this event
            indd1 = floor(intervals(j,1) - intervals(j,2));
            indd2 = floor(intervals(j,1) + intervals(j,2));
            %indices in cs_inds are mag cadence so find nearest ion
            %cadence point
            ind1 = find(abs(t-jt(indd1)) == min(abs(t-jt(indd1))),1);
            ind2 = find(abs(jt(indd2)-t) == min(abs(jt(indd2)-t)),1);
            dis_range = ind1:ind2;

            %extract data for this event
            jti = t(dis_range);
            vxi = vx(dis_range);
            vyi = vy(dis_range);
            vzi = vz(dis_range);
            pxi = px(dis_range);
            pyi = py(dis_range);
            pzi = pz(dis_range);
            vaxi = (1/1000)*(10^-9)*(vax(dis_range))./sqrt(mu*(1.67*10^(-27)));
            vayi = (1/1000)*(10^-9)*(vay(dis_range))./sqrt(mu*(1.67*10^(-27)));
            vazi = (1/1000)*(10^-9)*(vaz(dis_range))./sqrt(mu*(1.67*10^(-27)));

            cs_vx(cs) = nanmean(vxi);
            cs_vy(cs) = nanmean(vyi);
            cs_vz(cs) = nanmean(vzi);

%              figure; hold on
%              plot(temp./max(temp),'k-.')
%              plot(rho,'k')
%              plot(dis_range,ones(length(jti),1)*mean(rho),'g','LineWidth',2)
%              bounds = ginput(4);
%              press_balance_intervals(j,1) = floor(bounds(1,1));
%              press_balance_intervals(j,2) = ceil(bounds(2,1));
%              press_balance_intervals(j,3) = floor(bounds(3,1));
%              press_balance_intervals(j,4) = ceil(bounds(4,1));

            left = press_balance_intervals(j,1):press_balance_intervals(j,2);
            right = press_balance_intervals(j,3):press_balance_intervals(j,4);
            %left = ind1-5:ind1;
            %right = ind2:ind2+5;
            middle = ind1:ind2;
            %subplot(4,1,3);
            %plot(t(left),rho(left),'r','LineWidth',2)
            %plot(t(right),rho(right),'b','LineWidth',2)

            cs_rho_enh(cs,1) = mean(rho(left));
            cs_rho_enh(cs,2) = mean(rho(right));
            cs_rho_enh(cs,3) = mean(rho(middle));
            cs_temp_enh(cs,1) = mean(temp(left));
            cs_temp_enh(cs,2) = mean(temp(right));
            cs_temp_enh(cs,3) = mean(temp(middle));
            cs_p_enh(cs,1) = mean(p(left));
            cs_p_enh(cs,2) = mean(p(right));
            cs_p_enh(cs,3) = mean(p(middle));

            cs_Bmag(cs,1) = mean(bsqr(left));
            cs_Bmag(cs,2) = mean(bsqr(right));

            press_balls_std_T(cs,1) = std(temp(left));
            press_balls_std_T(cs,2) = std(temp(right));
            press_balls_std(cs,1) = std(p(left));
            press_balls_std(cs,2) = std(p(right));

            %subplot(4,1,1);
%            highlight whole event in black
            %plot(jti,pxi,'k','LineWidth',2)
            %plot(jti,pyi,'k','LineWidth',2)
            %plot(jti,pzi,'k','LineWidth',2)

            %do walen relation analysis
            [MVA1,MVA2,MVA3] = normal_coords(pxi,pyi,pzi);

            brotated_x = zeros(size(pxi));
            brotated_y = zeros(size(pxi));
            brotated_z = zeros(size(pxi));
            vrotated_x = zeros(size(pxi));
            vrotated_y = zeros(size(pxi));
            vrotated_z = zeros(size(pxi));
            varotated_x = zeros(size(pxi));
            varotated_y = zeros(size(pxi));
            varotated_z = zeros(size(pxi));
            for kle = 1:length(brotated_x)
                brotnew = [MVA1(:),MVA2(:),MVA3(:)]*[pxi(kle);pyi(kle);pzi(kle)];
                vrotnew = [MVA1(:),MVA2(:),MVA3(:)]*[vxi(kle);vyi(kle);vzi(kle)];
                varotnew = [MVA1(:),MVA2(:),MVA3(:)]*[vaxi(kle);vayi(kle);vazi(kle)];
                brotated_x(kle) = brotnew(1);
                brotated_y(kle) = brotnew(2);
                brotated_z(kle) = brotnew(3);
                vrotated_x(kle) = vrotnew(1);
                vrotated_y(kle) = vrotnew(2);
                vrotated_z(kle) = vrotnew(3);
                varotated_x(kle) = varotnew(1);
                varotated_y(kle) = varotnew(2);
                varotated_z(kle) = varotnew(3);
            end

%             figure
%             subplot(3,1,1); hold on
%             plot(brotated_x)
%             plot(brotated_y)
%             plot(brotated_z)
%             subplot(3,1,2); hold on
%             plot(vrotated_x-mean(vrotated_x))
%             plot(vrotated_y-mean(vrotated_y))
%             plot(vrotated_z-mean(vrotated_z))
%             subplot(3,1,3); hold on
%             plot(varotated_x)
%             plot(varotated_y)
%             plot(varotated_z)
            
            %[v_jumps2(cs,:),va_jumps2(cs,:),w_ind111,w_ind222] = walen(jti,vxi,vyi,vzi,vaxi,vayi,vazi,pxi,pyi,pzi,0.2);
            [v_jumps2(cs,:),va_jumps2(cs,:),w_ind111,w_ind222] = walen(jti,vrotated_x,vrotated_y,vrotated_z,varotated_x,varotated_y,varotated_z,brotated_x,brotated_y,brotated_z,0.2);
            cs_wal_len(cs) = w_ind222-w_ind111;
            %if w_ind111 > 0
            %    plot(jti(w_ind111:w_ind222),pxi(w_ind111:w_ind222),'c','LineWidth',3)
            %end

            %use single spacecraft and fast mag cadence for MVA analysis
            bxx = bx_int1(:,1);
            byy = by_int1(:,1);
            bzz = bz_int1(:,1);

            er = 0;
            wl_mva = 0;
            wind_lengg = floor((indd2-indd1)):floor((indd2-indd1)/5):(5*(indd2-indd1));
            %perform MVA analysis for progressively larger windows
            %until we get eigenvalues ratio > 10
            for wind_leng = wind_lengg
                if ((indd1-floor(wind_leng/2)) > 0 && (indd2+floor(wind_leng/2) < length(bxx)))
                    bxs = bxx(indd1-floor(wind_leng/2):indd2+floor(wind_leng/2));
                    bys = byy(indd1-floor(wind_leng/2):indd2+floor(wind_leng/2));
                    bzs = bzz(indd1-floor(wind_leng/2):indd2+floor(wind_leng/2));
                    [ndv1,ndv2,eigrat] = normal_dir_var(bxs,bys,bzs);
                    if eigrat > er
                        er = eigrat;
                        n1 = ndv1;
                        n2 = ndv2;
                        wl_mva = wind_leng;
                        if er > 10
                            break
                        end
                    end
                end
            end
            %find boundary normal by timing method for comparison
            %use this part to find all normals, else uncomment line after
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if sc == 4
%                 nt = normal_dir_timing4(jt(indd1:indd2),bx_int1(indd1:indd2,:),by_int1(indd1:indd2,:),bz_int1(indd1:indd2,:),...
%                     x_int1(indd1:indd2,:),y_int1(indd1:indd2,:),z_int1(indd1:indd2,:));
%             elseif sc == 3
%                 nt = normal_dir_timing3(jt(indd1:indd2),bx_int1(indd1:indd2,:),by_int1(indd1:indd2,:),bz_int1(indd1:indd2,:),...
%                     x_int1(indd1:indd2,:),y_int1(indd1:indd2,:),z_int1(indd1:indd2,:));
%             end
%             time_norms(j,:) = nt;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nt = time_norms(j,:); 
            
            %Determine B_|| and B_\perp in boundary normal coordinates 
            interval = ind1:ind2;

            dir_cs1 = ones(length(px),3); 
            dir_cs1(1,:) = [n1(1),n1(2),n1(3)];
            dir_cs1(2:end,1) = dir_cs1(2:end,1)*dir_cs1(1,1);
            dir_cs1(2:end,2) = dir_cs1(2:end,2)*dir_cs1(1,2);
            dir_cs1(2:end,3) = dir_cs1(2:end,3)*dir_cs1(1,3);

            dir_cs2 = ones(length(px),3); 
            dir_cs2(1,:) = [0,0,1];
            dir_cs2(2:end,1) = dir_cs2(2:end,1)*dir_cs2(1,1);
            dir_cs2(2:end,2) = dir_cs2(2:end,2)*dir_cs2(1,2);
            dir_cs2(2:end,3) = dir_cs2(2:end,3)*dir_cs2(1,3);

            bpar1vec = abs(sum([px,py,pz].*dir_cs1,2));
            bpar2vec = abs(sum([px,py,pz].*dir_cs2,2));
            minbpar1 = min(bpar1vec(interval));
            minbpar2 =  min(bpar2vec(interval));
            bpar1ind = find(bpar1vec == minbpar1);
            bpar2ind = find(bpar2vec == minbpar2);
            bpar1 = mean(bpar1vec((bpar1ind-6):(bpar1ind+6)));
            bpar2 = mean(bpar2vec((bpar2ind-6):(bpar2ind+6)));

            if (bpar1 > bpar2)
                dir_cs1 = dir_cs2;
                bpar1ind = bpar2ind;
                n1 = n2;
            end 
            bxd = cross([px,py,pz],dir_cs1);

            cs_Bpar(cs) = min([bpar1,bpar2]);
            cs_Bperp(cs) = mean(sqrt(bxd((bpar1ind-6):(bpar1ind+6),1).^2+bxd((bpar1ind-6):(bpar1ind+6),2).^2+bxd((bpar1ind-6):(bpar1ind+6),3).^2));
            cs_eigrats(cs) = er;
            cs_n1_nt(cs) = sum(n1.*nt);
            cs_norm_vx(cs) = acos(n1(1)/(n1(1)^2+n1(2)^2+n1(3)^2));
            cs_len(cs) = 2*intervals(j,2);

            cs_jumps_tot(cs) = sqrt((max(pxi)-min(pxi))^2+(max(pyi)-min(pyi))^2+(max(pzi)-min(pzi))^2); 
%             v_jumps_totx(cs) =  max(vxi)-min(vxi); 
%             v_jumps_toty(cs) =  max(vyi)-min(vyi); 
%             v_jumps_totz(cs) =  max(vzi)-min(vzi); 
%             va_jumps_totx(cs) =  max(vaxi)-min(vaxi); 
%             va_jumps_toty(cs) =  max(vayi)-min(vayi); 
%             va_jumps_totz(cs) =  max(vazi)-min(vazi); 
            cs_shear(cs) = (mean(px(left))*mean(px(right))+mean(py(left))*mean(py(right))+mean(pz(left))*mean(pz(right)))/...
                sqrt((mean(px(left))^2+mean(py(left))^2+mean(pz(left))^2)*(mean(px(right))^2+mean(py(right))^2+mean(pz(right))^2));
            v_jumps_totx(cs) =  max(vrotated_x)-min(vrotated_x); 
            v_jumps_toty(cs) =  max(vrotated_y)-min(vrotated_y); 
            v_jumps_totz(cs) =  max(vrotated_z)-min(vrotated_z); 
            va_jumps_totx(cs) =  max(varotated_x)-min(varotated_x); 
            va_jumps_toty(cs) =  max(varotated_y)-min(varotated_y); 
            va_jumps_totz(cs) =  max(varotated_z)-min(varotated_z); 
            cs_loc(cs,1) = x_int(1); cs_loc(cs,2) = y_int(1); cs_loc(cs,3) = z_int(1);
        end
       %saveas(gcf,['ion_cadence/discont1_',num2str(jjj),'.png'])
       %close all
     end
     %current_sheets(jjj).timing_normals = time_norms;
     %current_sheets(jjj).press_bal_ints = press_balance_intervals;
end
 save('cs_vx','cs_vx')
 save('cs_vy','cs_vy')
 save('cs_vz','cs_vz')
 save('va_jumps_totx','va_jumps_totx')
 save('va_jumps_toty','va_jumps_toty')
 save('va_jumps_totz','va_jumps_totz')
 save('v_jumps_totx','v_jumps_totx')
 save('v_jumps_toty','v_jumps_toty')
 save('v_jumps_totz','v_jumps_totz')
 save('v_jumps2','v_jumps2')
 save('va_jumps2','va_jumps2')
 save('cs_jumps_tot','cs_jumps_tot')
 save('tags','tags')
 save('cs_loc','cs_loc')    
 save('cs_Bpar','cs_Bpar')    
 save('cs_Bperp','cs_Bperp')    
 save('cs_eigrats','cs_eigrats')    
 save('cs_n1_nt','cs_n1_nt')   
 save('cs_rho_enh','cs_rho_enh')
 save('cs_p_enh','cs_p_enh') 
 save('cs_temp_enh','cs_temp_enh')
 save('cs_norm_vx','cs_norm_vx')   
 save('cs_wal_len','cs_wal_len')
 save('press_balls_std','press_balls_std')
 save('press_balls_std_T','press_balls_std_T')
 save('cs_Bmag','cs_Bmag')
 save('cs_shear','cs_shear')
save('cs_len','cs_len')
% open_boundaries
%  toc
 
