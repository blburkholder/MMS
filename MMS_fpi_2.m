tic
%evaluates walen relation at cadence of ion FPI measurements
mu = 4*pi*10^(-7);
%this file contains the start index and end index for current sheet
%intervals within the mms data files
load('current_sheets.mat')
%there are 454 current sheets within the 303 files
walegood = zeros(363,1);
cs_vx = zeros(363,1);
cs_vy = zeros(363,1);
cs_vz = zeros(363,1);
v_jumps_totx = zeros(363,1);
v_jumps_toty = zeros(363,1);
v_jumps_totz = zeros(363,1);
va_jumps_totx = zeros(363,1);
va_jumps_toty = zeros(363,1);
va_jumps_totz = zeros(363,1);
cs_jumps_tot = zeros(363,1);
v_jumps1 = zeros(363,3);
v_jumps2 = zeros(363,3);
va_jumps1 = zeros(363,3);
va_jumps2 = zeros(363,3);
cs_loc = zeros(363,3);
tags = zeros(363,1);
cs_eigrats = zeros(363,1);
cs_Bperp = zeros(363,1);
cs_Bpar = zeros(363,1);
cs_n1_nt = zeros(363,1);
cs_cond = zeros(363,1);
cs_rho_enh = zeros(363,3);
cs_temp_enh = zeros(363,3);
cs_norm_vx = zeros(363,1);
cs_wal_len = zeros(363,1);
cs_beta = zeros(363,1);
cs = 0;

%201 & 225 missing
for jjj = [1:200,202:224,226:303]
%for jjj = 251
     jjj
    BB = current_sheets(jjj).mag_data;
    vv = current_sheets(jjj).v_data;
    pospos = current_sheets(jjj).pos;
    rho_comp = current_sheets(jjj).rho_data;
    temp_comp = current_sheets(jjj).temp_data;
    etemp_comp = current_sheets(jjj).etemp_data;
    ti = current_sheets(jjj).ti;
    te = current_sheets(jjj).te;
    jt = current_sheets(jjj).jt;
    fpi_sc = current_sheets(jjj).fpi_sc;
    sc = current_sheets(jjj).sc;
    missing = current_sheets(jjj).m;

    bx_int1 = BB(:,:,1); by_int1 = BB(:,:,2); bz_int1 = BB(:,:,3);
    vx_comp = vv(:,:,1); vy_comp = vv(:,:,2); vz_comp = vv(:,:,3);
    x_int1 = pospos(:,:,1); y_int1 = pospos(:,:,2); z_int1 = pospos(:,:,3);
    
    %interpolating FPI data to same time stamps
     if fpi_sc > 1
        [vx_int1,vy_int1,vz_int1,it] =...
            interp_time_v(vx_comp,vy_comp,vz_comp,ti,fpi_sc);
        [rho_int1,~,~,et] =...
            interp_time_v(rho_comp,rho_comp,rho_comp,te,fpi_sc);
        [temp_int1,~,~,~] =...
            interp_time_v(temp_comp,temp_comp,temp_comp,ti,fpi_sc);
        [etemp_int1,~,~,~] =...
            interp_time_v(etemp_comp,etemp_comp,etemp_comp,te,fpi_sc);

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
        temp_int = zeros(length(t),fpi_sc);

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
            temp_int(:,h) = interp1(it,temp_int1(:,h),t)+interp1(et,etemp_int1(:,h),t);
        end

        %interpolating measurements from all available spacecraft to
        %center of formation
        if sc == 4 && fpi_sc == 4
            [px,py,pz] = nsimplex_4vec(bx_int,by_int,bz_int,x_int,y_int,z_int);
            [vx,vy,vz] = nsimplex_4vec(vx_int,vy_int,vz_int,x_int,y_int,z_int);
            rho = nsimplex_4(rho_int,x_int,y_int,z_int);
            temp = nsimplex_4(temp_int,x_int,y_int,z_int);

            [vax,vay,vaz] = nsimplex_4vec([bx_int(:,1)./sqrt(rho_int(:,1)),bx_int(:,2)./sqrt(rho_int(:,2)),bx_int(:,3)./sqrt(rho_int(:,3)),bx_int(:,4)./sqrt(rho_int(:,4))],...
                [by_int(:,1)./sqrt(rho_int(:,1)),by_int(:,2)./sqrt(rho_int(:,2)),by_int(:,3)./sqrt(rho_int(:,3)),by_int(:,4)./sqrt(rho_int(:,4))],...
                [bz_int(:,1)./sqrt(rho_int(:,1)),bz_int(:,2)./sqrt(rho_int(:,2)),bz_int(:,3)./sqrt(rho_int(:,3)),bz_int(:,4)./sqrt(rho_int(:,4))],x_int,y_int,z_int);
        elseif sc == 3 && fpi_sc == 3
            [px,py,pz] = nsimplex_3vec(bx_int,by_int,bz_int,x_int,y_int,z_int);  
            [vx,vy,vz] = nsimplex_3vec(vx_int,vy_int,vz_int,x_int,y_int,z_int);
            rho = nsimplex_3(rho_int,x_int,y_int,z_int);
            temp = nsimplex_3(temp_int,x_int,y_int,z_int);

            [vax,vay,vaz] = nsimplex_3vec([bx_int(:,1)./sqrt(rho_int(:,1)),bx_int(:,2)./sqrt(rho_int(:,2)),bx_int(:,3)./sqrt(rho_int(:,3))],...
                [by_int(:,1)./sqrt(rho_int(:,1)),by_int(:,2)./sqrt(rho_int(:,2)),by_int(:,3)./sqrt(rho_int(:,3))],...
                [bz_int(:,1)./sqrt(rho_int(:,1)),bz_int(:,2)./sqrt(rho_int(:,2)),bz_int(:,3)./sqrt(rho_int(:,3))],x_int,y_int,z_int);

        elseif sc == 4 && fpi_sc == 3
            [px,py,pz] = nsimplex_4vec(bx_int,by_int,bz_int,x_int,y_int,z_int);
            [vx,vy,vz] = nsimplex_3vec(vx_int,vy_int,vz_int,x_int,y_int,z_int);
            rho = nsimplex_3(rho_int,x_int,y_int,z_int);
            temp = nsimplex_3(temp_int,x_int,y_int,z_int);

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
            temp = nanmean(temp_int,2);

            fpi = find([1 2 3 4] ~= missing);
            vax = nanmean([bx_int(:,fpi(1))./sqrt(rho_int(:,1)),bx_int(:,fpi(2))./sqrt(rho_int(:,2))],2);
            vay = nanmean([by_int(:,fpi(1))./sqrt(rho_int(:,1)),by_int(:,fpi(2))./sqrt(rho_int(:,2))],2);
            vaz = nanmean([bz_int(:,fpi(1))./sqrt(rho_int(:,1)),bz_int(:,fpi(2))./sqrt(rho_int(:,2))],2);
        end
        % only do this when picking out current sheets
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %current_sheets(jjj).intervals = find_cs(jt,bx_int1(:,1),by_int1(:,1),bz_int1(:,1));
        %save('current_sheets','current_sheets')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        intervals = current_sheets(jjj).intervals;
        [csc,~] = size(intervals);

%         fpi = find([1 2 3 4] ~= missing);
%         px = bx_int(:,fpi(1)); py = by_int(:,fpi(1)); pz = bz_int(:,fpi(1));
%         vx = vx_int(:,1); vy = vy_int(:,1); vz = vz_int(:,1);
%         rho = rho_int(:,1); temp = temp_int(:,1);
%         vax = bx_int(:,fpi(1))./sqrt(rho_int(:,fpi(1)));
%         vay = by_int(:,fpi(1))./sqrt(rho_int(:,fpi(1)));
%         vaz = bz_int(:,fpi(1))./sqrt(rho_int(:,fpi(1)));


        %calculate de-hoffman teller velocity and subtract off
        [V_DHT] = get_DHT_frame(px,py,pz,vx,vy,vz);
        vx = vx - V_DHT(1);
        vy = vy - V_DHT(2);
        vz = vz - V_DHT(3);

%         %figure('Visible','off','rend','painters','pos',[10 10 1000 800])
%         figure('rend','painters','pos',[10 10 1000 800])
%         subplot(4,1,1); hold on;
%         plot_all_mms_B(t,bx_int,by_int,bz_int,sc)
%         lg = legend('B_x','B_y','B_z','Location','northeastoutside');
%         plot_help([datestr(todatenum(cdfepoch(t(1)))),' to ',...
%             datestr(todatenum(cdfepoch(t(end))))],t(1),t(end))
%         ylabel('nT')
%         subplot(4,1,2); hold on
%         plot(t,vx-nanmean(vx),'r')
%         plot(t,vy-nanmean(vy),'g','HandleVisibility','off')
%         plot(t,vz-nanmean(vz),'b','HandleVisibility','off')
%         plot(t,(1/1000)*(10^-9)*(vax-mean(vax))/sqrt(mu*(1.67*10^(-27))),'r--')
%         plot(t,(1/1000)*(10^-9)*(vay-mean(vay))/sqrt(mu*(1.67*10^(-27))),'g--','HandleVisibility','off')
%         plot(t,(1/1000)*(10^-9)*(vaz-mean(vaz))/sqrt(mu*(1.67*10^(-27))),'b--','HandleVisibility','off')
%         lg = legend('v_x','v_{Ax}','Location','northeastoutside');
%         ylabel('km/s')
%         plot_help('v and v_A',t(1),t(end))
%         subplot(4,1,3); hold on
%         plot(t,rho,'k')
%         plot_help('electron density',t(1),t(end))
%         legend('\rho^-','Location','northeastoutside')
%         ylabel('m^{-3}')
%         subplot(4,1,4)
%         plot(t,temp,'k')
%         plot_help('ion temperature',t(1),t(end))
%         legend('T','Location','northeastoutside')
%         ylabel('eV')
%         subplot(4,1,1) 
% %        fill([jti(w_ind111) jti(w_ind222) jti(w_ind222) jti(w_ind111)],[-10 -10  5 5],'c')
%         set(gca,'children',flipud(get(gca,'children')))

        %do analysis on each current sheet in this mms file
        for j = 1:csc
            cs = cs+1;
            tags(cs) = jjj;

            %get start and end indices for this event
            mid_point = intervals(j,1);
            indd1 = floor(intervals(j,1) - intervals(j,2));
            indd2 = floor(intervals(j,1) + intervals(j,2));
            %indices in cs_inds are for mag cadence so find nearest ion
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
            cs_rho_enh(cs,1) = mean((rho(ind1)-2):(rho(ind1)+2));
            cs_rho_enh(cs,2) = mean((rho(ind2)-2):(rho(ind2)+2));
            cs_rho_enh(cs,3) = mean(rho((floor((ind1+ind2)/2)-2):(floor((ind1+ind2)/2)+2)));
            cs_temp_enh(cs,1) = mean((temp(ind1)-2):(temp(ind1)+2));
            cs_temp_enh(cs,2) = mean((temp(ind2)-2):(temp(ind2)+2));
            cs_temp_enh(cs,3) = mean(temp((floor((ind1+ind2)/2)-2):(floor((ind1+ind2)/2)+2)));

            %highlight whole event in black
            %plot(jti,pxi,'k','LineWidth',2)
            %plot(jti,pyi,'k','LineWidth',2)
            %plot(jti,pzi,'k','LineWidth',2)

            %do walen relation analysis for event
            [v_jumps2(cs,:),va_jumps2(cs,:),walegood(cs),w_ind111,w_ind222] = walen(jti,vxi,vyi,vzi,vaxi,vayi,vazi,...
                pxi,pyi,pzi,0.2);
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
            if sc == 4
                [nt,cond] = normal_dir_timing4(jt(indd1:indd2),bx_int1(indd1:indd2,:),by_int1(indd1:indd2,:),bz_int1(indd1:indd2,:),...
                    x_int1(indd1:indd2,:),y_int1(indd1:indd2,:),z_int1(indd1:indd2,:));
            elseif sc == 3
                [nt,cond] = normal_dir_timing3(jt(indd1:indd2),bx_int1(indd1:indd2,:),by_int1(indd1:indd2,:),bz_int1(indd1:indd2,:),...
                    x_int1(indd1:indd2,:),y_int1(indd1:indd2,:),z_int1(indd1:indd2,:));
            end

            cs_norm_vx(cs) = acos(n1(1)/(n1(1)^2+n1(2)^2+n1(3)^2));

            %Determine B_|| and B_\perp in boundary normal coordinates 
            if w_ind111 == 0
                w_ind111 = ind1 + floor((ind2-ind1)/4); w_ind222 = ind2 - floor((ind2-ind1)/4);

                dir_cs1 = ones(length(w_ind111:w_ind222),3); 
                dir_cs1(1,:) = [n1(1),n1(2),n1(3)];
                dir_cs1(2:end,1) = dir_cs1(2:end,1)*dir_cs1(1,1);
                dir_cs1(2:end,2) = dir_cs1(2:end,2)*dir_cs1(1,2);
                dir_cs1(2:end,3) = dir_cs1(2:end,3)*dir_cs1(1,3);

                dir_cs2 = ones(length(w_ind111:w_ind222),3); 
                dir_cs2(1,:) = [n2(1),n2(2),n2(3)];
                dir_cs2(2:end,1) = dir_cs2(2:end,1)*dir_cs2(1,1);
                dir_cs2(2:end,2) = dir_cs2(2:end,2)*dir_cs2(1,2);
                dir_cs2(2:end,3) = dir_cs2(2:end,3)*dir_cs2(1,3);

                bpar1 = min(abs(sum([px(w_ind111:w_ind222),py(w_ind111:w_ind222),pz(w_ind111:w_ind222)].*dir_cs1,2)));
                bpar2 =  min(abs(sum([px(w_ind111:w_ind222),py(w_ind111:w_ind222),pz(w_ind111:w_ind222)].*dir_cs2,2)));
                if (bpar1 > bpar2)
                    dir_cs1 = dir_cs2;
                    n1 = n2;
                end
                bxd = cross([px(w_ind111:w_ind222),py(w_ind111:w_ind222),pz(w_ind111:w_ind222)],dir_cs1);
            else
                dir_cs1 = ones(length(w_ind111:w_ind222),3); 
                dir_cs1(1,:) = [n1(1),n1(2),n1(3)];
                dir_cs1(2:end,1) = dir_cs1(2:end,1)*dir_cs1(1,1);
                dir_cs1(2:end,2) = dir_cs1(2:end,2)*dir_cs1(1,2);
                dir_cs1(2:end,3) = dir_cs1(2:end,3)*dir_cs1(1,3);

                dir_cs2 = ones(length(w_ind111:w_ind222),3); 
                dir_cs2(1,:) = [n2(1),n2(2),n2(3)];
                dir_cs2(2:end,1) = dir_cs2(2:end,1)*dir_cs2(1,1);
                dir_cs2(2:end,2) = dir_cs2(2:end,2)*dir_cs2(1,2);
                dir_cs2(2:end,3) = dir_cs2(2:end,3)*dir_cs2(1,3);

                bpar1 = min(abs(sum([pxi(w_ind111:w_ind222),pyi(w_ind111:w_ind222),pzi(w_ind111:w_ind222)].*dir_cs1,2)));
                bpar2 =  min(abs(sum([pxi(w_ind111:w_ind222),pyi(w_ind111:w_ind222),pzi(w_ind111:w_ind222)].*dir_cs2,2)));
                if (bpar1 > bpar2)
                    dir_cs1 = dir_cs2;
                    n1 = n2;
                end
                bxd = cross([pxi(w_ind111:w_ind222),pyi(w_ind111:w_ind222),pzi(w_ind111:w_ind222)],dir_cs1);
            end

            cs_Bpar(cs) = min([bpar1,bpar2]);
            cs_Bperp(cs) = mean(sqrt(bxd(:,1).^2+bxd(:,2).^2+bxd(:,3).^2));
            cs_eigrats(cs) = er;
            cs_n1_nt(cs) = sum(n1.*nt');
            cs_cond(cs) = cond;
            cs_beta(cs) = (1.602e-19)*2*mu*mean(cs_rho_enh(cs,:))*mean(cs_temp_enh(cs,:))/((1e-18)*(mean(pxi.^2+pyi.^2+pzi.^2)));

            cs_jumps_tot(cs) = sqrt((max(pxi)-min(pxi))^2+(max(pyi)-min(pyi))^2+(max(pzi)-min(pzi))^2); 
            v_jumps_totx(cs) =  max(vxi)-min(vxi); 
            v_jumps_toty(cs) =  max(vyi)-min(vyi); 
            v_jumps_totz(cs) =  max(vzi)-min(vzi); 
            va_jumps_totx(cs) =  max(vaxi)-min(vaxi); 
            va_jumps_toty(cs) =  max(vayi)-min(vayi); 
            va_jumps_totz(cs) =  max(vazi)-min(vazi); 
            cs_loc(cs,1) = x_int(1); cs_loc(cs,2) = y_int(1); cs_loc(cs,3) = z_int(1);
        end
       %saveas(gcf,['ion_cadence/discont_',num2str(jjj),'.png'])
       %close all
     end
end
 save('walegood','walegood')
 save('cs_vx','cs_vx')
 save('cs_vy','cs_vy')
 save('cs_vz','cs_vz')
 save('va_jumps_totx','va_jumps_totx')
 save('va_jumps_toty','va_jumps_toty')
 save('va_jumps_totz','va_jumps_totz')
 save('v_jumps_totx','v_jumps_totx')
 save('v_jumps_toty','v_jumps_toty')
 save('v_jumps_totz','v_jumps_totz')
 save('v_jumps1','v_jumps1')
 save('v_jumps2','v_jumps2')
 save('va_jumps1','va_jumps1')
 save('va_jumps2','va_jumps2')
 save('cs_jumps_tot','cs_jumps_tot')
 save('tags','tags')
 save('cs_loc','cs_loc')    
 save('cs_Bpar','cs_Bpar')    
 save('cs_Bperp','cs_Bperp')    
 save('cs_eigrats','cs_eigrats')    
 save('cs_n1_nt','cs_n1_nt')
 save('cs_cond','cs_cond')    
 save('cs_rho_enh','cs_rho_enh')
 save('cs_temp_enh','cs_temp_enh')   
 save('cs_norm_vx','cs_norm_vx')   
 save('cs_wal_len','cs_wal_len')
 save('cs_beta','cs_beta')
 open_boundaries
toc
 
