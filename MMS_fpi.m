%evaluates walen relation at superfast mag cadence    
mu = 4*pi*10^(-7);

    fid = fopen('mms_datas.txt');
    A = fscanf(fid,'%s\n');
    load('cs_inds.mat')
    %cs_Bperp = zeros(455,1);
    %cs_Bpar = zeros(455,1);
    cs_walen1 = zeros(455,1);
    cs_walen2 = zeros(455,1);
    walegood1 = zeros(455,1);
    walegood2 = zeros(455,1);
    %cs_eigrats = zeros(455,1);
    cs_vx = zeros(455,1);
    cs_vy = zeros(455,1);
    cs_vz = zeros(455,1);
    v_jumps_totx = zeros(455,1);
    v_jumps_toty = zeros(455,1);
    v_jumps_totz = zeros(455,1);
    va_jumps_totx = zeros(455,1);
    va_jumps_toty = zeros(455,1);
    va_jumps_totz = zeros(455,1);
    cs_jumps_tot = zeros(455,1);
    v_jumps1 = zeros(455,3);
    v_jumps2 = zeros(455,3);
    va_jumps1 = zeros(455,3);
    va_jumps2 = zeros(455,3);
    %ion_beta = zeros(455,1);
    cs_loc = zeros(455,3);
    tags = zeros(455,1);
    cs = 0;
    skipthis = 0;

    %201 & 225 missing
   %for jjj = [1:200,202:224,226:303]
   for jjj = 251
        jjj
        datey = A((jjj-1)*34+1:jjj*34);
        date2 = datey(17:end-4);
        date1 = [date2(1:4),'/',date2(5:6),'/',date2(7:8)];
        path1 = '/home/computation/Documents/Sun/MMS/mms_data/mms/data/mms';
        e_path1 = ['/fpi/brst/l2/des-moms/',date1,'/mms'];
        e_path2 = ['_fpi_brst_l2_des-moms_',date2,'_v3.3.0.cdf'];
        i_path1 =  ['/fpi/brst/l2/dis-moms/',date1,'/mms'];
        i_path2 =  ['_fpi_brst_l2_dis-moms_',date2,'_v3.3.0.cdf'];
        B_path1 = ['/fgm/brst/l2/',date1,'/mms'];
        B_path2 =  ['_fgm_brst_l2_',date2];

        [j1,j2,bx_int1,by_int1,bz_int1,x_int1,y_int1,z_int1,jt,B_tag,sc] = cur_dens(path1,B_path1,B_path2);
        [wl,csc] = find_cs(jt,bx_int1(:,1),by_int1(:,1),bz_int1(:,1),jjj);

%         if sc == 4
%             jxm = j1(:,1); jym = j1(:,2); jzm = j1(:,3);
%         end
%         jxs = j2(:,1); jys = j2(:,2); jzs = j2(:,3);
% 
        B_path2 =  [B_path2,'_v5.',num2str(B_tag),'.0.cdf'];
        missing = 0;
        vx_comp = []; vy_comp = []; vz_comp = []; rho_comp = []; p_comp = []; temp_comp = []; temperr_comp = [];
        ion_dens = []; ti = []; te = [];

        fpi_sc = 0;
        for i = 1:4
            try
                bx = bx_int1(:,i);
            catch
                fprintf(['only ',num2str(sc), ' MMS FGM files available\n']);
            end
            try
                [electro_data,~] = spdfcdfread([path1,num2str(i),e_path1,num2str(i),e_path2]);
                [data,~] = spdfcdfread([path1,num2str(i),i_path1,num2str(i),i_path2]);
                fpi_sc = fpi_sc + 1;

                ti = data{1};
                te = electro_data{1};
                %ions = data{19};
                numberdensity = electro_data{23}; %cm^-3
                %ebulkv = electro_data{29} - electro_data{30};
                bulkv = data{25}; %GSE
%                 prestens = data{29}; %nPa
                temptens = data{33}; %
% 
%                 prees = zeros(length(t),1);
                 teemp = zeros(length(t),1);
                 for iii = 1:length(t)
%                     prees(iii) = trace(prestens(:,:,iii));
                     teemp(iii) = trace(temptens(:,:,iii));
                 end

                wind = 20+ceil(length(bulkv(:,1))/50);
                vx_comp(1:length(jt),fpi_sc) = smooth(interp1(ti,hamperl(fft(bulkv(:,1)),wind),jt),128);
                vy_comp(1:length(jt),fpi_sc) = smooth(interp1(ti,hamperl(fft(bulkv(:,2)),wind),jt),128);
                vz_comp(1:length(jt),fpi_sc) = smooth(interp1(ti,hamperl(fft(bulkv(:,3)),wind),jt),128);
                rho_comp(1:length(jt),fpi_sc) = (10^6)*smooth(interp1(te,numberdensity,jt),128);

                %ion_dens(1:length(jt),fpi_sc) = (10^6)*smooth(interp1(t,hamperl(fft(ions),wind),jt),128);
                %p_comp(1:length(jt),fpi_sc) = smooth(interp1(t,hamperl(fft(prees),wind),jt),128);
                temp_comp(1:length(jt),fpi_sc) = smooth(interp1(ti,hamperl(fft(teemp),wind),jt),128);
            catch
                fprintf(['no MMS',num2str(i), 'FPI\n']);
                missing = i;
            end
        end  

         if fpi_sc > 1
            bx_int = bx_int1; by_int = by_int1; bz_int = bz_int1;
            vx_int = vx_comp; vy_int = vy_comp; vz_int = vz_comp;
            x_int = x_int1; y_int = y_int1; z_int = z_int1;
            rho_int = rho_comp;
% %             jx = [];
% %             jy = [];
% %             jz = [];
            if sc == 4 && fpi_sc == 4
                [px,py,pz] = nsimplex_4vec(bx_int,by_int,bz_int,x_int,y_int,z_int);
                [vx,vy,vz] = nsimplex_4vec(vx_int,vy_int,vz_int,x_int,y_int,z_int);
                rho = nsimplex_4(rho_int,x_int,y_int,z_int);
                %ion_rho = nsimplex_4(ion_dens,x_int,y_int,z_int);
                %pres = nsimplex_4(p_comp,x_int,y_int,z_int);
                %temp = nsimplex_4(temp_comp,x_int,y_int,z_int);
                [vax,vay,vaz] = nsimplex_4vec([bx_int(:,1)./sqrt(rho_int(:,1)),bx_int(:,2)./sqrt(rho_int(:,2)),bx_int(:,3)./sqrt(rho_int(:,3)),bx_int(:,4)./sqrt(rho_int(:,4))],...
                    [by_int(:,1)./sqrt(rho_int(:,1)),by_int(:,2)./sqrt(rho_int(:,2)),by_int(:,3)./sqrt(rho_int(:,3)),by_int(:,4)./sqrt(rho_int(:,4))],...
                    [bz_int(:,1)./sqrt(rho_int(:,1)),bz_int(:,2)./sqrt(rho_int(:,2)),bz_int(:,3)./sqrt(rho_int(:,3)),bz_int(:,4)./sqrt(rho_int(:,4))],x_int,y_int,z_int);
%                 subplot(3,2,2); hold on;
%                 plot(jt,jxm,'r'); plot(jt,jym,'g'); plot(jt,jzm,'b');
%                 legend('jx','jy','jz')
%                 plot_help('current density curlometer',jt(1),jt(end));
%                 jx = jxm;
%                 jy = jym;
%                 jz = jzm;
            elseif sc == 3 && fpi_sc == 3
                [px,py,pz] = nsimplex_3vec(bx_int,by_int,bz_int,x_int,y_int,z_int);  
                [vx,vy,vz] = nsimplex_3vec(vx_int,vy_int,vz_int,x_int,y_int,z_int);
                rho = nsimplex_3(rho_int,x_int,y_int,z_int);
                %ion_rho = nsimplex_3(ion_dens,x_int,y_int,z_int);
                %pres = nsimplex_3(p_comp,x_int,y_int,z_int);
                %temp = nsimplex_3(temp_comp,x_int,y_int,z_int);
                [vax,vay,vaz] = nsimplex_3vec([bx_int(:,1)./sqrt(rho_int(:,1)),bx_int(:,2)./sqrt(rho_int(:,2)),bx_int(:,3)./sqrt(rho_int(:,3))],...
                    [by_int(:,1)./sqrt(rho_int(:,1)),by_int(:,2)./sqrt(rho_int(:,2)),by_int(:,3)./sqrt(rho_int(:,3))],...
                    [bz_int(:,1)./sqrt(rho_int(:,1)),bz_int(:,2)./sqrt(rho_int(:,2)),bz_int(:,3)./sqrt(rho_int(:,3))],x_int,y_int,z_int);
%                 jx = jxs;
%                 jy = jys;
%                 jz = jzs;
            elseif sc == 4 && fpi_sc == 3
                [px,py,pz] = nsimplex_4vec(bx_int,by_int,bz_int,x_int,y_int,z_int);
                [vx,vy,vz] = nsimplex_3vec(vx_int,vy_int,vz_int,x_int,y_int,z_int);
                rho = nsimplex_3(rho_int,x_int,y_int,z_int);
                %ion_rho = nsimplex_3(ion_dens,x_int,y_int,z_int);
                %pres = nsimplex_3(p_comp,x_int,y_int,z_int);
                %temp = nsimplex_3(temp_comp,x_int,y_int,z_int);
                fpi = find([1 2 3 4] ~= missing);
                [vax,vay,vaz] = nsimplex_3vec([bx_int(:,fpi(1))./sqrt(rho_int(:,1)),bx_int(:,fpi(2))./sqrt(rho_int(:,2)),bx_int(:,fpi(3))./sqrt(rho_int(:,3))],...
                    [by_int(:,fpi(1))./sqrt(rho_int(:,1)),by_int(:,fpi(2))./sqrt(rho_int(:,2)),by_int(:,fpi(3))./sqrt(rho_int(:,3))],...
                    [bz_int(:,fpi(1))./sqrt(rho_int(:,1)),bz_int(:,fpi(2))./sqrt(rho_int(:,2)),bz_int(:,fpi(3))./sqrt(rho_int(:,3))],x_int,y_int,z_int);
%                 subplot(3,2,2); hold on;
%                 plot(jt,jxm,'r'); plot(jt,jym,'g'); plot(jt,jzm,'b');
%                 legend('j_x','j_y','j_z')
%                 plot_help('current density curlometer',jt(1),jt(end));
%                 jx = jxm;
%                 jy = jym;
%                 jz = jzm;
            elseif sc == 4 && fpi_sc == 2
                [px,py,pz] = nsimplex_4vec(bx_int,by_int,bz_int,x_int,y_int,z_int);
                vx = nanmean(vx_int,2);
                vy = nanmean(vy_int,2);
                vz = nanmean(vz_int,2);
                rho = nanmean(rho_int,2);
                %ion_rho = nanmean(ion_dens,2);
                %pres = nanmean(p_comp,2);
                %temp = nanmean(temp_comp,2);
                fpi = find([1 2 3 4] ~= missing);
                vax = nanmean([bx_int(:,fpi(1))./sqrt(rho_int(:,1)),bx_int(:,fpi(2))./sqrt(rho_int(:,2))],2);
                vay = nanmean([by_int(:,fpi(1))./sqrt(rho_int(:,1)),by_int(:,fpi(2))./sqrt(rho_int(:,2))],2);
                vaz = nanmean([bz_int(:,fpi(1))./sqrt(rho_int(:,1)),bz_int(:,fpi(2))./sqrt(rho_int(:,2))],2);
            end

            [V_DHT] = get_DHT_frame(px,py,pz,vx,vy,vz);
            vx = vx - V_DHT(1);
            vy = vy - V_DHT(2);
            vz = vz - V_DHT(3);

            t = jt;

            %figure('Visible','off','rend','painters','pos',[10 10 1000 800])
            %figure('rend','painters','pos',[10 10 1000 800])
            figure
            subplot(3,1,1); hold on;
            plot_all_mms_B(t,bx_int,by_int,bz_int,sc)
            lg = legend('B_x','B_y','B_z','Location','northeastoutside');
            plot_help([datestr(todatenum(cdfepoch(t(1)))),' to ',...
                datestr(todatenum(cdfepoch(t(end))))],t(1),t(end))
            ylabel('nT')
            subplot(3,1,2); hold on
            plot(t,vx-nanmean(vx),'r')
            plot(t,vy-nanmean(vy),'g','HandleVisibility','off')
            plot(t,vz-nanmean(vz),'b','HandleVisibility','off')
            plot(t,(1/1000)*(10^-9)*(vax-mean(vax))/sqrt(mu*(1.67*10^(-27))),'r--')
            plot(t,(1/1000)*(10^-9)*(vay-mean(vay))/sqrt(mu*(1.67*10^(-27))),'g--','HandleVisibility','off')
            plot(t,(1/1000)*(10^-9)*(vaz-mean(vaz))/sqrt(mu*(1.67*10^(-27))),'b--','HandleVisibility','off')
            lg = legend('v_x','v_{Ax}','Location','northeastoutside');
            ylabel('km/s')
            plot_help('v and v_A',t(1),t(end))
            subplot(3,1,3); hold on
            plot(t,rho,'k')
            plot_help('electron density',t(1),t(end))
            legend('\rho^-','Location','northeastoutside')
            ylabel('m^{-3}')
            subplot(3,1,1)  
cs = 387;
            
%             maggg = sqrt(px.^2+py.^2+pz.^2);
%             dirb = [px./maggg,py./maggg,pz./maggg];
%     % 
%     %             subplot(3,2,5); hold on;
%     %             plot(jt,jxs,'r'); plot(jt,jys,'g'); plot(jt,jzs,'b');
%     %             legend('j_x','j_y','j_z')
%     %             plot_help('current density single spacecraft',jt(1),jt(end));
% 
%     %             subplot(3,2,4); hold on;
%     %             plot(jt,sum([jx,jy,jz].*dirb,2),'r')
%     %             jxb = cross([jx,jy,jz],dirb);
%     %             plot(jt,sqrt(jxb(:,1).^2+jxb(:,2).^2+jxb(:,3).^2),'g');
%     %             legend('j_{||}','j_\perp')
%     %             plot_help('j magnetic field aligned coordinates',jt(1),jt(end));
 
%                 dir_cs = ones(size(dirb(1:length(jt),:)));
%                 dir_cs(:,1) = dir_cs(:,1)*mva_dirx;
%                 dir_cs(:,2) = dir_cs(:,2)*mva_diry;
%                 dir_cs(:,3) = dir_cs(:,3)*mva_dirz;
%     %             subplot(3,2,6); hold on;
%     %             plot(jt,sum([jx,jy,jz].*dir_cs,2),'r')
%     %             jxb = cross([jx,jy,jz],dir_cs);
%     %             plot(jt,sqrt(jxb(:,1).^2+jxb(:,2).^2+jxb(:,3).^2),'g');
%     %             plot_help('j boundary normal coordinates',jt(1),jt(end));
%     %             legend('j_{||}','j_\perp')  

            for j = 1:length(csc)
                cs = cs+1;
                tags(cs) = jjj;

                ind1 = floor(cs_inds(cs,1));
                ind2 = floor(cs_inds(cs,2));

                dis_range = ind1:ind2;
                %t1 = jt(csc(j)-wl(j)); t2 = jt(csc(j)+wl(j));
                %dis_range = find(abs(t-t1) == min(abs(t-t1))):find(abs(t-t2) == min(abs(t-t2)));
                jti = t(dis_range);

                vxi = vx(dis_range);
                vyi = vy(dis_range);
                vzi = vz(dis_range);

                cs_vx(cs) = nanmean(vxi);
                cs_vy(cs) = nanmean(vyi);
                cs_vz(cs) = nanmean(vzi);

                pxi = px(dis_range);
                pyi = py(dis_range);
                pzi = pz(dis_range);
                %presi = pres(dis_range);

                vaxi = (1/1000)*(10^-9)*(vax(dis_range))./sqrt(mu*(1.67*10^(-27)));
                vayi = (1/1000)*(10^-9)*(vay(dis_range))./sqrt(mu*(1.67*10^(-27)));
                vazi = (1/1000)*(10^-9)*(vaz(dis_range))./sqrt(mu*(1.67*10^(-27)));


                [cs_walen1(cs),v_jumps1(cs,:),va_jumps1(cs,:),walegood1(cs),w_ind11,w_ind22] = walen(jti,vxi,vyi,vzi,vaxi,vayi,vazi,...
                    pxi,pyi,pzi,0.75,1.25,0);
                %if w_ind11 > 0
                %    plot(jti(w_ind11:w_ind22),pxi(w_ind11:w_ind22),'y','LineWidth',3)
                %end

                [cs_walen2(cs),v_jumps2(cs,:),va_jumps2(cs,:),walegood2(cs),w_ind111,w_ind222] = walen(jti,vxi,vyi,vzi,vaxi,vayi,vazi,...
                    pxi,pyi,pzi,0.9,1.1,0);
                %if w_ind111 > 0
                %    plot(jti(w_ind111:w_ind222),pxi(w_ind111:w_ind222),'c','LineWidth',3)
                %end
try
fill([jti(w_ind11) jti(w_ind22) jti(w_ind22) jti(w_ind11)],[-10 -10  5 5],'y')
fill([jti(w_ind111) jti(w_ind222) jti(w_ind222) jti(w_ind111)],[-10 -10  5 5],'c')

catch
end

                plot(jti,pxi,'k','LineWidth',3)
                plot(jti,pyi,'k','LineWidth',3)
                plot(jti,pzi,'k','LineWidth',3)

%                 er = 0;
%                 wl_mva = 0;
%                 if length(jt) > 10000
%                     wind_lengg = 1000:100:4200;
%                 else
%                     wind_lengg = 450:50:2000;
%                 end
%                 for wind_leng = wind_lengg
%                     if ((csc(j)-wind_leng/2) > 0 && (csc(j)+wind_leng/2 < length(px)))
%                         bxs = px(csc(j)-wind_leng/2:csc(j)+wind_leng/2);
%                         bys = px(csc(j)-wind_leng/2:csc(j)+wind_leng/2);
%                         bzs = px(csc(j)-wind_leng/2:csc(j)+wind_leng/2);
%                         [ndv1,ndv2,eigrat] = normal_dir_var(bxs,bys,bzs);
%                         if eigrat > er
%                             er = eigrat;
%                             n1 = ndv1;
%                             n2 = ndv2;
%                             wl_mva = wind_leng;
%                             if er > 10
%                                 break
%                             end
%                         end
%                     end
%                 end
                %scatter([jt(csc(j)-wl_mva/2) jt(csc(j)+wl_mva/2)],[0 0],'rp')

%                 figure; hold on
%                 plot((10^-12)*vax./sqrt(mu*(1.67*10^(-27))),'r')
%                 plot((10^-12)*vay./sqrt(mu*(1.67*10^(-27))),'g')
%                 plot((10^-12)*vaz./sqrt(mu*(1.67*10^(-27))),'b')
%                 plot(vx-mean(vx),'r--')   
%                 plot(vy-mean(vy),'g--')
%                 plot(vz-mean(vz),'b--')
%                 [cs_x,cs_y] = ginput(2);
%                 close 2
                %cs_inds(cs,:) = cs_x;

%                 if w_ind1 == 0
%                     dir_cs1 = ones(length(csc(j)-wl(j)/2:csc(j)+wl(j)/2),3); 
%                     dir_cs1(1,:) = [n1(1),n1(2),n1(3)];
%                     dir_cs1(2:end,1) = dir_cs1(2:end,1)*dir_cs1(1,1);
%                     dir_cs1(2:end,2) = dir_cs1(2:end,2)*dir_cs1(1,2);
%                     dir_cs1(2:end,3) = dir_cs1(2:end,3)*dir_cs1(1,3);
% 
%                     dir_cs2 = ones(length(csc(j)-wl(j)/2:csc(j)+wl(j)/2),3);
%                     dir_cs2(1,:) = [n2(1),n2(2),n2(3)];
%                     dir_cs2(2:end,1) = dir_cs2(2:end,1)*dir_cs2(1,1);
%                     dir_cs2(2:end,2) = dir_cs2(2:end,2)*dir_cs2(1,2);
%                     dir_cs2(2:end,3) = dir_cs2(2:end,3)*dir_cs2(1,3);
% 
%                     bpar1 = mean(abs(sum([px((csc(j)-wl(j)/2):(csc(j)+wl(j)/2)),...
%                         py((csc(j)-wl(j)/2):(csc(j)+wl(j)/2)),pz((csc(j)-wl(j)/2):(csc(j)+wl(j)/2))].*dir_cs1,2)));
%                     bpar2 =  mean(abs(sum([px((csc(j)-wl(j)/2):(csc(j)+wl(j)/2)),...
%                         py((csc(j)-wl(j)/2):(csc(j)+wl(j)/2)),pz((csc(j)-wl(j)/2):(csc(j)+wl(j)/2))].*dir_cs2,2)));
%                     if (bpar1 > bpar2)
%                         dir_cs1 = dir_cs2;
%                     end
%                     bxd = cross([px((csc(j)-wl(j)/2):(csc(j)+wl(j)/2)),py((csc(j)-wl(j)/2):(csc(j)+wl(j)/2)),...
%                        pz((csc(j)-wl(j)/2):(csc(j)+wl(j)/2))],dir_cs1);
%                else              
%                     dir_cs1 = ones(length(w_ind1:w_ind2),3); 
%                     dir_cs1(1,:) = [n1(1),n1(2),n1(3)];
%                     dir_cs1(2:end,1) = dir_cs1(2:end,1)*dir_cs1(1,1);
%                     dir_cs1(2:end,2) = dir_cs1(2:end,2)*dir_cs1(1,2);
%                     dir_cs1(2:end,3) = dir_cs1(2:end,3)*dir_cs1(1,3);
% 
%                     dir_cs2 = ones(length(w_ind1:w_ind2),3); 
%                     dir_cs2(1,:) = [n2(1),n2(2),n2(3)];
%                     dir_cs2(2:end,1) = dir_cs2(2:end,1)*dir_cs2(1,1);
%                     dir_cs2(2:end,2) = dir_cs2(2:end,2)*dir_cs2(1,2);
%                     dir_cs2(2:end,3) = dir_cs2(2:end,3)*dir_cs2(1,3);
% 
%                     bpar1 = mean(abs(sum([pxi(w_ind1:w_ind2),pyi(w_ind1:w_ind2),pzi(w_ind1:w_ind2)].*dir_cs1,2)));
%                     bpar2 =  mean(abs(sum([pxi(w_ind1:w_ind2),pyi(w_ind1:w_ind2),pzi(w_ind1:w_ind2)].*dir_cs2,2)));
%                     if (bpar1 > bpar2)
%                         dir_cs1 = dir_cs2;
%                     end
%                     bxd = cross([pxi(w_ind1:w_ind2),pyi(w_ind1:w_ind2),pzi(w_ind1:w_ind2)],dir_cs1);
%                end

%                 cs_Bpar(cs) = min([bpar1,bpar2]);
%                 cs_Bperp(cs) = mean(sqrt(bxd(:,1).^2+bxd(:,2).^2+bxd(:,3).^2));
%                 cs_eigrats(cs) = er;

%                 if w_ind1 == 0
%                     %2*mu_0*p / B^2 with p in nPa and B in nT thus only 1
%                     %factor of 10^-9 on the bottom
%                     ion_beta(cs) = 2*mu*mean(pres((csc(j)-wl(j)/2):(csc(j)+wl(j)/2)))/((cs_Bpar(cs)^2+cs_Bperp(cs)^2)*(10^-9));
%                 else
%                     ion_beta(cs) = 2*mu*mean(presi(w_ind1:w_ind2))/((cs_Bpar(cs)^2+cs_Bperp(cs)^2)*(10^-9));
%                 end 
                cs_jumps_tot(cs) = sqrt((px(ind1)-px(ind2))^2+(py(ind1)-py(ind2))^2+(pz(ind1)-pz(ind2))^2); 
                v_jumps_totx(cs) =  sqrt((vx(ind1)-vx(ind2))^2); 
                v_jumps_toty(cs) =  sqrt((vy(ind1)-vy(ind2))^2); 
                v_jumps_totz(cs) =  sqrt((vz(ind1)-vz(ind2))^2); 
                va_jumps_totx(cs) =  sqrt((vax(ind1)-vax(ind2))^2);  
                va_jumps_toty(cs) =  sqrt((vay(ind1)-vay(ind2))^2); 
                va_jumps_totz(cs) =  sqrt((vaz(ind1)-vaz(ind2))^2); 
                cs_loc(cs,1) = x_int(1); cs_loc(cs,2) = y_int(1); cs_loc(cs,3) = z_int(1);
            end
           saveas(1,['mag_cadence/discont_',num2str(jjj),'.png'])
           %close 1
        end
    end
%     save('mcs_walen1','cs_walen1')
%     save('mcs_walen2','cs_walen2')
%     save('mwalegood1','walegood1')
%     save('mwalegood2','walegood2')
%     save('mcs_vx','cs_vx')
%     save('mcs_vy','cs_vy')
%     save('mcs_vz','cs_vz')
%     save('mva_jumps_totx','va_jumps_totx')
%     save('mva_jumps_toty','va_jumps_toty')
%     save('mva_jumps_totz','va_jumps_totz')
%     save('mv_jumps_totx','v_jumps_totx')
%     save('mv_jumps_toty','v_jumps_toty')
%     save('mv_jumps_totz','v_jumps_totz')
%     save('mv_jumps1','v_jumps1')
%     save('mv_jumps2','v_jumps2')
%     save('mva_jumps1','va_jumps1')
%     save('mva_jumps2','va_jumps2')
%     save('mcs_jumps_tot','cs_jumps_tot')
%     save('mtags','tags')
%     save('mcs_loc','cs_loc')    
