function [] = MMS_fpi()
    mu = 4*pi*10^(-7);

    fid = fopen('mms_datas.txt');
    A = fscanf(fid,'%s\n');
    cs_Bperp = zeros(503,1);
    cs_Bpar = zeros(503,1);
    cs_walen1 = zeros(503,1);
    cs_walen2 = zeros(503,1);
    cs_walen3 = zeros(503,1);
    walegood1 = zeros(503,1);
    walegood2 = zeros(503,1);
    walegood3 = zeros(503,1);
    cs_eigrats = zeros(503,1);
    cs_vx = zeros(503,1);
    v_jumps_totx = zeros(503,1);
    v_jumps_toty = zeros(503,1);
    v_jumps_totz = zeros(503,1);
    va_jumps_totx = zeros(503,1);
    va_jumps_toty = zeros(503,1);
    va_jumps_totz = zeros(503,1);
    cs_jumps_tot = zeros(503,1);
    cs_jumps1 = zeros(503,1);
    cs_jumps2 = zeros(503,1);
    cs_jumps3 = zeros(503,1);
    cs_shear_tot = zeros(503,1);
    cs_shear1 = zeros(503,1);
    cs_shear2 = zeros(503,1);
    cs_shear3 = zeros(503,1);
    cs_shflow_stab = zeros(503,1);
    %ion_inlength = zeros(503,1);
    ion_beta = zeros(503,1);
    tags = zeros(503,1);
    cs = 0;
    skipthis = 0;

    %201 & 225 missing
   for jjj = [1:200,202:224,226:303]
   % for jjj = 1
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

        [j1,j2,bx_int,by_int,bz_int,x_int,y_int,z_int,jt,B_tag,sc] = cur_dens(path1,B_path1,B_path2);
        [wl,csc] = find_cs(jt,bx_int(:,1),by_int(:,1),bz_int(:,1),jjj);

%         if sc == 4
%             jxm = j1(:,1); jym = j1(:,2); jzm = j1(:,3);
%         end
%         jxs = j2(:,1); jys = j2(:,2); jzs = j2(:,3);
% 
        B_path2 =  [B_path2,'_v5.',num2str(B_tag),'.0.cdf'];
        missing = 0;
        vx_comp = []; vy_comp = []; vz_comp = []; rho_comp = []; p_comp = []; temp_comp = []; temperr_comp = [];
        ion_dens = [];

        fpi_sc = 0;
        for i = 1:4
            try
                bx = bx_int(:,i);
            catch
                fprintf(['only ',num2str(sc), ' MMS FGM files available\n']);
            end
            try
                [electro_data,~] = spdfcdfread([path1,num2str(i),e_path1,num2str(i),e_path2]);
                [data,~] = spdfcdfread([path1,num2str(i),i_path1,num2str(i),i_path2]);
                fpi_sc = fpi_sc + 1;
                 
                t = data{1};
                te = electro_data{1};
                ions = data{19};
                numberdensity = electro_data{23}; %cm^-3
                bulkv = data{25}; %GSE
                prestens = data{29}; %nPa
                temptens = data{33}; %

                prees = zeros(length(t),1);
                teemp = zeros(length(t),1);
                for iii = 1:length(t)
                    prees(iii) = trace(prestens(:,:,iii));
                    teemp(iii) = trace(temptens(:,:,iii));
                end

                wind = 20+ceil(length(bulkv(:,1))/50);
                vx_comp(1:length(jt),fpi_sc) = smooth(interp1(t,hamperl(fft(bulkv(:,1)),wind),jt),128);
                vy_comp(1:length(jt),fpi_sc) = smooth(interp1(t,hamperl(fft(bulkv(:,2)),wind),jt),128);
                vz_comp(1:length(jt),fpi_sc) = smooth(interp1(t,hamperl(fft(bulkv(:,3)),wind),jt),128);
                rho_comp(1:length(jt),fpi_sc) = (10^6)*smooth(interp1(te,numberdensity,jt),128);
                ion_dens(1:length(jt),fpi_sc) = (10^6)*smooth(interp1(t,hamperl(fft(ions),wind),jt),128);
                p_comp(1:length(jt),fpi_sc) = smooth(interp1(t,hamperl(fft(prees),wind),jt),128);
                temp_comp(1:length(jt),fpi_sc) = smooth(interp1(t,hamperl(fft(teemp),wind),jt),128);

            catch
                fprintf(['no MMS',num2str(i), 'FPI\n']);
                missing = i;
            end
        end
 
         if fpi_sc > 1
% %             jx = [];
% %             jy = [];
% %             jz = [];
            if sc == 4 && fpi_sc == 4
                [px,py,pz] = nsimplex_4vec(bx_int,by_int,bz_int,x_int,y_int,z_int);
                [vx,vy,vz] = nsimplex_4vec(vx_comp,vy_comp,vz_comp,x_int,y_int,z_int);
                rho = nsimplex_4(rho_comp,x_int,y_int,z_int);
                ion_rho = nsimplex_4(ion_dens,x_int,y_int,z_int);
                pres = nsimplex_4(p_comp,x_int,y_int,z_int);
                temp = nsimplex_4(temp_comp,x_int,y_int,z_int);
                [vax,vay,vaz] = nsimplex_4vec([bx_int(:,1)./sqrt(rho_comp(:,1)),bx_int(:,2)./sqrt(rho_comp(:,2)),bx_int(:,3)./sqrt(rho_comp(:,3)),bx_int(:,4)./sqrt(rho_comp(:,4))],...
                    [by_int(:,1)./sqrt(rho_comp(:,1)),by_int(:,2)./sqrt(rho_comp(:,2)),by_int(:,3)./sqrt(rho_comp(:,3)),by_int(:,4)./sqrt(rho_comp(:,4))],...
                    [bz_int(:,1)./sqrt(rho_comp(:,1)),bz_int(:,2)./sqrt(rho_comp(:,2)),bz_int(:,3)./sqrt(rho_comp(:,3)),bz_int(:,4)./sqrt(rho_comp(:,4))],x_int,y_int,z_int);
%                 subplot(3,2,2); hold on;
%                 plot(jt,jxm,'r'); plot(jt,jym,'g'); plot(jt,jzm,'b');
%                 legend('jx','jy','jz')
%                 plot_help('current density curlometer',jt(1),jt(end));
%                 jx = jxm;
%                 jy = jym;
%                 jz = jzm;
            elseif sc == 3 && fpi_sc == 3
                [px,py,pz] = nsimplex_3vec(bx_int,by_int,bz_int,x_int,y_int,z_int);  
                [vx,vy,vz] = nsimplex_3vec(vx_comp,vy_comp,vz_comp,x_int,y_int,z_int);
                rho = nsimplex_3(rho_comp,x_int,y_int,z_int);
                ion_rho = nsimplex_3(ion_dens,x_int,y_int,z_int);
                pres = nsimplex_3(p_comp,x_int,y_int,z_int);
                temp = nsimplex_3(temp_comp,x_int,y_int,z_int);
                [vax,vay,vaz] = nsimplex_3vec([bx_int(:,1)./sqrt(rho_comp(:,1)),bx_int(:,2)./sqrt(rho_comp(:,2)),bx_int(:,3)./sqrt(rho_comp(:,3))],...
                    [by_int(:,1)./sqrt(rho_comp(:,1)),by_int(:,2)./sqrt(rho_comp(:,2)),by_int(:,3)./sqrt(rho_comp(:,3))],...
                    [bz_int(:,1)./sqrt(rho_comp(:,1)),bz_int(:,2)./sqrt(rho_comp(:,2)),bz_int(:,3)./sqrt(rho_comp(:,3))],x_int,y_int,z_int);
%                 jx = jxs;
%                 jy = jys;
%                 jz = jzs;
            elseif sc == 4 && fpi_sc == 3
                [px,py,pz] = nsimplex_4vec(bx_int,by_int,bz_int,x_int,y_int,z_int);
                [vx,vy,vz] = nsimplex_3vec(vx_comp,vy_comp,vz_comp,x_int,y_int,z_int);
                rho = nsimplex_3(rho_comp,x_int,y_int,z_int);
                ion_rho = nsimplex_3(ion_dens,x_int,y_int,z_int);
                pres = nsimplex_3(p_comp,x_int,y_int,z_int);
                temp = nsimplex_3(temp_comp,x_int,y_int,z_int);
                fpi = find([1 2 3 4] ~= missing);
                [vax,vay,vaz] = nsimplex_3vec([bx_int(:,fpi(1))./sqrt(rho_comp(:,1)),bx_int(:,fpi(2))./sqrt(rho_comp(:,2)),bx_int(:,fpi(3))./sqrt(rho_comp(:,3))],...
                    [by_int(:,fpi(1))./sqrt(rho_comp(:,1)),by_int(:,fpi(2))./sqrt(rho_comp(:,2)),by_int(:,fpi(3))./sqrt(rho_comp(:,3))],...
                    [bz_int(:,fpi(1))./sqrt(rho_comp(:,1)),bz_int(:,fpi(2))./sqrt(rho_comp(:,2)),bz_int(:,fpi(3))./sqrt(rho_comp(:,3))],x_int,y_int,z_int);
%                 subplot(3,2,2); hold on;
%                 plot(jt,jxm,'r'); plot(jt,jym,'g'); plot(jt,jzm,'b');
%                 legend('j_x','j_y','j_z')
%                 plot_help('current density curlometer',jt(1),jt(end));
%                 jx = jxm;
%                 jy = jym;
%                 jz = jzm;
            elseif sc == 4 && fpi_sc == 2
                [px,py,pz] = nsimplex_4vec(bx_int,by_int,bz_int,x_int,y_int,z_int);
                vx = nanmean(vx_comp,2);
                vy = nanmean(vy_comp,2);
                vz = nanmean(vz_comp,2);
                rho = nanmean(rho_comp,2);
                ion_rho = nanmean(ion_dens,2);
                pres = nanmean(p_comp,2);
                temp = nanmean(temp_comp,2);
                fpi = find([1 2 3 4] ~= missing);
                vax = nanmean([bx_int(:,fpi(1))./sqrt(rho_comp(:,1)),bx_int(:,fpi(2))./sqrt(rho_comp(:,2))],2);
                vay = nanmean([by_int(:,fpi(1))./sqrt(rho_comp(:,1)),by_int(:,fpi(2))./sqrt(rho_comp(:,2))],2);
                vaz = nanmean([bz_int(:,fpi(1))./sqrt(rho_comp(:,1)),bz_int(:,fpi(2))./sqrt(rho_comp(:,2))],2);
            end
% 
            figure('Visible','off','rend','painters','pos',[10 10 1000 800])
            %figure
            subplot(4,1,1); hold on;
            plot_all_mms_B(jt,bx_int,by_int,bz_int,sc)
            %lg = legend('B_x','B_y','B_z','Location','northeast');
            %set(lg,'color','none');
            %set(lg, 'Box', 'off');
            plot_help(['B - num cs : ',num2str(length(csc))],jt(1),jt(end))
            %plot_help('B',jt(1),jt(end))
            subplot(4,1,2); hold on
            plot(jt,vx-nanmean(vx),'r')
            plot(jt,vy-nanmean(vy),'g','HandleVisibility','off')
            plot(jt,vz-nanmean(vz),'b','HandleVisibility','off')
            plot(jt,(1/1000)*(10^-9)*vax/sqrt(mu*(1.67*10^(-27))),'r--')
            plot(jt,(1/1000)*(10^-9)*vay/sqrt(mu*(1.67*10^(-27))),'g--','HandleVisibility','off')
            plot(jt,(1/1000)*(10^-9)*vaz/sqrt(mu*(1.67*10^(-27))),'b--','HandleVisibility','off')
            %lg = legend('v_x','v_{Ax}','Location','northeast');
            %set(lg,'color','none');
            %set(lg, 'Box', 'off');
            plot_help('v and v_A',jt(1),jt(end))
            subplot(4,1,3); hold on
            plot(jt,rho/max(rho),'k')
            plot(jt,ion_rho/max(ion_rho),'k--')
            plot(jt,pres/max(pres))
            plot(jt,temp/max(temp))
            plot_help_nice('normalized units',jt(1),jt(end))
            %lg = legend('n_e','n_i','p','T','Location','northeast');
            %set(lg,'color','none');
            %set(lg, 'Box', 'off');
            ylim([0 1])
            subplot(4,1,4)
            plot(jt,(10^-9)*pres+10^(-18)*(px.^2+py.^2+pz.^2)/(2*mu))
            plot_help('p+B^2/2\mu_0',jt(1),jt(end))
            subplot(4,1,1)  
            
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
                dis_range = csc(j)-wl(j):csc(j)+wl(j);
                jti = jt(dis_range);

                vxi = vx(dis_range);
                vyi = vy(dis_range);
                vzi = vz(dis_range);

                pxi = px(dis_range);
                pyi = py(dis_range);
                pzi = pz(dis_range);
                presi = pres(dis_range);

                vaxi = (1/1000)*(10^-9)*(vax(dis_range))./sqrt(mu*(1.67*10^(-27)));
                vayi = (1/1000)*(10^-9)*(vay(dis_range))./sqrt(mu*(1.67*10^(-27)));
                vazi = (1/1000)*(10^-9)*(vaz(dis_range))./sqrt(mu*(1.67*10^(-27)));

                scatter(jti,pxi,'k.')
                scatter(jti,pyi,'k.')
                scatter(jti,pzi,'k.')

                [cs_walen3(cs),cs_jumps3(cs),cs_shear3(cs),~,walegood3(cs),w_ind1,w_ind2] = walen(jti,vxi,vyi,vzi,vaxi,vayi,vazi,...
                    pxi,pyi,pzi,0.9,1.1);
                if w_ind1 > 0
                    scatter(jti(w_ind1:w_ind2),pxi(w_ind1:w_ind2),'m.')
                end
%                 for uuu = w_ind1:w_ind2
%                     h = plot([jti(uuu) jti(uuu)],[-10 10],'m');
%                     uistack(h,'bottom');
%                 end
                [cs_walen2(cs),cs_jumps2(cs),cs_shear2(cs),~,walegood2(cs),w_ind111,w_ind222] = walen(jti,vxi,vyi,vzi,vaxi,vayi,vazi,...
                    pxi,pyi,pzi,0.8,1.2);
                if w_ind111 > 0
                    scatter(jti(w_ind111:w_ind222),pxi(w_ind111:w_ind222),'c.')
                end
%                 for uuu = w_ind111:w_ind222
%                     h = plot([jti(uuu) jti(uuu)],[-10 10],'y');
%                     uistack(h,'bottom');
%                 end
                [cs_walen1(cs),cs_jumps1(cs),cs_shear1(cs),cs_vx(cs),walegood1(cs),w_ind11,w_ind22] = walen(jti,vxi,vyi,vzi,vaxi,vayi,vazi,...
                    pxi,pyi,pzi,0.75,1.25);
                if w_ind11 > 0
                    scatter(jti(w_ind11:w_ind22),pxi(w_ind11:w_ind22),'y.')
                end
%                 for uuu = w_ind11:w_ind22
%                     h = plot([jti(uuu) jti(uuu)],[-10 10],'c');
%                     uistack(h,'bottom');
%                 end

                er = 0;
                wl_mva = 0;
                if length(jt) > 10000
                    wind_lengg = 1000:100:4200;
                else
                    wind_lengg = 450:50:2000;
                end
                for wind_leng = wind_lengg
                    if ((csc(j)-wind_leng/2) > 0 && (csc(j)+wind_leng/2 < length(px)))
                        bxs = px(csc(j)-wind_leng/2:csc(j)+wind_leng/2);
                        bys = px(csc(j)-wind_leng/2:csc(j)+wind_leng/2);
                        bzs = px(csc(j)-wind_leng/2:csc(j)+wind_leng/2);
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
                scatter([jt(csc(j)-wl_mva/2) jt(csc(j)+wl_mva/2)],[0 0],'rp')

                if w_ind1 > 0
                    figure('Visible','off','rend','painters','pos',[10 10 900 700])
                    subplot(4,1,1); hold on;
                    plot_all_mms_B(jti,pxi,pyi,pzi,1)
                    scatter(jti(w_ind11:w_ind22),pxi(w_ind11:w_ind22),'c.')
                    scatter(jti(w_ind111:w_ind222),pxi(w_ind111:w_ind222),'y.')
                    scatter(jti(w_ind1:w_ind2),pxi(w_ind1:w_ind2),'m.')
                    lg = legend('B_x','B_y','B_z','Location','northeast');
                    set(lg,'color','none');
                    set(lg, 'Box', 'off');
                    plot_help(['B - tag : ',num2str(jjj)],jti(1),jti(end))
                    subplot(4,1,2); hold on
                    plot(jti,vxi-nanmean(vxi),'r')
                    plot(jti,vyi-nanmean(vyi),'g','HandleVisibility','off')
                    plot(jti,vzi-nanmean(vzi),'b','HandleVisibility','off')
                    plot(jti,(1/1000)*(10^-9)*vaxi/sqrt(mu*(1.67*10^(-27))),'r--')
                    plot(jti,(1/1000)*(10^-9)*vayi/sqrt(mu*(1.67*10^(-27))),'g--','HandleVisibility','off')
                    plot(jti,(1/1000)*(10^-9)*vazi/sqrt(mu*(1.67*10^(-27))),'b--','HandleVisibility','off')
                    lg = legend('v_x','v_{Ax}','Location','northeast');
                    set(lg,'color','none');
                    set(lg, 'Box', 'off');
                    plot_help('v and v_A',jti(1),jti(end))
                    subplot(4,1,3); hold on
                    plot(jti,rho(dis_range)/max(rho(dis_range)))
                    plot(jti,presi/max(presi))
                    plot(jti,temp(dis_range)/max(temp(dis_range)))
                    plot_help_nice('normalized units',jti(1),jti(end))
                    lg = legend('n','p','T','Location','northeast');
                    set(lg,'color','none');
                    set(lg, 'Box', 'off');
                    ylim([0 1])
                    subplot(4,1,4)
                    plot(jti,(10^-9)*presi+10^(-18)*(pxi.^2+pyi.^2+pzi.^2)/(2*mu))
                    plot_help('p+B^2/2\mu_0',jti(1),jti(end)) 
                    saveas(2,['cs_discont',num2str(cs),'.png'])
                    close 2
                end

                x1 = min(pxi);
                x2 = max(pxi);
                y1 = min(pyi);
                y2 = max(pyi);
                z1 = min(pzi);
                z2 = max(pzi);

                ind1x = find(pxi == x1,1);
                ind2x = find(pxi == x2,1);
                ind1y = find(pyi == y1,1);
                ind2y = find(pyi == y2,1);
                ind1z = find(pzi == z1,1);
                ind2z = find(pzi == z2,1);
                if abs(x1-x2) > abs(y1-y2) && abs(x1-x2) > abs(z1-z2)
                    ind1 = ind1x; ind2 = ind2x;
                elseif abs(y1-y2) > abs(x1-x2) && abs(y1-y2) > abs(z1-z2)
                    ind1 = ind1y; ind2 = ind2y;
                else
                    ind1 = ind1z; ind2 = ind2z;
                end

                if w_ind1 == 0
                    dir_cs1 = ones(length(csc(j)-wl(j)/2:csc(j)+wl(j)/2),3); 
                    dir_cs1(1,:) = [n1(1),n1(2),n1(3)];
                    dir_cs1(2:end,1) = dir_cs1(2:end,1)*dir_cs1(1,1);
                    dir_cs1(2:end,2) = dir_cs1(2:end,2)*dir_cs1(1,2);
                    dir_cs1(2:end,3) = dir_cs1(2:end,3)*dir_cs1(1,3);

                    dir_cs2 = ones(length(csc(j)-wl(j)/2:csc(j)+wl(j)/2),3);
                    dir_cs2(1,:) = [n2(1),n2(2),n2(3)];
                    dir_cs2(2:end,1) = dir_cs2(2:end,1)*dir_cs2(1,1);
                    dir_cs2(2:end,2) = dir_cs2(2:end,2)*dir_cs2(1,2);
                    dir_cs2(2:end,3) = dir_cs2(2:end,3)*dir_cs2(1,3);

                    bpar1 = mean(abs(sum([px((csc(j)-wl(j)/2):(csc(j)+wl(j)/2)),...
                        py((csc(j)-wl(j)/2):(csc(j)+wl(j)/2)),pz((csc(j)-wl(j)/2):(csc(j)+wl(j)/2))].*dir_cs1,2)));
                    bpar2 =  mean(abs(sum([px((csc(j)-wl(j)/2):(csc(j)+wl(j)/2)),...
                        py((csc(j)-wl(j)/2):(csc(j)+wl(j)/2)),pz((csc(j)-wl(j)/2):(csc(j)+wl(j)/2))].*dir_cs2,2)));
                    if (bpar1 > bpar2)
                        dir_cs1 = dir_cs2;
                    end
                    bxd = cross([px((csc(j)-wl(j)/2):(csc(j)+wl(j)/2)),py((csc(j)-wl(j)/2):(csc(j)+wl(j)/2)),...
                       pz((csc(j)-wl(j)/2):(csc(j)+wl(j)/2))],dir_cs1);

                    [~,~,rv1,rv2] = riemann_problemator(pxi(ind1),pxi(ind2),pyi(ind1),pyi(ind2),pzi(ind1),pzi(ind2),...
                        vxi(ind1)-mean(vxi),vyi(ind1)-mean(vyi),vzi(ind1)-mean(vzi),vxi(ind2)-mean(vxi),vyi(ind2)-mean(vyi),vzi(ind2)-mean(vzi));

                    cs_shflow_stab(cs) = abs(rv1(3) - rv2(3))/(sqrt(mean(vaxi(ind1:ind2).^2)+mean(vayi(ind1:ind2).^2)+...
                        mean(vazi(ind1:ind2).^2)));
 
                else              
                    dir_cs1 = ones(length(w_ind1:w_ind2),3); 
                    dir_cs1(1,:) = [n1(1),n1(2),n1(3)];
                    dir_cs1(2:end,1) = dir_cs1(2:end,1)*dir_cs1(1,1);
                    dir_cs1(2:end,2) = dir_cs1(2:end,2)*dir_cs1(1,2);
                    dir_cs1(2:end,3) = dir_cs1(2:end,3)*dir_cs1(1,3);

                    dir_cs2 = ones(length(w_ind1:w_ind2),3); 
                    dir_cs2(1,:) = [n2(1),n2(2),n2(3)];
                    dir_cs2(2:end,1) = dir_cs2(2:end,1)*dir_cs2(1,1);
                    dir_cs2(2:end,2) = dir_cs2(2:end,2)*dir_cs2(1,2);
                    dir_cs2(2:end,3) = dir_cs2(2:end,3)*dir_cs2(1,3);

                    bpar1 = mean(abs(sum([pxi(w_ind1:w_ind2),pyi(w_ind1:w_ind2),pzi(w_ind1:w_ind2)].*dir_cs1,2)));
                    bpar2 =  mean(abs(sum([pxi(w_ind1:w_ind2),pyi(w_ind1:w_ind2),pzi(w_ind1:w_ind2)].*dir_cs2,2)));
                    if (bpar1 > bpar2)
                        dir_cs1 = dir_cs2;
                    end
                    bxd = cross([pxi(w_ind1:w_ind2),pyi(w_ind1:w_ind2),pzi(w_ind1:w_ind2)],dir_cs1);

                    [~,~,rv1,rv2] = riemann_problemator(pxi(w_ind1),pxi(w_ind2),pyi(w_ind1),pyi(w_ind2),pzi(w_ind1),pzi(w_ind2),...
                        vxi(w_ind1)-mean(vxi),vyi(w_ind1)-mean(vyi),vzi(w_ind1)-mean(vzi),vxi(w_ind2)-mean(vxi),vyi(w_ind2)-mean(vyi),vzi(w_ind2)-mean(vzi));

                    cs_shflow_stab(cs) = abs(rv1(3) - rv2(3))/sqrt(mean(vaxi(w_ind1:w_ind2).^2)+mean(vayi(w_ind1:w_ind2).^2)+mean(vazi(w_ind1:w_ind2).^2));
                end

                cs_Bpar(cs) = min([bpar1,bpar2]);
                cs_Bperp(cs) = mean(sqrt(bxd(:,1).^2+bxd(:,2).^2+bxd(:,3).^2));
                cs_eigrats(cs) = er;

                if w_ind1 == 0
                    %2*mu_0*p / B^2 with p in nPa and B in nT thus only 1
                    %factor of 10^-9 on the bottom
                    ion_beta(cs) = 2*mu*mean(pres((csc(j)-wl(j)/2):(csc(j)+wl(j)/2)))/((cs_Bpar(cs)^2+cs_Bperp(cs)^2)*(10^-9));
                else
                    ion_beta(cs) = 2*mu*mean(presi(w_ind1:w_ind2))/((cs_Bpar(cs)^2+cs_Bperp(cs)^2)*(10^-9));
                end

                cs_jumps_tot(cs) = sqrt((pxi(ind1)-pxi(ind2))^2+(pyi(ind1)-pyi(ind2))^2+(pzi(ind1)-pzi(ind2))^2); 
                v_jumps_totx(cs) =  sqrt((vxi(ind1)-vxi(ind2))^2); 
                v_jumps_toty(cs) =  sqrt((vyi(ind1)-vyi(ind2))^2); 
                v_jumps_totz(cs) =  sqrt((vzi(ind1)-vzi(ind2))^2); 
                va_jumps_totx(cs) =  sqrt((vaxi(ind1)-vaxi(ind2))^2); 
                va_jumps_toty(cs) =  sqrt((vayi(ind1)-vayi(ind2))^2); 
                va_jumps_totz(cs) =  sqrt((vazi(ind1)-vazi(ind2))^2); 
                cs_shear_tot(cs) = (pxi(ind1)*pxi(ind2) + pyi(ind1)*pyi(ind2) +...
                    pzi(ind1)*pzi(ind2))/sqrt((pxi(ind1)^2 + pyi(ind1)^2 +...
                    pzi(ind1)^2)*(pxi(ind2)^2 + pyi(ind2)^2 + pzi(ind2)^2));
            end
          saveas(1,['discont_',num2str(jjj),'.png'])
          close 1
        end
    end
    save('cs_Bpar','cs_Bpar')
    save('cs_Bperp','cs_Bperp')
    save('cs_walen1','cs_walen1')
    save('cs_walen2','cs_walen2')
    save('cs_walen3','cs_walen3')
    save('walegood1','walegood1')
    save('walegood2','walegood2')
    save('walegood3','walegood3')
    save('cs_eigrats','cs_eigrats')
    save('cs_vx','cs_vx')
    save('va_jumps_totx','va_jumps_totx')
    save('va_jumps_toty','va_jumps_toty')
    save('va_jumps_totz','va_jumps_totz')
    save('v_jumps_totx','v_jumps_totx')
    save('v_jumps_toty','v_jumps_toty')
    save('v_jumps_totz','v_jumps_totz')
    save('cs_jumps1','cs_jumps1')
    save('cs_jumps2','cs_jumps2')
    save('cs_jumps3','cs_jumps3')
    save('cs_jumps_tot','cs_jumps_tot')
    save('cs_shear_tot','cs_shear_tot')
    save('cs_shear1','cs_shear1')
    save('cs_shear2','cs_shear2')
    save('cs_shear3','cs_shear3')
    save('cs_shflow_stab','cs_shflow_stab')
    save('ion_beta','ion_beta')
    save('tags','tags')
end

function [] = plot_help(tit,miin,maax)
    title(tit);
    ax1 = gca;
    ax1.XTick = [];
    xlim([miin,maax]);
end

function [] = plot_help_nice(tit,miin,maax)
    title(tit);
    ax2 = gca;
    ax2.XTick = [miin,maax];
    ax2.XTickLabel = datestr(todatenum([cdfepoch(miin),cdfepoch(maax)]));
    xlim([miin,maax]);
end

function [] = plot_all_mms_B(jt,bx_int,by_int,bz_int,sc)
    for i = 1:sc
        plot(jt,bx_int(:,i),'r')
        plot(jt,by_int(:,i),'g')
        plot(jt,bz_int(:,i),'b')
    end
end

  
    
