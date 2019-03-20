%this file contains strings of all the mms data files
fid = fopen('mms_datas.txt');
A = fscanf(fid,'%s\n');

%current_sheets(303) = struct();
load('current_sheets.mat')

for jjj = [1:200,202:224,226:303]
    jjj
    %determine the path of the jjj-th mms data file
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

    %this function gets B and GSE position for all available spacecraft
    %[curlometer,single spacecraft,bx,by,bz,x,y,z, time vector, thing for path, number available spacecraft]
    [j1,j2,bx_int1,by_int1,bz_int1,x_int1,y_int1,z_int1,jt,B_tag,sc] = cur_dens(path1,B_path1,B_path2);

    B_path2 =  [B_path2,'_v5.',num2str(B_tag),'.0.cdf'];
    missing = 0;

    %not sure about lengths of vectors and #sc so cant preallocate :(
    vx_comp = []; vy_comp = []; vz_comp = []; rho_comp = []; p_comp = []; temp_comp = []; etemp_comp = [];
    ion_dens = []; ti = []; te = [];

    %getting FPI data
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

            ti(1:length(data{1}),fpi_sc) = data{1};
            te(1:length(electro_data{1}),fpi_sc) = electro_data{1};
            numberdensity = electro_data{23}; %cm^-3
            bulkv = data{25}; %GSE

            etemptens = electro_data{36};
            temptens = data{33}; 
            teemp = zeros(length(bulkv(:,1)),1);
            eteemp = zeros(length(numberdensity),1);
             for iii = 1:length(bulkv(:,1))
                 teemp(iii) = trace(temptens(:,:,iii));
             end
             for iii = 1:length(numberdensity)
                 eteemp(iii) = trace(etemptens(:,:,iii));
             end

            wind = 20+ceil(length(bulkv(:,1))/50);
            vx_comp(1:length(bulkv(:,1)),fpi_sc) = smooth(hamperl(fft(bulkv(:,1)),wind),12);
            vy_comp(1:length(bulkv(:,2)),fpi_sc) = smooth(hamperl(fft(bulkv(:,2)),wind),24);
            vz_comp(1:length(bulkv(:,3)),fpi_sc) = smooth(hamperl(fft(bulkv(:,3)),wind),12);
            rho_comp(1:length(electro_data{23}),fpi_sc) = (10^6)*smooth(numberdensity,30);
            temp_comp(1:length(teemp),fpi_sc) = smooth(hamperl(fft(teemp),wind),12);
            etemp_comp(1:length(eteemp),fpi_sc) = smooth((eteemp),30);

        catch
            fprintf(['no MMS',num2str(i), 'FPI\n']);
            missing = i;
        end
    end  

    if fpi_sc > 0
%     vx  = bulkv(:,1);
%     vx_clean = vx_comp(:,end);
%     vy  = bulkv(:,2);
%     vy_clean = vy_comp(:,end);
%     vz  = bulkv(:,3);
%     vz_clean = vz_comp(:,end);
%     temp  = teemp;
%     temp_clean = temp_comp(:,end);

%     figure('Visible','off');
%     subplot(2,2,1); plot(vx); hold on; plot(vx_clean(1:end-1)); title('v_x')
%     subplot(2,2,2); plot(vy); hold on; plot(vy_clean(1:end-1)); title('v_y')
%     subplot(2,2,3);  plot(vz); hold on; plot(vz_clean(1:end-1)); title('v_z')
%     subplot(2,2,4); plot(temp); hold on; plot(temp_clean); title('temperature')
%     saveas(gcf,['cleanage',num2str(jjj),'.png'])

    B = repmat(bx_int1,1,1,3);
    B(:,:,2) = by_int1; B(:,:,3) = bz_int1;

    V = repmat(vx_comp,1,1,3);
    V(:,:,2) = vy_comp; V(:,:,3) = vz_comp;

    pos = repmat(x_int1,1,1,3);
    pos(:,:,2) = y_int1; pos(:,:,3) = z_int1;

    current_sheets(jjj).mag_data = B;
    current_sheets(jjj).v_data = V;
    current_sheets(jjj).rho_data = rho_comp;
    current_sheets(jjj).temp_data = temp_comp;
    current_sheets(jjj).etemp_data = etemp_comp;
    current_sheets(jjj).pos = pos;
    current_sheets(jjj).ti = ti;
    current_sheets(jjj).te = te;
    current_sheets(jjj).jt = jt;
    current_sheets(jjj).fpi_sc = fpi_sc;
    current_sheets(jjj).sc = sc;
    current_sheets(jjj).m = missing;
    end
end
save('current_sheets','current_sheets')
