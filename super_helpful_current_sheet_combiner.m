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
    vx_comp = []; vy_comp = []; vz_comp = []; rho_comp = []; p_comp = []; temppar_comp = []; etemppar_comp = [];
    ion_dens = []; ti = []; te = []; tempperp_comp = []; etempperp_comp = []; p_comp = []; ep_comp = [];

    %getting FPI data
    fpi_sc = 0;
    for i = 1:4
        try
            bx = bx_int1(:,i);
        catch
            fprintf(['only ',num2str(sc), ' MMS FGM files available\n']);
        end
        try
            [electro_data,infos] = spdfcdfread([path1,num2str(i),e_path1,num2str(i),e_path2]);
            [data,infso] = spdfcdfread([path1,num2str(i),i_path1,num2str(i),i_path2]);
            fpi_sc = fpi_sc + 1;

            ti(1:length(data{1}),fpi_sc) = data{1};
            te(1:length(electro_data{1}),fpi_sc) = electro_data{1};
            numberdensity = electro_data{23}; %cm^-3
            densi = data{19};
            bulkv = data{25}; %GSE
            % pay attention to this error flag
            % https://lasp.colorado.edu/galaxy/pages/viewpage.action?pageId=37618954
            % erreflag = electro_data{4};
    
            etempperp = electro_data{46};
            etemppar = electro_data{45};
            tempperp = data{41};
            temppar = data{40}; 

            eep = electro_data{33};
            ep = zeros(length(electro_data{23}),1);
            for jk = 1:length(electro_data{23})
                ep(jk) = trace(eep(:,:,jk));
            end
            pp = data{29}; 
            p = zeros(length(pp),1);
            for jk = 1:length(pp)
                p(jk) = trace(pp(:,:,jk));
            end

            wind = 20+ceil(length(bulkv(:,1))/50);
            vx_comp(1:length(bulkv(:,1)),fpi_sc) = smooth(hamperl(fft(hamperl(fft(bulkv(:,1)),wind)),wind),6);
            vy_comp(1:length(bulkv(:,2)),fpi_sc) = smooth(hamperl(fft(hamperl(fft(bulkv(:,2)),wind)),wind),6);
            vz_comp(1:length(bulkv(:,3)),fpi_sc) = smooth(hamperl(fft(hamperl(fft(bulkv(:,3)),wind)),wind),6);
            rho_comp(1:length(electro_data{23}),fpi_sc) = (10^6)*smooth(numberdensity,30);

            temppar_comp(1:length(temppar),fpi_sc) = smooth(hamperl(fft(temppar),wind),12);
            etemppar_comp(1:length(etemppar),fpi_sc) = smooth((etemppar),30);
            tempperp_comp(1:length(tempperp),fpi_sc) = smooth(hamperl(fft(tempperp),wind),12);
            etempperp_comp(1:length(etempperp),fpi_sc) = smooth((etempperp),30);

            p_comp(1:length(densi),fpi_sc) = smooth(hamperl(fft(hamperl(fft(p),wind)),wind),12);
            ep_comp(1:length(electro_data{23}),fpi_sc) = smooth(ep,30);

%             p_comp(1:length(densi),fpi_sc) = smooth(hamperl(fft(densi.*(2*tempperp+temppar)/3),wind),12);
%             ep_comp(1:length(electro_data{23}),fpi_sc) = smooth(numberdensity.*(2*etempperp+etemppar)/3,30);
       catch
            fprintf(['no MMS',num2str(i), 'FPI\n']);
            missing = i;
       end
    end  

    if fpi_sc > 0
    
    B = repmat(bx_int1,1,1,3);
    B(:,:,2) = by_int1; B(:,:,3) = bz_int1;

    V = repmat(vx_comp,1,1,3);
    V(:,:,2) = vy_comp; V(:,:,3) = vz_comp;

    pos = repmat(x_int1,1,1,3);
    pos(:,:,2) = y_int1; pos(:,:,3) = z_int1;

    current_sheets(jjj).mag_data = B;
    current_sheets(jjj).v_data = V;
    current_sheets(jjj).rho_data = rho_comp;

    current_sheets(jjj).temp_data = temppar_comp;
    current_sheets(jjj).etemp_data = etemppar_comp;
    current_sheets(jjj).tempperp_data = tempperp_comp;
    current_sheets(jjj).etempperp_data = etempperp_comp;
    current_sheets(jjj).p_data = p_comp;
    current_sheets(jjj).ep_data = ep_comp;

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
