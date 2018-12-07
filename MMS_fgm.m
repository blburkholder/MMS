function [x_comp,y_comp,z_comp,t_pos,bx_comp,by_comp,bz_comp,t_B,B_tag,sc] = MMS_fgm(pathpath,B_path1,B_path2)
    sc = 0;
    span = 64;
    for i = 1:4
        for j = 1:300
            try
                B_path22 = [B_path2,'_v5.',num2str(j),'.0.cdf'];
                [data,iinfso] = spdfcdfread([pathpath,num2str(i),B_path1,num2str(i),B_path22]);
                B_tag = j;
                sc = sc + 1;

                B_time = data{1};
                B_GSE = data{2};
                %B_GSM = data{3};
                %B_DMPA = data{4};
                %B_BCS = data{5};

                pos_time = data{7};
                pos_GSE = data{8}; %position every 30 secs?
                %pos_GSM = data{9};
                
                bx1 = smooth(B_GSE(:,1),span);
                by1 = smooth(B_GSE(:,2),span);
                bz1 = smooth(B_GSE(:,3),span);

                bx_comp(1:length(B_GSE(:,1)),sc) = bx1;
                by_comp(1:length(B_GSE(:,2)),sc) = by1;
                bz_comp(1:length(B_GSE(:,3)),sc) = bz1;
                t_B(1:length(B_time),sc) = B_time;
                t_pos(1:length(pos_time),sc) = pos_time;
                x_comp(1:length(pos_GSE(:,1)),sc) = pos_GSE(:,1);    
                y_comp(1:length(pos_GSE(:,2)),sc) = pos_GSE(:,2);    
                z_comp(1:length(pos_GSE(:,3)),sc) = pos_GSE(:,3);
                break;
            catch MException

            end
        end
    end
end




