function [interpedx,interpedy,interpedz] = nsimplex_3vec(datx,daty,datz,x_comp_int,y_comp_int,z_comp_int)
        mx = mean(x_comp_int,2);
        my = mean(y_comp_int,2);
        mz = mean(z_comp_int,2);

        x12 = x_comp_int(:,1) - x_comp_int(:,2);
        x23 = x_comp_int(:,2) - x_comp_int(:,3);
        x13 = x_comp_int(:,1) - x_comp_int(:,3);
        y12 = y_comp_int(:,1) - y_comp_int(:,2);
        y23 = y_comp_int(:,2) - y_comp_int(:,3);
        y13 = y_comp_int(:,1) - y_comp_int(:,3);
        z12 = z_comp_int(:,1) - z_comp_int(:,2);
        z23 = z_comp_int(:,2) - z_comp_int(:,3);
        z13 = z_comp_int(:,1) - z_comp_int(:,3);

        l12 = sqrt(x12.^2+y12.^2+z12.^2);
        l23 = sqrt(x23.^2+y23.^2+z23.^2);
        l13 = sqrt(x13.^2+y13.^2+z13.^2);
        s = (l12+l23+l13)/2;
        A_abc = sqrt(s.*(s-l12).*(s-l23).*(s-l13));

        x1p = x_comp_int(:,1) - mx;
        x2p = x_comp_int(:,2) - mx;
        x3p = x_comp_int(:,3) - mx;

        y1p = y_comp_int(:,1) - my;
        y2p = y_comp_int(:,2) - my;
        y3p = y_comp_int(:,3) - my; 

        z1p = z_comp_int(:,1) - mz;
        z2p = z_comp_int(:,2) - mz;
        z3p = z_comp_int(:,3) - mz;

        l1p = sqrt(x1p.^2+y1p.^2+z1p.^2);
        l2p = sqrt(x2p.^2+y2p.^2+z2p.^2);
        l3p = sqrt(x3p.^2+y3p.^2+z3p.^2);

        s = (l1p+l2p+l12)/2;
        A_abp = sqrt(s.*(s-l1p).*(s-l2p).*(s-l12));
        s = (l1p+l3p+l13)/2;
        A_acp = sqrt(s.*(s-l1p).*(s-l3p).*(s-l13));
        s = (l2p+l3p+l23)/2;
        A_bcp = sqrt(s.*(s-l2p).*(s-l3p).*(s-l23));

        interpedx = (A_abp(:).*datx(:,3)+A_acp(:).*datx(:,2)+A_bcp(:).*datx(:,1))./A_abc;
        interpedy = (A_abp(:).*daty(:,3)+A_acp(:).*daty(:,2)+A_bcp(:).*daty(:,1))./A_abc;
        interpedz = (A_abp(:).*datz(:,3)+A_acp(:).*datz(:,2)+A_bcp(:).*datz(:,1))./A_abc;
end