function [interpedx,interpedy,interpedz] = nsimplex_4vec(datx,daty,datz,x_comp_int,y_comp_int,z_comp_int)
        mx = mean(x_comp_int,2);
        my = mean(y_comp_int,2);
        mz = mean(z_comp_int,2);

        xa = x_comp_int(:,1) - x_comp_int(:,4);
        xb = x_comp_int(:,2) - x_comp_int(:,4);
        xc = x_comp_int(:,3) - x_comp_int(:,4);
        ya = y_comp_int(:,1) - y_comp_int(:,4);
        yb = y_comp_int(:,2) - y_comp_int(:,4);
        yc = y_comp_int(:,3) - y_comp_int(:,4);
        za = z_comp_int(:,1) - z_comp_int(:,4);
        zb = z_comp_int(:,2) - z_comp_int(:,4);
        zc = z_comp_int(:,3) - z_comp_int(:,4);
        v_abcd = abs(xa.*(yb.*zc-zb.*yc) - ya.*(xb.*zc-zb.*xc)+za.*(xb.*yc-yb.*xc))/6;

        xa = x_comp_int(:,1) - mx;
        xb = x_comp_int(:,2) - mx;
        xc = x_comp_int(:,3) - mx;
        xd = x_comp_int(:,4) - mx;
        ya = y_comp_int(:,1) - my;
        yb = y_comp_int(:,2) - my;
        yc = y_comp_int(:,3) - my; 
        yd = y_comp_int(:,4) - my;
        za = z_comp_int(:,1) - mz;
        zb = z_comp_int(:,2) - mz;
        zc = z_comp_int(:,3) - mz;
        zd = z_comp_int(:,4) - mz;

        v_abcp = abs(xa.*(yb.*zc-zb.*yc) - ya.*(xb.*zc-zb.*xc)+za.*(xb.*yc-yb.*xc))/6;
        v_abdp = abs(xa.*(yb.*zd-zb.*yd) - ya.*(xb.*zd-zb.*xd)+za.*(xb.*yd-yb.*xd))/6;
        v_acdp = abs(xa.*(yc.*zd-zc.*yd) - ya.*(xc.*zd-zc.*xd)+za.*(xc.*yd-yc.*xd))/6;
        v_bcdp = abs(xb.*(yc.*zd-zc.*yd) - yb.*(xc.*zd-zc.*xd)+zb.*(xc.*yd-yc.*xd))/6;

        interpedx = (datx(:,4).*v_abcp + datx(:,3).*v_abdp + datx(:,2).*v_acdp + datx(:,1).*v_bcdp)./v_abcd;
        interpedy = (daty(:,4).*v_abcp + daty(:,3).*v_abdp + daty(:,2).*v_acdp + daty(:,1).*v_bcdp)./v_abcd;
        interpedz = (datz(:,4).*v_abcp + datz(:,3).*v_abdp + datz(:,2).*v_acdp + datz(:,1).*v_bcdp)./v_abcd;
end