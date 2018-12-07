function [x_comp_int,y_comp_int,z_comp_int,bx_comp_int,by_comp_int,bz_comp_int,tt] = interp_time(bx_comp,by_comp,bz_comp,x_comp,y_comp,z_comp,t_B,t_pos,sc)

    bx1 = bx_comp(:,1); bx1 = bx1(bx1 ~= 0 & ~isnan(bx1));
    bx2 = bx_comp(:,2); bx2 = bx2(bx2 ~= 0 & ~isnan(bx2));
    bx3 = bx_comp(:,3); bx3 = bx3(bx3 ~= 0 & ~isnan(bx3));
    by1 = by_comp(:,1); by1 = by1(by1 ~= 0 & ~isnan(by1));
    by2 = by_comp(:,2); by2 = by2(by2 ~= 0 & ~isnan(by2));
    by3 = by_comp(:,3); by3 = by3(by3 ~= 0 & ~isnan(by3));
    bz1 = bz_comp(:,1); bz1 = bz1(bz1 ~= 0 & ~isnan(bz1));
    bz2 = bz_comp(:,2); bz2 = bz2(bz2 ~= 0 & ~isnan(bz2));
    bz3 = bz_comp(:,3); bz3 = bz3(bz3 ~= 0 & ~isnan(bz3));

    px1 = x_comp(:,1); px1 = px1(px1 ~= 0 & ~isnan(px1));
    px2 = x_comp(:,2); px2 = px2(px2 ~= 0 & ~isnan(px2));
    px3 = x_comp(:,3); px3 = px3(px3 ~= 0 & ~isnan(px3));
    py1 = y_comp(:,1); py1 = py1(py1 ~= 0 & ~isnan(py1));
    py2 = y_comp(:,2); py2 = py2(py2 ~= 0 & ~isnan(py2));
    py3 = y_comp(:,3); py3 = py3(py3 ~= 0 & ~isnan(py3));
    pz1 = z_comp(:,1); pz1 = pz1(pz1 ~= 0 & ~isnan(pz1));
    pz2 = z_comp(:,2); pz2 = pz2(pz2 ~= 0 & ~isnan(pz2));
    pz3 = z_comp(:,3); pz3 = pz3(pz3 ~= 0 & ~isnan(pz3));

    t1 = t_B(:,1); t1 = t1(t1 ~= 0  & ~isnan(t1));
    t2 = t_B(:,2); t2 = t2(t2 ~= 0  & ~isnan(t2));
    t3 = t_B(:,3); t3 = t3(t3 ~= 0  & ~isnan(t3));

    p1 = t_pos(:,1); p1 = p1(p1 ~= 0  & ~isnan(p1));
    p2 = t_pos(:,2); p2 = p2(p2 ~= 0  & ~isnan(p2));
    p3 = t_pos(:,3); p3 = p3(p3 ~= 0  & ~isnan(p3));

    if sc == 4
        bx4 = bx_comp(:,4); bx4 = bx4(bx4 ~= 0 & ~isnan(bx4));
        by4 = by_comp(:,4); by4 = by4(by4 ~= 0 & ~isnan(by4));
        bz4 = bz_comp(:,4); bz4 = bz4(bz4 ~= 0 & ~isnan(bz4));
        px4 = x_comp(:,4); px4 = px4(px4 ~= 0 & ~isnan(px4));
        py4 = y_comp(:,4); py4 = py4(py4 ~= 0 & ~isnan(py4));
        pz4 = z_comp(:,4); pz4 = pz4(pz4 ~= 0 & ~isnan(pz4));
        t4 = t_B(:,4); t4 = t4(t4 ~= 0  & ~isnan(t4));
        p4 = t_pos(:,4); p4 = p4(p4 ~= 0  & ~isnan(p4));

        ind_1 = find([t1(end),t2(end),t3(end),t4(end)] - min([t1(end),t2(end),t3(end),t4(end)]) == 0);
        ind_2 = find([t1(1),t2(1),t3(1),t4(1)] - max([t1(1),t2(1),t3(1),t4(1)]) == 0);

        tt1 = []; tt2 = []; tt3 = []; tt4 = [];
        if ind_1 == 1
            tt2 = t2(t2<t1(end)); tt3 = t3(t3<t1(end)); tt4 = t4(t4<t1(end));
            c1 = t1; tt1 = t1;
        elseif ind_1 == 2
            tt1 = t1(t1<t2(end)); tt3 = t3(t3<t2(end)); tt4 = t4(t4<t2(end));
            c1 = t2; tt2 = t2;
        elseif ind_1 == 3
            tt2 = t2(t2<t3(end)); tt1 = t1(t1<t3(end)); tt4 = t4(t4<t3(end));
            c1 = t3; tt3 = t3;
        elseif ind_1 == 4
            tt2 = t2(t2<t4(end)); tt3 = t3(t3<t4(end)); tt1 = t1(t1<t4(end));
            c1 = t4; tt4 = t4;
        end
        if ind_2 == 1
            c2 = tt1;
        elseif ind_2 == 2
            c2 = tt2;
        elseif ind_2 == 3
            c2 = tt3;
        elseif ind_2 == 4
            c2 = tt4;
        end

        ind = find(c1 > c2(1));
        tt = c1(ind(1):end);

        pll = length(tt);
        x_comp_int = zeros(pll,4);
        y_comp_int = zeros(pll,4);
        z_comp_int = zeros(pll,4);
        bx_comp_int = zeros(pll,4);
        by_comp_int = zeros(pll,4);
        bz_comp_int = zeros(pll,4);

        for k = 1:4
            if k == 1
                t = t1; p = p1; bx = bx1; by = by1; bz = bz1; px = px1; py = py1; pz = pz1;
            elseif k == 2
                t = t2; p = p2; bx = bx2; by = by2; bz = bz2; px = px2; py = py2; pz = pz2;
            elseif k == 3
                t = t3; p = p3; bx = bx3; by = by3; bz = bz3; px = px3; py = py3; pz = pz3;
            else
                t = t4; p = p4; bx = bx4; by = by4; bz = bz4; px = px4; py = py4; pz = pz4;
            end
            x_comp_int(:,k) = interp1(p,px,tt);
            y_comp_int(:,k) = interp1(p,py,tt);
            z_comp_int(:,k) = interp1(p,pz,tt);
            bx_comp_int(:,k) = interp1(t,bx,tt);
            by_comp_int(:,k) = interp1(t,by,tt);  
            bz_comp_int(:,k) = interp1(t,bz,tt);
        end
    else
        ind_1 = find([t1(end),t2(end),t3(end)] - min([t1(end),t2(end),t3(end)]) == 0);
        ind_2 = find([t1(1),t2(1),t3(1)] - max([t1(1),t2(1),t3(1)]) == 0);

        tt1 = []; tt2 = []; tt3 = [];
        if ind_1 == 1
            tt2 = t2(t2<t1(end)); tt3 = t3(t3<t1(end));
            c1 = t1; tt1 = t1;
        elseif ind_1 == 2
            tt1 = t1(t1<t2(end)); tt3 = t3(t3<t2(end));
            c1 = t2; tt2 = t2;
        elseif ind_1 == 3
            tt2 = t2(t2<t3(end)); tt1 = t1(t1<t3(end));
            c1 = t3; tt3 = t3;
        end
        if ind_2 == 1
            c2 = tt1;
        elseif ind_2 == 2
            c2 = tt2;
        elseif ind_2 == 3
            c2 = tt3;
        end

        ind = find(c1 > c2(1));
        tt = c1(ind(1):end);

        pll = length(tt);
        x_comp_int = zeros(pll,3);
        y_comp_int = zeros(pll,3);
        z_comp_int = zeros(pll,3);
        bx_comp_int = zeros(pll,3);
        by_comp_int = zeros(pll,3);
        bz_comp_int = zeros(pll,3);

        for k = 1:3
            if k == 1
                t = t1; p = p1; bx = bx1; by = by1; bz = bz1; px = px1; py = py1; pz = pz1;
            elseif k == 2
                t = t2; p = p2; bx = bx2; by = by2; bz = bz2; px = px2; py = py2; pz = pz2;
            elseif k == 3
                t = t3; p = p3; bx = bx3; by = by3; bz = bz3; px = px3; py = py3; pz = pz3;
            end
            x_comp_int(:,k) = interp1(p,px,tt);
            y_comp_int(:,k) = interp1(p,py,tt);
            z_comp_int(:,k) = interp1(p,pz,tt);
            bx_comp_int(:,k) = interp1(t,bx,tt);
            by_comp_int(:,k) = interp1(t,by,tt);  
            bz_comp_int(:,k) = interp1(t,bz,tt);
        end
    end
end