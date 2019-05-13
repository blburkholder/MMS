function [bx_comp_int,by_comp_int,bz_comp_int,tt] = interp_time_v(bx_comp,by_comp,bz_comp,t_B,sc)
    t1 = t_B(:,1); 
    bx1 = bx_comp(:,1); bx1 = bx1(t1 ~= 0 & ~isnan(bx1));
    by1 = by_comp(:,1); by1 = by1(t1 ~= 0 & ~isnan(by1));
    bz1 = bz_comp(:,1); bz1 = bz1(t1 ~= 0 & ~isnan(bz1));
    t1 = t1(t1 ~= 0  & ~isnan(t1));

    if sc == 4
        t2 = t_B(:,2); 
        t3 = t_B(:,3); 
        t4 = t_B(:,4); 
        bx2 = bx_comp(:,2); bx2 = bx2(t2 ~= 0 & ~isnan(bx2));
        by2 = by_comp(:,2); by2 = by2(t2 ~= 0 & ~isnan(by2));
        bz2 = bz_comp(:,2); bz2 = bz2(t2 ~= 0 & ~isnan(bz2));
        bx3 = bx_comp(:,3); bx3 = bx3(t3 ~= 0 & ~isnan(bx3));
        by3 = by_comp(:,3); by3 = by3(t3 ~= 0 & ~isnan(by3));
        bz3 = bz_comp(:,3); bz3 = bz3(t3 ~= 0 & ~isnan(bz3));
        bx4 = bx_comp(:,4); bx4 = bx4(t4 ~= 0 & ~isnan(bx4));
        by4 = by_comp(:,4); by4 = by4(t4 ~= 0 & ~isnan(by4));
        bz4 = bz_comp(:,4); bz4 = bz4(t4 ~= 0 & ~isnan(bz4));
        t2 = t2(t2 ~= 0  & ~isnan(t2));
        t3 = t3(t3 ~= 0  & ~isnan(t3));
        t4 = t4(t4 ~= 0  & ~isnan(t4));

        ind_1 = find([t1(end),t2(end),t3(end),t4(end)] - min([t1(end),t2(end),t3(end),t4(end)]) == 0,1);
        ind_2 = find([t1(1),t2(1),t3(1),t4(1)] - max([t1(1),t2(1),t3(1),t4(1)]) == 0,1);

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
        bx_comp_int = zeros(pll,4);
        by_comp_int = zeros(pll,4);
        bz_comp_int = zeros(pll,4);

        for k = 1:4
            if k == 1
                t = t1; bx = bx1; by = by1; bz = bz1; 
            elseif k == 2
                t = t2; bx = bx2; by = by2; bz = bz2;
            elseif k == 3
                t = t3; bx = bx3; by = by3; bz = bz3;
            else
                t = t4; bx = bx4; by = by4; bz = bz4;
            end
            bx_comp_int(:,k) = interp1(t,bx,tt);
            by_comp_int(:,k) = interp1(t,by,tt);  
            bz_comp_int(:,k) = interp1(t,bz,tt);
        end
    elseif sc == 3
        t2 = t_B(:,2); 
        t3 = t_B(:,3); 
        bx2 = bx_comp(:,2); bx2 = bx2(t2 ~= 0 & ~isnan(bx2));
        by2 = by_comp(:,2); by2 = by2(t2 ~= 0 & ~isnan(by2));
        bz2 = bz_comp(:,2); bz2 = bz2(t2 ~= 0 & ~isnan(bz2));
        bx3 = bx_comp(:,3); bx3 = bx3(t3 ~= 0 & ~isnan(bx3));
        by3 = by_comp(:,3); by3 = by3(t3 ~= 0 & ~isnan(by3));
        bz3 = bz_comp(:,3); bz3 = bz3(t3 ~= 0 & ~isnan(bz3));
        t2 = t2(t2 ~= 0  & ~isnan(t2));
        t3 = t3(t3 ~= 0  & ~isnan(t3));

        ind_1 = find([t1(end),t2(end),t3(end)] - min([t1(end),t2(end),t3(end)]) == 0,1);
        ind_2 = find([t1(1),t2(1),t3(1)] - max([t1(1),t2(1),t3(1)]) == 0,1);

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
        bx_comp_int = zeros(pll,3);
        by_comp_int = zeros(pll,3);
        bz_comp_int = zeros(pll,3);

        for k = 1:3
            if k == 1
                t = t1; bx = bx1; by = by1; bz = bz1;
            elseif k == 2
                t = t2; bx = bx2; by = by2; bz = bz2;
            elseif k == 3
                t = t3; bx = bx3; by = by3; bz = bz3;
            end
            bx_comp_int(:,k) = interp1(t,bx,tt);
            by_comp_int(:,k) = interp1(t,by,tt);  
            bz_comp_int(:,k) = interp1(t,bz,tt);
        end
    elseif sc == 2
        t2 = t_B(:,2);
        bx2 = bx_comp(:,2); bx2 = bx2(t2 ~= 0 & ~isnan(bx2));
        by2 = by_comp(:,2); by2 = by2(t2 ~= 0 & ~isnan(by2));
        bz2 = bz_comp(:,2); bz2 = bz2(t2 ~= 0 & ~isnan(bz2));
        t2 = t2(t2 ~= 0  & ~isnan(t2));

        ind_1 = find([t1(end),t2(end)] - min([t1(end),t2(end)]) == 0,1);
        ind_2 = find([t1(1),t2(1)] - max([t1(1),t2(1)]) == 0,1);

        tt1 = []; tt2 = [];
        if ind_1 == 1
            tt2 = t2(t2<t1(end));
            c1 = t1; tt1 = t1;
        elseif ind_1 == 2
            tt1 = t1(t1<t2(end));
            c1 = t2; tt2 = t2;
        end
        if ind_2 == 1
            c2 = tt1;
        elseif ind_2 == 2
            c2 = tt2;
        end

        ind = find(c1 > c2(1));
        tt = c1(ind(1):end);

        pll = length(tt);
        bx_comp_int = zeros(pll,2);
        by_comp_int = zeros(pll,2);
        bz_comp_int = zeros(pll,2);

        for k = 1:2
            if k == 1
                t = t1; bx = bx1; by = by1; bz = bz1;
            elseif k == 2
                t = t2; bx = bx2; by = by2; bz = bz2;
            end
            bx_comp_int(:,k) = interp1(t,bx,tt);
            by_comp_int(:,k) = interp1(t,by,tt);  
            bz_comp_int(:,k) = interp1(t,bz,tt);
        end
    elseif sc == 1
        bx_comp_int = bx1;
        by_comp_int = by1;
        bz_comp_int = bz1;
        tt = t1;
    end
end