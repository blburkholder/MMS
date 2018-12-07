function [dvdva,jump,shear,cs_vx,rsqr,ind1,ind2] = walen(jt,vx,vy,vz,vax,vay,vaz,bx,by,bz,wll,wlh)  
    %wlh = walen limit high
    %wll = walen limit low
    ind1 = 0;
    ind2 = 0;

    px = zeros(length(vx),1);
    max_windl = 32;
    wind_inc = 4;
    vxd = zeros(max_windl/wind_inc,1);
    vyd = zeros(max_windl/wind_inc,1);
    vzd = zeros(max_windl/wind_inc,1);
    vaxd = zeros(max_windl/wind_inc,1);
    vayd = zeros(max_windl/wind_inc,1);
    vazd = zeros(max_windl/wind_inc,1);

%     figure
%     plot(jt,vx-mean(vx),'r'); hold on; plot(jt,vy-mean(vy),'g'); plot(jt,vz-mean(vz),'b');
%     plot(jt,vax-mean(vax),'r--'); plot(jt,vay-mean(vay),'g--'); plot(jt,vaz-mean(vaz),'b--');
    start = max_windl+1;
    first = 0;
    pmax = 0;
    here = 0;
    heremax = 0;
    rsqr = 0;
    bjump = 0;
    bjumpmax = 0;
    while start <= length(vx) - max_windl
        for wind = wind_inc:wind_inc:max_windl
            vxd(wind/wind_inc) = abs(vx(start-wind) - vx(start+wind)); 
            vyd(wind/wind_inc) = abs(vy(start-wind) - vy(start+wind));     
            vzd(wind/wind_inc) = abs(vz(start-wind) - vz(start+wind));     
            vaxd(wind/wind_inc) = abs(vax(start-wind) - vax(start+wind));     
            vayd(wind/wind_inc) = abs(vay(start-wind) - vay(start+wind));     
            vazd(wind/wind_inc) = abs(vaz(start-wind) - vaz(start+wind));
        end

        p = polyfit([vxd;vyd;vzd],[vaxd;vayd;vazd],1);
        f = polyval(p,[vxd;vyd;vzd]);
        [goodness ,~] = rsquare([vaxd;vayd;vazd],f);

        px(start) = p(1);
        if (p(1) >= wll) && (p(1) <= wlh)
            first = first + 1;
            if first == 1
                here = start;
                rsqr = 0;
            end
            start = start+max_windl;
            %if first > pmax
            %    pmax = first;
            %    heremax = here;
            rsqr = rsqr + goodness;
            bjump = (bx(here)-bx(start))^2+(by(here)-by(start))^2+(bz(here)-bz(start))^2;
            %end
        else
            if bjump > bjumpmax
                bjumpmax = bjump;
                heremax = here;
                pmax = first;
                rsqr = rsqr/pmax;
            end
            first = 0;
            start = start+1;
            bjump = 0;
        end
    end
    if heremax ~= 0
        ind1 = heremax-max_windl/2;
        ind2 = heremax+(pmax+1)*max_windl/2;
        jump = sqrt(max([(bx(ind1)-bx(ind2))^2,(by(ind1)-by(ind2))^2,(bz(ind1)-bz(ind2))^2]));
        shear = (bx(ind1)*bx(ind2) + by(ind1)*by(ind2) + bz(ind1)*bz(ind2)) / ...
                    sqrt((bx(ind1)^2 + by(ind1)^2 + bz(ind1)^2)*(bx(ind2)^2 + by(ind2)^2 + bz(ind2)^2));
    else
        jump = 0;
        shear =0;
    end
    dvdva = pmax;
    cs_vx = nanmean(vx);
    %rsqr = rsqr/pmax;
end