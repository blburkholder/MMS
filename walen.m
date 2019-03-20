function [vjump,vajump,rsqr,ind1,ind2] = walen(jt,vx,vy,vz,vax,vay,vaz,bx,by,bz,eppps)  
    wlh = 1+eppps;%walen limit high
    wll = 1-eppps;%walen limit low
    ind1 = 0;
    ind2 = 0;
    wind_inc = 1;

    max_windl = 2;
    vxd = zeros(6,1);
    vyd = zeros(6,1);
    vzd = zeros(6,1);
    vaxd = zeros(6,1);
    vayd = zeros(6,1);
    vazd = zeros(6,1);

%     figure
%     plot(jt,vx-mean(vx),'r'); hold on; plot(jt,vy-mean(vy),'g'); plot(jt,vz-mean(vz),'b');
%     plot(jt,vax-mean(vax),'r--'); plot(jt,vay-mean(vay),'g--'); plot(jt,vaz-mean(vaz),'b--');
    start = max_windl+1; % index within data stream
    first = 0; % number of consecutive windows satisfying walen relation
    pmax = 0; % number of consecutive windows for the interval which is chosen
    here = 0; % start index of window(s) satisfying walen relation
    heremax = 0; % start index of window(s) for the interval which is chosen
    rsqr = 0; % R^2 value for interval satisfying walen relation
    bjump = 0; % magnetic field jump across interval satisfying walen relation
    bjumpmax = 0; % magnetic field jump across interval which is chosen
    while start <= length(vx) - max_windl
        for wind = wind_inc:wind_inc:max_windl
            vxd(wind/wind_inc) = abs(vx(start-wind) - vx(start+wind)); 
            vyd(wind/wind_inc) = abs(vy(start-wind) - vy(start+wind));     
            vzd(wind/wind_inc) = abs(vz(start-wind) - vz(start+wind));     
            vaxd(wind/wind_inc) = abs(vax(start-wind) - vax(start+wind));     
            vayd(wind/wind_inc) = abs(vay(start-wind) - vay(start+wind));     
            vazd(wind/wind_inc) = abs(vaz(start-wind) - vaz(start+wind));
        end

        vxd(3) = abs(vx(start) - vx(start+1)); vxd(4) = abs(vx(start) - vx(start-1));
        vxd(5) = abs(vx(start) - vx(start+2)); vxd(6) = abs(vx(start) - vx(start-2));
        vyd(3) = abs(vy(start) - vy(start+1)); vyd(4) = abs(vy(start) - vy(start-1));
        vyd(5) = abs(vy(start) - vy(start+2)); vyd(6) = abs(vy(start) - vy(start-2));     
        vzd(3) = abs(vz(start) - vz(start+1)); vzd(4) = abs(vz(start) - vz(start-1));
        vzd(5) = abs(vz(start) - vz(start+2)); vzd(6) = abs(vz(start) - vz(start-2));     
        vaxd(3) = abs(vax(start) - vax(start+1)); vaxd(4) = abs(vax(start) - vax(start-1));
        vaxd(5) = abs(vax(start) - vax(start+2)); vaxd(6) = abs(vax(start) - vax(start-2));
        vayd(3) = abs(vay(start) - vay(start+1)); vayd(4) = abs(vay(start) - vay(start-1));
        vayd(5) = abs(vay(start) - vay(start+2)); vayd(6) = abs(vay(start) - vay(start-2));    
        vazd(3) = abs(vaz(start) - vaz(start+1)); vazd(4) = abs(vaz(start) - vaz(start-1));
        vazd(5) = abs(vaz(start) - vaz(start+2)); vazd(6) = abs(vaz(start) - vaz(start-2));

        pxx = [vaxd;vayd;vazd]\[vxd;vyd;vzd];
        %pyy = vayd\vyd;
        %pzz = vazd\vzd;
        mslope = pxx*[vxd;vyd;vzd];
        %fy = pyy*vyd;
        %fz = pzz*vzd;       
        [goodness ,~] = rsquare([vaxd;vayd;vazd],mslope);
        %[gy ,~] = rsquare(vayd,fy);
        %[gz ,~] = rsquare(vazd,fz);
        %goodness = (gx+gy+gz)/3;
        %mslope = (pxx(1)+pyy(1)+pzz(1))/3;

        %if (pxx(1) >= wll) && (pxx(1) <= wlh) && (pzz(1) >= wll) && (pzz(1) <= wlh) && (mslope >= wll) && (mslope <= wlh)
        if (pxx(1) >= wll) && (pxx(1) <= wlh) && goodness > 0.5
            first = first + 1;
            if first == 1
                here = start;
                rsqr = 0;
                %figure; hold on
            end
            bjump = (bx(here-max_windl)-bx(start+max_windl))^2+(by(here-max_windl)-...
                by(start+max_windl))^2+(bz(here-max_windl)-bz(start+max_windl))^2;
            rsqr = rsqr + goodness;
            %scatter(vxd,vaxd,'r.')
            %scatter(vyd,vayd,'g.')
            %scatter(vzd,vazd,'b.') 
            %title(['#wind - ',num2str(first),' R^2 = ',num2str(rsqr/first),' mslope = ',num2str(pxx(1))])
        else
            if bjump > bjumpmax
                bjumpmax = bjump;
                heremax = here;
                pmax = first;
                rsqr = rsqr/pmax;
            end
            first = 0;
            bjump = 0;
        end
        start = start+1;
    end
    if heremax ~= 0
        ind1 = heremax-max_windl;
        ind2 = heremax+pmax+max_windl;
        vjump(1) = abs(vx(ind1)-vx(ind2));
        vjump(2) = abs(vy(ind1)-vy(ind2));
        vjump(3) = abs(vz(ind1)-vz(ind2));
        vajump(1) = abs(vax(ind1)-vax(ind2));
        vajump(2) = abs(vay(ind1)-vay(ind2));
        vajump(3) = abs(vaz(ind1)-vaz(ind2));
    else
        vjump(1) = 0; vjump(2) = 0; vjump(3) = 0;
        vajump(1) = 0; vajump(2) = 0; vajump(3) = 0;
    end
end