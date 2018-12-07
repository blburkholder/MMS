function [wind_lengs,cs_centers] = find_cs(jt,bx,by,bz,jjj)

    ink_inds = [15,16,26,28,31,43,45,80,94,96,98,99,100,105,...
            107,112,114,115,117,119,128,133,146,148,150,156,163,166,...
            171,177,179,183,185,186,187,192,210,212,216,217,224,...
            234,235,241,256,260,268,269,271,277,278,279,283,284];
    start_wind = 100;   
    wind = start_wind;
    xf = zeros(length(bx),1);
    yf = zeros(length(by),1);
    zf = zeros(length(bz),1);

%     figure('visible','off')
%     %figure
%     subplot(3,1,1); hold on;
%     plot(jt,bx,'r');plot(jt,by,'g');plot(jt,bz,'b');
%     xlim([1 length(bx)])
%     subplot(3,1,2); hold on;
%     xlim([1 length(bx)])

    while wind < length(bx)/4
        for i = 1:length(bx)-wind
            bxs = bx(i:i+wind);
            bys = by(i:i+wind);
            bzs = bz(i:i+wind);

            B = sqrt(mean(bxs)^2+mean(bys)^2+mean(bzs)^2);

            bx1 = bxs(1:ceil(length(bxs)/2));
            by1 = bys(1:ceil(length(bxs)/2));
            bz1 = bzs(1:ceil(length(bxs)/2));
            bx2 = bxs(floor(length(bxs)/2):end);
            by2 = bys(floor(length(bxs)/2):end);
            bz2 = bzs(floor(length(bxs)/2):end);
            xf(i+wind/2) = xf(i+wind/2) + abs(mean(bx1)-mean(bx2))/B;
            yf(i+wind/2) = yf(i+wind/2) + abs(mean(by1)-mean(by2))/B;
            zf(i+wind/2) = zf(i+wind/2) + abs(mean(bz1)-mean(bz2))/B;
        end
        wind = wind+100;
    end
%     plot(jt,xf,'r')
%     plot(jt,yf,'g')
%     plot(jt,zf,'b')
%     plot(jt,(xf+yf+zf)/3,'k')

    jumpfuck = (xf+yf+zf)/3;
    [peaksx,loc] = findpeaks(jumpfuck);
%     plot(jt(loc),peaksx,'kv')
    
    [peaksx,px] = sort(peaksx,'descend');
    loc = loc(px);

    csss = hand_determined_cs(jjj);
    wind_lengs = zeros(sum(csss~=0),1);
    cs_centers = zeros(sum(csss~=0),1);

    this_guy = 1;
    for z = csss(csss ~= 0)
%         subplot(3,1,2)
%         plot(jt(loc(z)),peaksx(z),'kv','MarkerFaceColor','k')
        wx_help = 0;
        wy_help = 0;
        wz_help = 0;
        for wll = 100:50:500
            difx = abs(bx(loc(z)-wll/2) - bx(loc(z)+wll/2));
            if difx > wx_help
                wx_help = difx;
                wlx = wll;
            else
                break;
            end
        end

        for wll = 100:50:500
            dify = abs(by(loc(z)-wll/2) - by(loc(z)+wll/2));
            if dify > wy_help
                wy_help = dify;
                wly = wll;
            else
                break;
            end
        end

        for wll = 100:50:500
            difz = abs(bz(loc(z)-wll/2) - bz(loc(z)+wll/2));
            if difz > wz_help
                wz_help = difz;
                wlz = wll;
            else
                break;
            end
        end
        wl = max([wlx,wly,wlz]);
        if (sum(ink_inds == jjj) > 0)
            wl = wl*2;
        end    
%         subplot(3,1,1)
%         plot(jt((loc(z)-wl):(loc(z)+wl)),bx((loc(z)-wl):(loc(z)+wl)),'k.')
%         plot(jt((loc(z)-wl):(loc(z)+wl)),by((loc(z)-wl):(loc(z)+wl)),'k.')
%         plot(jt((loc(z)-wl):(loc(z)+wl)),bz((loc(z)-wl):(loc(z)+wl)),'k.')
        cs_centers(this_guy) = loc(z);
        wind_lengs(this_guy) = wl;
        this_guy = this_guy + 1;
    end
end
            

