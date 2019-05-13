function [intervals] = find_cs(jt,bx,by,bz)

    start_wind = 100;   
    wind = start_wind;
    xf = zeros(length(bx),1);
    yf = zeros(length(by),1);
    zf = zeros(length(bz),1);

    resfrac = ceil(length(bx)/2000);
    inds = 1:length(bx);
    inds = inds(1:resfrac:end);

    %figure('visible','off')
    figure('rend','painters','pos',[10 10 1000 800])
    subplot(3,1,1); hold on;
    plot(inds,bx(1:resfrac:end),'r');plot(inds,by(1:resfrac:end),'g');plot(inds,bz(1:resfrac:end),'b');
    %xlim([1 length(bx)])
    subplot(3,1,2); hold on;
    %xlim([1 length(bx)])

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
    plot(inds,xf(1:resfrac:end),'r')
    plot(inds,yf(1:resfrac:end),'g')
    plot(inds,zf(1:resfrac:end),'b')

    jumpfuck = (xf+yf+zf)/3;
    subplot(3,1,3); hold on
    plot(inds,jumpfuck(1:resfrac:end),'k')
    [peaksx,loc] = findpeaks(jumpfuck);
    plot(loc,peaksx,'kv')

    [~,so] = sort(peaksx,'descend');

    intervals = [];
    accepted = 0;
    more_cs = 1;
    ints = 0;
    if input('do it?')
        while more_cs
            while ~accepted
                x = input('which fucking peak?');
                plot([loc(so(x)),loc(so(x))],[0,1],'k')
                accepted = input('do you fucking accept?');
            end
            accepted = 0;
            ints = ints + 1;
            while ~accepted
                window = ginput(1);
                plot([2*loc(so(x))-window(1),2*loc(so(x))-window(1)],[0,1],'r')
                plot([window(1),window(1)],[0,1],'r')
                accepted = input('do you fucking accept?');
            end
            intervals(ints,1) = loc(so(x));
            intervals(ints,2) = abs(window(1) - loc(so(x)));
            more_cs = input('any more fucking current sheets?');
            accepted = 0;
        end
    end
end
            

