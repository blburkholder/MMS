function [hi] = hamperl(dfft,wind)
    h = zeros(size(dfft));
    halfsies = ceil(length(dfft)/2);
    std_mult = 1;
    %figure; hold on
    %loglog(abs(dfft),'k')
    for i = 1:halfsies
        h(i) = dfft(i); 
        if i > floor(wind/2)
            if i > (halfsies-floor(wind/2))
                rr = dfft(halfsies-wind:halfsies);
            else
                rr = dfft(i-floor(wind/2):i+floor(wind/2));                
            end
            mr = median(rr);
            stdr = std(rr);

            if abs(dfft(i)) - abs(mr) > std_mult*stdr
                h(i) = mr;
            end
        end
    end
    %h(end-50:end) = 10;
    if halfsies == length(dfft)/2;
        h(end-halfsies+2:end) = flipud(h(2:halfsies));
        h(halfsies+1) = h(halfsies);
    else
        h(end-halfsies+2:end) = flipud(h(2:halfsies));
        %h(halfsies+1) = h(halfsies);
    end
    %loglog(abs(h),'r')
    %set(gca,'XScale','log')
    %set(gca,'YScale','log')
    hi = ifft(h,'symmetric');
    %loglog(abs(h(1:131)),'k')
end