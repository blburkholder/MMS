function [b11,b22,v11,v22] = riemann_problemator(b1x,b2x,b1y,b2y,b1z,b2z,v1x,v1y,v1z,v2x,v2y,v2z)
    cz = 1;
    cy = cz*(b2x/b1x*b1z-b2z)/(b2y - b2x/b1x*b1y);
    cx = -(cy*b1y+cz*b1z)/b1x;
    
    dz = 1;
    dy = dz*((cx/(b1x-b2x))*(b1z-b2z)-cz)/(cy-(cx/(b1x-b2x))*(b1z-b2z));
    dx = -(dy*(b1y-b2y)+dz*(b1z-b2z))/(b1x-b2x);

    ez = 1;
    ey = ez*(dx/cx*cz-dz)/(dy-dx/cx*cy);
    ex = -(ey*cy+ez*cz)/cx;

    xp = [cx;cy;cz]/norm([cx,cy,cz]);
    yp = [dx;dy;dz]/norm([dx,dy,dz]);
    zp = [ex;ey;ez]/norm([ex,ey,ez]);

    b1 = [[b1x,b1y,b1z]*xp;[b1x,b1y,b1z]*yp;[b1x,b1y,b1z]*zp];
    b2 = [[b2x,b2y,b2z]*xp;[b2x,b2y,b2z]*yp;[b2x,b2y,b2z]*zp];

    v11 = [[v1x,v1y,v1z]*xp;[v1x,v1y,v1z]*yp;[v1x,v1y,v1z]*zp];
    v22 = [[v2x,v2y,v2z]*xp;[v2x,v2y,v2z]*yp;[v2x,v2y,v2z]*zp];

    b11 = b1/max(abs([b1(1) b1(2) b1(3) b2(1) b2(2) b2(3)]));
    b22 = b2/max(abs([b1(1) b1(2) b1(3) b2(1) b2(2) b2(3)]));