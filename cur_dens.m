%problems:
%25,54,74,103
function [j1,j2,bx_comp_int,by_comp_int,bz_comp_int,x_comp_int,y_comp_int,z_comp_int,tt,B_tag,sc] = cur_dens(pathpath,B_path1,B_path2)
    [x_comp,y_comp,z_comp,t_pos,bx_comp,by_comp,bz_comp,t_B,B_tag,sc] = MMS_fgm(pathpath,B_path1,B_path2);
    [x_comp_int,y_comp_int,z_comp_int,bx_comp_int,by_comp_int,bz_comp_int,tt] =...
        interp_time_B(bx_comp,by_comp,bz_comp,x_comp,y_comp,z_comp,t_B,t_pos,sc);

    j2 = [];
      
     if sc == 4  
%         [j1x,j1y,j1z] = cur_dens_curl(x_comp_int(:,1),y_comp_int(:,1),z_comp_int(:,1),bx_comp_int(:,1),by_comp_int(:,1),bz_comp_int(:,1));
%         [j2x,j2y,j2z] = cur_dens_curl(x_comp_int(:,2),y_comp_int(:,2),z_comp_int(:,2),bx_comp_int(:,2),by_comp_int(:,2),bz_comp_int(:,2));
%         [j3x,j3y,j3z] = cur_dens_curl(x_comp_int(:,3),y_comp_int(:,3),z_comp_int(:,3),bx_comp_int(:,3),by_comp_int(:,3),bz_comp_int(:,3));
%         [j4x,j4y,j4z] = cur_dens_curl(x_comp_int(:,4),y_comp_int(:,4),z_comp_int(:,4),bx_comp_int(:,4),by_comp_int(:,4),bz_comp_int(:,4));
% 
%         [jx2,jy2,jz2] = nsimplex_4([j1x,j2x,j3x,j4x],[j1y,j2y,j3y,j4y],[j1z,j2z,j3z,j4z],x_comp_int,y_comp_int,z_comp_int);
%         j2 = [jx2(:),jy2(:),jz2(:)];
% 
         d_R2 = abs([x_comp_int(:,4),y_comp_int(:,4),z_comp_int(:,4)] - [x_comp_int(:,2),y_comp_int(:,2),z_comp_int(:,2)]);
         d_R3 = abs([x_comp_int(:,4),y_comp_int(:,4),z_comp_int(:,4)] - [x_comp_int(:,3),y_comp_int(:,3),z_comp_int(:,3)]);
         d_R4 = abs([x_comp_int(:,4),y_comp_int(:,4),z_comp_int(:,4)] - [x_comp_int(:,1),y_comp_int(:,1),z_comp_int(:,1)]);
         d_B2 = abs([bx_comp_int(:,4),by_comp_int(:,4),bz_comp_int(:,4)] - [bx_comp_int(:,2),by_comp_int(:,2),bz_comp_int(:,2)]);
         d_B3 = abs([bx_comp_int(:,4),by_comp_int(:,4),bz_comp_int(:,4)] - [bx_comp_int(:,3),by_comp_int(:,1),bz_comp_int(:,3)]);
         d_B4 = abs([bx_comp_int(:,4),by_comp_int(:,4),bz_comp_int(:,4)] - [bx_comp_int(:,1),by_comp_int(:,3),bz_comp_int(:,1)]);
 
         cdr23 = cross(d_R2,d_R3);
         cdr34 = cross(d_R3,d_R4);
         cdr42 = cross(d_R4,d_R2);
 
         s23 = sqrt(sum(cdr23.^2,2));
         s34 = sqrt(sum(cdr34.^2,2));
         s42 = sqrt(sum(cdr42.^2,2));
 
         j23 = (sum(d_B2.*d_R3,2) - sum(d_B3.*d_R2,2))./s23;
         j34 = (sum(d_B3.*d_R4,2) - sum(d_B4.*d_R3,2))./s34;
         j42 = (sum(d_B4.*d_R2,2) - sum(d_B2.*d_R4,2))./s42;
 
         e23 = cdr23./[s23,s23,s23];
         e34 = cdr34./[s34,s34,s34];
         e42 = cdr42./[s42,s42,s42];
 
         e1 = zeros(length(bx_comp_int),3);
         e2 = zeros(length(bx_comp_int),3);
         e3 = zeros(length(bx_comp_int),3);
         e1(:,1) = 1; e2(:,2) = 1; e3(:,3) = 1;
 
         jx1 = j23.*sum(e23.*e1,2) + j34.*sum(e34.*e1,2) + j42.*sum(e42.*e1,2);
         jy1 = j23.*sum(e23.*e2,2) + j34.*sum(e34.*e2,2) + j42.*sum(e42.*e2,2);
         jz1 = j23.*sum(e23.*e3,2) + j34.*sum(e34.*e3,2) + j42.*sum(e42.*e3,2);
         j1 = [jx1(:),jy1(:),jz1(:)];
% 
% %          j_spaceframe_square = zeros(length(j23),1);
% %         for s = 1:length(j23)
% %             j_spaceframe_square(s) = gendotprod([j23(s),j34(s),j42(s)],[j23(s),j34(s),j42(s)],e23(s,:),e34(s,:),e42(s,:));
% %         end
     else
%         [j1x,j1y,j1z] = cur_dens_curl(x_comp_int(:,1),y_comp_int(:,1),z_comp_int(:,1),bx_comp_int(:,1),by_comp_int(:,1),bz_comp_int(:,1));
%         [j2x,j2y,j2z] = cur_dens_curl(x_comp_int(:,2),y_comp_int(:,2),z_comp_int(:,2),bx_comp_int(:,2),by_comp_int(:,2),bz_comp_int(:,2));
%         [j3x,j3y,j3z] = cur_dens_curl(x_comp_int(:,3),y_comp_int(:,3),z_comp_int(:,3),bx_comp_int(:,3),by_comp_int(:,3),bz_comp_int(:,3));
% 
%         [jx2,jy2,jz2] = nsimplex_3([j1x,j2x,j3x],[j1y,j2y,j3y],[j1z,j2z,j3z],x_comp_int,y_comp_int,z_comp_int);
%         j2 = [jx2(:),jy2(:),jz2(:)];
         j1 = [];
     end
end

