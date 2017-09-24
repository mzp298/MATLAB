function [ stress ] = VM_Stress(sx, sy, sz, txy, tyz, tzx)
%   This routine computes the von Mises stresses from the applied stresses
%
%   sx = stress in X direction
%   sy = stress in Y direction
%   sz = stress in Z direction
%   txy = shear on a plane normal to x in a y direction
%   txy = shear on a plane normal to y in a z direction
%   txy = shear on a plane normal to z in a z direction
%
%  The von Mises stress is returned to the calling function
%
   p1 = (sx - sy)^2 + (sy - sz)^2 + (sz - sx)^2;
   p2 = 6 * (txy^2 + tyz^2 + tzx^2);
   stress = sqrt((p1 + p2)/2);
end

