% function [theta,phi] = cart_to_sphere(x,y,z)
%     theta = acos(z);
%     phi = atan2(y,x);
%     
%     phi(phi < 0) = phi(phi < 0) + 2*pi;
% end

function [theta,phi] = cart_to_sphere(x)
    theta = acos(x(:,3));
    phi = atan2(x(:,2),x(:,1));
    
    phi(phi < 0) = phi(phi < 0) + 2*pi;    
end