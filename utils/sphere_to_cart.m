function pt = sphere_to_cart(theta,phi)
    x = sin(theta).*cos(phi);
    y = sin(theta).*sin(phi);
    z = cos(theta);
    
    pt = [x,y,z];
    pt = normr(pt).';
end
