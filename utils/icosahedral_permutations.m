function perms = icosahedral_permutations(P)  
    Rset = icosahedral_rotations();
   
    perms = zeros(size(P,1),size(Rset,3));

    for i = 1:size(Rset,3)
        [~,I] = max(P * Rset(:,:,i)' * P');

        if (length(I) == length(unique(I)))
            perms(:,i) = I;
        else
            error("Rotations must be symmetries of the icosahedral mesh.");
        end
    end
end