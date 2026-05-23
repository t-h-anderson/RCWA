classdef Silver < rcwa.materials.Material
    % Silver permittivity from the original silver.mat tabulated data.

    methods
        function eps = permittivity(~, nmlambda)
            % The original `Silver('nmlambda', ...)` function returns the
            % complex refractive index (n + i k). BuildEps squared it to
            % get permittivity; do the same.
            % Use feval to avoid the in-package class name shadowing the
            % path function.
            n = feval('Silver', 'nmlambda', nmlambda);
            eps = n.^2;
        end
    end
end
