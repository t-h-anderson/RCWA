classdef FTO < rcwa.materials.Material
    % Fluorine-doped tin oxide (FTO) refractive index from FTOnk1200 table.

    methods
        function eps = permittivity(~, nmlambda)
            n = feval('FTOnk', nmlambda);
            eps = n.^2;
        end
    end
end
