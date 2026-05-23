classdef Constant < rcwa.materials.Material
    % Wavelength-independent permittivity (a placeholder for "epsm = 1+3i").

    properties
        eps_r (1,1) double
    end

    methods
        function obj = Constant(eps_r)
            obj.eps_r = eps_r;
        end

        function eps = permittivity(obj, nmlambda)
            eps = obj.eps_r * ones(size(nmlambda));
        end
    end
end
