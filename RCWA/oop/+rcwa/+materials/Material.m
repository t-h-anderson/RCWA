classdef (Abstract) Material
    % Abstract dispersive-material interface.
    %
    % Subclasses return the complex relative permittivity epsilon_r at a
    % set of wavelengths in nm.

    methods (Abstract)
        eps_r = permittivity(obj, nmlambda)
    end

    methods
        function n = refractiveIndex(obj, nmlambda)
            % Complex refractive index n = sqrt(eps_r). Convention: branch
            % with non-negative real part (MATLAB sqrt default for complex).
            n = sqrt(obj.permittivity(nmlambda));
        end
    end
end
