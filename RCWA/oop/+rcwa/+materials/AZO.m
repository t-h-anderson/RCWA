classdef AZO < rcwa.materials.Material
    % Aluminium-doped zinc oxide (AZO) permittivity.
    %
    % The original AZO.m returns the complex refractive index (n+ik); the
    % original BuildEps stored that *directly* into loc.epsd / loc.epsw
    % without squaring. We preserve that behaviour here (see BUGS.md
    % BUG-007 -- it is almost certainly wrong but we leave it faithful).
    %
    % Pass squared = true to instead get the physically correct
    % eps = (n+ik)^2.

    properties
        squared (1,1) logical = false
    end

    methods
        function obj = AZO(varargin)
            for i = 1:2:numel(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
        end

        function eps = permittivity(obj, nmlambda)
            % feval avoids the in-package class name shadowing the
            % path function of the same name.
            n = feval('AZO', nmlambda);
            if obj.squared
                eps = n.^2;
            else
                eps = n;        % faithful to BuildEps (BUG-007)
            end
        end
    end
end
