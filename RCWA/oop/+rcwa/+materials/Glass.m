classdef Glass < rcwa.materials.Material
    % Dispersive glass refractive index table (450-1300 nm).
    %
    % Wraps the tabulated function `glass.m` from the original code.
    % By default reproduces the BuildEps convention:
    %     eps = (n + 1e-3 i)^2
    % i.e. a small imaginary part is added to the index before squaring
    % to give the glass a tiny absorption. Pass `absorption = 0` to
    % recover the pure real n^2.

    properties
        absorption (1,1) double = 1e-3
    end

    methods
        function obj = Glass(varargin)
            for i = 1:2:numel(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
        end

        function eps = permittivity(obj, nmlambda)
            n = glass('nmlambda', nmlambda);
            eps = (n + obj.absorption * 1i).^2;
        end
    end
end
