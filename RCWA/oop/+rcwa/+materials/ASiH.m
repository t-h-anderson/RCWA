classdef ASiH < rcwa.materials.Material
    % Hydrogenated amorphous silicon (a-Si:H) permittivity via the Cody
    % band-edge / Lorentz-oscillator model.
    %
    % The bandgap is itself spatially varying through the i-layer; the
    % `permittivity(nmlambda)` method requires an eVEg vector to be set.

    properties
        eVEg     (:,1) double = 1.6       % Bandgap profile per z-slice (eV)
        optimum  (1,1) double = 0         % 0 thesis material, 1 idealised
        epsinf   (1,1) double = 1
    end

    methods
        function obj = ASiH(varargin)
            for i = 1:2:numel(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
        end

        function eps = permittivity(obj, nmlambda)
            eps = aSiHGC(nmlambda, obj.eVEg.', obj.optimum, obj.epsinf);
        end
    end
end
