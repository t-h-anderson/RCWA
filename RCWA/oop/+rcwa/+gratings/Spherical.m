classdef Spherical < rcwa.gratings.Grating
    % Spherical (circular cross-section) grating sunk into the mirror.

    properties
        nmda (1,1) double = 0   % grating-region thickness needed for geometry
    end

    methods
        function obj = Spherical(varargin)
            obj@rcwa.gratings.Grating(varargin{:});
        end

        function g = relief(obj, nmx)
            idx = sign(abs(nmx) - 0.5 * obj.zeta * obj.nmLx);
            g = zeros(size(nmx));
            mask = (idx == -1);
            c = 0.5 * obj.zeta * obj.nmLx;
            a = 0.5 * (c^2/obj.nmLg - obj.nmda);
            r = a + obj.nmLg;
            g(mask) = sqrt(r^2 - nmx(mask).^2) - a;
        end

        function eps_slice = epsAt(obj, nmx, nmz, nmLm, epsm, epsd, ~)
            eps_slice = obj.epsAtCommon(nmx, nmz, nmLm, epsm, epsd);
        end
    end
end
