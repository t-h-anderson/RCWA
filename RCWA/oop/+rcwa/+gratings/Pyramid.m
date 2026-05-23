classdef Pyramid < rcwa.gratings.Grating
    % Triangular (pyramid cross-section) grating.

    methods
        function obj = Pyramid(varargin)
            obj@rcwa.gratings.Grating(varargin{:});
        end

        function g = relief(obj, nmx)
            idx = sign(abs(nmx) - 0.5 * obj.zeta * obj.nmLx);
            g = zeros(size(nmx));
            mask = (idx == -1);
            m = obj.nmLg / (obj.nmLx * obj.zeta * 0.5);
            g(mask) = obj.nmLg - m * abs(nmx(mask));
        end

        function eps_slice = epsAt(obj, nmx, nmz, nmLm, epsm, epsd, ~)
            eps_slice = obj.epsAtCommon(nmx, nmz, nmLm, epsm, epsd);
        end
    end
end
