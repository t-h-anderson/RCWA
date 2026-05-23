classdef Sinusoidal < rcwa.gratings.Grating
    % Sinusoidal grating relief.

    methods
        function obj = Sinusoidal(varargin)
            obj@rcwa.gratings.Grating(varargin{:});
        end

        function g = relief(obj, nmx)
            idx = sign(abs(nmx) - 0.5 * obj.zeta * obj.nmLx);
            g = zeros(size(nmx));
            mask = (idx == -1);
            g(mask) = obj.nmLg * sin( ...
                pi * (nmx(mask) + 0.5 * obj.zeta * obj.nmLx) ...
                ./ (obj.zeta * obj.nmLx));
        end

        function eps_slice = epsAt(obj, nmx, nmz, nmLm, epsm, epsd, ~)
            eps_slice = obj.epsAtCommon(nmx, nmz, nmLm, epsm, epsd);
        end
    end
end
