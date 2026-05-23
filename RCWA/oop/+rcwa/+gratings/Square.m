classdef Square < rcwa.gratings.Grating
    % Rectangular (square) grating.
    %
    % The relief is `nmLg` between |x| <= zeta*Lx/2 and 0 elsewhere.

    methods
        function obj = Square(varargin)
            obj@rcwa.gratings.Grating(varargin{:});
        end

        function g = relief(obj, nmx)
            % Faithful port of the original gx.m for type==0.
            idx = sign(abs(nmx) - 0.5 * obj.zeta * obj.nmLx);
            g = zeros(size(nmx));
            g(idx ==  1) = 0;
            g(idx ==  0) = obj.nmLg;
            g(idx == -1) = obj.nmLg;
        end

        function eps_slice = epsAt(obj, nmx, nmz, nmLm, epsm, epsd, ~)
            eps_slice = obj.epsAtCommon(nmx, nmz, nmLm, epsm, epsd);
        end

        function f = fourierExplicit(obj, n, epsm, epsd)
            % Closed-form Fourier coefficients of the rectangular profile
            % (used in the grating region only). `n` is a vector of mode
            % indices; n == 0 must be handled by the caller (average).
            f = (epsm - epsd) .* sin(n * pi * obj.zeta) ./ (pi * n);
        end
    end
end
