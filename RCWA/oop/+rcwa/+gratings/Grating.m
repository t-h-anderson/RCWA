classdef (Abstract) Grating
    % Abstract base for grating relief profiles z = g(x).
    %
    % Subclasses define the profile shape used in the grating region. All
    % subclasses take a period Lx (nm), duty cycle zeta in (0,1] and a
    % peak-to-trough relief height Lg (nm).

    properties
        nmLx (1,1) double {mustBePositive} = 400
        zeta (1,1) double {mustBePositive} = 0.5
        nmLg (1,1) double = 0
    end

    methods (Abstract)
        g = relief(obj, nmx);
        % Permittivity slice at depth nmz, sampled at nmx.
        eps_slice = epsAt(obj, nmx, nmz, nmLm, epsm, epsd, nmda);
    end

    methods (Static)
        function obj = fromConfig(cfg)
            % Factory: pick the right subclass from cfg.type.
            switch cfg.type
                case 0
                    obj = rcwa.gratings.Square('nmLx', cfg.nmLx, ...
                        'zeta', cfg.zeta, 'nmLg', cfg.nmLg);
                case 1
                    obj = rcwa.gratings.Sinusoidal('nmLx', cfg.nmLx, ...
                        'zeta', cfg.zeta, 'nmLg', cfg.nmLg);
                case 2
                    obj = rcwa.gratings.Pyramid('nmLx', cfg.nmLx, ...
                        'zeta', cfg.zeta, 'nmLg', cfg.nmLg);
                otherwise
                    obj = rcwa.gratings.Spherical('nmLx', cfg.nmLx, ...
                        'zeta', cfg.zeta, 'nmLg', cfg.nmLg, ...
                        'nmda', cfg.nmda);
            end
        end
    end

    methods
        function obj = Grating(varargin)
            for i = 1:2:numel(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
        end

        function eps_slice = epsAtCommon(obj, nmx, nmz, nmLm, epsm, epsd)
            % Default implementation shared by Square/Sinusoidal/Pyramid:
            % half-thickness step function compared with the relief curve.
            nmgval = obj.relief(nmx);
            idx = sign(nmz - nmLm - nmgval);
            eps_slice = zeros(size(idx));
            eps_slice(idx ==  1) = epsd;
            eps_slice(idx ==  0) = (epsd + epsm) / 2;
            eps_slice(idx == -1) = epsm;
        end
    end
end
