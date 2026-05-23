classdef Config
    % Simulation configuration.
    %
    % Replaces the `loc` struct used throughout the original code.
    % Fields are validated on construction; unknown fields raise an error
    % rather than silently being ignored (the original ApplyVarargin
    % erroneously errored only on the second varargin pass, never the
    % struct-merge pass).

    properties
        % Computation
        parforArg  = Inf;       % Number of parallel workers (Inf = default pool)
        plotting   = 0;         % Produce plots (0/1)
        epscalc    = 0;         % Eps exact (0) or Fourier (1)
        pol        = 1;         % Polarisation (0 s, 1 p)
        degtheta0  = 0;         % Min incident angle (deg)
        degtheta1  = 0;         % Max incident angle (deg)
        ntheta     = 1;         % Number of angles
        Nt         = 2;         % Number of Fourier modes (per side)
        nmlambda0  = 300;       % Min wavelength (nm)
        nmlambda1  = 1240/1.6;  % Max wavelength (nm)
        nlambda    = 10;        % Number of wavelengths
        Nx         = 50;        % x samples for reconstruction

        % Grating
        nmLx   = 400;
        zeta   = 0.5;
        type   = 0;             % 0 rect, 1 sine, 2 pyramid, 3 sphere
        nmda   = 0;             % grating region thickness
        nmLg   = 0;             % grating relief
        Ng     = 1;
        nmLm   = 100;           % mirror thickness
        Nm     = 1;
        gmat   = 1;             % 1 silver, else constant
        epsm   = 1 + 3i;
        dmat   = 2;             % 1 glass, 2 AZO, else constant
        epsd   = 11.559 + 0.204i;
        NFFT   = 900;
        ftype  = 1;             % 0 FFT, 1 explicit

        % Junction
        nmLn   = 15;
        nmLi   = 200;
        nmLp   = 15;
        Nz     = 46;
        epsJ   = 11.559 + 0.204i;
        material = 1;           % 0 const, 1 aSiH(GC), 2 Faryad

        % Bandgap profile
        Eg0    = 1.6;
        A      = 0;
        kappa  = 0;
        phi    = 0;
        alpha  = 0;

        % Window
        nmLw   = 75;
        epsW   = 4 + 1e-6i;
        Nw     = 1;
        wmat   = 2;             % 1 glass, 2 AZO, else constant

        % Air
        nmLair = 1000;
        Nair   = 1;
        nsa    = 1 + 1e-6i;     % Refractive index of containing medium
    end

    methods
        function obj = Config(varargin)
            % Construct a Config. Accepts either a struct or name/value pairs,
            % or a leading struct followed by name/value pairs.
            if nargin == 0
                return
            end

            args = varargin;

            % Handle leading struct or Config
            if isstruct(args{1}) || isa(args{1}, 'rcwa.Config')
                src = args{1};
                args = args(2:end);
                names = fieldnames(src);
                for i = 1:numel(names)
                    if isprop(obj, names{i})
                        obj.(names{i}) = src.(names{i});
                    end
                end
            end

            if mod(numel(args), 2) ~= 0
                error('rcwa:Config:badArgs', ...
                    'Config(name, value, ...) pairs expected.');
            end
            for i = 1:2:numel(args)
                name = args{i};
                if ~ischar(name) && ~(isstring(name) && isscalar(name))
                    error('rcwa:Config:badArgs', ...
                        'Property names must be strings.');
                end
                name = char(name);
                if ~isprop(obj, name)
                    error('rcwa:Config:unknownProp', ...
                        '%s is not a recognised parameter name', name);
                end
                obj.(name) = args{i+1};
            end
        end

        function obj = verify(obj)
            % Apply consistency rules from the original VerifyLoc.m.
            % Original wrote loc.Lp/Li/Ln instead of loc.nmLp/nmLi/nmLn (bug
            % BUG-008); here we zero the correct fields.
            if obj.Nair == 0, obj.nmLair = 0; end
            if obj.Nw   == 0, obj.nmLw   = 0; end
            if obj.Nz   == 0
                obj.nmLp = 0; obj.nmLi = 0; obj.nmLn = 0;
            end
            if obj.Ng   == 0
                obj.nmda = 0; obj.nmLg = 0;
            end
            if obj.Nm   == 0, obj.nmLm = 0; end
            if obj.nmda < obj.nmLg
                obj.nmda = obj.nmLg;
            end
        end

        function s = toStruct(obj)
            % Convert to plain struct for interop with the original *.m files.
            names = properties(obj);
            s = struct();
            for i = 1:numel(names)
                s.(names{i}) = obj.(names{i});
            end
        end

    end

    methods (Static)
        function obj = withMisha()
            % Preset matching MishaLoc.m.
            obj = rcwa.Config( ...
                'parforArg', 0, 'plotting', 1, 'pol', 1, ...
                'Nt', 14, 'nmlambda0', 562, 'nmlambda1', 562, ...
                'nlambda', 1, 'Nx', 300, 'nmLx', 500, ...
                'type', 0, 'nmda', 0, 'nmLg', 40, 'Ng', 20, ...
                'nmLm', 150, 'Nm', 20, 'gmat', 1, 'epsm', 0, ...
                'dmat', 0, 'NFFT', 900, 'ftype', 1, ...
                'nmLn', 0, 'nmLi', 275, 'nmLp', 0, 'Nz', 100, ...
                'epsJ', 11.559 + 0.204i, 'material', 0, ...
                'nmLw', 75, 'Nw', 20, 'wmat', 0, ...
                'nmLair', 1000, 'Nair', 100, 'nsa', 1 + 1e-6i);
        end
    end
end
