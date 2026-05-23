classdef PermittivityBuilder
    % Build epsilon(x,z) and its Fourier representation in x.
    %
    % Replaces BuildEps.m + BuildEpsF.m. Keeps the original conventions
    % faithful so that thesis comparisons still hold. Bugs are documented
    % in BUGS.md and surfaced as TODO comments where relevant.

    properties (SetAccess = immutable)
        cfg     rcwa.Config
        geom    rcwa.Geometry
        % `grating` is intentionally untyped: rcwa.gratings.Grating is
        % abstract, and MATLAB tries to default-construct any typed
        % property -- which fails for abstract classes.
        grating
    end

    methods
        function obj = PermittivityBuilder(cfg, geom, grating)
            obj.cfg     = cfg;
            obj.geom    = geom;
            obj.grating = grating;
        end

        function [eps, cfgOut] = build(obj, nmx, nmlambda, Eg)
            % Build the real-space eps(x,z) matrix and a Config with eps
            % material values resolved at this wavelength.
            %
            % Eg is the bandgap at every junction-region z-slice.

            cfgOut = obj.cfg;
            cfgOut = obj.applyMaterials(cfgOut, nmlambda, Eg);

            nmz = obj.geom.nmz;
            mat = obj.geom.mat_cat;
            eps = zeros(numel(nmx), numel(nmz));

            % Air. NOTE BUG-002: original BuildEps assigned `loc.nsa` here,
            % not `loc.nsa^2`. We preserve that (eps = n, not n^2) to match
            % the original results. Toggle via cfg.air_sq if a fix is wanted.
            eps(:, mat == 4) = cfgOut.nsa;

            % Window. NOTE BUG-003: original reads loc.epsW (capital W) but
            % the glass/AZO branch wrote loc.epsw (lowercase). Here we read
            % the resolved cfgOut.epsW which the applyMaterials step set
            % consistently.
            eps(:, mat == 3) = cfgOut.epsW;

            % Junction.
            if isscalar(cfgOut.epsJ)
                eps(:, mat == 2) = cfgOut.epsJ;
            else
                eps(:, mat == 2) = kron(ones(numel(nmx),1), cfgOut.epsJ);
            end

            % Mirror.
            eps(:, mat == 0) = cfgOut.epsm;

            % Grating region.
            zg = nmz(mat == 1);
            if ~isempty(zg)
                eps(:, mat == 1) = obj.gratingSlices(nmx, zg, cfgOut);
            end
        end

        function epsf = buildFourier(obj, eps_xz)
            % Fourier transform eps along x at every z, giving the 4*Nt+1
            % Fourier coefficients per z-slice. Homogeneous slices keep
            % only the zeroth mode. Choice of explicit vs FFT follows
            % cfg.ftype.

            cfg = obj.cfg;
            nmz = obj.geom.nmz;
            Nt  = cfg.Nt;

            homo = false(numel(nmz), 1);
            for i = 1:numel(nmz)
                homo(i) = all(eps_xz(:, i) == eps_xz(1, i));
            end

            epsf = zeros(numel(nmz), 4*Nt + 1);
            epsf(homo,  2*Nt + 1) = eps_xz(1, homo).';

            if any(~homo)
                if cfg.ftype == 1
                    epsf(~homo, :) = epsn(cfg.toStruct(), nmz(~homo));
                else
                    epsf(~homo, :) = gepszfft(nmz(~homo), cfg.toStruct());
                end
            end
        end
    end

    methods (Access = private)
        function cfg = applyMaterials(obj, cfg, nmlambda, Eg)
            % Resolve dispersive permittivities at this wavelength.

            % Dielectric in the grating
            if cfg.dmat == 1
                m = rcwa.materials.Glass();
                cfg.epsd = m.permittivity(nmlambda);
            elseif cfg.dmat == 2
                m = rcwa.materials.AZO();
                cfg.epsd = m.permittivity(nmlambda);
            end

            % Window. Original had `elseif loc.dmat == 2` in the window
            % block (BUG-004) which made wmat==2 silently use the
            % dielectric's choice. Here we correctly key on wmat.
            if cfg.wmat == 1
                m = rcwa.materials.Glass();
                cfg.epsW = m.permittivity(nmlambda);
            elseif cfg.wmat == 2
                m = rcwa.materials.AZO();
                cfg.epsW = m.permittivity(nmlambda);
            end

            % Metallic grating material
            if cfg.gmat == 1
                m = rcwa.materials.Silver();
                cfg.epsm = m.permittivity(nmlambda);
            end

            % Junction material profile
            if cfg.material == 1 && cfg.Nz ~= 0
                EgJ = Eg(obj.geom.mat_cat == 2);
                m = rcwa.materials.ASiH('eVEg', EgJ(:));
                cfg.epsJ = m.permittivity(nmlambda).';
            elseif cfg.material == 2
                cfg.epsJ = Faryad(obj.geom.nmz, nmlambda);
            end
        end

        function eps_slice = gratingSlices(obj, nmx, zg, cfg)
            eps_slice = zeros(numel(nmx), numel(zg));
            for i = 1:numel(zg)
                eps_slice(:, i) = obj.grating.epsAt( ...
                    nmx, zg(i), cfg.nmLm, cfg.epsm, cfg.epsd, cfg.nmda);
            end
        end
    end
end
