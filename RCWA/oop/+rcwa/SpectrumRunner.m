classdef SpectrumRunner
    % Sweep the RCWA solver over a wavelength / angle grid.
    %
    % Replaces RCWAsetup.m. Splits the original monolith into clear
    % phases:
    %   1. setup    - build geometry, bandgap profile, grating
    %   2. forward  - per (lambda, theta) call the solver, store fields
    %   3. absorb   - reconstruct fields in real space, integrate Q, G
    %   4. summary  - optical short-circuit current density (mA/cm^2)
    %
    % Constructed with a Config; call `run()` to execute.

    properties (SetAccess = immutable)
        cfg     rcwa.Config
        geom    rcwa.Geometry
        % `grating` is intentionally untyped: rcwa.gratings.Grating is
        % abstract, and MATLAB tries to default-construct any typed
        % property -- which fails for abstract classes.
        grating
        bandgap rcwa.BandgapProfile
        builder rcwa.PermittivityBuilder
        solver  rcwa.Solver
    end

    methods
        function obj = SpectrumRunner(cfg)
            cfg = cfg.verify();
            obj.cfg     = cfg;
            obj.geom    = rcwa.Geometry(cfg);
            obj.grating = rcwa.gratings.Grating.fromConfig(cfg);
            obj.bandgap = rcwa.BandgapProfile.fromConfig(cfg);
            obj.builder = rcwa.PermittivityBuilder(cfg, obj.geom, obj.grating);
            obj.solver  = rcwa.Solver(cfg);
        end

        function out = run(obj)
            cfg  = obj.cfg;
            geom = obj.geom;
            Nt   = cfg.Nt;

            nmlambda = linspace(cfg.nmlambda0, cfg.nmlambda1, cfg.nlambda);
            inmk0    = 2*pi ./ nmlambda;
            radtheta = rcwa.Units.rad_from_deg( ...
                linspace(cfg.degtheta0, cfg.degtheta1, cfg.ntheta));

            nmx = linspace(-0.5*cfg.nmLx, 0.5*cfg.nmLx, cfg.Nx);

            Eg = obj.junctionBandgap();

            Nz = numel(geom.nmz);
            Efl  = zeros(cfg.nlambda, cfg.ntheta, 3, Nz, 2*Nt+1);
            epsl = zeros(cfg.nlambda, cfg.Nx, Nz);
            epsf = zeros(cfg.nlambda, Nz, 4*Nt+1);
            r0 = zeros(cfg.nlambda, cfg.ntheta);
            r1 = zeros(cfg.nlambda, cfg.ntheta);
            t0 = zeros(cfg.nlambda, cfg.ntheta);
            t1 = zeros(cfg.nlambda, cfg.ntheta);

            % --- Phase 2: forward solve ---
            for ll = 1:cfg.nlambda
                [eps_xz, ~] = obj.builder.build(nmx, nmlambda(ll), Eg);
                epsl(ll, :, :) = eps_xz;
                epsf(ll, :, :) = obj.builder.buildFourier(eps_xz);

                for tt = 1:cfg.ntheta
                    res = obj.solver.solve(inmk0(ll), radtheta(tt), ...
                            squeeze(epsf(ll, :, :)), geom.nmdz);
                    Efl(ll, tt, :, :, :) = res.E;
                    if cfg.plotting == 1
                        [r0(ll,tt), r1(ll,tt), t0(ll,tt), t1(ll,tt)] = ...
                            rcwa.Specular.compute(cfg, res.R, ...
                                res.Tn{end}, inmk0(ll), radtheta(tt));
                    end
                end
            end

            % --- Phase 3: absorption + generation ---
            absResult = obj.absorption(nmx, nmlambda, Efl, epsl, Eg);

            % --- Phase 4: summary ---
            nmdJ = obj.junctionSliceThickness();
            nmdx = cfg.nmLx / (cfg.Nx - 1);
            JOpt = rcwa.Constants.q * sum(absResult.G(:)) ...
                   * rcwa.Units.m_from_nm(nmdJ) * nmdx / (10 * cfg.nmLx);

            out = struct( ...
                'nmlambda', nmlambda, 'radtheta', radtheta, 'nmx', nmx, ...
                'nmz', geom.nmz, 'mat_cat', geom.mat_cat, ...
                'Efl', Efl, 'epsl', epsl, 'epsf', epsf, ...
                'r0', r0, 'r1', r1, 't0', t0, 't1', t1, ...
                'Gx', absResult.Gx, 'G', absResult.G, 'Ql', absResult.Ql, ...
                'die_absl', absResult.die_absl, ...
                'metal_absl', absResult.metal_absl, ...
                'mAicm2JOpt', JOpt);
        end

        function nmdJ = junctionSliceThickness(obj)
            cfg = obj.cfg;
            % NOTE BUG-009: original RCWAsetup divides by cfg.Nz without
            % guarding Nz==0, which gives a divide-by-zero. We guard.
            if cfg.Nz == 0
                nmdJ = 0;
            else
                nmdJ = (cfg.nmLp + cfg.nmLi + cfg.nmLn) / cfg.Nz;
            end
        end
    end

    methods (Access = private)
        function Eg = junctionBandgap(obj)
            % NOTE BUG-001 acknowledgement: the EgProfile original treated
            % the i-layer offset wrong. We delegate to BandgapProfile which
            % defaults to legacy=true to reproduce the original numbers.
            cfg  = obj.cfg;
            geom = obj.geom;
            Eg = zeros(numel(geom.nmz), 1);
            zJ = geom.nmz(geom.mat_cat == 2) - geom.nmgtop;
            if isempty(zJ), return; end
            Eg_vals = obj.bandgap.evaluate(zJ);
            % NOTE BUG-010: original used fliplr on the row Eg_vals before
            % assigning into a column slice -- works by coincidence of
            % MATLAB's auto-orient. We use flip(Eg_vals) and explicit reshape.
            Eg(geom.mat_cat == 2) = flip(Eg_vals(:));
        end

        function out = absorption(obj, nmx, nmlambda, Efl, epsl, Eg)
            % Reconstruct fields in real space, build the local
            % absorption density Q(x,z), then the carrier generation rate
            % G(x,z) within the junction layers.

            cfg  = obj.cfg;
            geom = obj.geom;
            Nz   = numel(geom.nmz);
            spec = rcwa.Spectrum();
            W2inmim2Sl = spec.spectralIrradiance(nmlambda);
            Enorm = sqrt(2 * rcwa.Constants.eta0);

            nmdx = cfg.nmLx / (cfg.Nx - 1);
            nmdJ = obj.junctionSliceThickness();
            nmdg = obj.sliceThickness(cfg.nmda,  cfg.Ng);
            nmdm = obj.sliceThickness(cfg.nmLm,  cfg.Nm);
            nmdlambda = (cfg.nmlambda1 - cfg.nmlambda0) / cfg.nlambda;
            nmdlambda = max(nmdlambda, 1);

            nJunction = sum(geom.mat_cat == 2);
            Ql = zeros(cfg.nlambda, cfg.Nx, Nz);
            Q  = zeros(cfg.Nx, Nz);
            Gl = zeros(cfg.nlambda, cfg.Nx, nJunction);
            G  = zeros(cfg.Nx, nJunction);
            die_absl   = zeros(cfg.nlambda, cfg.ntheta);
            metal_absl = zeros(cfg.nlambda, cfg.ntheta);

            for ll = 1:cfg.nlambda
                Q_temp = zeros(1, cfg.Nx, Nz);
                for tt = 1:cfg.ntheta
                    Ef = squeeze(Efl(ll, tt, :, :, :));
                    E_temp = obj.reconstructField(nmx, Ef, Enorm);

                    if cfg.epscalc == 1
                        eps_temp = zeros(cfg.Nx, Nz);
                        epsf_ll = squeeze(obj.builder.buildFourier( ...
                            squeeze(epsl(ll, :, :))));
                        for zz = 1:Nz
                            eps_temp(:, zz) = rcwa.InverseFourier.reconstruct( ...
                                nmx', epsf_ll(zz, :), cfg.nmLx, 0, cfg.Nt);
                        end
                        epsl(ll, :, :) = eps_temp; %#ok<NASGU>
                    end

                    normE2 = abs(E_temp(1,:,:)).^2 + abs(E_temp(2,:,:)).^2 ...
                             + abs(E_temp(3,:,:)).^2;
                    normE2 = squeeze(normE2);

                    Qtt = (pi / (rcwa.Constants.eta0 * ...
                                 rcwa.Units.m_from_nm(nmlambda(ll)))) ...
                          * abs(imag(epsl(ll, :, :))) ...
                          .* reshape(normE2, 1, cfg.Nx, Nz);

                    [die_abst, metal_abst] = obj.layerSums(Qtt, geom.mat_cat, ...
                        cfg, nmdJ, nmdg, nmdm, nmx);
                    die_absl(ll, tt)   = die_abst;
                    metal_absl(ll, tt) = metal_abst;
                    Q_temp = Q_temp + Qtt;
                end
                Q_temp = Q_temp / cfg.ntheta;
                Ql(ll, :, :) = Ql(ll, :, :) + Q_temp;
                Q = Q + squeeze(Ql(ll, :, :));

                if nJunction > 0
                    EgJ = Eg(geom.mat_cat == 2);
                    above = sign(rcwa.Units.eV_from_nm(nmlambda(ll)) - EgJ);
                    gamma = rcwa.Constants.h * rcwa.Constants.c / ...
                            rcwa.Units.m_from_nm(nmlambda(ll));
                    % MATLAB strips trailing singleton dimensions, so
                    % zeros(1, Nx, 1) is actually a 1xNx matrix. squeeze
                    % then yields a row vector that broadcasts against
                    % the (Nx, 1) accumulator G into (Nx, Nx). Force the
                    % shape explicitly through reshape to keep the
                    % accumulation unambiguous for any nJunction.
                    QJ = reshape(Q_temp(1, :, geom.mat_cat == 2), ...
                                 cfg.Nx, nJunction);
                    Gl_this = (W2inmim2Sl(ll) / gamma) * QJ;
                    mask = reshape(above == 1, 1, nJunction);
                    Gl_this = bsxfun(@times, Gl_this, mask);
                    Gl(ll, :, :) = reshape(Gl_this, 1, cfg.Nx, nJunction);
                    G = G + nmdlambda * Gl_this;
                end
            end

            % G integrated over x (per junction z-slice)
            if nJunction > 0
                nmgtop = geom.nmgtop;
                nmLz   = cfg.nmLp + cfg.nmLi + cfg.nmLn;
                zCoord = nmLz - (geom.nmz(geom.mat_cat == 2) - nmgtop);
                xIntG  = sum(G, 1) * nmdx / cfg.nmLx;
                % Force columns to avoid vertcat shape surprises when
                % nJunction == 1 (sum collapses the row to a scalar with
                % a different effective orientation than nmz(...)).
                Gx = [zCoord(:), xIntG(:)];
            else
                Gx = zeros(0, 2);
            end

            out = struct('Ql', Ql, 'Q', Q, 'Gl', Gl, 'G', G, ...
                'die_absl', die_absl, 'metal_absl', metal_absl, 'Gx', Gx);
        end

        function E_temp = reconstructField(obj, nmx, Ef, Enorm)
            cfg = obj.cfg;
            Nz  = numel(obj.geom.nmz);
            E_temp = zeros(3, cfg.Nx, Nz);
            if cfg.pol == 0
                for zz = 1:Nz
                    E_temp(2, :, zz) = Enorm * rcwa.InverseFourier.reconstruct( ...
                        nmx', Ef(2, zz, :), cfg.nmLx, 0, cfg.Nt);
                end
            else
                for zz = 1:Nz
                    E_temp(1, :, zz) = Enorm * rcwa.InverseFourier.reconstruct( ...
                        nmx', Ef(1, zz, :), cfg.nmLx, 0, cfg.Nt);
                    E_temp(3, :, zz) = Enorm * rcwa.InverseFourier.reconstruct( ...
                        nmx', Ef(3, zz, :), cfg.nmLx, 0, cfg.Nt);
                end
            end
        end

        function [die_abst, metal_abst] = layerSums(~, Qtt, mat_cat, ...
                cfg, nmdJ, nmdg, nmdm, nmx)
            mid = numel(nmx);
            die_abst   = sum(sum(Qtt(1, :, mat_cat == 2))) ...
                         * rcwa.Units.m_from_nm(nmdJ) / cfg.Nx ...
                       + sum(sum(Qtt(1, ceil(mid/2):mid, mat_cat == 1))) ...
                         * rcwa.Units.m_from_nm(nmdg) / cfg.Nx;
            metal_abst = sum(sum(Qtt(1, :, mat_cat == 0))) ...
                         * rcwa.Units.m_from_nm(nmdm) / cfg.Nx ...
                       + sum(sum(Qtt(1, 1:floor(mid/2), mat_cat == 1))) ...
                         * rcwa.Units.m_from_nm(nmdg) / cfg.Nx;
        end
    end

    methods (Static, Access = private)
        function dz = sliceThickness(L, N)
            if N > 0, dz = L / N; else, dz = 0; end
        end
    end
end
