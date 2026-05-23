classdef PermittivityBuilderTest < matlab.unittest.TestCase
    methods (Test)
        function homogeneousSlicesHaveOnlyZerothFourierMode(tc)
            % A purely homogeneous stack should land all permittivity into
            % the zeroth Fourier mode (column 2*Nt+1) and zero elsewhere.
            cfg = rcwa.Config('Nt', 3, 'Nx', 21, 'Ng', 0, ...
                              'nmda', 0, 'nmLg', 0, ...
                              'gmat', 0, 'dmat', 0, 'wmat', 0, ...
                              'material', 0, 'Nz', 2, 'Nair', 1, 'Nw', 1);
            geom    = rcwa.Geometry(cfg.verify());
            grating = rcwa.gratings.Grating.fromConfig(cfg);
            pb      = rcwa.PermittivityBuilder(cfg.verify(), geom, grating);
            nmx = linspace(-cfg.nmLx/2, cfg.nmLx/2, cfg.Nx);
            Eg  = zeros(numel(geom.nmz), 1);
            [eps, ~] = pb.build(nmx, 500, Eg);
            epsf = pb.buildFourier(eps);

            zeroth = 2*cfg.Nt + 1;
            other_cols = [1:zeroth-1, zeroth+1:size(epsf,2)];
            tc.verifyTrue(all(epsf(:, other_cols) == 0, 'all'), ...
                'Non-zeroth Fourier modes should vanish for homogeneous stack');
        end

        function gratingZerothModeIsZetaAverage(tc)
            % For a square grating with metal/dielectric, the zeroth Fourier
            % mode is zeta*epsm + (1-zeta)*epsd.
            cfg = rcwa.Config('Nt', 2, 'Nx', 51, 'type', 0, 'zeta', 0.5, ...
                              'Ng', 4, 'nmda', 40, 'nmLg', 40, ...
                              'gmat', 0, 'dmat', 0, 'wmat', 0, ...
                              'material', 0, 'Nz', 1, 'Nair', 1, 'Nw', 1, ...
                              'epsm', 2+0i, 'epsd', 6+0i, 'ftype', 1);
            cfg = cfg.verify();
            geom = rcwa.Geometry(cfg);
            grat = rcwa.gratings.Grating.fromConfig(cfg);
            pb   = rcwa.PermittivityBuilder(cfg, geom, grat);
            nmx = linspace(-cfg.nmLx/2, cfg.nmLx/2, cfg.Nx);
            Eg = zeros(numel(geom.nmz), 1);
            [eps, ~] = pb.build(nmx, 500, Eg);
            epsf = pb.buildFourier(eps);

            % Find a grating slice
            idx = find(geom.mat_cat == 1, 1);
            avg = cfg.zeta * cfg.epsm + (1 - cfg.zeta) * cfg.epsd;
            tc.verifyEqual(epsf(idx, 2*cfg.Nt + 1), avg, 'AbsTol', 1e-9);
        end

        function fourierMatchesLegacy(tc)
            % Cross-check OOP buildFourier against gepszfft for a grating
            % slice with ftype = 0 (FFT path).
            cfg = rcwa.Config('Nt', 3, 'Nx', 51, 'type', 0, 'zeta', 0.4, ...
                              'Ng', 2, 'nmda', 50, 'nmLg', 50, ...
                              'gmat', 0, 'dmat', 0, 'wmat', 0, ...
                              'material', 0, 'Nz', 1, 'Nair', 1, 'Nw', 1, ...
                              'epsm', 1+2i, 'epsd', 4+0i, ...
                              'ftype', 0, 'NFFT', 200);
            cfg = cfg.verify();
            geom = rcwa.Geometry(cfg);
            grat = rcwa.gratings.Grating.fromConfig(cfg);
            pb   = rcwa.PermittivityBuilder(cfg, geom, grat);
            nmx = linspace(-cfg.nmLx/2, cfg.nmLx/2, cfg.Nx);
            Eg = zeros(numel(geom.nmz), 1);
            [eps, ~] = pb.build(nmx, 500, Eg);
            epsf_oop = pb.buildFourier(eps);

            % Legacy gepszfft on the grating slices only
            grat_z = geom.nmz(geom.mat_cat == 1);
            epsf_leg = gepszfft(grat_z, cfg.toStruct());

            % Compare rows aligned to the grating region
            tc.verifyEqual(epsf_oop(geom.mat_cat == 1, :), epsf_leg, ...
                'AbsTol', 1e-9);
        end
    end
end
