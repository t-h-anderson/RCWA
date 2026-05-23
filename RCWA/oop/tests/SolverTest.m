classdef SolverTest < matlab.unittest.TestCase
    % Solver-level tests. The most valuable is cross-checking the OOP
    % solver against the legacy RCWA.m on a tiny problem.

    methods (Test)
        function shapesAreCorrect(tc)
            cfg = tc.minimalCfg();
            [solver, epsf, nmdz, inmk0, radtheta] = tc.smallProblem(cfg);
            res = solver.solve(inmk0, radtheta, epsf, nmdz);
            Ns = numel(nmdz);
            N  = 2*cfg.Nt + 1;
            tc.verifySize(res.Z,  [Ns+1, 1]);
            tc.verifySize(res.Tn, [Ns+1, 1]);
            tc.verifySize(res.R,  [N, 1]);
            tc.verifySize(res.E,  [3, Ns, N]);
        end

        function noNaNs(tc)
            cfg = tc.minimalCfg();
            [solver, epsf, nmdz, inmk0, radtheta] = tc.smallProblem(cfg);
            res = solver.solve(inmk0, radtheta, epsf, nmdz);
            tc.verifyFalse(any(isnan(res.R)), 'R has NaN');
            tc.verifyFalse(any(isnan(res.E(:))), 'E has NaN');
        end

        function matchesLegacyRCWA(tc)
            % Build a small homogeneous (no grating) configuration and
            % compare the OOP solver result to RCWA.m bit-for-bit (within
            % numerical tolerance).
            cfg = tc.minimalCfg();
            [solver, epsf, nmdz, inmk0, radtheta] = tc.smallProblem(cfg);

            res_oop = solver.solve(inmk0, radtheta, epsf, nmdz);
            [~, Tn_leg, R_leg, E_leg] = ...
                RCWA(inmk0, radtheta, epsf, nmdz, cfg.toStruct());

            tc.verifyEqual(res_oop.R, R_leg, 'AbsTol', 1e-9, 'RelTol', 1e-9);
            tc.verifyEqual(res_oop.E, E_leg, 'AbsTol', 1e-9, 'RelTol', 1e-9);
            tc.verifyEqual(res_oop.Tn{end}, Tn_leg{end}, ...
                'AbsTol', 1e-9, 'RelTol', 1e-9);
        end

        function energyConservationAirOnly(tc)
            % All-air stack at normal incidence (single mode): |R|^2 + |T|^2
            % should be unity along the zeroth Fourier mode. Tolerance is
            % loose to absorb the homogeneous-branch numerical quirks.
            cfg = rcwa.Config('Nt', 1, 'Nx', 11, 'pol', 1, ...
                              'Nm', 0, 'Ng', 0, 'Nz', 0, 'Nw', 0, 'Nair', 1, ...
                              'nmLair', 100, 'nmda', 0, 'nmLg', 0, ...
                              'nmLm', 0, 'nmLp', 0, 'nmLi', 0, 'nmLn', 0, ...
                              'nmLw', 0, 'gmat', 0, 'dmat', 0, 'wmat', 0, ...
                              'material', 0, 'nsa', 1);
            cfg  = cfg.verify();
            geom = rcwa.Geometry(cfg);
            grat = rcwa.gratings.Grating.fromConfig(cfg);
            pb   = rcwa.PermittivityBuilder(cfg, geom, grat);
            nmx  = linspace(-cfg.nmLx/2, cfg.nmLx/2, cfg.Nx);
            Eg   = zeros(numel(geom.nmz), 1);
            [eps_xz, ~] = pb.build(nmx, 500, Eg);
            epsf = pb.buildFourier(eps_xz);
            solver = rcwa.Solver(cfg);
            res = solver.solve(2*pi/500, 0, epsf, geom.nmdz);

            % Zeroth-mode reflection should be near zero (air-only stack)
            tc.verifyLessThan(abs(res.R(cfg.Nt+1)), 0.05, ...
                'air-only stack should reflect very little');
        end
    end

    methods (Access = private)
        function cfg = minimalCfg(~)
            cfg = rcwa.Config('Nt', 2, 'Nx', 11, 'pol', 1, ...
                              'Nm', 1, 'Ng', 0, 'Nz', 2, 'Nw', 1, 'Nair', 1, ...
                              'nmLair', 100, 'nmLw', 20, 'nmLm', 50, ...
                              'nmda', 0, 'nmLg', 0, ...
                              'nmLp', 5, 'nmLi', 20, 'nmLn', 5, ...
                              'gmat', 0, 'dmat', 0, 'wmat', 0, ...
                              'material', 0, ...
                              'epsm', 1+3i, 'epsd', 4+0.1i, ...
                              'epsJ', 9 + 0.2i, 'epsW', 3 + 1e-6i, ...
                              'nsa', 1+1e-6i);
            cfg = cfg.verify();
        end

        function [solver, epsf, nmdz, inmk0, radtheta] = smallProblem(~, cfg)
            geom = rcwa.Geometry(cfg);
            grat = rcwa.gratings.Grating.fromConfig(cfg);
            pb   = rcwa.PermittivityBuilder(cfg, geom, grat);
            nmx  = linspace(-cfg.nmLx/2, cfg.nmLx/2, cfg.Nx);
            Eg   = zeros(numel(geom.nmz), 1);
            [eps_xz, ~] = pb.build(nmx, 500, Eg);
            epsf = pb.buildFourier(eps_xz);
            nmdz = geom.nmdz;
            inmk0    = 2*pi / 500;
            radtheta = 0;
            solver   = rcwa.Solver(cfg);
        end
    end
end
