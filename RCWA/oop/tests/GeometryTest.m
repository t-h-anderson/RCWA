classdef GeometryTest < matlab.unittest.TestCase
    methods (Test)
        function matchesLegacyBuildMaterial(tc)
            % Geometry must reproduce BuildMaterial.m bit-for-bit for the
            % default loc.
            cfg = rcwa.Config('Nm', 3, 'Ng', 5, 'Nz', 10, 'Nw', 2, ...
                              'Nair', 4, 'nmda', 50, 'nmLg', 0);
            g   = rcwa.Geometry(cfg);
            [mat_cat, nmz, nmdz] = BuildMaterial(cfg.toStruct());
            tc.verifyEqual(g.nmz,     nmz,     'AbsTol', 1e-9);
            tc.verifyEqual(g.nmdz,    nmdz,    'AbsTol', 1e-9);
            tc.verifyEqual(g.mat_cat, mat_cat);
        end

        function nmzIsDescending(tc)
            cfg = rcwa.Config();
            g   = rcwa.Geometry(cfg);
            tc.verifyTrue(issorted(g.nmz, 'descend'), ...
                'nmz should run top-of-cell to bottom (decreasing z)');
        end

        function regionCounts(tc)
            cfg = rcwa.Config('Nm', 2, 'Ng', 0, 'Nz', 7, 'Nw', 3, 'Nair', 4);
            cfg.nmda = 0; cfg.nmLg = 0;
            g = rcwa.Geometry(cfg);
            tc.verifyEqual(sum(g.mat_cat == 0), 2);
            tc.verifyEqual(sum(g.mat_cat == 1), 0);
            tc.verifyEqual(sum(g.mat_cat == 2), 7);
            tc.verifyEqual(sum(g.mat_cat == 3), 3);
            tc.verifyEqual(sum(g.mat_cat == 4), 4);
        end

        function totalThicknessSumsCorrectly(tc)
            cfg = rcwa.Config('Nm', 2, 'Ng', 4, 'Nz', 5, 'Nw', 3, 'Nair', 1, ...
                              'nmda', 40, 'nmLg', 20);
            g = rcwa.Geometry(cfg);
            total = cfg.nmLm + cfg.nmda + cfg.nmLp + cfg.nmLi + cfg.nmLn ...
                    + cfg.nmLw + cfg.nmLair;
            tc.verifyEqual(sum(g.nmdz), total, 'RelTol', 1e-12);
        end

        function noRegionWhenZeroSlices(tc)
            cfg = rcwa.Config('Nm', 0, 'Ng', 0, 'Nz', 0, 'Nw', 0, 'Nair', 1);
            cfg = cfg.verify();
            g = rcwa.Geometry(cfg);
            tc.verifyEqual(sum(g.mat_cat == 0), 0);
            tc.verifyEqual(sum(g.mat_cat == 4), 1);
        end
    end
end
