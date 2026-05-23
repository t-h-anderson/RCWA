classdef ConfigTest < matlab.unittest.TestCase
    methods (Test)
        function defaultsMatchOriginalLoc(tc)
            % The OOP defaults are copied from DefaultLoc.m. Verify the
            % public surface matches the legacy struct field-for-field.
            cfg = rcwa.Config();
            loc = DefaultLoc();
            cs  = cfg.toStruct();
            fns = fieldnames(loc);
            for i = 1:numel(fns)
                tc.verifyTrue(isfield(cs, fns{i}), ...
                    sprintf('Config missing field %s', fns{i}));
                tc.verifyEqual(cs.(fns{i}), loc.(fns{i}), ...
                    sprintf('Config.%s differs from DefaultLoc', fns{i}));
            end
        end

        function namedOverrides(tc)
            cfg = rcwa.Config('Nt', 7, 'pol', 0, 'nmLx', 250);
            tc.verifyEqual(cfg.Nt, 7);
            tc.verifyEqual(cfg.pol, 0);
            tc.verifyEqual(cfg.nmLx, 250);
        end

        function rejectsUnknownProperty(tc)
            tc.verifyError(@() rcwa.Config('not_a_real_prop', 1), ...
                'rcwa:Config:unknownProp');
        end

        function rejectsOddPairs(tc)
            tc.verifyError(@() rcwa.Config('Nt'), 'rcwa:Config:badArgs');
        end

        function copyFromStruct(tc)
            s = struct('Nt', 9, 'pol', 0, 'unknown_thing', 'x');
            cfg = rcwa.Config(s);
            tc.verifyEqual(cfg.Nt, 9);
            tc.verifyEqual(cfg.pol, 0);
        end

        function verifyZeroesOutCorrectFields(tc)
            % BUG-008 regression: verify() must zero nmLp/nmLi/nmLn (not
            % the non-existent Lp/Li/Ln that the original wrote).
            cfg = rcwa.Config('Nz', 0, 'nmLp', 10, 'nmLi', 20, 'nmLn', 30);
            cfg = cfg.verify();
            tc.verifyEqual(cfg.nmLp, 0);
            tc.verifyEqual(cfg.nmLi, 0);
            tc.verifyEqual(cfg.nmLn, 0);
        end

        function verifyClampsNmda(tc)
            cfg = rcwa.Config('nmda', 0, 'nmLg', 30);
            cfg = cfg.verify();
            tc.verifyEqual(cfg.nmda, 30);
        end
    end
end
