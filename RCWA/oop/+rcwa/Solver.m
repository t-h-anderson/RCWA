classdef Solver
    % Stable RCWA solver for a 1D-periodic stack.
    %
    % Re-implementation of RCWA.m in OOP form. The numerics are kept
    % bit-equivalent (or as close as practical) to the original; bugs are
    % only documented, not silently fixed -- see BUGS.md.
    %
    % Algorithm reference: Mackay, Polo & Lakhtakia,
    % "Electromagnetic Surface Waves" (EMSW), pp. 75-77.
    %
    % Usage
    % -----
    %   solver = rcwa.Solver(cfg);
    %   result = solver.solve(inmk0, radtheta, epsf, nmdz);
    %   result.R, result.Tn, result.Z, result.E
    %
    % result.E is laid out as (component, z-slice, mode), matching the
    % original `E` array.

    properties (SetAccess = immutable)
        cfg rcwa.Config
    end

    methods
        function obj = Solver(cfg)
            obj.cfg = cfg;
        end

        function result = solve(obj, inmk0, radtheta, epsf, nmdz)
            cfg = obj.cfg;
            Ns  = numel(nmdz);
            Nt  = cfg.Nt;
            twoNt1 = 2*Nt + 1;
            fourNt2 = 4*Nt + 2;

            % --- Incident / reflected Y matrices ---
            inmk0x = cfg.nsa*inmk0*sin(radtheta) + ...
                     (-Nt:Nt) * 2*pi / cfg.nmLx;
            inmk0z = sqrt(cfg.nsa^2 * inmk0^2 - inmk0x.^2);

            [Y_inc, Y_ref] = obj.boundaryY(cfg, inmk0, inmk0z);

            Z = cell(Ns+1, 1);
            Z{Ns+1} = Y_inc;

            % --- Per-slice eigen-decomposition pass ---
            Vn   = cell(Ns+1, 1);
            Gn   = cell(Ns+1, 1);
            XU   = cell(Ns+1, 1);

            NullMat = zeros(twoNt1);
            inmKX   = diag(inmk0x);
            epst_temp = []; invepst_temp = [];
            epsz_prev = NaN;

            for n = Ns:-1:1
                epsz = epsf(n, :);

                if ~isequal(epsz_prev, epsz) || n == Ns
                    epsz_prev = epsz;
                    [epst_temp, invepst_temp, V, inmG] = ...
                        obj.layerDecomp(cfg, inmk0, inmk0x, inmKX, ...
                                        NullMat, epsz);
                end

                if obj.isHomogeneous(cfg, epsz)
                    % Original code calls invV*Z explicitly for the
                    % homogeneous branch. layerDecomp already inverted V
                    % analytically and returned it through the sorted V
                    % (we now share the V\\Z path safely).
                    X = V \ Z{n+1};
                else
                    X = V \ Z{n+1};
                end

                Vn{n+1} = V;
                Gn{n+1} = inmG;

                G_U = inmG(1:twoNt1);
                G_L = inmG(twoNt1+1:fourNt2);
                expG_U = diag(exp( 1i * nmdz(n) * G_U));
                expG_L = diag(exp(-1i * nmdz(n) * G_L));

                XU{n+1} = X(1:twoNt1, :);
                U = expG_L * ( X(twoNt1+1:fourNt2, :) * (XU{n+1} \ expG_U) );

                Z{n} = V * [eye(twoNt1); U];
            end

            % --- Solve the boundary system for [T0; R] ---
            A = zeros(twoNt1, 1);
            A(Nt+1) = 1;
            T0R = [Z{1}, -Y_ref] \ (Y_inc * A);
            T0  = T0R(1:twoNt1, :);
            R   = T0R(twoNt1+1:fourNt2, :);

            % --- Propagate the transmission vectors forward ---
            Tn = cell(Ns+1, 1);
            Tn{1} = T0;
            for n = 1:Ns
                G_U = Gn{n+1}(1:twoNt1);
                Tn{n+1} = XU{n+1} \ (diag(exp(1i * nmdz(n) * G_U)) * Tn{n});
            end

            % --- Fields per slice ---
            E = obj.assembleFields(cfg, Ns, Z, Tn, epsf, inmk0, inmKX);

            result = struct('Z', {Z}, 'Tn', {Tn}, 'R', R, 'E', E, ...
                            'inmk0x', inmk0x, 'inmk0z', inmk0z);
        end
    end

    methods (Static, Access = private)
        function tf = isHomogeneous(cfg, epsz)
            % Original test: at least 4*Nt of the 4*Nt+1 Fourier modes
            % are exactly zero, i.e. the slice is laterally uniform.
            tf = sum(epsz == 0) >= 4*cfg.Nt;
        end

        function [Y_inc, Y_ref] = boundaryY(cfg, inmk0, inmk0z)
            N = 2*cfg.Nt + 1;
            if cfg.pol == 1
                Ye_inc = -diag(inmk0z / (cfg.nsa * inmk0));
                Ye_ref = -Ye_inc;
                Yh_inc = -diag(ones(N, 1));
                Yh_ref =  Yh_inc;
            else
                Ye_inc =  cfg.nsa * diag(ones(N, 1));
                Ye_ref =  Ye_inc;
                Yh_inc = -diag(inmk0z / (cfg.nsa * inmk0));
                Yh_ref = -Yh_inc;
            end
            Y_inc = [Ye_inc; Yh_inc];
            Y_ref = [Ye_ref; Yh_ref];
        end

        function [epst, invepst, V, inmG] = layerDecomp(cfg, inmk0, ...
                inmk0x, inmKX, NullMat, epsz)
            % Eigen-decomposition of the slice's wave operator. Two
            % branches: laterally homogeneous (analytic) vs grating
            % (numerical eig).
            Nt = cfg.Nt; N = 2*Nt + 1;

            if rcwa.Solver.isHomogeneous(cfg, epsz)
                % NOTE BUG-005: original computes
                %   epst_temp = epsz(2N+1) * diag(2N+1)
                % where diag(scalar) returns the scalar (a 1x1 matrix),
                % not eye(N). We use eye(N), the obviously intended form.
                epst = epsz(N) * eye(N);
                invepst = (1/epsz(N)) * eye(N);

                if cfg.pol == 1
                    % The factorised form below looks dimensionally odd
                    % but p_bl * p_tr = sqrt(eps*k0^2 - kx^2), the
                    % correct eigenvalue. The trick keeps the eigenvector
                    % matrix V well-scaled for the analytic inverse.
                    p_tr = sqrt(inmk0 - inmk0x.*inmk0x/(epsz(N)*inmk0));
                    p_bl = sqrt(inmk0 * epsz(N));
                else
                    p_tr = sqrt(-inmk0);
                    p_bl = sqrt(-inmk0*epsz(N) + (1./inmk0)*inmk0x.*inmk0x);
                end

                pd_p = diag(p_tr ./ p_bl);
                p_pd = diag(p_bl ./ p_tr);

                V    = [pd_p, -pd_p; eye(N), eye(N)];
                inmG = [p_bl * p_tr, -p_bl * p_tr].';
                inmG = inmG(:);

                [V, order] = sortmat(V, inmG);
                inmG = inmG(order);
            else
                epst = toeplitz(epsz(N:-1:1), epsz(N:4*Nt+1));
                invepst = inv(epst);

                if cfg.pol == 1
                    P14 = inmk0*eye(N) - (1/inmk0)*inmKX*invepst*inmKX;
                    P41 = inmk0 * epst;
                    P = [NullMat, P14; P41, NullMat];
                else
                    P23 = -inmk0*eye(N);
                    P32 = -inmk0*epst + (1./inmk0)*inmKX*inmKX;
                    P = [NullMat, P23; P32, NullMat];
                end

                [V, GD] = eig(P);
                inmG = diag(GD);
                [V, order] = sortmat(V, inmG);
                inmG = inmG(order);
            end
        end

        function E = assembleFields(cfg, Ns, Z, Tn, epsf, inmk0, inmKX)
            Nt = cfg.Nt; N = 2*Nt + 1;
            e_x = zeros(Ns, N);
            e_y = zeros(Ns, N);
            e_z = zeros(Ns, N);
            h_y = zeros(Ns, N);

            for n = 1:Ns
                f = (Z{n+1} * Tn{n+1}).';
                if cfg.pol == 1
                    e_x(n, :) = f(1:N);
                    h_y(n, :) = f(N+1:4*Nt+2);
                    epsz = epsf(n, :);
                    if rcwa.Solver.isHomogeneous(cfg, epsz)
                        invepst = (1/epsz(N)) * eye(N);
                    else
                        invepst = inv(toeplitz(epsz(N:-1:1), epsz(N:4*Nt+1)));
                    end
                    e_z(n, :) = -(1/inmk0) * invepst * (inmKX * h_y(n, :).');
                else
                    e_y(n, :) = f(1:N);
                end
            end

            E = zeros(3, Ns, N);
            E(1, :, :) = e_x;
            E(2, :, :) = e_y;
            E(3, :, :) = e_z;
        end
    end
end
