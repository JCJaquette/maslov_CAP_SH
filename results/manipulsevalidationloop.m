%% Parameter sweep: params.order and params.pulse.order
Case_number = 1;

% Define parameter ranges
order_vals       = [35, 40, 45, 50];
pulse_order_vals = [2^8, 2^9, 2^10, 2^11];

% Flags
BOOL_load_bndl    = 0;
BOOL_load_pulse   = 0;
BOOL_save_pulse   = 0;
BOOL_load_Euminus = 1;
BOOL_save_Euminus = 0;
bndl_BOOL.save_data = 0;

% Base parameters (fixed)
params_base.rho      = .99;
params_base.tol      = 4e-14;
params_base.bd_scale = .2;
params_base.new      = 1.01;
params_base.Eu.order = 2^9;
params_base.mu       = 0.1;
params_base.scale    = .3;
params_base.xi       = 0;    %change this by case
params_base.nu       = 1.6;
params_base.lambda   = 0;
params_base.isIntval = 0;


results = struct();
run_idx = 0;
total_runs = length(order_vals) * length(pulse_order_vals);
fprintf('\n===== Starting parameter sweep: %d total runs =====\n\n', total_runs);

%% Main sweep loop
for i = 1:length(order_vals)

    % Build params with current manifold order
    params            = params_base;
    params.order      = order_vals(i);
    params.mfld.order = params.order;
    if params.isIntval
        params.mu = intval(num2str(params_base.mu));
        params.nu = intval('1.6');
    end

    %% Compute bundle once per manifold order
    fprintf('=== Computing bundle for order=%d ===\n', params.order)
    bundle_ok = false;
    try
        if BOOL_load_bndl
            data_str = "data_bndl_nu_1p6_mu_0p1";
            load(data_str)
            disp('Loaded bundles and manifolds')
        else
            bndl_BOOL.plot       = 0;
            bndl_BOOL.save_image = 0;
            bndl_BOOL.Lminus     = 1;
            bndl_BOOL.stable     = 1;
            [mflds, bndl] = get_all_bundles(params, bndl_BOOL);
        end
        bundle_ok = true;
    catch ME
        warning('  Bundle FAILED for order=%d: %s', params.order, ME.message)
    end

    if ~bundle_ok
        % Log a failed entry for every pulse order in this manifold order
        for j = 1:length(pulse_order_vals)
            run_idx = run_idx + 1;
            results(run_idx).order       = params.order;
            results(run_idx).pulse_order = pulse_order_vals(j);
            results(run_idx).verif       = NaN;
            results(run_idx).r           = NaN;
            results(run_idx).file        = '';
            results(run_idx).error       = ['Bundle failed: ', ME.message];
        end
        fprintf('  Skipping all pulse runs for order=%d\n\n', params.order)
        continue  % skip to next manifold order
    end

    %%
    for j = 1:length(pulse_order_vals)

        run_idx = run_idx + 1;
        params.pulse.order = pulse_order_vals(j);

        fprintf('--- Run %d/%d | order=%d | pulse.order=%d ---\n', ...
            run_idx, total_runs, params.order, params.pulse.order)

        try
            disp('  Computing pulse seed...')
            seed = get_newton_seed(params, mflds);

            disp('  Refining with Newton...')
            pulse4D = refine_cheb_orbit(seed, mflds, params);

            disp('  Validating pulse...')
            [verif, pulse4D.r, vali_data] = verify_homoclinic_orbit( ...
                params, mflds, pulse4D, params.new);

            file_str = ['save_MO', int2str(params.order), ...
                        '_MSc',    num2str(params.scale),  ...
                        '_PO',     int2str(params.pulse.order)];
            z          = flattenstruct(vali_data, '');
            data_Table = struct2table(z, 'AsArray', true);
            save(file_str, 'data_Table')
            fprintf('  Saved: %s\n', file_str)

            results(run_idx).order       = params.order;
            results(run_idx).pulse_order = params.pulse.order;
            results(run_idx).verif       = verif;
            results(run_idx).r           = pulse4D.r;
            results(run_idx).file        = file_str;
            results(run_idx).error       = '';

        catch ME
            warning('  Run %d FAILED: %s', run_idx, ME.message)
            results(run_idx).order       = params.order;
            results(run_idx).pulse_order = params.pulse.order;
            results(run_idx).verif       = NaN;
            results(run_idx).r           = NaN;
            results(run_idx).file        = '';
            results(run_idx).error       = ME.message;
        end

        fprintf('\n')
    end
end