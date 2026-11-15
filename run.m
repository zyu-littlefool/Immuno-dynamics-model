function Results = run(config)
%   group contains burst, sustained, and burst-sustained release modes (45 simulations
%   in total). No retry, qualification screening, folder renaming,
%   collection, or downstream Excel filtering is performed.
%   Optional configuration fields:
%     targetRoot     - Output directory for the selected run folders.
%     helperDir      - Directory containing the required helper .m files.
%     clusterProfile - Parallel cluster profile (default: 'local').
%     nWorkers       - Number of parallel workers (default: 1).
%     seedGroupIDs   - One-based seed groups to run (default: 1:15).
%                      Example: 15 runs the original m = 43, 44, 45 only.
%   Before running, save this file as run.m. The required helper files are
%   banlance_fast.m, Copy_of_New_xiaoshi.m, Copy_4_of_GC.m, and
%   get_Shannon_2.m.
%   This configuration is fixed as: philox4x32_10 + threefry4x64_20, seedOffset = 10400000,
%   globalSeed = 200018, and runPrefix = run.

if nargin < 1 || isempty(config)
    config = struct();
end
config = apply_run_defaults(config);
clc

if ~exist(config.targetRoot, 'dir')
    mkdir(config.targetRoot);
end

helperNames = {
    'banlance_fast.m'
    'Copy_of_New_xiaoshi.m'
    'Copy_4_of_GC.m'
    'get_Shannon_2.m'
};
attachedFiles = cellfun(@(name) fullfile(config.helperDir, name), ...
    helperNames, 'UniformOutput', false);
missingMask = ~cellfun(@(path) exist(path, 'file') == 2, attachedFiles);
if any(missingMask)
    error('Missing helper file(s) in helperDir:\n%s', ...
        strjoin(string(attachedFiles(missingMask)), newline));
end
addpath(config.helperDir);

% Each selected seed group contains three release modes. The original
% run indices are preserved so that seed group 15 produces run43-run45.
nReleaseModes = 3;
mList = [];
for seedGroupID = config.seedGroupIDs
    mList = [mList, (seedGroupID - 1) * nReleaseModes + ...
        (1:nReleaseModes)]; %#ok<AGROW>
end
nSimulations = numel(mList);
runPrefix = 'run';
streamType1 = 'philox4x32_10';
streamType2 = 'threefry4x64_20';
seedOffset = 10400000;
runSeed = 200018;

fprintf('Output root: %s\n', config.targetRoot);
fprintf('Seed groups: %s\n', mat2str(config.seedGroupIDs));
fprintf('Original run indices: %s\n', mat2str(mList));
fprintf('Workers: %d\n', config.nWorkers);
targetRoot = config.targetRoot;

rng(runSeed, 'twister');

%% Model parameters
d_Ag = 1.5;
d_Ag2 = 1.5;
Ag_can = 5e4;
miu = 1.2;
nita = 0.1;
decay_T = 0.23;
alpha_CD4 = 1.5;
decay_base = 0.23;
decay_TFH_can = 0.8;


k_deposit = 24;
d_IC = 0.15;
k_Ig = 0.008;
d_Ig = 0.025;

h = 0.01;
dt = 0.01;
t_max = 8000;
N_T0 = 1200;

%% Initial affinity distributions
N = 1000;
p_less = 0.3;
p_zero = 0.5;
p_neg = 0.19;
p_pos = 0.01;

n_less = round(p_less*N);
n_zero = round(p_zero*N);
n_neg = round(p_neg*N);
n_pos = round(p_pos*N);
n_over_dom = round(0.19*n_pos);
n_over_sub = round(0.72*n_pos);

mu = -0.2;
sigma = 0.1;
n = 190;
r = mu + sigma * randn(1, n);
r = r';

part_1 = -100*ones(1, n_less);
part_2 = zeros(1, n_zero);
part_3 = r(randperm(n_neg))';

mu_2 = 0.6;       % WT affinity change of the dominant epitope
sigma_2 = 0.1;
r_2 = positive_normal(mu_2, sigma_2, n_pos);
part_wt_pos = r_2(randperm(n_pos));

mu_5 = 0.4;       % WT affinity change of the subdominant epitope
sigma_5 = 0.1;
r_5 = positive_normal(mu_5, sigma_5, n_pos);
part_wt_pos_sub = r_5(randperm(n_pos));

mu_3 = 0.1;       % Variant affinity change of the dominant epitope
sigma_3 = 0.1;
r_3 = positive_normal(mu_3, sigma_3, n_pos);
part_var_pos_dom = r_3(randperm(n_pos));
part_var_pos_dom_1 = part_var_pos_dom(1:n_over_dom);
part_var_pos_dom_2 = part_var_pos_dom(n_over_dom+1:n_pos);

mu_4 = 0.4;       % Variant affinity change of the subdominant epitope
sigma_4 = 0.1;
r_4 = positive_normal(mu_4, sigma_4, n_pos);
part_var_pos_sub = r_4(randperm(n_pos));
part_var_pos_sub_1 = part_var_pos_sub(1:n_over_sub);
part_var_pos_sub_2 = part_var_pos_sub(n_over_sub+1:n_pos);

E_wt = [part_wt_pos, part_2, part_3];
E_wt_sub = [part_wt_pos_sub, part_2, part_3];
E_i = [part_1, E_wt];
E_i_2 = [part_1, E_wt_sub];

E_var_dom = [part_var_pos_dom_2, part_2, part_3];
E_var_dom = E_var_dom(randperm(length(E_var_dom)));
E_j = [part_1, part_var_pos_dom_1, E_var_dom];

E_var_sub = [part_var_pos_sub_2, part_2, part_3];
E_var_sub = E_var_sub(randperm(length(E_var_sub)));
E_j_2 = [part_1, part_var_pos_sub_1, E_var_sub];

xE_3 = E_i';
xE_4 = E_j';
xE = E_i_2';
xE_2 = E_j_2';

target_folder_base = parallel.pool.Constant(targetRoot);
D = parallel.pool.DataQueue;
afterEach(D, @(msg) fprintf('[run] %s', msg));


%% Result preallocation
Results(nSimulations) = struct( ...
    'workerFolder', "", ...
    'main_file', "", ...
    'Var_Names', [], ...
    'Data_Matrix', [], ...
    'geneTables', [], ...
    'geneFiles', [], ...
    'm', [], ...
    'seedGroup', [], ...
    'release_mod', [], ...
    'mode', "" ...
);



parfor job = 1:nSimulations
    m = mList(job);


    %% Initialize reproducible random streams
        N_rand = t_max+1;
        stateCell = cell(1, N_rand);
        stateCell_2 = cell(1, N_rand);
        for k = 1:N_rand
            stateCell{k} = cell(1,9);
            stateCell_2{k} = cell(1,9);
        end

        seedGroupForStream = floor((m-0.1)/3);

        for k = 1:N_rand
            for q = 1:9
                stateIndex = (k-1)*9 + q;
                stateSeed1 = seedOffset + seedGroupForStream + ...
                    stateIndex;
                stateSeed2 = seedOffset + 1e5 + ...
                    seedGroupForStream + stateIndex;
                tmpStream = RandStream(streamType1, ...
                    'Seed', stateSeed1);
                stateCell{k}{q} = tmpStream.State;           % Store the stream state
                tmpStream_2 = RandStream(streamType2, ...
                    'Seed', stateSeed2);
                stateCell_2{k}{q} = tmpStream_2.State;
            end
        end

        seedBase = seedOffset + seedGroupForStream;

        stream_1 = RandStream(streamType1, 'Seed', seedBase + 1);
        stream_2 = RandStream(streamType2, ...
            'Seed', seedBase + 1e5 + 1);


        release_mod = mod(m, 3);

        if release_mod == 1       % burst
            k_release_1 = 10;
            k_release_2 = 10;
            t_trans = 0;
        elseif release_mod == 2   % sustained
            k_release_1 = 50;
            k_release_2 = 50;
            t_trans = 0;
        else                      % burst-sustained (release_mod == 0)
            k_release_1 = 10;
            k_release_2 = 50;
            t_trans = 500;
        end

    Ag_remain = zeros(1,t_max+1);
    Ag = zeros(1,t_max+1);
    Ag_decay = zeros(1,t_max+1);
    pC = zeros(1, t_max+1);
    pC_2 = zeros(1, t_max+1);

    T = zeros(1,t_max+1);
    T_2 = zeros(1,t_max+1);
    T_FH = zeros(1,t_max+1);
    T_FH_2 = zeros(1,t_max+1);
    decay_TFH = zeros(1,t_max+1);
    decay_TFH_2 = zeros(1,t_max+1);



    B_act = zeros(1, t_max+1);
    B_act_2 = zeros(1, t_max+1);


    IC = zeros(1,t_max+1);
    IC_2 = zeros(1,t_max+1);
    IC_FDC = zeros(1,t_max+1);
    IC_FDC_2 = zeros(1,t_max+1);
    Ig = zeros(1,t_max+1);
    Ig_2 = zeros(1,t_max+1);
    Ig_titer = zeros(1, t_max+1);
    Ig_titer_2 = zeros(1, t_max+1);


    Ka = zeros(1, t_max+1);
    Ka_2 = zeros(1, t_max+1);
    Ka_sub = zeros(1, t_max+1);
    Ka_2_sub = zeros(1, t_max+1);
    Ka_avg_i = zeros(1, t_max+1);
    Ka_sub_avg_i = zeros(1, t_max+1);
    Ka2_avg_i = zeros(1, t_max+1);
    Ka2_sub_avg_i = zeros(1, t_max+1);


    Ka_PC = zeros(1, t_max+1);
    Ka_PC_2 = zeros(1, t_max+1);
    Ka_PC_sub = zeros(1, t_max+1);
    Ka_PC_2_sub = zeros(1, t_max+1);
    Ka_BM = zeros(1, t_max+1);
    Ka_BM_2 = zeros(1, t_max+1);
    PC = zeros(1,t_max+1);
    PC_2 = zeros(1,t_max+1);
    BM = zeros(1,t_max+1);
    BM_2 = zeros(1,t_max+1);
    C = zeros(1,t_max+1);
    C_2 = zeros(1,t_max+1);

    E_WT_aver = zeros(1, t_max+1);
    E_Var_aver = zeros(1, t_max+1);
    E_WT_sub_aver = zeros(1,t_max+1);
    E_Var_sub_aver = zeros(1, t_max+1);

    E_WT_aver_2 = zeros(1, t_max+1);
    E_Var_aver_2 = zeros(1, t_max+1);
    E_WT_sub_aver_2 = zeros(1,t_max+1);
    E_Var_sub_aver_2 = zeros(1, t_max+1);

    GC_E_WT_aver = zeros(1,t_max+1);
    GC_E_Var_aver = zeros(1, t_max+1);
    GC_E_WT_sub_aver = zeros(1,t_max+1);
    GC_E_Var_sub_aver = zeros(1, t_max+1);

    GC_E_WT_aver_2 = zeros(1,t_max+1);
    GC_E_Var_aver_2 = zeros(1, t_max+1);
    GC_E_WT_sub_aver_2 = zeros(1,t_max+1);
    GC_E_Var_sub_aver_2 = zeros(1, t_max+1);

    all_GC_E_WT_aver = zeros(1,t_max+1);
    all_GC_E_Var_aver = zeros(1, t_max+1);
    all_GC_E_WT_sub_aver = zeros(1,t_max+1);
    all_GC_E_Var_sub_aver = zeros(1, t_max+1);

    all_GC_E_WT_aver_2 = zeros(1,t_max+1);
    all_GC_E_Var_aver_2 = zeros(1, t_max+1);
    all_GC_E_WT_sub_aver_2 = zeros(1,t_max+1);
    all_GC_E_Var_sub_aver_2 = zeros(1, t_max+1);


    Div_GC = zeros(1, t_max+1);
    Div_GC_2 = zeros(1, t_max+1);
    Div_PC = zeros(1, t_max+1);
    Div_PC_2 = zeros(1, t_max+1);

    unique_gene = cell(t_max+1, 4);


    %%
    %对于显性表位和隐性表位采用1和2，分开进行突变和gc反应
    %对于原型株和变异株采用i和j，仅有E有i和j的区分，但是他们是对应于同一个gene的

    N_naive_1 = 150;
    N_naive_2 = 50;
    state_cell_i_1 = zeros(1,N_naive_1);
    state_cell_i_2 = zeros(1,N_naive_2);

    E_cell_j_1 = zeros(1, N_naive_1)+6; % dominant to variant
    E_cell_j_2 = zeros(1, N_naive_2)+6; % subdominant to variant
    E_cell_i_2 = zeros(1,N_naive_2)+3;  % subdominant to WT
    E_cell_i_1 = zeros(1, N_naive_1)+6; % dominant to WT

    E_cell_j_naive_1 = zeros(1, 1500)+6;
    E_cell_j_naive_2 = zeros(1, 500)+6;


    E_cell_naive_1 = zeros(1, 1500)+6;
    E_cell_naive_2 = zeros(1,500)+3;

    gen_cell_i_1 = zeros(N_naive_1,N);
    gen_cell_i_2 = zeros(N_naive_2,N);
    gen_num_1 = ones(1, N_naive_1);
    gen_num_2 = ones(1, N_naive_2);
    del_num_1 = ones(1, N_naive_1);
    del_num_2 = ones(1, N_naive_2);
    div_num_1 = ones(1, N_naive_1);
    div_num_2 = ones(1, N_naive_2);


    N_GCB = zeros(1,t_max+1);
    N_GCB_2 = zeros(1,t_max+1);
    N_GC = ones(1,t_max+1);

    N_GCB_enhance = zeros(1,t_max+1);
    N_GCB_2_enhance= zeros(1,t_max+1);
    PC_enhance= zeros(1,t_max+1);
    PC_2_enhance= zeros(1,t_max+1);

    SHM_1 = zeros(1, t_max+1);
    SHM_2 = zeros(1, t_max+1);
    SHM_1_pos = zeros(1, t_max+1);
    SHM_2_pos = zeros(1, t_max+1);
    SHM_1_pos_j = zeros(1, t_max+1);
    SHM_2_pos_j = zeros(1, t_max+1);

    Ag_remain(1) = 3;
    Ag(1) = 0;
    T(1) = 28;
    T_FH(1) = 0;
    T_2(1) = 28;
    T_FH_2(1) = 0;
    Ig(1) = 0;
    Ig_2(1) = 0;
    IC(1) = 0;
    IC_2(1) = 0;
    Ka(1) = 1e-3;
    Ka_2(1) = 1e-3;
    Ka_sub(1) = 1e-6;
    Ka_2_sub(1) = 1e-3;
    IC_FDC(1) = 0;
    IC_FDC_2(1) = 0;
    Ka_PC(1) = 1e-3;
    Ka_PC_2(1) = 1e-3;
    Ka_PC_sub(1) = 1e-6;
    Ka_PC_2_sub(1) = 1e-3;
    PC(1) = 0;
    PC_2(1) = 0;


    Ag_over_1 = 0;
    Ag_over_2 = 0;

    X_prev = zeros(4, 1);

    % 记录主模型每一步释放后的抗原，供 CD8T/TEFF 使用。
    CD8_AntigenSignal = zeros(1, t_max + 1);
    CD8_AntigenSignal(1) = Ag(1);


    for i = 1:t_max  %随着时间进行


        if mod(i,500)==0                    % 每 500 步汇报一次
                send(D, sprintf('worker%d: inner %d/%d\n', m, i, t_max));
         end


    if i < (t_trans+1)
        k_release_ag = Ag_remain(1)/k_release_1;
    else
        k_release_ag = Ag_remain(1)/k_release_2;
    end

    %抗原释放
    if Ag_remain(i)>0
        Ag_remain(i+1) = Ag_remain(i) - h*k_release_ag;
        New_Ag = Ag(i) + h*k_release_ag;
        Ag_pre = IC(i) + New_Ag;
    else
        New_Ag = Ag(i);
        Ag_pre = IC(i) + New_Ag;
        Ag_over_1 = 1;
        Ag_over_2 = 1;
    end


    %抗原衰减

    New_Ag_decay = Ag_decay(i) + h*d_Ag*New_Ag;
    Ag_pre_2 = New_Ag_decay + IC_2(i);


     % 新增：只记录，不修改主模型的抗原。
    CD8_AntigenSignal(i+1) = Ag_pre + Ag_pre_2;

    pC(i+1) = Ag_pre * Ag_can;
    pC_2(i+1) = Ag_pre_2* Ag_can;

    C_dom = min(pC(i+1),5e40);
    C_sub = min(pC_2(i+1),5e40);


    pC(i+1) = pC(i+1) - h * miu * pC(i+1);
    pC_2(i+1) = pC_2(i+1) - h * miu * pC_2(i+1);

    T(i+1) = T(i) + h*(alpha_CD4 * T(i) * C_dom/(T(i)  +C_dom))-h*nita*T(i);
    T_2(i+1) = T_2(i) + h*(alpha_CD4 * T_2(i) * C_sub/ (T_2(i) + C_sub))-h*nita*T_2(i);
    T_FH(i+1) = T_FH(i) + h*(nita*T(i));
    T_FH_2(i+1) = T_FH_2(i) + h*(nita*T_2(i));

    decay_TFH(i) = decay_base + decay_TFH_can * (1 - Ag_pre/(Ag_pre+0.01));
    decay_TFH_2(i) = decay_base + decay_TFH_can * (1 - Ag_pre_2/(Ag_pre_2+0.01));

    T(i+1) = T(i+1) - h*decay_T*T(i+1);
    T_FH(i+1) = T_FH(i+1) - h*decay_TFH(i)*T_FH(i+1);
    T_2(i+1) = T_2(i+1) - h*decay_T*T_2(i+1);
    T_FH_2(i+1) = T_FH_2(i+1) - h*decay_TFH_2(i)*T_FH_2(i+1);


    %% 抗体反应

    New_Ig = h*(k_Ig*PC(i));
    New_Ig_2 = h*(k_Ig*PC_2(i));

    % 1. 合并抗体并更新平均亲和力 (加权平均)
    Ig_tot = Ig(i) + New_Ig;
    % 更新 Ka_tot (主) 和 Ka_sub_tot (亚)\

    if Ig_tot > 1e-18 % 防止除以0
        Ka_avg = (Ig(i)*Ka_avg_i(i) + New_Ig*Ka_PC(i)) / Ig_tot;
        Ka_sub_avg = (Ig(i)*Ka_sub_avg_i(i) + New_Ig*Ka_PC_sub(i)) / Ig_tot;
    else
        Ka_avg = Ka(1);
        Ka_sub_avg = Ka_sub(1);
    end

    Ka_avg_i(i+1) = Ka_avg;
    Ka_sub_avg_i(i+1) = Ka_sub_avg;

    Ig2_tot = Ig_2(i) + New_Ig_2;
    if Ig2_tot > 1e-18
        Ka2_avg = (Ig_2(i)*Ka2_avg_i(i) + New_Ig_2*Ka_PC_2(i)) / Ig2_tot;
        Ka2_sub_avg = (Ig_2(i)*Ka2_sub_avg_i(i) + New_Ig_2*Ka_PC_2_sub(i)) / Ig2_tot;
    else
        Ka2_avg = Ka_2(1);
        Ka2_sub_avg = Ka_2_sub(1);
    end
    Ka2_avg_i(i+1) = Ka2_avg;
    Ka2_sub_avg_i(i+1) = Ka2_sub_avg;

    % 2. 求解四元平衡
    % 传入初值 [x11; x12; x21; x22] 实现热启动
    [x11, x12, x21, x22] = banlance_fast(New_Ag, New_Ag_decay, ...
                            Ig_tot, Ig2_tot, ...
                            Ka_avg, Ka_sub_avg, ...
                            Ka2_avg, Ka2_sub_avg, ...
                            X_prev);

    % 3. 更新状态
    % 抗体池减去结合掉的总量
    Ig(i+1) = (Ig_tot - x11 - x12) * (1 - h*d_Ig);
    Ig_2(i+1) = (Ig2_tot - x21 - x22) * (1 - h*d_Ig);

    % 抗原池更新
    Ag_free = New_Ag - x11 - x21;
    Ag_decay_free = New_Ag_decay - x12 - x22;
    Ag(i+1) = Ag_free * (1 - h*d_Ag);
    Ag_decay(i+1) = Ag_decay_free * (1 - h*d_Ag2 );

    % 保存解供下一轮热启动
    X_prev = [x11; x12; x21; x22];

    % 记录免疫复合物产量 (用于后续信号传递)
    New_IC_1 = IC(i)+ x11 + x12;
    New_IC_2 = IC_2(i) + x21 + x22;


    New_IC_FDC_1 = IC_FDC(i) + h*k_deposit * New_IC_1;
    New_IC_FDC_2 = IC_FDC_2(i) + h*k_deposit * New_IC_2;

    IC_FDC(i+1) = New_IC_FDC_1 - h*d_IC*New_IC_FDC_1;
    IC_FDC_2(i+1) = New_IC_FDC_2 - h*d_IC*New_IC_FDC_2;

    IC(i+1) = New_IC_1 - h* (k_deposit + d_IC) * New_IC_1;
    IC_2(i+1) = New_IC_2 - h* (k_deposit + d_IC) * New_IC_2;

    Ig_titer(i+1) = Ig(i+1)*Ka_avg;
    Ig_titer_2(i+1) = Ig_2(i+1) * Ka2_sub_avg;


    [state_cell_i_1, E_cell_i_1, E_cell_j_1, gen_cell_i_1, gen_num_1, del_num_1, div_num_1] = ...
    Copy_of_New_xiaoshi(state_cell_i_1, gen_num_1, gen_cell_i_1, E_cell_i_1, E_cell_j_1, del_num_1, div_num_1, stream_1, stateCell{i});
    [state_cell_i_2, E_cell_i_2, E_cell_j_2, gen_cell_i_2, gen_num_2, del_num_2, div_num_2] = ...
    Copy_of_New_xiaoshi(state_cell_i_2, gen_num_2, gen_cell_i_2, E_cell_i_2, E_cell_j_2, del_num_2, div_num_2, stream_2, stateCell_2{i});


    %% GC反应

    C(i+1) = 0.01*(Ag(i+1)+IC(i+1))+(IC_FDC(i+1));
    C_2(i+1) = 0.01*(Ag_decay(i+1)+IC_2(i+1))+(IC_FDC_2(i+1));


    N_GC(i+1) = max(round((T_FH(i+1)+T_FH_2(i+1))/N_T0),1);
    if i > 1
        if N_GC(i+1) > N_GC(i)
            addGC = N_GC(i+1) - N_GC(i);
            addB_1 = addGC * length(E_cell_naive_1);
            addB_2 = addGC * length(E_cell_naive_2);
        else
            addGC = 0;
            addB_1 = 0;
            addB_2 = 0;
        end
    else
        addGC = 0;
        addB_1 = 0;
        addB_2 = 0;
    end

    % dominate
    state_cell_i_1 = [state_cell_i_1, zeros(1,addB_1)];
    E_cell_i_1 = [E_cell_i_1, repmat(E_cell_naive_1, 1, addGC)]; % WT
    E_cell_j_1 = [E_cell_j_1, repmat(E_cell_j_naive_1, 1, addGC)]; % Variant
    gen_cell_i_1 = [gen_cell_i_1; zeros(addB_1, N)];
    gen_num_1 = [gen_num_1, ones(1, addB_1)];
    del_num_1 = [del_num_1, ones(1, addB_1)];
    div_num_1 = [div_num_1, zeros(1, addB_1)];

    % subdominate
    state_cell_i_2 = [state_cell_i_2, zeros(1,addB_2)];
    E_cell_i_2 = [E_cell_i_2, repmat(E_cell_naive_2, 1, addGC)]; % WT
    E_cell_j_2 = [E_cell_j_2, repmat(E_cell_j_naive_2, 1, addGC)]; % Variant
    gen_cell_i_2 = [gen_cell_i_2; zeros(addB_2, N)];
    gen_num_2 = [gen_num_2, ones(1, addB_2)];
    del_num_2 = [del_num_2, ones(1, addB_2)];
    div_num_2 = [div_num_2, zeros(1, addB_2)];


    [state_cell_i_1, E_cell_i_1, E_cell_j_1 , gen_cell_i_1, gen_num_1,del_num_1, SHM_1(i+1), SHM_1_pos(i+1), SHM_1_pos_j(i+1),div_num_1] = ...
        Copy_4_of_GC(state_cell_i_1, E_cell_i_1, E_cell_j_1, gen_cell_i_1, gen_num_1,del_num_1, C(i+1), N_GC(i+1), T(i+1), T_FH(i+1), xE_3, xE_4, SHM_1(i), SHM_1_pos(i), SHM_1_pos_j(i),div_num_1,stream_1,stateCell{i},Ag_over_1);

    [state_cell_i_2, E_cell_i_2, E_cell_j_2, gen_cell_i_2, gen_num_2, del_num_2, SHM_2(i+1), SHM_2_pos(i+1), SHM_2_pos_j(i+1),div_num_2] = ...
         Copy_4_of_GC(state_cell_i_2, E_cell_i_2, E_cell_j_2, gen_cell_i_2, gen_num_2, del_num_2, C_2(i+1), N_GC(i+1), T_2(i+1), T_FH_2(i+1), xE, xE_2, SHM_2(i), SHM_2_pos(i), SHM_2_pos_j(i), div_num_2,stream_2,stateCell_2{i},Ag_over_2);


    gc_state_1 = state_cell_i_1 == 2 | ...
             state_cell_i_1 == 10 | ...
             state_cell_i_1 == 12;

    gc_state_2 = state_cell_i_2 == 2 | ...
                 state_cell_i_2 == 10 | ...
                 state_cell_i_2 == 12;

    B_act(i+1) = length(state_cell_i_1(state_cell_i_1 == 1));
    B_act_2(i+1) = length(state_cell_i_2(state_cell_i_2 == 1));


    PC(i+1) = length(find(state_cell_i_1 == 5));
    PC_2(i+1) = length(find(state_cell_i_2 == 5));

    if PC(i+1) ~= 0
        Ka_PC(i+1) = mean(10.^(E_cell_i_1(state_cell_i_1 == 5)-9));
        Ka_PC_sub(i+1) = mean(10.^(E_cell_j_1(state_cell_i_1   == 5)-9));
        if PC_2(i+1) ~= 0
            Ka_PC_2(i+1) = mean(10.^(E_cell_i_2(state_cell_i_2 == 5)-9));
            Ka_PC_2_sub(i+1) = mean(10.^(E_cell_j_2(state_cell_i_2 == 5)-9));
        else
            Ka_PC_2(i+1) = Ka_PC_2(i);
            Ka_PC_2_sub(i+1) = Ka_PC_2_sub(i);
        end
    else
        Ka_PC(i+1) = Ka_PC(i);
        Ka_PC_sub(i+1) = Ka_PC_sub(i);
        if PC_2(i+1) ~= 0
            Ka_PC_2(i+1) = mean(10.^(E_cell_i_2(state_cell_i_2 == 5)-9));
            Ka_PC_2_sub(i+1) = mean(10.^(E_cell_j_2(state_cell_i_2 == 5)-9));
        else
            Ka_PC_2(i+1) = Ka_PC_2(i);
            Ka_PC_2_sub(i+1) = Ka_PC_2_sub(i);
        end
    end


    E_WT_aver_2(i+1) = log10(mean(10.^(E_cell_i_1(state_cell_i_1 == 5))));
    E_WT_sub_aver_2(i+1) = log10(mean(10.^(E_cell_i_2(state_cell_i_2 == 5))));
    E_Var_aver_2(i+1) =  log10(mean(10.^(E_cell_j_1(state_cell_i_1 == 5))));
    E_Var_sub_aver_2(i+1) = log10(mean(10.^(E_cell_j_2(state_cell_i_2 == 5))));

    all_GC_E_WT_aver_2(i+1) = log10(mean(10.^(E_cell_i_1(state_cell_i_1 == 2))));
    all_GC_E_WT_sub_aver_2(i+1) = log10(mean(10.^(E_cell_i_2(state_cell_i_2 == 2))));
    all_GC_E_Var_aver_2(i+1) = log10(mean(10.^(E_cell_j_1(state_cell_i_1 == 2))));
    all_GC_E_Var_sub_aver_2(i+1) = log10(mean(10.^(E_cell_j_2(state_cell_i_2 == 2))));


    % GC B细胞总数
    N_GCB(i+1) = nnz(gc_state_1);
    N_GCB_2(i+1) = nnz(gc_state_2);

    % 亲和力增强的GC B细胞
    N_GCB_enhance(i+1) = nnz(gc_state_1 & E_cell_i_1 > 6);
    N_GCB_2_enhance(i+1) = nnz(gc_state_2 & E_cell_j_2 > 6);
    PC_enhance(i+1) =  length(find(state_cell_i_1 == 5 & E_cell_i_1 > 6));
    PC_2_enhance(i+1) = length(find(state_cell_i_2 == 5 & E_cell_j_2 > 6));


    if N_GCB(i+1) ~= 0
        [Div_GC(i+1),unique_gene{i+1,1}] = get_Shannon_2(gen_cell_i_1(gc_state_1,:), xE_3, xE_4);
    end
    if N_GCB_2(i+1) ~= 0
        [Div_GC_2(i+1),unique_gene{i+1,2}] = get_Shannon_2(gen_cell_i_2(gc_state_2,:), xE, xE_2);
    end

    if  PC(i+1) ~= 0
        [Div_PC(i+1),unique_gene{i+1,3}] = get_Shannon_2(gen_cell_i_1(state_cell_i_1 == 5,:), xE_3, xE_4);
    end
    if PC_2(i+1) ~= 0
        [Div_PC_2(i+1),unique_gene{i+1,4}] = get_Shannon_2(gen_cell_i_2(state_cell_i_2 == 5,:), xE, xE_2);
    end


    end

     %% Sample outputs and align the time axis


    output = 0:20:t_max;
    output = output + 1; % 转换为 MATLAB 从 1 开始的索引
    time_out = (output - 1) * dt; % dt 是你的步长，确保时间轴计算准确
    time_out = time_out';

    % 统一提取抽样数据
    Ag_remain_out = Ag_remain(output)';
    Ag_out = Ag(output)';
    Ag_decay_out = Ag_decay(output)';
    C_out = C(output)';
    C_2_out = C_2(output)';
    IC_FDC_out = IC_FDC(output)';
    IC_FDC_2_out = IC_FDC_2(output)';
    C_dom_out = pC(output)';
    C_sub_out = pC_2(output)';
    N_GCB_out = N_GCB(output)';
    N_GCB_2_out = N_GCB_2(output)';
    T_out = T(output)';
    T_2_out = T_2(output)';
    T_FH_out = T_FH(output)';
    T_FH_2_out = T_FH_2(output)';
    Ig_titer_out = Ig_titer(output)';
    Ig_titer_2_out = Ig_titer_2(output)';
    B_act_out = B_act(output)';
    B_act_2_out = B_act_2(output)';

    E_WT_aver_2_out = E_WT_aver_2(output)';
    E_WT_sub_aver_2_out = E_WT_sub_aver_2(output)';
    E_Var_aver_2_out = E_Var_aver_2(output)';
    E_Var_sub_aver_2_out = E_Var_sub_aver_2(output)';

    all_GC_E_WT_aver_2_out = all_GC_E_WT_aver_2(output)';
    all_GC_E_WT_sub_aver_2_out = all_GC_E_WT_sub_aver_2(output)';
    all_GC_E_Var_aver_2_out = all_GC_E_Var_aver_2(output)';
    all_GC_E_Var_sub_aver_2_out = all_GC_E_Var_sub_aver_2(output)';


    SHM_1_out = SHM_1(output)';
    SHM_2_out = SHM_2(output)';
    SHM_1_pos_out = SHM_1_pos(output)';
    SHM_2_pos_out = SHM_2_pos(output)';
    SHM_1_pos_j_out = SHM_1_pos_j(output)';
    SHM_2_pos_j_out = SHM_2_pos_j(output)';
    PC_out = PC(output)';
    PC_2_out = PC_2(output)';
    BM_out = BM(output)';
    BM_2_out = BM_2(output)';
    Ig_out = Ig(output)';
    Ig_2_out = Ig_2(output)';
    Div_PC_out = Div_PC(output)';
    Div_PC_2_out = Div_PC_2(output)';
    Div_GC_out = Div_GC(output)';
    Div_GC_2_out = Div_GC_2(output)';

    Ka_avg_i_out = Ka_avg_i(output)';
    Ka_sub_avg_i_out =  Ka_sub_avg_i(output)';
    Ka2_avg_i_out = Ka2_avg_i(output)';
    Ka2_sub_avg_i_out =  Ka2_sub_avg_i(output)';

    N_GCB_ALL = N_GCB_out + N_GCB_2_out;
    T_ALL = T_out + T_2_out;
    T_FH_ALL = T_FH_out + T_FH_2_out;
    B_act_ALL = B_act_out + B_act_2_out;
    PC_ALL = PC_out + PC_2_out;

    N_GCB_enhance_out = N_GCB_enhance(output)';
    N_GCB_2_enhance_out = N_GCB_2_enhance(output)';
    PC_enhance_out = PC_enhance(output)';
    PC_2_enhance_out = PC_2_enhance(output)';


    % Gene snapshot indices for subdominant GC and PC populations.
    [~, rel_idx2] = max(N_GCB_2);
    [~, rel_idx4] = max(PC_2);

    %% Construct the main output matrix
    Var_Names = ["time","Ag_remain","Ag_decay","Ag","C","C_2","IC_FDC","IC_FDC_2","C_dom","C_sub","N_GCB_all","T_all","T_FH_all","B_Act_all","PC_all",...
                 "N_GCB","N_GCB_2","T","T_2","T_FH","T_FH_2","B_act","B_act_2",...
                 "Ig","Ig_2","Ig_titer","Ig_titer_2","Ka_avg","Ka_sub_avg","Ka2_avg","Ka2_sub_avg",...
                 "log_E_WT_avg","log_E_WT_sub_avg","log_E_Var_avg","log_E_Var_sub_avg",...
                 "log_GC_WT_avg","log_GC_WT_sub_avg","log_GC_Var_avg","log_GC_Var_sub_avg",...
                 "GC_+E","GCSUB_+E","PC+E","PCSUB+E","SHM_1","SHM_2",...
                 "PC","PC_2","Div_PC","Div_PC2","Div_GC","Div_GC2"];

    Data_Matrix = [time_out, Ag_remain_out, Ag_decay_out, Ag_out,C_out, C_2_out,IC_FDC_out,IC_FDC_2_out,C_dom_out,C_sub_out, N_GCB_ALL,T_ALL, T_FH_ALL,B_act_ALL, PC_ALL,...
                   N_GCB_out, N_GCB_2_out, T_out, T_2_out, T_FH_out, T_FH_2_out, B_act_out, B_act_2_out, ...
                   Ig_out, Ig_2_out, Ig_titer_out, Ig_titer_2_out, Ka_avg_i_out, Ka_sub_avg_i_out, Ka2_avg_i_out, Ka2_sub_avg_i_out, ...
                   E_WT_aver_2_out, E_WT_sub_aver_2_out, E_Var_aver_2_out, E_Var_sub_aver_2_out,...
                   all_GC_E_WT_aver_2_out, all_GC_E_WT_sub_aver_2_out, all_GC_E_Var_aver_2_out, all_GC_E_Var_sub_aver_2_out, ...
                   N_GCB_enhance_out, N_GCB_2_enhance_out,PC_enhance_out, PC_2_enhance_out,SHM_1_out, SHM_2_out, ...
                   PC_out, PC_2_out, Div_PC_out, Div_PC_2_out, Div_GC_out, Div_GC_2_out];

    %% Write run outputs

    switch  release_mod
        case 0
            mode = "burst-sustained";
        case 1
            mode = "burst";
        case 2
            mode = "sustained";
    end


    workerFolder = fullfile(target_folder_base.Value, ...
        sprintf('%s%02d', runPrefix, m));

    if ~exist(workerFolder, 'dir')
        mkdir(workerFolder);
    end

   AntigenSignalRaw = CD8_AntigenSignal(:);
   sampleIndex = output(:);
   signalSource = "Ag_pre_plus_Ag_pre_2";

   cacheFile = fullfile(workerFolder, 'CD8_input_cache.mat');

   save(cacheFile, ...
       'AntigenSignalRaw', ...  % Full-resolution CD8 antigen input
       'Ag', ...                % Main-model free antigen trajectory
       'dt', ...                % Time step
       'sampleIndex', ...       % Indices used for sampled outputs
       'mode', ...              % Release mode
       'signalSource', ...      % Signal definition for provenance
       '-v7');

    % 主文件名包含 m，防止并行任务互相覆盖。
    main_file = fullfile(workerFolder, sprintf('data_main_%s_seed%d_m%02d.csv', mode, seedBase, m));

    % --- 收集 gene 表以及对应文件名 ---
    gene_types = ["GC_dom", "GC_sub", "PC_dom", "PC_sub"];

    % 关键点：rel_idx2 对应第 2 列 GC_sub；rel_idx4 对应第 4 列 PC_sub
    number_raw = [rel_idx2, rel_idx4];
    key_cols   = [2, 4];

    % 额外关键点
    extra_points = [1401, 2101, 2801];

    % 预分配：关键点 2 个 + extra_points 每个 4 类
    max_items = numel(number_raw) + numel(extra_points) * 4;
    geneTables = cell(max_items, 1);
    geneFiles  = cell(max_items, 1);
    valid_count = 0;

    %% 1) 先收集 rel_idx2 / rel_idx4 这两个关键点
    for kk = 1:numel(number_raw)

        curr_idx = number_raw(kk);
        curr_col = key_cols(kk);

        if curr_idx <= size(unique_gene, 1) && curr_col <= size(unique_gene, 2)

            curr_gene_data = unique_gene{curr_idx, curr_col};

            if ~isempty(curr_gene_data) && isstruct(curr_gene_data)

                T_gene = table((1:curr_gene_data.seednum)', ...
                    curr_gene_data.count(:), ...
                    curr_gene_data.frequen(:), ...
                    curr_gene_data.E_chenge_wt(:), ...
                    curr_gene_data.E_chenge_var(:), ...
                    curr_gene_data.char(:), ...
                    'VariableNames', {'PatternID', 'Count', 'Frequency', ...
                                      'E_change_wt', 'E_change_var', 'char'});

                valid_count = valid_count + 1;
                geneTables{valid_count, 1} = T_gene;

                geneFiles{valid_count, 1} = fullfile(workerFolder, ...
                    sprintf('Gene_%s_at_%d_relmod%d_seed%d_mIter%02d.csv', ...
                    gene_types(curr_col), curr_idx, release_mod, seedBase, m));
            end
        end
    end

    %% 2) 再收集 extra_points
    for ep = extra_points

        if ep <= size(unique_gene, 1)

            for k = 1:4

                if k <= size(unique_gene, 2)

                    curr_gene_data = unique_gene{ep, k};

                    if ~isempty(curr_gene_data) && isstruct(curr_gene_data)

                        T_gene = table((1:curr_gene_data.seednum)', ...
                            curr_gene_data.count(:), ...
                            curr_gene_data.frequen(:), ...
                            curr_gene_data.E_chenge_wt(:), ...
                            curr_gene_data.E_chenge_var(:), ...
                            curr_gene_data.char(:), ...
                            'VariableNames', {'PatternID', 'Count', 'Frequency', ...
                                              'E_change_wt', 'E_change_var', 'char'});

                        valid_count = valid_count + 1;
                        geneTables{valid_count, 1} = T_gene;

                        geneFiles{valid_count, 1} = fullfile(workerFolder, ...
                            sprintf('Gene_%s_at_%d_%s_relmod%d_seed%d_mIter%02d.csv', ...
                            char(mode), ep, gene_types(k), release_mod, seedBase, m));
                    end
                end
            end
        end
    end

    % 去掉空项
    geneTables = geneTables(1:valid_count);
    geneFiles  = geneFiles(1:valid_count);

    % --- 写入 Results(m) ---
    Results(job).workerFolder = string(workerFolder);
    Results(job).main_file = string(main_file);
    Results(job).Var_Names = Var_Names;
    Results(job).Data_Matrix = Data_Matrix;
    Results(job).geneTables = geneTables;
    Results(job).geneFiles = string(geneFiles);
    Results(job).m = m;
    Results(job).seedGroup = floor((m - 0.1) / 3) + 1;
    Results(job).release_mod = release_mod;
    Results(job).mode = mode;

    tf = Results(job).workerFolder;

    if strlength(tf) == 0
        continue
    end

    if ~exist(tf, 'dir')
        mkdir(tf);
    end

    % 1) 写主矩阵 CSV
    header = cellstr(Results(job).Var_Names(:)');  % 强制横向表头
    main_file = Results(job).main_file;

    writecell(header, main_file);
    writematrix(Results(job).Data_Matrix, main_file, 'WriteMode', 'append');

    % 2) 写 gene CSV
    geneFiles = Results(job).geneFiles;
    geneTables = Results(job).geneTables;

    for j = 1:numel(geneTables)
        try
            writetable(geneTables{j}, geneFiles(j));
        catch ME
            warning("Failed writing %s : %s", char(geneFiles(j)), ME.message);
        end
    end


end % parfor

clear target_folder_base

fprintf(['Finished: all %d simulations were exported ' ...
    'without screening to %s.\n'], nSimulations, config.targetRoot);
end

function config = apply_run_defaults(config)
% Apply and validate run-configuration settings without changing the model.

defaults = struct();
defaults.targetRoot = fullfile(fileparts(mfilename('fullpath')), 'output');
defaults.helperDir = fileparts(mfilename('fullpath'));
defaults.clusterProfile = 'local';
defaults.nWorkers = 1;
defaults.seedGroupIDs = 1:15;

defaultNames = fieldnames(defaults);
for k = 1:numel(defaultNames)
    name = defaultNames{k};
    if ~isfield(config, name) || isempty(config.(name))
        config.(name) = defaults.(name);
    end
end

config.targetRoot = char(config.targetRoot);
config.helperDir = char(config.helperDir);
config.clusterProfile = char(config.clusterProfile);

validateattributes(config.seedGroupIDs, {'numeric'}, ...
    {'nonempty', 'vector', 'integer', '>=', 1, '<=', 15});
config.seedGroupIDs = unique(config.seedGroupIDs(:).');
validateattributes(config.nWorkers, {'numeric'}, ...
    {'nonempty', 'scalar', 'integer', '>=', 1});

localCluster = parcluster(config.clusterProfile);
if config.nWorkers > localCluster.NumWorkers
    error(['The requested %d workers exceed the %d workers supported by ' ...
        'profile "%s". Reduce config.nWorkers or adjust the profile.'], ...
        config.nWorkers, localCluster.NumWorkers, config.clusterProfile);
end

existingPool = gcp('nocreate');
if isempty(existingPool)
    parpool(config.clusterProfile, config.nWorkers);
else
    if existingPool.NumWorkers ~= config.nWorkers
        fprintf(['Using the existing parallel pool with %d workers ' ...
            'instead of the requested %d workers.\n'], ...
            existingPool.NumWorkers, config.nWorkers);
    end
    config.nWorkers = existingPool.NumWorkers;
end
end

function x = positive_normal(mu, sigma, n)
    x = mu + sigma * randn(1, n);

    invalid = (x <= 0);
    while any(invalid)
        x(invalid) = mu + sigma * randn(1, nnz(invalid));
        invalid = (x <= 0);
    end
end