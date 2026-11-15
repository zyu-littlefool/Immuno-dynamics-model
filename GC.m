function [state_cell_i, E_cell_i, E_cell_j, gen_cell_i, gen_num, del_num, SHM,SHM_pos, SHM_pos_j, div_num] = ...
    Copy_4_of_GC(state_cell_i, E_cell_i, E_cell_j, gen_cell_i, gen_num, del_num, C, N_GC, T, T_FH, xE, xE_2, SHM, SHM_pos, SHM_pos_j, div_num, stream_1, stateCell,Ag_over)

    d_GC = 0.5;
    d_Bact = 0.5; %match 3 day to 5%
    % d_PC = 0.17;
    dt = 0.01;
    E0 = 6;
    C0 = 0.2;
    beta_max = 4;
    p1 = 0.1;
    p2 = 0.1;
    % p3 = 0.3;
    % K = 0.5;

    
    ActB_decay_start = 100;      % delay before decay starts, counted by del_num

    %% activation
    stream_1.State = stateCell{1};
    E_cell_active = max(E_cell_i, E_cell_j);
    A = sqrt(10.^(min(E_cell_active,10) - E0)*(C/C0));
    Pr_active = min(A,1);
    rr = rand(stream_1, 1, length(A));
    match = (state_cell_i == 0) & (rr < Pr_active);
    state_cell_i(match)=1; % get antigen
    del_num(match) = 0;


    %% enter GC
    N_activated = length(find(state_cell_i == 1));
    A_mean = mean(A(state_cell_i == 1));
    lamda = ((T/(N_activated*A_mean)).*A)./(1 + (T/(N_activated*A_mean)).*A);
    Pr_enter = 1-exp(-dt*lamda);
    stream_1.State = stateCell{2};
    rr = rand(stream_1, 1, length(A));
    match = (state_cell_i == 1) & (rr < Pr_enter);
    state_cell_i(match) = 2;  % in gc
    % cells that successfully enter GC no longer need the ActB waiting timer
    del_num(match) = 0;
    gen_num(match) = 0;


    actB_id = find(state_cell_i == 1);

    if ~isempty(actB_id)

    del_num(actB_id) = del_num(actB_id) + 1;

    actB_age = del_num(actB_id);
    decay_id = actB_id(actB_age >= ActB_decay_start);
        if ~isempty(decay_id)
            stream_1.State = stateCell{9};
            rr_decay = rand(stream_1, 1, length(decay_id));
            Pr_ActB_decay = 1 - exp(-d_Bact*dt);
            decay_match = decay_id(rr_decay < Pr_ActB_decay);
            if Ag_over == 1
                state_cell_i(decay_match) = [];
                del_num(decay_match) = [];
                A(decay_match) = [];
                E_cell_i(decay_match) = [];
                E_cell_j(decay_match) = [];
                gen_cell_i((decay_match),:) = [];
                gen_num(decay_match) = [];
                div_num(decay_match) = [];
            elseif Ag_over == 0        
                state_cell_i(decay_match) = 0;
                del_num(decay_match) = 0;
                gen_num(decay_match) = 1;
            end
        end

    end

    %% positive select
    N_activated = length(find(state_cell_i == 2));
    A_mean = mean(A(state_cell_i == 2));
    beta = beta_max*(T_FH.*A/(N_GC*N_activated*A_mean))./(1 + (T_FH.*A/(N_GC*N_activated*A_mean)));
    Pr_effec = 1-exp(-dt*beta);
    stream_1.State = stateCell{3};
    rr = rand(stream_1, 1, length(A));
    match = (state_cell_i == 2) & (rr < Pr_effec);
    state_cell_i(match) = 3; % be selected

    %% proliferation & differariation

    stream_1.State = stateCell{4};
    rr = rand(stream_1,1, length(A));
    match = (state_cell_i  == 3) & (rr < p1);
    state_cell_i(match) = 12; % ready to leave GC
    gen_num(match) = 0;
    gen_num(state_cell_i == 12) = gen_num(state_cell_i == 12) + 1;
    match = (state_cell_i == 12) & (gen_num > 25); % divition time to PC or BM; dt = 0.01
    state_cell_i(match) = 4; %leave GC
    match_2 = (state_cell_i == 3) & (rr >= p1); %proliferation & mutation
    state_cell_i(match_2) = 11; % ready to mutate
    gen_num(match_2) = 0;
    gen_num(state_cell_i == 11) = gen_num(state_cell_i == 11) + 1;
    match = (state_cell_i == 11); %
    state_cell_i(match) = 7;
    gen_num(state_cell_i == 7) = 0;

    pmut = 0.6-div_num*0.4/5;
    stream_1.State = stateCell{5};
    rr = rand(stream_1,1, length(state_cell_i));
    state_cell_i((state_cell_i == 7) & (rr < pmut)) = 9; % mutate
    state_cell_i((state_cell_i == 7) & (rr >= pmut)) = 8; % double
    div_num((state_cell_i == 9)) = div_num((state_cell_i == 9)) + 1;

    %%  proliferation
    state_cell_new = [8*ones(1,length(find(state_cell_i == 8))), 9*ones(1,length(find(state_cell_i == 9)))];
    E_cell_new_1 = [E_cell_i(state_cell_i == 8), E_cell_i(state_cell_i == 9)];
    E_cell_new_2 = [E_cell_j(state_cell_i == 8), E_cell_j(state_cell_i == 9)];
    gen_cell_new = [gen_cell_i(state_cell_i == 8, :); gen_cell_i(state_cell_i == 9, :)];
    gen_num_new = ones(1, length(state_cell_new));
    del_num_new = ones(1, length(state_cell_new));
    div_num_new = [div_num(state_cell_i == 8), div_num(state_cell_i == 9)];
    state_cell_i = [state_cell_i, state_cell_new];
    E_cell_i = [E_cell_i, E_cell_new_1];
    E_cell_j = [E_cell_j, E_cell_new_2];
    gen_cell_i = [gen_cell_i; gen_cell_new];
    del_num = [del_num, del_num_new];
    gen_num = [gen_num, gen_num_new];
    div_num = [div_num, div_num_new];

    
    %% differariation

    stream_1.State = stateCell{6};
    rr = rand(stream_1,1, length(state_cell_i));
    div_num(state_cell_i == 4) = 1;
    match = (state_cell_i == 4) & (rr < p2);
    state_cell_i(match) = 5; %PC
    state_cell_i(state_cell_i == 4) = 6; % BM


    %% mutation

    mut_id = find(state_cell_i == 9);
    if ~isempty(mut_id)
        gen_cell_new = gen_cell_i(mut_id,:);
        stream_1.State = stateCell{7};
        rr = 1000*rand(stream_1,1,length(mut_id));
        wz = 1 + floor(rr);
        row = 1:numel(wz);
        idx_change = sub2ind(size(gen_cell_new), row, wz);
        if gen_cell_new(idx_change)  == 0
            gen_cell_new(idx_change) = 1;
        else
            gen_cell_new(idx_change) = 0;
        end
        gen_cell_i(mut_id,:) = gen_cell_new;
        E_cell_i(mut_id) = E_cell_i(mut_id) + (gen_cell_i(mut_id,:) * xE)';
        E_cell_j(mut_id) = E_cell_j(mut_id) + (gen_cell_i(mut_id,:) * xE_2)';
        changeE = (gen_cell_i(mut_id, :) * xE)';
        changeE_j = (gen_cell_i(mut_id,:) * xE_2)';
        mut_death = mut_id(changeE < -50);
        E_cell_i(mut_death) = [];
        E_cell_j(mut_death) = [];
        gen_cell_i((mut_death),:) = [];
        gen_num(mut_death) = [];
        del_num(mut_death) = [];
        state_cell_i(mut_death) = [];
        div_num(mut_death) = [];
        changeE(changeE < -50) = [];
        changeE_j(changeE < -50) = [];
        new_SHM_pos = length(find(changeE > 0));
        new_SHM_pos_j = length(find(changeE_j > 0));
        new_SHM = length(changeE);
        SHM = SHM + new_SHM;
        SHM_pos = SHM_pos + new_SHM_pos;
        SHM_pos_j = SHM_pos_j + new_SHM_pos_j;
    end

    % back
    state_cell_i(div_num >= 7) = 4;
    state_cell_i(state_cell_i == 8 | state_cell_i == 9) = 10; % ready to back to LZ
    gen_num(state_cell_i == 10) = gen_num(state_cell_i == 10) + 1;
    match_exit = (state_cell_i == 10) & (gen_num > 25); 
    match_exit_id = (match_exit == 1);
    state_cell_i(match_exit_id) = 0;
    gen_num(match_exit_id) = 1;
    
     %GC decay

        Pr1 = 1 - exp(-d_GC*dt);
        stream_1.State = stateCell{8};
        rr = rand(stream_1, 1, length(state_cell_i));
        match = (state_cell_i == 2) & (rr < Pr1);
        E_cell_i(match) = [];
        E_cell_j(match) = [];
        gen_cell_i((match),:) = [];
        gen_num(match) = [];
        del_num(match) = [];
        div_num(match) = [];
        state_cell_i(match) = [];
    
end

