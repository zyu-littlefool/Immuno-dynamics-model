function [state_cell_i, E_cell_wt, E_cell_var, gen_cell_i, gen_num_i, del_num_i, div_num_1] = ...
    Copy_of_New_xiaoshi(state_cell_i, gen_num, gen_cell_i, E_cell_i, E_cell_j,  del_num, div_num_1, stream, stateCell)
    d_PC = 0.17;
    dt = 0.01;
    stream.State = stateCell{9};
    Pr2 = 1 - exp(-d_PC*dt);
    rr = rand(stream, 1,length(state_cell_i));
    match = (state_cell_i == 5) & (rr < Pr2);
    E_cell_i(match) = [];
    E_cell_j(match) = [];
    gen_cell_i((match),:) = [];
    gen_num(match) = [];
    del_num(match) = [];
    div_num_1(match) = [];
    state_cell_i(match) = [];
    gen_num_i = gen_num;
    del_num_i = del_num;
    E_cell_wt = E_cell_i;
    E_cell_var = E_cell_j;
end

