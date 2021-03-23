% decodes a vector u -> matrix

function u = state_dec(u, N_y)
u = (vec2mat(u,N_y)).';