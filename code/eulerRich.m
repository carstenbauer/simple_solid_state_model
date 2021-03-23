% Euler-Richardson Solver for integration of equations of motion

function [t, y_total] = eulerRich(dy,times,y0)

y_total(1,:) = y0;
t = times';

for tcur = times(1:end-1);
    idx = find(times==tcur);
    tau = times(idx+1) - times(idx);
    ycur = y_total(end,:);

    k1 = tau*dy(tcur, ycur);   
    k2 = tau*dy((tcur+tau/2) , (ycur+k1.'/2));
    y = ycur + k2.';
    y_total = [y_total; y];
end