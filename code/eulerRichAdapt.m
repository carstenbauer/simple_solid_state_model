% Euler-Richardson Adapative Solver for integration of equations of motion

function [t, y] = eulerRichAdapt(y0,dy,t0,dt,t1,RelTol)

t = t0; %set start time
tc = t0; %set current time to start time
y = y0; %set start value


while (tc < t1)

    while 1

        [t_dump1 y_big] = eulerRich(dy,tc:dt:(tc+dt),y(end,:)); %one step Euler-Richardson
        [t_dump2 y_small] = eulerRich(dy,tc:(dt./2):(tc+dt),y(end,:)); %two half-steps Euler-Richardson

        e = abs(y_big(2,:) - y_small(3,:)); %difference between one/two steps

        if (norm(e) > RelTol * norm(y_small(3,:))) %if not accurate enough
            dt = dt ./2; %half step, try again
        elseif (norm(e) < RelTol ./4 * norm(y_small(3,:))) %if too accurate
            t  = [t; (tc+dt)]; %update time
            y = [y; y_small(3,:)]; %update value
            tc = tc+dt; %progress current time
            dt = 2 .* dt; %double step, continue
            break
        else %if everything is alright
            t = [t; (tc+dt)]; %update time
            tc = tc+dt; %progress current time
            y = [y; y_small(3,:)]; %update value, continue
            break
        end
    end   
end

if (t(end)~=t1) %if upper time limit is not reached exactly
    if (t(end-1)<t1)
        y(end,:) = y(end-1,:) + (y(end,:)-y(end-1,:))/(t(end)-t1)*(t1-t(end-1));
    end
end