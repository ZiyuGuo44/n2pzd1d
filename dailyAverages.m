function [mX,tspan] = dailyAverages(X,dt)

[msize,nsize] = size(X);
window = 1/dt; %one points per day.

mX = [];
stdX = [];
jspan = [1:window:nsize];

for j = jspan
    J = j:j+(window-1);
    Xj = X(:,J);
    mXj = mean(Xj,2);
    mX  = [mX,mXj]; %OUPUT.
end

%TIME OUTPUT =============================================================
t0 = 1*dt; %t0 = dt
tmax = nsize*dt; %tmax = 360
dtspan = window*dt;
tspan = [t0:dtspan:tmax];

return