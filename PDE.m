clc
close all
clear all

k.Bdivision = 0.2; % bacterial growth rate per bacterium per macrophage
k.Bdeath = 0.0; % intrinsic bacterial death rate
k.Mdeath = 0.1; % macrophage death rate per bacterium
k.Mformation = 0.1; % rate of macrophage formation
k.Mrecruitment = 0.3; % rate of macrophage recruitment to the infection site
k.Aconc = 0.5; % antibiotic concentration (antibiotic-dependent bacterial death rate)
k.Aperiod = 6; % periodicity of antibiotic regime
k.immuneResponse = 0; % increase in bacterial death rate over time due to the host immunity

k.diffB = 0.4; % bacterial diffusion rate
k.diffM = 0.1; % macrophagal diffusion rate


k.N = 50; % size of the grid

tStart = 0;
tEnd = 10;

% initial bacteria concentrations
minB = 0;
maxB = 1;

% initial macrophage concentrations
minM = 0;
maxM = 10;

% starting conditions
initB = minB + (maxB-minB)*rand(1,k.N^2);
initB = rescale(exp(1).^(30*initB)); % exponential distribution so that very few grid cells have many bacteria, the number you multiply initB by controls the skewness
initM = minM + (maxM-minM)*rand(1,k.N^2);


optionsODE = odeset('AbsTol', 1e-5.*ones(1,k.N^2*2), 'RelTol', 1e-5);
[T, Out] = ode45(@(t,y)tissuePeriodic(t,y,k),[tStart,tEnd],[initB,initM], optionsODE);

%%
% plotting
filename = strcat('animation-',datestr(now,'yymmDD-HHMMSS'), ".gif");
t = 1;
dataB = reshape(Out(t,1:k.N^2),k.N,k.N);
h = heatmap(dataB);
len = (length(T));
for t = 1:len
    if mod(t,10) == 0
        strcat(num2str(100*t/len), '%')  % displays progress
    end
    dataB = reshape(Out(t,1:k.N^2),k.N,k.N);
    h.ColorData = dataB;
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if t == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,...
        'DelayTime',0.02);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append',...
        'DelayTime',0.02);
    end
end
%%
% plotting
filename2 = strcat('animation2-',datestr(now,'yymmDD-HHMMSS'), ".gif");
t = 1;
dataB = reshape(Out(t,1*k.N^2+1:2*k.N^2),k.N,k.N);
hh = heatmap(dataB);
len = (length(T));
for t = 1:len
    if mod(t,10) == 0
        strcat(num2str(100*t/len), '%')  % displays progress
    end
    dataB = reshape(Out(t,1*k.N^2+1:2*k.N^2),k.N,k.N);
    hh.ColorData = dataB;
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if t == 1
        imwrite(imind,cm,filename2,'gif', 'Loopcount',inf,...
        'DelayTime',0.02);
    else
        imwrite(imind,cm,filename2,'gif','WriteMode','append',...
        'DelayTime',0.02);
    end
end

%%
% functions


function  dY = tissuePeriodic(t,Y,k)
dY = zeros(2*k.N^2,1);
for i=1:k.N
    for j = 1:k.N
        if i==1 
            if j==1 %corner
                dY((i-1)*k.N+j) = Bdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k,t) + k.diffB*(-4*Y((i-1)*k.N+j)+Y(i*k.N+j)+Y((i-1)*k.N+j+1)+Y(k.N)+Y(k.N^2-k.N+1));
            elseif j==k.N %corner
                dY((i-1)*k.N+j) = Bdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k,t) + k.diffB*(-4*Y((i-1)*k.N+j)+Y(i*k.N+j)+Y((i-1)*k.N+j-1));
            else %edge
                dY((i-1)*k.N+j) = Bdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k,t) + k.diffB*(-4*Y((i-1)*k.N+j)+Y(i*k.N+j)+Y((i-1)*k.N+j+1)+Y((i-1)*k.N+j-1)+Y(k.N*(k.N-1)+j));
            end
        elseif i == k.N 
            if j==1 %corner
                dY((i-1)*k.N+j) = Bdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k,t) + k.diffB*(-4*Y((i-1)*k.N+j)+Y((i-2)*k.N+j)+Y((i-1)*k.N+j+1)+Y(1)+Y(k.N^2));
            elseif j==k.N %corner
                dY((i-1)*k.N+j) = Bdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k,t) + k.diffB*(-4*Y((i-1)*k.N+j)+Y((i-2)*k.N+j)+Y((i-1)*k.N+j-1)+Y(k.N)+Y(k.N^2+1-k.N));
            else %edge
                dY((i-1)*k.N+j) = Bdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k,t) + k.diffB*(-4*Y((i-1)*k.N+j)+Y((i-1)*k.N+j+1)+Y((i-2)*k.N+j)+Y((i-1)*k.N+j-1)+Y(j));
            end
        else
            if j==1 %edge
                dY((i-1)*k.N+j) = Bdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k,t) + k.diffB*(-4*Y((i-1)*k.N+j)+Y(i*k.N+j)+Y((i-1)*k.N+j+1)+Y((i-2)*k.N+j)+Y(i*k.N));
            elseif j==k.N %edge
                dY((i-1)*k.N+j) = Bdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k,t) + k.diffB*(-4*Y((i-1)*k.N+j)+Y(i*k.N+j)+Y((i-2)*k.N+j)+Y((i-1)*k.N+j-1)+Y(k.N*(i-1)+1));
            else % inside
                dY((i-1)*k.N+j) = Bdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k,t) + k.diffB*(-4*Y((i-1)*k.N+j)+Y(i*k.N+j)+Y((i-1)*k.N+j+1)+Y((i-2)*k.N+j)+Y((i-1)*k.N+j-1));
            end
        end
    end
end
for i=1:k.N
    for j = 1:k.N
        if i==1 
            if j==1 %corner
                dY(k.N^2+(i-1)*k.N+j) = Mdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k) + k.diffM*(-2*Y(k.N^2+(i-1)*k.N+j)+Y(k.N^2+i*k.N+j)+Y(k.N^2+(i-1)*k.N+j+1));
            elseif j==k.N %corner
                dY(k.N^2+(i-1)*k.N+j) = Mdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k) + k.diffM*(-2*Y(k.N^2+(i-1)*k.N+j)+Y(k.N^2+i*k.N+j)+Y(k.N^2+(i-1)*k.N+j-1));
            else %edge
                dY(k.N^2+(i-1)*k.N+j) = Mdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k) + k.diffM*(-3*Y(k.N^2+(i-1)*k.N+j)+Y(k.N^2+i*k.N+j)+Y(k.N^2+(i-1)*k.N+j+1)+Y(k.N^2+(i-1)*k.N+j-1));
            end
        elseif i == k.N 
            if j==1 %corner
                dY(k.N^2+(i-1)*k.N+j) = Mdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k) + k.diffM*(-2*Y(k.N^2+(i-1)*k.N+j)+Y(k.N^2+(i-2)*k.N+j)+Y(k.N^2+(i-1)*k.N+j+1));
            elseif j==k.N %corner
                dY(k.N^2+(i-1)*k.N+j) = Mdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k) + k.diffM*(-2*Y(k.N^2+(i-1)*k.N+j)+Y(k.N^2+(i-2)*k.N+j)+Y(k.N^2+(i-1)*k.N+j-1));
            else %edge
                dY(k.N^2+(i-1)*k.N+j) = Mdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k) + k.diffM*(-3*Y(k.N^2+(i-1)*k.N+j)+Y(k.N^2+(i-1)*k.N+j+1)+Y(k.N^2+(i-2)*k.N+j)+Y(k.N^2+(i-1)*k.N+j-1));
            end
        else
            if j==1 %edge
                dY(k.N^2+(i-1)*k.N+j) = Mdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k) + k.diffM*(-3*Y(k.N^2+(i-1)*k.N+j)+Y(k.N^2+i*k.N+j)+Y(k.N^2+(i-1)*k.N+j+1)+Y(k.N^2+(i-2)*k.N+j));
            elseif j==k.N %edge
                dY(k.N^2+(i-1)*k.N+j) = Mdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k) + k.diffM*(-3*Y(k.N^2+(i-1)*k.N+j)+Y(k.N^2+i*k.N+j)+Y(k.N^2+(i-2)*k.N+j)+Y(k.N^2+(i-1)*k.N+j-1));
            else % inside
                dY(k.N^2+(i-1)*k.N+j) = Mdynamics(Y((i-1)*k.N+j),Y(k.N^2+(i-1)*k.N+j),k) + k.diffM*(-4*Y(k.N^2+(i-1)*k.N+j)+Y(k.N^2+i*k.N+j)+Y(k.N^2+(i-1)*k.N+j+1)+Y(k.N^2+(i-2)*k.N+j)+Y(k.N^2+(i-1)*k.N+j-1));
            end
        end
    end
end


end



%%

function dB = Bdynamics(B,M,k,t)
if (B <= 0)
dB = k.Bdivision*B*M;
    if (M<=0)
        dB = 0;
    end  
else
% dB = k.Bdivision*B*M - sqrt(B)*k.Aconc*(-cos(t*k.Aperiod*pi)+1) - k.Bdeath*B*(1+t*k.immuneResponse);
dB = k.Bdivision*B*M - sqrt(B)*k.Aconc*(-cos(t*k.Aperiod*pi)+1);
end
end
function dM = Mdynamics(B,M,k)

if (M <= 0)
dM = k.Mformation + k.Mrecruitment*B;
else  
dM = -k.Mdeath*B/M + k.Mformation + k.Mrecruitment*B;
end
end
