function Surface = genSurface(Support)
%clear, clc
% Support =     [8.5203    8.5315    8.5792    8.5766    8.6202    8.1842    8.1741    8.5673    8.5673    8.8287     8.0
%                1.2984    1.3534    1.3032    1.3057    1.2340    1.1687    1.4311    1.4311    1.3406    1.3378     1.3
%                6.8674    6.8675    6.8616    6.8533    6.8514    6.8514    6.8514    6.8514    6.8425    6.8426     6.8];

% Support = [7.8997    7.8284    8.1646    7.6752    8.0187    8.0187    8.0723    8.1487    8.1487    8.1487    8.0977    7.9611    7.9611    7.9109    7.6449    7.6449    7.8159    8.1284    8.1284    7.8761    7.9911    7.9911    7.7290    7.7290    7.7290    7.7290    7.5033    7.6604    7.4490    7.5152    7.5152    7.5152    7.5152    7.6766    7.7952    7.6947    7.8234    7.8234    8.1482    7.5784    7.5784
%            1.1742    1.1742    1.1742    1.1742    1.1742    1.3026    1.3026    1.2180    1.2180    1.0614    1.0807    1.1240    0.9480    0.8779    0.8429    1.0169    1.0169    1.0169    1.0169    1.0169    1.3169    1.1227    1.1227    0.9193    0.9193    0.9193    0.9193    0.9193    0.9193    0.9193    0.7560    0.9281    0.7527    0.7527    0.7527    0.7527    0.9785    0.8328    1.2215    0.9961    0.9961
%            3.9435    3.9401    3.9501    3.9748    3.9575    3.9619    3.9660    3.9597    3.9597    3.9678    3.9294    3.9316    3.9288    3.9159    3.8636    3.8922    3.8983    3.9569    3.9569    3.9020    4.0367    3.9549    3.9505    3.9234   -0.0051    3.8802    3.8573    3.8867    3.8387    3.8560    3.8269    3.8499    3.8398    3.8547    3.9198    3.8575    3.9081    3.8957    3.9550    3.8923    3.8923];

% Support =    [5.6283    5.6768    5.6076    5.6074    5.5566    5.5566    5.5566    5.5566    5.5647    5.5647    5.5647    5.6976     5.6565    5.5310    5.3305    5.2610    5.3775    5.4815    5.5898    5.6431    5.6431    5.2422    5.3875    5.5447    5.5097    5.5632    5.5632    5.4753    5.5768    5.4988    5.6170    5.5760    5.3185    5.3185    5.3965    5.3965    5.3149    5.4529    5.4529    5.4529
%               1.3822    1.3822    1.5068    1.7544    1.9568    1.7031    1.9296    1.7088    1.5042    1.3039    1.2356    1.2356     1.4420    1.6016    1.6016    1.6016    1.4483    1.3522    1.3419    1.1888    1.3625    1.7490    1.7490    1.7490    1.4420    1.6016    1.6016    1.6016    1.4483    1.3522    1.3419    1.1888    1.3625    1.7490    1.7490    1.7490    1.6660    1.7025    1.5958    1.5097
%               19.0022   18.9981   19.0198   19.0390   19.0643   19.0477   19.0708   19.0480   19.0189   18.9878   18.9774   18.9797    19.0175   19.0517   18.9995   18.9704   18.9852   18.9852   19.0126   18.9687   19.0045   18.9962   19.0508   19.0593   19.0742   19.0807   19.0975   19.0892   19.0892   19.0895   19.0709   19.0542   19.0125   18.9937   19.0170   19.0548   19.0194   19.0718   19.0551   19.0321];

interval = 0.02;
Power = 2;
% find 4 corners
[Xmin,Xminpos] = min(Support(1,:));
[Ymin,Yminpos] = min(Support(2,:));
[Xmax,Xmaxpos] = max(Support(1,:));
[Ymax,Ymaxpos] = max(Support(2,:));

Xmat = Xmin:interval:Xmax; % Surface X values
Ymat =  Ymin:interval:Ymax; %Surface Y values
Surface = zeros(3,length(Ymat)*length(Xmat)); % initialize surface
W = zeros(length(Support),1); % Initiallize Weights
% Populate surface X values
for q = 1:length(Ymat)
    for p = 1:length(Xmat)
        Surface(1, p + (q-1)*length(Xmat)) = Xmat(p);
    end
end
% Populate surface Y values
for q = 1:length(Ymat)
    for p = 1:length(Xmat)
        Surface(2, p + (q-1)*length(Xmat)) = Ymat(q);
    end
end
% Populate surface Z values
for m = 1:length(Surface)
    for n = 1:length(Support)
        W(n) = 1/((Support(1,n) - Surface(1,m))^Power + (Support(2,n) - Surface(2,m))^Power); % Weight by distance from point
        Surface(3,m) = Surface(3,m) + Support(3,n) * W(n); %Sum of z value times weights
    end
    Surface(3,m) = Surface(3,m) / sum(W); % Divide by Weights
end

% set up 4 corner points
P(:,1) = Support(:,Xminpos);
P(:,2) = Support(:,Ymaxpos);
P(:,3) = Support(:,Xmaxpos);
P(:,4) = Support(:,Yminpos);
P(:,5) = P(:,1); %Complete the circle

t = 1;
while t == 1 % interate until all suport points are included in shape
    for n = 1:length(P)-1
        for m = 1:length(Support)
            % d = (Ax + By + C) / sqrt(A^2 + B^2)
            % (y1 – y2)x + (x2 – x1)y + (x1y2 – x2y1) = 0
            A(n) = P(2,n) - P(2,n+1);
            B(n) = P(1,n+1) - P(1,n);
            C(n) = P(1,n)*P(2,n+1) - P(1,n+1)*P(2,n);
            D(n,m) = (A(n)*Support(1,m) + B(n)*Support(2,m) + C(n)) / sqrt(A(n)^2 + B(n)^2); % calculate distance from line
        end
    end
    
    if max(max(D)) > 0.01 %if there is a point outside the area
        [Row,Col] = find(D > 0.01); %find the points outside the area
        RowMin = min(Row); %find the first point row
        [Jubjub, ColMin] = max(D(RowMin,:));
        % ColMin = Col(RowIndex); %find the first point Col
        for r = length(P):-1:RowMin+1 %Shift other points to make room 
            P(:,r+1) = P(:,r);
        end
        P(:,RowMin+1) = Support(:,ColMin); %Add new point
    else 
        t = 0; % else end loop
    end

%     clf
%     figure(1)
%     scatter(Support(1,:), Support(2,:))
%     hold on
%     for v = 1:length(P)-1
%     fplot(@(x) -1.*A(v).*x./B(v) - C(v)./B(v))
%     end
%     xlim([7.3 8.3])
%     ylim([0 2])
end 



for n = 1:length(P)-1
    for m = 1:length(Surface)
        % d = (Ax + By + C) / sqrt(A^2 + B^2)
        % (y1 – y2)x + (x2 – x1)y + (x1y2 – x2y1) = 0
        A(n) = P(2,n) - P(2,n+1);
        B(n) = P(1,n+1) - P(1,n);
        C(n) = P(1,n)*P(2,n+1) - P(1,n+1)*P(2,n);
        D(n,m) = (A(n)*Surface(1,m) + B(n)*Surface(2,m) + C(n)) / sqrt(A(n)^2 + B(n)^2);
        if D(n,m) > 0
            Surface(:,m) = 100; % Set surface to 100 if outside of the region
        end
    end
end
Surface(:,any(Surface==100)) = []; % removes columns if there is an element of 100

%get rid of any NAN
[rNAN, cNAN] = find(isnan(Surface));
if length(rNAN) == 0;
    wNAN = 0;
else
    wNAN = 1;
end
while wNAN == 1
    [rNAN, cNAN] = find(isnan(Surface));
    Surface(rNAN(1), cNAN(1)) = Surface(rNAN(1), cNAN(1) - 1);
    [rNAN, cNAN] = find(isnan(Surface));
    if length(rNAN) == 0;
        wNAN = 0;
    end
end

% figure()
% plot(Surface(1,414:431), Surface(3, 414:431))

% adds z stack to surface
zSurface = []; zSurf = []; 
z = [-0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03];
for n = 1:length(Surface)
    for m = 1:length(z)
        zSurf(:, m) = Surface(:, n); % Assigns the Surface position to all locations on zSurf
        zSurf(3, m) = zSurf(3,m) + z(m); % modifies the z position of z
        zSurface(:, (n*5) + m - 5) = zSurf(:, m); %combines zSurf into zsurface hopefully at the correct indices
    end
end
Surface = zSurface;
end