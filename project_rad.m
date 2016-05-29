% function eE = project_2(as)

%% as = 0 or -20 %%

%% initial setting

as = 0;
CA = -180 : 180; % Crank Angle (deg)
Bo = 83e-03; % Bore (m)
St = 91.4e-03; % Stroke (m)
n = 6; % Cylinder number
Vd = (pi*(Bo^2)*St*n)/4; % displacement (m^2)
rc = 10; % compression ratio
Vc = Vd/(rc-1); % clearance volume
l = 155e-03; % connecting rod length (m)
a = l/3; % crankshaft radius (assumption that a/l = 1/3)
R = 2500; % angular velocity (RPM) (Maximum)
w = R*(60); % angular velocity (deg/s)
T1 = 298; % initial temp (K)
P1 = 101.325; % initial pressure (kPa)
s = a*(cosd(CA))+sqrt(l.^2-((a.*sind(CA)).^2)); % piston displacement
% volume displacement
V = zeros(1,length(CA));
for i=1:length(CA)
    V(i) = Vc + (Vd/2)*((2*l/s(i))+1-cosd(CA(i))-sqrt(((2*l/s(i))^.2)-(sind(CA(i)).^2)));
end
% k = 1.35; % isentropic coefficient (1.35)
% C = P1*((Vd+Vc)^k); % isentropic constant for volumetric calculation
% Pi = C*(V.^-k); % Pressure (isentropic)

% piston velocity %
pv = a*w*sind(CA)-w*a.^2*sind(CA).*cosd(CA)/sqrt(l.^2-a.^2*(sind(CA)).^2);
pa = - a*w^2*cosd(CA) - (a^2*w^2*cosd(2*CA))/l; % piston acceleration
pm = 0.1; % piston mass (kg)
IF = -pm*pa; % inertia force (N)
Tq = IF.*sind(CA)*a; % torque (Nm)

fm = 2e-2; % mass of fuel (0.02kg)
shru = 1.31; % specific heat ratio (unburned)
shrb = 1.21; % specific heat ratio (burned)
Te = 2394.5; % ignition temperature (K)

% calculate heat income (molecular mass of methane = 16)
% note that exhaust gas temperature is ignition temperature
Pe = P1*(Te/T1); % Pressure of products
eM = (-74850-(8.314*T1)); % enthalpy of react methane (kJ/kmol)
eOr = 2*(-8.314*T1); % enthalpy of react oxygen (kJ/kmol)
eOp = 2*(83174-8682-(8.314*Te)); % enthalpy of product oxygen (kJ/kmol)
eNr = 7.52*(-8.314*T1); % enthalpy of react nitrogen (kJ/kmol)
eNp = 7.52*(79320-8669-(8.314*Te)); % enthalpy of product nitrogen (kJ/kmol)
eH = 2*(-241820+103508-9904-(8.314*Te)); % enthalpy of product vapor (kJ/kmol)
eC = (-393520+125152-9364)-(8.314*Te); % enthalpy of product carbon dioxide (kJ/kmol)

qin = (eM+eOr+eNr)-(eOp+eNp+eH+eC); % [kJ/kmol]
Qin = (fm/16)*qin; % [kJ]

%% calculating data section by section

% mass fraction; efficiency factor(a)=5, form factor(n)=3, interval = 40
xb = zeros(1,length(CA));
dQ = zeros(1,length(CA));
Q = 0;
for i=1:length(CA)
    if CA(i) <= as
        xb(i) = 0;
        dQ(i) = 0;
    else
        xb(i) = 1-exp(-5*((CA(i)-as)/40).^3);
        % diff of heat by radian angle
        dQ(i) = Qin*((15/40)*(1-xb(i))*((CA(i)-as)/40).^2)*(180/pi);
        Q = Q + dQ(i)*(pi/180);
    end
end

% diff of Volume

dV = zeros(1,length(CA));
for i=1:length(CA)
    dV(i) =(Vd/2)*sind(CA(i))*(1+(cosd(CA(i))*...
        ((((2*l)/s(i))^2)-(sind(CA(i)))^2)^-0.5));
end

Ve = (fm*0.5812*Te)/Pe;

dP = zeros(1,length(CA));
P = zeros(1,length(CA));
P(1) = P1;
for i=1:length(CA)
    if CA(i) < as % before combustion
        dP(i) = (-shru*(P(i)/V(i)))*dV(i)+((shru-1)/V(i))*dQ(i);
        P(i) = P(1)*(V(1)/V(i)).^(shru);
    elseif CA(i) >= as && CA(i) <= as+40 % combustion
        dP(i) = (-shrb*(P(i)/V(i)))*dV(i)+((shru-1)/V(i))*dQ(i);
        P(i+1) = P(i) + dP(i)*(pi/180); % pressure that concerned burning
    else % after combustion
        dP(i) = (-shrb*(P(i)/V(i)))*dV(i)+((shru-1)/V(i))*dQ(i);
        P(i) = Pe*(Ve/V(i)).^(shrb);
    end
%     if i ~= length(CA)
%         P(i+1) = P(i) + dP(i)*(pi/180); % pressure that concerned burning
%     end
end

% Temperature than concerns combustion

Tb = zeros(1,length(CA));
Tu = zeros(1,length(CA));

% for i=1:length(CA)
%     if CA(i) < as % before combustion (zero index is 181)
%     T(i) = T1*(V(1)/V(i))^(shru-1); 
%     elseif CA(i) >= as && CA(i) <= as+40 % combustion
%     T(i) = Te*(xb(i)/(V(i)/V(181-as)))^(shrb-1); 
%     else % after combustion
%     T(i) = Te*(V(181-as)/V(i))^(shrb-1);
%     end
% end


for i=1:length(CA)
    if CA(i) < as % before combustion (zero index is 181)
        Tu(i) = T1*(V(1)/V(i)).^(shru-1); % assuming that this process is isentropic
        Tb(i) = 0;
    elseif CA(i) >= as && CA(i) <= as+40 % combustion
        Tu(i) = (P(i)*V(i))/(fm*0.287);
        Tb(i) = (P(i)*V(i))/(fm*0.287);
    else % after combustion
        Tu(i) = 0;
        Tb(i) = Te*(Ve/V(i))^(shrb-1);
    end
end

% engine work

dW = zeros(1,length(CA)-1);
W = 0;
for i=1:length(CA)
    dW(i) = P(i)*dV(i);
    W = W + dW(i)*(pi/180);
end

eE = W/Q; % engine efficiency

%% plotting w/ print the efficiency

% xlabel('Crank Angle (¡Æ)')
% ylabel('Pressure (kPa)')
% plot(CA, P)
subplot(2,2,1)
plot(CA, xb)
subplot(2,2,2)
plot(CA, dQ)
figure(2)
plot(CA, P)
hold on
plot(CA, Tu)
plot(CA, Tb)
fprintf('engine efficiency : %4.3f\n',eE);