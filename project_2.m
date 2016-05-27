%% initial setting

CA = -180 : 180; % Crank Angle (deg)
Bo = 83e-03; % Bore (m)
St = 91.4e-03; % Stroke (m)
n = 6; % Cylinder number
Vd = (pi*(Bo^2)*St*n)/4; % displacement (cc)
rc = 16.8; % compression ratio
Vc = Vd/(rc-1); % clearance volume
l = 155e-03; % connecting rod length (m)
a = l/3; % crankshaft radius (assumption that a/l = 1/3)
R = 2500; % angular velocity (RPM) (Maximum)
w = R*(60); % angular velocity (deg/s)
T1 = 25+273; % initial temp (K)
P1 = 101.325; % initial pressure (kPa)
s = a*(cosd(CA))+sqrt(l.^2-((a.*sind(CA)).^2)); % piston displacement
V = Vc + ((pi*Bo^2)/4)*(l+a-s); % volume displacement
k = 1.35; % isentropic coefficient
C = P1*((Vd+Vc)^k); % isentropic constant for volumetric calcuration
Pi = C*(V.^-k); % Pressure (isentropic)
% piston velocity %
pv = a*w*sind(CA)-w*a.^2*sind(CA).*cosd(CA)/sqrt(l.^2-a.^2*(sind(CA)).^2);
pa = - a*w^2*cosd(CA) - (a^2*w^2*cosd(2*CA))/l; % piston acceleration
pm = 0.1; % piston mass (kg)
IF = -pm*pa; % inertia force (N)
Tq = IF.*sind(CA)*a; % torque (Nm)
as = 0; % or -20
Qin = 2000;%55625; % [kJ/kg]
shru = 1.31; % specific heat ratio (unburned)
shrb = 1.21; % specific heat ratio (burned)

%% calculating data section by section

% burn fraction; efficiency factor(a)=5, form factor(n)=3, interval = 40
xb = zeros(length(CA),1);
dQ = zeros(length(CA),1);
for i=1:length(CA)
    if CA(i) <= as
    xb(i) = 0;
    dQ(i) = 0;
    else
    xb(i) = 1-exp(-5*((CA(i)-as)/40).^3);
    % diff of heat by radian angle
    dQ(i) = Qin*((15/40/(pi/180))*(1-xb(i))*((CA(i)-as)/40).^2);
    end
end

% dQ = Qin*(15/40)*(1-xb)*(((CA-as)/40).^2);

% diff of Volume
dV = zeros(length(CA),1);
for i=1:length(CA)
dV(i) = (Vd/2)*sind(CA(i))*(1+(cosd(CA(i))*((((2*l)/s(i))^2)-(sind(CA(i)))^2)^-0.5));
end

% diff of Pressure

dP = zeros(length(CA)-1,1);
P = zeros(length(CA)-1,1);
P(1) = P1;
for i=1:length(CA)-1
    dP(i) = (-k*P(i)/V(i))*dV(i)+((k-1)/V(i))*dQ(i);
    P(i+1) = P(i) + dP(i)*pi/180; % pressure that concerned burning
end

for i=1:length(CA)-1
   
end

%% plotting

% subplot(3,2,1)
% plot(CA, T)
% xlabel('Crank Angle (¡Æ)')
% ylabel('Temperature (K)')
% subplot(3,2,2)
plot(CA, P)
xlabel('Crank Angle (¡Æ)')
ylabel('Pressure (kPa)')
% subplot(3,2,3)
% plot(CA, pv)
% xlabel('Crank Angle (¡Æ)')
% ylabel('Piston Velocity (m/s)')
% subplot(3,2,4)
% plot(CA, pa)
% xlabel('Crank Angle (¡Æ)')
% ylabel('Piston Acceleration (m/s^2)')
% subplot(3,2,5)
% plot(CA, IF)
% xlabel('Crank Angle (¡Æ)')
% ylabel('Inertia Force (N)')
% subplot(3,2,6)
% plot(CA, Tq)
% xlabel('Crank Angle (¡Æ)')
% ylabel('Torque (N*m)')