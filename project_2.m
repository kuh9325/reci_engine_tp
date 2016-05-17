CA = -360 : 360; %Crank Angle (deg)
Bo = 83e-03; %Bore (m)
St = 91.4e-03; %Stroke (m)
n = 6; %Cylinder number
Vd = (pi*(Bo^2)*St*n)/4; %displacement (cc)
rc = 16.8; %compression ratio
Vc = Vd/(rc-1); %clearance volume
l = 155e-03; %connecting rod length (m)
a = l/3; %crankshaft radius (assumption that a/l = 1/3)
R = 2500; %angular velocity (RPM) (Maximum)
w = R*(60); %angular velocity (deg/s)
T1 = 25+273; %initial temp (K)
P1 = 101.325; %initial pressure (kPa)
s = a*(cosd(CA))+sqrt(l.^2-((a.*sind(CA)).^2)); %piston displacement
V = Vc + ((pi*Bo^2)/4)*(l+a-s); %volume displacement
C = P1*((Vd+Vc)^1.35); %k=1.35 (isentropic constant)
P = C*(V.^-1.35); %Pressure
T = T1*(Vd./V).^0.35; %Temperature
% piston velocity %
pv = a*w*sind(CA)-w*a.^2*sind(CA).*cosd(CA)/sqrt(l.^2-a.^2*(sind(CA)).^2);
pa = - a*w^2*cosd(CA) - (a^2*w^2*cosd(2*CA))/l; %piston acceleration
pm = 0.1; %piston mass (kg)
IF = -pm*pa; %inertia force (N)
Tq = IF.*sind(CA)*a; %torque (Nm)
subplot(3,2,1)
plot(CA, T)
xlabel('Crank Angle (¡Æ)')
ylabel('Temperature (K)')
subplot(3,2,2)
plot(CA, P)
xlabel('Crank Angle (¡Æ)')
ylabel('Pressure (kPa)')
subplot(3,2,3)
plot(CA, pv)
xlabel('Crank Angle (¡Æ)')
ylabel('Piston Velocity (m/s)')
subplot(3,2,4)
plot(CA, pa)
xlabel('Crank Angle (¡Æ)')
ylabel('Piston Acceleration (m/s^2)')
subplot(3,2,5)
plot(CA, IF)
xlabel('Crank Angle (¡Æ)')
ylabel('Inertia Force (N)')
subplot(3,2,6)
plot(CA, Tq)
xlabel('Crank Angle (¡Æ)')
ylabel('Torque (N*m)')