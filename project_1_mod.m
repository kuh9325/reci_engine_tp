CA = -360 : 360; %Crank Angle (deg)
syms t
Bo = 83; %Bore (mm)
St = 91.4; %Stroke (mm)
n = 6; %Cylinder number
Vd = (pi*(Bo^2)*St*6)/4; %displacement (cc)
r = 16.8; %compression ratio
Vc = Vd/(r-1); %clearance volume
l = 155e-03; %connecting rod length (mm)
a = l/3; %crankshaft radius (assumption that a/l = 1/3)
T1 = 25+273; %initial temp (K)
P1 = 101.325; %initial pressure (kPa)
s = a*(cos(t*(pi/180)))+(l*sqrt(1-((a/l).^2)*((cos(t*(pi/180))).^2))); %piston displacement
V = Vc + ((pi*Bo^2)/4)*(l+a-s); %volume displacement
t3 = solve('a*(cos(t*(pi/180)))+(l*sqrt(1-((a/l)^2)*((cos(t*(pi/180)))^2))) = (-Vc*(4/((Bo^2)*pi))+(l+a))',t);
% crank angle @ point 3
if CA >= -180 && t <= 0 %process 1-2
    C = P1*(Vd^1.35); %k=1.35
    T = ((V^-0.35)*C)/0.287; %R=0.287
    P = C*(V^-1.35);
elseif CA >= 0 && t <= t3 %process 2-3
    T = T1*(r^0.4)*2; %cut-off ratio = 2
    P = P1*(r^1.4);
elseif t >= t3 && t <= 180 %process 3-4
    T = (T1*(r^0.4)*2)*(2/r)^0.4;
    P = (P1*(r^1.4))*((T1*(r^0.4)*2)*(2/r)^0.4)/(T1*(r^0.4)*2)*(2/r);
else %process 4-1
    T = T1;
    P = P1;
end
ezplot(T,[-360 360])
pv = diff(s,t); %piston velocity
pa = diff(pv,t); %piston acceleration
pm = 0.1; %piston mass
IF = -pm*pa; %inertia force