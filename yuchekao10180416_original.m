% Parameters
U=1500;
Rg=1.987;
Goi=3.98e7;
Goii=4.81e11;
Kgi=2.55e5;
Kgii=5.51e5;
Tm=212+273.15;
No=63;
k=1/54978;

t = 0:0.02:5;
ntsteps = numel(t);
T = 170 - 22 * t;

% Pre-allocate arrays rather than resizing in the loops
[G, r, R, v, N] = deal(zeros(1, ntsteps));

% Growth rate
for i = 1:ntsteps
    if T(i)>=120
G(i)=(Goi*exp(-U/(Rg*(T(i)-30)))*exp(-Kgi*(T(i)+273.15+Tm)/(2*(T(i)+273.15)*(T(i)+273.15)*(Tm-T(i)-273.15))))/60; 
else 
G(i)=(Goii*exp(-U/(Rg*(T(i)-30)))*exp(-Kgii*(T(i)+273.15+Tm)/(2*(T(i)+273.15)*(T(i)+273.15)*(Tm-T(i)-273.15))))/60; 
   end
r(i)=G(i).*0.02;
R(i)=sum(r(1:i));
v(i)=((4/3)*pi)*R(i).^3;
end

% PLA nucleation 
 for i=1:ntsteps
     N(i)=No*exp(-(1/k)/((T(i)+273.15)*(Tm-T(i)-273.15)));
 end

% Rt(i).cumsum = cumsum(r(i:end)) = Rt(i)=struct('cumsum',cumsum(r(i:end)))
for i = 1:ntsteps
    Rt1(i).cumsum = cumsum(r(i:end));
    V1(i).cumsum=((4/3)*pi*Rt1(i).cumsum.^3);
    Vt1(i).cumsum=((4/3)*pi*round(N(i))*Rt1(i).cumsum.^3);
end

% Calculate total volume, Vtot.
Vtot1 = 0;
for i = 1:ntsteps 
    V1 = Vt1(i).cumsum(end);
    Vtot1 = Vtot1 + V1;
end

% calculate crystals volumes in each column (from Vtot2 to Vtot251) do not
% need to plot them
for i = 1:ntsteps
    Rt2(i).cumsum = cumsum(r(i:end-1));
    Vt2(i).cumsum=((4/3)*pi*round(N(i))*Rt2(i).cumsum.^3);
    Rt3(i).cumsum = cumsum(r(i:end-2));
    Vt3(i).cumsum=((4/3)*pi*round(N(i))*Rt3(i).cumsum.^3);
    Rt4(i).cumsum = cumsum(r(i:end-3));
    Vt4(i).cumsum=((4/3)*pi*round(N(i))*Rt4(i).cumsum.^3);
    Rt5(i).cumsum = cumsum(r(i:end-4));
    Vt5(i).cumsum=((4/3)*pi*round(N(i))*Rt5(i).cumsum.^3);
    Rt6(i).cumsum = cumsum(r(i:end-5));
    Vt6(i).cumsum=((4/3)*pi*round(N(i))*Rt6(i).cumsum.^3);
    Rt7(i).cumsum = cumsum(r(i:end-6));
    Vt7(i).cumsum=((4/3)*pi*round(N(i))*Rt7(i).cumsum.^3);
    Rt8(i).cumsum = cumsum(r(i:end-7));
    Vt8(i).cumsum=((4/3)*pi*round(N(i))*Rt8(i).cumsum.^3);
    Rt9(i).cumsum = cumsum(r(i:end-8));
    Vt9(i).cumsum=((4/3)*pi*round(N(i))*Rt9(i).cumsum.^3);
    Rt10(i).cumsum = cumsum(r(i:end-9));
    Vt10(i).cumsum=((4/3)*pi*round(N(i))*Rt10(i).cumsum.^3);
    Rt11(i).cumsum = cumsum(r(i:end-10));
    Vt11(i).cumsum=((4/3)*pi*round(N(i))*Rt11(i).cumsum.^3);
    Rt12(i).cumsum = cumsum(r(i:end-11));
    Vt12(i).cumsum=((4/3)*pi*round(N(i))*Rt12(i).cumsum.^3);
    .....................
    Rt251(i).cumsum = cumsum(r(i:end-250));
    Vt251(i).cumsum=((4/3)*pi*round(N(i))*Rt251(i).cumsum.^3);
end

% calculate total volumes in each column (from Vtot2 to Vtot251)
Vtot2 = 0;
for i = 1:ntsteps 
    V2 = Vt2(i).cumsum(end);
    Vtot2 = Vtot2 + V2;
end

Vtot3 = 0;
for i = 1:ntsteps 
    V3 = Vt3(i).cumsum(end);
    Vtot3 = Vtot + V3;
end

Vtot4 = 0;
for i = 1:ntsteps 
    V4 = Vt4(i).cumsum(end);
    Vtot4 = Vtot + V4;
end
% and go on



% Plotting
figure(1);
subplot(1,3,1);plot(t,G,'LineWidth',3);
xlabel('Time (second)');
ylabel('Growth rate (µm s^-1)');
subplot(1,3,2);plot(T,G,'LineWidth',3);
xlabel('Temperature (°C)');
ylabel('Growth rate (µm s^-1)');
axis([60,170,0,0.08]);set(gca,'xtick',60:10:170);set(gca,'xdir','reverse')
subplot(1,3,3);plot(t,r,'LineWidth',3);
xlabel('Time (second)');
ylabel('Radius (µm)');

figure(2);
subplot(1,2,1);plot(t,R,'LineWidth',3);
xlabel('Time (second)');
ylabel('Radius (µm)');
set(gca,'ytick',0:0.01:0.18);
subplot(1,2,2);plot(t,v,'LineWidth',3);
xlabel('Time (second)');
ylabel('Spherical volume (µm^3)');
set(gca,'ytick',0:0.001:0.025);

figure(3);
hold on
for i = 1:numel(Rt1)
    plot(t(i:end), Rt1(i).cumsum,'LineWidth',1.25);
end
xlabel('Time (second)');
ylabel('Radius (µm)');
set(gca,'ytick',0:0.01:0.18);

figure(4);
hold on
for i=1:numel(V1)
    plot(t(i:end), V1(i).cumsum,'LineWidth',1.25);
end
xlabel('Time (second)');
ylabel('Spherical volume (µm^3)');
set(gca,'ytick',0:0.002:0.024);

figure(5);
subplot(1,2,1);plot(t,N,'LineWidth',3);
xlabel('Time (second)');
ylabel('Number of spherulities');
axis([0,5,0,25]);set(gca,'xtick',0:0.5:5);set(gca,'ytick',0:2:24);
subplot(1,2,2);plot(T,N,'LineWidth',3);
xlabel('Temperature (°C)');
ylabel('Number of spherulities');
axis([60,170,0,25]);set(gca,'xtick',60:10:170);set(gca,'xdir','reverse');set(gca,'ytick',0:2:24)

figure(6);
hold on
for i=1:numel(Vt1)
    plot(t(i:end),Vt1(i).cumsum,'LineWidth',1.25);
end
xlabel('Time (second)');
ylabel('Spherical volume (µm^3)');
set(gca,'ytick',0:0.01:0.15);

figure(7);
plot(t(1:end),Vt1(1).cumsum,t(26:end),Vt1(26).cumsum,t(51:end),Vt1(51).cumsum,t(101:end),Vt1(101).cumsum,t(151:end),Vt1(151).cumsum,t(201:end),Vt1(201).cumsum,'LineWidth',3);
xlabel('Time (second)');
ylabel('Spherical volume (µm^3)');
legend('from 0s --> 5s','from 0.5s --> 5s','from 1s --> 5s','from 2s --> 5s','from 3s --> 5s','from 4s --> 5s')




