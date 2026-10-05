

hbar = 1;
kB = 1;



energy = 0.2;
Delta = 1;

H = [-energy/2, -Delta/2;
     -Delta/2,   energy/2];

HS = H;

sigmaz = [1, 0;
          0,-1];



rho_nought = [1, 0;
              0, 0];



Nb = 15;

a = diag(sqrt(1:Nb-1),1);

adag = a';

omega = 1.2;

HB = hbar*omega*(adag*a);



Isystem = eye(2);
Ibath = eye(Nb);



HS_full = kron(HS,Ibath);

HB_full = kron(Isystem,HB);



g = 0.2;

HI = g*kron(sigmaz,(a + adag));



Htotal = HS_full + HB_full + HI;



T = 1;

rhoB = expm(-HB/(kB*T));
rhoB = rhoB/trace(rhoB);

rhoTotal0 = kron(rho_nought,rhoB);

rhoVec0 = rhoTotal0(:);

dim = size(Htotal,1);

tspan = 0:0.01:100;

[t,rhoVec] = ode45(@(t,rhoVec) ...
    densityEquation(t,rhoVec,Htotal,hbar,dim), ...
    tspan,rhoVec0);

n = length(t);



rho1 = zeros(n,1);
rho2 = zeros(n,1);
coherence = zeros(n,1);
purity = zeros(n,1);



for k = 1:n

    rhoTotal = reshape(rhoVec(k,:),dim,dim);

    rhoSystem = zeros(2,2);

    for i = 1:2
        for j = 1:2

            rows = (i-1)*Nb + (1:Nb);
            cols = (j-1)*Nb + (1:Nb);

            rhoSystem(i,j) = trace(rhoTotal(rows,cols));

        end
    end

    rho1(k) = real(rhoSystem(1,1));
    rho2(k) = real(rhoSystem(2,2));

    coherence(k) = abs(rhoSystem(1,2));

    purity(k) = real(trace(rhoSystem*rhoSystem));

end


figure

plot(t,rho1,'LineWidth',1.5)
hold on
plot(t,rho2,'LineWidth',1.5)

xlabel('Time')
ylabel('Population')
legend('State 1','State 2')
title('TLS Populations')
grid on



figure

plot(t,coherence,'LineWidth',1.5)

xlabel('Time')
ylabel('|rho_{12}|')
title('TLS Coherence')
grid on


figure

plot(t,purity,'LineWidth',1.5)

xlabel('Time')
ylabel('Tr(\rho^2)')
title('TLS Purity')
grid on



function drhoVec = densityEquation(~,rhoVec,Htotal,hbar,dim)

    rhoTotal = reshape(rhoVec,dim,dim);

    drhoTotal = -(1i/hbar) * ...
        (Htotal*rhoTotal - rhoTotal*Htotal);

    drhoVec = drhoTotal(:);

end