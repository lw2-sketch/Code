close all

kb = 1;
hbar = 1;

T = 20;
gamma = 0.3;

N = 30;

a = diag(sqrt(1:N-1),1);

adag = a';


m = 1;
w = 1;

x = sqrt(hbar/(2*m*w)) * (a + adag);

p = -1i * sqrt((m * hbar * w)/2) * (a - adag);

L = (sqrt((4*m*gamma*kb*T)/((hbar)^2)) * x) + 1i*(sqrt(gamma/(4*m*kb*T)) * p);

H = hbar*w*(adag*a + 0.5*eye(N));
Heff = H + ((gamma/2)*(x*p + p*x));

psi0 = zeros(N,1);

psi0(1) = 1/sqrt(2);
psi0(2) = 1/sqrt(2);

rhonought = psi0*psi0';

rhoVecnought = rhonought(:);

timespan = 0:0.01:20;


[t,rhoVec] = ode45(@(t,rhoVec) ...
    CL_Lindblad(timespan,rhoVec,Heff,L,hbar,N), ...
    timespan,rhoVecnought);


nt = length(t);

P0 = zeros(nt,1);
P1 = zeros(nt,1);
Prest = zeros(nt,1);
Psubspace = zeros(nt,1);


coherence01 = zeros(nt,1);


purity = zeros(nt,1);


traceRho = zeros(nt,1);


for k = 1:nt

    
    rho = reshape(rhoVec(k,:),N,N);

    
    P0(k) = real(rho(1,1));

    P1(k) = real(rho(2,2));

    Prest(k) = sum(real(diag(rho(3:end,3:end))));

   
    Psubspace(k) = P0(k) + P1(k);    

    coherence01(k) = abs(rho(1,2));

    purity(k) = real(trace(rho*rho));

    traceRho(k) = real(trace(rho));

end


figure

plot(t,P0,'LineWidth',1.5)
hold on

plot(t,P1,'LineWidth',1.5)

plot(t,Prest,'LineWidth',1.5)

xlabel('Time')
ylabel('Population')

legend('|0>','|1>','Other states','Location','best')

title('Population Dynamics')

grid on

ylim([0 1])



figure

plot(t,coherence01,'LineWidth',1.5)

xlabel('Time')
ylabel('|\rho_{01}|')

title('Coherence Between |0> and |1>')

grid on

ylim([0 0.5])



figure

plot(t,purity,'LineWidth',1.5)

xlabel('Time')
ylabel('Tr(\rho^2)')

title('Purity')

grid on

ylim([0 1])

figure

plot(t,Psubspace,'LineWidth',1.5)

xlabel('Time')
ylabel('P_0 + P_1')

title('Population Remaining in Initial Two-State Subspace')

grid on

ylim([0 1])

figure

plot(t,traceRho,'LineWidth',1.5)

xlabel('Time')
ylabel('Tr(\rho)')

title('Trace Check')

grid on

ylim([0.99 1.01])


function drhoVec = CL_Lindblad(~,rhoVec,Heff,L,hbar,N)

    rho = reshape(rhoVec,N,N);

    % hamiltonian
    HamiltonianTerm = -(1i/hbar)*(Heff*rho - rho*Heff);

    % dissipator
    dissipator = L*rho*L' ...
        - 0.5*(L'*L*rho + rho*L'*L);

    % Total evolution
    drho = HamiltonianTerm + dissipator;

    drhoVec = drho(:);

end


