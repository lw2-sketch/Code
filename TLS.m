hbar = 1;

energy = 0.2;   
Delta = 1;      

H = [-energy/2, -Delta/2;
     -Delta/2,   energy/2];


sigmaz = [1,0;0,-1];
gamma = 0.1;





rho_nought = [1, 0; 0, 0];


rhovec_nought = rho_nought(:);

tspan = 0:0.01:100;


[t, rhovec] = ode45(@(t,rhovec) ...
    densityEquation(t,rhovec,H,hbar,gamma,sigmaz), ...
    tspan,rhovec_nought);
n = length(t);


rho1 = zeros(n,1);
rho2 = zeros(n,1);
coherence = zeros(n,1);


for i = 1:n

    rho = reshape(rhovec(i,:),2,2);

    rho1(i) = real(rho(1,1));
    rho2(i) = real(rho(2,2));

    coherence(i) = abs(rho(1,2));

end


figure

plot(t,rho1,'LineWidth',1.5)
hold on
plot(t,rho2,'LineWidth',1.5)

xlabel('Time')
ylabel('Population')
legend('State 1','State 2')
grid on


figure

plot(t,coherence,'LineWidth',1.5)

xlabel('Time')
ylabel('|rho_{12}|')
title('Coherence')
grid on



function drhovec = densityEquation(~,rhovec,H,hbar,gamma,sigmaz)
    rho = reshape(rhovec,2,2);
    hamiltonianiterm = -(1i/hbar)*(H*rho - rho*H);
    dephasingterm = (gamma/2)*(sigmaz*rho*sigmaz - rho);
    drho = hamiltonianiterm + dephasingterm;
    
    drhovec = drho(:);


    

end





