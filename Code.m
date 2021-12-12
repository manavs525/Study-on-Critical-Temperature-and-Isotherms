A1 = 7.28166;   %for carbon disulphide
B1 = 1446.937;
C1 = 227.6;
R = 8.314;
Tc = 588;
Pc = 600;
P = 0;
v = 0;
w = -2.24;
Q = 0;

i=0;
sigma = 1+sqrt(2);
epsilon = 1-sqrt(2);
alpha = 0;
a = 0.45724*R^2*Tc^2/Pc;
b = 0.0778*R*Tc/Pc;
for T=100:10:800
    Psat = 10^(A1-B1/(T-273.15+C1));
    Psat = Psat*0.133322;
    P = Psat;
    while Q<=0.0001
    %w = -1-log10(P/Pc);
    
    
    
    alpha = (1+(0.37464+1.54226*w-0.26992*w^2)*(1-(T/Tc)^0.5))^2;
    A = alpha*a*P/(R^2*T^2);
    B = b*P/(R^2*T^2);
    k1 = 1;
    k2 = B-1;
    k3 = A-2*B-3*B^2;
    k4 = B^3+B^2-A*B;
    k = [k1 k2 k3 k4];
    Z = roots(k);
    Zliq = min(Z);
    Zvap = max(Z);
    vliq = min(Z)*R*T/P;
    vvap = max(Z)*R*T/P;
    rholiq = 1/vliq;
    rhovap = 1/vvap;
    LHS = -log(1-rholiq*b)-(a*alpha/(2*sqrt(2)*b*R*T))*log((1+sigma*rholiq*b)/(1+epsilon*rholiq*b))+Zliq-1-log(Zliq);
    RHS = -log(1-rhovap*b)-(a*alpha/(2*sqrt(2)*b*R*T))*log((1+sigma*rhovap*b)/(1+epsilon*rhovap*b))+Zvap-1-log(Zvap);
    Q = abs(LHS-RHS);
    
    P = P-0.01;
    
    
    end 
        
    
    
    
    
    j=1;
    for P=50:10:5000
        pressure(j) = P;
    A = alpha*a*P/(R^2*T^2);
    B = b*P/(R^2*T^2);
    k1 = 1;
    k2 = B-1;
    k3 = A-2*B-3*B^2;
    k4 = B^3+B^2-A*B;
    k = [k1 k2 k3 k4];
    Z = roots(k);
    Zliq = min(Z);
    Zvap = max(Z);
    vliq = min(Z)*R*T/P;
    vvap = max(Z)*R*T/P;
    if P < Psat
            v = vvap;
        else 
            v = vliq;            
        end
        volume(j) = v;
    
    j = j+1;
    end
    plot(volume,pressure);
    grid on;
    hold on;
    
end
