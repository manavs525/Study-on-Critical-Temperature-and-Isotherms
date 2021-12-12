clc;
%Defining Constants value in SI Units
R = 8.314;
Tc = 552;
Pc = 7900000;
Omega = 0.1107;
A = 6.94219;
B = 1169.11;
C = 241.59;
Psat = 300000;           %Initial guess Psat for finding Psat for corresponding temperatures

%Arrays to store values of Volume of liquid, vapor and Saturation Pressure
%for corresponding temperature
Vliq = zeros(26,1);
Vvap = zeros(26,1);
PSAT = zeros(26,1);
VOLR = zeros(7,1);
PssR = zeros(7,1);
c = 1;

%for loop to calculate Psat at different temperatures till Tc
for T = 300:10:550
    alpha = 1+(0.37464 + 1.54226*Omega - 0.26992*Omega*Omega)*(1-sqrt(T/Tc));
    a = (0.45724*alpha*(R*Tc)^2)/Pc;
    b = 0.0775*R*Tc/Pc;
    
    %Antoine equation is not being used as the guess value provided by it
    %is giving imaginary values
    %Psat = 10^(A - B/(C+T-273.15));
    %Psat = (Psat*101.325)/760;
    
    
    flag = 0;
    while flag == 0 
        c1 = Psat;
        c2 = (Psat*b - R*T);
        c3 = (a - 2*b*R*T -3*b*b*Psat);
        c4 = (Psat*(b^3) + R*T*b*b - a*b);
        p = [c1 c2 c3 c4];
        V = roots(p);
        if (isreal(V)) == 1          
            Vg = max(V);
            Vl = 10000;
            for i = 1:3
                if V(i)<Vl && V(i)>0
                    Vl = V(i);
                end
            end 
            %fugacity(similar to chemical potential) is calculated for
            %liquid and vapour
            LHS = Psat*exp(Value(Vl,a,b,R,T,Psat));
            RHS = Psat*exp(Value(Vg,a,b,R,T,Psat));
            if abs(LHS - RHS) < 1000
                flag = 1;
            elseif RHS > LHS
                Psat = Psat - 1000;
            else
                Psat = Psat + 1000;
            end
        else
            %if the roots are imaginary
            Psat = Psat + 100000;
        end
    end
    Vliq(c,1) = Vl;
    Vvap(c,1) = Vg;
    PSAT(c,1) = Psat;
    c = c+1;
end

%plot for dome
plot(Vliq,PSAT,'color','blue','LineWidth',2.5);
hold on;
plot(Vvap,PSAT,'color','blue','LineWidth',2.5);
hold on;

%plot for lines inside the dome
for i = 1:26
    x(1) = Vvap(i,1);
    x(2) = Vliq(i,1);
    y = [PSAT(i,1) PSAT(i,1)];
    plot(x,y,'color','red');
    hold on;
end

%plot for lines at T<Tc
T1 = [300:10:550];
for i = 1:4:26
    Ps1 = PSAT(i,1);
    j = 1;
    while(Ps1 <= 1.25*Pc)
        C1 = Ps1;
        C2 = (Ps1*b - R*(T1(i)));
        C3 = (a - 2*b*R*(T1(i)) -3*b*b*Ps1);
        C4 = (Ps1*(b^3) + R*(T1(i))*b*b - a*b);
        p = [C1 C2 C3 C4];
        r = roots(p);
        for k =1:3
            if isreal(r(k))
                VOLL(j) = r(k);
            end
        end
        PssL(j) = Ps1;
        Ps1 = Ps1 + 1000;
        j = j+1;
    end
    plot(VOLL,PssL,'color','red');
    hold on;
end

%plot for Right side of dome
T2 = [300:10:550];
for i = 1:4:26
    Ps2 = PSAT(i,1);
    j = 1;
    while(Ps2 >= 250000)
        C1 = Ps2;
        C2 = (Ps2*b - R*(T2(i)));
        C3 = (a - 2*b*R*(T2(i)) -3*b*b*Ps2);
        C4 = (Ps2*(b^3) + R*(T2(i))*b*b - a*b);
        p = [C1 C2 C3 C4];
        r = roots(p);
        r = real(r);
        VOLR(j) = r(1);
        PssR(j) = Ps2;
        Ps2 = Ps2 - 1000;
        j = j+1;
    end
    plot(VOLR,PssR,'color','red');
    hold on;
end

%plot for T>Tc
for i = 1:0.25:3
    j = 1;
    for Ps = (1.5*Pc):-1000:50000
        C1 = Ps;
        C2 = (Ps*b - R*i*Tc);
        C3 = (a - 2*b*R*i*Tc -3*b*b*Ps);
        C4 = (Ps*(b^3) + R*i*Tc*b*b - a*b);
        p = [C1 C2 C3 C4];
        r = roots(p);
        r = real(r);
        VOLT(j) = r(1);
        PssT(j) = Ps;
        j = j+1;
    end
    plot(VOLT,PssT,'color','red');
    xlim([-0.001 0.01]);
    ylim([0 10^7]);
    title('P-V Isotherm for Carbon Disulfide(Peng Robinson EOS)');
    xlabel('Volume(m^3)');
    ylabel('Pressure(Pa)');
    hold on;
end

%function to calculate chemical potential
function Y = Value(V,a,b,R,T,P)
    A1 = log(abs(1 - b/V));
    B1 = a*log(abs((V + 2.414*b)/(V - 0.414*b)))/(2.828*b*R*T);
    C1 = ((V*P)/(R*T) - 1) - log(abs((V*P)/(R*T)));
    Y = -A1 -B1 +C1;
end 
