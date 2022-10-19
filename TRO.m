function [Qtro, titulo, hint, ho, dhevap, storex, storeh, storeBo, storeCo, storePhi, storeZ] = TRO2(vazao, tempCond, tempOleo, frequencia, estimInicialCoefConvec, numDivisoesTRO, suctionLineEnthalpy)
format shortg
%Dados de entrada
pi = 3.1415;
g = 9.81; %acel gravidade
D = 0.00635; %Diametro do tubo
Ttubo = tempCond;
R = 0.01*0.5; %Raio do rotor
L = 0.25; %comprimento do tubo

Tw = 0.5*(tempOleo + Ttubo);
Tfilme = 0.5*(tempOleo + Tw);
Tfilmeo = 0;

%Inicializacao
x = zeros(2,1); %Vetor solucao

%Propriedades de vapor e liquido
hlcond=(py.CoolProp.CoolProp.PropsSI('H','T',tempCond,'Q',0,'Propane'));
hl=(py.CoolProp.CoolProp.PropsSI('H','T',Ttubo,'Q',0,'Propane'));
hv=(py.CoolProp.CoolProp.PropsSI('H','T',Ttubo,'Q',1,'Propane'));
hlv = hv - hl;
rhov=(py.CoolProp.CoolProp.PropsSI('D','T',Ttubo,'Q',1,'Propane'));
rhol=(py.CoolProp.CoolProp.PropsSI('D','T',Ttubo,'Q',0,'Propane'));
Prl=(py.CoolProp.CoolProp.PropsSI('PRANDTL','T',Ttubo,'Q',0,'Propane'));
kl=(py.CoolProp.CoolProp.PropsSI('conductivity','T',Ttubo,'Q',0,'Propane'));
mul=(py.CoolProp.CoolProp.PropsSI('viscosity','T',Ttubo,'Q',0,'Propane'));
muv=(py.CoolProp.CoolProp.PropsSI('viscosity','T',Ttubo,'Q',1,'Propane'));
sigma=(py.CoolProp.CoolProp.PropsSI('I','T',Ttubo,'Q',0,'Propane'));
cpl=(py.CoolProp.CoolProp.PropsSI('C','T',Ttubo,'Q',0,'Propane'));

vl=1/rhol;
vv=1/rhov;

hmono=0.023*((4*vazao/(pi*D*mul))^0.8)*(Prl^0.4)*kl/D;
Fr = ((4*vazao/(pi*D^2))^2)/(g*D*rhol^2); %numero de Froude do liquido

iter = 0;

div = numDivisoesTRO;
storex = zeros(div,1);
storeh = zeros(div,1);
storeBo = zeros(div,1);
storeCo = zeros(div,1);
storePhi = zeros(div,1);

dL = L/div; %Comprimento de cada trecho

while (abs(Tfilme - Tfilmeo) > 0.01) && (iter < 10)
Tfilmeo = Tfilme;
iter = iter + 1;

%Propriedades do oleo (Mermond, 1999)
DO=1.019*1000+1.994*0.1*Tfilme-(1.318/1000)*Tfilme^2;%934.138 - 0.424401*(tempOleo-273.15) - 0.000986513*(tempOleo-273.15)^2;%852; %densidade do oleo
KO=0.128+0.17*10^(-3)*(60-Tfilme+273.15);%0.124082 + 0.000619921*(tempOleo-273.15) - 0.0000202753*(tempOleo-273.15)^2 + 2.16727E-07*(tempOleo-273.15)^3;%0.138; %Condutividade do oleo
MUO=((9.252*10^(-11))*exp((3.985*1000)/(Tfilme)))*DO;%(28.6029 - 0.73264*(tempOleo-273.15) + 0.00693052*(tempOleo-273.15)^2 - 0.0000222494*(tempOleo-273.15)^3)/1000; %viscosidade em Pa.s
CPO = 4.186*((0.388+0.00045*(1.8*Tfilme+32))/sqrt(DO/998.5))*1000;
PRO=MUO*CPO/KO;%499.3; %Prandtl do oleo
VO=MUO/DO; %Viscosidade cinematica
uo=2*pi*frequencia*R; %velocidade do oleo

%Propriedades avaliadas na temperatura da parede
KOsup=0.128+0.17*10^(-3)*(60-Tw+273.15);%0.124082 + 0.000619921*(tempOleo-273.15) - 0.0000202753*(tempOleo-273.15)^2 + 2.16727E-07*(tempOleo-273.15)^3;%0.138; %Condutividade do oleo
MUOsup=((9.252*10^(-11))*exp((3.985*1000)/(Tw)))*DO;%(28.6029 - 0.73264*(tempOleo-273.15) + 0.00693052*(tempOleo-273.15)^2 - 0.0000222494*(tempOleo-273.15)^3)/1000; %viscosidade em Pa.s
CPOsup = 4.186*((0.388+0.00045*(1.8*Tw+32))/sqrt(DO/998.5))*1000;
PROsup=MUOsup*CPOsup/KOsup; %Prandtl do oleo

%Calculo do h do oleo 
REO=(uo*D)/(VO); %Reynolds
NUO=0.26*(REO^0.6)*(PRO^0.36)*(PRO/PROsup)^0.25; %Zhukauskas
%NUO=0.683*(REO^0.466)*(PRO^0.33); %Hilpert
ho=(NUO*KO)/D; %Coeficiente do oleo

z = 0; %cota da entrada do tubo

% Estimativa inicial e numero maximo de iteracoes
xo = (hlcond - hl)/(hv - hl); %Titulo na entrada do TRO
x(1) = xo + 0.001;
x(2) = hmono;

for j = 1:div
    iterx1 = 0;
    x1old = 0;
    z = z + dL;
%Calcular o novo título de entrada no compresor x1 com base em Shah1982
while ((abs(x(1)-x1old) > 0.00001) && (iterx1 < 2000))
    x1old = x(1);
    iterx1 = iterx1 + 1; %contador
    Co = (((1-x(1))/x(1))^0.8)*((rhov/rhol)^0.5); %Numero de convecao
    
    if Fr > 0.04
        Const = Co;
    else
        Const = 0.38*Co*Fr^(-0.3);
    end
    
    phiCB = 1.8*Const^(-0.8);
    
    iterx2 = 0;
    x2old = 0;
    while (abs(x(2)-x2old)>0.1)&& (iterx2 < 2000)
    x2old = x(2);
    iterx2 = iterx2 + 1;
    Bo = ((x(2)*ho/(x(2)+ho))*(tempOleo-Ttubo))/(4*vazao*hlv/(pi*D^2));
    
    if Bo > 11*10^(-4)
        Fator = 14.7;
    else
        Fator = 15.43;
    end
        if Const > 1.0
            if Bo > 0.3*10^(-4)
                phiNB = 230*Bo^0.5;
            else
                phiNB = 1+46*Bo^0.5;
            end
                if phiNB > phiCB
                    phi = phiNB;
                else
                    phi = phiCB;
                end
        elseif (Const > 0.1) && (Const < 1.0)
            phiBS = Fator*(Bo^0.5)*exp(2.74*Const^(-0.1));
                if phiBS > phiCB
                    phi = phiBS;
                else
                    phi = phiCB;
                end
        else
            phiBS = Fator*(Bo^0.5)*exp(2.47*Const^(-0.15));
                if phiBS > phiCB
                    phi = phiBS;
                else
                    phi = phiCB;
                end
        end
        
        x(2) = phi*hmono; 
    end
        x(1) = xo + (x(2)*ho/(x(2)+ho))*pi*D*dL*(tempOleo-Ttubo)/(vazao*(hv-hl));
end
xo = x(1);
storex(j,1)=x(1);
storeh(j,1)=x(2);
storeBo(j,1)=Bo;
storeCo(j,1)=Co;
storePhi(j,1)=phi;
storeZ(j,1)=z;

end
hsai=(py.CoolProp.CoolProp.PropsSI('H','T',Ttubo,'Q',x(1),'Propane'));
Qtro = vazao*(hsai-hlcond);
Tw = tempOleo - Qtro/(ho*pi*D*L);
Tfilme = 0.5*(tempOleo + Tw);
end
titulo = x(1);
hint = x(2);
dhevap = 1-((hsai - suctionLineEnthalpy)/(hlcond - suctionLineEnthalpy));


