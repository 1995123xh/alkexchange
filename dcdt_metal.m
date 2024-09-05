function dcdt_metal = dcdt_metal(t,c, Mbalance, rList, EIE, sf,k, Temp)
% function for calculating the time-derivative of concentration for each
% isotopomer in the alkexchange model

% bond dissociation energies, primary, secondary and tertiary, in kJ/mol
BDE=[100.7 98.2 95.7]*4.18; 
BDE_methane=105*4.18;
DBDE=BDE-BDE_methane;
r=zeros(length(Mbalance),1);  %reaction rates
aH=3;

for i=1:length(r)
    deg=rList(3,i);
    if deg~=0
        BDEfactor=exp(-sf*DBDE(deg)*1000/8.314/Temp);
        if rList(2, i) == 29 % ethane related reactions
            % Higher BDE for methane based on Gribov 2003
            BDEfactor=BDEfactor*exp(sf*-6.4*1000/8.314/Temp); 
        end
    else
        BDEfactor=1;
    end
    if i<length(r)*1/4+1   % adsorption reactions,H
        r(i)=k*rList(4,i)*BDEfactor*c(rList(1,i));
    elseif i<length(r)*1/2+1   % adsorption reactions,D
        r(i)=k*rList(4,i)*BDEfactor*c(rList(1,i))/aH*EIE(rList(1,i));
    elseif i<length(r)*3/4+1  %desorption reactions, H
        r(i)=k*1e-7*BDEfactor*c(rList(2,i))*c(end-1);
    else %desorption reactions, D
        r(i)=k*1e-7*BDEfactor*c(rList(2,i))*c(end)/aH;
    end
end

dcdt_metal=Mbalance*r;

% display(t(end))
end




