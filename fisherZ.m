function Ztr = fisher(Rho);


Ztr =   triu(.5*log((1+Rho)./(1-Rho)),1);
Ztr = Ztr+tril(Ztr',-1);

end
