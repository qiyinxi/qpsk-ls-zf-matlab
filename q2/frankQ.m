function [FrankZ]=frankQ(U)
 
for q = 0:(sqrt(U)-1)
    for p = 0:(sqrt(U)-1)
        F(p+q*sqrt(U)+1) = 2*pi*p*q/sqrt(U);
    end
end
 
I = cos(F);
Q = sin(F);
FrankZ = I + 1i*Q;
end
