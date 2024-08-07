function Out = Demod_AC(I1,I2,I3)
A = (I1 - I2).^2;
B = (I2 - I3).^2;
C = (I3 - I1).^2; 
Out = (sqrt(2)/3) * sqrt(A + B + C);
end