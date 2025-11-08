% PA 7

n2 = 2;
[M2, T2, S2, G2] = Nint(1,3,@func82,2);
n4 =4;
[M4, T4, S4, G4] = Nint(1,3,@func82,4);
n8 = 8;
[M8, T8, S8, G8] = Nint(1,3,@func82,8);

exact = 52;

EM2 = abs(exact - M2);
EM4 = abs(exact - M4);
EM8 = abs(exact - M8);

loglog([n2, n4, n8], [EM2, EM4, EM8])

Mm = (log(EM8) - log(EM4));