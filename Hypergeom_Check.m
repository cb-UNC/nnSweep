n       =   0.666;
m           =   1;
gamma       =   1;
etaO        =   10;
etaI        =   0.01;
delta       =   etaO-etaI;
f           =   (m*gamma)^n;

a   =   1;
b   =   2/n;
c   =   (n+1)/n;
d   =   -f;

H = real(hypergeom([a,b],c,d) );

display(H)