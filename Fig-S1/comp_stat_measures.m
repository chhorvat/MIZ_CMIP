clear

z = [-82.5946  -44.8103  -45.5426  -30.3671  -71.6475  -61.2510  -50.0263  -15.7408];
y = [-75.5140  -75.2553  -72.6858];


[tval(1),pval(1),ci(1,:),stat1] = ttest2(z,y,'vartype','equal')

[tval(2),pval(2),ci(2,:),stat2] = ttest2(z,y,'vartype','unequal')


mA = mean(z); 
mB = mean(y); 

sA = std(z); 
sB = std(y); 

nA = length(z); 
nB = length(y); 

num = mA-mB; 

denom_Welch = sqrt(sA^2/nA + sB^2/nB); 

denom_Stud = sqrt(((nA-1)*sA^2 + (nB-1)*sB^2)/(nA+nB-2))*(sqrt((1/nA) + (1/nB))); 

Welch_t = num./denom_Welch; 
Stud_t = num./denom_Stud; 
