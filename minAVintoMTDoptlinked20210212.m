%first want to define and calculate the unswitched system%


residualoptarray= [];
residualoptmat= [];
residualoptcell= {};

define_constants;
mpc=loadcase(case14_20210127);
results = rundcpf(mpc);

%defining topology matrix from matpower%

B=full(makeBdc(mpc));
var=0.0001;
W=(1/var)*eye(20);
D = diag([B(1,2) B(1,5)	B(2,3)	B(2,4)	B(2,5)	B(3,4)	B(4,5)	B(4,7)	B(4,9)	B(5,6)	B(6,11)	B(6,12)	B(6,13)	B(7,8)	B(7,9)	B(9,10)	B(9,14)	B(10,11)	B(12,13)	B(13,14)]);
I=full(makeIncidence(mpc));

%H is the topology matrix%

H = D*I;

%getting the power flows from matpower%

z=results.branch(:, PF);    




%define line capacity%
    
cap = transpose([200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200])

amat = [];
acell = {};
ztestmat = [];



for requiredattackvectors = 1:1 
for branch = 1:20
amat = [];
cmat = [];
%gives us the min required change to load shedd%

co = (abs(cap)-abs(z)).*sign(z)

%change voltage angle Vthi such that we get the needed diff%
zol = z;
zol(branch) = z(branch)+ co(branch) 

c0 = transpose([0,0,0,0,0,0,0,0,0,0,0,0,0,0]);

for busnumber = 1:14
    

ubc = [0,0,0,0,0,0,0,0,0,0,0,0,0,0];

lbc = [0,0,0,0,0,0,0,0,0,0,0,0,0,0];

ubc(busnumber)=100 
lbc(busnumber)=-100

mincfun = @(c)(sum(abs((zol(branch)-z(branch)) - subsref(H*c, struct('type', '()', 'subs', {{branch}})))));
c = fmincon(mincfun,c0,[],[],[],[],lbc,ubc);



cmat = [cmat,c];



for quicketest = 1:1

ztest = z+H*c 
ztestresult=abs(ztest(branch))>=abs(zol(branch))-0.05
a = H*c



if ztestresult==1;
    
acell{busnumber,branch} = a;

amat = [amat,a];



else
end

end



end

    

end

end






for branch = 1:20    
for busnumber = 1:14
    
     a = acell{busnumber,branch}    
     
     if isempty(a)
     else
     
    
    
for MTDupperlimit = 1:0.01:1.15

fun = @(X2opt0)(1-sum(abs((a) -diag(X2opt0(:,:))*I(:,2:14)*(inv((transpose(diag(X2opt0(:,:))*I(:,2:14))*inv(W)*diag(X2opt0(:,:))*I(:,2:14)))*(transpose(diag(X2opt0(:,:))*I(:,2:14))*inv(W))*(a)))));

%starting point no MTD%

X2opt0 = [B(1,2) B(1,5)	B(2,3)	B(2,4)	B(2,5)	B(3,4)	B(4,5)	B(4,7)	B(4,9)	B(5,6)	B(6,11)	B(6,12)	B(6,13)	B(7,8)	B(7,9)	B(9,10)	B(9,14)	B(10,11)	B(12,13)	B(13,14)];
lb = X2opt0*MTDupperlimit;
ub = X2opt0;

    
x0 = [B(1,2) B(1,5)	B(2,3)	B(2,4)	B(2,5)	B(3,4)	B(4,5)	B(4,7)	B(4,9)	B(5,6)	B(6,11)	B(6,12)	B(6,13)	B(7,8)	B(7,9)	B(9,10)	B(9,14)	B(10,11)	B(12,13)	B(13,14)];

X2opt0 = fmincon(fun,x0,[],[],[],[],lb,ub);

residualnormalDVnorm = sum(abs((a) -diag(X2opt0(:,:))*I(:,2:14)*(inv((transpose(diag(X2opt0(:,:))*I(:,2:14))*inv(W)*diag(X2opt0(:,:))*I(:,2:14)))*(transpose(diag(X2opt0(:,:))*I(:,2:14))*inv(W))*((a)))));

residualoptarray= [residualoptarray;residualnormalDVnorm];

end

residualoptcell{busnumber,branch} =residualoptarray
residualoptarray = [];

     end
     
    

end

end
