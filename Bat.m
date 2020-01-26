format long
inp = csvread('musa.csv');
n = size(inp);
n = n(1);
t = inp(:,2);
tn = max(t);
ti = sum(t);
x = 15;
display('Traditional calculation: ');
nbeta = normalcalculation(n,t,x);
integT= 0:10:1000;
upto = size(integT);
for i=1:upto(2)
rel(i) = reliabilitycalculation(nbeta(1),nbeta(2),integT(i));
end
plot(integT,rel,'r');
xlabel('time');
ylabel('reliability');
title('Plot of Reliability Calculation');
display('################################');
display('Solving using BAT Algorithm');
bbeta = batToSolveBETA1(n,t,x)
integT= 0:10:1000;
upto = size(integT);
for i=1:upto(2)
rel1(i) = reliabilitycalculation(bbeta(1),bbeta(2),integT(i));
end
hold on
plot(integT,rel1,'b');
xlabel('time');
ylabel('reliability');
title('Plot of Reliability Calculation');
display('################################');
display('Solving using E-BAT Algorithm');
bbeta = batToSolveBETA2(n,t,x)
integT= 0:10:1000;
upto = size(integT);
for i=1:upto(2)
rel2(i) = reliabilitycalculation(bbeta(1),bbeta(2),integT(i));
end
hold on
plot(integT,rel2,'g');
xlabel('time');
ylabel('reliability');
title('Plot of Reliability Calculation');
format long
result = [integT;rel;rel1;rel2;rel-rel1;rel1-rel];
csvwrite('result.csv',result.');
display('#########NRMSE calculation##############')
nrmse = sqrt(sum(square(rel1 - rel))/sum(square(rel1)));
display('BAt vs traditional')
display(abs(nrmse))
nrmse = sqrt(sum(square(rel2 - rel))/sum(square(rel2)));
display('EBAt vs traditional')
display(abs(nrmse))
nrmse = sqrt(sum(square(rel2 - rel1))/sum(square(rel2)));
display('EBAt vs BAt')
display(abs(nrmse))

function beta = batToSolveBETA1(an,at,ax)
mn=an;
mt = at;
mtn=max(at);
mx = ax;
mti = sum(mt);
para=[40 800 0.5 0.5];
n=para(1);
N_gen=para(2);
A=para(3);
r=para(4);
Qmin=0;
Qmax=2;
N_iter=0;
d=1;
Lb=-2*ones(1,d);
Ub=2*ones(1,d);
Q=zeros(n,1);
v=zeros(n,d);
for i=1:n,
Sol(i,:)=Lb+(Ub-Lb).*rand(1,d);
Fitness(i)=Fun(Sol(i,:));
end
[fmin,I]=min(Fitness);
best=Sol(I,:);
for t=1:N_gen,
for i=1:n,
Q(i)=Qmin+(Qmin-Qmax)*rand;
v(i,:)=v(i,:)+(Sol(i,:)-best)*Q(i);
S(i,:)=Sol(i,:)+v(i,:);
Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
if rand>r
S(i,:)=best+0.001*randn(1,d);
end
Fnew=Fun(S(i,:));
if (Fnew<=Fitness(i)) & (rand<A) ,
Sol(i,:)=S(i,:);
Fitness(i)=Fnew;
end
if Fnew<=fmin,
best=S(i,:);
fmin=Fnew;
end
end
N_iter=N_iter+n;
end
beta(1) = best;
beta(2) = BetaZero(best);
function s=simplebounds(s,Lb,Ub)
ns_tmp=s;
I=ns_tmp<Lb;
ns_tmp(I)=Lb(I);
J=ns_tmp>Ub;
ns_tmp(J)=Ub(J);
s=ns_tmp;
end
function z=Fun(var)
z = mn/var - (mn*(mtn+mx))/(exp(var*(mtn+mx))-1) - mti;
z = abs(z);
end
function F = BetaZero(betaOne)
F = mn/(1-exp(betaOne*-1*(mtn+mx)));
end
end
function beta = batToSolveBETA2(an,at,ax)
mn=an ;
mt = at;
mtn=max(at);
mx = ax;
mti = sum(mt);
para=[20 800 0.5 0.5];
n=para(1);
N_gen=para(2);
A=para(3);
r=para(4);
Qmin=0;
Qmax=2;
N_iter=0;
d=1;
Lb=-2*ones(1,d);
Ub=2*ones(1,d);
Q=zeros(n,1);
v=zeros(n,d);
for i=1:n,
Sol(i,:)=Lb+(Ub-Lb).*rand(1,d);
Fitness(i)=Fun(Sol(i,:));
end
[fmin,I]=min(Fitness);
best=Sol(I,:);
w=rand(1,1);
l1 = 0.6;
l2 = 1-l1;
k = ceil(rand(1,1)*n);
for t=1:N_gen,
for i=1:n,
Q(i)=Qmin+(Qmin-Qmax)*rand;
v(i,:)=w*v(i,:)+(Sol(i,:)-best)*Q(i)*l1+(Sol(i,:)-Sol(k,:))*Q(i)*l2;
S(i,:)=Sol(i,:)+v(i,:);
Sol(i,:)=simplebounds(Sol(i,:),Lb,Ub);
if rand>r
S(i,:)=best+0.001*randn(1,d);
end
Fnew=Fun(S(i,:));
if (Fnew<=Fitness(i)) & (rand<A) ,
Sol(i,:)=S(i,:);
Fitness(i)=Fnew;
end
if Fnew<=fmin,
best=S(i,:);
fmin=Fnew;
end
end
N_iter=N_iter+n;
end
beta(1) = best;
beta(2) = BetaZero(best);
function s=simplebounds(s,Lb,Ub)
ns_tmp=s;
I=ns_tmp<Lb;
ns_tmp(I)=Lb(I);
J=ns_tmp>Ub;
ns_tmp(J)=Ub(J);
s=ns_tmp;
end
function z=Fun(var)
z = mn/var - (mn*(mtn+mx))/(exp(var*(mtn+mx))-1) - mti;
z = abs(z);
end
function F = BetaZero(betaOne)
F = mn/(1-exp(betaOne*-1*(mtn+mx)));
end
end
function answer = normalcalculation(an,at,ax)
n = 0;
t = 0;
tn = 0;
x = 0;
ti = 0;
function assignvalues
n = an;
t = at;
tn= max(at);
x = ax;
ti = sum(t);
display('normalcalculation: values are assignied')
end
function F = likeLIhood(b0,b1)
te1 = b0^n * b1^n;
te2=0;
for i= 1:n
te2 = te2+ exp(-b1*t(i));
end
te3 = exp(-b0*(1-exp(-b1*(tn+x))));
F = te1*te2*te3;
end
function F = BetaZero(betaOne)
F = n/(1-exp(betaOne*-1*(tn+x)));
end
function F = BetaOne(var)
F = n/var - (n*(tn+x))/(exp(var*(tn+x))-1) - ti;
end
function solveIT
display('normalcalculation: solving the values')
x1 = 1;
[x1,fval]= fsolve(@BetaOne,x1);
beta1 = x1
beta0 = BetaZero(beta1)
answer(1) = beta1;
answer(2) = beta0;
end
assignvalues()
solveIT()
end
function rel = reliabilitycalculation(b1,b0,at)
myfun = @(x) b0*b1*exp(-b1*x);
area=integral(myfun,at,at+5);
rel = exp(-area);
end
t=csvread('musa.csv');
xlag = lagmatrix(t(:,2),[0 1 2 3 4 5]);
xlag = xlag(6:101,:);
xlag = [xlag(:,6) xlag(:,5) xlag(:,4) xlag(:,3) xlag(:,2) xlag(:,1) ];
csvwrite('lagcheck.csv',xlag);
nrmse=0;
nrmse1=0;
for lop = 1:5
format long
inp = xlag(:,lop);
n = size(inp);
n = n(1);
t = inp;
tn = max(t);
ti = sum(t);
x = 15;
display('Traditional calculation: ');
display('goes to normal calculation');
nbeta = normalcalculation(n,t,x);
integT= 0:10:1000;
upto = size(integT);
for i=1:upto(2)
rel(i) = reliabilitycalculation(nbeta(1),nbeta(2),integT(i));
end
plot(integT,rel,'r');
xlabel('time');
ylabel('reliability');
title('Plot of Reliability Calculation');
display('################################');
display('Solving using BAT Algorithm');
bbeta = batToSolveBETA1(n,t,x)
integT= 0:10:1000;
upto = size(integT);
for i=1:upto(2)
rel1(i) = reliabilitycalculation(bbeta(1),bbeta(2),integT(i));
end
hold on
%subplot(2,1,2);
plot(integT,rel1,'b');
%legend('normal','bat');
xlabel('time');
ylabel('reliability');
title('Plot of Reliability Calculation');
display('################################');
display('Solving using E-BAT Algorithm');
bbeta = batToSolveBETA2(n,t,x)
integT= 0:10:1000;
upto = size(integT);
for i=1:upto(2)
rel2(i) = reliabilitycalculation(bbeta(1),bbeta(2),integT(i));
end
hold on
plot(integT,rel2,'g');
xlabel('time');
ylabel('reliability');
title('Plot of Reliability Calculation');
format long
result = [integT;rel;rel1;rel2];
csvwrite(strcat('result',int2str(lop),'.csv'),result.');
display(strcat('dumped into: result',int2str(lop),'.csv'));
display('#########NRMSE calculation##############')
nrmse(lop) = abs(sqrt(sum(square(rel - rel2))/sum(square(rel))));
display('EBAt vs tradi')
display(abs(nrmse))
nrmse1(lop) = abs(sqrt(sum(square(rel - rel1))/sum(square(rel))));
display('Bat vs tradi')
display(abs(nrmse1))
display('>>>>>>>>>>>>&&&&&&&&&&&&&<<<<<<<<<<<<<<<<<<<')
end