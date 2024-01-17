function params=InitializeFromArt(Art,patternNo,listAmps,Tmax,var0,xfix,xfixind,permutation)
%Gonzalo Mena, 03/2016
load arrayPositions512

params.global.options=optimoptions('fminunc','Algorithm','trust-region','GradObj','on');

params.global.x0=[2.6766    2.6729    1.5639    2.5233    1.9566  -20.7433    0.3893 32.0274];
%params.global.x0=[0.7550 2.2292 0.5780 -2.6504 0.0269 -20.1574 -1.5853 29.9947];
params.global.positions=positions;
params.global.maxIter=5;
params.global.thresEI=30;
params.global.sortData=0;
x0=params.global.x0;
params.global.useStimElec = 0;
options=params.global.options;
stimElecs=patternNo;



    times=[1:Tmax]';
    
    
for j=1:length(times)
for i=1:length(times)

Dif1(i,j)=abs(i-j);
end
end
Dif1=Dif1/(max(max(Dif1)));
Difs{1}=Dif1;    
   
for i=1:512
for j=1:512

Dif2(i,j)=norm(positions(i,:)-positions(j,:),2);
end
end

params.global.Dif1 = Dif1;

 params.global.Dif2 = Dif2;;
  

[theta,rho] = cart2pol(positions(:,1)-positions(stimElecs,1),positions(:,2)-positions(stimElecs,2));

ind=setdiff([1:512],find(rho==0));

rho=rho(ind);


for j=1:length(listAmps)
for i=1:length(listAmps)
Dif3(i,j)=abs(listAmps(i)-listAmps(j));
end
end
Difs{2}=Dif2(ind,ind)/max(max(Dif2(ind,ind)));
Difs{3}=Dif3;
Difs{3}=Difs{3}/max(max(Difs{3}));

Diags{1}=[1:Tmax]'/Tmax;
Diags{2}=rho/max(rho);
Diags{3}=listAmps'/max(listAmps);
     

f1=@(Art,x)logDetKronShift(Art(:,ind,:),[x xfix -100 -100 log(var0)],Difs,setdiff([1:10],union([8 9],xfixind)),[1 1 1],Diags,[3 3 3],permutation);
%f1=@(Art,x)logDetKronShift(Art(:,ind,:),[x(1:7) -100 -100 x(8) log(var0)],Difs,setdiff([1:10],[8 9]),[1 1 1],Diags,[3 3 3],[1:10]); 
g1=@(x)f1(Art,x);
pt(permutation)=1:length(permutation);

x01=[x0(1:7) -100 -100 x0(8)];
x01=x01(pt);
x01=x01(1:8-length(xfix));
 x1 = fminunc(g1,x01,options);
 %x=[x1(1:7) -100 -100 x1(end)];
x=[x1 xfix -100 -100];
x=x(permutation);
    
factp=[0 3 6 9];

for k=1:3
[Ker KerD]=evalKernels(Difs{k},Diags{k},x(factp(k)+1:factp(k+1)),[1 1 1]);
[a b]=eig(Ker);
Q{k}=a';
Qt{k}=a;
dL{k}=diag(b);
Kers{k}=Ker;
end

 
 params.arrayInfo.x= x;
 params.arrayInfo.patternNo =patternNo;
 params.arrayInfo.stimElecs= stimElecs;
 params.patternInfo.Difs=Difs;
 params.patternInfo.Diags=Diags;
 params.patternInfo.Kers=Kers;