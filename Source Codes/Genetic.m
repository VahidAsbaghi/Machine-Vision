function [pop,I,M_M,table1]=Genetic(im_inp,kk,table,ci,ci1,method,M_M)

if method==1
    taw=2;
    beta=1.2;
elseif method==2
    taw=0.4;
    beta=0.1;
else
    taw=2;
    beta=0.5;
end
pop_size=40;
MR=0.25;
CR=0.8;

temp=zeros(1,M_M);

pop=zeros(pop_size,2*(M_M-1)+1);
for i=1:pop_size
    temp(1,2:M_M)=randi([1,254],[1,M_M-1]);
    temp=sort(temp);
    tempt=temp(2:M_M)-temp(1:M_M-1);
    pop(i,1:M_M-1)=tempt;
end

for i=1:pop_size
    temp1(1,2:M_M)=randi([1,10000*taw],[1,M_M-1])./100000;
    temp1=sort(temp1);
    tempt1=temp1(2:M_M)-temp1(1:M_M-1);
    pop(i,M_M:2*(M_M-1))=tempt1;%1/(M_M-1);
    pop(i,2*(M_M-1)+1)=0;
end
fit_nes=zeros(pop_size,1);
out_taw=zeros(pop_size,1);

for j=1:30
    for i=1:pop_size
        chr_inp=pop(i,:);
        [fit_nes(i),~,~]=interpol_pix(im_inp,M_M,kk,chr_inp,0,ci,ci1,method);
        pop(i,2*M_M-1)=fit_nes(i);
    end
    selpop=selectn(fit_nes,CR,pop_size,pop);
    crospop=crosov(selpop,M_M);
    mupop=mutat(crospop,MR,M_M,beta);
    newpop=replace(pop,mupop,fit_nes,pop_size,CR);
    pop=newpop;
end
for i=1:pop_size
    chr_inp=pop(i,:);
    [fit_nes(i),out_taw(i),~]=interpol_pix(im_inp,M_M,kk,chr_inp,0,ci,ci1,method);
    pop(i,2*M_M-1)=fit_nes(i);
end
[~,I]=max(fit_nes);

table(ci+1,ci1+1)=out_taw(I);
table1=table;
end
