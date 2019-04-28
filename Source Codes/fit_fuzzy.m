function out_taw=fit_fuzzy(chr_inp,ci,ci1,M)
rule_b=zeros(M*M,5);
t=0;
for i=1:M
    for j=1:M
        t=t+1;
        for k=1:M
            if j==i+M-k
                rule_b(t,:)=[i,j,M-k+1,1,1];
            elseif i==j+M-k
                rule_b(t,:)=[i,j,M-k+1,1,1];
            end
        end
    end
end
a1=newfis('lenfis','mamdani','min','max','min',...
    'max','centroid');
a1=addvar(a1,'input','ci',[0 255]);
a1=addvar(a1,'input','ci1',[0 255]);
a1=addvar(a1,'output','taw',[0 sum(chr_inp(M:2*(M-1)))]);
j=0;
t=0;
while j<M-1
    j=j+1;
    if (j==1)
        a1=addmf(a1,'input',1,'ci','trimf',[0 0 chr_inp(j)]);
        a1=addmf(a1,'input',2,'ci1','trimf',[0 0 chr_inp(j)]);
    else
        a1=addmf(a1,'input',1,'ci','trimf',[t t+chr_inp(j-1) t+chr_inp(j-1)+chr_inp(j)]);
        a1=addmf(a1,'input',2,'ci1','trimf',[t t+chr_inp(j-1) t+chr_inp(j-1)+chr_inp(j)]);
        t=t+chr_inp(j-1);
    end
end
a1=addmf(a1,'input',1,'ci','trapmf',[t t+chr_inp(j) 255 255]);
a1=addmf(a1,'input',2,'ci1','trapmf',[t t+chr_inp(j) 255 255]);
j=0;
t=0;
while j<M-1
    j=j+1;
    if (j==1)
        a1=addmf(a1,'output',1,'taw','trimf',[0 0 chr_inp(j+M-1)]);
    else
        a1=addmf(a1,'output',1,'taw','trimf',[t t+chr_inp(j+M-2) chr_inp(j+M-2)+t+chr_inp(j+M-1)]);
        t=t+chr_inp(j+M-2);
    end
end
a1=addmf(a1,'output',1,'taw','trimf',[t t+chr_inp(j+M-1) t+chr_inp(j+M-1)]);
a1=addrule(a1,rule_b);
[out_taw]=evalfis(double([ci;ci1]),a1);
end