function crospop=crosov(selpop,M)
%% ************************Crossover Function******************************

% two point crossover do:
%1-product two random number in length of choromosomes boundry
%2-divide two consecutive choromosomes by those two number
%3-displace each part of choromosomes by each other
% to save ascendent order of genes reordered those ascendent

%% ************************************************************************
len=length(selpop(:,1));
len1=length(selpop(1,:));
i=0;
newpop=zeros(len,len1);
l=len;
while i<len
    if l~=1
        p=randi([1,len1-2],1,1);
        i=i+1;
        newpop(i,:)=[selpop(i+1,1:p(1)) selpop(i,p(1)+1:len1)];
        newpop(i+1,:)=[selpop(i,1:p(1)) selpop(i+1,p(1)+1:len1)];
        newpop(i,len1)=0;
        newpop(i+1,len1)=0;
        if (sum(newpop(i,1:M-1))>=254)
            newpop(i,:)=selpop(i,:);
            newpop(i,len1)=selpop(i,len1);
        end
        if (sum(newpop(i+1,1:M-1))>=254)
            newpop(i+1,:)=selpop(i+1,:);
            newpop(i+1,len1)=selpop(i+1,len1);
        end        
        i=i+1;
    else
        i=i+1;
        newpop(i,:)=selpop(i,:);
    end
    l=l-2;
end
crospop=newpop;

%**************************************************************************
%******************************End Function********************************
end
%**************************************************************************
%**************************************************************************