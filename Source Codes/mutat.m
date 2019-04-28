function mupop=mutat(crospop,Pm,M,beta)
%% **********************Mutation Process**********************************

%1- for each coromosome: generate a random double number in [0,1]
%2- if 'random number'< Pm then do:
%3- replace 5% of genes with new random number


%% ************************************************************************
len=length(crospop(:,1));
len1=length(crospop(1,:));
len2=1;
for i=1:len
    ra=randi([1 1000])/1000;
    if (ra<Pm)
        ra1=randperm(len1-1);
        ra2=ra1(1:len2);
        if (ra2(1)<=M-1)
            s=floor(sum(crospop(i,1:M-1))-crospop(i,ra2(1)));
            if (s<=253)
            crospop(i,ra2(1))=randi([1,254-s]);
            crospop(i,len1)=0;
            end
        else
            crospop(i,ra2(1))=randi([1,10000*beta])/10000;
            crospop(i,len1)=0;
        end
    end
end
mupop=crospop;

%**************************************************************************
%***************************End Function***********************************
end
%**************************************************************************
%**************************************************************************