function newpop=replace(Pop,mupop,fitval,popsize,Pc)
%% *********************Replacement Process********************************

%select 20% of best choromosomes from last population
%select Pc*Population numbers from mate process(crossover and mutation)
%reminder select randomly from last population

%% ************************************************************************
relite=0.2;
rrand=1-(relite+Pc);

ra=randi([1,popsize],1,round(rrand*popsize));
newpop=Pop;
newpop(1:round(rrand*popsize),:)=(Pop(ra,:));
for i=1:relite*popsize
    [~,I]=max(fitval);
    rar=round(rrand*popsize);
    rar1=rar+i;
    newpop(rar1,:)=Pop(I,:);
    fitval(I)=0;
end
newpop(rrand*popsize+relite*popsize+1:popsize,:)=mupop;

%**************************************************************************
%********************************End Function******************************
end
%**************************************************************************
%**************************************************************************