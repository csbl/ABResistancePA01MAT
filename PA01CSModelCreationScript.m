%Example Script that uses the function createCSModelTiger to create the CS
%models and save the models as .mat cobra models for analysis in python.
%%

clear all
load tmodel.mat;

assession = cellstr(['GDS4244 ';'GDS3572 ';'GSE90620';'GSE65870';'GSE30021']);
width = 10;
nsamples = [12,6,12,6,9];
normalFact = 75;

[names,expressor,ave,sd,names2] = dataProcessor(assession,nsamples,width,tmodel,normalFact);
%{
%%
parameterResults = zeros(8^4,6);
parameters = [linspace(30,80,4);linspace(0,90,4);linspace(60,150,4);linspace(10,50,4);linspace(10,50,4)]

%%
m =1;
oldscore = 0;
for a = parameters(1,:)
    for b = parameters(2,:)
        for c = parameters(3,:)
            for d = parameters(4,:)
                for e = parameters(5,:)
                    z = [a,b,c,d,e];
                    j = 1;
                    states = cell(5,1);
                    genes = cell(5,1); 
                    allScores = [];

                    for i = 1:length(nsamples)

                        disp(assession(i))

                        [CSCobraModels,errorCount,names3,states{j},genes{j},totalScore] = createCSModelTiger(names{i},expressor{i},nsamples(i),tmodel,sd,z,'width',width,'percentile',ave,'ThreshNames',names2);

                        allScores = [allScores , totalScore];
                        for k = 1:nsamples(i)
                            tempModel = CSCobraModels(k);
                            fil = sprintf('%s_%d',assession{i},k);
                            %save(strcat(fil,'.mat'),'tempModel')
                            file = fopen(strcat(fil,'.txt'),'w');
                            for l = 1:length(states{j})
                                fprintf(file,'%s %d\n',genes{j}{l,k},round(states{j}(l,k)));
                            end
                            %fprintf(file,'%s  %d\n',genes{j}{:,k},round(states{j}(:,k)));
                            fclose(file);
                        end

                        j= j+1;
                    end

                    score = testSum(allScores,states);
                    if score < oldscore || e == parameters(5,length(parameters(5,:)))
                        parameterResults(m,1) = max([score,oldscore]);
                        parameterResults(m,2:end) = [a,b,c,d,e];
                        m = m +1; 
                        break;
                    end
                    oldscore = score;
                end  
            end   
        end   
    end  
end

%%
%}
%save('gimmeSpaceRes.mat','parameterResults')
load gimmeSpaceRes.mat
[n,index] = max(parameterResults(:,1));
z = parameterResults(index,2:6);
j = 1;
states = cell(5,1);
genes = cell(5,1); 
allScores = [];

for i = 1:length(nsamples)

    disp(assession(i))

    [CSCobraModels,errorCount,names3,states{j},genes{j},totalScore] = createCSModelTiger(names{i},expressor{i},nsamples(i),tmodel,sd,z,'width',width,'percentile',ave,'ThreshNames',names2);

    allScores = [allScores , totalScore];
    for k = 1:nsamples(i)
        tempModel = CSCobraModels(k);
        fil = sprintf('%s_%d',assession{i},k);
        %save(strcat(fil,'.mat'),'tempModel')
        file = fopen(strcat(fil,'.txt'),'w');
        for l = 1:length(states{j})
            fprintf(file,'%s %d\n',genes{j}{l,k},round(states{j}(l,k)));
        end
        %fprintf(file,'%s  %d\n',genes{j}{:,k},round(states{j}(:,k)));
        fclose(file);
    end

    j= j+1;
end

score = testSum(allScores,states);
                    
   
%%

controls = [0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1];
assession = cellstr(['GDS4244 ';'GDS3572 ';'GSE90620';'GSE65870';'GSE30021']);
media = [1,0,1,2,1];
types = ["MH","LB","PBM"];

numE = sum(controls);
numC = length(controls) - numE;
k = 1;
genieE = [];
genieC = [];
genieEM = [];
genieCM = [];


%Creates list of all genes that are deleted in experimental groups
for j = 1:length(genes)

    A = states{j} == 0;
    [numrows , numcols] = size(states{j});
    for i = 1 : numcols
       if(controls(k) == 1)
           genieEM = [genieEM ; strcat(genes{j}(A(:,i)),[string(media(j))])];
           genieE = [genieE ; genes{j}(A(:,i))];
       else
           genieCM = [genieCM ; strcat(genes{j}(A(:,i)),string(media(j)))];
           genieC = [genieC ; genes{j}(A(:,i))]; 
       end       
       k = k +1;
    end
end
hitsEM = zeros(length(tmodel.genes),3);
hitsCM = zeros(length(tmodel.genes),3);

for j = 1:3
   for i = 1:length(tmodel.genes)
      hitsEM(i,j) = sum(sum(genieEM==strcat(tmodel.genes(i),string(j-1))));
      hitsCM(i,j) = sum(sum(genieCM==strcat(tmodel.genes(i),string(j-1))));
   end
end

A = sum(hitsEM > 0,2) == 3;
geneDeletionsEM = string(tmodel.genes(logical(A)));
hitsEM2 = zeros(sum(A),3);
for i = 1:3
    temp = hitsEM(:,i);
    hitsEM2(:,i) = temp(A);
end

A = sum(hitsCM > 0,2) == 2;
geneDeletionsCM = string(tmodel.genes(logical(A)));
hitsCM2 = zeros(sum(A),3);
for i = 1:2
    temp = hitsCM(:,i);
    hitsCM2(:,i) = temp(A);
end 


%compress list
hitsE = zeros(length(tmodel.genes),1);
hitsC = zeros(length(tmodel.genes),1);
for i = 1 : length(tmodel.genes)
    hitsE(i) = sum(sum(genieE == tmodel.genes(i)));
    hitsC(i) = sum(sum(genieC == tmodel.genes(i)));
end
geneDeletionsE = tmodel.genes(hitsE > 0);%Names of gene deletions
hitsE = hitsE(hitsE > 0); %get number
geneDeletionsC = tmodel.genes(hitsC > 0);
hitsC = hitsC(hitsC>0);

summary{1} = hitsC;
summary{2} = geneDeletionsC;
summary{3} = hitsE;
summary{4} = geneDeletionsE;

file = [fopen('deletedGenesC.csv','w'),fopen('deletedGenesE.csv','w')];
minim = min([length(geneDeletionsC),length(geneDeletionsE)]);
[maxim,in] = max([length(geneDeletionsC),length(geneDeletionsE)]);
diff = maxim - minim;
for i = 1:minim
    for j = 1:length(file)
        fprintf(file(j),'%s,%d\n',summary{(j-1)*2+2}{i},summary{(j-1)*2+1}(i));
    end
end

for i = minim+1:maxim
    fprintf(file(in),'%s,%d\n',summary{2+2*(in-1)}{i},summary{1+2*(in-1)}(i));
end

    
geneDel = {geneDeletionsCM,geneDeletionsEM};
hits = {hitsCM2,hitsEM2};

file = [fopen('deletedGenesCM.csv','w'),fopen('deletedGenesEM.csv','w')];
minim = min([length(geneDeletionsCM),length(geneDeletionsEM)]);
[maxim,in] = max([length(geneDeletionsCM),length(geneDeletionsEM)]);
diff = maxim - minim;
%for k = 1:minim
    for j = 1:length(file)
        fprintf(file(j),'%s,%s,%s,%s\n',[geneDel{j},string(hits{j})]');
    end
%end

%for k = minim+1:maxim
   % fprintf(file(in),'%s,%s,%s,%s\n',[geneDel{in},string(hits{in})]');

    
fclose('all')
