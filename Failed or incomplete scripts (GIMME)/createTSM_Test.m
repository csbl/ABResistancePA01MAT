%function [geneExpression] = testRun()
%Test file for using createTissueSpecificModel (does not work)
clear all
close all

load Data;
%initCobraToolbox
%changeCobraSolver('gurobi');

len = length(Data);
k=1;
cutoff = 75; %percentile 
for i = 1:len  %extracts relevant genes, this may need to be removed...
   if startsWith(Data(i,1),'PA') == 1 
        temp = Data(i,1);
        Data2(k,1)=extractAfter(temp,3);
        Data2(k,2) = Data(i,2);
        k = k+1;
   end
end

geneExpression = struct;
%geneExpression.Locus = double(cell2mat(Data2(:,1)));
geneExpression.Data = cell2mat(Data2(:,2));%convert types
geneExpression.Locus = cellstr(Data2(:,1));


minim = min(geneExpression.Data);
maxim = max(geneExpression.Data)+ -1 *minim;  %find extrema
geneExpression.Data = (geneExpression.Data - minim) ./ maxim;
threshold = prctile(geneExpression.Data,cutoff);
geneExpression.Data = double(geneExpression.Data > threshold); %transform to binary

Model2 = load('myModel.mat');
Model2 = Model2.exported_model;

k=1;
i=1;

%while i <= length(Model2.genes)  
  % if startsWith(Model2.genes(i),'PA') == 1 
        %temp = Model2.genes(i);
       % Model2.genes(k)=extractAfter(temp,3);
       % while startsWith(Model2.genes(k),'PA') == 1 
        %    temp = Model2.genes(k);
        %    Model2.genes(k)=extractAfter(temp,3);
        %    print('Entered double loop')
       % end
     %   k = k+1;
        %disp(k)
   %else
   %     disp(Model2.genes(i))
   %     Model2.genes(i) = [];
  % end
  % i= i+1;
   
   %disp(i)
%end

%if startsWith(Model2.genes(length(Model2.genes)),'PA') == 1 
       % temp = Model2.genes(length(Model2.genes));
        %Model2.genes(k)=extractAfter(temp,3);
%end
        
%Model2.genes = cellstr(Model2.genes)

%save('GeneExpresssion.mat','geneExpression')

%scatter(1:length(geneExpression.Data),geneExpression.Data) %visual verification

%model2 = load('myModel.mat');

%model2 = model2.exported_model;
%initCobraToolbox
model.genes = cellstr(Model2.genes)
model.genes = cell
overlayModel = createTissueSpecificModel(Model2,geneExpression);



%end





