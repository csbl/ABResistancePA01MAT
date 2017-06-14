function [geneNames,Expressor] = GSETXTDataProcessor(fileName,width)

geoData = geoseriesread(fileName);

geneNames = string();

for i = 1:length(geoData.Data.RowNames)
   temp = char(geoData.Data.RowNames(i));
   geneNames(i,1) = temp(1,1:6);
end

condition = 1;
numAddP = 0;
numAdd = 0;

while(condition)
    for i =2:length(geneNames)
        if startsWith(geneNames(i-1),'PA') == 1 && PADownStream(geneNames,i,width) == 1 && startsWith(geneNames(i),'PA') ==0 
            value = str2num(temp(3:length(temp)))+1;
            value = sprintf('%04d',value);
            geneNames(i,1) = string(['P','A',value]); 
            numAdd = numAdd+1;
        end
    end
    if(numAdd - numAddP == 0)
        condition = 0;
    end
    numAddP = numAdd;
    numAdd = 0;
end

Expressor = double(geoData.Data);

end

function [result] = PADownStream(geneNames,i,width)

if (i + width) > length(geneNames)
    result = max(startsWith(geneNames(i+1:length(geneNames)),'PA'));
else 
    result = max(startsWith(geneNames(i+1:i+width),'PA'));

end

if result == 1
    result = true;
else
    result = false;
end
end