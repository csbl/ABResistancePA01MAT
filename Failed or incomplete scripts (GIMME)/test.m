    for k = 1:nsamples(i)
        tempModel = CSCobraModels(k);
        fil = sprintf('%s_%d',assession(i),k);
        save(strcat(fil,'.mat'),'tempModel');
        file = fopen(strcat(fil,'.txt'),'w');
        fprintf('%s  %d\n',genes{j}{:,k},round(states{j}(:,k)))
        fclose(file);
    end
    