start_tiger();

model = load('myModel.mat');
model = model.exported_model;


tmodel = cobra_to_tiger(model);
save('tmodel.mat','tmodel')