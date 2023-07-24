function prediction(pInput,pOutput)

rng('default')

K = [1:10];

data = load(sprintf('/nas/longleaf/home/mrcole/Encore/results/pca/%s',pInput));

N_subs = size(data.score,1);
y = [false(N_subs/2,1); true(N_subs/2,1)];

N_REPEATS = 100;

%%

results_mean = zeros(1,length(K));
results_std = zeros(1,length(K));

for k = 1:length(K)
    X = data.score(:,1:K(k));
    cv = cvpartition(y,'KFold',10,'Stratify',true);    

    clasif_error = zeros(1,N_REPEATS);

    for i = 1:N_REPEATS
        cv = repartition(cv);
        idx = 1;

        prediction = zeros(1,N_subs);
        truth = zeros(1,N_subs);

        for j = 1:10
            idxTrain = training(cv,j);
            idxTest = ~idxTrain;
            XTrain = X(idxTrain,:);
            yTrain = y(idxTrain);
            XTest = X(idxTest,:);
            yTest = y(idxTest);        
            
%             [B,FitInfo] = lassoglm(XTrain,yTrain,'binomial','CV',10,'MaxIter',1e6);
%             idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
%             B0 = FitInfo.Intercept(idxLambdaMinDeviance);
%             coef = [B0; B(:,idxLambdaMinDeviance)];
%
%             yhat = glmval(coef,XTest,'logit');                        
        
            mdl = glmfit(XTrain,yTrain,'binomial');
            yhat = glmval(mdl,XTest,'logit');

            len = length(yTest);
            prediction(idx:(idx+len-1)) = yhat;
            truth(idx:(idx+len-1)) = yTest;

            idx = idx + len;
        end 
    
        [u,v,~,auc] = perfcurve(truth,prediction,true);

        clasif_error(i) = auc;
        roc_coords{i,k} = [u,v];
    end    

    disp('RESULT')
    disp(K(k))
    disp(mean(clasif_error))
    results_mean(k) = mean(clasif_error);
    results_std(k) = std(clasif_error); 
end

%save(sprintf('/nas/longleaf/home/mrcole/ScratchingTheSurface/results/%s.mat', pOutput),'results_std','results_mean','roc_coords')

end

% plot(squeeze(roc_coords(:,1,:,1)),squeeze(roc_coords(:,2,:,1)),'color',[0.2 0.5 0.9 0.05]); 
% hold on; 
% plot(mean(squeeze(roc_coords(:,1,:,1)),2),mean(squeeze(roc_coords(:,2,:,1)),2),'Linewidth',2,'Color','blue'); 
% hold off