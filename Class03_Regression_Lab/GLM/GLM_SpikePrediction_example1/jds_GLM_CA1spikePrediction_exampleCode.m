function jds_GLM_CA1spikePrediction_exampleCode(allanim_CA1resp, allanim_PFCmatrix)

allanim_pg = []; %Compile actual prediction gain across all animals
allanim_pg_s = []; %Compile shuffled prediction gain across all animals

nanim=1;
for a = 1:nanim
    disp(['Processing animal # ' num2str(a)])
    %set animal data
    CA1_resp = allanim_CA1resp{a};
    PFCmatrix = allanim_PFCmatrix{a};
    hpnum = length(CA1_resp);
    
    %GLM with n-fold cross validation
    mse_CA1 = []; %Error compiling for actual data
    s_mse_CA1 = []; %Error compiling for shuffled data
    for i = 1:hpnum
        spk_cnt_str = num2str(CA1_resp{i});
        K = 5;
        cv = cvpartition(spk_cnt_str, 'kfold',K);
        disp(['Cell number ' num2str(i) ' out of ' num2str(hpnum)])% ' - step ' num2str(steps)])
        mse = zeros(K,1);
        shuf_mse = zeros(K,1);
        allshufs = [];
        yhats = [];
        for k=1:K
            % training/testing indices for this fold
            trainIdx = cv.training(k);
            testIdx = cv.test(k);
            
            % train GLM model
            CA1mat = CA1_resp{i};
            CA1mat2 = CA1mat(trainIdx);
            mdl = fitglm(PFCmatrix(trainIdx,:), CA1mat2,'linear','Distribution', 'poisson'); %Shuffle PFCmat2 to determine shuffled data to get error bars
            
            % predict regression output
            Y_hat = predict(mdl, PFCmatrix(testIdx,:));
            
            %Do shuffling
            for s = 1:5000
                shuf = Y_hat(randperm(length(Y_hat)));
                shuf_err(s) = mean(abs(CA1mat(testIdx) - shuf));
            end
            
            % compute mean absolute error
            mse(k) = mean(abs(CA1mat(testIdx) - Y_hat));
            shuf_mse(k) = mean(shuf_err);
            allshufs = [allshufs; shuf_err'];
            yhats = [yhats; (CA1mat(testIdx) - Y_hat)];
        end
        mse_CA1{i}.mse = mse;
        mse_CA1{i}.shuf_mse = shuf_mse;
        
        %SHUFFLE ORIG DATA
        cv2 = cvpartition(spk_cnt_str, 'kfold',K);
        allshufs_s = [];
        yhats_s = [];
        for kk=1:K
            trainIdx2 = cv2.training(kk);
            testIdx2 = cv2.test(kk);
            for p = 1 %Shuffle model just once
                % train GLM model
                %shuffling the order of PFC ripple events
                PFCmatrixshuf = PFCmatrix(randperm(length(PFCmatrix(:,1))),:);
                CA1mat = CA1_resp{i};
                CA1mat2 = CA1mat(trainIdx2);
                mdl2 = fitglm(PFCmatrixshuf(trainIdx2,:), CA1mat2,'Distribution', 'poisson'); %Shuffle PFCmat2 to determine shuffled data to get error bars
                
                % predict regression output
                Y_hat_s = predict(mdl2, PFCmatrix(testIdx2,:));
                
                %Do shuffling
                for s = 1:5000
                    shuf_shuf = Y_hat_s(randperm(length(Y_hat_s)));
                    shuf_shuf_err(s) = mean(abs(CA1mat(testIdx2) - shuf_shuf));
                end
                
                mse_s(kk) = mean(abs(CA1mat(testIdx2) - Y_hat_s));
                shuf_mse_s(kk) = mean(shuf_shuf_err);
                allshufs_s = [allshufs_s; shuf_shuf_err'];
                yhats_s = [yhats_s; (CA1mat(testIdx2) - Y_hat_s)];
            end
        end
        s_mse_CA1{i}.mse_s = mse_s;
        s_mse_CA1{i}.shuf_mse_s = shuf_mse_s;
    end
    
    %Calculate prediction gain for actual data
    pred_gain = [];
    for i = 1:length(mse_CA1)
        mse = mean(mse_CA1{i}.mse);
        shuf_mse = mean(mse_CA1{i}.shuf_mse);
        pg_log = log10(shuf_mse/mse);
        pred_gain = [pred_gain pg_log];
    end
    allanim_pg = [allanim_pg; pred_gain'];
    
    %Now go through all shuffled models and calculate prediction gain
    pred_gain_s = [];
    for i = 1:length(s_mse_CA1)
        mse_s = mean(s_mse_CA1{i}.mse_s);
        shuf_mse_s = mean(s_mse_CA1{i}.shuf_mse_s);
        pg_log_s = log10(shuf_mse_s/mse_s);
        pred_gain_s = [pred_gain_s pg_log_s];
    end
    allanim_pg_s = [allanim_pg_s; pred_gain_s'];
    
    disp(['Animal # ' num2str(a) ' processing complete'])
end

%%
%Not a native matlab function
% permutationTest(allanim_pg, allanim_pg_s, 10000, 'plotresult', 1,'sidedness','larger');
%%
mean_data = mean(allanim_pg);
mean_shuf = mean(allanim_pg_s);
data_sem = std(allanim_pg)/sqrt(size(allanim_pg,1));
shuf_sem = std(allanim_pg_s)/sqrt(size(allanim_pg_s,1));
figure
combdata = [mean_data mean_shuf];
combsem = [data_sem shuf_sem];
bar(combdata,'k'); hold on;
er = errorbar([1 2],combdata,combsem);
er.Color = 'k';
er.LineStyle = 'none';
er.LineWidth = 2;
xticklabels({'Data','Shuffled'})

keyboard

