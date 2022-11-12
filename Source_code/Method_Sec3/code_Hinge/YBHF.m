function [traindata,trainlabel,testdata,testlabel] =  YBHF(SDPDataPositive,SDPDataNegative);

        [SnPositive,SmPositive] = size(SDPDataPositive);
        arraypositive=randperm(SnPositive);

        testset_positive = SDPDataPositive(arraypositive(1:ceil(1*SnPositive/10)),:);
        trainset_positive  = SDPDataPositive(arraypositive(ceil(1*SnPositive/10)+1:SnPositive),:);

        [SnNegative,SmNegative] = size(SDPDataNegative);
        arraynegative=randperm(SnNegative);

        testset_negative=SDPDataNegative(arraynegative(1:ceil(1*SnNegative/10)),:);
        trainset_negative=SDPDataNegative(arraynegative(ceil(1*SnNegative/10)+1:SnNegative),:);
        
        %训练数据和训练标签
        trainset1 = [trainset_positive; trainset_negative];
        rowrank_train = randperm(size(trainset1, 1)); 
        trainset = trainset1(rowrank_train, :);
        traindata = trainset(:,1:SmPositive-1);
%         for i = 1:size(traindata,2)
%             traindata(:,i)=traindata(:,i)/norm(traindata(:,i));%按列归一化
%         end        %traindata = zscore(traindata);%归一化
        trainlabel= trainset(:,end);
        %trainlabel';
        %测试数据和测试标签
        testset1 = [testset_positive; testset_negative];
        rowrank_test = randperm(size(testset1, 1)); 
        testset = testset1(rowrank_test, :);
        testdata = testset(:,1:size(testset,2)-1);
%         for i = 1:size(testdata,2)
%             testdata(:,i)=testdata(:,i)/norm(testdata(:,i));%按列归一化 %testdata = zscore(testdata);%归一化
%         end
        testlabel= testset(:,end);