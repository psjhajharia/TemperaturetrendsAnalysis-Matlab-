mean_temp=mean()


abc=temp_data(1,1,:)
abc=squezz(temp_data(1,1,:))
sum(1:432)
mean(1:432);
abc=squeeze(temp_data(1,1,:));
meantemp=mean(abc)
abc1=abc-meantemp
test1=zeros(432)
test1=zeros(432,1)
for i=1:432
a(i)=i
end
for i=1:432
test1(i)=i;
end
meandate=mean(test1)
test2=test1-meandate
(sum(abc1.*test2))
new3=(sum(abc1.*test2))^2
new4=sqrt(sum((abc1).^2)*sum((test2).^2))
new3/new4
b=regress(abc);
b=regress(test1,abc);
[p,S] = polyfit(test1,abc,1)
mdl = LinearModel.fit(abc)
mdl = LinearModel.fit(test1,abc);
for i=1:432
a(i)=i
end
for i=1:432
end
for i=1:432
test1(i)=i/10;
end
mdl = LinearModel.fit(test1,abc);