list=dir('C:\Users\havm2701\Desktop\New folder (4)\textdata\');
for i=3:length(list)+2
name = list(i).name;    
a=load(name);
name = strrep(name,'.txt','');
truth = a(3,:);

healthyidx=find(truth<.5);
tumoridx=find(.5<truth & truth<=1.5);
edemaidx=find(truth>1.5);
t1 = a(1,:);
t2 = a(2,:);
t1c = t1 - mean(t1);
t2c = t2 - mean(t2);

%t1=downsample(a(1,:),10);
%t2=downsample(a(2,:),10);
%truth=downsample(a(3,:),10);

h=figure(1);
plot(downsample(t1c(healthyidx),20) ,downsample(t2c(healthyidx),20),'.b');
hold on;

plot(downsample(t1c(edemaidx),20),downsample(t2c(edemaidx),20),'.g');
plot(downsample(t1c(tumoridx),20),downsample(t2c(tumoridx),20),'.r');
axis([-1000 2000 -1000 2000])
saveas(h,['C:\Users\havm2701\Dropbox\PhD\Plots\',name,'ms.jpg']);
hold off
h=figure(2);
plot(downsample(t1(healthyidx),20),downsample(t2(healthyidx),20),'.b');
hold on
plot(downsample(t1(tumoridx),20),downsample(t2(tumoridx),20),'.r');
plot(downsample(t1(edemaidx),20),downsample(t2(edemaidx),20),'.g');
axis([-200 2000 -200 2000])
saveas(h,['C:\Users\havm2701\Dropbox\PhD\Plots\',name,'.jpg']);
hold off;

end
