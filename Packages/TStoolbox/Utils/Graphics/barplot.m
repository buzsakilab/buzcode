function [hb, he] = barplot(data,gca)
    
N = length(data);

mu = zeros(1,N);
er = zeros(1,N);

for ii=1:N
    mu(ii)=nanmean( data{ii}); 
    er(ii)=nansem( data{ii}); 
end

bar(mu); %,'b','FaceColor',[0.5,0.5,1],'LineWidth',1.5)
hold on
for ii=1:N
  plot([ii,ii],[mu(ii)-er(ii),mu(ii)+er(ii)],'-k','LineWidth',4)
end
