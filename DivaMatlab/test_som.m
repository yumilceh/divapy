load soundtest
arts = Arts2;
[aud, som, outline] = diva_synth(arts','audsom');


clearvars -except som arts
som2=som*0;

for i=1:80;
    [~, som_, ~] = diva_synth(arts(i,:)','audsom');
    som2(:,i)=som_;
end

%If diva_syth is called with 3 output arguments returns a different result
%than when it is called with 4 output arguments, the difference is that
%when the function is called with 3 output arguments it uses forward
%computation with previously computed information.

% plot(aud')
% figure
% plot(som')
% figure
% plot(arts)
x=som-som2