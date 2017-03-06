clc
close all
% art = 0.7 * ones(13,1);
% [Aud,Som,Outline,af] = diva_synth(art);


art = [3 * ones(13,1), 0.7 * ones(13,1)];
art(11:13,1) = 1;
art(11:13,2) = 1;

art_ = repmat(art', [40,1]);

[Aud, af] = diva_synth(art_','sound');
plot(Aud)
sound(Aud, 11025)
figure
plot(af)

% art = 0.2 * ones(2,13);
% art(1,11:13) = 1;
% art(2,11:13) = 1;
% [Aud, af] = diva_synth(art','sound');
% plot(Aud)
% sound(Aud, 11025)