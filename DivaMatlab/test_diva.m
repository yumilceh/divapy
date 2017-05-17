%
%20th March
%

clear
clc
% art = [-2.82586455, 0.6511439, 0.41705491, -0.16636287, 1.52612695, -0.65089551,...
%         2.71660497, -2.8833028, -1.46663703, -2.10624563, 0.21783021, 0.74918289,...
%         0.48187699];
% arts = repmat(art,60,1);
load soundtest
arts = Arts2;

[s,a,b,c] = diva_synth([arts, arts]','audsom');

% sound(s,11025)
% figure,
% plot(s,'r')
% hold on
% plot(Aud2,'b')
% legend('Matlab', 'Python')