function [Aud,Som,Outline,af,filt]=diva_synth(Art,option)
if nargin<2, if size(Art,2)>1, option='sound'; else option='audsom'; end; end
% Art(1:10) vocaltract shape params
% Art(11:13) F0/P/V params

Art=tanh(Art);%max(-1,min(1,Art));
switch(lower(option))
    case 'explicit' 
        [Aud,Som,Outline,af,d]=diva_synth_sample(Art);
        filt=a2h(max(0,af),d,1000,11025);
    case 'sound' % outputs soundwave associated with sequence of articulatory states
        [Aud, Som]=diva_synth_sound(Art);   %Aud is the sound wave and Som is the area function
    case 'af'
        Aud = diva_synth_af(Art);
           
    case 'audsom' % outputs auditory/somatosensory representation associated with a given articulatory state
        ndata=size(Art,2);
        if ndata>1
            Aud=cell(1,ndata);Som=cell(1,ndata);Outline=cell(1,ndata);af=cell(1,ndata);
            for n1=1:size(Art,2),
                [Aud{n1},Som{n1},Outline{n1}]=diva_synth_sample(Art(:,n1));
            end
            Aud=cat(2,Aud{:});
            Som=cat(2,Som{:});
            Outline=cat(2,Outline{:});
             if nargout>3,
                 laf=cellfun('length',af);
                 mlaf=max(laf);
                 af=cellfun(@(x)x(min(numel(x),1:mlaf)),af,'UniformOutput',false);
                 af=cat(2,af{:});
                 display(af);
             end
        else
            if nargout>3
                [Aud,Som,Outline,af]=diva_synth_sample(Art);
            elseif nargout>1
                [Aud,Som,Outline]=diva_synth_sample(Art);
            else
                Aud=diva_synth_sample(Art);
            end
        end
end
end

% synthesizes sound
function [s af]=diva_synth_sound(Art)
% global DIVA_x;
persistent vt;
if isempty(vt)
    [filepath,filename]=fileparts(mfilename);
    load(fullfile(filepath,[filename,'.mat']),'vt');
end

synth=struct('fs',11025,'update_fs',200); %Modify sample frequency
synth.f0=120;
synth.samplesperperiod=ceil(synth.fs/synth.f0);
synth.glottalsource=glotlf(0,(0:1/synth.samplesperperiod:1-1/synth.samplesperperiod)');
synth.f=[0,1];
synth.filt=[0,0];
synth.pressure=0;
%synth.modulation=1;
synth.voicing=1;
synth.pressurebuildup=0;
synth.pressure0=0;
synth.sample=zeros(synth.samplesperperiod,1);
synth.k1=1;
synth.numberofperiods=1;
synth.samplesoutput=0;

vt.idx=1:10;
vt.pressure=0;
vt.f0=120;
vt.closed=0;
vt.closure_time=0;
vt.closure_position=0;
vt.opening_time=0;

voices=struct('F0',{120,340},'size',{1,.7});
opt.voices=1;

ndata=size(Art,2);
dt=.005;
time=0;
s=zeros(ceil((ndata+1)*dt*synth.fs),1);
while time<(ndata+1)*dt;
    % sample articulatory parameters
    t0=floor(time/dt);
    display(t0)
    t1=(time-t0*dt)/dt;
    [nill,nill,nill,af1,d]=diva_synth_sample(Art(:,min(ndata,1+t0)));
    [nill,nill,nill,af2,d]=diva_synth_sample(Art(:,min(ndata,2+t0)));
    naf1=numel(af1);naf2=numel(af2);
    if naf2<naf1,af2(end+(1:naf1-naf2))=af2(end); end
    if naf1<naf2,af1(end+(1:naf2-naf1))=af1(end); end
    af=af1*(1-t1)+af2*t1;
    FPV=max(-1,min(1, Art(end-2:end,min(ndata,1+t0))*(1-t1)+Art(end-2:end,min(ndata,2+t0))*t1 ));
    vt.voicing=(1+tanh(3*FPV(3)))/2;
    vt.pressure=FPV(2);
    vt.pressure0=vt.pressure>.01;
    vt.f0=100+20*FPV(1);
    
    af0=max(0,af);
    k=.025;af0(af0>0&af0<k)=k;
    minaf=min(af);
    minaf0=min(af0);
    vt.af=af;
%      display(af);
%      DIVA_x.af_sample=DIVA_x.af_sample+1;
%      DIVA_x.af(:,DIVA_x.af_sample)=af;
%    tracks place of articulation
    if minaf0==0, 
        release=0;
        vt.opening_time=0; vt.closure_time=vt.closure_time+1;
        vt.closure_position=find(af0==0,1,'last');
        if ~vt.closed, closure=vt.closure_position; else closure=0; end;
        vt.closed=1;
    else
        if vt.closed, release=vt.closure_position; release_closure_time=vt.closure_time; else release=0; end;
        if (vt.pressure0&&~synth.pressure0) vt.opening_time=0; end;
        vt.opening_time=vt.opening_time+1;
        vt.closure_time=0;
        [nill,vt.closure_position]=min(af);
        closure=0;
        vt.closed=0;
    end
%     display(vt.closed);
    if release>0  af=max(k,af);minaf=max(k,minaf);minaf0=max(k,minaf0); end
    
    if release>0, 
                    vt.f0=(.95+.0*rand)*voices(opt.voices).F0; %MODIFIED 0.1 Noise
                    synth.pressure=0;%modulation=0; 
    elseif  (vt.pressure0&&~synth.pressure0) 
                    vt.f0=(.95+.0*rand)*voices(opt.voices).F0;  %MODIFIED 0.1 Noise
                    synth.pressure=vt.pressure; synth.f0=1.25*vt.f0; 
                    synth.pressure=1;%synth.modulation=1; 
    elseif  (~vt.pressure0&&synth.pressure0&&~vt.closed), synth.pressure=synth.pressure/10;
    end
    
    % computes glottal source
    synth.samplesperperiod=ceil(synth.fs/synth.f0);
    pp=[.6,.2-.1*synth.voicing,.25];%10+.15*max(0,min(1,1-vt.opening_time/100))];
    synth.glottalsource=10*.25*glotlf(0,(0:1/synth.samplesperperiod:1-1/synth.samplesperperiod)',pp)+10*.025*synth.k1*glotlf(1,(0:1/synth.samplesperperiod:1-1/synth.samplesperperiod)',pp);
    numberofperiods=synth.numberofperiods;
        
    % computes vocal tract filter
    [synth.filt,synth.f,synth.filt_closure]=a2h(af0,d,synth.samplesperperiod,synth.fs,vt.closure_position,minaf0);
    synth.filt=2*synth.filt/max(eps,synth.filt(1));
    synth.filt(1)=0;
    synth.filt_closure=2*synth.filt_closure/max(eps,synth.filt_closure(1));
    synth.filt_closure(1)=0;
    
    % computes sound signal
    w=linspace(0,1,synth.samplesperperiod)';
    if release>0,%&&synth.pressure>.01,
        u=synth.voicing*1*.010*(synth.pressure+20*synth.pressurebuildup)*synth.glottalsource + (1-synth.voicing)*1*.010*(synth.pressure+20*synth.pressurebuildup)*0*randn(synth.samplesperperiod,1);
%         if release_closure_time<40
%             u=1*.010*synth.pressure*synth.glottalsource;%.*(0.25+.025*randn(synth.samplesperperiod,1)); % vocal tract filter
%         else
%             u=1*.010*(synth.pressure+synth.pressurebuildup)*randn(synth.samplesperperiod,1);
%         end
        v0=real(ifft(fft(u).*synth.filt_closure));
        numberofperiods=numberofperiods-1;
        synth.pressure=synth.pressure/10;
        vnew=v0(1:synth.samplesperperiod);
        v0=(1-w).*synth.sample(ceil(numel(synth.sample)*(1:synth.samplesperperiod)/synth.samplesperperiod))+w.*vnew;
        synth.sample=vnew;        
    else v0=[]; end
    if numberofperiods>0,
        %u=0.25*synth.modulation*synth.pressure*synth.glottalsource.*(1+.1*randn(synth.samplesperperiod,1)); % vocal tract filter
        u=0.25*synth.pressure*synth.glottalsource.*(1+.0*randn(synth.samplesperperiod,1)); % vocal tract filter
        u=(synth.voicing*u+(1-synth.voicing)*.0*synth.pressure*randn(synth.samplesperperiod,1));
        if minaf0>0&&minaf0<=k, u=minaf/k*u+(1-minaf/k)*.0*synth.pressure*randn(synth.samplesperperiod,1); end
        v=real(ifft(fft(u).*synth.filt));
        
        vnew=v(1:synth.samplesperperiod);
        v=(1-w).*synth.sample(ceil(numel(synth.sample)*(1:synth.samplesperperiod)'/synth.samplesperperiod))+w.*vnew;
        synth.sample=vnew;
        
        if numberofperiods>1
            v=cat(1,v,repmat(vnew,[numberofperiods-1,1]));
        end
    else v=[]; end
    
    display(time)
    display(size(v0))
    display(size(v))
    v=cat(1,v0,v);
    v=v+.000*randn(size(v));
    v=(1-exp(-v))./(1+exp(-v));
    s(synth.samplesoutput+(1:numel(v)))=v;
    time=time+numel(v)/synth.fs;
    synth.samplesoutput=synth.samplesoutput+numel(v);
    
    % computes f0/amp/voicing/pressurebuildup modulation
    synth.pressure0=vt.pressure0;
    alpha=min(1,(.1)*synth.numberofperiods);
    beta=100/synth.numberofperiods;
    synth.pressure=synth.pressure+alpha*(vt.pressure*(max(1,1.5-vt.opening_time/beta))-synth.pressure);
    alpha=min(1,.5*synth.numberofperiods);beta=100/synth.numberofperiods;
    synth.f0=synth.f0+2*sqrt(alpha)*0*randn+alpha*(vt.f0*max(1,1.25-vt.opening_time/beta)-synth.f0);%147;%120;
    synth.voicing=max(0,min(1, synth.voicing+.5*(vt.voicing-synth.voicing) ));
    %synth.modulation=max(0,min(1, synth.modulation+.1*(2*(vt.pressure>0&&minaf>-k)-1) ));
    alpha=min(1,.1*synth.numberofperiods);
    synth.pressurebuildup=max(0,min(1, synth.pressurebuildup+alpha*(2*(vt.pressure>0&minaf<0)-1) ));
    synth.numberofperiods=max(1,numberofperiods);
end
s=s(1:ceil(synth.fs*ndata*dt));
end

% computes auditory/somatosensory representations
% Art(1:10) vocaltract shape params
% Art(11:13) F0/P/V params
% Aud(1:4) F0-F3 pitch&formants
% Som(1:6) place of articulation (~ from pharyngeal to labial closure)
% Som(7:8) P/V params (pressure,voicing)
function [Aud,Som,Outline,af,d]=diva_synth_sample(Art)
persistent vt fmfit;
if isempty(vt)
    [filepath,filename]=fileparts(mfilename);
    load(fullfile(filepath,[filename,'.mat']),'vt','fmfit');
end
% computes vocal tract configuration
idx=1:10;
x=vt.Scale(idx).*Art(idx);
Outline=vt.Average+vt.Base(:,idx)*x;
% computes somatosensory output (explicitly from vocal tract configuration)
Som=zeros(8,1);
if nargout>3
    [a,b,sc,af,d]=xy2ab(Outline);
    Som(1:6)=max(-1,min(1, -tanh(1*sc) ));
    Som(7:8)=Art(end-1:end);
%         [a,b,sc]=xy2ab(Outline);
end
% computes auditory/somatosensory output (through previously computed forward fit)
Aud=zeros(4,1);
if ~isempty(fmfit)
    Aud(1)=100+50*Art(end-2);
    dx=bsxfun(@minus,Art(idx)',fmfit.mu);
    p=-sum((dx*fmfit.iSigma).*dx,2)/2;
    p=fmfit.p.*exp(p-max(p));
    p=p/sum(p);
    px=p*[Art(idx)',1];
    Aud(2:4)=fmfit.beta_fmt*px(:);
end
if ~isempty(fmfit)&&nargout>1&&nargout<=3,
    Som(1:6)=fmfit.beta_som*px(:);
    Som(7:8)=Art(end-1:end);
end
end

% computes area function
function [a,b,sc,af,d]=xy2ab(x,y)
persistent ab_alpha ab_beta;
if isempty(ab_alpha),
    %     ab_alpha=[1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.79,1.78,1.78,1.78,1.77,1.77,1.76,1.75,1.74,1.73,1.72,1.71,1.69,1.67,1.65,1.63,1.60,1.58,1.55,1.52,1.49,1.46,1.43,1.40,1.37,1.34,1.31,1.29,1.26,1.24,1.22,1.20,1.19,1.17,1.16,1.16,1.15,1.15,1.15,1.15,1.16,1.17,1.17,1.19,1.20,1.21,1.23,1.24,1.25,1.27,1.28,1.30,1.31,1.32,1.34,1.35,1.36,1.37,1.37,1.38,1.38,1.38,1.39,1.39,1.39,1.39,1.39,1.39,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.37,1.37,1.37,1.37,1.37,1.36,1.36,1.36,1.36,1.36,1.35,1.35,1.35,1.35,1.35,1.35,1.36,1.36,1.37,1.37,1.38,1.39,1.40,1.41,1.43,1.44,1.46,1.47,1.49,1.51,1.53,1.55,1.57,1.60,1.62,1.64,1.66,1.69,1.71,1.73,1.75,1.77,1.79,1.82,1.84,1.87,1.90,1.94,1.98,2.02,2.07,2.13,2.19,2.25,2.32,2.40,2.48,2.56,2.65,2.75,2.85,2.95,3.05,3.16,3.27,3.37,3.48,3.59,3.69,3.79,3.89,3.99,4.07,4.16,4.24,4.31,4.38,4.44,4.49,4.54,4.58,4.62,4.64,4.67,4.69,4.70,4.71,4.71,4.72,4.72,4.72,4.72,4.72,4.72,4.72,4.72,4.72,4.72,4.72,4.72,4.72,4.72,4.72,4.72,4.72]';
    %     ab_beta=[1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.38,1.39,1.39,1.39,1.40,1.40,1.40,1.41,1.42,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.51,1.52,1.53,1.54,1.55,1.55,1.56,1.56,1.57,1.57,1.57,1.56,1.56,1.55,1.54,1.53,1.52,1.50,1.48,1.46,1.44,1.42,1.40,1.38,1.36,1.33,1.31,1.29,1.26,1.24,1.22,1.20,1.19,1.17,1.15,1.14,1.13,1.12,1.11,1.10,1.10,1.09,1.09,1.09,1.09,1.10,1.10,1.11,1.12,1.12,1.13,1.14,1.15,1.17,1.18,1.19,1.21,1.22,1.24,1.25,1.27,1.29,1.30,1.32,1.34,1.35,1.37,1.38,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.46,1.47,1.47,1.47,1.47,1.46,1.46,1.45,1.45,1.44,1.43,1.42,1.41,1.40,1.38,1.37,1.36,1.35,1.34,1.33,1.31,1.30,1.29,1.28,1.27,1.27,1.26,1.26,1.26,1.27,1.27,1.28,1.30,1.32,1.34,1.36,1.39,1.42,1.46,1.50,1.54,1.58,1.62,1.67,1.72,1.77,1.82,1.86,1.91,1.96,2.01,2.06,2.10,2.14,2.19,2.22,2.26,2.29,2.32,2.35,2.38,2.40,2.42,2.43,2.45,2.46,2.46,2.47,2.47,2.48,2.48,2.48,2.48,2.48,2.48,2.48,2.48,2.48,2.48,2.48,2.48,2.48,2.48,2.48,2.48,2.48,2.48]';
    %     ab_alpha=[0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.26,0.26,0.26,0.27,0.28,0.29,0.30,0.31,0.33,0.34,0.36,0.38,0.40,0.42,0.45,0.47,0.50,0.53,0.55,0.58,0.61,0.64,0.67,0.70,0.72,0.75,0.78,0.80,0.83,0.85,0.87,0.89,0.91,0.92,0.94,0.95,0.96,0.97,0.98,0.99,0.99,0.99,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.01,1.01,1.02,1.03,1.04,1.05,1.06,1.08,1.10,1.12,1.15,1.17,1.20,1.23,1.26,1.30,1.33,1.37,1.40,1.44,1.48,1.52,1.56,1.60,1.63,1.67,1.70,1.74,1.77,1.81,1.84,1.88,1.92,1.95,1.99,2.04,2.08,2.13,2.17,2.23,2.28,2.34,2.40,2.46,2.52,2.59,2.66,2.73,2.81,2.88,2.96,3.04,3.12,3.19,3.27,3.34,3.41,3.48,3.54,3.60,3.66,3.71,3.76,3.80,3.84,3.87,3.90,3.93,3.95,3.96,3.98,3.98,3.99,4.00,4.00,4.00,4.00,4.00,4.00,4.00,4.00,4.00,4.00,4.00,4.00,4.00,4.00,4.00,4.00,4.00,4.00]';
    %     ab_beta=[1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.25,1.26,1.26,1.26,1.27,1.28,1.29,1.30,1.31,1.33,1.34,1.36,1.38,1.40,1.42,1.45,1.47,1.50,1.53,1.55,1.58,1.61,1.64,1.67,1.70,1.72,1.75,1.78,1.80,1.83,1.85,1.87,1.89,1.91,1.92,1.94,1.95,1.96,1.97,1.98,1.99,1.99,1.99,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00]';
    %     alpha=[1.79, 1.34, 0.73, 1.39, 1.34, 1.92, 4.72];
    %     beta=[1.38, 1.62, 1.81, 1.08, 1.51, 1.20, 2.48];
    %alpha=[.2, 1, 1, 1, 1, 1, 2.5];
    %beta=[.125, 1.25, 1.25, 1.25, 1.25, 1.25, 2.5];
    amax=220;
    alpha=[1, 1, 1, 1, 1, 1, 1];
    beta=[.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25];
    idx={1:60,61:70,71:80,81:120,121:150,151:190,191:amax};
    ab_alpha=zeros(amax,1);
    ab_beta=zeros(amax,1);
    for n1=1:numel(idx),
        ab_alpha(idx{n1})=alpha(n1);
        ab_beta(idx{n1})=beta(n1);
    end; 
    h=hanning(51)/sum(hanning(51));
    ab_alpha=convn(ab_alpha([ones(1,25),1:end,end+zeros(1,25)]),h,'valid');
    ab_beta=convn(ab_beta([ones(1,25),1:end,end+zeros(1,25)]),h,'valid');
    %sprintf('%0.2f,',ab_alpha(:)')
    %sprintf('%0.2f,',ab_beta(:)')
end


if nargin==1, x=exp(-1i*pi/12)*x; y=imag(x);x=real(x); end
% Grid
x0=45;%90;
y0=-100;%-60;
r=60;%30;
k=pi*r/2;
d=.75/10;%unitstocm

a=zeros(size(x));
b=zeros(size(x));
i1=y<y0;
i2=x<x0;
i3=y>=y0&x>=x0;

% a,b: "linearized" coordinates along vocal tract
a(i1)=y(i1)-y0;
b(i1)=x(i1)-x0;
a(i2)=k+x0-x(i2);
b(i2)=y(i2)-y0;
z=x(i3)-x0+1i*(y(i3)-y0);
a(i3)=r*angle(z);
b(i3)=abs(z);

if nargout>2,
    % tube area
    
    olips=30:45;
    ilips=257:302;
    owall=45:164;
    iwall=164+10:257;
    oall=30:164;
    iall=164:302;
    xmin=-20;ymin=-160;%xmin=0;ymin=-140;
    amin=ymin-y0;amax=ceil((x0-xmin+k-amin));
    
    fact=3;
    wallab1=accumarray(max(1,min(fact*9, ceil(fact*9*(a(oall)-amin)/amax))),b(oall),[fact*9,1],@min,nan);
    wallab2=accumarray(max(1,min(fact*9, ceil(fact*9*(a(iwall)-amin)/amax))),b(iwall),[fact*9,1],@max,nan);
    lipsab1=min(b(olips));
    lipsab2=max(b(ilips));
    mind=min(reshape(wallab1(fact*2+1:fact*8)-wallab2(fact*2+1:fact*8),[fact,6]),[],1);
    sc=d*[mind(1:4)';min(mind(5:6));lipsab1-lipsab2];
    
    %xx=x;yy=y;[x,y]=meshgrid(linspace(min(x),max(x),1000),linspace(min(y),max(y),1000));
    %id=max(1,min(fact*9, ceil(fact*9*(a-amin)/amax)));id(id<=fact*2)=max(id(:));
    %idx=find(isnan(xx));xx(idx)=150;yy(idx)=[yy(idx(1)-1);yy(idx(end)+1)];
    %clf;imagesc(x(1,:),y(:,1),ceil(id/3)); hold on; patch(xx,yy,'k','linewidth',5,'facecolor','w','edgecolor',.75*[1,1,1]); set(gca,'xdir','reverse','ydir','normal'); axis equal; colormap([copper;1,1,1]); set(gcf,'color','w'); axis off;
    if nargout>3
        w=2;
        ab1=accumarray(max(1,min(amax, round((a(oall)-amin)))),b(oall),[amax,1],@min,nan);
        ab2=accumarray(max(1,min(amax, round((a(iall)-amin)))),b(iall),[amax,1],@max,nan);
        for n1=1:w,ab1(2:end-1)=min(min(ab1(2:end-1),ab1(1:end-2)),ab1(3:end));end
        for n1=1:w,ab2(2:end-1)=max(max(ab2(2:end-1),ab2(1:end-2)),ab2(3:end));end
        i=ab1>0&ab2>0;
        af=d*(ab1(i)-ab2(i));
        idx=find(af>0,1); % source position
        % af: area function
        af=min(0,af)+ab_alpha(i).*max(0,af).^ab_beta(i);
        af=af(idx:end);
    end
end
end

% computes vocal tract filter
function [H,f,Hc]=a2h(a,l,n,fs,closure,mina)

if nargin<6 || isempty(mina), mina=min(a,[],1); end
if nargin<5 || isempty(closure), closure=0; end
if nargin<4 || isempty(fs), fs=11025; end
if sum(size(a)>1)<=1,a=a(:);end
if sum(size(l)>1)<=1,l=l(:);end

c=34326;% speed of sound (cm/s)
[NL,ML]=size(l);
[N,M]=size(a);
m=ceil(n/2)+1;
f=fs*(0:ceil(n/2))'/n;
t=l/c;
Rrad=.9*exp(-(abs(f)/4e3).^2);% reflection at lips (low pass)
H=zeros(m,M);%H=zeros(n,M);
Hc=zeros(m,M);
if mina==0, a=max(.05,a); end
%k=0.995;%.999;
for nM=1:M,
    if 1,%mina(nM)>0,
        coswt=cos(2*pi*f*t(:,min(ML,nM))');
        sinwt=1i*sin(2*pi*f*t(:,min(ML,nM))');
        R=[.9;(a(2:N,nM)-a(1:N-1,nM))./max(eps,a(2:N,nM)+a(1:N-1,nM))];
        U=coswt+sinwt;
        V=coswt-sinwt;
        if 1,
            h1=ones(m,1);                                               % signal at glottis
            h2=zeros(m,1);
            for nN=1:N-1,
                RnN=-R(nN);
                u=h1+RnN*h2; v=h2+RnN*h1;    % reflection
                if closure==nN, Hc(:,nM)=u-v; end
                if NL==1, h1=U.*u;h2=V.*v;        % delay
                else      h1=U(:,nN).*u;h2=V(:,nN).*v; end
                %h(:,1)=h(:,1)/k;h(:,2)=h(:,2)*k;
%                 display(sum(h1))
%                 display(sum(h2))
%                 display(RnN)
            end
            
            u=h1-Rrad.*h2; %v=h2-Rrad.*h1;    % reflection
            h=u;              %h(:,2)=v;
            if closure>=N, Hc(:,nM)=u-(h2-Rrad.*h1); end
        elseif 1,
            h=[ ones(m,1), zeros(m,1) ];                                 % signal at glottis
            for nN=1:N-1,
                u=h(:,1)-R(nN).*h(:,2); v=-R(nN).*h(:,1)+h(:,2);    % reflection
                if closure==nN, Hc(:,nM)=u; end
                if NL==1, h(:,1)=U.*u;h(:,2)=V.*v;        % delay
                else      h(:,1)=U(:,nN).*u;h(:,2)=V(:,nN).*v; end
                %h(:,1)=h(:,1)/k;h(:,2)=h(:,2)*k;
            end
            u=h(:,1)-Rrad.*h(:,2); %v=-Rrad.*h(:,1)+h(:,2);    % reflection
            h(:,1)=u;              %h(:,2)=v;
        else
            h=[ ones(m,1), -Rrad ];                                 % signal at lips
            for nN=N:-1:1,
                if NL==1, h(:,1)=U.*h(:,1);h(:,2)=V.*h(:,2);        % delay
                else      h(:,1)=U(:,nN).*h(:,1);h(:,2)=V(:,nN).*h(:,2); end
                %h(:,1)=h(:,1)/k;h(:,2)=h(:,2)*k;
                u=h(:,1)-R(nN).*h(:,2); v=-R(nN).*h(:,1)+h(:,2);    % reflection
                h(:,1)=u;h(:,2)=v;
                if closure==nN, Hc(:,nM)=(1+Rrad).*prod(1+R(nN:N))./h(:,1); end
            end
        end
        H(:,nM)=(1+Rrad).*prod(1+R)./h(:,1);
        if closure>0, Hc(:,nM)=(1+Rrad).*prod(1+R(closure+1:N)).*Hc(:,nM)./h(:,1); end
    end
end
H=cat(1,H,conj(H(1+(n-m:-1:1),:)));
Hc=cat(1,Hc,conj(Hc(1+(n-m:-1:1),:)));
f=cat(1,f,-f(1+(n-m:-1:1)));
if mina==0, H=0*H; end
end

% returns area function
function af=diva_synth_af(Art)
% global DIVA_x;
persistent vt;
if isempty(vt)
    [filepath,filename]=fileparts(mfilename);
    load(fullfile(filepath,[filename,'.mat']),'vt');
end

synth=struct('fs',11025,'update_fs',200); 
synth.f0=120;
synth.samplesperperiod=ceil(synth.fs/synth.f0);
synth.glottalsource=glotlf(0,(0:1/synth.samplesperperiod:1-1/synth.samplesperperiod)');
synth.f=[0,1];
synth.filt=[0,0];
synth.pressure=0;
%synth.modulation=1;
synth.voicing=1;
synth.pressurebuildup=0;
synth.pressure0=0;
synth.sample=zeros(synth.samplesperperiod,1);
synth.k1=1;
synth.numberofperiods=1;
synth.samplesoutput=0;

vt.idx=1:10;
vt.pressure=0;
vt.f0=120;
vt.closed=0;
vt.closure_time=0;
vt.closure_position=0;
vt.opening_time=0;

voices=struct('F0',{120,340},'size',{1,.7});
opt.voices=1;

ndata=size(Art,2);
dt=.005;
time=0;
s=zeros(ceil((ndata+1)*dt*synth.fs),1);
while time<(ndata+1)*dt;
    % sample articulatory parameters
    t0=floor(time/dt);
    t1=(time-t0*dt)/dt;
    [nill,nill,nill,af1,d]=diva_synth_sample(Art(:,min(ndata,1+t0)));
    [nill,nill,nill,af2,d]=diva_synth_sample(Art(:,min(ndata,2+t0)));
    naf1=numel(af1);naf2=numel(af2);
    if naf2<naf1,af2(end+(1:naf1-naf2))=af2(end); end
    if naf1<naf2,af1(end+(1:naf2-naf1))=af1(end); end
    af=af1*(1-t1)+af2*t1;
    FPV=max(-1,min(1, Art(end-2:end,min(ndata,1+t0))*(1-t1)+Art(end-2:end,min(ndata,2+t0))*t1 ));
    vt.voicing=(1+tanh(3*FPV(3)))/2;
    vt.pressure=FPV(2);
    vt.pressure0=vt.pressure>.01;
    vt.f0=100+20*FPV(1);
    
    af0=max(0,af);
    k=.025;af0(af0>0&af0<k)=k;
    minaf=min(af);
    minaf0=min(af0);
    vt.af=af;
%      display(af);
%      DIVA_x.af_sample=DIVA_x.af_sample+1;
%      DIVA_x.af(:,DIVA_x.af_sample)=af;
%    tracks place of articulation
    if minaf0==0, 
        release=0;
        vt.opening_time=0; vt.closure_time=vt.closure_time+1;
        vt.closure_position=find(af0==0,1,'last');
        if ~vt.closed, closure=vt.closure_position; else closure=0; end;
        vt.closed=1;
    else
        if vt.closed, release=vt.closure_position; release_closure_time=vt.closure_time; else release=0; end;
        if (vt.pressure0&&~synth.pressure0) vt.opening_time=0; end;
        vt.opening_time=vt.opening_time+1;
        vt.closure_time=0;
        [nill,vt.closure_position]=min(af);
        closure=0;
        vt.closed=0;
    end
%     display(vt.closed);
    if release>0  af=max(k,af);minaf=max(k,minaf);minaf0=max(k,minaf0); end
    
    if release>0, 
                    vt.f0=(.95+.0*rand)*voices(opt.voices).F0;
                    synth.pressure=0;%modulation=0; 
    elseif  (vt.pressure0&&~synth.pressure0) 
                    vt.f0=(.95+.0*rand)*voices(opt.voices).F0;
                    synth.pressure=vt.pressure; synth.f0=1.25*vt.f0; 
                    synth.pressure=1;%synth.modulation=1; 
    elseif  (~vt.pressure0&&synth.pressure0&&~vt.closed), synth.pressure=synth.pressure/10;
    end
    
    % computes glottal source
    synth.samplesperperiod=ceil(synth.fs/synth.f0);
    pp=[.6,.2-.1*synth.voicing,.25];%10+.15*max(0,min(1,1-vt.opening_time/100))];
    synth.glottalsource=10*.25*glotlf(0,(0:1/synth.samplesperperiod:1-1/synth.samplesperperiod)',pp)+10*.025*synth.k1*glotlf(1,(0:1/synth.samplesperperiod:1-1/synth.samplesperperiod)',pp);
    numberofperiods=synth.numberofperiods;
        
    % computes vocal tract filter
    [synth.filt,synth.f,synth.filt_closure]=a2h(af0,d,synth.samplesperperiod,synth.fs,vt.closure_position,minaf0);
    synth.filt=2*synth.filt/max(eps,synth.filt(1));
    synth.filt(1)=0;
    synth.filt_closure=2*synth.filt_closure/max(eps,synth.filt_closure(1));
    synth.filt_closure(1)=0;
    
    % computes sound signal
    w=linspace(0,1,synth.samplesperperiod)';
    if release>0,%&&synth.pressure>.01,
        u=synth.voicing*1*.010*(synth.pressure+20*synth.pressurebuildup)*synth.glottalsource + (1-synth.voicing)*1*.00*(synth.pressure+20*synth.pressurebuildup)*randn(synth.samplesperperiod,1);
%         if release_closure_time<40
%             u=1*.010*synth.pressure*synth.glottalsource;%.*(0.25+.025*randn(synth.samplesperperiod,1)); % vocal tract filter
%         else
%             u=1*.010*(synth.pressure+synth.pressurebuildup)*randn(synth.samplesperperiod,1);
%         end
        v0=real(ifft(fft(u).*synth.filt_closure));
        numberofperiods=numberofperiods-1;
        synth.pressure=synth.pressure/10;
        vnew=v0(1:synth.samplesperperiod);
        v0=(1-w).*synth.sample(ceil(numel(synth.sample)*(1:synth.samplesperperiod)/synth.samplesperperiod))+w.*vnew;
        synth.sample=vnew;        
    else v0=[]; end
    if numberofperiods>0,
        %u=0.25*synth.modulation*synth.pressure*synth.glottalsource.*(1+.1*randn(synth.samplesperperiod,1)); % vocal tract filter
        u=0.25*synth.pressure*synth.glottalsource.*(1+.0*randn(synth.samplesperperiod,1)); % vocal tract filter
        u=(synth.voicing*u+(1-synth.voicing)*.0*synth.pressure*randn(synth.samplesperperiod,1));
        if minaf0>0&&minaf0<=k, u=minaf/k*u+(1-minaf/k)*.0*synth.pressure*randn(synth.samplesperperiod,1); end
        v=real(ifft(fft(u).*synth.filt));
        
        vnew=v(1:synth.samplesperperiod);
        v=(1-w).*synth.sample(ceil(numel(synth.sample)*(1:synth.samplesperperiod)'/synth.samplesperperiod))+w.*vnew;
        synth.sample=vnew;
        
        if numberofperiods>1
            v=cat(1,v,repmat(vnew,[numberofperiods-1,1]));
        end
    else v=[]; end
    v=cat(1,v0,v);
    v=v+.000*randn(size(v));
    v=(1-exp(-v))./(1+exp(-v));
    s(synth.samplesoutput+(1:numel(v))) =v;
    time=time+numel(v)/synth.fs;
    synth.samplesoutput=synth.samplesoutput+numel(v);
    
    % computes f0/amp/voicing/pressurebuildup modulation
    synth.pressure0=vt.pressure0;
    alpha=min(1,(.1)*synth.numberofperiods);beta=100/synth.numberofperiods;
    synth.pressure=synth.pressure+alpha*(vt.pressure*(max(1,1.5-vt.opening_time/beta))-synth.pressure);
    alpha=min(1,.5*synth.numberofperiods);beta=100/synth.numberofperiods;
    synth.f0=synth.f0 + ... 
                   0*sqrt(alpha)*randn + ...
                       alpha * (vt.f0*max(1,1.25-vt.opening_time/beta)-synth.f0);%147;%120;
    synth.voicing=max(0,min(1, synth.voicing+.5*(vt.voicing-synth.voicing) ));
    %synth.modulation=max(0,min(1, synth.modulation+.1*(2*(vt.pressure>0&&minaf>-k)-1) ));
    alpha=min(1,.1*synth.numberofperiods);
    synth.pressurebuildup=max(0,min(1, synth.pressurebuildup+alpha*(2*(vt.pressure>0&minaf<0)-1) ));
    synth.numberofperiods=max(1,numberofperiods);
end
s=s(1:ceil(synth.fs*ndata*dt));
end

