function [ ampbins,phasebins,sig2powerskew,sig2prefangle,phaseamphist ] = PhaseAmpCouplingByAmp( sig1phase,sig1amp,sig2amp,numbins )
%PhaseAmpCouplingByAmp(sig1phase,sig1amp,sig2amp,ampbins,numphasebins)
%calculates phase amplitude coupling between the phase of signal 1 and the
%amplitude of signal 2, with respect to the amplitude of signal 1.
%
%INPUTS
%   sig1phase
%   sig1amp
%   sig2amp
%   ampbins
%   numbins
%   Detailed explanation goes here
%
%OUTPUTS
%
%
%DLevenstein Fall 2016
%THIS FUNCTION NEEDS DOCUMENTATION, I/O IMPROVEMENTS
%%

phasebins = linspace(-pi,pi,numbins+1);
phasebins=phasebins(1:end-1)+diff(phasebins(1:2));
ampbins = linspace(-2.5,2.5,numbins);
numampbins = length(ampbins);

%%
sig1amp = zscore(sig1amp);
sig2amp = zscore(sig2amp);

%%
sig1binpower = interp1(ampbins,ampbins,sig1amp,'nearest');
sig1binphase = interp1(phasebins,phasebins,sig1phase,'nearest');


powerhist = zeros(numampbins,numbins);
sig2prefangle = zeros(numampbins,1);
sig2powerskew = zeros(numampbins,1);
for bb = 1:numampbins
    for bbb = 1:numbins
        ampwintimes = sig2amp(sig1binpower==ampbins(bb) & sig1binphase==phasebins(bbb));
        phaseamphist(bb,bbb) = mean(ampwintimes);
    end
    sig2powerskew(bb) = mean(sig2amp(sig1binpower==ampbins(bb)).*exp(1i.*sig1phase(sig1binpower==ampbins(bb))));
    sig2prefangle(bb) = angle(sig2powerskew(bb));
    sig2powerskew(bb) = abs(sig2powerskew(bb));
end

%% Figure

rwbcolormap = makeColorMap([0 0 0.8],[1 1 1],[0.8 0 0]);
plotx = linspace(-pi,3*pi,100);
figure
subplot(2,2,1)
hold on
    imagesc(phasebins,ampbins,phaseamphist)
    imagesc(phasebins+2*pi,ampbins,phaseamphist)
    plot(sig2prefangle,ampbins,'.k')
    plot(sig2prefangle+2*pi,ampbins,'.k')
    plot(plotx,cos(plotx),'k')
    colormap(gca,rwbcolormap)
    axis xy
    axis tight
    ColorbarWithAxis([-0.5 0.5],['Mean Amp.'])
    caxis([-0.5 0.5])
  %  xlim([-pi 3*pi]);ylim(ampbins([1 end]))
    xlabel('Signal 1 Phase');ylabel('Signal 1 Amp (Z)')
subplot(4,2,2)
    plot(ampbins,sig2powerskew,'k','LineWidth',1)
    xlabel('Signal 1 Amp. (Z)');
    ylabel('Phase-Amp. Modulation (mrl)')
    axis tight
subplot(4,2,6)
    histogram(sig1amp,ampbins)
    xlabel('Signal 1 Amp. (Z)');
    ylabel('Occupancy')
    axis tight
    title('Signal1/2 Amp. Distributions')
    
subplot(4,2,8)
    histogram(sig2amp,ampbins)
    xlabel('Signal 2 Amp. (Z)');
    ylabel('Occupancy')
    axis tight 

subplot(2,2,3)
    gasphist = hist3([sig1amp,sig2amp],{ampbins,ampbins});
    imagesc(ampbins,ampbins,gasphist)
    axis xy
    xlabel('Signal 1 Power');ylabel('Signal 2 Power')
end

