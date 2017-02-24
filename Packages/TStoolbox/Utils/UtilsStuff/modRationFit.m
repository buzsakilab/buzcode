function [modRatio, prefPh, fitCurve, H]= modRatioFit(phase)


[H,X] = hist(phase,100);
[estimates, model] = fitPhaseHist(2*pi*X,H);

[sse, fitCurve] = model(estimates);
prefPh = estimates(3);
modRatio = estimates(2)/estimates(1);

end


function [estimates, model] = fitPhaseHist(xdata, ydata)

start_point = rand(1, 3);
model = @fitFct;
estimates = fminsearch(model, start_point);

    function [sse, FittedCurve] = fitFct(params)
        A = params(1);
	B = params(2);
        phi0 = params(3);
        FittedCurve = A + B.*cos(xdata-phi0);
        ErrorVector = FittedCurve - ydata;
        sse = sum(ErrorVector .^ 2);
    end

end