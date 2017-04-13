N = 100;
parfor_progress(N);
parfor i=1:N
    pause(rand); % Replace with real code
    parfor_progress;
end
parfor_progress(0);