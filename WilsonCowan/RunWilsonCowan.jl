include("functions/WC_InitSetup.jl")
ntrial = 1
R_trial = zeros(140,140,ntrial)
for jj = 1:ntrial

    global WCS = setup()

    WCS.modelparams.ρ=0.15
    WCS.solver_opts.tWindow = 2000
    WCS.stimparams.stim = false
    WCS.stimparams.Tstim = [750,800]
    WCS.stimparams.stimStr = -1.

    @time wilsoncowan_windows_run()


    N=140
    BOLD = WCS.bold_out[:,800:end]
    R = zeros(N,N)
    for i = 1:N
        for j = i+1:N
            R[i,j] = cor(BOLD[i,:],BOLD[j,:])
            R[j,i] = R[i,j]
        end
    end

    R_trial[:,:,jj] = R
end

heatmap(R_trial[:,:,1])