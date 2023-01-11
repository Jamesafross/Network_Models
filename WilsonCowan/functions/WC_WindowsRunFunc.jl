function wilsoncowan_windows_run()
    @unpack modelparams,balloonparams,network,stimparams,helpers,solver_opts = WCS
 
    @unpack cEE,cEI,cIE,cII,τE,τI,τx,Pext,θE,θI,β,η,σ,τISP,ρ = modelparams
    @unpack W,lags,N = network

    @unpack stimWindow,stimNodes,stimStr,Tstim = stimparams
    @unpack wmat,d,current_window = helpers

    @unpack tWindow, nWindow = solver_opts
    BOLD_saveat = collect(0:1.0:tWindow)
    size_out = length(BOLD_saveat)
    BOLD_out = zeros(N,size_out,nWindow)
    
    for j = 1:nWindow
    
        WCS.helpers.current_window = j
       
        println("working on window . . . ",j)
        if j == 1
            WCS.IC = rand(5N)
            WCS.IC[3N+1:4N] .= WCS.modelparams.cIE
            WCS.IC[4N+1:5N] .= WCS.modelparams.cEE
            hparams = WCS.IC
            h = h1 
        else
            WCS.IC = sol[:,end]
            iStart = findfirst(sol.t .> tWindow - 1.2)
            u_hist = make_uhist(sol.t[iStarst:end] .- sol.t[end],sol[:,iStart:end])
            hparams = u_hist
            h = h2
        end

        tspan = (0.0,tWindow)
        
        p = hparams
        
        if WCS.solver_opts.use_noise_process == true
            prob = SDDEProblem(WC, dW, WCS.IC, h, tspan, p,noise=WCS.solver_opts.noise_process)
        else
            prob = SDDEProblem(WC, dW, WCS.IC, h, tspan, p)
        end
       
        println("solving...")
        global sol = solve(prob,RKMil(),maxiters = 1e20,save_noise=true)

        save("$WORKDIR/_helpers/noise_processes/noise_process$(j).jld","noise_process$(j)",sol.W)
       
        BalloonIn= make_In(sol.t,sol[1:N,:])
        tspanB = (sol.t[1],sol.t[end])
        balloonParams = balloonparams,BalloonIn,N
        if j == 1
            b0 =  cat(zeros(N),ones(3N),dims=1)
        else
            b0 = endBM
        end
        println(". . . Running Balloon Model")

        out,endBM = runBalloon(b0,balloonParams,tspanB,BOLD_saveat,N)
        
        if j == 1
            WCS.bold_out = out
        else
            WCS.bold_out = cat(WCS.BOLD_out,out,dims=2)
        end

    end

end


