function WC(du,u,h,p,t)
    hparams = p

    @unpack modelparams,network,stimparams,helpers,solver_opts = WCS
 
    @unpack cEE,cEI,cIE,cII,τE,τI,τx,Pext,Qext,θE,θI,β,η,σ,τISP,ρ = modelparams
    @unpack W,lags,N = network

    @unpack stimWindow,stimNodes,stimStr,Tstim = stimparams
    @unpack wmat,d,current_window = helpers


    make_hist_mat2_threads!(h,W,u,hparams,N,lags,t,wmat)
    
    sum!(d,wmat)

    @inbounds Threads.@threads for i = 1:N

        E = u[i]
        I = u[i+N]
        s =  stim(t,i,stimparams,current_window) 
        du[i] = (1/τE)*(-E + f(u[i+4N]*E -  u[i+3N]*I + s+ u[i+2N]+ Pext + (η)*d[i],β,θE))
        du[i+N] =(1/τI)*( -I + f(cEI*E - cII*I+u[i+2N] + Qext,β,θI) )
        du[i+2N] = (-1/τx)*u[i+2N]
        du[i+3N] = 0.01*(I*(E-ρ))
        du[i+4N] = 0.0001*(E*(E - h(hparams,t-1.;idxs=i)))
    end
end

function dW(du,u,h,p,t)
    @unpack modelparams,network= WCS
    @unpack σ = modelparams
    @unpack N = network

    for i = 1:N
        du[i+2N] = σ
    end
end