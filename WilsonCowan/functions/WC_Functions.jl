function adapt_global_coupling(hparams,N::Int64,W::Matrix{Float64},lags::Matrix{Float64},h,t::Float64,u::Vector{Float64})
    @inbounds for ii = 1:N
        @inbounds for jj = 1:N
            if W[jj,ii]  > 0.0
                if lags[jj,ii] == 0.0
                    W[ii,jj]+= 0.001*u[ii]*(u[jj] - h(hparams,t-1.0;idxs=jj))
                else
                    W[ii,jj]+= 0.001*h(hparams,t-lags[jj,ii];idxs=ii)*(u[jj] - h(hparams,t-1.0;idxs=jj))
                end


                if W[ii,jj] > SC[ii,jj] + SC[ii,jj]*0.05
                    W[ii,jj] = SC[ii,jj] + SC[ii,jj]*0.05
                elseif W[ii,jj] < SC[ii,jj] - SC[ii,jj]*0.05
                    W[ii,jj] = SC[ii,jj] - SC[ii,jj]*0.05
                end

        
                #W[ii,jj] += W[ii,jj]*(1-exp(-abs(SC[ii,jj] - W[ii,jj])/0.01))*dw
            end
        end
     

    end

    @inbounds for k=1:N #Maintaining symmetry in the weights between regions
        @inbounds for l = k:N
                 W[k,l] = (W[k,l]+W[l,k])/2
                 W[l,k] = W[k,l]
        end
    end

    W .= W./maximum(W)
    

    return W
  
end


f(x::Float64,β::Float64,θ::Float64) = 1/(1+exp(-β*(x-θ)))

invg(x::Float64,widthp::Real) = -exp((-(x^2))/widthp) + 1
function Θ(x)
    if x > 0.
        return 1.
    else 
        return 0.
    end
end






function adapt_local_func(h,hparams,t,weights,WCp,rE,rI,i,N,c)
    @unpack cEE,cEI,cIE,cII,τE,τI,τx,Pext,θE,θI,β,η,σ = WCp
    @unpack cEEv,cIEv,cEIv,cIIv,cSUM = weights

    hEE = rE - h(hparams,t-1.0;idxs = i)
    hIE = rI - h(hparams,t-1.0;idxs = i+N)
    hEI = rE - h(hparams,t-1.0;idxs = i)
   
    #cEEv[i] = (cEEv[i] + c*rE*(hEE))
    #cIEv[i] = (cIEv[i] + c*rE*(hIE))
    #cIEv[i] = (cIEv[i] + 0.0001*rI*(rE - 0.3))
    #cEIv[i] = (cEIv[i] + c*rI*(hEI))

    #cIIv[i] = (cIIv[i] + c*rI*(rI - h(hparams,t-1.0;idxs = i+N)))
    

    cEEv[i], cIEv[i],cEIv[i], cIIv[i] = cSUM*[cEEv[i], cIEv[i], cEIv[i], cIIv[i]]/(cEEv[i] + cIEv[i] + cEIv[i] + cIIv[i])
 

return cEEv[i],cIEv[i],cEIv[i],cIIv[i]
end



	

function getModelFC(BOLD_TRIALS,nTrials,tstart,step,size)
    modelFC = []
    for i = 1:nTrials
        if i == 1
            modelFC = get_FC(BOLD_TRIALS[:,:,i],tstart,step,size)/nTrials
        else
            modelFC += get_FC(BOLD_TRIALS[:,:,i],tstart,step,size)/nTrials
        end
    end
    return modelFC
end


function remove_zeros_W!(W,N)
    indx_zeros = 0 
    break_cond = false
    for i = 1:size(W,3)
        if W[:,:,i] == zeros(N,N)
            indx_zeros = i
            break
        end
    end

    return W[:,:,1:indx_zeros-1]
end