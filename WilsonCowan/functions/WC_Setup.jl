function setup()
   #load data and make struct & dist matrices
   
   println("setting up Wilson Cowan Model")

   c = 13000
   constant_delay = 0.005

   SC,dist,lags,N,FC,missingROIs = networksetup(c,constant_delay)
   
   wmat = zeros(N,N)
   d=zeros(N)

   network_params = networkParameters(SC,dist,round.(lags,digits=5),N)



   WCS = WCstruct(WCparams(Pext = 0.31,η=0.089,σ=0.001),
                  network_params,
                  balloonModelParameters(),
                  zeros(4N),
                  StimParams(),
                  SolverOptions(),
                  Helpers(wmat,d,1),
                  zeros(N)
   )
  

   return WCS
end


