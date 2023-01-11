@with_kw mutable struct WCparams{R}
    cEE::R = 3.5
    cEI::R  = 3.75
    cIE::R  = 2.5
    cII::R  = 0.
    θE::R  = 1.
    θI::R  = 1.
    β::R = 4.
    η::R  = 0.12
    σ::R  = 0.005
    τE::R  = 0.01
    τI::R  =0.02
    τx::R = .01
    Pext::R  =0.315
    Qext::R = 0.0
    τISP::R=0.25
    ρ::R=0.30
end


mutable struct networkParameters
    W::Matrix{Float64}
    dist::Matrix{Float64}
    lags::Matrix{Float64}
    N::Int64
end

mutable struct Helpers
    wmat::Array{Float64}
    d::Array{Float64}
    current_window::Real
end

@with_kw mutable struct PlasticityOptions
    plast_start_win::Real = 1 
    plast_time::Real = 0.01
end

@with_kw mutable struct weights
    cEE::Real = 3.5
    cEI::Real = 3.75
    cIE::Real = 2.5
    cII::Real = 0.
end



@with_kw mutable struct StimParams
    stim::Bool = false
    stimWindow::Real = 1
    stimNodes::Vector{Real} = [39]
    stimStr::Real = 0.
    Tstim::Vector{Real} = [30,60]
end

@with_kw mutable struct SolverOptions
    tWindow::Real = 6000.
    nWindow::Real = 1
    noise_process = WienerProcess(0.0,0.0,0.0)
    use_noise_process::Bool = false
end


@with_kw mutable struct weightsSave
    cEEv
    cIEv
    cEIv
    cIIv
    Wv
    count
end

mutable struct WCstruct
    modelparams::WCparams
    network::networkParameters
    balloonparams::balloonModelParameters
    IC::Vector{Float64}
    stimparams::StimParams
    solver_opts::SolverOptions
    helpers::Helpers
    bold_out
end




	
