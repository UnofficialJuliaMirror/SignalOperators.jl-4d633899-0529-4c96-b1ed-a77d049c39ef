using Random

# helpers
astuple(x::Number) = (x,)
astuple(x::Tuple) = x
astuple(x) = error("Function must return number or tuple of numbers.")
ntuple_T(::Type{<:NTuple{<:Any,T}}) where T = T
ntuple_N(::Type{<:NTuple{N}}) where N = N

# signals can be generated by functions of time
struct SignalFunction{Fn,Fr,El,T,Fs} <: AbstractSignal{El}
    fn::Fn
    first::El
    ω::Fr
    ϕ::Float64
    samplerate::Fs
    function SignalFunction(fn::Fn,first::El,ω::Fr,ϕ,
        sr::Fs=missing) where {Fn,El,Fr,Fs}

        new{Fn,Fr,El,ntuple_T(El),Fs}(fn,first,ω,ϕ,sr)
    end
end
SignalTrait(::Type{<:SignalFunction{<:Any,<:Any,El,Fs}}) where {El,Fs} = 
    IsSignal{ntuple_T(El),Fs,Nothing}()
nchannels(x::SignalFunction) = ntuple_N(typeof(x.first))
nsamples(x::SignalFunction) = nothing
samplerate(x::SignalFunction) = x.samplerate

@Base.propagate_inbounds function sampleat!(result,
    x::SignalFunction{Fn,Fr},::IsSignal,i,j) where {Fn,Fr}

    t = j/x.samplerate
    if Fn <: typeof(sin)
        if Fr <: Missing
            writesink(result,i,@. sinpi(2*(t+x.ϕ)))
        else
            writesink(result,i,@. sinpi(2*(t*x.ω + x.ϕ)))
        end
    else
        if Fr <: Missing
            writesink(result,i,@. x.fn(t + x.ϕ))
        else
            writesink(result,i,@. x.fn(2π*(t*x.ω + x.ϕ)))
        end
    end
end
tosamplerate(x::SignalFunction,::IsSignal,::ComputedSignal,fs) = 
    SignalFunction(x.fn,x.first,x.ω,x.ϕ,coalesce(fs,x.fs))

function signal(fn::Function,samplerate::Union{Missing,Number}=missing;
    ω=missing,frequency=ω,ϕ=0,phase=ϕ)

    SignalFunction(fn,astuple(fn(0)),inHz(ω),
        inradians(Float64,ϕ)/2π,
        inHz(Float64,samplerate))
end

struct RandFn{R}
    rng::R
end
signal(x::typeof(randn),fs::Union{Missing,Number}=missing;rng=Random.GLOBAL_RNG) =
    SignalFunction(RandFn(rng),(randn(rng),),missing,0.0,inHz(Float64,fs))
@Base.propagate_inbounds function sampleat!(result,
    x::SignalFunction{<:RandFn},::IsSignal,i,j)

    writesink(result,i,randn(x.fn.rng))
end