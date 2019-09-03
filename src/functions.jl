
# helpers
astuple(x::Number) = (x,)
astuple(x::Tuple) = x
astuple(x) = error("Function must return number or tuple of numbers.")
ntuple_T(x::NTuple{<:Any,T}) where T = T
ntuple_N(x::NTuple{N}) where N = T

# signals can be generated by functions of time
struct SignalFunction{Fn,Fr,El,Fs} <: AbstractSignal
    fn::Fn
    first::El
    ω::Fr
    ϕ::Float64
    samplerate::Fs
    offset::Int
    SignalFunction(fn,ω,ϕ,sr=missing,offset=0) = new(fn,ω,ϕ,sr,offset)
end
SignalTrait(::Type{<:SignalFunction{<:Any,<:Any,El,Fs}}) where {El,Fs} = 
    IsSignal{ntuple_T(El),Fs,Nothing}()
nchannels(x::SignalFunction) = ntuple_N(typeof(x.first))
nsamples(x::SignalFunction) = nothing
samplerate(x::SignalFunction) = x.samplerate

function Base.iterate(x::SignalFunction{Fn,Fr},i=x.offset) where {Fn,Fr}
    if iszero(i)
        x.first, i+1
    elseif Fn <: typeof(sin) && 
        t = i/x.samplerate
        if Fr <: Missing
            (sinpi(2*(t+x.ϕ)),), i+1
        else
            (sinpi(2*(t*x.ω + x.ϕ)),), i+1
        end
    else
        t = i/x.samplerate
        if Fr <: Missing
            astuple(x.fn(t + x.ϕ)), i+1
        else
            astuple(x.fn(2π*(t*x.ω + x.ϕ)),), i+1
        end
    end
end
function Base.Iterators.drop(x::SignalFunction,n::Int)
    SignalFunction(x.fn,x.first,x.ω,x.ϕ,x.samplerate,x.offset+n)
end
tosamplerate(x::SignalFunction,::IsSignal,::ComputedSignal,fs=missing) = 
    SignalFunction(x.fn,x.first,x.ω,x.ϕ,coalesce(fs,x.fs),x.offset)

function signal(fn::Function,samplerate=missing;
    ω=missing,frequency=ω,ϕ=0,phase=ϕ)

    SignalFunction(fn,astuple(fn(0)),inHz(ω),
        inradians(Float64,ϕ)/2π,
        inHz(Float64,samplerate),0)
end

signal(x::typeof(randn),fs=missing;rng=Random.GLOBAL_RNG) =
    SignalFunction(x,(randn(rng),),fs,missing,inHz(Float64,fs),0)
Base.iterate(x::SignalFunction{typeof(randn)},i=blank) = (randn(),), blank