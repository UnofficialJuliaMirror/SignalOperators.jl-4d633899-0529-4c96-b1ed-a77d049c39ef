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
SignalTrait(::Type{<:SignalFunction{<:Any,<:Any,<:Any,T,Fs}}) where {T,Fs} = 
    IsSignal{T,Fs,InfiniteLength}()
nchannels(x::SignalFunction) = ntuple_N(typeof(x.first))
nsamples(x::SignalFunction) = inflen
samplerate(x::SignalFunction) = x.samplerate
EvalTrait(x::SignalFunction) = ComputedSignal()

function Base.show(io::IO, ::MIME"text/plain",x::SignalFunction) 
    if ismissing(x.ω) && iszero(x.ϕ)
        write(io,string(x.fn))
        show_fs(io,x)
    else
        write(io,"signal(")
        write(io,string(x.fn))
        write(io,",")
        !ismissing(x.ω) && write(io,"ω=",string(x.ω))
        !iszero(x.ϕ) && write(io,"ϕ=",string(x.ϕ))
        write(io,")")
        show_fs(io,x)
    end
end

@Base.propagate_inbounds sampleat!(result,x::SignalFunction,i,j,check) =
    writesink!(result,i,x.fn(2π*((t/x.samplerate*x.ω + x.ϕ) % 1)))

@Base.propagate_inbounds sampleat!(result,
    x::SignalFunction{<:Any,Missing},i,j,check) =
    writesink!(result,i,x.fn(j/x.samplerate + x.ϕ))

@Base.propagate_inbounds sampleat!(result,
    x::SignalFunction{typeof(sin)},i,j,check) =
    writesink!(result,i,sinpi(2*(j/x.samplerate*x.ω + x.ϕ)))

@Base.propagate_inbounds sampleat!(result,
    x::SignalFunction{typeof(sin),Missing},i,j,check) =
    writesink!(result,i,sinpi(2*(j/x.samplerate + x.ϕ)))

tosamplerate(x::SignalFunction,::IsSignal,::ComputedSignal,fs;blocksize) = 
    SignalFunction(x.fn,x.first,x.ω,x.ϕ,coalesce(inHz(Float64,fs),x.samplerate))

abstract type Functor
end

"""
## Functions

    signal(fn,[samplerate];[ω/frequency],[ϕ/phase])

Functions can define infinite length signals of known or unknown sample rate.
The function `fn` can either return a number, or, for multi-channel signals,
a tuple of values. 

The input to `fn` is either a phase value or a time value. If passed a
frequency (using either the ω or frequency keyword), the input to `fn` will
be a phase value in radians, ranging from 0 to 2π. If no frequency is
specified the value passed to `fn` is the time in seconds. Specifying phase
(by the ϕ or phase keyword) will first add that value to the input before
passing it to `fn`. The phase is assumed to be in units of radians (but you
can also pass degrees by using `°` or a unit of time (e.g. `s` for seconds)).

"""
function signal(fn::Union{Function,Functor},
    samplerate::Union{Missing,Number}=missing;
    ω=missing,frequency=ω,ϕ=0,phase=ϕ)

    SignalFunction(fn,astuple(fn(0)),inHz(ω),
        inradians(Float64,ϕ,ω)/2π,
        inHz(Float64,samplerate))
end

struct RandFn{R}
    rng::R
end
"""

If `fn == randn` no frequency or phase can be specified. Instead there is a
single keyword argument, `rng`, which allows you to specify the random number
generator; `rng` defaults to `Random.GLOBAL_RNG`.

"""
signal(x::typeof(randn),fs::Union{Missing,Number}=missing;rng=Random.GLOBAL_RNG) =
    SignalFunction(RandFn(rng),(randn(rng),),missing,0.0,inHz(Float64,fs))
@Base.propagate_inbounds function sampleat!(result,
    x::SignalFunction{<:RandFn},i,j,check)

    writesink!(result,i,randn(x.fn.rng))
end