
# signals can be generated by functions of time
struct SignalFunction{Fn,Fr,El} <: AbstractSignal
    fn::Fn
    first::El
    ω::Fr
    ϕ::Float64
    samplerate::Float64
    offset::Int
end
wrap(x::Number) = (x,)
wrap(x::Tuple) = x
wrap(x) = error("Function must return number or tuple of numbers.")

SignalFunction(fn,ω,ϕ,samplerate) = SignalFunction(fn,ω,ϕ,samplerate,0)
SignalTrait(x::SignalFunction{<:Any,<:Any,T}) where T =
    IsSignal{T}(x.samplerate)
function Base.iterate(x::SignalFunction,i=x.offset) 
    if iszero(i)
        x.first, i+1
    else
        t = i/x.samplerate
        wrap(x.fn(2π*(t*x.ω + x.ϕ)),), i+1
    end
end
function Base.iterate(x::SignalFunction{typeof(sin),Float64},i=x.offset) 
    if iszero(i)
        x.first, i+1
    else
        t = i/x.samplerate
        (sinpi(2*(t*x.ω + x.ϕ)),), i+1
    end
end
function Base.iterate(x::SignalFunction{<:Any,Missing},i=x.offset)
    if iszero(i)
        x.first, i+1
    else
        t = i/x.samplerate
        wrap(x.fn(t + x.ϕ)), i+1
    end
end
function Base.Iterators.drop(x::SignalFunction,n)
    SignalFunction(x.fn,x.first,x.ω,x.ϕ,x.samplerate,x.n+n)
end

Base.Iterators.IteratorEltype(::Type{<:SignalFunction}) = Iterators.HasEltype()
Base.eltype(::Type{<:SignalFunction{Fn,Fr,El}}) where {Fn,Fr,El} = El
Base.Iterators.IteratorSize(::Type{<:SignalFunction}) = Iterators.IsInfinite()

function signal(fn::Function,samplerate;
    ω=missing,frequency=ω,ϕ=0,phase=ϕ)

    SignalFunction(fn,wrap(fn(0)),inHz(ω),
        Float64(inradians(ϕ)/2π),
        Float64(inHz(samplerate)),0)
end

# TODO: handle missing sample rates more broadly for functions
signal(x::typeof(randn),fs=nothing) = SignalFunction(x,(rand(),),nothing,0.0,0.0)
function Base.iterate(x::SignalFunction{typeof(randn)},i=0)
    (randn(),), 0
end

signal(fn::typeof(zero),x) = signal(fn,x,SignalTrait(x))
function signal(fn::typeof(zero),x,::IsSignal)
    signal(zero(signal_eltype(x)),samplerate(x))
end
signal(fn::typeof(zero),x,::Nothing) = error("Value is not a signal: $x")

signal(fn::typeof(one),x) = signal(fn,x,SignalTrait(x))
function signal(fn::typeof(one),x,::IsSignal)
    signal(one(signal_eltype(x)),samplerate(x))
end
signal(fn::typeof(one),x,::Nothing) = error("Value is not a signal: $x")