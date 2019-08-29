using Statistics

zero_helper(sig) = zero_helper(eltype(sig),length(sig))
zero_helper(x::Type{<:NTuple{M,T}},N) where {M,T} = zeros(T,N,M)


"""
    sink(signal)

Converts any signal into an array, with time as the rows and channels as
the columns.
"""
sink(x::AbstractArray) = x
sink(x) = sink(x,SignalTrait(x))
function sink(x,s::IsSignal) 
    smp = samples(x)
    sink(x,s,smp,Iterators.IteratorSize(x))
end
sink(x, ::Nothing) = error("Don't know how to interpret value as an array: $x")
function sink(xs,::IsSignal,smp,::Iterators.HasLength)
    result = zero_helper(smp)
    samples_to_result!(result,smp)
end
function samples_to_result!(result,smp)
    for (i,x) in enumerate(smp)
        result[i,:] .= x
    end
    result
end
function sink(x,::IsSignal,smp,::Iterators.IsInfinite)
    error("Cannot store infinite signal in an array. (Use `until`?)")
end

abstract type WrappedSignal{T} <: AbstractSignal
end

"""
    childsignal(x)

Retrieve the signal wrapped by x of type `WrappedSignal`
"""
function childsignal
end
samplerate(x::WrappedSignal) = samplerate(childsignal(x))
SignalTrait(x::WrappedSignal) = SignalTrait(childsignal(x))

Base.Iterators.IteratorEltype(::Type{<:WrappedSignal}) = Iterators.HasEltype()
Base.eltype(::Type{<:WrappedSignal{T}}) where T = eltype(T)
Base.Iterators.IteratorSize(::Type{<:WrappedSignal{T}}) where T =
    Iterators.IteratorSize(T)
Base.length(x::WrappedSignal) = length(childsignal(x))