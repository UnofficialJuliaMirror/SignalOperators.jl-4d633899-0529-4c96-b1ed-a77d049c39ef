export append, prepend, pad, cycle, mirror, lastsample

################################################################################
# appending signals

struct AppendSignals{Si,Sis,T,L} <: WrappedSignal{Si,T}
    signals::Sis
    len::L
end
SignalTrait(x::Type{T}) where {Si,T <: AppendSignals{Si}} =
    SignalTrait(x,SignalTrait(Si))
function SignalTrait(x::Type{<:AppendSignals{Si,Rst,T,L}},
        ::IsSignal{T,Fs}) where {Si,Rst,T,L,Fs}
    IsSignal{T,Fs,L}()
end
child(x::AppendSignals) = x.signals[1]
nsamples(x::AppendSignals) = x.len
duration(x::AppendSignals) = sum(duration.(x.signals))

"""
    append(x,y,...)

Append a series of signals, one after the other.
"""
append(y) = x -> append(x,y)
"""
    prepend(x,y,...)

Prepend the series of signals: `prepend(xs...)` is equivalent to
`append(reverse(xs)...)`.

"""
prepend(x) = y -> append(x,y)
prepend(x,y,rest...) = prepend(reverse((x,y,rest...)...))

function append(xs...)
    if any(isinf ∘ nsamples,xs[1:end-1])
        error("Cannot append to the end of an infinite signal")
    end
    xs = uniform(xs,channels=true)

    El = promote_type(channel_eltype.(xs)...)
    xs = map(xs) do x
        if channel_eltype(x) != El
            toeltype(x,El)
        else
            x
        end
    end

    len = sum(nsamples,xs)
    AppendSignals{typeof(xs[1]),typeof(xs),El,typeof(len)}(xs, len)
end
tosamplerate(x::AppendSignals,s::IsSignal{<:Any,<:Number},c::ComputedSignal,fs;blocksize) =
    append(tosamplerate.(x.signals,fs;blocksize=blocksize)...)
tosamplerate(x::AppendSignals,s::IsSignal{<:Any,Missing},__ignore__,fs;
    blocksize) = append(tosamplerate.(x.signals,fs;blocksize=blocksize)...)

struct AppendChunk{S,C} <: AbstractChunk
    signal::S
    child::C
    init::Ref{Any}
    k::Int
end
nsamples(x::AppendChunk) = nsamples(x.child)
@Base.propagate_inbounds sample(x::AppendChunk,i) = sample(x.child,i)

function initchunk(x::AppendSignals)
    k = 1
    chunk = initchunk(x.signals[k])
    K = length(x.signals)
    while k < K && isnothing(chunk)
        chunk = initchunk(x.signals[k])
    end
    AppendChunk(x.signals[k],chunk,Ref(nothing),k)
end

function nextchunklen(x::AppendSignals,chunk::AppendChunk,maxlen,skip)
    curlen = nextchunklen(chunk.signal,chunk.child)

    if curlen > 0
        chunk.init[] = chunk.k, chunk.child
        curlen
    else
        k = chunk.k
        K = length(x.signals)
        nextchild = chunk.child
        while k < K && curlen == 0
            k += 1
            nextchild = initchunk(x.signals[k])
            curlen = nextchunklen(x.signals[k],nextchild,maxlen,skip)
        end
        chunk.init[] = (k,nextchild)
        curlen
    end
end

function nextchunk(x::AppendSignals,chunk::AppendChunk,maxlen,skip)
    len = nextchunklen(x,chunk,maxlen,skip)

    k,childchunk = chunk.init[]
    childchunk = nextchunk(x.signals[k],childchunk,maxlen,skip)

    if !isnothing(childchunk)
        AppendChunk(x.signals[k],childchunk,Ref{Any}(nothing),k)
    end
end

Base.show(io::IO,::MIME"text/plain",x::AppendSignals) = pprint(io,x)
function PrettyPrinting.tile(x::AppendSignals)
    if length(x.signals) == 2
        child = signaltile(x.signals[1])
        operate = literal("append(") * signaltile(x.signals[2]) * literal(")") |
            literal("append(") / indent(4) * signaltile(x.signals[2]) / literal(")")
        tilepipe(child,operate)
    else
        list_layout(map(signaltile,x.signals),prefix="append",sep=",",sep_brk=",")
    end
end
signaltile(x::AppendSignals) = PrettyPrinting.tile(x)

################################################################################
# padding
struct PaddedSignal{S,T} <: WrappedSignal{S,T}
    signal::S
    pad::T
end
SignalTrait(x::Type{T}) where {S,T <: PaddedSignal{S}} =
    SignalTrait(x,SignalTrait(S))
SignalTrait(x::Type{<:PaddedSignal},::IsSignal{T,Fs}) where {T,Fs} =
    IsSignal{T,Fs,InfiniteLength}()
nsamples(x::PaddedSignal) = inflen
duration(x::PaddedSignal) = inflen
tosamplerate(x::PaddedSignal,s::IsSignal{<:Any,<:Number},c::ComputedSignal,fs;blocksize) =
    PaddedSignal(tosamplerate(x.signal,fs,blocksize=blocksize),x.pad)
tosamplerate(x::PaddedSignal,s::IsSignal{<:Any,Missing},__ignore__,fs;
    blocksize) = PaddedSignal(tosamplerate(x.signal,fs;blocksize=blocksize),x.pad)

"""

    pad(x,padding)

Create a signal that appends an infinite number of values, `padding`, to `x`.
The value `padding` can be:

- a number
- a tuple or vector
- a type function: a one argument function of the `channel_eltype` of `x`
- a value function: a one argument function of the signal `x` for which
    `SignalOperators.valuefunction(padding) == true`.
- an indexing function: a three argument function following the same type
  signature as `getindex` for two dimensional arrays.

If the signal is already infinitely long (e.g. a previoulsy padded signal),
`pad` has no effect.

If `padding` is a number it is used as the value for all samples and channels
past the end of `x`.

If `padding` is a tuple or vector it is the value for all samples past the end
of `x`.

If `padding` is a type function it is passed the [`channel_eltype`](@ref) of
the signal and the resulting value is used as the value for all samples past
the end of `x`. Examples include `zero` and `one`

If `padding` is a value function it is passed `x` just before padding during
`sink` begins and it should return a tuple of `channel_eltype(x)` values.
This value is repeated for the remaining samples. It is generally only useful
when x is an AbstractArray.

If `padding` is an indexing function (it accepts 3 arguments) it will be used
to retrieve samples from the signal `x` assuming it conforms to the
`AbstractArray` interface, with the first index being samples and the second
channels. If the sample index goes past the bounds of the array, it should be
transformed to an index within the range of that array. Note that such
padding functions only work on signals that are also AbstractArray objects.
You can always generate an array from a given signal by first passing it
through `sink` or `sink!`.

## See also

[`cycle`](@ref)
[`mirror`](@ref)
[`lastsample`](@ref)
[`valuefunction`](@ref)
"""
pad(p) = x -> pad(x,p)
function pad(x,p)
    x = signal(x)
    isinf(nsamples(x)) ? x : PaddedSignal(x,p)
end

"""
    lastsample

When passed as an argument to `pad`, allows padding using the last sample of a
signal. You cannot use this function in other contexts, and it will normally
throw an error. See [`pad`](@ref).
"""
lastsample(x) = error("Must be passed as argument to `pad`.")

"""
    SignalOperators.valuefunction(fn)

Returns true if `fn` should be treated as a value function. See
[`pad`](@ref). If you wish your own function to be a value function, you can
do this as follows.

    SignalOperators.valuefunction(::typeof(myfun)) = true

"""
valuefunction(x) = false
valuefunction(::typeof(lastsample)) = true

"""
    cycle(x,i,j)

An indexing function which wraps index i using mod, thus
repeating the signal when i > size(x,1). It can be passed as the second
argument to [`pad`](@ref).
"""
@Base.propagate_inbounds cycle(x,i,j) = x[(i-1)%end+1,j]

"""
    mirror(x,i,j)

An indexing function which mirrors the indices when i > size(x,1). This means
that past the end of the signal x, the signal first repeats with samples in
reverse order, then repeats in the original order, so on and so forth. It
can be passed as the second argument to  [`pad`](@ref).
"""
@Base.propagate_inbounds function mirror(x,i,j)
    function helper(i,N)
       count,remainder = divrem(i-1,N)
       iseven(count) ? remainder+1 : N-remainder
    end
    x[helper(i,end),j]
end

usepad(x::PaddedSignal,chunk) = usepad(x,SignalTrait(x),chunk)
usepad(x::PaddedSignal,s::IsSignal,chunk) = usepad(x,s,x.pad,chunk)
usepad(x::PaddedSignal,s::IsSignal{T},p::Number,chunk) where T =
    Fill(convert(T,p),nchannels(x.signal))
function usepad(x::PaddedSignal,s::IsSignal{T},
    p::Union{Array,Tuple},chunk) where T

    map(x -> convert(T,x),p)
end
usepad(x::PaddedSignal,s::IsSignal,::typeof(lastsample),chunk) =
    sample(chunk,nsamples(chunk))
function usepad(x::PaddedSignal,s::IsSignal{T},fn::Function,chunk) where T
    nargs = map(x -> x.nargs - 1, methods(fn).ms)
    if 3 ∈ nargs
        if indexable(x.signal)
            i -> fn(x.signal,i,:)
        else
            io = IOBuffer()
            show(io,MIME("text/plain"),x)
            sig_string = String(take!(io))
            error("Attemped to specify an indexing pad function for the ",
                  "following signal, which is not known to support ",
                  "`getindex`.\n",sig_string)
        end
    elseif 1 ∈ nargs
        if valuefunction(fn)
            fn(x.signal)
        else
            Fill(fn(T),nchannels(x.signal))
        end
    else
        error("Pad function ($fn) must take 1 or 3 arguments. ",
              "Refering `pad` help.")
    end
end

child(x::PaddedSignal) = x.signal

struct UsePad
end
const use_pad = UsePad()

struct PadChunk{P,C} <: AbstractChunk
    pad::P
    child_or_len::C
    offset::Int
end
child(x::PadChunk{<:Nothing}) = x.child_or_len
child(x::PadChunk) = nothing
nsamples(x::PadChunk{<:Nothing}) = nsamples(child(x))
nsamples(x::PadChunk) = x.child_or_len

sample(x::PadChunk{<:Nothing},i) = sample(child(x),i)
sample(x::PadChunk{<:Function},i) = x.pad(i + x.offset)
sample(x::PadChunk,i) = x.pad

function initchunk(x::PaddedSignal)
    chunk = initchunk(child(x))
    if nextchunklen(child(x),chunk,inflen,false) == 0
        PadChunk(usepad(x,chunk),0)
    else
        PadChunk(nothing,chunk,0)
    end
end

function nextchunklen(x::PaddedSignal,chunk::PadChunk{<:Nothing},maxlen,skip)
    clen = nextchunklen(child(x),child(chunk),maxlen,skip)
    clen == 0 ? inflen : clen
end
nextchunklen(x::PaddedSignal,::PadChunk,maxlen,skip) = min(maxlen,inflen)
function nextchunk(x::PaddedSignal,chunk::PadChunk{<:Nothing},maxlen,skip)
    len = nextchunklen(x,chunk,maxlen,skip)
    childchunk = nextchunk(child(x),child(chunk),len,skip)
    if isnothing(childchunk)
        PadChunk(usepad(x,chunk),len,nsamples(chunk) + chunk.offset)
    else
        PadChunk(nothing,childchunk,nsamples(chunk) + chunk.offset)
    end
end
function nextchunk(x::PaddedSignal,chunk::PadChunk,maxlen,skip)
    len = nextchunklen(x,chunk,maxlen,skip)
    PadChunk(chunk.pad,len,nsamples(chunk) + chunk.offset)
end

Base.show(io::IO,::MIME"text/plain",x::PaddedSignal) = pprint(io,x)
function PrettyPrinting.tile(x::PaddedSignal)
    child = signaltile(x.signal)
    operate = literal(string("pad(",x.pad,")"))
    tilepipe(child,operate)
end
signaltile(x::PaddedSignal) = PrettyPrinting.tile(x)
