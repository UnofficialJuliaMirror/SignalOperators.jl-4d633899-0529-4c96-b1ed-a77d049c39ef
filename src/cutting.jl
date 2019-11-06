export until, after

################################################################################
# cutting signals

struct CutApply{Si,Tm,K,T} <: WrappedSignal{Si,T}
    signal::Si
    time::Tm
end
CutApply(signal::T,time,fn) where T = CutApply(signal,SignalTrait(T),time,fn)
CutApply(signal::Si,::IsSignal{T},time::Tm,kind::K) where {Si,Tm,K,T} =
    CutApply{Si,Tm,K,T}(signal,time)

SignalTrait(::Type{T}) where {Si,T <: CutApply{Si}} =
    SignalTrait(T,SignalTrait(Si))
function SignalTrait(::Type{<:CutApply{Si,Tm,K}},::IsSignal{T,Fs,L}) where
    {Si,Tm,K,T,Fs,L}

    if Fs <: Missing
        IsSignal{T,Missing,Missing}()
    elseif K <: Val{:until}
        IsSignal{T,Float64,Int}()
    elseif K <: Val{:after}
        IsSignal{T,Float64,L}()
    else
        error("Unexpected cut apply type $K")
    end
end

child(x::CutApply) = x.signal
resolvelen(x::CutApply) = insamples(Int,maybeseconds(x.time),samplerate(x))

const UntilApply{S,T} = CutApply{S,T,Val{:until}}
const AfterApply{S,T} = CutApply{S,T,Val{:after}}

"""
    until(x,time)

Create a signal of all samples of `x` up until and including `time`.
"""
until(time) = x -> until(x,time)
until(x,time) = CutApply(signal(x),time,Val{:until}())

"""
    after(x,time)

Create a signal of all samples of `x` after `time`.
"""
after(time) = x -> after(x,time)
after(x,time) = CutApply(signal(x),time,Val{:after}())

Base.show(io::IO,::MIME"text/plain",x::CutApply) = pprint(io,x)
function PrettyPrinting.tile(x::CutApply)
    operate = literal(string(cutname(x),"(",(x.time),")"))
    tilepipe(signaltile(x.signal),operate)
end
signaltile(x::CutApply) = PrettyPrinting.tile(x)

cutname(x::UntilApply) = "until"
cutname(x::AfterApply) = "after"

nsamples(x::UntilApply) = min(nsamples(x.signal),resolvelen(x))
duration(x::UntilApply) =
    min(duration(x.signal),inseconds(Float64,maybeseconds(x.time),samplerate(x)))

nsamples(x::AfterApply) = max(0,nsamples(x.signal) - resolvelen(x))
duration(x::AfterApply) =
    max(0,duration(x.signal) - inseconds(Float64,maybeseconds(x.time),samplerate(x)))

EvalTrait(x::AfterApply) = DataSignal()
function tosamplerate(x::UntilApply,s::IsSignal{<:Any,<:Number},c::ComputedSignal,fs;blocksize)
    CutApply(tosamplerate(child(x),fs;blocksize=blocksize),x.time,
        Val{:until}())
end
function tosamplerate(x::CutApply{<:Any,<:Any,K},s::IsSignal{<:Any,Missing},
    __ignore__,fs; blocksize) where K

    CutApply(tosamplerate(child(x),fs;blocksize=blocksize),x.time,K())
end

struct CutChunk{C} <: AbstractChunk
    n::Int
    child::C
end
child(x::CutChunk) = x.child

initchunk(x::UntilApply) = CutChunk(resolvelen(x),initchunk(child(x)))
function initchunk(x::AfterApply)
    skipped = 0
    len = resolvelen(x)
    chunk = initchunk(child(x))
    while !isnothing(chunk) && skipped < len
        chunk = nextchunk(x,chunk,len - skipped,true)
        skipped += nsamples(chunk)
    end
    @assert skipped == len
    CutChunk(len,chunk)
end

maxchunklen(x::UntilApply,chunk::CutChunk) =
    min(chunk.n,maxchunklen(child(x),child(chunk)))
function nextchunk(x::UntilApply,chunk::CutChunk,maxlen,skip)
    len = min(maxlen,maxchunklen(x,chunk))
    childchunk = nextchunk(x,child(chunk),len,skip)
    CutChunk(x.n - len,childchunk)
end

nextchunk(x::CutApply,chunk::CutChunk,maxlen,skip) =
    CutChunk(nextchunk(x,child(chunk),maxlen,skip))
nsamples(x::CutChunk) = nsamples(child(x))
@Base.propagate_inbounds sample(x::CutChunk,i) = sample(child(x),i)