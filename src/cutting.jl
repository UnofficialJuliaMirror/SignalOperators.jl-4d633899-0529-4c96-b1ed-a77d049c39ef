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

function nextchunk(x::AfterApply,maxlen,skip)
    len = resolvelen(x)
    childchunk = nextchunk(child(x),len,true)
    skipped = nsamples(childchunk)
    while !isnothing(childchunk) && skipped < len
        skipped += nsamples(childchunk)
        childchunk = nextchunk(child(x),min(maxlen,len - skipped),true,
            childchunk)
        isnothing(childchunk) && break
    end
    if skipped < len
        io = IOBuffer()
        signalshow(io,child(x))
        sig_string = String(take!(io))

        error("Signal is too short to skip $(maybeseconds(x.time)): ",
            sig_string)
    end
    @assert skipped == len
    nextchunk(x,maxlen,skip,CutChunk(0,childchunk))
end
function nextchunk(x::AfterApply,maxlen,skip,chunk::CutChunk)
    childchunk = nextchunk(child(x),maxlen,skip,child(chunk))
    if !isnothing(childchunk)
        CutChunk(0,childchunk)
    end
end
nextchunk(x::AfterApply,maxlen,skip,chunk::CutChunk{Nothing}) = nothing

initchunk(x::UntilApply) = CutChunk(resolvelen(x),nothing)
function nextchunk(x::UntilApply,len,skip,chunk::CutChunk=initchunk(x))
    nextlen = chunk.n - nsamples(chunk)
    if nextlen > 0
        childchunk = !isnothing(child(chunk)) ?
            nextchunk(child(x),min(nextlen,len),skip,child(chunk)) :
            nextchunk(child(x),min(nextlen,len),skip)
        if !isnothing(childchunk)
            CutChunk(nextlen,childchunk)
        end
    end
end

nsamples(x::CutChunk) = nsamples(child(x))
nsamples(x::CutChunk{Nothing}) = 0
@Base.propagate_inbounds sample(x::CutApply,chunk::CutChunk,i) =
    sample(child(x),child(chunk),i)