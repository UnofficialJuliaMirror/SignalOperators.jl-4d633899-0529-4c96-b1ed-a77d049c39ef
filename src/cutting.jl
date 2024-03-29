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

struct AfterCheckpoint{S,C} <: AbstractCheckpoint{S}
    diff::Int
    child::C
end
checkindex(c::AfterCheckpoint) = checkindex(c.child)-c.diff
# No long offset, use 'stopat' index
function atcheckpoint(x::AfterApply,offset::Number,stopat::Int)
    n = resolvelen(x)
    childcheck = atcheckpoint(child(x),offset+n,stopat+n)
    if !isnothing(childcheck)
        AfterCheckpoint{typeof(x),typeof(childcheck)}(n,childcheck)
    end
end

function atcheckpoint(x::S,check::AfterCheckpoint{S},stopat::Int) where
    S <: AfterApply

    n = resolvelen(x)
    childcheck = atcheckpoint(child(x),check.child,stopat+n)
    if !isnothing(childcheck)
        AfterCheckpoint{S,typeof(childcheck)}(n,childcheck)
    end
end


atcheckpoint(x::UntilApply,offset::Number,stopat) =
    atcheckpoint(child(x),offset,stopat)
atcheckpoint(x::UntilApply,offset::AbstractCheckpoint,stopat) =
    atcheckpoint(child(x),offset,stopat)
@Base.propagate_inbounds function sampleat!(result,x::AfterApply,i,j,check)
    sampleat!(result,x.signal,i,j+check.diff,check.child)
end

@Base.propagate_inbounds function sampleat!(result,x::UntilApply,i,j,check)
    sampleat!(result,x.signal,i,j,check)
end