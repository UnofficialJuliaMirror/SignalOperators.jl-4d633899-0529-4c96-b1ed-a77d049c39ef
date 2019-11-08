struct NumberSignal{T,S,DB} <: AbstractSignal{T}
    val::T
    samplerate::S
end

NumberSignal(x::T,sr::Fs;dB=false) where {T,Fs} = NumberSignal{T,Fs,dB}(x,sr)
function Base.show(io::IO, ::MIME"text/plain", x::NumberSignal{<:Any,<:Any,true})
    show(io,MIME("text/plain"), uconvertrp(Units.dB, x.val))
    show_fs(io,x)
end
function Base.show(io::IO, ::MIME"text/plain", x::NumberSignal{<:Any,<:Any,false})
    show(io, MIME("text/plain"), x.val)
    show_fs(io,x)
end

"""

## Numbers

Numbers can be treated as infinite length, constant signals of unknown
sample rate.

"""
signal(val::Number,::Nothing,fs) = NumberSignal(val,inHz(Float64,fs))
signal(val::Unitful.Gain{<:Any,<:Any,T},::Nothing,fs) where T =
    NumberSignal(float(T)(uconvertrp(NoUnits,val)),inHz(Float64,fs),dB=true)

SignalTrait(::Type{<:NumberSignal{T,S}}) where {T,S} = IsSignal{T,S,InfiniteLength}()

nchannels(x::NumberSignal) = 1
nsamples(x::NumberSignal) = inflen
samplerate(x::NumberSignal) = x.samplerate

tosamplerate(x::NumberSignal{<:Any,<:Any,DB},::IsSignal,::ComputedSignal,
    fs=missing;blocksize) where DB = NumberSignal(x.val,fs,dB=DB)

struct NumberChunk{T}
    val::T
    len::Int
end
initchunk(x::NumberSignal) = NumberChunk(x.val,0)
nextchunk(x::NumberSignal,chunk,maxlen,skip) = NumberChunk(x.val,maxlen)
maxchunklen(x::NumberSignal,chunk,maxlen,skip) = min(maxlen,inflen)
nsamples(chunk::NumberChunk) = chunk.len
sample(chunk::NumberChunk,i) = chunk.val
