
export rampon, rampoff, ramp, fadeto, sinramp

sinramp(x) = sinpi(0.5x)

struct RampSignal{D,S,Tm,Fn,T} <: WrappedSignal{S,T}
    signal::S
    time::Tm
    fn::Fn
end
function RampSignal(D,signal::S,time::Tm,fn::Fn) where {S,Tm,Fn}

    T = channel_eltype(signal)
    RampSignal{D,S,Tm,Fn,float(T)}(signal,time,fn)
end

SignalTrait(::Type{T}) where {S,T <: RampSignal{<:Any,S}} =
    SignalTrait(T,SignalTrait(S))
function SignalTrait(::Type{<:RampSignal{D,S,Tm,Fn,T}},::IsSignal{<:Any,Fs,L}) where
    {D,S,Tm,Fn,T,Fs,L}

    IsSignal{T,Fs,L}()
end

child(x::RampSignal) = x.signal
resolvelen(x::RampSignal) = max(1,insamples(Int,maybeseconds(x.time),samplerate(x)))

function tosamplerate(
    x::RampSignal{D},
    s::IsSignal{<:Any,<:Number},
    c::ComputedSignal,fs;blocksize) where D

    RampSignal(D,tosamplerate(child(x),fs;blocksize=blocksize),x.time,x.fn)
end
function tosamplerate(
    x::RampSignal{D},
    s::IsSignal{<:Any,Missing},
    __ignore__,fs; blocksize) where D

    RampSignal(D,tosamplerate(child(x),fs;blocksize=blocksize),x.time,x.fn)
end

struct RampChunk{Fn,T} <: AbstractChunk
    ramp::Fn
    marker::Int
    stop::Int
    offset::Int
    len::Int
end
RampChunk(x,fn,marker,stop,offset,len) =
    RampChunk{typeof(fn),float(channel_eltype(x))}(fn,marker,stop,offset,len)
nsamples(x::RampChunk) = x.len

sample(x::RampSignal,chunk::RampChunk{Nothing,T},i) where T =
    Fill(one(T),nchannels(x))
function sample(x::RampSignal{:on},chunk::RampChunk,i)
    ramplen = chunk.marker
    rampval = chunk.ramp((i-1) / ramplen)
    Fill(rampval,nchannels(x))
end
function sample(x::RampSignal{:off},chunk::RampChunk,i)
    startramp = chunk.marker
    rampval = chunk.stop > startramp ?
        rampval = chunk.ramp(1-(i - startramp)/(chunk.stop - startramp)) :
        rampval = chunk.ramp(1)
    Fill(rampval,nchannels(x))
end

initchunk(x::RampSignal{:on}) =
    RampChunk{:on}(x,x.fn,resolvelen(x),nsamples(x),0,0)
initchunk(x::RampSignal{:off}) =
    RampChunk{:off}(x,x.fn,nsamples(x) - resolvelen(x),nsamples(x),0,0)

nextchunk(x::RampSignal,maxlen,skip,chunk) =
    nextchunk(x,maxlen,skip,initchunk(x))
function nextchunk(x::RampSignal{:on},maxlen,skip,chunk::RampChunk)
    len = min(maxlen,chunk.marker - chunk.offset)
    offset = chunk.offset + chunk.len
    if offset < chunk.marker
        RampChunk{:on}(x,x.fn,chunk.marker,chunk.stop,offset,len)
    else
        RampChunk{:on}(x,nothing,chunk.marker,chunk.stop,offset,len)
    end
end

function nextchunk(x::RampSignal{:on},maxlen,skip,chunk::RampChunk{Nothing})
    len = min(maxlen,chunk.stop - chunk.offset)
    offset = chunk.offset + chunk.len
    RampChunk{:on}(x,nothing,chunk.marker,chunk.stop,offset,len)
end

function nextchunk(x::RampSignal{:off},maxlen,skip,chunk::RampChunk{Nothing})
    len = min(maxlen,chunk.marker - chunk.offset)
    offset = chunk.offset + chunk.len
    if offset < chunk.marker
        RampChunk{:off}(x,nothing,chunk.marker,chunk.stop,offset,len)
    else
        RampChunk{:off}(x,x.fn,chunk.marker,chunk.stop,offset,len)
    end
end

function nextchunk(x::RampSignal{:off},maxlen,skip,chunk::RampChunk)
    len = min(maxlen,chunk.stop - chunk.offset)
    offset = chunk.offset + chunk.len
    RampChunk{:off}(x,x.fn,chunk.marker,chunk.stop,offset,len)
end

function Base.show(io::IO, ::MIME"text/plain",x::RampSignal{D}) where D
    if x.fn isa typeof(sinramp)
        if D == :on
            write(io,"rampon_fn(",string(x.time),")")
        elseif D == :off
            write(io,"rampoff_fn(",string(x.time),")")
        else
            error("Reached unexpected code")
        end
    else
        if D == :on
            write(io,"rampon_fn(",string(x.time),",",string(x.fn),")")
        elseif D == :off
            write(io,"rampoff_fn(",string(x.time),",",string(x.fn),")")
        else
            error("Reached unexpected code")
        end
    end
end

"""

    rampon(x,[len=10ms],[fn=x -> sinpi(0.5x)])

Ramp the onset of a signal, smoothly transitioning from 0 to full amplitude
over the course of `len` seconds.

The function determines the shape of the ramp and should be non-decreasing
with a range of [0,1] over the domain [0,1]. It should map over the entire
range: that is `fn(0) == 0` and `fn(1) == 1`.

Both `len` and `fn` are optional arguments: either one or both can be
specified, though `len` must occur before `fn` if present.

"""
rampon(fun::Function) = rampon(10ms,fun)
rampon(len::Number=10ms,fun::Function=sinramp) = x -> rampon(x,len,fun)
function rampon(x,len::Number=10ms,fun::Function=sinramp)
    x = signal(x)
    x |> amplify(RampSignal(:on,x,len,fun))
end

"""

    rampoff(x,[len=10ms],[fn=x -> sinpi(0.5x)])

Ramp the offset of a signal, smoothly transitioning from full amplitude to 0
amplitude over the course of `len` seconds.

The function determines the shape of the ramp and should be non-decreasing
with a range of [0,1] over the domain [0,1]. It should map over the entire
range: that is `fn(0) == 0` and `fn(1) == 1`.

Both `len` and `fn` are optional arguments: either one or both can be
specified, though `len` must occur before `fn` if present.

"""
rampoff(fun::Function) = rampoff(10ms,fun)
rampoff(len::Number=10ms,fun::Function=sinramp) = x -> rampoff(x,len,fun)
function rampoff(x,len::Number=10ms,fun::Function=sinramp)
    x = signal(x)
    x |> amplify(RampSignal(:off,x,len,fun))
end

"""

    ramp(x,[len=10ms],[fn=x -> sinpi(0.5x)])

Ramp the onset and offset of a signal, smoothly transitioning from 0 to full
amplitude over the course of `len` seconds at the start and from full to 0
amplitude over the course of `len` seconds.

The function determines the shape of the ramp and should be non-decreasing
with a range of [0,1] over the domain [0,1]. It should map over the entire
range: that is `fn(0) == 0` and `fn(1) == 1`.

Both `len` and `fn` are optional arguments: either one or both can be
specified, though `len` must occur before `fn` if present.

"""
ramp(fun::Function) = ramp(10ms,fun)
ramp(len::Number=10ms,fun::Function=sinramp) = x -> ramp(x,len,fun)
function ramp(x,len::Number=10ms,fun::Function=sinramp)
    x = signal(x)
    x |> rampon(len,fun) |> rampoff(len,fun)
end

"""

    fadeto(x,y,[len=10ms],[fn=x->sinpi(0.5x)])

Append x to y, with a smooth transition lasting `len` seconds fading from
`x` to `y` (so the total length is `duration(x) + duration(y) - len`).

This fade is accomplished with a [`rampoff`](@ref) of `x` and a
[`rampon`](@ref) for `y`. `fn` should be non-decreasing with a range of [0,1]
over the domain [0,1]. It should map over the entire range: that is
`fn(0) == 0` and `fn(1) == 1`.

Both `len` and `fn` are optional arguments: either one or both can be
specified, though `len` must occur before `fn` if present.

"""
fadeto(y,fun::Function) = fadeto(y,10ms,fun)
fadeto(y,len::Number=10ms,fun::Function=sinramp) = x -> fadeto(x,y,len,fun)
function fadeto(x,y,len::Number=10ms,fun::Function=sinramp)
    x,y = uniform((x,y))
    x = signal(x)
    if ismissing(samplerate(x))
        error("Unknown sample rate is not supported by `fadeto`.")
    end
    n = insamples(Int,maybeseconds(len),samplerate(x))
    silence = signal(zero(channel_eltype(y))) |> until((nsamples(x) - n)*samples)
    x |> rampoff(len,fun) |> mix(
        y |> rampon(len,fun) |> prepend(silence))
end
