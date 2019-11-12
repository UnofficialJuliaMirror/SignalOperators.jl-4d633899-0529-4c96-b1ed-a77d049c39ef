
"""
    sink([signal],[to=AxisArray];duration,samplerate)

Creates a given type of object (`to`) from a signal. By default it is an
`AxisArray` with time as the rows and channels as the columns. If a filename
is specified for `to`, the signal is written to the given file. If given a
type (e.g. `Array`) the signal is written to a value of that type.

# Sample Rate

The sample rate does not need to be specified, it will use either the sample
rate of `signal` or a default sample rate (which raises a warning). If
specified, the given sample rate is passed to [`signal`](@ref) when coercing
the input to a signal.

# Duration

You can limit the output of the given signal to the specified duration. If
this duration exceedes the duration of the passed signal an error will be
thrown.

"""
sink(to::Type=AxisArray;kwds...) = x -> sink(x,to;kwds...)
function sink(x::T,::Type{A}=AxisArray;
        duration=missing,
        samplerate=SignalOperators.samplerate(x)) where {T,A}

    if ismissing(samplerate) && ismissing(SignalOperators.samplerate(x))
        @warn("No sample rate was specified, defaulting to 44.1 kHz.")
        samplerate = 44.1kHz
    end
    x = signal(x,samplerate)
    duration = coalesce(duration,nsamples(x)*samples)

    if isinf(duration)
        error("Cannot store infinite signal. Specify a finite duration when ",
            "calling `sink`.")
    end

    n = insamples(Int,maybeseconds(duration),SignalOperators.samplerate(x))
    if n > nsamples(x)
        error("Requested signal duration is too long for passed signal: $x.")
    end

    sink(x,SignalTrait(x),n,A)
end

function sink(x,sig::IsSignal{El},len::Number,::Type{<:Array}) where El
    result = Array{El}(undef,len,nchannels(x))
    sink!(result,x)
end
function sink(x,sig::IsSignal{El},len,::Type{<:AxisArray}) where El
    result = sink(x,sig,len,Array)
    times = Axis{:time}(range(0s,length=size(result,1),step=float(s/samplerate(x))))
    channels = Axis{:channel}(1:nchannels(x))
    AxisArray(result,times,channels)
end
sink(x, ::IsSignal, ::Nothing, ::Type) = error("Don't know how to interpret value as a signal: $x")

"""
    sink!(array,x;[samplerate])

Write `size(array,1)` samples of signal `x` to `array`. If no sample rate has
been specified for `x` you can specify it now, using `samplerate` (it will
default to 44.1kHz).

"""
sink!(result::Union{AbstractVector,AbstractMatrix};kwds...) =
    x -> sink!(result,x;kwds...)
function sink!(result::Union{AbstractVector,AbstractMatrix},x;
    samplerate=SignalOperators.samplerate(x))

    if ismissing(samplerate) && ismissing(SignalOperators.samplerate(x))
        @warn("No sample rate was specified, defaulting to 44.1 kHz.")
        samplerate = 44.1kHz
    end
    x = signal(x,samplerate)

    if nsamples(x) < size(result,1)
        error("Signal is too short to fill buffer of length $(size(result,1)).")
    end
    x = tochannels(x,size(result,2))

    sink!(result,x,SignalTrait(x))
    result
end

abstract type AbstractChunk; end
function sample
end
function nextchunk
end

fold(x) = zip(x,Iterators.drop(x,1))
sink!(result,x,sig::IsSignal) =
    sink!(result,x,sig,nextchunk(x,size(result,1),false))
function sink!(result,x,::IsSignal,chunk)
    written = 0
    while !isnothing(chunk) && written < size(result,1)
        @assert nsamples(chunk) > 0
        sink_helper!(result,written,x,chunk)
        written += nsamples(chunk)
        maxlen = size(result,1)-written
        chunk = nextchunk(x,maxlen,false,chunk)
    end
    @assert written == size(result,1)

    chunk
end

@noinline function sink_helper!(result,written,x,chunk)
    @inbounds @simd for i in 1:nsamples(chunk)
        writesink!(result,i+written,sample(x,chunk,i))
    end
end

@Base.propagate_inbounds function writesink!(result::AbstractArray,i,v)
    for ch in 1:length(v)
        result[i,ch] = v[ch]
    end
    v
end
