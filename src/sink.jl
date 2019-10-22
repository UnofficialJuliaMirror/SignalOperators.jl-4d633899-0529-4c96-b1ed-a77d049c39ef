
"""
    sink([signal],[to=AxisArray];duration,samplerate)

Creates a given type of object (`to`) from a signal. By default it is an
`AxisArray` with time as the rows and channels as the columns.

# Keyword arguments

## Sample Rate

The sample rate does not need to be specified, it will use either the sample
rate of `signal` or a default sample rate (which raises a warning). If
specified, the given sample rate is passed to [`signal`](@ref) when coercing
the input to a signal.

## Duration

You can limit the output of the given signal to the specified duration. If
this duration exceedes the duration of the passed signal an error will be
thrown.

# Values for `to`

## Type

If `to` is a type (e.g. `Array`) the signal is written to a value of that
type.

"""
sink(to::Type=AxisArray;kwds...) = x -> sink(x,to;kwds...)
sink(x;kwds...) = sink(x,AxisArray;kwds...)
function sink(x,::Type{<:AxisArray};kwds...)
    result = sink(x,Array;kwds...)
    times = Axis{:time}(range(0s,length=size(result,1),step=float(s/samplerate(x))))
    channels = Axis{:channel}(1:nchannels(x))
    AxisArray(result,times,channels)
end
function sink(x,::Type{<:Array{El}};kwds...) where El
    x,n = process_sink_params(x;kwds...)
    result = Array{El}(undef,n,nchannels(x))
    sink!(result,x)
end

process_sink_params(x;duration=missing,
    samplerate=SignalOperators.samplerate(x))
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

    x,n
end

"""

## Filename

If `to` is a string, it is assumed to describe the name of a file to which
the signal will be written. You will need to call `import` or `using` on an
appropriate backend for writing to the given file type.

Available backends include the following pacakges
- [WAV](https://codecov.io/gh/haberdashPI/SignalOperators.jl/src/master/src/WAV.jl)
- [LibSndFile](https://github.com/JuliaAudio/LibSndFile.jl)

"""
sink(to::String;kwds...) = x -> sink(x,to;kwds...)
function sink(x,to::String;kwds...)
    x,n = process_sink_params(x;kwds...)
    save_signal(filetype(to),to,x,n)
end
function save_signal(::Val{T},filename,x,len) where T
    error("No backend loaded for file of type $T. Refer to the documentation ",
          "of `signal` to find a list of available backends.")
end

"""
    sink!(array,x;[samplerate],[offset])

Write `size(array,1)` samples of signal `x` to `array`, starting from the sample after
`offset`. If no sample rate has been specified for `x` you can specify it
now, using `samplerate` (it will default to 44.1kHz).

"""
sink!(result::Union{AbstractVector,AbstractMatrix};kwds...) =
    x -> sink!(result,x;kwds...)
function sink!(result::Union{AbstractVector,AbstractMatrix},x;
    samplerate=SignalOperators.samplerate(x),offset=0)

    if ismissing(samplerate) && ismissing(SignalOperators.samplerate(x))
        @warn("No sample rate was specified, defaulting to 44.1 kHz.")
        samplerate = 44.1kHz
    end
    x = signal(x,samplerate)

    if nsamples(x)-offset < size(result,1)
        error("Signal is too short to fill buffer of length $(size(result,1)).")
    end
    x = tochannels(x,size(result,2))

    sink!(result,x,SignalTrait(x),offset)
end

abstract type AbstractCheckpoint{S}
end
struct EmptyCheckpoint{S} <: AbstractCheckpoint{S}
    n::Int
end
checkindex(x::EmptyCheckpoint) = x.n
child(x::Number) = x
checkindex(x::Number) = x+1

atcheckpoint(x,offset::Number,stopat) =
    offset ≥ stopat ? nothing : EmptyCheckpoint{typeof(x)}(1+offset)
atcheckpoint(x::S,check::AbstractCheckpoint{S},stopat) where S =
    checkindex(check) == stopat+1 ? nothing : EmptyCheckpoint{S}(stopat+1)
atcheckpoint(x,check) =
    error("Internal error: signal inconsistent with checkpoint")

fold(x) = zip(x,Iterators.drop(x,1))
sink!(result,x,sig::IsSignal,offset::Number) =
    sink!(result,x,sig,offset,atcheckpoint(x,offset,size(result,1)+offset))
function sink!(result,x,sig::IsSignal,offset,check)
    written = 0
    while !isnothing(check) && written < size(result,1)
        next = atcheckpoint(x,check,offset+size(result,1))
        isnothing(next) && break

        len = checkindex(next) - checkindex(check)
        @assert len > 0

        sink_helper!(result,offset,x,check,len)
        written += len
        check = next
    end
    @assert written == size(result,1)
    result
end

@noinline function sink_helper!(result,offset,x,check,len)
    @inbounds @simd for i in checkindex(check):(checkindex(check)+len-1)
        sampleat!(result,x,i-offset,i,check)
    end
end

@Base.propagate_inbounds function writesink!(result::AbstractArray,i,v)
    for ch in 1:length(v)
        result[i,ch] = v[ch]
    end
    v
end
