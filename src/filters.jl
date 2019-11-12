export lowpass, highpass, bandpass, bandstop, normpower, filtersignal

const default_blocksize = 2^12

struct FilterFn{D,M,A}
    design::D
    method::M
    args::A
end
(fn::FilterFn)(fs) =
    digitalfilter(fn.design(inHz.(fn.args)...,fs=inHz(fs)),fn.method)
filterfn(design,method,args...) = FilterFn(design,method,args)

function nyquist_check(x,hz)
    if !ismissing(samplerate(x)) && inHz(hz) ≥ 0.5samplerate(x)
        error("The frequency $(hz) cannot be represented at a sampling rate ",
              "of $(samplerate(x)) Hz. Increase the sampling rate or lower ",
              "the frequency.")
    end
end

"""
    lowpass(x,low;[order=5],[method=Butterworth(order)],[blocksize])

Apply a lowpass filter to x at the given cutoff frequency (`low`).
See [`filtersignal`](@ref) for details on `blocksize`.
"""
lowpass(low;kwds...) = x->lowpass(x,low;kwds...)
function lowpass(x,low;order=5,method=Butterworth(order),
    blocksize=default_blocksize)

    nyquist_check(x,low)
    filtersignal(x, filterfn(Lowpass,method,low), blocksize=blocksize)
end

"""
    highpass(x,high;[order=5],[method=Butterworth(order)],[blocksize])

Apply a highpass filter to x at the given cutoff frequency (`low`).
See [`filtersignal`](@ref) for details on `blocksize`.
"""
highpass(high;kwds...) = x->highpass(x,high;kwds...)
function highpass(x,high;order=5,method=Butterworth(order),
    blocksize=default_blocksize)

    nyquist_check(x,high)
    filtersignal(x, filterfn(Highpass,method,high),blocksize=blocksize)
end

"""
    bandpass(x,low,high;[order=5],[method=Butterworth(order)],[blocksize])

Apply a bandpass filter to x at the given cutoff frequencies (`low` and `high`).
See [`filtersignal`](@ref) for details on `blocksize`.
"""
bandpass(low,high;kwds...) = x->bandpass(x,low,high;kwds...)
function bandpass(x,low,high;order=5,method=Butterworth(order),
    blocksize=default_blocksize)

    nyquist_check(x,low)
    nyquist_check(x,high)
    filtersignal(x, filterfn(Bandpass,method,low,high),blocksize=blocksize)
end

"""
    bandstop(x,low,high;[order=5],[method=Butterworth(order)],[blocksize])

Apply a bandstop filter to x at the given cutoff frequencies (`low` and `high`).
See [`filtersignal`](@ref) for details on `blocksize`.
"""
bandstop(low,high;kwds...) = x->bandstop(x,low,high;kwds...)
function bandstop(x,low,high;order=5,method=Butterworth(order),
    blocksize=default_blocksize)

    nyquist_check(x,low)
    nyquist_check(x,high)
    filtersignal(x, filterfn(Bandstop,method,low,high),blocksize=blocksize)
end

"""
    filtersignal(x,h;[blocksize])

Apply the given filter `h` (from [`DSP`](https://github.com/JuliaDSP/DSP.jl))
to signal `x`.

## Blocksize

Blocksize determines the size of the buffer used when computing intermediate
values of the filter. It defaults to 4096. It need not normally be adjusted.

"""
filtersignal(h;blocksize=default_blocksize) =
    x -> filtersignal(x,h;blocksize=blocksize)
filtersignal(x,fn::Union{FilterFn,Function};kwds...) =
    filtersignal(x,SignalTrait(x),fn;kwds...)
filtersignal(x,h;kwds...) =
    filtersignal(x,SignalTrait(x),RawFilterFn(h);kwds...)
function filtersignal(x,::Nothing,args...;kwds...)
    filtersignal(signal(x),args...;kwds...)
end

struct RawFilterFn{H}
    h::H
end
(fn::RawFilterFn)(fs) = deepcopy(fn.h)

resolve_filter(x) = DSP.Filters.DF2TFilter(x)
resolve_filter(x::FIRFilter) = x
function filtersignal(x::Si,s::IsSignal,fn;blocksize,newfs=samplerate(x)) where {Si}
    FilteredSignal(x,fn,blocksize,newfs)
end
struct FilteredSignal{T,Si,Fn,Fs} <: WrappedSignal{Si,T}
    signal::Si
    fn::Fn
    blocksize::Int
    samplerate::Fs
end
function FilteredSignal(signal::Si,fn::Fn,blocksize::Number,newfs::Fs) where {Si,Fn,Fs}
    T = float(channel_eltype(signal))
    FilteredSignal{T,Si,Fn,Fs}(signal,fn,Int(blocksize),newfs)
end
SignalTrait(x::Type{T}) where {S,T <: FilteredSignal{<:Any,S}} =
    SignalTrait(x,SignalTrait(S))
SignalTrait(x::Type{<:FilteredSignal{T}},::IsSignal{<:Any,Fs,L}) where {T,Fs,L} =
    IsSignal{T,Fs,L}()
child(x::FilteredSignal) = x.signal
samplerate(x::FilteredSignal) = x.samplerate
EvalTrait(x::FilteredSignal) = ComputedSignal()

Base.show(io::IO,::MIME"text/plain",x::FilteredSignal) = pprint(io,x)
function PrettyPrinting.tile(x::FilteredSignal)
    child = signaltile(x.signal)
    operate = literal(filterstring(x.fn))
    tilepipe(child,operate)
end
signaltile(x::FilteredSignal) = PrettyPrinting.tile(x)
filterstring(fn::FilterFn) =
    string(filterstring(fn.design),"(",join(string.(fn.args),","),")")
filterstring(fn::Function) = string("filtersignal(",string(fn),")")
function filtertring(fn::RawFilterFn)
    io = IOBuffer()
    show(IOContext(io,:displaysize=>(1,30),:limit=>true),
        MIME("text/plain"),x)
    string("filtersignal(",String(take!(io)),")")
end
filterstring(::Type{<:Lowpass}) = "lowpass"
filterstring(::Type{<:Highpass}) = "highpass"
filterstring(::Type{<:Bandpass}) = "bandpass"
filterstring(::Type{<:Bandstop}) = "bandstop"
filterstring(x) = string(x)

function tosamplerate(x::FilteredSignal,s::IsSignal{<:Any,<:Number},::ComputedSignal,fs;
blocksize)
    # is this a non-resampling filter?
    if samplerate(x) == samplerate(x.signal)
        FilteredSignal(tosamplerate(x.signal,fs,blocksize=blocksize),
            x.fn,x.blocksize,fs)
    else
        tosamplerate(x.signal,s,DataSignal(),fs,blocksize=blocksize)
    end
end
function tosamplerate(x::FilteredSignal,::IsSignal{<:Any,Missing},__ignore__,fs;
        blocksize)
    FilteredSignal(tosamplerate(x.signal,fs,blocksize=blocksize),
        x.fn,x.blocksize,fs)
end

function nsamples(x::FilteredSignal)
    if ismissing(samplerate(x.signal))
        missing
    elseif samplerate(x) == samplerate(x.signal)
        nsamples(x.signal)
    else
        ceil(Int,nsamples(x.signal)*samplerate(x)/samplerate(x.signal))
    end
end

struct FilterChunk{H,S,T,C}
    len::Int
    last_output_index::Int

    last_input_offset::Int
    last_output_offset::Int

    hs::Vector{H}
    input::Matrix{S}
    output::Matrix{T}

    child::C
end
init_length(x::FilteredSignal) = min(nsamples(x),x.blocksize)
init_length(x::FilteredSignal{<:Any,<:Any,<:ResamplerFn}) =
    trunc(Int,min(nsamples(x),x.blocksize) / x.fn.ratio)
function FilterChunk(x::FilteredSignal)
    hs = [resolve_filter(x.fn(samplerate(x))) for _ in 1:nchannels(x.signal)]
    len = init_length(x)
    input = Array{channel_eltype(x.signal)}(undef,len,nchannels(x))
    output = Array{channel_eltype(x)}(undef,x.blocksize,nchannels(x))

    FilterChunk(0,size(output,1), 0,0, hs,input,output,nothing)
end
nsamples(x::FilterChunk) = x.len
@Base.propagate_inbounds sample(::FilteredSignal,x::FilterChunk,i) =
    view(x.state.output,i+x.last_output_index,:)

inputlength(x,n) = n
outputlength(x,n) = n
inputlength(x::DSP.Filters.Filter,n) = DSP.inputlength(x,n)
outputlength(x::DSP.Filters.Filter,n) = DSP.outputlength(x,n)

function nextchunk(x::FilteredSignal,maxlen,skip,
    chunk::FilterChunk=FilterChunk(x))

    last_output_index = chunk.last_output_index + chunk.len

    if last_output_index < size(chunk.output,1)
        len = min(maxlen, size(chunk.output) - last_output_index)
        FilterChunk(len, last_output_index, chunk.last_input_offset,
            chunk.last_output_offset, chunk.hs, chunk.input, chunk.output,
            chunk.child)
    else
        len = min(maxlen,size(chunk.output,1))

        psig = pad(x.signal,zero)
        childchunk = !isnothing(child(chunk)) ? child(chunk) :
            nextchunk(psig,size(chunk.input,1),false)
        childchunk = sink!(chunk.input,psig,SignalTrait(psig),childchunk)
        last_input_offset = chunk.last_input_offset + size(chunk.input,1)

        # filter the input to the output buffer
        out_len = outputlength(chunk.hs[1],size(chunk.input,1))
        for ch in 1:size(chunk.output,2)
            filt!(view(chunk.output,1:out_len,ch),chunk.hs[ch],
                    view(chunk.input,:,ch))
        end
        last_output_offset = chunk.last_output_offset + out_len
        last_output_index = 0

        FilterChunk(len, last_output_index, last_input_offset,
            last_output_offset, chunk.hs, chunk.input, chunk.output, childchunk)
    end
end

# TODO: create an online version of normpower?
# TODO: this should be excuted lazzily to allow for unkonwn samplerates
struct NormedSignal{Si,T} <: WrappedSignal{Si,T}
    signal::Si
end
child(x::NormedSignal) = x.signal
nsamples(x::NormedSignal) = nsamples(x.signal)
NormedSignal(x::Si) where Si = NormedSignal{Si,float(channel_eltype(Si))}(x)
SignalTrait(x::Type{T}) where {S,T <: NormedSignal{S}} =
    SignalTrait(x,SignalTrait(S))
SignalTrait(x::Type{<:NormedSignal{<:Any,T}},::IsSignal{<:Any,Fs,L}) where {T,Fs,L} =
    IsSignal{T,Fs,L}()
function tosamplerate(x::NormedSignal,s::IsSignal{<:Any,<:Number},
    ::ComputedSignal,fs;blocksize)

    NormedSignal(tosamplerate(x.signal,fs,blocksize=blocksize))
end
function tosamplerate(x::NormedSignal,::IsSignal{<:Any,Missing},
    __ignore__,fs;blocksize)

    NormedSignal(tosamplerate(x.signal,fs,blocksize=blocksize))
end

struct NormedChunk{A}
    offset::Int
    len::Int
    vals::A
end
nsamples(x::NormedChunk) = x.len
@Base.propagate_inbounds sample(x::NormedChunk,i) = view(x.vals,i,:)

function initchunk(x::NormedSignal)
    if isinf(nsamples(x))
        error("Cannot normalize an infinite-length signal. Please ",
              "use `until` to take a prefix of the signal")
    end
    vals = Array{channel_eltype(x)}(undef,nsamples(x),nchannels(x))
    sink!(vals, x.signal)

    rms = sqrt(mean(x -> float(x)^2,vals))
    vals ./= rms

    S,V = typeof(x), typeof(vals)
    NormedChunk(0,0,vals)
end

function nextchunk(x::NormedSignal,maxlen,skip,chunk::NormedChunk)
    len = min(maxlen,nsamples(x) - chunk.offset)
    NormedChunk(chunk.offset + chunk.len, len, chunk.vals)
end

"""
    normpower(x)

Return a signal with normalized power. That is, divide all samples by the
root-mean-squared value of the entire signal.

"""
function normpower(x)
    x = signal(x)
    NormedSignal{typeof(x),float(channel_eltype(typeof(x)))}(signal(x))
end

Base.show(io::IO,::MIME"text/plain",x::NormedSignal) = pprint(io,x)
function PrettyPrinting.tile(x::NormedSignal)
    tilepipe(signaltile(x.signal),literal("normpower"))
end
signaltile(x::NormedSignal) = PrettyPrinting.tile(x)