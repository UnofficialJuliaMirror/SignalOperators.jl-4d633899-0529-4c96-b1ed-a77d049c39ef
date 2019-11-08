export sink
using AxisArrays

errordim() = error("To treat an array as a signal it must have 1 or 2 dimensions")

# signals can be arrays with some metadata
"""
## Arrays

Arrays can be treated as signals. The first dimension is time, the second
channels.

[`AxisArrays`](https://github.com/JuliaArrays/AxisArrays.jl), if they have an
axis labeled `time` and one or zero additional axes, can be treated as a
signal. The time dimension must be represented using on object with the `step`
function defined (e.g. any `AbstractRange`).

"""
function signal(x::AbstractArray{<:Any,N},::IsSignal,
    fs::Union{Missing,Number}=missing) where N

    if N == 1
        ismissing(fs) && return x
        times = range(0s,length=size(x,1),step=float(s/inHz(fs)))
        AxisArray(x,Axis{:time}(times))
    elseif N == 2
        ismissing(fs) && return x
        times = range(0s,length=size(x,1),step=float(s/inHz(fs)))
        channels = 1:size(x,2)
        AxisArray(x,Axis{:time}(times),Axis{:channel}(channels))
    else
        errordim()
    end
end
function signal(x::AxisArray,::IsSignal,fs::Union{Missing,Number}=missing)
    if !isconsistent(fs,samplerate(x))
        error("Signal expected to have sample rate of $(inHz(fs)) Hz.")
    else
        x
    end
end

tosamplerate(x::AbstractArray,::IsSignal{<:Any,Missing},::DataSignal,fs::Number;blocksize) =
    signal(x,fs)
tosamplerate(x::AxisArray,s::IsSignal{<:Any,<:Number},::DataSignal,fs::Number;blocksize) =
    __tosamplerate__(x,s,fs,blocksize)


function SignalTrait(::Type{A}) where{T,N,A<:AbstractArray{T,N}}
    if N ∈ [1,2]
        if A <: AxisArray
            IsSignal{T,Float64,Int}()
        else
            IsSignal{T,Missing,Int}()
        end
    else
        error("Array must have 1 or 2 dimensions to be treated as a signal.")
    end
end

nsamples(x::AxisArray) = length(AxisArrays.axes(x,Axis{:time}))
nsamples(x::AbstractVecOrMat) = size(x,1)

function nchannels(x::AxisArray)
    chdim = axisdim(x,Axis{:time}) == 1 ? 2 : 1
    size(x,chdim)
end
nchannels(x::AbstractVecOrMat) = size(x,2)
function samplerate(x::AxisArray)
    times = axisvalues(AxisArrays.axes(x,Axis{:time}))[1]
    inHz(1/step(times))
end
samplerate(x::AbstractVecOrMat) = missing

const WithAxes{Tu} = AxisArray{<:Any,<:Any,<:Any,Tu}
const AxTimeD1 = Union{
    WithAxes{<:Tuple{Axis{:time}}},
    WithAxes{<:Tuple{Axis{:time},<:Any}}}
const AxTimeD2 = WithAxes{<:Tuple{<:Any,Axis{:time}}}
const AxTime = Union{AxTimeD1,AxTimeD2}

struct ArrayChunk{A} <: AbstractChunk
    offset::Int
    ar::A
end
nsamples(x::ArrayChunk) = size(x.ar,1)
sample(x::ArrayChunk,i) = view(x.ar,i,:)

maxchunklen(x::AxTime,chunk::ArrayChunk,maxlen,skip) =
    min(maxlen,size(x,timedim(x)) - chunk.offset)
timedim(x::AxTimeD1) = 1
timedim(x::AxTimeD2) = 2

initchunk(x::AxTime) = ArrayChunk(0,1:0)
function nextchunk(x::AxTime,chunk::ArrayChunk,maxlen,skip)
    offset = chunk.offset + nsamples(chunk)
    if offset < nsamples(x)
        len = maxchunklen(x,chunk,maxlen,skip)
        indices = offset .+ (1:len)
        array_chunk_helper(x,offset,indices)
    else
        nothing
    end
end
array_chunk_helper(x::AxTimeD1,offset,indices) =
    ArrayChunk(offset, view(x,indices,:))
array_chunk_helper(x::AxTimeD2,offset,indices) =
    ArrayChunk(offset, view(x,:,indices))

function signaltile(x)
    io = IOBuffer()
    signalshow(io,x)
    literal(String(take!(io)))
end
function signalshow(io,x::AbstractArray)
    show(IOContext(io,:displaysize=>(1,30),:limit=>true),
        MIME("text/plain"),x)
    show_fs(io,x)
end
function signalshow(io,x::AxisArray)
    signalshow(io,x.data)
    show_fs(io,x)
end