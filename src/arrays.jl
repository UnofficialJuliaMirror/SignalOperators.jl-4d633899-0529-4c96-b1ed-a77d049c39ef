export sink
using AxisArrays

errordim() = error("To treat an array as a signal it must have 1 or 2 dimensions")

# signals can be arrays with some metadata
function signal(x::AbstractArray{<:Any,N},
    fs::Union{Missing,Number}=missing) where N

    if N == 1
        ismissing(fs) && return x
        times = range(0s,length=size(x,1),step=s/inHz(fs))
        AxisArray(x,Axis{:time}(times))
    elseif N == 2
        ismissing(fs) && return x
        times = range(0s,length=size(x,1),step=s/inHz(fs))
        channels = 1:size(x,2)
        AxisArray(x,Axis{:time}(times),Axis{:channel}(channels))
    else
        errordim()
    end
end

function signal(x::AxisArray,fs::Union{Missing,Number}=missing)
    times = axisvalues(AxisArrays.axes(x,Axis{:time}))[1]
    !isconsistent(fs,1/step(times))
    x
end

function SignalTrait(::Type{A}) where{T,N,A<:AbstractArray{T,N}}
    if N ∈ [1,2]
        if A isa AxisArray
            IsSignal{T,Float64,Int}()
        else
            IsSignal{T,Missing,Int}()
        end
    else
        errordim()
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
@Base.propagate_inbounds function sampleat!(result,x::AxTimeD1,
    ::IsSignal,i::Number,j::Number,check)

    writesink(result,i,view(x,j,:))
end
@Base.propagate_inbounds function sampleat!(result,x::AxTimeD2,
    ::IsSignal,i::Number,j::Number,check)

    writesink(result,i,view(x,:,j))
end