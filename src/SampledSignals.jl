using .SampledSignals: SampleBuf
function signal(x::SampleBuf,::IsSignal,fs::Union{Missing,Number}=missing)
    if !isconsistent(fs,samplerate(x))
        error("Signal expected to have sample rate of $fs Hz.")
    else
        x
    end
end
SignalTrait(::Type{<:SampleBuf{T}}) where T = IsSignal{T,Float64,Int}()
nsamples(x::SampleBuf) = size(x,1)
nchannels(x::SampleBuf) = size(x,2)
samplerate(x::SampleBuf) = SampledSignals.samplerate(x)

@Base.propagate_inbounds function sampleat!(result,x::SampleBuf,i,j,check)
    writesink!(result,i,view(x,j,:))
end

function sink(x,::Type{<:SampleBuf{El}};kwds...) where El
    x,n = process_sink_params(x;kwds...)
    result = SampleBuf{El}(n,nchannels(x))
    sink!(result,x)
end