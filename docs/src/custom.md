# Custom Signals

Many signals can be readily created by passing a function to [`signal`](@ref) or by using [`mapsignal`](@ref). However, in a few cases it may be preferable to define a new subtype of `AbstractSignal`. This allows for the most flexibility in terms of how the signal will behave within a chain of signal operators.

The interface for `AbstractSignal` is built on the idea of "chunks" which are sequentially extracted from a signal. Chunks are short buffers of samples created from the signal. The contract of these chunks is that each sample from a chunk is guaranted to be read once or never (the samples may never be read if the signal is passed through [`after`](@ref) or [`until`](@ref)). Importantly the ordering of when samples are read from a chunk is *not* guaranteed.

Within [`sink`](@ref) and [`sink!](@ref), the funtion [`nextchunk`](@ref) is called repeatedly until it returns nothing. Given a signal and the last chunk (if it exists), it returns the next chunk of samples.

The chunks returned by [`nextchunk`](@ref) must implement two methods.

* [`nsamples`](@ref) Like a signal, each chunk has a fixed number of samples. Unlike signals, this cannot be an infinite value.
* [`sample`](@ref) Individuals samples of the chunk can be accessed by their index within the chunk (falling in the range of `1:nsamples(chunk)`).

This apporach to signals allows for a very tight inner loop where `sample` is called. Any slow code (e.g. branching or type unstable) should occur during `nextchunk`; `nsamples` and `sample` should both be fast, type-stable functions.

To make your new type of signal work with the rest of the operators you
probably also want to define an easy way to construct that signal by defining
a method of [`signal`](@ref). It should take two arguments: the first should
be of a type that can be used to construct your new signal type and the
second should be defined as `samplerate::Union{Missing,Number}=missing` which
should allow the user to specify the expected sample rate of the signal. You
can take additional arguments after these two, typically as keyword
arguments. If you follow this contract for `signal`, then any of the
operators defined by `SignalOperators` should be able to operate over your
signal seamlessly.