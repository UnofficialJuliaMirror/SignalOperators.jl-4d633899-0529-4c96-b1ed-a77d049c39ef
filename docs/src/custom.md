# Custom Signals

Many signals can be readily created by passing a function to [`signal`](@ref) or by using [`mapsignal`](@ref). However, in a few cases it may be preferable to define a new subtype of `AbstractSignal`. This allows for the most flexibility in terms of how the signal will behave within a chain of signal operators.

The interface for an `AbstractSignal` includes the following functions which are called during [`sink!`](@ref) and [`sink`](@ref).

* [`initchunk`](@ref) - initializes an empty chunk to be fed into `nextchunk`, no samples will be read from the chunk.
* [`maxchunklen`](@ref) - gets the maximum permitted sample count of the next `chunk`
* [`nextchunk`](@ref) - given a signal and the last chunk, return the next chunk of samples or `nothing` if there are no samples left for the given signal.

This interface is built on the idea of "chunks" which are sequentially extracted from a signal. Chunks are short buffers of samples created from the signal. The contract of these chunks is that each sample from a chunk is guaranted to be read once or never (the samples may never be read if the signal is passed through [`after`](@ref) or [`until`](@ref)). Importantly the ordering of when samples are read from a chunk is *not* guaranteed.

The idea is that `sink!` first calls `initchunk` and then sequentially requests each chunk from a signal, using `nextchunk` to request a chunk of a given sample length. The number of samples in the next chunk is no greater than `maxchunklen` and the returned chunk is always exactly as long as the requested length so long as it is less than this maximum length. The value of `maxchunklen` is used, for instance, by [`mapsignal`](@ref) to make sure the returned chunk for each signal being operated on for a multi-argument function is of the same length.

The chunks returned by `nextchunk` must implement two methods.

* [`nsamples`](@ref) Like a signal, each chunk has a fixed number of samples. Unlike signals, this cannot be an infinite value.
* [`sample`](@ref) Individuals samples of the chunk can be accessed by their index within the chunk (falling in the range of `1:nsamples(chunk)`).

This apporach to signals allows for a very tight inner loop (or parlallized code) where `sample` is called. Any slow code (e.g. branching or type unstable) should occur during `initchunk` and `nextchunk`; `nsamples` and `sample` should all be fast, type-stable functions. (In addition, `maxchunklen` should also ideally be fast, as it can be called multiple times).

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