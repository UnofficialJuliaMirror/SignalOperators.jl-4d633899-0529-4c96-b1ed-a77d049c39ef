
SignalOperators is a package that aims to provide a clean interface for generating and manipulating regularly sampled signals, typically sounds.

# Key concepts

There are several important concepts employed across the public interface. Let's step through one of the examples from the homepage (and README.md), which demonstrates most of these concepts.

```julia
sound1 = signal(sin,ω=1kHz) |> until(5s) |> ramp |> normpower |> amplify(-20dB)
```
This example creates a 1 kHz pure-tone (sine wave) that lasts 5 serconds. It's amplitude is 20 dB lower than a signal with unit 1 power.

There are a few things going on here: piping, the use of units, infinite length signals, and unspecified sample rates 

## Piping

All of the functions implemented in SignalOperators can be piped. This means
that instead of passing the first argument to a signal operator, you can pipe it using `|>`. For example, the two statements below have the same meaning.

```julia
sound1 = signal(sin,ω=1kHz) |> until(5s)
sound1 = until(signal(sin,ω=1kHz),5s)
```

The use of piping makes it easier to read the sequence of operations that
are performed on the signal.

## Units

In any place where a signal operator needs a time or a frequency, it can be specified using units. If units are not specified, time is assumed to be in seconds, and frequency in Hertz. So, for example, the following two statements
are equivalent.

```julia
sound1 = signal(sin,ω=1kHz)
sound1 = signal(sin,ω=1000)
```

To make use of the unit variables you must call `using SignalOperators.Units`.
The units defined are `Hz`, `s` `kHz` `ms` `°`, `rad` and `dB`. You can just include the ones you want using e.g. `using SignalOperators.Units: Hz`, or you can include more by adding the [`Unitful`](https://github.com/PainterQubits/Unitful.jl) package to your project and adding the desired units from there.

Note that the output of all functions to inspect a signal (e.g. `duration`, `samplerate`) are `Float64` values in the default unit (seconds or Hertz).

### decibels

You can pass an amplification value as unitless or a unitful value in `dB`,
a unitless value is not assumed to be in decibels. Instead it's assumed to be the actual ratio by which you wish to multiply the signal. E.g. `amplify(x,2)`
will make x twice as loud.

## Infinite lengths

Some of the ways you can define a signal lead to an infinite length signal.
You cannot store an infinite signal. It is represented as a function of some
kind. Operations on signals are generally lazy, meaning the samples of the
signal aren't evaluated until necessary. To allow actual data to be created
from a signal, you have to specify the length, using [`until`](@ref). For example,
when using `signal(sin)`, the signal is an infinite length sine wave. That's
why, in the example above we use [`until`](@ref) to specify the length, like so:

```julia
signal(sin,ω=1kHz) |> until(5s)
```

Infinite lengths are represented as the value `inflen`. This has overloaded definitions of various operators to play nicely with ordering, arithmetic etc...

## Unspecified sample rates

You may notice that the above signal has no defined sample rate. Such a signal
is defined by a function, and can be sampled at whatever rate you desire. If you add a signal to the chain of operations that does have a defined sample rate, the unspecified sample rate will be resolved to that same rate (see Signal promotion, below). If there is no defined sample rate by the time you call [`sink`](@ref) you can specify it then.

## Sinking

Once you have defined a signal, you can create some concrete sequence of samples from it. This is done using [`sink`](@ref). The resulting value is itself a signal, so you can pass this to other signal operators. The function [`sink`](@ref) is also used to create a file. Sink must consume a finite-length signal. To store the five second signal in the above example to "example.wav" we could write the following.

```julia
sound1 |> sink("example.wav",samplerate=44.1kHz)
```

In this case, since `sound1` had no defined sample rate, we have to tell sink what the sample rate is. If the signal already had a sample rate defined, we could just call `sink("example.wav")`.

## Signal promotion

A final concept, which may not be as obvious from the examples, is the use
of automatic signal promotion. When multiple signals are passed to the same
operator, and they have a different element type (e.g. `Float32` vs
`Float64`), different number of channels, or different sample rate, the signals
are first converted to the highest fidelity format and then operated on. This allows for a relatively seamless chain of operations where you don't have to worry about the specific format of the signal, and you won't loose information about your signals unless you explicitly request a lower fidelity signal format
(e.g. using [`tochannels`](@ref) or [`tosamplerate`](@ref)).

# Signal generation

There are three basic types that can be interpreted as signals: numbers, arrays and functions. Internally the function [`signal`](@ref) is called on any object passed to a signal operator, you can call this function yourself if you want to specify the exact sample rate by which you want to interpret the signal.

## Numbers

A number is treated as an infinite length signal, with unknown sample rate. 

```julia
1 |> until(1s) |> sink(samplerate=10Hz) == ones(10)
```

## Arrays

A standard array is treated as a finite signal, with unknown sample rate.

```julia
rand(10,2) |> sink(samplerate=10Hz) |> duration == 1
```

An `AxisArray` is treated as a finite signal with a known sample rate (and is the default output of [`sink`](@ref))

```julia
using AxisArrays
x = AxisArray(rand(10,1),Axis{:time}(range(0,1,length=10)))
samplerate(x) == 10
```

## Functions

A function is can be treated as an infinite signal. It should take a single argument which is the time in radians (*not* seconds). Radians are used to make defining the frequency of the signal for typical functions (e.g. sin), which is defind using the keywored argument `ω` or `frequency`. See [`signal`](@ref)'s documentation for more detail's.

```julia
signal(sin,ω=1kHz) |> duration |> isinf == true
```

An exception to this is `randn`. It can be used directly as a signal with unknown sample rate. 

```julia
randn |> duration == isinf
```

# Signal inspection

You can examine the properties of a signal using [`nsamples`](@ref), [`nchannels`](@ref), [`samplerate`](@ref), and [`duration`](@ref).

# Signal manipulation

There are several categories of signal manipulation: extending, cutting, filtering, ramping, and mapping.

## Extending

You can extend a signal using [`pad`](@ref) or [`append`](@ref). A padded signal becomes infinite and ends with the specified value, usually `one` or `zero`. You can append to or more signals (or [`preprend`](@ref)) so the occur one after another.

```julia
pad(x,zero) |> duration |> isinf == true
append(x,y) |> duration == duration(x) + duration(y)
```

## Cutting
You can cut signals apart, removing either the end of the signal ([`until`](@ref)) or the beginning ([`after`](@ref)). The operations are exact compliments of one another.

```julia
append(until(x,2s),after(x,2s)) |> nsamples == nsamples(x)
```

## Filtering
You can filter signals, removing undesired frequencies using [`lowpass`](@ref), [`highpass`](@ref), [`bandpass`](@ref), [`bandstop`](@ref) and [`filtersignal`](@ref). The latter allows the use of any arbitrary filter defined using `DSP`. 

```julia
signal(randn) |> lowpass(20Hz)
```

Note that if you use `using DSP` you will have to also call `dB = SignalOperators.Units.dB` if you want to make use of the proper meaning of `dB` for `SignalOperators`: `DSP` also defined a value for `dB`.

A unusual filter is [`normpower`](@ref): it computes the root mean squared power of the signal and then normalizes each sample by that value.

## Ramping
A ramp allows for smooth transition from 0 amplitude to the full amplitude of the signal. It is useful for avoid clicks in the onset or offset of a sound. For example, pure-tones are typically ramped when presented.

```julia
signal(sin,ω=2kHz) |> until(5s) |> ramp
```

You can ramp only the start of a signal ([`rampon`](@ref)), or the end of it ([`rampoff`](@ref)) and you can use ramps to create a smooth transition between two signals ([`fadeto`](@ref)). 

## Mapping

Probably the most powerful operator is [`mapsignal`](@ref). It works a lot like `map` but automatically promotes the signals, as with all operators, *and* it pads the end of the signal appropriately, so different length signals can be combined. The output is always the length of the longest *finite*-length signal.

```julia
a = signal(sin,ω=2kHz) |> until(2s)
b = signal(sin,ω=1kHz) |> until(3s)
a_minus_b = mapsignal(-,a,b)
```

The function `mapsignal` cannot, itself be piped, due to ambigity in the arugments, but shortcuts for this function have been provided for addition [`mix`](@ref)] and multiplication [`amplify`](@ref)], the two most common operations, and these two shortcuts have piped versions available.

```julia
a_plus_b = mix(a,b)
a_times_b = amplify(a,b)
```

You can add or select out channels using [`addchannel`](@ref) and [`channel`](@ref), which are defined in terms of calls to [`mapsignal`](@ref). These use a variant of `mapsignal` where the keyword `across_channels` is defined as true (see `mapsignal`'s documentation for details).