using .LibSndFile

for fmt in LibSndFile.supported_formats
    if fmt != DataFormat{:WAV}
    @eval function load_signal(::$fmt,filename,fs=missing)
        signal(load(filename),fs)
    end
    @eval function save_signal(::$fmt,filename,x,len)
        data = sink(x,duration=len*samples)
        save(filename,data)
    end
end