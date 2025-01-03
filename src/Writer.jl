function write2csv(file_prefix)
    println("writing " * file_prefix * ".csv")
    open(file_prefix * ".csv", w) do file
    end
end