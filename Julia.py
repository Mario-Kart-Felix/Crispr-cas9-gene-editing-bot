Make transcript by sample
using CSV

using DataFrames

tp_tr_sa = DataFrame()

for di in readdir(joinpath("../output/", ex, "psuedoalign/"), join = true)

    if !occursin("DS_Store", di)

        pa = string(joinpath(di, "abundance.tsv"))

        tpm = DataFrame(CSV.File(pa, delim = "	"))[:, [:target_id, :tpm]]

        sa = last(splitdir(di))

        tpm = rename!(tpm, :tpm => sa)

        if isempty(tp_tr_sa)

            tp_tr_sa = tpm

        else

            println(ncol(tp_tr_sa))

            id = ncol(tp_tr_sa) + 1

            insertcols!(tp_tr_sa, id, sa => tpm[:, sa])

        end

    end

end

tp_tr_sa = rename!(tp_tr_sa, :target_id => :id)

println(first(tp_tr_sa))

CSV.write(joinpath("../output/", ex, "transcript_x_sample.tsv"), tp_tr_sa)
Make gene by sample
using Statistics

tr_ge = DataFrame(CSV.File(pat, delim = "	"))

tr_ge = rename!(tr_ge, Dict("Transcript stable ID version" => :id, "Gene name" => :gene))


# Map transcript to gene name

tp_trge_sa = sort!(innerjoin(tp_tr_sa, tr_ge, on = :id), :gene)

tp_ge__sa = select!(tp_trge_sa, [n for n in names(tp_trge_sa) if n != "id"])


# Save the mean tpm for each gene

gr = groupby(tp_ge__sa, :gene)

sa_ = [n for n in names(tp_trge_sa) if n != "gene"]

tp_ge_sa = DataFrame()

for sa in sa_

    ge_sa = combine(gr, sa => sum)

    if isempty(tp_ge_sa)

        append!(tp_ge_sa, ge_sa)

    else

        tp_ge_sa = innerjoin(tp_ge_sa, ge_sa, on = :gene)

    end

end

println(tp_ge_sa[1:5, :])


# Save gene by sample

CSV.write(joinpath("../output/"
