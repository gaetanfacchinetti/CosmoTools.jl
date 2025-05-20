push!(LOAD_PATH,"../src/")

using Documenter, CosmoTools

makedocs(sitename="CosmoTools.jl",
 pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "Public API" => ["background.md", "transfer.md"],
    ],)