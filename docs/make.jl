using Documenter
using FT_Observer

makedocs(
    sitename = "FT_Observer.jl",
    modules = [FT_Observer],
    pages = [
        "Home" => "index.md",
        "Feedback 1" => "feedback1.md",
        "Feedback 2" => "feedback2.md",
        "Feedback 3" => "feedback3.md",
    ],
)
deploydocs(
    repo = "github.com/AnasXbouali/FT_Observer.git",
    devbranch = "main",
)
