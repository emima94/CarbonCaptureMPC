# Restore toml to only Jump related stuff and util packages (no SciML fluff)

using Pkg
Pkg.activate(".")

Pkg.add(["JuMP", 
        "Ipopt", 
        "CSV",
        "DataFrames", 
        "JSON", 
        "Revise", 
        "Plots", 
        "LaTeXStrings", 
        "Measures"])