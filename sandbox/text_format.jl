using Plots
using LaTeXStrings
using Printf

lambda_F = 123.4
number_str = @sprintf("%.0e", lambda_F)
title = L"\mathrm{efwe: }\lambda_F = %$number_str"

#title = L"\mathrm{Opt.\ control},\ \lambda_F = %$(@sprintf("%.1e", lambda_F))"


#title = L"\lambda_F = %$(@sprintf("%.1e", lambda_F))"

plt = plot(1:10, rand(10), title=title)
savefig(plt, "sandbox/test_plot.pdf")