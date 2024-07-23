using FFTW
using Plots
using DSP 
using LinearAlgebra
using ProgressMeter

plotly()

folder = "D:/Marlone/Profilometrie/120424/"
h_total1 = load( joinpath(folder, "h_x750.jld2"), "h_profile")

h = 1e-3
E = 70e3
nu = 0
g = 9.81
rho = 1e3
T_ = 9
B = h^3 * E / (12 * (1 - nu^2))

lambda_ = 5e-3:1e-3:0.6
k_th = 2 * π ./ lambda_

function relation(T, k, h, g, B, rho)
    return sqrt.(g * k + T .* k.^3 ./ rho + B .* k.^5 ./ rho) ./ (2 .* π)
end

f_th = relation(T_, k_th, h, g, B, rho)
f_th0 = relation(0, k_th, h, g, B, rho)

h_bis = diff(h_total1, dims=1)

Nx, Ny = size(h_bis)

window_x = hanning(Nx)
window_y = hanning(Ny)
window = window_x .* window_y'

faq = 120 # acquisition frequency
pix_par_m = 8 / 3e-3

lmin = 1 / pix_par_m
kmax = 2 * π / lmin

dk = kmax / Ny

f = (-Int(floor(Nx / 2)):Int(floor(Nx / 2))) .* (faq / Nx)
kabs = (-Int(floor(Ny / 2)):Int(floor(Ny / 2 - 1))) .* dk

Kabs = repeat(kabs, 1, length(f))'
F = repeat(f', length(kabs), 1)'

# Apply window to the FFT result
h_totalw = h_bis .* window

fft2 = fft(h_totalw)

plotfft = log.(abs.(fftshift(fft2)))

x1, x2 = 625, 665


plot(heatmap(x=Kabs[:, x1:x2], y=F[:, x1:x2], z=plotfft[:, x1:x2]))

plot!(-k_th, f_th, line=:dash, color=:white)
plot!(-k_th, f_th0, color=:white)

#region 

x_mid = 645

plotfft_tl = plotfft[11998:end, x1:x_mid]
plotfft_tr = plotfft[11998:end, x_mid:x2]
plotfft_bl = plotfft[1:11998, x1:x_mid]
plotfft_br = plotfft[1:11998, x_mid:x2]

plotfft_t = plotfft_tl[:, end:-1:1] .+ plotfft_tr
plotfft_b = plotfft_bl[:, end:-1:1] .+ plotfft_br

plotfft_total = (plotfft_tl[2:end, end:-1:1] .+ plotfft_br[end:-1:1, :]) / 2

heatmap(Kabs[11999:end, x_mid:x2], F[11999:end, x_mid:x2], plotfft_total, 
        aspect_ratio=:equal, color=:viridis, 
        xlims=(12.2, 250), ylims=(0, 60), 
        xlabel="k [m⁻¹]", ylabel="f [Hz]", 
        title="T=$T_ N/m",
        clims=(0.5, 2.0),
        c=:both)

plot!(k_th, f_th, line=:dash, color=:white)
plot!(k_th, f_th0, color=:white)

#endregion
