include("continue.jl")
# Refine the changing length-six point before increasing chi; keep its source.
continue_cylinder("results/zigzag_L6_w1_chi512_star_projected0_seed7112.jls",512;cutoff=1e-12,noise_first=[1e-8,1e-9,0])
for L in [4,6,10]
    path=L==6 ? "results/zigzag_L6_w1_chi512_star_seed7112.jls" : "results/zigzag_L$(L)_w1_chi512_star_projected0_seed7112.jls"
    continue_cylinder(path,1024;cutoff=1e-12,noise_first=[1e-8,1e-9,0])
end
