% generates multiple PSF with erf profile
function PSF = MultiplePSFGenerator(options)

PSF = cell(options.PSF.number, 1);
sq = sqrt(options.PSF.number);
for jk = 1:sq
    for ik = 1:sq
        PSF{sq*(jk-1)+ik} = generatePSF(options.PSF.params(ik, :), ...
            options.PSF.params(jk, :), options.PSF.sz);
    end
end

end