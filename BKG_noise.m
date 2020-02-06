% adds noise that is reciprucal to distance from center
function BKG = BKG_noise(options)
    Sz = options.gridSz;
    BKG = zeros(Sz);
    for ik = 1:Sz
        for jk = 1:Sz
            if sqrt((ik - Sz/2).^2 + (jk - Sz/2).^2) > options.BeamStopSz
                BKG(ik, jk) = options.BKG.A.*(sqrt((ik - Sz/2).^2 + (jk - Sz/2).^2)).^options.BKG.r;
            end
        end
    end
end