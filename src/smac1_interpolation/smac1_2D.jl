global const sunm = 1.989e33
global const pc   = 3.0856780e18

"""
   cic_mapping_2D( Pos, HSML, 
                  M, Rho, 
                  Bin_Q, Weights;
                  param::mappingParameters, kernel::AbstractSPHKernel,
                  show_progress::Bool=false )

Underlying function to map SPH data to a 2D grid.
"""
function smac1_mapping_2D( Pos, HSML, 
                        M, Rho, 
                        Bin_Q, Weights;
                        param::mappingParameters, kernel::AbstractSPHKernel,
                        show_progress::Bool=false )

    N = size(M,1)  # number of particles
    
    # max number of pixels over which the particle can be distributed
    N_distr = param.Npixels[1] * param.Npixels[2]


    wi_tab = zeros(N_distr)
    d1_tab = zeros(N_distr)
    d2_tab = zeros(N_distr)

    image = zeros(Float64, N_distr, 2)

    if param.periodic
        k_start = 0
    else
        k_start = 7
    end

    if show_progress
        P = Progress(N)
    end

    # loop over all particles
    @inbounds for p = 1:N

        bin_q = Float64(Bin_Q[p])

        if bin_q == 0.0
            continue
        end

        _pos = Pos[:,p]
        weight = Weights[p]
        hsml = HSML[p]
        hsml_inv = 1/hsml
        rho = Rho[p] / sunm * pc^3
        dbin = param.pixelSideLength

        cen = param.center

        x1 = _pos[1]-cen[1] 
        y1 = _pos[2]-cen[2]

        xp1 = (x1-hsml) / dbin
        xp2 = (x1+hsml) / dbin
        yp1 = (y1-hsml) / dbin
        yp2 = (y1+hsml) / dbin

        ip1 = max(floor(Integer, xp1)+1, 1)
        ip2 = min(floor(Integer, xp2)+1, param.Npixels[1])
        jp1 = max(floor(Integer, yp1)+1, 1)
        jp2 = min(floor(Integer, yp2)+1, param.Npixels[2])

        area = (xp2-xp1)*(yp2-yp1) * dbin^2
        d3   = (M[p]/rho) / area / dbin

        wit1    = 0.
        wit2    = 0.
        wit2tot = 0.
        witd    = 0.
        jcount  = 1

        @inbounds for ii = ip1:ip2, jj = jp1:jp2
            dimin = max(xp1, ii-1.0)
            djmin = max(yp1, jj-1.0)

            dimax = max(xp2, Real(ii))
            djmax = max(yp2, Real(jj))

            d1 = dimax - dimin # [pix]
            d2 = djmax - djmin  

            xx = x1 - (ii - 0.5)*dbin
            yy = y1 - (jj - 0.5)*dbin

            dist = √(xx^2 + yy^2) * hsml_inv

            wi = 𝒲(kernel, dist, 1.0)

            wi_tab[jcount] = wi
            d1_tab[jcount] = d1
            d2_tab[jcount] = d2

            jcount += 1

            wit1 += wi * d1 * d2
            witd += d1 * d2

            if wi > 0
                wit2 += 1
            end

            wit2tot += 1
        end # end loop over pixels

        fak = 1.0/(ip2-ip1)/(jp1-jp2)

        if wit1 > 0
            fak = wit2/wit1
        else
            wit2 = wit2tot
            fak = wit2tot / witd
        end

        fak_hsml=area/(wit2*dbin*dbin)

        jcount=1
        @inbounds for ii = ip1:ip2, jj = jp1:jp2

            wi     = wi_tab[jcount]
            d1     = d1_tab[jcount]
            d2     = d2_tab[jcount]
            jcount += 1

            if wit1 <= 0
                wi = 1.0
            end
            wi = wi * fak * fak_hsml

            idx = calculate_index(ii, jj, param.Npixels[1])

            image[idx,1] += bin_q * weight * d1 * d2 * d3 * wi
            image[idx,2] += weight * d1 * d2 * d3 * wi  
        end

         # update for ProgressMeter
        if show_progress
            next!(P)
        end
    end # p

    return image

end # function