    subroutine distaz(sta, sto, epa, epo, dk, dd, daze, dazs)
    real*8 :: sta, sto, epa, epo
    real*8 :: dk, dd, daze, dazs
    
    real*8 :: pi, rad
    real*8 :: sa, ea, ssa, csa, so, eo, sea, cea, ces, ses
    
    pi = 4.0*atan(1.0)
    rad = pi/180.0
    
    sa = atan(0.993270*tan(sta*rad))
    ea = atan(0.993270*tan(epa*rad))
    ssa = sin(sa)
    csa = cos(sa)
    so = sto*rad
    eo = epo*rad
    sea = sin(ea)
    cea = cos(ea)
    ces = cos(eo-so)
    ses = sin(eo-so)
    
    if (abs(sa-ea) < 1e-9) then
        if (abs(sto-epo) < 1e-9) then
            dk = 0.0
            dd = 0.0
            daze = 0.0
            dazs = 0.0
            return
        end if
    end if
    
    if (sta == 90.0) then
        if (epa == 90.0) then
            dk = 0.0
            dd = 0.0
            daze =0.0
            dazs = 0.0
            return
        end if
    end if
    
    if (sta == -90.0) then
        if (epa == -90.0) then
            dk = 0.0
            dd = 0.0
            daze =0.0
            dazs = 0.0
            return
        end if
    end if
    
    dd = ssa*sea + csa*cea*ces
    if (dd /= 0.0) dd = atan(sqrt(1.0-dd*dd)/dd)
    if (dd == 0.0) dd = pi/2.0
    if (dd < 0.0) dd = dd + pi
    dd = dd/rad
    dk = dd*111.19
    
    dazs = atan2(-ses,(ssa/csa*cea-sea*ces))
    daze = atan2(ses,(sea*csa/cea-ssa*ces))
    dazs = dazs/rad
    daze = daze/rad
    if (dazs < 0.0) dazs = dazs + 360.0
    if (daze < 0.0) daze = daze + 360.0
    
    return
    end subroutine distaz