# -*- coding: utf-8 -*-
"""
autovar_core.py
===============
Python port (2D, no drill-hole id, no anchors / geolocated variogram, no
"attempt to improve variogram", general anisotropy) of the Fortran program
AutoVar v1.003 (Luis Davila, 2022-2024), written for the GeostatisticsLab
QGIS plugin.  Only numpy is required.

The routines keep the names and the order of operations of the Fortran
sources so that the results can be compared one to one:

    datspac_horz  <- datspac_xy.f90
    omnipair      <- omni.f90
    moi / jacobi  <- moi.f90
    dblemodsortem <- sortem.f90
    _corrc        <- expvar.f90 (corrc / varexp, general variogram)
    varmodelfit   <- varmodel.f90 (random-perturbation fit, MT19937, seed 53783)

Behaviours that differ on purpose from the Fortran code are marked "DIFF";
they are summarised in AUTOVAR_DIFFERENCES at the end of this file.
"""

import math
import random as _pyrandom

import numpy as np

AUTOVAR_CORE_VERSION = 2   # 2: corrected algorithm + `legacy` mode

F32 = np.float32
DEG2RAD = math.pi / 180.0
RAD2DEG = 180.0 / math.pi
EPS_MAIN = float(F32(1.0e-6))     # EPSLON of Autovar.f90 / expvar.f90 (REAL*4)
EPS_DAT = float(F32(1.0e-2))      # EPSLON of datspac_xy.f90 (REAL*4)
SMALLDBLE = 1.0e-6


class AutoVarError(Exception):
    pass


class AutoVarCancelled(Exception):
    pass


# ---------------------------------------------------------------------------
# sortem.f90 :: dblemodsortem  (Singleton quicksort, same tie-breaking)
# ---------------------------------------------------------------------------
def dblemodsortem(a, carried=()):
    """Sort `a` ascending applying the same permutation to the arrays in
    `carried`.  Returns (a_sorted, [carried_sorted, ...])."""
    n = len(a)
    A = [0.0] + [float(t) for t in a]
    C = [[0.0] + [float(t) for t in c] for c in carried]
    lt = [0] * 80
    ut = [0] * 80
    j = n
    m = 1
    i = 1
    while True:
        d = j - i - 1
        if d >= 0:
            if d == 0:                                   # label 90
                if A[i] > A[j]:
                    A[i], A[j] = A[j], A[i]
                    for c in C:
                        c[i], c[j] = c[j], c[i]
            else:                                        # label 15
                p = (j + i) // 2
                ta = A[p]
                A[p] = A[i]
                tc = []
                for c in C:
                    tc.append(c[p])
                    c[p] = c[i]
                q = j
                k = i
                while True:                              # label 20
                    k += 1
                    if k > q:
                        break                            # -> 60
                    if A[k] <= ta:
                        continue
                    found = False
                    while True:                          # label 30
                        if A[q] < ta:
                            found = True                 # -> 40
                            break
                        q -= 1
                        if q > k:
                            continue
                        break                            # -> 50
                    if not found:
                        q = k - 1                        # label 50
                        break
                    A[k], A[q] = A[q], A[k]              # label 40
                    for c in C:
                        c[k], c[q] = c[q], c[k]
                    q -= 1
                A[i] = A[q]                              # label 60
                A[q] = ta
                for c, t in zip(C, tc):
                    c[i] = c[q]
                    c[q] = t
                if 2 * q <= i + j:                       # label 70
                    lt[m] = q + 1
                    ut[m] = j
                    j = q - 1
                else:
                    lt[m] = i
                    ut[m] = q - 1
                    i = q + 1
                m += 1
                continue                                 # label 10
        m -= 1                                           # label 100
        if m <= 0:
            break
        i = lt[m]
        j = ut[m]
    return np.array(A[1:]), [np.array(c[1:]) for c in C]


# ---------------------------------------------------------------------------
# ACORN generator used by datspac_horz for de-spiking (common block /iaco/)
# ---------------------------------------------------------------------------
class _Acorn:
    MAXINT = 2 ** 30

    def __init__(self):
        self.ixv = [0] * 13            # static common block -> zeros

    def next(self):
        ixv = self.ixv
        for i in range(12):
            ixv[i + 1] = ixv[i + 1] + ixv[i]
            if ixv[i + 1] >= self.MAXINT:
                ixv[i + 1] -= self.MAXINT
        return ixv[12] / self.MAXINT


# ---------------------------------------------------------------------------
# datspac_xy.f90 :: datspac_horz
# ---------------------------------------------------------------------------
def datspac_horz(xyz, acorn, ndt=None, legacy=False):
    """Average of: minimum / average / maximum distance to the (up to) 4 nearest
    neighbours.  DIFF: in the Fortran source `ndt` is uninitialised when there is
    no drill-hole id column; here ndt = n (unless given)."""
    xyz = np.asarray(xyz, dtype=float)
    n = len(xyz)
    if ndt is None:
        ndt = n
    minspac2 = 1.0
    maxspac = math.sqrt((xyz[:, 0].max() - xyz[:, 0].min()) ** 2 +
                        (xyz[:, 1].max() - xyz[:, 1].min()) ** 2 +
                        (xyz[:, 2].max() - xyz[:, 2].min()) ** 2)
    maxspac2 = maxspac * maxspac
    ndmax = 4

    acorn.ixv[0] = 69069
    for _ in range(1000):
        acorn.next()

    xloc = xyz[:, 0].copy()
    yloc = xyz[:, 1].copy()
    if legacy:                         # quirk of the original code (removed in the corrected version)
        if xloc[0] < 0:
            xloc[0] = xyz[:, 0].min()
        if yloc[0] < 0:
            yloc[0] = xyz[:, 1].min()

    distmin = np.zeros(n)
    avgdist = np.zeros(n)
    distmax = np.zeros(n)
    tmpdist = np.zeros(n)              # persistent: stale value read at tmpdist(na+1)
    xs = xloc[:ndt]
    ys = yloc[:ndt]
    for idat in range(ndt):
        dx = xs - xs[idat]
        dy = ys - ys[idat]
        h2 = dx * dx + dy * dy
        ok = (h2 <= maxspac2) & (h2 >= minspac2)
        ok[idat] = False
        hh = np.sqrt(h2[ok])
        na = hh.size
        if na == 0:
            continue
        tmpdist[:na] = np.sort(hh)
        if legacy:                     # stale tmpdist(na+1) of the original code
            keep = (tmpdist[1:na + 1] - tmpdist[:na]) >= 1.0
        else:                          # the last distance has no successor: always kept
            keep = np.ones(na, bool)
            keep[:-1] = (tmpdist[1:na] - tmpdist[:na - 1]) >= 1.0
        sel = tmpdist[:na][keep][:ndmax]
        nb = sel.size
        if nb == 0:
            continue
        distmin[idat] = sel[0]
        distmax[idat] = sel[nb - 1]
        avgdist[idat] = sel.sum() / nb

    distmin[:ndt] = np.sort(distmin[:ndt])
    avgdist[:ndt] = np.sort(avgdist[:ndt])
    distmax[:ndt] = np.sort(distmax[:ndt])
    for i in range(ndt - 1, 0, -1):                  # de-spike ties
        if distmin[i - 1] == distmin[i]:
            rnd = acorn.next()
            distmin[i] = distmin[i] + rnd * EPS_DAT if rnd > 0.5 else distmin[i] - rnd * EPS_DAT
        if avgdist[i - 1] == avgdist[i]:
            rnd = acorn.next()
            avgdist[i] = avgdist[i] + rnd * EPS_DAT if rnd > 0.5 else avgdist[i] - rnd * EPS_DAT
        if distmax[i - 1] == distmax[i]:
            rnd = acorn.next()
            distmax[i] = distmax[i] + rnd * EPS_DAT if rnd > 0.5 else distmax[i] - rnd * EPS_DAT
    return (distmin[:ndt].sum() / ndt, avgdist[:ndt].sum() / ndt,
            distmax[:ndt].sum() / ndt)


# ---------------------------------------------------------------------------
# general.f90 :: get_max_dist (2D)
# ---------------------------------------------------------------------------
def get_max_dist(xyz, block=512):
    P = np.asarray(xyz, float)[:, :2]
    best = 0.0
    for s in range(0, len(P), block):
        B = P[s:s + block]
        d2 = ((B[:, None, :] - P[None, :, :]) ** 2).sum(axis=2)
        best = max(best, float(d2.max()))
    return math.sqrt(best)


# ---------------------------------------------------------------------------
# omni.f90 :: omnipair
# ---------------------------------------------------------------------------
def omnipair(x, y, z, v, maxdist, legacy=False):
    """All ordered pairs (i tail, j head) within 0.2*maxdist, with the
    pseudo-correlation (v_i v_j - m^2)/s^2."""
    n = len(x)
    srch = 0.2 * maxdist
    av = float(np.sum(v)) / max(n, 1)
    ss = float(np.sum(v * v)) / max(n, 1) - av * av
    if ss <= 0.0:
        raise AutoVarError("The variable has zero variance.")
    idx = np.arange(1, n + 1, dtype=float)
    xs, (ys, zs, ids) = dblemodsortem(x, [y, z, idx])
    ids = ids.astype(int) - 1
    vs = v[ids]
    last = n - 1 if legacy else n     # original code: j runs 1..n-1 (last sorted point never a head/tail)
    O_h, O_dx, O_dy, O_dz, O_sq = [], [], [], [], []
    for i in range(n - 1 if legacy else n):
        lo = int(np.searchsorted(xs, xs[i] - srch, side='left'))
        hi = min(int(np.searchsorted(xs, xs[i] + srch, side='right')), last)
        if hi <= lo:
            continue
        j = np.arange(lo, hi)
        j = j[j != i]
        if j.size == 0:
            continue
        h = np.sqrt((xs[i] - xs[j]) ** 2 + (ys[i] - ys[j]) ** 2 + (zs[i] - zs[j]) ** 2)
        m = h <= srch
        if not m.any():
            continue
        j = j[m]
        O_h.append(h[m])
        O_dx.append(xs[j] - xs[i])
        O_dy.append(ys[j] - ys[i])
        O_dz.append(zs[j] - zs[i])
        O_sq.append((vs[i] * vs[j] - av ** 2) / ss)
    if not O_h:
        raise AutoVarError("No pairs inside the search radius.")
    return (np.concatenate(O_h), np.concatenate(O_dx), np.concatenate(O_dy),
            np.concatenate(O_dz), np.concatenate(O_sq))


# ---------------------------------------------------------------------------
# moi.f90 :: jacobi (Numerical Recipes; internal scalars are REAL*4)
# ---------------------------------------------------------------------------
def jacobi(a_in):
    a = np.array(a_in, dtype=float)
    n = 3
    v = np.eye(3)
    b = [F32(a[i, i]) for i in range(n)]
    d = [float(b[i]) for i in range(n)]
    z = [F32(0.0)] * n
    one = F32(1.0)
    nrot = 0
    for it in range(1, 101):
        sm = F32(0.0)
        for ip in range(n - 1):
            for iq in range(ip + 1, n):
                sm = F32(float(sm) + abs(a[ip, iq]))
        if sm == 0.0:
            return np.array(d), v, nrot
        tresh = F32(F32(F32(0.2) * sm) / F32(n * n)) if it < 4 else F32(0.0)
        for ip in range(n - 1):
            for iq in range(ip + 1, n):
                g = F32(100.0 * abs(a[ip, iq]))
                if (it > 4 and abs(d[ip]) + float(g) == abs(d[ip]) and
                        abs(d[iq]) + float(g) == abs(d[iq])):
                    a[ip, iq] = 0.0
                elif abs(a[ip, iq]) > float(tresh):
                    h = F32(d[iq] - d[ip])
                    if F32(abs(h) + g) == abs(h):
                        t = F32(a[ip, iq] / float(h))
                    else:
                        theta = F32(float(F32(F32(0.5) * h)) / a[ip, iq])
                        t = F32(one / F32(abs(theta) + np.sqrt(F32(one + F32(theta * theta)))))
                        if theta < 0:
                            t = F32(-t)
                    c = F32(one / np.sqrt(F32(one + F32(t * t))))
                    s = F32(t * c)
                    tau = F32(s / F32(one + c))
                    h = F32(float(t) * a[ip, iq])
                    z[ip] = F32(z[ip] - h)
                    z[iq] = F32(z[iq] + h)
                    d[ip] = d[ip] - float(h)
                    d[iq] = d[iq] + float(h)
                    a[ip, iq] = 0.0

                    def rot(M, r1, c1, r2, c2):
                        gg = F32(M[r1, c1])
                        hh = F32(M[r2, c2])
                        M[r1, c1] = float(F32(gg - F32(s * F32(hh + F32(gg * tau)))))
                        M[r2, c2] = float(F32(hh + F32(s * F32(gg - F32(hh * tau)))))
                    for j in range(0, ip):
                        rot(a, j, ip, j, iq)
                    for j in range(ip + 1, iq):
                        rot(a, ip, j, j, iq)
                    for j in range(iq + 1, n):
                        rot(a, ip, j, iq, j)
                    for j in range(n):
                        rot(v, j, ip, j, iq)
                    nrot += 1
        for ip in range(n):
            b[ip] = F32(b[ip] + z[ip])
            d[ip] = float(b[ip])
            z[ip] = F32(0.0)
    raise AutoVarError("Too many iterations in jacobi")


def _jacobi_double(a_in):
    """Jacobi eigen-decomposition in double precision (same algorithm, REAL*8)."""
    a = np.array(a_in, dtype=float)
    n = 3
    v = np.eye(3)
    b = [a[i, i] for i in range(n)]
    d = list(b)
    z = [0.0] * n
    for it in range(1, 101):
        sm = sum(abs(a[ip, iq]) for ip in range(n - 1) for iq in range(ip + 1, n))
        if sm == 0.0:
            return np.array(d), v, 0
        tresh = 0.2 * sm / (n * n) if it < 4 else 0.0
        for ip in range(n - 1):
            for iq in range(ip + 1, n):
                g = 100.0 * abs(a[ip, iq])
                if it > 4 and abs(d[ip]) + g == abs(d[ip]) and abs(d[iq]) + g == abs(d[iq]):
                    a[ip, iq] = 0.0
                elif abs(a[ip, iq]) > tresh:
                    h = d[iq] - d[ip]
                    if abs(h) + g == abs(h):
                        t = a[ip, iq] / h
                    else:
                        theta = 0.5 * h / a[ip, iq]
                        t = 1.0 / (abs(theta) + math.sqrt(1.0 + theta * theta))
                        if theta < 0:
                            t = -t
                    c = 1.0 / math.sqrt(1.0 + t * t)
                    s = t * c
                    tau = s / (1.0 + c)
                    h = t * a[ip, iq]
                    z[ip] -= h
                    z[iq] += h
                    d[ip] -= h
                    d[iq] += h
                    a[ip, iq] = 0.0

                    def rot(M, r1, c1, r2, c2):
                        gg, hh = M[r1, c1], M[r2, c2]
                        M[r1, c1] = gg - s * (hh + gg * tau)
                        M[r2, c2] = hh + s * (gg - hh * tau)
                    for j in range(0, ip):
                        rot(a, j, ip, j, iq)
                    for j in range(ip + 1, iq):
                        rot(a, ip, j, j, iq)
                    for j in range(iq + 1, n):
                        rot(a, ip, j, iq, j)
                    for j in range(n):
                        rot(v, j, ip, j, iq)
        for ip in range(n):
            b[ip] += z[ip]
            d[ip] = b[ip]
            z[ip] = 0.0
    raise AutoVarError("Too many iterations in jacobi")


# ---------------------------------------------------------------------------
# moi.f90 :: moi
# ---------------------------------------------------------------------------
def moi(sq, dist, dx, dy, dz, legacy=False):
    ss = len(sq)
    tmp = np.sort(sq)
    perc_max = int(F32(0.99) * F32(ss)) if legacy else int(0.99 * ss)
    if ss % 99 == 0:
        max_corr = tmp[perc_max - 1]
    else:
        max_corr = (tmp[perc_max - 1] + tmp[perc_max]) / 2.0
    nd = int(np.sum(sq <= max_corr))
    if nd < 2:
        raise AutoVarError("Too few data for calculating the MOI matrix.")
    # DIFF (Fortran bug reproduced): the program copies sqdfmoi(nd), dxmoi(nd), ...
    # (the FIRST nd elements) instead of the selected ones, and the sum runs to nd-1.
    if legacy:
        val = sq[:nd]
        k = slice(0, nd - 1)
        p = dist[:nd][k] ** 2.5
        ddx, ddy, ddz = dx[:nd][k], dy[:nd][k], dz[:nd][k]
        vv = val[k]
    else:                                             # pairs below the 99th percentile, all of them
        sel = sq <= max_corr
        p = dist[sel] ** 2.5
        ddx, ddy, ddz = dx[sel], dy[sel], dz[sel]
        vv = sq[sel]
    M = np.zeros((3, 3))
    M[1, 0] = -np.sum(vv * ddx * ddy / p)
    M[2, 0] = -np.sum(vv * ddx * ddz / p)
    M[2, 1] = -np.sum(vv * ddy * ddz / p)
    M[0, 0] = np.sum(vv * (ddy ** 2 + ddz ** 2) / p)
    M[1, 1] = np.sum(vv * (ddx ** 2 + ddz ** 2) / p)
    M[2, 2] = np.sum(vv * (ddx ** 2 + ddy ** 2) / p)
    M[0, 1], M[0, 2], M[1, 2] = M[1, 0], M[2, 0], M[2, 1]
    eigval, eigvec, _ = jacobi(M) if legacy else _jacobi_double(M)
    return M, eigvec, eigval, float(max_corr)


# ---------------------------------------------------------------------------
# expvar.f90 :: corrc / varexp (general variogram, disseminated, 2D)
# ---------------------------------------------------------------------------
def _corrc(P, v, lagsize, nlag, vspc, dir_, block=256):
    """P: (n,3) rotated coordinates (major, minor, third).  dir_ 1 = along the major
    axis (band on the minor coordinate), 2 = along the minor axis."""
    n = len(P)
    azmtol = 22.5 * DEG2RAD
    diptol = 22.5 * DEG2RAD
    azmtoldist = 0.5 * lagsize / math.tan(azmtol)
    diptoldist = 0.5 * vspc / math.tan(diptol)
    lagtol = 0.5 * lagsize
    ia, ib = (0, 1) if dir_ == 1 else (1, 0)
    L = nlag + 1
    num = np.zeros(L); dis = np.zeros(L); cov = np.zeros(L)
    mu = np.zeros(L); muh = np.zeros(L); su = np.zeros(L); suh = np.zeros(L)
    cols = np.arange(n)
    for s in range(0, n, block):
        e = min(s + block, n)
        I = np.arange(s, e)
        d_lat = P[None, :, ib] - P[I, ib][:, None]
        d_z = P[None, :, 2] - P[I, 2][:, None]
        d_a = np.abs(P[None, :, ia] - P[I, ia][:, None])
        ok = (np.abs(d_lat) <= lagsize) & (np.abs(d_z) <= vspc)       # disseminated
        ok &= (I[:, None] != cols[None, :])
        ang = np.abs(np.arctan(d_lat / np.maximum(d_a, SMALLDBLE)))
        ok &= ~((d_a < azmtoldist) & (ang > azmtol))
        angz = np.abs(np.arctan(d_z / np.maximum(d_a, SMALLDBLE)))
        ok &= ~((d_a < diptoldist) & (angz > diptol))
        ii, jj = np.nonzero(ok)
        if ii.size == 0:
            continue
        ig = ii + s
        h = np.sqrt((P[ig, 0] - P[jj, 0]) ** 2 + (P[ig, 1] - P[jj, 1]) ** 2 +
                    (P[ig, 2] - P[jj, 2]) ** 2)
        xd = d_a[ii, jj]
        vi = v[ig]
        vj = v[jj]
        l0 = np.floor(h / lagsize + 0.5).astype(int)
        for cand in (l0 - 1, l0, l0 + 1):      # a pair on the bin edge belongs to both
            m = ((cand >= 0) & (cand <= nlag) & (h >= cand * lagsize - lagtol) &
                 (h <= cand * lagsize + lagtol))
            if not m.any():
                continue
            c = cand[m]
            num += np.bincount(c, minlength=L)
            dis += np.bincount(c, weights=xd[m], minlength=L)
            cov += np.bincount(c, weights=(vi * vj)[m], minlength=L)
            mu += np.bincount(c, weights=vi[m], minlength=L)
            muh += np.bincount(c, weights=vj[m], minlength=L)
            su += np.bincount(c, weights=(vi * vi)[m], minlength=L)
            suh += np.bincount(c, weights=(vj * vj)[m], minlength=L)
    return num, dis, cov, mu, muh, su, suh


def _finish_lags(num, dis, cov, mu, muh, su, suh):
    """Pseudo-correlogram per lag -> 1 - rho(h); lags with < 5 pairs are dropped."""
    out_num, out_dis, out_val = [], [], []
    for l in range(len(num)):
        if num[l] < 5:
            continue
        n_ = num[l]
        c = cov[l] / n_ - (mu[l] / n_) * (muh[l] / n_)
        s1 = max(su[l] / n_ - (mu[l] / n_) ** 2, EPS_MAIN)
        s2 = max(suh[l] / n_ - (muh[l] / n_) ** 2, EPS_MAIN)
        rho = c / math.sqrt(s1 * s2)
        d = dis[l] / n_
        if d > 0.0:
            out_num.append(n_)
            out_dis.append(d)
            out_val.append(1.0 - rho)
    return np.array(out_num), np.array(out_dis), np.array(out_val)


# ---------------------------------------------------------------------------
# varmodel.f90 :: Mersenne Twister, fixsills, varmodelfit
# ---------------------------------------------------------------------------
def _init_genrand(seed):
    """init_genrand(seed) of random.f90 -> python Random whose getrandbits(32)
    returns the same tempered 32-bit words as the Fortran grnd()."""
    mt = [0] * 624
    mt[0] = seed & 0xFFFFFFFF
    latest = mt[0]
    for i in range(1, 624):
        latest = (latest ^ (latest >> 30)) & 0xFFFFFFFF
        latest = (latest * 1812433253 + i) & 0xFFFFFFFF
        mt[i] = latest
    r = _pyrandom.Random()
    r.setstate((3, tuple(mt) + (624,), None))
    return r


def _fixsills(cc, lo=0.01, hi=1.0e21, sill=1.0):
    """fixsills for: nugget fixed at 0, total sill fixed, structure sills in [lo, hi]."""
    rtol = 1.0e-6 * sill
    niter = 0
    while True:
        tot = 0.0
        for i in range(len(cc)):
            cc[i] = min(max(cc[i], lo), hi)
            tot += cc[i]
        if abs(tot - sill) < rtol:
            break
        if tot > 0.0:
            f = sill / tot
            for i in range(len(cc)):
                cc[i] *= f
        niter += 1
        if niter > 100:
            break
    for i in range(len(cc)):
        cc[i] = min(max(cc[i], lo), hi)
    return cc


def _sph(h, a):
    hr = h / max(a, SMALLDBLE)
    return np.where(hr < 1.0, 1.5 * hr - 0.5 * hr ** 3, 1.0)


def varmodelfit(lagdist, value, numpairs, wts_in, is_major, ahmax_b, ahmin_b,
                maxiter=2000000, rseed=53783, minpairs=10, progress=None, nst=3, legacy=False):
    """Random-perturbation fit: `nst` spherical structures, nugget 0, total sill 1,
    ranges in [bounds].  Inputs are the full variogram table (major lags first, then
    minor): lag distance, value, number of pairs, weight (0 = not used), is_major."""
    value1 = np.minimum(value, 1.0)                  # values above the sill are clipped to 1
    sel = wts_in > 0.0
    npoint = int(sel.sum())
    vl = lagdist[sel]
    vv = value1[sel]
    vmaj = is_major[sel]
    # DIFF (Fortran bug reproduced): the unfiltered array varnumpairs is indexed with the
    # filtered index; it only matters when some lag has weight 0.
    npairs_idx = numpairs[:npoint] if legacy else numpairs[sel]
    w = np.where(npairs_idx < max(minpairs, 1), 0.0, 1.0)
    if int((npairs_idx >= max(minpairs, 1)).sum()) < 2:
        raise AutoVarError("Too few variogram values to optimise the model.")
    cc = _fixsills([0.0] * nst)
    maxlag_maj = max([vl[i] for i in range(npoint) if vmaj[i]] + [SMALLDBLE])
    maxlag_min = max([vl[i] for i in range(npoint) if not vmaj[i]] + [SMALLDBLE])
    if maxlag_maj <= SMALLDBLE:
        maxlag_maj = float(vl.max())
    if maxlag_min <= SMALLDBLE:
        maxlag_min = float(vl.max())
    ahmax, ahmin = [], []
    for ist in range(1, nst + 1):
        frac = float(F32(ist) / F32(nst))
        ahmax.append(min(max(frac * maxlag_maj, ahmax_b[0]), ahmax_b[1]))
        ahmin.append(min(max(frac * maxlag_min, ahmin_b[0]), ahmin_b[1]))

    def mkrow(amax, amin):
        r = np.empty(npoint)
        r[vmaj] = _sph(vl[vmaj], amax)
        r[~vmaj] = _sph(vl[~vmaj], amin)
        return r

    S = np.vstack([mkrow(ahmax[k], ahmin[k]) for k in range(nst)])
    rowpar = [(ahmax[k], ahmin[k]) for k in range(nst)]

    def objective(c):
        r = np.dot(c, S) - vv
        return float(np.dot(w, r * r))

    old_obj = start_obj = objective(cc)
    old_cc, old_amax, old_amin = list(cc), list(ahmax), list(ahmin)
    old_sill = sum(cc)
    rng = _init_genrand(rseed)
    grb = rng.getrandbits
    inv = 1.0 / 4294967295.0
    lo_max, hi_max = ahmax_b
    lo_min, hi_min = ahmin_b
    report = max(maxiter // 400, 1)
    floor = math.floor
    for it in range(1, maxiter + 1):
        if grb(32) * inv < 0.5:                                  # perturb a sill
            ip = int(floor(grb(32) * inv * nst + 1))
            shift = 0.15 * (grb(32) * inv - 0.5) * old_sill
            cc[ip - 1] += shift
            _fixsills(cc)
        else:                                                    # perturb a range
            ist = int(floor(grb(32) * inv * nst + 1)) - 1
            ip = int(floor(grb(32) * inv * 3 + 1))
            shift = 0.15 * (grb(32) * inv - 0.5)
            if ip == 3:
                pass                                             # vertical range (0 in 2D)
            elif ip == 2:
                ahmin[ist] = ahmin[ist] + shift * max(ahmin)
                if ist < nst - 1 and ahmin[ist] > ahmin[ist + 1]:
                    ahmin[ist] = ahmin[ist + 1]
                ahmin[ist] = min(max(ahmin[ist], lo_min), hi_min)
            else:
                ahmax[ist] = ahmax[ist] + shift * max(ahmax)
                if ist < nst - 1 and ahmax[ist] > ahmax[ist + 1]:
                    ahmax[ist] = ahmax[ist + 1]
                ahmax[ist] = min(max(ahmax[ist], lo_max), hi_max)
            if ahmin[ist] > ahmax[ist]:
                ahmin[ist] = ahmax[ist]
            for k in range(nst):
                if ahmin[k] > ahmax[k]:
                    ahmin[k] = ahmax[k]
            for k in range(nst):
                if rowpar[k][0] != ahmax[k] or rowpar[k][1] != ahmin[k]:
                    S[k] = mkrow(ahmax[k], ahmin[k])
                    rowpar[k] = (ahmax[k], ahmin[k])
        obj = objective(cc)
        if obj < old_obj:
            old_obj = obj
            old_cc[:] = cc
            old_amax[:] = ahmax
            old_amin[:] = ahmin
            old_sill = sum(cc)
        else:
            cc[:] = old_cc
            if ahmax != old_amax or ahmin != old_amin:
                ahmax[:] = old_amax
                ahmin[:] = old_amin
                for k in range(nst):
                    if rowpar[k][0] != ahmax[k] or rowpar[k][1] != ahmin[k]:
                        S[k] = mkrow(ahmax[k], ahmin[k])
                        rowpar[k] = (ahmax[k], ahmin[k])
        if progress is not None and it % report == 0:
            progress(it, maxiter, old_obj)
    cc = list(old_cc)
    _fixsills(cc)
    return dict(cc=cc, ahmax=list(old_amax), ahmin=list(old_amin),
                start_obj=start_obj, final_obj=old_obj, npoint=npoint)


# ---------------------------------------------------------------------------
# Model evaluation helpers (what the plugin plots)
# ---------------------------------------------------------------------------
def model_gamma(h, cc, ranges, c0=0.0):
    """Sum of spherical structures (sill cc_i, range a_i) evaluated at distances h
    in one of the two principal directions (use ahmax or ahmin as `ranges`)."""
    h = np.asarray(h, float)
    g = np.zeros_like(h) + 0.0
    for c, a in zip(cc, ranges):
        g += c * _sph(h, a)
    return g + np.where(h > 0, c0, 0.0)


# ---------------------------------------------------------------------------
# Driver: Autovar.f90 (2D / disseminated / general variogram)
# ---------------------------------------------------------------------------
def run_autovar(x, y, v, maxiter=2000000, rseed=53783, progress=None, log=None,
                legacy_projection=None, ndt=None, float32_input=False, nst=3,
                lag_override=None, legacy=False):
    """Run AutoVar on the data (x, y, v).  Returns a dictionary with:
       exp_major / exp_minor : dict(dist, num, val)  experimental "variogram" (1-rho)
       model                 : dict(nst, cc, ahmax, ahmin, azimuth, c0)
       az_major, az_minor, angle_gslib, eigenvectors, eigenvalues, lags, objective ..."""
    if legacy_projection is None:
        legacy_projection = legacy
    def _log(msg):
        if log is not None:
            log(msg)

    x = np.asarray(x, float).copy()
    y = np.asarray(y, float).copy()
    v = np.asarray(v, float).copy()
    good = np.isfinite(x) & np.isfinite(y) & np.isfinite(v)
    x, y, v = x[good], y[good], v[good]
    if float32_input:                                # the Fortran program reads REAL*4
        x, y, v = (x.astype(np.float32).astype(float), y.astype(np.float32).astype(float),
                   v.astype(np.float32).astype(float))
    nd = len(x)
    if nd < 30:
        raise AutoVarError("At least 30 data are needed.")
    z = np.ones(nd)                                  # 2D: Z is a constant
    res = dict(nd=nd)
    res['mean'] = float(v.mean())
    res['variance'] = float(np.mean(v * v) - v.mean() ** 2)

    acorn = _Acorn()
    xyz = np.column_stack([x, y, z])
    avgmin, avgavg, avgmax = datspac_horz(xyz, acorn, ndt, legacy)
    res['spacing'] = (avgmin, avgavg, avgmax)
    _log("Data spacing (min, avg, max): %.3f  %.3f  %.3f" % (avgmin, avgavg, avgmax))

    # --- omnidirectional pairs inside 0.2 * maximum distance -----------------------------
    max_dist = get_max_dist(xyz)
    h, dxp, dyp, dzp, sq = omnipair(x, y, z, v, max_dist, legacy)
    total_pairs = len(h)
    sq = np.where(sq < 0.0, 0.0, sq)
    pos = sq > 0.0
    h, dxp, dyp, dzp, sq = h[pos], dxp[pos], dyp[pos], dzp[pos], sq[pos]
    res['npairs_omni'] = (total_pairs, int(pos.sum()))
    _log("Omnidirectional pairs: %d (positive pseudo-correlation: %d)" % res['npairs_omni'])

    # --- moment of inertia -> principal directions ---------------------------------------
    M, eigvec, eigval, p99 = moi(sq, h, dxp, dyp, dzp, legacy)
    _, (idxf,) = dblemodsortem(eigval.copy(), [np.array([1.0, 2.0, 3.0])])
    idx = idxf.astype(int) - 1
    e = np.zeros((3, 3))
    e[0, 0], e[1, 0], e[2, 0] = eigvec[0, idx[0]], eigvec[1, idx[0]], 0.0
    e[0, 1], e[1, 1], e[2, 1] = eigvec[0, idx[1]], eigvec[1, idx[1]], 0.0
    e[0, 2], e[1, 2], e[2, 2] = eigvec[0, idx[2]], eigvec[1, idx[2]], 1.0
    res['moi'] = M
    res['eigenvectors'] = e
    res['eigenvalues'] = eigval
    try:
        a1 = (math.pi / 2.0 - math.atan(e[1, 0] / e[0, 0])) * RAD2DEG
    except ZeroDivisionError:
        a1 = 0.0
    az_major = math.atan2(e[0, 0], e[1, 0]) * RAD2DEG
    az_minor = math.atan2(e[0, 1], e[1, 1]) * RAD2DEG
    res['angle_gslib'] = a1
    res['az_major'] = az_major
    res['az_minor'] = az_minor
    _log("Principal directions (azimuth): major %.4f  minor %.4f" % (az_major, az_minor))

    lrot = np.array([[e[0, 0], e[1, 0], e[2, 0]],
                     [e[0, 1], e[1, 1], e[2, 1]],
                     [e[0, 2], e[1, 2], e[2, 2]]])          # rows = eigenvectors

    # --- extents and lags (data are re-ordered by successive sorts, as in Autovar.f90) ---
    il = int(F32(1.0) + F32(0.04) * F32(nd))
    ih = nd - il
    xs, (ys, zs, vs) = dblemodsortem(x, [y, z, v])
    deltax = xs[ih - 1] - xs[il - 1]
    ys, (xs, zs, vs) = dblemodsortem(ys, [xs, zs, vs])
    deltay = ys[ih - 1] - ys[il - 1]
    zs, (xs, ys, vs) = dblemodsortem(zs, [xs, ys, vs])
    deltaz = zs[ih - 1] - zs[il - 1]
    usex = deltax > EPS_MAIN
    usey = deltay > EPS_MAIN
    usez = deltaz > EPS_MAIN
    dmajco = lrot[0, 0] * (xs if usex else 0.0) + lrot[0, 1] * (ys if usey else 0.0) + \
        lrot[0, 2] * (zs if usez else 0.0)
    dminco = lrot[1, 0] * (xs if usex else 0.0) + lrot[1, 1] * (ys if usey else 0.0) + \
        lrot[1, 2] * (zs if usez else 0.0)
    dterco = lrot[2, 0] * (xs if usex else 0.0) + lrot[2, 1] * (ys if usey else 0.0) + \
        lrot[2, 2] * (zs if usez else 0.0)
    dmajco = np.broadcast_to(dmajco, (nd,)).copy()
    dminco = np.broadcast_to(dminco, (nd,)).copy()
    dterco = np.broadcast_to(dterco, (nd,)).copy()
    xyz_anis = np.column_stack([dmajco, dminco, dterco])
    a_min, a_avg, a_max = datspac_horz(xyz_anis, acorn, ndt, legacy)
    dmajlag = max(EPS_MAIN, abs(a_avg))
    dminlag = max(EPS_MAIN, abs(a_min))
    dterlag = 1.0
    if lag_override is not None:             # diagnostic use only
        dmajlag, dminlag = lag_override
    max_maj = 0.5 * abs(dmajco.max() - dmajco.min())
    max_min = 0.5 * abs(dminco.max() - dminco.min())
    max_min = min(max_min, max_maj)
    nmaj_lag = int(max_maj / dmajlag) + 1
    nmin_lag = int(max_min / dminlag) + 1
    res['lags'] = dict(major=dmajlag, minor=dminlag, extent_major=max_maj, extent_minor=max_min,
                       n_major=nmaj_lag, n_minor=nmin_lag)
    _log("Lag size major / minor: %.3f / %.3f  (extent %.1f / %.1f)" %
         (dmajlag, dminlag, max_maj, max_min))

    # --- experimental variograms in the principal directions (original data order) -------
    if legacy_projection:        # xyzrot = matmul(xyz, lrotmat) of the Fortran source
        P1 = x * lrot[0, 0] + y * lrot[1, 0]
        P2 = x * lrot[0, 1] + y * lrot[1, 1]
    else:                        # true projection onto the eigenvectors
        P1 = x * lrot[0, 0] + y * lrot[0, 1]
        P2 = x * lrot[1, 0] + y * lrot[1, 1]
    P = np.column_stack([P1, P2, np.full(nd, 0.5)])
    vspc = 1.0
    maj = _finish_lags(*_corrc(P, v, dmajlag, nmaj_lag, vspc, 1))
    mnr = _finish_lags(*_corrc(P, v, dminlag, nmin_lag, vspc, 2))
    if len(maj[0]) < 3 or len(mnr[0]) < 3:
        raise AutoVarError("Too few lags with at least 5 pairs in a principal direction.")
    wmaj = np.ones(len(maj[0]))
    wmin = np.ones(len(mnr[0]))
    if maj[2][0] > maj[2][1] and maj[0][0] < maj[0][1]:
        wmaj[0] = 0.0
    if mnr[2][0] > mnr[2][1] and mnr[0][0] < mnr[0][1]:
        wmin[0] = 0.0
    dmin_max = mnr[1].max()
    for i in range(len(mnr[0]) - 1, 0, -1):         # only the minor direction (as in the source)
        if mnr[2][i] < 1.0 and mnr[1][i] > 0.5 * dmin_max:
            if mnr[2][i] < mnr[2][i - 1]:
                wmin[i] = 0.0
    if maj[1][0] < mnr[1][0]:
        wmin[0] = 0.0
        wmaj[0] = 0.0
    res['exp_major'] = dict(num=maj[0], dist=maj[1], val=maj[2], wt=wmaj)
    res['exp_minor'] = dict(num=mnr[0], dist=mnr[1], val=mnr[2], wt=wmin)

    # --- automatic fit ---------------------------------------------------------------------
    ahmax_b = (0.75 * maj[1].min(), 4.0 * maj[1].max())
    ahmin_b = (0.75 * mnr[1].min(), 4.0 * mnr[1].max())
    lagdist = np.concatenate([maj[1], mnr[1]])
    value = np.concatenate([maj[2], mnr[2]])
    numpairs = np.concatenate([maj[0], mnr[0]]).astype(int)
    wts = np.concatenate([wmaj, wmin])
    is_major = np.concatenate([np.ones(len(maj[1]), bool), np.zeros(len(mnr[1]), bool)])
    fit = varmodelfit(lagdist, value, numpairs, wts, is_major, ahmax_b, ahmin_b,
                      maxiter=maxiter, rseed=rseed, progress=progress, nst=nst, legacy=legacy)
    res['fit'] = fit
    res['model'] = dict(nst=nst, c0=0.0, cc=fit['cc'], ahmax=fit['ahmax'], ahmin=fit['ahmin'],
                        azimuth=a1, azimuth_major=az_major, azimuth_minor=az_minor)
    res['range_bounds'] = dict(major=ahmax_b, minor=ahmin_b)
    return res


# ---------------------------------------------------------------------------
# Output in the same layout as the Fortran program (to compare file by file)
# ---------------------------------------------------------------------------
def format_experimental(res):
    out = ["Experimental variogram", "6", "Variogram Index", "Lag Distance",
           "Number of Pairs", "Variogram Value", "Inferred Azimuth", "Inferred Dip"]
    for idx, key, az in ((1, 'exp_major', res['az_major']), (2, 'exp_minor', res['az_minor'])):
        e = res[key]
        for d, n, vl in zip(e['dist'], e['num'], e['val']):
            out.append("%1d%16.4f%10d%14.4f%14.4f%14.4f" % (idx, d, int(n), vl, az, 0.0))
    return "\n".join(out) + "\n"


def format_gslib_model(res):
    m = res['model']
    a1 = res['angle_gslib']
    out = ["%3d%10.2f    -nst, nugget effect" % (m['nst'], 0.0)]
    for c, amx, amn in zip(m['cc'], m['ahmax'], m['ahmin']):
        out.append("%3d%10.2f%10.2f%10.2f%10.2f    -it,cc,azm,dip,tilt" % (1, c, a1, 0.0, 0.0))
        out.append("%14.2f%14.2f%14.2f    -a_hmax, a_hmin, a_vert (ranges)" % (amx, amn, 0.0))
    return "\n".join(out) + "\n"


def format_model_points(res):
    m = res['model']
    out = ["Variogram model points", "7", "Variogram Index", "Lag Distance",
           "Number of Pairs", "Variogram Value", "Variogram Number",
           "Inferred Azimuth", "Inferred Dip"]
    for idx, rng, az in ((1, m['ahmax'], res['az_major']), (2, m['ahmin'], res['az_minor'])):
        sep = 1.1 * max(m['ahmax'] if idx == 1 else m['ahmin']) / 1000.0
        hh = sep * np.arange(1, 1001)
        g = model_gamma(hh, m['cc'], rng)
        for d, gv in zip(hh, g):
            out.append("%6d%18.6f%10d%18.6f%6d%18.6f%18.6f" % (idx, d, 1, gv, 1, az, 0.0))
    return "\n".join(out) + "\n"


AUTOVAR_DIFFERENCES = [
    "ndt (number of data scanned by datspac_horz) is uninitialised in the Fortran source when no "
    "drill-hole id is used; the port uses ndt = number of data.",
    "The Fortran program reads the data as REAL*4 and runs the Jacobi diagonalisation in single "
    "precision; the port uses double precision except inside jacobi (emulated REAL*4).",
    "Projection onto the principal axes in the experimental variogram: the Fortran code multiplies "
    "xyz by the transposed eigenvector matrix (valid only if the matrix is symmetric); the port "
    "projects onto the eigenvectors (identical when the matrix is symmetric).",
    "Fortran bugs reproduced on purpose: moi uses the first nd elements instead of the selected "
    "ones; varnumpairs is indexed with the filtered index in varmodelfit.",
    "Fit objective evaluated in double precision (Fortran: GSLIB cova33 in REAL*4); with 2,000,000 "
    "accept-only-improvement steps the random path is identical, the final values can differ by "
    "the single-precision noise.",
]
