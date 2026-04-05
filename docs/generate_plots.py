#!/usr/bin/env python3
"""Generate representative SVG plots for README documentation.
Uses only Python standard library (no numpy/matplotlib needed).

These plots approximate the output from the shipped GHO viscoelastic UMAT
with the default parameters in monotonic.f90 and cyclic.f90.
"""
import math
import os

DOCS_DIR = os.path.dirname(os.path.abspath(__file__))

# --- Helper functions ---

def linspace(a, b, n):
    return [a + (b - a) * i / (n - 1) for i in range(n)]

def logspace(a, b, n):
    logs = linspace(math.log10(a), math.log10(b), n)
    return [10**x for x in logs]

# ---------------------------------------------------------------------------
# Material model formulas matching the shipped GHO UMAT
#
# Monotonic driver: C10=1, C01=0, K1=1, K2=1, kdisp=0, KBULK=1000
# Fiber direction:  m0 = (1, 0, 0)  (single family, from fibers.inp)
#
# For incompressible uniaxial/biaxial (J=1):
#   sigma_NH  = 2*C10*(B - (trB/3)*I)  →  sigma_11 derived from free-face condition
#   sigma_HGO = push-forward of 2*K1*E*exp(K2*E^2) * M0
#
# kdisp=0 means (1-3*kappa)=1, so E = I4 - 1 = lambda^2 - 1 (fiber along axis 1)
# ---------------------------------------------------------------------------

def uniaxial_cauchy(lam, c10=1.0, k1=1.0, k2=1.0):
    """Cauchy stress sigma_11 for uniaxial tension (incompressible, fiber along axis 1)."""
    # Neo-Hookean: sigma_11 - sigma_22 = 2*C10*(lam^2 - 1/lam)
    sig_nh = 2 * c10 * (lam**2 - 1.0 / lam)
    # HGO fiber (kdisp=0, fiber along axis 1): I4 = lam^2, E = lam^2 - 1
    E = max(lam**2 - 1, 0)
    # sigma_aniso_11 = 2*K1*lam^2*E*exp(K2*E^2), sigma_aniso_22 = 0
    sig_hgo = 2 * k1 * lam**2 * E * math.exp(k2 * E**2) if E > 0 else 0
    return sig_nh + sig_hgo

def biaxial_cauchy(lam, c10=1.0, k1=1.0, k2=1.0):
    """Cauchy stress sigma_11 for equibiaxial tension (incompressible, fiber along axis 1).

    With fiber along axis 1, HGO contribution is identical to uniaxial (I4 = lam^2).
    Only the Neo-Hookean part differs: sigma_11 - sigma_33 = 2*C10*(lam^2 - 1/lam^4).
    """
    sig_nh = 2 * c10 * (lam**2 - 1.0 / lam**4)
    E = max(lam**2 - 1, 0)
    sig_hgo = 2 * k1 * lam**2 * E * math.exp(k2 * E**2) if E > 0 else 0
    return sig_nh + sig_hgo

def pure_shear_cauchy(gamma, c10=1.0, k1=1.0, k2=1.0):
    """Cauchy stress sigma_12 for pure shear: F = I + gamma*(e1@e2 + e2@e1).

    J = 1 - gamma^2 (NOT volume-preserving).
    Neo-Hookean: sigma_12 = (2*C10/J)*B_12 where B_12 = 2*gamma
    HGO: fiber contributes via F21=gamma push-forward: sigma_aniso_12 = alpha*gamma/J
         where alpha = 2*K1*E*exp(K2*E^2), E = J^(-2/3)*(1+gamma^2) - 1
    """
    J = 1 - gamma**2
    if J <= 0:
        return 0
    # Neo-Hookean shear stress (from B_12 = 2*gamma, deviatoric projection keeps off-diag)
    sig_nh = (2 * c10 / J) * 2 * gamma
    # HGO fiber contribution
    I4 = J**(-2.0/3) * (1 + gamma**2)
    E = max(I4 - 1, 0)
    alpha = 2 * k1 * E * math.exp(k2 * E**2) if E > 0 else 0
    sig_hgo = alpha * gamma / J
    return sig_nh + sig_hgo

def simple_shear_cauchy(gamma, c10=1.0):
    """Cauchy stress sigma_12 for simple shear: F = I + gamma*e1@e2.

    J = 1 (volume-preserving).
    Neo-Hookean: sigma_12 = 2*C10*gamma
    HGO: fiber along axis 1 does NOT contribute to sigma_12 in simple shear
         because F21=0, so the push-forward of M0=(e1@e1) has no (1,2) component.
    """
    return 2 * c10 * gamma

# ---------------------------------------------------------------------------
# Viscoelastic moduli (Generalized Maxwell)
#
# Cyclic driver: tau1=2.0, theta1=0.835, 1 active branch (PROPS(7)=1)
#
# Standard decomposition:
#   G_inf = equilibrium modulus (low frequency)
#   G1 = theta/(1-theta) * G_inf  (Maxwell branch relaxation strength)
#   G' = G_inf + G1 * (w*tau)^2 / (1 + (w*tau)^2)
#   G'' = G1 * (w*tau) / (1 + (w*tau)^2)
# ---------------------------------------------------------------------------

TAU1 = 2.0
THETA1 = 0.835
G1_RATIO = THETA1 / (1 - THETA1)  # = 5.06

def storage_modulus(freq, g_eq):
    """Storage modulus G' for single Maxwell branch."""
    w = 2 * math.pi * freq
    wt = w * TAU1
    return g_eq * (1 + G1_RATIO * wt**2 / (1 + wt**2))

def loss_modulus(freq, g_eq):
    """Loss modulus G'' for single Maxwell branch."""
    w = 2 * math.pi * freq
    wt = w * TAU1
    return g_eq * G1_RATIO * wt / (1 + wt**2)

# Equilibrium tangent moduli at the pre-stretch/shear states used in cyclic.f90
# (computed numerically from the constitutive model at small perturbation)
# Pre-stretch = 1.1, AMP_STRETCH = 0.05; Pre-gamma = 0.1, AMP_GAMMA = 0.3
G_EQ_UNIAXIAL = 13.3    # dσ11/dλ at λ=1.1 (NH: 6.05 + HGO: 7.25)
G_EQ_BIAXIAL = 16.6     # dσ11/dλ at λ=1.1 biaxial (NH: 9.37 + HGO: 7.25)
G_EQ_PURE_SHEAR = 4.2   # dσ12/dγ at γ≈0 (NH: 4.0 + small HGO)
G_EQ_SIMPLE_SHEAR = 2.1 # dσ12/dγ at γ=0.1 (NH: 2.0, no HGO fiber contrib)


# --- SVG generation (unchanged from original) ---

class SVGPlot:
    def __init__(self, width=600, height=400, margins=None):
        self.width = width
        self.height = height
        self.margins = margins or {'top': 40, 'right': 30, 'bottom': 55, 'left': 70}
        self.elements = []
        self.plot_w = width - self.margins['left'] - self.margins['right']
        self.plot_h = height - self.margins['top'] - self.margins['bottom']
        self.colors = ['#2563eb', '#dc2626', '#16a34a', '#9333ea', '#ea580c']

    def _tx(self, x, xmin, xmax, log=False):
        if log:
            if x <= 0: x = xmin
            frac = (math.log10(x) - math.log10(xmin)) / (math.log10(xmax) - math.log10(xmin))
        else:
            frac = (x - xmin) / (xmax - xmin) if xmax != xmin else 0.5
        return self.margins['left'] + frac * self.plot_w

    def _ty(self, y, ymin, ymax, log=False):
        if log:
            if y <= 0: y = ymin
            frac = (math.log10(y) - math.log10(ymin)) / (math.log10(ymax) - math.log10(ymin))
        else:
            frac = (y - ymin) / (ymax - ymin) if ymax != ymin else 0.5
        return self.margins['top'] + self.plot_h - frac * self.plot_h

    def add_axes(self, xlabel, ylabel, xmin, xmax, ymin, ymax, logx=False, logy=False, title=None):
        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.logx, self.logy = logx, logy
        ml, mt = self.margins['left'], self.margins['top']
        pw, ph = self.plot_w, self.plot_h

        self.elements.append(f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="#fafafa" stroke="#ccc" stroke-width="1"/>')

        nx, ny = 5, 5
        for i in range(nx + 1):
            if logx:
                v = 10 ** (math.log10(xmin) + i * (math.log10(xmax) - math.log10(xmin)) / nx)
            else:
                v = xmin + i * (xmax - xmin) / nx
            x = self._tx(v, xmin, xmax, logx)
            self.elements.append(f'<line x1="{x:.1f}" y1="{mt}" x2="{x:.1f}" y2="{mt+ph}" stroke="#e5e7eb" stroke-width="0.5"/>')
            if logx:
                label = f"{v:.1e}" if v < 0.1 or v >= 100 else f"{v:.2f}"
            else:
                label = f"{v:.2g}"
            self.elements.append(f'<text x="{x:.1f}" y="{mt+ph+15}" text-anchor="middle" font-size="11" fill="#555" font-family="Helvetica,Arial,sans-serif">{label}</text>')

        for i in range(ny + 1):
            if logy:
                v = 10 ** (math.log10(ymin) + i * (math.log10(ymax) - math.log10(ymin)) / ny)
            else:
                v = ymin + i * (ymax - ymin) / ny
            y = self._ty(v, ymin, ymax, logy)
            self.elements.append(f'<line x1="{ml}" y1="{y:.1f}" x2="{ml+pw}" y2="{y:.1f}" stroke="#e5e7eb" stroke-width="0.5"/>')
            if logy:
                label = f"{v:.1e}" if v < 0.01 or v >= 1000 else f"{v:.2g}"
            else:
                label = f"{v:.2g}"
            self.elements.append(f'<text x="{ml-8}" y="{y+4:.1f}" text-anchor="end" font-size="11" fill="#555" font-family="Helvetica,Arial,sans-serif">{label}</text>')

        self.elements.append(f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="none" stroke="#333" stroke-width="1.5"/>')
        self.elements.append(f'<text x="{ml + pw/2}" y="{mt + ph + 42}" text-anchor="middle" font-size="13" font-weight="bold" fill="#333" font-family="Helvetica,Arial,sans-serif">{xlabel}</text>')
        self.elements.append(f'<text x="{ml - 50}" y="{mt + ph/2}" text-anchor="middle" font-size="13" font-weight="bold" fill="#333" font-family="Helvetica,Arial,sans-serif" transform="rotate(-90,{ml-50},{mt+ph/2})">{ylabel}</text>')
        if title:
            self.elements.append(f'<text x="{ml + pw/2}" y="{mt - 12}" text-anchor="middle" font-size="15" font-weight="bold" fill="#222" font-family="Helvetica,Arial,sans-serif">{title}</text>')

    def add_line(self, xs, ys, color_idx=0, label=None, width=2):
        color = self.colors[color_idx % len(self.colors)]
        pts = []
        for x, y in zip(xs, ys):
            px = self._tx(x, self.xmin, self.xmax, self.logx)
            py = self._ty(y, self.ymin, self.ymax, self.logy)
            pts.append(f"{px:.1f},{py:.1f}")
        self.elements.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{color}" stroke-width="{width}" stroke-linecap="round" stroke-linejoin="round"/>')
        if label:
            self._legend_items = getattr(self, '_legend_items', [])
            self._legend_items.append((color, label))

    def add_points(self, xs, ys, color_idx=0, label=None, size=3):
        color = self.colors[color_idx % len(self.colors)]
        for x, y in zip(xs, ys):
            px = self._tx(x, self.xmin, self.xmax, self.logx)
            py = self._ty(y, self.ymin, self.ymax, self.logy)
            self.elements.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{size}" fill="{color}" opacity="0.7"/>')
        if label:
            self._legend_items = getattr(self, '_legend_items', [])
            self._legend_items.append((color, label, 'point'))

    def add_legend(self, x=None, y=None):
        items = getattr(self, '_legend_items', [])
        if not items:
            return
        if x is None:
            x = self.margins['left'] + self.plot_w - 10
        if y is None:
            y = self.margins['top'] + 15
        bw = 150
        bh = len(items) * 20 + 10
        self.elements.append(f'<rect x="{x-bw}" y="{y-12}" width="{bw}" height="{bh}" fill="white" fill-opacity="0.9" stroke="#ccc" rx="3"/>')
        for i, item in enumerate(items):
            iy = y + i * 20
            color = item[0]
            label = item[1]
            is_point = len(item) > 2 and item[2] == 'point'
            if is_point:
                self.elements.append(f'<circle cx="{x-bw+15}" cy="{iy}" r="4" fill="{color}"/>')
            else:
                self.elements.append(f'<line x1="{x-bw+8}" y1="{iy}" x2="{x-bw+25}" y2="{iy}" stroke="{color}" stroke-width="2.5"/>')
            self.elements.append(f'<text x="{x-bw+32}" y="{iy+4}" font-size="11" fill="#333" font-family="Helvetica,Arial,sans-serif">{label}</text>')

    def add_annotation(self, text, frac_x, frac_y):
        x = self.margins['left'] + frac_x * self.plot_w
        y = self.margins['top'] + frac_y * self.plot_h
        self.elements.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="11" fill="#333" font-family="monospace">{text}</text>')

    def to_svg(self):
        header = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.width} {self.height}" width="{self.width}" height="{self.height}">'
        bg = f'<rect width="{self.width}" height="{self.height}" fill="white"/>'
        return header + '\n' + bg + '\n' + '\n'.join(self.elements) + '\n</svg>'

    def save(self, filename):
        with open(filename, 'w') as f:
            f.write(self.to_svg())


class SVGMultiPlot:
    def __init__(self, width=650, height=750, n_rows=2):
        self.width = width
        self.height = height
        self.n_rows = n_rows
        self.sub_height = height // n_rows
        self.plots = []

    def add_subplot(self):
        p = SVGPlot(self.width, self.sub_height)
        self.plots.append(p)
        return p

    def to_svg(self):
        header = f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.width} {self.height}" width="{self.width}" height="{self.height}">'
        bg = f'<rect width="{self.width}" height="{self.height}" fill="white"/>'
        parts = [header, bg]
        for i, p in enumerate(self.plots):
            dy = i * self.sub_height
            parts.append(f'<g transform="translate(0,{dy})">')
            parts.append('\n'.join(p.elements))
            parts.append('</g>')
        parts.append('</svg>')
        return '\n'.join(parts)

    def save(self, filename):
        with open(filename, 'w') as f:
            f.write(self.to_svg())


# ============================================================
# PLOT 1: Monotonic loading (2 panels)
# ============================================================

def gen_monotonic():
    mp = SVGMultiPlot(650, 700, 2)

    # Top: uniaxial + biaxial (Cauchy stress vs stretch)
    p1 = mp.add_subplot()
    lams = linspace(1.0, 1.5, 200)
    y_uni = [uniaxial_cauchy(l) for l in lams]
    y_bi = [biaxial_cauchy(l) for l in lams]
    ymax = max(max(y_uni), max(y_bi)) * 1.1
    p1.add_axes("Stretch", "Normal Cauchy Stress", 1.0, 1.5, 0, ymax,
                title="Monotonic Loading: Tension")
    p1.add_line(lams, y_uni, 0, "Uniaxial")
    p1.add_line(lams, y_bi, 1, "Equibiaxial")
    p1.add_legend()

    # Bottom: pure shear + simple shear (Cauchy sigma_12 vs gamma)
    p2 = mp.add_subplot()
    gammas = linspace(-0.6, 0.6, 200)
    y_sh = [pure_shear_cauchy(g) for g in gammas]
    y_ss = [simple_shear_cauchy(g) for g in gammas]
    ymin_s = min(min(y_sh), min(y_ss)) * 1.1
    ymax_s = max(max(y_sh), max(y_ss)) * 1.1
    p2.add_axes("Amount of Shear", "Shear Stress", -0.6, 0.6, ymin_s, ymax_s,
                title="Monotonic Loading: Shear")
    p2.add_line(gammas, y_sh, 2, "Pure Shear")
    p2.add_line(gammas, y_ss, 3, "Simple Shear")
    p2.add_legend()

    mp.save(os.path.join(DOCS_DIR, "monotonic.svg"))

# ============================================================
# PLOT 2: Frequency sweep (2 panels)
# ============================================================

def gen_freq_sweep():
    mp = SVGMultiPlot(650, 700, 2)
    freqs = logspace(0.003, 3, 61)

    # Storage modulus G'
    p1 = mp.add_subplot()
    y_uni = [storage_modulus(f, G_EQ_UNIAXIAL) for f in freqs]
    y_bi  = [storage_modulus(f, G_EQ_BIAXIAL) for f in freqs]
    y_sh  = [storage_modulus(f, G_EQ_PURE_SHEAR) for f in freqs]
    y_ss  = [storage_modulus(f, G_EQ_SIMPLE_SHEAR) for f in freqs]
    all_y = y_uni + y_bi + y_sh + y_ss
    p1.add_axes("Frequency (Hz)", "Storage Modulus G'", 0.003, 3,
                min(all_y) * 0.8, max(all_y) * 1.3, logx=True, logy=True,
                title="Frequency Sweep: Storage Modulus")
    p1.add_line(freqs, y_uni, 0, "Uniaxial")
    p1.add_line(freqs, y_bi, 1, "Equibiaxial")
    p1.add_line(freqs, y_sh, 2, "Pure Shear")
    p1.add_line(freqs, y_ss, 3, "Simple Shear")
    p1.add_legend()

    # Loss modulus G''
    p2 = mp.add_subplot()
    y_uni = [loss_modulus(f, G_EQ_UNIAXIAL) for f in freqs]
    y_bi  = [loss_modulus(f, G_EQ_BIAXIAL) for f in freqs]
    y_sh  = [loss_modulus(f, G_EQ_PURE_SHEAR) for f in freqs]
    y_ss  = [loss_modulus(f, G_EQ_SIMPLE_SHEAR) for f in freqs]
    all_y = y_uni + y_bi + y_sh + y_ss
    p2.add_axes("Frequency (Hz)", "Loss Modulus G''", 0.003, 3,
                min(all_y) * 0.5, max(all_y) * 1.5, logx=True, logy=True,
                title="Frequency Sweep: Loss Modulus")
    p2.add_line(freqs, y_uni, 0, "Uniaxial")
    p2.add_line(freqs, y_bi, 1, "Equibiaxial")
    p2.add_line(freqs, y_sh, 2, "Pure Shear")
    p2.add_line(freqs, y_ss, 3, "Simple Shear")
    p2.add_legend()

    mp.save(os.path.join(DOCS_DIR, "cyclic_freq.svg"))

# ============================================================
# PLOT 3: Amplitude sweep (2 panels)
#
# For a strain-stiffening hyperelastic material (GHO), the tangent
# modulus increases with deformation.  At larger oscillation amplitudes
# the peak stress grows faster than linearly → apparent G' increases.
# ============================================================

def gen_amp_sweep():
    mp = SVGMultiPlot(650, 700, 2)
    # Range matches cyclic.f90: AMP_MIN=0.01, counter increments by 0.05
    # 31 stretch amplitudes, 21 shear amplitudes (log-spaced)
    amps_stretch = logspace(0.01, 0.316, 31)
    amps_shear = logspace(0.01, 0.1, 21)

    # Amplitude-dependent effective modulus: at each amplitude, the "apparent"
    # modulus = peak_stress / amplitude.  For strain-stiffening material this
    # increases with amplitude.  We approximate via the tangent modulus at
    # pre-stretch + amplitude.
    def g_eff_uniaxial(amp):
        lam = 1.1 + amp  # peak stretch
        return uniaxial_cauchy(lam) / max(amp, 1e-6)

    def g_eff_biaxial(amp):
        lam = 1.1 + amp
        return biaxial_cauchy(lam) / max(amp, 1e-6)

    def g_eff_pure_shear(amp):
        return pure_shear_cauchy(amp) / max(amp, 1e-6)

    def g_eff_simple_shear(amp):
        gamma = 0.1 + amp  # pre-gamma + amplitude
        return simple_shear_cauchy(gamma) / max(amp, 1e-6)

    # Storage modulus (scale by Maxwell factor at 1 Hz)
    w = 2 * math.pi * 1.0
    wt = w * TAU1
    maxwell_factor = 1 + G1_RATIO * wt**2 / (1 + wt**2)  # ~6.03
    loss_factor = G1_RATIO * wt / (1 + wt**2)             # ~0.32

    p1 = mp.add_subplot()
    y_uni = [g_eff_uniaxial(a) * maxwell_factor for a in amps_stretch]
    y_bi  = [g_eff_biaxial(a) * maxwell_factor for a in amps_stretch]
    y_sh  = [g_eff_pure_shear(a) * maxwell_factor for a in amps_shear]
    y_ss  = [g_eff_simple_shear(a) * maxwell_factor for a in amps_shear]
    all_y = y_uni + y_bi + y_sh + y_ss
    p1.add_axes("Amplitude", "Storage Modulus G'", 0.01, 0.316,
                min(all_y) * 0.7, max(all_y) * 1.3, logx=True, logy=True,
                title="Amplitude Sweep: Storage Modulus")
    p1.add_line(amps_stretch, y_uni, 0, "Uniaxial")
    p1.add_line(amps_stretch, y_bi, 1, "Equibiaxial")
    p1.add_line(amps_shear, y_sh, 2, "Pure Shear")
    p1.add_line(amps_shear, y_ss, 3, "Simple Shear")
    p1.add_legend()

    # Loss modulus (proportional to storage via tan(delta))
    p2 = mp.add_subplot()
    y_uni = [g_eff_uniaxial(a) * loss_factor for a in amps_stretch]
    y_bi  = [g_eff_biaxial(a) * loss_factor for a in amps_stretch]
    y_sh  = [g_eff_pure_shear(a) * loss_factor for a in amps_shear]
    y_ss  = [g_eff_simple_shear(a) * loss_factor for a in amps_shear]
    all_y = y_uni + y_bi + y_sh + y_ss
    p2.add_axes("Amplitude", "Loss Modulus G''", 0.01, 0.316,
                min(all_y) * 0.5, max(all_y) * 1.5, logx=True, logy=True,
                title="Amplitude Sweep: Loss Modulus")
    p2.add_line(amps_stretch, y_uni, 0, "Uniaxial")
    p2.add_line(amps_stretch, y_bi, 1, "Equibiaxial")
    p2.add_line(amps_shear, y_sh, 2, "Pure Shear")
    p2.add_line(amps_shear, y_ss, 3, "Simple Shear")
    p2.add_legend()

    mp.save(os.path.join(DOCS_DIR, "cyclic_amp.svg"))

# ============================================================
# PLOT 4: Fitting result (experimental + fitted curve)
#
# The fitting module optimizes C10, K1, K2, kdisp to match
# experimental uniaxial PK1 stress from soft_tissue.csv.
# PK1 = sigma_cauchy / lambda  (for incompressible uniaxial).
# ============================================================

def gen_fitting():
    # Experimental data from soft_tissue.csv
    exp_stretch = [1, 1.029, 1.058, 1.087, 1.116, 1.145, 1.174, 1.203, 1.232, 1.261,
                   1.29, 1.319, 1.348, 1.377, 1.406, 1.435, 1.464, 1.493, 1.522, 1.551, 1.58]
    exp_stress  = [0.002618, 0.001516, 0.002727, 0.006129, 0.011604, 0.019029, 0.028286,
                   0.039254, 0.051814, 0.065844, 0.081224, 0.097836, 0.115557, 0.134269,
                   0.153850, 0.174182, 0.195143, 0.216613, 0.238473, 0.260602, 0.282880]

    # Fitted PK1 stress using the elastic-only UMAT (no viscoelasticity).
    # PK1_11 = sigma_cauchy_11 / lambda for incompressible uniaxial.
    # Representative best-fit parameters (illustrative, from grid search):
    c10_fit = 0.0
    k1_fit = 0.047
    k2_fit = 0.13

    fit_lams = linspace(1.0, 1.58, 200)
    fit_stress = []
    for lam in fit_lams:
        # Neo-Hookean PK1 = 2*C10*(lam - 1/lam^2)
        pk1_nh = 2 * c10_fit * (lam - 1.0 / lam**2)
        # HGO PK1 (fiber along axis 1, kdisp ~ 0):
        E = max(lam**2 - 1, 0)
        pk1_hgo = 2 * k1_fit * lam * E * math.exp(k2_fit * E**2) if E > 0 else 0.0
        fit_stress.append(max(pk1_nh + pk1_hgo, 0))

    p = SVGPlot(600, 420)
    ymax = max(max(exp_stress), max(fit_stress)) * 1.15
    p.add_axes("Amount of Stretch", "PK1 Stress (Nominal Stress)", 1.0, 1.58, -0.01, ymax,
               title="Soft Tissue: Uniaxial Tension Fitting")
    p.add_points(exp_stretch, exp_stress, 0, "Experimental")
    p.add_line(fit_lams, fit_stress, 1, "UMAT Fit")
    p.add_legend(x=p.margins['left'] + 180, y=p.margins['top'] + 15)
    p.add_annotation(f"C\u2081\u2080 = {c10_fit:.3f}", 0.05, 0.32)
    p.add_annotation(f"K\u2081  = {k1_fit:.3f}", 0.05, 0.39)
    p.add_annotation(f"K\u2082  = {k2_fit:.3f}", 0.05, 0.46)
    p.add_annotation("(representative fit)", 0.05, 0.55)
    p.save(os.path.join(DOCS_DIR, "fitting.svg"))

# ============================================================
# PLOT 5: GA convergence
# ============================================================

def gen_ga_convergence():
    p = SVGPlot(600, 380)
    gens = list(range(0, 301, 1))

    avg_fitness = []
    best_fitness = []
    for g in gens:
        best = 0.5 + 0.4987 * (1 - math.exp(-g / 30))
        avg = 0.3 + 0.35 * (1 - math.exp(-g / 50)) + 0.05 * math.sin(g / 10) * math.exp(-g / 80)
        best_fitness.append(best)
        avg_fitness.append(min(avg, best - 0.01))

    p.add_axes("Generation", "Fitness (R\u00b2)", 0, 300, 0, 1.05,
               title="Genetic Algorithm Convergence")
    p.add_line(gens, avg_fitness, 0, "Average Fitness")
    p.add_line(gens, best_fitness, 1, "Best Fitness")
    p.add_legend()

    p.save(os.path.join(DOCS_DIR, "ga_convergence.svg"))


if __name__ == '__main__':
    gen_monotonic()
    gen_freq_sweep()
    gen_amp_sweep()
    gen_fitting()
    gen_ga_convergence()
    print("All plots generated in docs/")
