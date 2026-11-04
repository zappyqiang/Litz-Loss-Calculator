"""2D FEM check of the analytical Litz model (NGSolve).
Run:  pip install ngsolve scipy numpy ; python fem_check.py 2   (argument = number of hex rings: 2 -> 19 strands, 3 -> 37)
Hex-packed bundle, every strand forced to carry I/N (ideal transposition).
Compares FEM AC resistance with the analytical skin + internal proximity model."""
import sys, numpy as np
from ngsolve import *
from netgen.occ import *

MU0 = 4e-7*np.pi; RHO = 1.724e-8; SIG = 1/RHO

def hex_centers(rings, pitch):
    pts = [(0.0, 0.0)]
    for n in range(1, rings+1):
        for k in range(6*n):
            seg, i = divmod(k, n)
            a0 = np.pi/3*seg; a1 = np.pi/3*(seg+1)
            p0 = np.array([np.cos(a0), np.sin(a0)])*n; p1 = np.array([np.cos(a1), np.sin(a1)])*n
            pts.append(tuple((p0 + (p1-p0)*i/n)*pitch))
    return pts

def build(centers, d):
    N = len(centers); r = d/2
    Rair = 25*max(max(np.hypot(*c) for c in centers), r) + 10*d
    outer = WorkPlane().Circle(Rair).Face(); outer.faces.name = "air"; outer.edges.name = "outer"
    shapes = [outer]
    for i, (x, y) in enumerate(centers):
        c = WorkPlane(Axes((x, y, 0), Z, X)).Circle(r).Face(); c.faces.name = f"s{i}"; c.faces.maxh = d/7
        shapes.append(c)
    # glue: strands sit inside the air disc
    g = Glue(shapes)
    mesh = Mesh(OCCGeometry(g, dim=2).GenerateMesh(maxh=Rair/5, grading=0.3))
    mesh.Curve(3)
    return mesh, N, r

def solve(mesh, N, r, f, Itot_rms=1.0):
    w = 2*np.pi*f
    Is = Itot_rms*np.sqrt(2)/N   # peak per strand
    H1c = H1(mesh, order=2, complex=True, dirichlet="outer")
    fes = FESpace([H1c] + [NumberSpace(mesh, complex=True)]*N)
    X, Y = fes.TnT()
    A, v = X[0], Y[0]
    a = BilinearForm(fes, symmetric=False)
    a += (1/MU0)*grad(A)*grad(v)*dx
    L = LinearForm(fes)
    for i in range(N):
        U, lam = X[i+1], Y[i+1]
        reg = mesh.Materials(f"s{i}")
        a += (1j*w*SIG*A*v - SIG*U*v)*dx(definedon=reg)
        a += (-1j*w*SIG*A*lam + SIG*U*lam)*dx(definedon=reg)   # constraint test fn
        L += (Is/(np.pi*r*r))*lam*dx(definedon=reg)
    a.Assemble(); L.Assemble()
    gf = GridFunction(fes)
    gf.vec.data = a.mat.Inverse(fes.FreeDofs()) * L.vec
    Ps = 0
    for i in range(N):
        reg = mesh.Materials(f"s{i}")
        J = -1j*w*SIG*gf.components[0] + SIG*gf.components[i+1]
        Ps += Integrate(((J*Conj(J)).real)/(2*SIG), mesh, definedon=reg)
    return Ps/(Itot_rms**2)   # R_ac per metre (ohm/m)


from scipy.special import iv

def analytic_rac_ratio(d, N, Db, f, T=20.0):
    """Same equations as the web calculator: skin + internal proximity, Rac/Rdc."""
    rho = RHO*(1+0.00393*(T-20)); sig = 1/rho; r = d/2; w = 2*np.pi*f
    delta = np.sqrt(2*rho/(w*MU0)); g = (1+1j)/delta; z = g*r
    FR = (z*iv(0, z)/(2*iv(1, z))).real
    C = 2*MU0/(g*iv(0, z))
    dI1 = iv(0, z) - iv(1, z)/z
    Ppk = np.pi*r*w/(2*MU0)*(1j*abs(C)**2*iv(1, z)*np.conj(g)*np.conj(dI1)).real
    kp = 2*Ppk                                        # W/m per (A/m rms)^2
    Hi2 = 1/(8*np.pi**2*(Db/2)**2)                    # for I = 1 A rms
    Rdc = rho/(N*np.pi*r**2)
    return FR + N*kp*Hi2/Rdc, Rdc

if __name__ == "__main__":
    d = 0.1e-3; rings = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    pitch = 1.15*d
    centers = hex_centers(rings, pitch); N = len(centers)
    mesh, N, r = build(centers, d)
    Acu = N*np.pi*d*d/4
    Db_eq = 2*np.sqrt(N*np.sqrt(3)/2*pitch**2/np.pi)    # equal-area bundle diameter
    print(f"N={N}, pitch={pitch*1e3:.3f} mm, Db_eq={Db_eq*1e3:.3f} mm, fill={Acu/(np.pi*Db_eq**2/4):.3f}, mesh elements={mesh.ne}")
    for f in [1e4, 1e5, 3e5, 1e6, 2e6, 4e6]:
        Rfem = solve(mesh, N, r, f)
        ratio, Rdc = analytic_rac_ratio(d, N, Db_eq, f)
        print(f"f={f:9.0f} Hz  Rac/Rdc  FEM={Rfem/Rdc:.5f}  analytic={ratio:.5f}  diff={100*(ratio/(Rfem/Rdc)-1):+.2f}%")
