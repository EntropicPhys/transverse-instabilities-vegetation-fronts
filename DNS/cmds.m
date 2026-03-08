import numpy as np
import dedalus.public as d3
import logging
from scipy.optimize import fsolve
from mpi4py import MPI

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

comm  = MPI.COMM_WORLD
rank  = comm.rank
nproc = comm.size

def eval_global_scalar(expr):
    """Evaluate a Dedalus scalar expression and return a Python float on all MPI ranks."""
    fld = expr.evaluate()
    arr = fld['g']
    local = float(arr.flat[0]) if arr.size else 0.0
    return comm.allreduce(local, op=MPI.SUM)
    
def tail_fraction(obj,scales):
    field = obj.evaluate() if hasattr(obj, "evaluate") else obj
    field.change_scales(1)
    field.require_coeff_space()
    tmp = field.copy()
    tmp.low_pass_filter(scales=scales)
    diff = field-tmp
    return eval_global_scalar(integ(diff*diff)/(integ(field*field)+1e-30))


# ----------------------------
# Simulation Params
# ----------------------------
Lx, Ly = 280, 600
Nx, Ny = 2*280,2*600

# dimensional params
PP=144.5; Lam=0.33/2; Nw=6.4; Nh=7;  Rw=0.73; Rh=0.3; Ga=1.9
A=20; Q=0.065;  K=0.8; M=3.55; sx=0; sy=0; 
Db=1; Dw=350; Dh=3000; EE = 1.58; f  = 0.4
# non-dimensional params
nw=Nw/M; nh=Nh/M; a=A/M; q=Q/K; et=EE*K; pp=PP*Lam/M**2; ga=Ga*K/M; dw=Dw/Db; dh=Dh/Db

# homogeneous states
def hss(x):
    b, w, h = x
    I  = a*(b+q*f)/(b+q)
    Lw = nw/(1+Rw*b)
    Lh = nh/(1+Rh*b)
    return (
        w*b*(1-b)*(1+et*b)**2 - b,
        I*h - Lw*w - ga*w*b*(1+et*b)**2,
        pp - I*h - Lh*h
    )

b00=0.65; w00=0.645
h00=pp/(a*(b00+q*f)/(b00+q)+nw/(1+Rw*b00))
b0, w0, h0 = fsolve(hss, [b00, w00, h00])

dealias = 1.1
stop_sim_time = 1e1
timestepper = d3.RK443
timestep = 1.3e-1
dtype = np.float64

# ----------------------------
# Bases + Distributor
# ----------------------------
coords = d3.CartesianCoordinates('x','y')
mesh = (nproc,) if nproc > 1 else None
dist = d3.Distributor(coords, dtype=dtype, mesh=mesh) if mesh else d3.Distributor(coords, dtype=dtype)
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(0, Lx), dealias=dealias)
ybasis = d3.ChebyshevT(coords['y'], size=Ny, bounds=(0, Ly), dealias=dealias)

# Helper: integrate over the whole domain
integ = lambda A: d3.Integrate(A, ('x','y'))
Area  = Lx * Ly
# ----------------------------
# Fields
# ----------------------------
b = dist.Field(name='b', bases=(xbasis, ybasis))
w = dist.Field(name='w', bases=(xbasis, ybasis))
h = dist.Field(name='h', bases=(xbasis, ybasis))
u = dist.VectorField(coords, name='u', bases=(xbasis, ybasis))
p = dist.Field(name='p')   # scalar
v = dist.Field(name='v')   # scalar

tau_bx  = dist.Field(name='tau_bx', bases=xbasis)
tau_bxx = dist.Field(name='tau_bxx', bases=xbasis)
tau_wx  = dist.Field(name='tau_wx', bases=xbasis)
tau_wxx = dist.Field(name='tau_wxx', bases=xbasis)
tau_hx  = dist.Field(name='tau_hx', bases=xbasis)
tau_hxx = dist.Field(name='tau_hxx', bases=xbasis)

# ----------------------------
# Substitutions / Operators
# ----------------------------
x, y   = dist.local_grids(xbasis, ybasis)
ex, ey = coords.unit_vector_fields(dist)
lift_basis = ybasis.derivative_basis(1)
lift = lambda A: d3.Lift(A, lift_basis, -1)
grad_b = d3.grad(b) + ey*lift(tau_bx)
grad_w = d3.grad(w) + ey*lift(tau_wx)
grad_h = d3.grad(h) + ey*lift(tau_hx)
by = ey @ grad_b
wy = ey @ grad_w
hy = ey @ grad_h
bLap = d3.div(grad_b) + lift(tau_bxx)
wLap = d3.div(grad_w) + lift(tau_wxx)
hLap = d3.div(grad_h) + lift(tau_hxx)
nodalb = w*b*(1-b)*(1+et*b)**2 - b
nodalw = h*a*(b+q*f)/(b+q) - nw*w/(1+Rw*b) - ga*w*b*(1+et*b)**2
nodalh = p - h*a*(b+q*f)/(b+q) - nh*h/(1+Rh*b)

in1 = (
    by*(nodalb + bLap) +
    wy*(nodalw + wLap) +
    hy*(nodalh + hLap)
)
in2 = by**2 + wy**2 + hy**2

# # x-derivative
bx = ex @ grad_b
bx_rms = integ(bx*bx)
b_total = integ(b)/(Lx*Ly)

L_b = -b+bLap
N_b = w*b*(1-b)*(1+et*b)**2+v*by
L_w = wLap
N_w = h*a*(b+q*f)/(b+q) - nw*w/(1+Rw*b) - ga*w*b*(1+et*b)**2 + v*wy
L_h = hLap
N_h = -h*a*(b+q*f)/(b+q) - nh*h/(1+Rh*b) + v*hy

# ----------------------------
# Stationarity proxy: RHS norms (time-derivative residuals)
# ----------------------------
dbdt = -b + bLap + w*b*(1-b)*(1+et*b)**2 + v*by
dwdt =  dw*wLap + h*a*(b+q*f)/(b+q) - nw*w/(1+Rw*b) - ga*w*b*(1+et*b)**2 + v*wy
dhdt =  p + dh*hLap - h*a*(b+q*f)/(b+q) - nh*h/(1+Rh*b) + v*hy

R2 = integ(dbdt*dbdt + dwdt*dwdt + dhdt*dhdt) / Area
R  = R2**0.5

# Initialize v consistently (compute scalar from expression)
v0 = eval_global_scalar(v_expr)
if v['g'].size:
    v['g'][...] = v0

# ----------------------------
# Problem
# ----------------------------
problem = d3.IVP(
    [b, w, h, p, v, tau_bx, tau_bxx, tau_wx, tau_wxx, tau_hx, tau_hxx],
    namespace=locals()
)

problem.add_equation("dt(b) + b - bLap = w*b*(1-b)*(1+et*b)**2 + v*by")
problem.add_equation("dt(w) - dw*wLap = h*a*(b+q*f)/(b+q) - nw*w/(1+Rw*b) - ga*w*b*(1+et*b)**2 + v*wy")
problem.add_equation("dt(h) - p - dh*hLap = -h*a*(b+q*f)/(b+q) - nh*h/(1+Rh*b) + v*hy")
problem.add_equation("dt(p) = 0")
problem.add_equation("v = -integ(in1)/integ(in2)")

# Neumann BCs via tau/lift formulation
problem.add_equation("ey@grad_b(y=0)  = 0")
problem.add_equation("ey@grad_b(y=Ly) = 0")
problem.add_equation("ey@grad_w(y=0)  = 0")
problem.add_equation("ey@grad_w(y=Ly) = 0")
problem.add_equation("ey@grad_h(y=0)  = 0")
problem.add_equation("ey@grad_h(y=Ly) = 0")

# ----------------------------
# Solver
# ----------------------------
solver = problem.build_solver(timestepper)

# ----------------------------
# Stationarity criteria (RESIDUAL-based)  <<<< CHANGED
# ----------------------------
X_iter = 1000                 # check cadence in iterations

# ----------------------------
# Analysis: snapshots
# ----------------------------
outdir_var = f"bwh_FI_f_{float(f):.3f}_eta_{float(et):.3f}_var_2"
outdir_var = outdir_var.replace('.','p')
snapshots = solver.evaluator.add_file_handler(outdir_var, iter=X_iter, 
                                              max_writes=1,mode='overwrite')
snapshots.add_task(b, name='b')
snapshots.add_task(w, name='w')
snapshots.add_task(h, name='h')
snapshots.add_task(v, name='v')
snapshots.add_task(p, name='p')


# ----------------------------
# Main loop  
# ----------------------------
try:
    if rank == 0:
        logger.info("Starting main loop")
        logger.info(f"Using nproc={nproc}, mesh={mesh}")
        logger.info(f"RUN: f={f:.6g} | et=EE*K={et:.6g}" )

    while solver.proceed:

        v_now = eval_global_scalar(v)
        u['g'] = v_now
        dt = CFL.compute_timestep()
        solver.step(dt)

        if (solver.iteration - 1) % X_iter == 0:
            # Compute residuals
            R_now    = eval_global_scalar(R)
            bx_rms_now = eval_global_scalar(bx_rms)
            L_b_tail = tail_fraction(L_b_f,scales=(1.0,2/3))
            N_b_tail = tail_fraction(N_b_f,scales=(1.0,2/3))
            L_w_tail = tail_fraction(L_w_f,scales=(1.0,2/3))
            N_w_tail = tail_fraction(N_w_f,scales=(1.0,2/3))
            L_h_tail = tail_fraction(L_h_f,scales=(1.0,2/3))
            N_h_tail = tail_fraction(N_h_f,scales=(1.0,2/3))

            if rank == 0:
                logger.info(f"Iter={solver.iteration}, t={solver.sim_time:.3e}, dt={dt:.3e},"
                            f"R={R_now:.3e}, CFL={Ly/abs(Ny*v_now):.3e}, v={v_now:.3e}, bx_rms={bx_rms_now:.3e}") # 
                logger.info(f"tails:"
                            f"L_b[y]={L_b_tail:.2e}, N_b[y]={N_b_tail:.2e} |"
                            f"L_w[y]={L_w_tail:.2e}, N_w[y]={N_w_tail:.2e} |"
                            f"L_h[y]={L_h_tail:.2e}, N_h[y]={N_h_tail:.2e}" 
                            )


finally:
    solver.log_stats()



