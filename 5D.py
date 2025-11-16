from vpython import *
#Web VPython 3.2
# GlowScript 3.2 VPython
# 5D collision visualiser: H2O (molecule1) + H2 (molecule2)
from vpython import *

# ----------------- Parameters / geometry -----------------
scale = 0.7

# H2O geometry (molecule1)
OH = 0.5 * scale
HOH_angle = radians(104.5)
half_angle = HOH_angle / 2

# H2 geometry (molecule2)
HH = 0.74 * scale   # H-H bond length

# Degrees of freedom (initial)
th1 = 0.0    # θ for molecule1 (H2O)
phi1 = 0.0   # φ for molecule1 (tilt)
th2 = 0.0    # θ for molecule2 (H2)
phi2 = 0.0   # φ for molecule2 (tilt)
R = 2.5      # center-of-mass separation along +x

# scan limits
Rmin, Rmax = 1.0, 2.5
dR = -0.005

# animation toggles
scan_R = True
rotate_th1 = True
rotate_phi1 = True
rotate_th2 = True
rotate_phi2 = True
atom_approach = True

speed = 120

# ----------------- Scene -----------------
scene = canvas(width=1300, height=700, background=color.white)
scene.center = vector(0,0,0)
scene.forward = vector(-1, -0.25, -0.25)

soft_gray = vector(0.65,0.65,0.65)

def dotted_axis(start, direction, n=24, radius=0.008, col=soft_gray):
    for i in range(1, n+1):
        frac = i/n
        sphere(pos=start + frac*direction, radius=radius, color=col)
        sphere(pos=start - frac*direction, radius=radius, color=col)

# global dotted axes
dotted_axis(vector(0,0,0), vector(2.5,0,0), n=40)
dotted_axis(vector(0,0,0), vector(0,1,0), n=28)
dotted_axis(vector(0,0,0), vector(0,0,1), n=28)

label(pos=vector(2.6,0,0), text="+z", box=False, color=soft_gray)
label(pos=vector(0,1.1,0), text="+y", box=False, color=soft_gray)
label(pos=vector(0,0,1.1), text="+x", box=False, color=soft_gray)

# ----------------- Build molecule 1: H2O at origin -----------------
O_pos = vector(0,0,0)   # molecule1 COM ~ O atom (approx)
H1_local = vector(OH*cos(half_angle),  OH*sin(half_angle), 0)
H2_local = vector(OH*cos(half_angle), -OH*sin(half_angle), 0)

O1 = sphere(pos=O_pos, radius=0.12*scale, color=vector(0.2,0.2,0.8))
H1_1 = sphere(pos=O_pos + H1_local, radius=0.08*scale, color=color.white)
H1_2 = sphere(pos=O_pos + H2_local, radius=0.08*scale, color=color.white)
bond1a = cylinder(pos=O1.pos, axis=H1_1.pos - O1.pos, radius=0.02*scale, color=color.gray(0.5))
bond1b = cylinder(pos=O1.pos, axis=H1_2.pos - O1.pos, radius=0.02*scale, color=color.gray(0.5))

# body-frame arrows for molecule1
m1_x = arrow(pos=O1.pos, axis=vector(0.6*scale,0,0), shaftwidth=0.02, color=color.blue)
m1_y = arrow(pos=O1.pos, axis=vector(0,0.6*scale,0), shaftwidth=0.02, color=color.green)
m1_z = arrow(pos=O1.pos, axis=vector(0,0,0.6*scale), shaftwidth=0.02, color=color.orange)

# ----------------- Build molecule 2: H2 at +R along x -----------------
# Place H2 COM at x = R (we'll place two H atoms symmetrically around COM along local x)
def h2_local_positions(bond):
    half = bond/2.0
    return [vector(-half, 0, 0), vector(half, 0, 0)]

H2_local_positions = h2_local_positions(HH)

# initial COM for molecule2
COM2 = vector(R, 0, 0)
H2_1 = sphere(pos=COM2 + H2_local_positions[0], radius=0.09*scale, color=color.white)
H2_2 = sphere(pos=COM2 + H2_local_positions[1], radius=0.09*scale, color=color.white)
bond2 = cylinder(pos=H2_1.pos, axis=H2_2.pos - H2_1.pos, radius=0.02*scale, color=color.gray(0.4))

# body-frame arrows for molecule2 (attached to COM2)
m2_x = arrow(pos=COM2, axis=vector(0.5*scale,0,0), shaftwidth=0.02, color=color.blue)
m2_y = arrow(pos=COM2, axis=vector(0,0.5*scale,0), shaftwidth=0.02, color=color.green)
m2_z = arrow(pos=COM2, axis=vector(0,0,0.5*scale), shaftwidth=0.02, color=color.orange)

# incoming atom visual (optional: a separate single atom could be used; here molecule2 is the collider)
# we keep a small trail for molecule2 COM to help visualize motion
COM2_trail_ball = sphere(pos=COM2, radius=0.01, color=color.red, make_trail=True, retain=200)
COM2_trail_ball.visible = False  # hide point but keep trail if desired

# ----------------- R connector -----------------
R_line = cylinder(pos=O1.pos, axis=COM2 - O1.pos, radius=0.01*scale, color=color.black)

# ----------------- Labels (use line=False to avoid the leader line artifact) -----------------
R_label = label(pos=(O1.pos+COM2)/2 + vector(0,0.15*scale,0),
                text="R = {:.2f} Å".format(R),
                height=12, box=False, color=color.black, line=False)

th1_label = label(pos=O1.pos + vector(0.6*scale, 0.6*scale, 0),
                  text="th1 = {:.1f}°".format(degrees(th1)), height=14, box=False, color=color.red, line=False)

ph1_label = label(pos=O1.pos + vector(0.0, 0.9*scale, 0.5*scale),
                  text="ph1 = {:.1f}°".format(degrees(phi1)), height=14, box=False, color=color.green, line=False)

th2_label = label(pos=COM2 + vector(0.6*scale, 0.6*scale, 0),
                  text="th2 = {:.1f}°".format(degrees(th2)), height=14, box=False, color=color.red, line=False)

ph2_label = label(pos=COM2 + vector(0.0, 0.9*scale, 0.5*scale),
                  text="ph2 = {:.1f}°".format(degrees(phi2)), height=14, box=False, color=color.green, line=False)

# ----------------- Arcs for both molecules -----------------
# theta arcs (xy plane), phi arcs (xz plane), anchored at each COM
m1_theta_arc = curve(radius=0.012, color=color.red)
m1_phi_arc   = curve(radius=0.012, color=color.green)
m2_theta_arc = curve(radius=0.012, color=color.red)
m2_phi_arc   = curve(radius=0.012, color=color.green)

def make_theta_arc_at(center, ang, radius=0.6*scale, n=48):
    pts = []
    for i in range(n+1):
        a = ang * i / n
        pts.append(center + vector(radius*cos(a), radius*sin(a), 0))
    return pts

def make_phi_arc_at(center, ang, radius=0.5*scale, n=48):
    pts = []
    for i in range(n+1):
        a = ang * i / n
        pts.append(center + vector(radius*cos(a), 0, radius*sin(a)))
    return pts

# ----------------- Rotation helpers -----------------
def rot_z(v, ang):
    return vector(v.x*cos(ang) - v.y*sin(ang),
                  v.x*sin(ang) + v.y*cos(ang),
                  v.z)

def rot_y(v, ang):
    return vector(v.x*cos(ang) + v.z*sin(ang),
                  v.y,
                 -v.x*sin(ang) + v.z*cos(ang))

def apply_rotation(v, th, ph):
    v1 = rot_z(v, th)
    v2 = rot_y(v1, ph)
    return v2

# ----------------- Animation loop -----------------
frame = 0
while True:
    rate(speed)
    frame += 1

    # auto-rotate toggles
    if rotate_th1: th1 += 0.009
    if rotate_phi1: phi1 += 0.005
    if rotate_th2: th2 += 0.012
    if rotate_phi2: phi2 += 0.007

    # scan R
    if scan_R:
        R += dR
        if R < Rmin or R > Rmax:
            dR = -dR

    # update COM2 and its visuals
    COM2 = vector(R, 0, 0)
    # rotate H2 local positions by th2/phi2 then place around COM2
    H2_rot0 = apply_rotation(H2_local_positions[0], th2, phi2)
    H2_rot1 = apply_rotation(H2_local_positions[1], th2, phi2)
    H2_1.pos = COM2 + H2_rot0
    H2_2.pos = COM2 + H2_rot1
    bond2.pos = H2_1.pos
    bond2.axis = H2_2.pos - H2_1.pos

    # update molecule1 H positions by th1/phi1
    H1_1.pos = O1.pos + apply_rotation(H1_local, th1, phi1)
    H1_2.pos = O1.pos + apply_rotation(H2_local, th1, phi1)
    bond1a.pos = O1.pos
    bond1a.axis = H1_1.pos - O1.pos
    bond1b.pos = O1.pos
    bond1b.axis = H1_2.pos - O1.pos

    # update body arrows (frame) for both molecules
    m1_x.pos = O1.pos; m1_x.axis = apply_rotation(vector(0.6*scale,0,0), th1, phi1)
    m1_y.pos = O1.pos; m1_y.axis = apply_rotation(vector(0,0.6*scale,0), th1, phi1)
    m1_z.pos = O1.pos; m1_z.axis = apply_rotation(vector(0,0,0.6*scale), th1, phi1)

    m2_x.pos = COM2; m2_x.axis = apply_rotation(vector(0.5*scale,0,0), th2, phi2)
    m2_y.pos = COM2; m2_y.axis = apply_rotation(vector(0,0.5*scale,0), th2, phi2)
    m2_z.pos = COM2; m2_z.axis = apply_rotation(vector(0,0,0.5*scale), th2, phi2)

    # update R connector & label
    R_line.pos = O1.pos
    R_line.axis = COM2 - O1.pos
    R_label.pos = (O1.pos + COM2)/2 + vector(0, 0.13*scale, 0)
    R_label.text = "R = {:.2f} Å".format(R)

    # update labels (use simple ascii names to avoid parser artifacts; line=False to remove leader)
    th1_label.pos = O1.pos + vector(0.6*scale, 0.6*scale, 0)
    th1_label.text = "θ1 = {:.1f}°".format(degrees(th1)%360)
    ph1_label.pos = O1.pos + vector(0.0, 0.9*scale, 0.5*scale)
    ph1_label.text = "ϕ1 = {:.1f}°".format(degrees(phi1)%360)

    th2_label.pos = COM2 + vector(0.6*scale, 0.6*scale, 0)
    th2_label.text = "θ2 = {:.1f}°".format(degrees(th2)%360)
    ph2_label.pos = COM2 + vector(0.0, 0.9*scale, 0.5*scale)
    ph2_label.text = "ϕ2 = {:.1f}°".format(degrees(phi2)%360)

    # update arcs (anchored at COMs)
    m1_theta_pts = make_theta_arc_at(O1.pos, th1%(2*pi), radius=0.6*scale)
    m1_theta_arc.clear(); m1_theta_arc.append(pos=m1_theta_pts)
    m1_phi_pts = make_phi_arc_at(O1.pos, phi1%(2*pi), radius=0.5*scale)
    m1_phi_arc.clear(); m1_phi_arc.append(pos=m1_phi_pts)

    m2_theta_pts = make_theta_arc_at(COM2, th2%(2*pi), radius=0.45*scale)
    m2_theta_arc.clear(); m2_theta_arc.append(pos=m2_theta_pts)
    m2_phi_pts = make_phi_arc_at(COM2, phi2%(2*pi), radius=0.4*scale)
    m2_phi_arc.clear(); m2_phi_arc.append(pos=m2_phi_pts)

    # optional: small visible marker for COM2 (useful when R changes)
    COM2_trail_ball.pos = COM2
    COM2_trail_ball.visible = False

    # clear trails periodically
    if frame % 800 == 0:
        for obj in (H2_1, H2_2):
            obj.clear_trail()
