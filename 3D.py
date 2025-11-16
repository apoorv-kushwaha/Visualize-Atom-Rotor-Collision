from vpython import *
#Web VPython 3.2
# GlowScript 3.2 VPython
# H2O - Atom 3D collision visualiser with dotted axes, R label, and rotation arcs

from vpython import *

# ----------------- Parameters -----------------
OH = 0.5         # OH bond length (scaled smaller)
HOH_angle = radians(104.5)
half_angle = HOH_angle/2

theta = 0.0
phi   = 0.0
R     = 2.5
Rmin, Rmax = 1.0, 2.5
dR = -0.01

rotate_theta = True
rotate_phi   = True
scan_R       = True

# ----------------- Scene -----------------
scene = canvas(width=1500, height=700, background=color.white)
scene.center = vector(0,0,0)
scene.forward = vector(-1,-0.2,-0.4)

soft_gray = vector(0.6,0.6,0.6)

def dotted_axis(start, direction, n=20, radius=0.01, col=soft_gray):
    for i in range(1, n+1):
        frac = i/n
        posp = start + frac*direction
        posm = start - frac*direction
        sphere(pos=posp, radius=radius, color=col)
        sphere(pos=posm, radius=radius, color=col)

# Dotted global axes
dotted_axis(vector(0,0,0), vector(1,0,0), n=20)  # x
dotted_axis(vector(0,0,0), vector(0,1,0), n=20)  # y
dotted_axis(vector(0,0,0), vector(0,0,1), n=20)  # z

label(pos=vector(2.6,0,0), text="+z", box=False, color=soft_gray)
label(pos=vector(0,1.1,0), text="+y", box=False, color=soft_gray)
label(pos=vector(0,0,1.1), text="+x", box=False, color=soft_gray)

label(pos=vector(-1.1,0,0), text="-z", box=False, color=soft_gray)
label(pos=vector(0,-1.1,0), text="-y", box=False, color=soft_gray)
label(pos=vector(0,0,-1.1), text="-x", box=False, color=soft_gray)

# ----------------- Molecule -----------------
O_atom = sphere(pos=vector(0,0,0), radius=0.12, color=color.blue)
H1_local = vector(OH*cos(half_angle),  OH*sin(half_angle), 0)
H2_local = vector(OH*cos(half_angle), -OH*sin(half_angle), 0)

H1 = sphere(pos=H1_local, radius=0.08, color=color.white)
H2 = sphere(pos=H2_local, radius=0.08, color=color.white)

bond1 = cylinder(pos=O_atom.pos, axis=H1.pos - O_atom.pos, radius=0.02, color=color.gray(0.5))
bond2 = cylinder(pos=O_atom.pos, axis=H2.pos - O_atom.pos, radius=0.02, color=color.gray(0.5))

# Incoming atom
atom = sphere(pos=vector(R,0,0), radius=0.1, color=color.red, make_trail=True)

R_line = cylinder(pos=O_atom.pos, axis=atom.pos - O_atom.pos, radius=0.01, color=color.black)

# ----------------- Labels -----------------
R_label = label(pos=atom.pos/2, text="R = {:.2f}".format(R),
                xoffset=20, yoffset=20, height=12, box=False, color=color.black)

#theta_label = label(pos=vector(0.6,0.6,0), text="R = {:.2f}".format(R),
#                xoffset=20, yoffset=20, height=12, box=False, color=color.black)
                
theta_label = label(pos=vector(0.6,0.6,0), 
                    text="θ = {:.2f}°".format(degrees(theta)), 
                    height=20, box=False, color=color.red)

phi_label   = label(pos=vector(0,0.4,0.6), 
                    text=" ϕ = {:.2f}°".format(degrees(phi)), 
                    height=20, box=False, color=color.green)




# ----------------- Arcs -----------------
theta_arc = curve(radius=0.01, color=color.red)    # θ rotation around z
phi_arc   = curve(radius=0.01, color=color.green)  # φ rotation around y

def update_theta_arc(c, ang, radius=0.6, n=40):
    """θ arc in xy-plane around z-axis"""
    pts = []
    for i in range(n+1):
        frac = ang*i/n
        pts.append(vector(radius*cos(frac), radius*sin(frac), 0))
    c.clear()
    c.append(pos=pts)

def update_phi_arc(c, ang, radius=0.5, n=40):
    """φ arc in xz-plane around y-axis"""
    pts = []
    for i in range(n+1):
        frac = ang*i/n
        pts.append(vector(radius*cos(frac), 0, radius*sin(frac)))
    c.clear()
    c.append(pos=pts)
# ----------------- Rotation helper -----------------
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

# ----------------- Animation -----------------
while True:
    rate(60)

    if rotate_theta:
        theta += 0.02
    if rotate_phi:
        phi += 0.01

    if scan_R:
        R += dR
        if R < Rmin or R > Rmax:
            dR = -dR
        atom.pos = vector(R,0,0)

    # Rotate H positions
    H1_rot = apply_rotation(H1_local, theta, phi)
    H2_rot = apply_rotation(H2_local, theta, phi)

    H1.pos = O_atom.pos + H1_rot
    H2.pos = O_atom.pos + H2_rot

    bond1.pos = O_atom.pos
    bond1.axis = H1.pos - O_atom.pos
    bond2.pos = O_atom.pos
    bond2.axis = H2.pos - O_atom.pos

    # Update R line & label
    R_line.axis = atom.pos - O_atom.pos
    R_label.pos = atom.pos/2
    R_label.text = "R = {:.2f}".format(R)

    # Update θ & φ labels
    theta_label.text = "θ  = {:.1f}°".format(degrees(theta)%360)
    phi_label.text   = "ϕ = {:.1f}°".format(degrees(phi)%360)
    # Update arcs
    update_theta_arc(theta_arc, theta%(2*pi))
    update_phi_arc(phi_arc, phi%(2*pi))
