from vpython import *
#Web VPython 3.2
from vpython import *


######################## --- Input Parameters ---###############################
L1 = 1.5   # Rotor 1 bond length
L2 = 0.7   # Rotor 2 bond length
R  = 2.0   # Distance between COM1 and COM2

# Initial angles (radians)
th1_init = 0
th2_init = 0
phi_init = 0   # tilt

# Enable/disable rotations
rotate_th1 = True
rotate_th2 = True
rotate_phi = False
scan_R     = False   

# R scanning setup
R = 5.0
R_min, R_max = 1.0, 5.0
dR = -0.01   # step size (negative for decreasing R)

speed = 60  # FPS
############################## ---  END ---#####################################




# --- Canvas ---
#scene = canvas(width=700, height=700, background=color.white)
#scene.camera.pos = vector(0,0,5)
#scene.userspin = False           # no rotation
# Camera setup


# color define
soft_gray = vector(0.7, 0.7, 0.7)   # light grey, closer to neutral

# --- Generate heading based on initial values ---
heading_text = f"R scan ({R_max}-{R_min}) | " if scan_R else f"R fixed ({R}) | "
heading_text += "θ₁ rotating | " if rotate_th1 else f"θ₁ fixed ({th1_init* 180 / pi}°) | "
heading_text += "θ₂ rotating | " if rotate_th2 else f"θ₂ fixed ({th2_init* 180 / pi}°)| "
heading_text += "ϕ rotating" if rotate_phi else f"ϕ fixed ({phi_init* 180 / pi}°)"

# --- Canvas with heading ---
scene = canvas(
    width=1500, height=700, background=color.white,
    caption=f"<b>{heading_text}</b>\n\n"
)
scene.center = vector(1,0,0)   # focus point for rotation (default 0,0,0)
scene.userspin = True    # Enable spin/orbit
scene.userzoom = True    # keep zoom with right-drag
scene.userpan  = True    # Allow pan with shift+drag

# --- COM positions ---
COM1 = vector(0, 0, 0)
COM2 = vector(R, 0, 0)

label(pos=COM1, text="COM₁", xoffset=20, yoffset=-20, box=False, height=12)
COM2_label = label(pos=COM2, text="COM₂", xoffset=-20, yoffset=-20, box=False, height=12)

# --- Reference axes setup ---
axis_len = 1.5
x_axis_len = R + 2

def dotted_axis(start, direction, color, n_dots=30, radius=0.01):
    for sign in [+1, -1]:
        for i in range(1, n_dots+1):
            pos = start + sign*(i/n_dots)*direction
            sphere(pos=pos, radius=radius, color=color)

# --- Global x-axis ---
dotted_axis((COM1+COM2)/2, vector(x_axis_len/2,0,0), soft_gray, n_dots=40, radius=0.008)
label(pos=COM2 + vector(1.25,0,0), text="+z", box=False, color=soft_gray)
label(pos=COM1 - vector(1.25,0,0), text="−z", box=False, color=soft_gray)

# --- Local y, z axes at COM1 ---
dotted_axis(COM1, vector(0,axis_len,0), soft_gray, radius=0.008)
dotted_axis(COM1, vector(0,0,axis_len), soft_gray, radius=0.008)
label(pos=COM1 + vector(0,axis_len,0), text="+y₁", box=False, color=soft_gray)
label(pos=COM1 - vector(0,axis_len,0), text="−y₁", box=False, color=soft_gray)
label(pos=COM1 + vector(0,0,axis_len), text="+x₁", box=False, color=soft_gray)
label(pos=COM1 - vector(0,0,axis_len), text="−x₁", box=False, color=soft_gray)

# --- Local y, z axes at COM2 ---
y2_pos = vector(0, axis_len, 0)
z2_pos = vector(0, 0, axis_len)

# store dotted axis points as curve objects (instead of spheres, easier to move)
def dotted_axis_objects(start, direction, color, n_dots=30, radius=0.01):
    """Return a list of spheres representing a dotted axis."""
    objs = []
    for sign in [+1, -1]:
        for i in range(1, n_dots+1):
            pos = start + sign*(i/n_dots)*direction
            objs.append(sphere(pos=pos, radius=radius, color=color))
    return objs

# Local y, z axes at COM2 (as object lists)
axis_len = 1.5
y2_axis_objs = dotted_axis_objects(COM2, vector(0, axis_len, 0), soft_gray, radius=0.008)
z2_axis_objs = dotted_axis_objects(COM2, vector(0, 0, axis_len), soft_gray, radius=0.008)

# Labels (stored for later updates)
y2_label_p = label(pos=COM2 + vector(0, axis_len, 0), text="+y₂", box=False, color=soft_gray)
y2_label_m = label(pos=COM2 - vector(0, axis_len, 0), text="−y₂", box=False, color=soft_gray)
x2_label_p = label(pos=COM2 + vector(0, 0, axis_len), text="+x₂", box=False, color=soft_gray)
x2_label_m = label(pos=COM2 - vector(0, 0, axis_len), text="−x₂", box=False, color=soft_gray)


# --- Rotor 1 setup ---
r1_atoms = [vector(-L1/2, 0, 0), vector(L1/2, 0, 0)]
rotor1 = [sphere(pos=COM1 + r, radius=0.1, color=color.blue) for r in r1_atoms]
bond1  = cylinder(pos=rotor1[0].pos, axis=rotor1[1].pos-rotor1[0].pos, radius=0.02, color=color.blue)

# --- Rotor 2 setup ---
r2_atoms = [vector(-L2/2, 0, 0), vector(L2/2, 0, 0)]
rotor2 = [sphere(pos=COM2 + r, radius=0.08, color=color.green) for r in r2_atoms]
bond2  = cylinder(pos=rotor2[0].pos, axis=rotor2[1].pos-rotor2[0].pos, radius=0.015, color=color.red)

# --- Angle arcs and labels ---
arc1 = curve(color=color.yellow, radius=0.02)
arc2 = curve(color=color.orange, radius=0.02)
phi_arc = curve(color=color.purple, radius=0.02)
R_line = cylinder(pos=COM1, axis=COM2-COM1, radius=0.01, color=vector(0.5,0.5,0.5))

theta1_label = label(pos=vector(0,0,0), text="θ₁", box=False, color=color.blue)
theta2_label = label(pos=COM2, text="θ₂", box=False, color=color.red)
phi_label    = label(pos=COM2, text="ϕ",  box=False, color=color.purple)

R_label = label(pos=(COM1+COM2)/2, text="R", xoffset=0, yoffset=20, box=False, color=color.black)

# --- Simulation variables ---
th1 = th1_init
th2 = th2_init
phi = phi_init

circle_path = curve(color=color.cyan, radius=0.01)

while True:
    rate(speed)

    # --- Scan R if enabled ---
    if scan_R:
        R += dR
        if R <= R_min or R >= R_max:
            dR *= -1   # reverse direction

        COM2 = vector(R, 0, 0)  # update COM2

        # --- Update dotted axes at COM2 ---
        def update_dotted_axis(objs, origin, direction):
            n = len(objs)//2
            for i in range(n):
                frac = (i+1)/n
                objs[i].pos      = origin + frac*direction      # positive side
                objs[i+n].pos    = origin - frac*direction      # negative side
        
        update_dotted_axis(y2_axis_objs, COM2, vector(0, axis_len, 0))
        update_dotted_axis(z2_axis_objs, COM2, vector(0, 0, axis_len))
        
        # --- Update labels ---
        y2_label_p.pos = COM2 + vector(0, axis_len, 0)
        y2_label_m.pos = COM2 - vector(0, axis_len, 0)
        x2_label_p.pos = COM2 + vector(0, 0, axis_len)
        x2_label_m.pos = COM2 - vector(0, 0, axis_len)

        R_line.axis = COM2 - COM1
        R_label.pos = (COM1 + COM2)/2
        
        # --- When COM2 changes ---
        COM2_label.pos = COM2
        
        
    if rotate_th1: th1 += 0.02
    if rotate_th2: th2 += 0.02
    if rotate_phi: phi += 0.02

    # Rotor 1 rotation
    for i, r in enumerate(r1_atoms):
        x = r.x*cos(th1) - r.y*sin(th1)
        y = r.x*sin(th1) + r.y*cos(th1)
        rotor1[i].pos = COM1 + vector(x, y, r.z)
    bond1.pos = rotor1[0].pos
    bond1.axis = rotor1[1].pos - rotor1[0].pos

    # Rotor 2 rotation (now relative to updated COM2)
    for i, r in enumerate(r2_atoms):
        x = r.x*cos(th2) - r.y*sin(th2)
        y = r.x*sin(th2) + r.y*cos(th2)
        z = r.z
        y_new = y*cos(phi) - z*sin(phi)
        z_new = y*sin(phi) + z*cos(phi)
        rotor2[i].pos = COM2 + vector(x, y_new, z_new)
    bond2.pos = rotor2[0].pos
    bond2.axis = rotor2[1].pos - rotor2[0].pos

    # --- Circle path update ---
    circle_radius = L2/2
    circle_points = []
    N = 100
    
    circle_path.clear()
    for i in range(N+1):
        t = 2*pi*i/N
        x = circle_radius*cos(t)
        y = circle_radius*sin(t)*cos(phi)
        z = circle_radius*sin(t)*sin(phi)
        circle_points.append(COM2 + vector(x,y,z))
    circle_path = curve(pos=circle_points, color=color.cyan, radius=0.01)

    # --- Update arcs (same as before) ---
    n_pts = 20
    arc1.clear()
    for i in range(n_pts):
        t = (th1 % (2*pi)) * i/n_pts
        arc1.append(COM1 + vector(0.3*cos(t), 0.3*sin(t), 0))
    theta1_label.pos = COM1 + vector(1*cos(th1), 1*sin(th1), 0)

    arc2.clear()
    for i in range(n_pts):
        t = (th2 % (2*pi)) * i/n_pts
        arc2.append(COM2 + vector(0.3*cos(t), 0.3*sin(t), 0))
    theta2_label.pos = COM2 + vector(0.6*cos(th2), 0.6*sin(th2), 0)

    phi_arc.clear()
    for i in range(n_pts):
        t = (phi % (2*pi)) * i/n_pts
        phi_arc.append(COM2 + vector(0.0, 0.3*cos(t), 0.3*sin(t)))
    phi_label.pos = COM2 + vector(0.3, 0.3*cos(phi), 0.3*sin(phi))

