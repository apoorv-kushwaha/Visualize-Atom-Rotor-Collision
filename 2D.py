from vpython import *
#Web VPython 3.2
from vpython import *

# --------------------------------------------------------------------
# LAYOUT: Two canvases side-by-side + centered text
# --------------------------------------------------------------------
scene1 = canvas(width=700, height=200, background=color.white)
scene2 = canvas(width=700, height=500, background=color.white)
scene3 = canvas(width=700, height=50, background=color.white)

# --------------------------------------------------------------------
# 1) LEFT ANIMATION — ROTOR ROTATING (Blue Diatom + Angle)
# --------------------------------------------------------------------
R = 2

# rotor atoms
r_atoms = [vector(-0.5,0,0), vector(0.5,0,0)]
rotor = [sphere(canvas=scene1, pos=r, radius=0.1, color=color.blue) for r in r_atoms]

# fixed red atom
fixed_atom = vector(R,0,0)
atom_fixed = sphere(canvas=scene1, pos=fixed_atom, radius=0.12, color=color.red)

# labels
label(canvas=scene1, pos=fixed_atom, text="Atom", xoffset=20, yoffset=20, box=False)
label(canvas=scene1, pos=vector(0,0,0), text="COM", xoffset=-30, yoffset=-30, box=False)

# double arrow R
shaft1 = cylinder(canvas=scene1, pos=vector(0,0,0), axis=fixed_atom, radius=0.01)
arrow1a = cone(canvas=scene1, pos=fixed_atom*0.98, axis=fixed_atom*0.02, radius=0.05)
arrow1b = cone(canvas=scene1, pos=vector(0,0,0), axis=-fixed_atom*0.02, radius=0.05)
label(canvas=scene1, pos=fixed_atom/2, text="R", box=False, height=12)

# bond between rotor atoms
bond1 = cylinder(canvas=scene1, pos=r_atoms[0], axis=r_atoms[1]-r_atoms[0], radius=0.01, color=color.blue)

# angle arc
arc1 = curve(canvas=scene1, color=color.yellow, radius=0.02)
theta_label1 = label(canvas=scene1, pos=vector(0.7,0.3,0), text="θ", box=False, height=18, opacity=0)

# --------------------------------------------------------------------
# 2) RIGHT ANIMATION — RED ATOM ORBITING (COM fixed)
# --------------------------------------------------------------------
# rotor fixed
rotor2 = [sphere(canvas=scene2, pos=r, radius=0.1, color=color.blue) for r in r_atoms]
bond2 = cylinder(canvas=scene2, pos=r_atoms[0], axis=r_atoms[1]-r_atoms[0], radius=0.01, color=color.blue)

# orbiting red atom
atom2 = sphere(canvas=scene2, pos=vector(R,0,0), radius=0.12, color=color.red)
#label(canvas=scene2, pos=atom2.pos, text="Atom", xoffset=20, yoffset=20, box=False)
label(canvas=scene2, pos=vector(0,0,0), text="COM", xoffset=-30, yoffset=-30, box=False)

R_label2 = label(canvas=scene2, pos=atom2.pos/2, text="R", box=False, height=12)

# R arrow
shaft2 = cylinder(canvas=scene2, pos=vector(0,0,0), axis=atom2.pos, radius=0.01)
arrow2a = cone(canvas=scene2, pos=atom2.pos*0.98, axis=atom2.pos*0.02, radius=0.05)
arrow2b = cone(canvas=scene2, pos=vector(0,0,0), axis=-atom2.pos*0.02, radius=0.05)

# arc for theta
arc2 = curve(canvas=scene2, color=color.yellow, radius=0.02)
theta_label2 = label(canvas=scene2, pos=vector(0.7,0.3,0), text="θ", box=False, height=18, opacity=0)
arcR = 0.6


# --------------------------------------------------------------------
# 1) Center ANIMATION 
# --------------------------------------------------------------------

# Center text below both
label(canvas=scene3,pos=vector(0,1,0), text="Both pictures are equivalent",
      height=20, box=False, color=color.black)


# --------------------------------------------------------------------
# RUN BOTH ANIMATIONS
# --------------------------------------------------------------------
theta = 0
while True:
    rate(24)
    theta += 0.03
    # -----------------------------------
    # LEFT animation: Rotate rotor atoms
    # -----------------------------------
    rpos = []
    for r in r_atoms:
        newpos = vector(r.x*cos(theta) - r.y*sin(theta),
                        r.x*sin(theta) + r.y*cos(theta), 0)
        rpos.append(newpos)
    rotor[0].pos = rpos[0]
    rotor[1].pos = rpos[1]
    bond1.pos = rpos[0]
    bond1.axis = rpos[1] - rpos[0]

    # arc update
    th = theta % (2*pi)
    arc1.clear()
    for t in [i*th/25 for i in range(26)]:
        arc1.append(vector(0.6*cos(t), 0.6*sin(t), 0))
    theta_label1.pos = vector(0.75*cos(th*0.95), 0.75*sin(th*0.95), 0)
    # Update θ label at middle of arc        
    # Update blue bond line to follow atoms
    # -----------------------------------
    # RIGHT animation: Orbiting atom
    # -----------------------------------
    atom2.pos = vector(R*cos(theta), R*sin(theta), 0)
    R_label2.pos = atom2.pos/2

    shaft2.axis = atom2.pos
    arrow2a.pos = atom2.pos*0.98
    arrow2a.axis = atom2.pos*0.02
    arrow2b.axis = -atom2.pos*0.02

    # short arc
    arc2.clear()
    for t in [i*th/25 for i in range(26)]:
        arc2.append(vector(arcR*cos(t), arcR*sin(t), 0))
        
    theta_label2.pos = vector(arcR*1.3*cos(th*0.95), arcR*1.3*sin(th*0.95), 0)
