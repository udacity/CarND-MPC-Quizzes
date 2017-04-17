1. Forces Quiz (Type Radio):

Question:

What forces does a dynamic model take into account?

A [ ] Tire forces, longitudinal and lateral forces, gravity, air resistance, drag.

B [ ] Mass, vehicle width, vehicle length and other vehicle characteristics.

C [ ] A and B.

Answer:

A

The answers in B aren't forces but vehicle characteristics. Kinematic models often take these characteristics into account but forces and accelerations are specific to dynamic models.

If they selected B

Are these forces? Could these be in a kinematic model?

if they selected C

Are the vehicle characteristics forces that affect the motion of the vehicle?

2. Kinematic Dynamics Quiz (Might change this to the Rajamani equations)

Question (Type Calculate Stuff):

x_{k+1} = x_{k} + v_{k} * cos(psi_{k} + sa_{k}) * dt
y_{k+1} = y_{k} + v_{k} * sin(psi_{k} + sa_{k}) * dt
psi_{k+1} = psi_{k} + v_{k} * sin(sa_{k}) * dt
v_{k+1} = v_{k} + a_{k} * dt

Given certain values ? calculate the future values.

Answer:

whatever stuff

3. Error minimization quiz (Type Checkbox)

Question:

What errors do we want to minimize?

A [] Distance of vehicle from trajectory

B [] Difference of vehicle orientation and trajectory orientation

C [] The speed of the vehicle

Answer:

A and B

If A is checked. While we might want to set the speed to a specific value, it's not something we'd like to minimize.

4. Acutator Constraints Quiz

Q: Why would we want to constrain actuators?

[ ] Adhere to physics.

[ ] Mimic limitations of test vehicle.

A:

Both (A and B)

5. Curvature Quiz

TODO: Flesh out this quiz more

Q: When will curvature give provide a significant advantage over straight line assumptions?

[ ] When travelling at high speeds.

[ ] 

[ ] Always.

A: A or C

C is correct, but, in some situations curvature offers more benefit than others.

