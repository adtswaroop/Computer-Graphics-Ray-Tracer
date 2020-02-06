Building a Ray Tracer.

![](Output%20Images/000.jpg)

Step 1: Uniformly send out rays from the camera location. 
Since the camera does not have to move, assume that its location is (0,0,0). 
Implemented backwards ray tracing where rays are sent from the camera, one ray per pixel. 

Step 2: Write the intersection code. 

Step 3: Implement the illumination equation
I = lightColor * (kd * (L dot N) + ks * (R dot V) ^ sh)      
L: light source vector
N: Normal Vector
R: Reflected Ray Vector
V: Eye position Vector

Step 4: Create still images showing off ray tracer.

