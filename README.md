# hw_1_glut_shader

First homework in 2110514 Realtime Computer Graphics and Physics Simulation at Chulalongkorn University,
almost the same as the one in CS184: Introduction to Computer Graphics - Spring 2006

## Compilation
In Mac OS, you could compile using this command
`g++ -w -O3 -o main main.cpp  -L/System/Library/Frameworks -framework GLUT -framework OpenGL`

## Usage
`./main [PARAMS]`, where PARAMS are set of,
* -ka r g b 
This is the ambient color coefficients of the sphere material. The parameters r g b are numbers between 0 and 1 inclusive. 
* -kd r g b 
This is the diffuse color coefficients of the sphere material. The parameters r g b are numbers between 0 and 1 inclusive. 
* -ks r g b 
This is the specular color coefficients of the sphere material. The parameters r g b are numbers between 0 and 1 inclusive. • -sp v 
This is the power coefficient on the specular term. It is a number between 0 and max_float. 
* -pl x y z r g b 
This adds a point light to the scene. The x y z values are the location of the light. The r g b values are it's color. Note that the x y z values are relative to the sphere. That is, the center of the sphere is at the origin and the radius of the sphere defines one unit of length. The Y direction is UP, the X direction is to the right on the screen, and the Z direction is "in your face." The r g b value are between 0 and max_float, NOT between 0 and 1 (that is, the r g b values encode the brightness of the light). 
* -dl x y z r g b 
This adds a directional light to the scene. The x y z values are the direction that the light points in. The r g b values are it's color. See -pl for coordinate system notes. There may be up to 5 point lights and 5 directional lights (10 total) in a scene. The r g b values of 1.0 should be mapped to a display values of 255. 
* -toon x
This will render the sphere with toon shader where x is the dividend for color variable (recommended start value: 100000)

Example: 
./main -ka 0.2 0.2 0.2 -kd 1 1 1 -ks 1 1 1 -sp 200 -pl 5 5 5 1 0 0 -pl -5 5 5 0 1 0 -pl -5 -5 5 0 0 1 -pl 0 1 0 1 1 1



