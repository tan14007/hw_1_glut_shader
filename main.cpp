// Simple OpenGL example for CS184 F06 by Nuttapong Chentanez, modified from sample code for CS184 on Sp06
// Modified for Realtime-CG class

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#define OSX

#ifdef _WIN32
#	include "windows.h"
#else
#	include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <time.h>
#include <math.h>
#define M_PI 3.14159265

#include "algebra3.h"

#ifdef _WIN32
static DWORD lastTime;
#else
static struct timeval lastTime;
#endif

#define PI 3.14159265
#define EPSILON 0.2
#define MINI_EPSILON 0.000001

using namespace std;

//****************************************************
// Some Classes
//****************************************************

class Viewport;

class Viewport {
public:
	int w, h; // width and height
};

class Material{
public:
	vec3 ka; // Ambient color
	vec3 kd; // Diffuse color
	vec3 ks; // Specular color
	float sp; // Power coefficient of specular
	bool toonShader; // Toon Shader flag
	float toonResolution;

	Material() : ka(0.0f), kd(0.0f), ks(0.0f), sp(0.0f) {
	}
};

class Light{
public:
	enum LIGHT_TYPE{POINT_LIGHT, DIRECTIONAL_LIGHT};

	vec3 posDir;  // Position (Point light) or Direction (Directional light)
	vec3 color;   // Color of the light
	LIGHT_TYPE type;

	Light() : posDir(0.0f), color(0.0f), type(POINT_LIGHT) {
	}
};

// Material and lights
Material material;
vector<Light> lights;

//****************************************************
// Global Variables
//****************************************************
Viewport	viewport;
int 		drawX = 0;
int 		drawY = 0;

void initScene(){
	if(material.toonShader)
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // White background for toon shader
	else
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f); // Black background for normal case

	glViewport (0,0,viewport.w,viewport.h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0,viewport.w, 0, viewport.h);
}

float min(float x, float y)
{
  return (x < y ? x : y);
}

float max(float x, float y)
{
  return (x > y ? x : y);
}


//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
	viewport.w = w;
	viewport.h = h;

	glViewport (0,0,viewport.w,viewport.h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, viewport.w, 0, viewport.h);

	drawX = (int)(viewport.w*0.5f);
	drawY = (int)(viewport.h*0.5f);


}

void setPixel(int x, int y, GLfloat r, GLfloat g, GLfloat b) {
	glColor3f(r, g, b);
	glVertex2f(x+0.5, y+0.5);
}

vec3 elementWiseMultiply(vec3 a, vec3 b) {
	vec3 c = vec3(0,0,0);

	c[0] = a[0] * b[0];
	c[1] = a[1] * b[1];
	c[2] = a[2] * b[2];
	return c;
}


vec3 computeShadedColor(vec3 pos) {

	// TODO: Your shading code mostly go here
	vec3 normal = vec3(pos);
	normal.normalize();
	vec3 color = vec3(0.0f,0.0f,0.0f); //Default black

	vec3 viewVector = vec3(0,0,1);
	
	for (int i=0; i<lights.size(); i++) {
		vec3 lightColor = vec3(lights[i].color);
		vec3 lightVector = vec3(0,0,0);
		
		if(lights[i].type == Light::POINT_LIGHT)
			lightVector = (lights[i].posDir) - pos;
		else
			lightVector = lights[i].posDir;

		lightVector.normalize();
		
		vec3 reflectionVector = -lightVector + 2*(lightVector*pos)*pos;
		
		reflectionVector.normalize();
		//Ambient term
		color += elementWiseMultiply(lightColor,material.ka);
		//Diffuse term
		color += max(elementWiseMultiply(material.kd,lightColor)*(lightVector*pos), 0.0f);
		//Specular term
		color += elementWiseMultiply(material.ks,lightColor)*pow(max(reflectionVector*viewVector, 0.0f), material.sp);
		
		if (material.toonShader && (viewVector*normal < EPSILON)) {
			color = vec3(0, 0, 0);
		}
	}

	// printf("Original: %f %f %f\n",color.r,color.g,color.b);

	if(material.toonShader){

		double rp = color.r;
		double gp = color.g;
		double bp = color.b;

		double cmax = max(rp, max(gp, bp));
		double cmin = min(rp, min(gp, bp));
		double lambda = cmax-cmin;

		// printf("C: %f %f %f\n",cmax,cmin,lambda);

		double h = 0;
		if(lambda < MINI_EPSILON) h = 0;
		else if (fabs(cmax - rp) < MINI_EPSILON) h = 60*(fmod((gp-bp)/lambda, 6.0f));
		else if (fabs(cmax - gp) < MINI_EPSILON) h = 60*(((bp-rp)/lambda) + 2);
		else if (fabs(cmax - bp) < MINI_EPSILON) h = 60*(((rp-gp)/lambda) + 4);

		h = fmod(h + 360.0f, 360.0f);

		// printf("Diff: %f %f %f\n", (gp-bp)/lambda, (bp-rp)/lambda, (rp-gp)/lambda);

		double l = (cmax+cmin)/2;

		double s = 0;
		if(cmax > MINI_EPSILON || l < 1) s = lambda/(1-fabs(2*l-1));

		
		// printf("HSV(Before): %f %f %f\n",h,s,v);

		s = floor(s/material.toonResolution*1000000)*material.toonResolution/1000000;
		// l = floor(l/material.toonResolution*1000000)*material.toonResolution/1000000;

		// printf("HSV(Divided): %f %f %f\n",h,s,v);

		double c = (1-fabs(2*l-1))*s;
		double x = c*((1-fabs(fmod(h/60,2.0f)-1)));
		double m = l-c/2;

		if(h >= 0 && h< 60){
			rp = c;
			gp = x;
			bp = 0;
		}
		else if(h >= 60 && h< 120){
			rp = x;
			gp = c;
			bp = 0;
		}
		else if(h >= 120 && h< 180){
			rp = 0;
			gp = c;
			bp = x;
		}
		else if(h >= 180 && h< 240){
			rp = 0;
			gp = x;
			bp = c;
		}
		else if(h >= 240 && h< 300){
			rp = x;
			gp = 0;
			bp = c;
		}
		else if(h >= 300 && h< 360){
			rp = c;
			gp = 0;
			bp = x;
		}

// 		if(fabs(rp+m - color.r) > MINI_EPSILON*10 || fabs(gp+m - color.g) > MINI_EPSILON*10 || fabs(bp+m - color.b) > MINI_EPSILON*10){
// 			printf("HSL: %lf %lf %lf\n", h, s, l);
// printf("Original: %lf %lf %lf\n",color.r,color.g,color.b);
// printf("Color Final: %lf %lf %lf\n", rp+m, gp+m, bp+m);
// 		}

		color.r = (rp+m);
		color.g = (gp+m);
		color.b = (bp+m);

		// printf("Color Final: %f %f %f\n", rp+m, gp+m, bp+m);
	}
	
	return color;
}

//****************************************************
// function that does the actual drawing of stuff
//***************************************************
void myDisplay() {

	glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer

	glMatrixMode(GL_MODELVIEW);					// indicate we are specifying camera transformations
	glLoadIdentity();							// make sure transformation is "zero'd"


	int drawRadius = min(viewport.w, viewport.h)/2 - 10;  // Make it almost fit the entire window
	float idrawRadius = 1.0f / drawRadius;
	// Start drawing sphere
	glBegin(GL_POINTS);

	for (int i = -drawRadius; i <= drawRadius; i++) {
		int width = floor(sqrt((float)(drawRadius*drawRadius-i*i)));
		for (int j = -width; j <= width; j++) {

			// Calculate the x, y, z of the surface of the sphere
			float x = j * idrawRadius;
			float y = i * idrawRadius;
			float z = sqrtf(1.0f - x*x - y*y);
			vec3 pos(x,y,z); // Position on the surface of the sphere

			vec3 col = computeShadedColor(pos);

			// Set the red pixel
			setPixel(drawX + j, drawY + i, col.r, col.g, col.b);
		}
	}
	glEnd();

	glFlush();
	glutSwapBuffers();					// swap buffers (we earlier set double buffer)
}


//****************************************************
// for updating the position of the circle
//****************************************************

void myFrameMove() {
	float dt;
	// Compute the time elapsed since the last time the scence is redrawn
#ifdef _WIN32
	DWORD currentTime = GetTickCount();
	dt = (float)(currentTime - lastTime)*0.001f;
#else
	timeval currentTime;
	gettimeofday(&currentTime, NULL);
	dt = (float)((currentTime.tv_sec - lastTime.tv_sec) + 1e-6*(currentTime.tv_usec - lastTime.tv_usec));
#endif

	// Store the time
	lastTime = currentTime;
	glutPostRedisplay();
}


void parseArguments(int argc, char* argv[]) {
	int i = 1;
	while (i < argc) {
		if (strcmp(argv[i], "-ka") == 0) {
			// Ambient color
			material.ka.r = (float)atof(argv[i+1]);
			material.ka.g = (float)atof(argv[i+2]);
			material.ka.b = (float)atof(argv[i+3]);
			i+=4;
		} else
		if (strcmp(argv[i], "-kd") == 0) {
			// Diffuse color
			material.kd.r = (float)atof(argv[i+1]);
			material.kd.g = (float)atof(argv[i+2]);
			material.kd.b = (float)atof(argv[i+3]);
			i+=4;
		} else
		if (strcmp(argv[i], "-ks") == 0) {
			// Specular color
			material.ks.r = (float)atof(argv[i+1]);
			material.ks.g = (float)atof(argv[i+2]);
			material.ks.b = (float)atof(argv[i+3]);
			i+=4;
		} else
		if (strcmp(argv[i], "-sp") == 0) {
			// Specular power
			material.sp = (float)atof(argv[i+1]);
			i+=2;
		} else
		if ((strcmp(argv[i], "-pl") == 0) || (strcmp(argv[i], "-dl") == 0)){
			Light light;
			// Specular color
			light.posDir.x = (float)atof(argv[i+1]);
			light.posDir.y = (float)atof(argv[i+2]);
			light.posDir.z = (float)atof(argv[i+3]);
			light.color.r = (float)atof(argv[i+4]);
			light.color.g = (float)atof(argv[i+5]);
			light.color.b = (float)atof(argv[i+6]);
			if (strcmp(argv[i], "-pl") == 0) {
				// Point
				light.type = Light::POINT_LIGHT;
			} else {
				// Directional
				light.type = Light::DIRECTIONAL_LIGHT;
			}
			lights.push_back(light);
			i+=7;
		} else
		if (strcmp(argv[i], "-toon") == 0) {
			// Toon Shader
			material.toonShader = true;
			material.toonResolution = (float)atof(argv[i+1]);
			i+=2;
		} 
	}
}

//****************************************************
// the usual stuff, nothing exciting here
//****************************************************
int main(int argc, char *argv[]) {

	parseArguments(argc, argv);

  	//This initializes glut
  	glutInit(&argc, argv);

  	//This tells glut to use a double-buffered window with red, green, and blue channels
  	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

  	// Initalize theviewport size
  	viewport.w = 400;
  	viewport.h = 400;

  	//The size and position of the window
  	glutInitWindowSize(viewport.w, viewport.h);
  	glutInitWindowPosition(0,0);
  	glutCreateWindow(argv[0]);

   	// Initialize timer variable
	#ifdef _WIN32
	lastTime = GetTickCount();
	#else
	gettimeofday(&lastTime, NULL);
	#endif

  	initScene();							// quick function to set up scene

  	glutDisplayFunc(myDisplay);					// function to run when its time to draw something
  	glutReshapeFunc(myReshape);					// function to run when the window gets resized
  	glutIdleFunc(myFrameMove);

  	glutMainLoop();							// infinite loop that will keep drawing and resizing and whatever else

  	return 0;
}








