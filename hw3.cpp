/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: <Your name here>
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include <glm/gtc/type_ptr.hpp>

#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

#define EPSILON -1e-5
int counter = 10000;

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

void getVertex(Triangle t, glm::vec3 *triangleA, glm::vec3 *triangleB, glm::vec3 *triangleC)
{
	*triangleA = { t.v[0].position[0], t.v[0].position[1], t.v[0].position[2] };
	*triangleB = { t.v[1].position[0], t.v[1].position[1], t.v[1].position[2] };
	*triangleC = { t.v[2].position[0], t.v[2].position[1], t.v[2].position[2] };
}

void getEdges(Triangle t, glm::vec3 *BminusA, glm::vec3 *CminusA, glm::vec3 *CminusB, glm::vec3 *AminusC)
{
	glm::vec3 triangleA, triangleB, triangleC;
	getVertex(t, &triangleA, &triangleB, &triangleC);

	*BminusA = { triangleB.x - triangleA.x , triangleB.y - triangleA.y , triangleB.z - triangleA.z };
	*CminusA = { triangleC.x - triangleA.x , triangleC.y - triangleA.y , triangleC.z - triangleA.z };
	*CminusB = { triangleC.x - triangleB.x , triangleC.y - triangleB.y , triangleC.z - triangleB.z };
	*AminusC = { triangleA.x - triangleC.x , triangleA.y - triangleC.y , triangleA.z - triangleC.z };
}

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;
int num_rays = 0;

bool supersample = false, softshadows = false;

double t0, t1, t;
double alpha, beta, gamma;
double r = 0.0, g = 0.0, b = 0.0;
double tanfov = std::tan((fov / 2.0) * (3.1415 / 180.0));
double aspectRatio = (double)WIDTH / (double)HEIGHT;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

double clampValues(double value)
{
	double clampedVal = NULL;
	if (value > 1.0)
		clampedVal = 1.0;
	else if (value < 0.0)
		clampedVal = 0.0;
	else
		clampedVal = value;
	return clampedVal;
}


void computeBarrycentricCoords(Triangle triangle, glm::vec3 inPt, double *alpha, double *beta)
{
	glm::vec3 triangleA, triangleB, triangleC;
	getVertex(triangle, &triangleA, &triangleB, &triangleC);
	
	glm::vec3  BminusA, CminusA, CminusB, AminusC;
	getEdges(triangle, &BminusA, &CminusA, &CminusB, &AminusC);

	//populate the normal (B-A) x (C-A)
	glm::vec3 n = glm::cross(BminusA, CminusA);

	glm::vec3 nX = { 1.0, 0.0, 0.0 };
	glm::vec3 nY = { 0.0, 1.0, 0.0 };
	glm::vec3 nZ = { 0.0, 0.0, 1.0 };

	double XY = glm::dot(n, nZ);
	double YZ = glm::dot(n, nX);
	double ZX = glm::dot(n, nY);

	glm::vec2 C, C0, C1, C2;
	double areaCC0C1, areaCC1C2, areaCC0C2, areaC0C1C2;

	//find the plane
	if (XY > ZX && XY > YZ) {
		//projection Plane XY
		C = { inPt.x, inPt.y };
		C0 = { triangle.v[0].position[0], triangle.v[0].position[1] };
		C1 = { triangle.v[1].position[0], triangle.v[1].position[1] };
		C2 = { triangle.v[2].position[0], triangle.v[2].position[1] };

	}
	else if (ZX > YZ && ZX > XY)
	{
		//projection Plane ZX
		C = { inPt.x, inPt.y };
		C0 = { triangle.v[0].position[2], triangle.v[0].position[0] };
		C1 = { triangle.v[1].position[2], triangle.v[1].position[0] };
		C2 = { triangle.v[2].position[2], triangle.v[2].position[0] };
	}
	else
	{
		//projection Plane YZ
		C = { inPt.x, inPt.y };
		C0 = { triangle.v[0].position[1], triangle.v[0].position[2] };
		C1 = { triangle.v[1].position[1], triangle.v[1].position[2] };
		C2 = { triangle.v[2].position[1], triangle.v[2].position[2] };
	}

	//area = (1/2) ((bx – ax)(cy – ay) – (cx – ax) (by – ay))
	areaCC0C1 = 0.5 * (((C1.x - C0.x) * (C.y - C0.y)) - ((C.x - C0.x) * (C1.y - C0.y)));
	areaCC1C2 = 0.5 * (((C1.x - C.x)  * (C2.y - C.y)) - ((C2.x - C.x) * (C1.y - C.y)));
	areaCC0C2 = 0.5 * (((C.x - C0.x)  * (C2.y - C0.y)) - ((C2.x - C0.x) * (C.y - C0.y)));
	areaC0C1C2 = 0.5 * (((C1.x - C0.x) * (C2.y - C0.y)) - ((C2.x - C0.x) * (C1.y - C0.y)));

	*alpha = areaCC1C2 / areaC0C1C2;
	*beta = areaCC0C2 / areaC0C1C2;
}


class Ray
{
private:
	glm::vec3 rayOrigin;
	glm::vec3 rayDirection;

public:
	Ray() {}
	Ray(glm::vec3 o, glm::vec3 d) : rayOrigin(o), rayDirection(d) {}

	bool sphereIntersectionTest(const Sphere &sphere, glm::vec3 &iPt)
	{
		glm::vec3 spherePosition = { sphere.position[0], sphere.position[1], sphere.position[2] };
		//glm::vec3 raySphereDistance = rayOrigin - spherePosition;

		double a = glm::dot(rayDirection, rayDirection);
		double b = 2 * (rayDirection.x * (rayOrigin.x - spherePosition.x) +
			rayDirection.y * (rayOrigin.y - spherePosition.y) +
			rayDirection.z * (rayOrigin.z - spherePosition.z));
		double c = (rayOrigin.x - spherePosition.x) * (rayOrigin.x - spherePosition.x) +
			(rayOrigin.y - spherePosition.y) * (rayOrigin.y - spherePosition.y) +
			(rayOrigin.z - spherePosition.z) * (rayOrigin.z - spherePosition.z) -
			(sphere.radius * sphere.radius);

		double discriminant = b * b - 4 * c;

		//if discriminant < 0 then no intersection
		if (discriminant < 0)
			return false;

		//if discriminant = 0 then t0 = t1
		else if (discriminant == 0)
		{
			t0 = -b / 2.0;
			t1 = t0;
		}

		//if discriminant > 0 then calculate t0 and t1
		else
		{
			t0 = (-b + std::sqrt(discriminant)) / 2.0;
			t1 = (-b - std::sqrt(discriminant)) / 2.0;
		}

		//check if both t0 and t1 > 0, find the minimum
		if (t0 > 0 && t1 > 0)
		{
			if (t1 < t0)
				t = t1;
			else
				t = t0;
		}

		else if (t0 > 0 && t1 < 0)
		{
			t = t0;
		}
		else if (t0 < 0 && t1 > 0)
		{
			t = t1;
		}
		//if both are negative, no intersection
		else
			return false;

		//calculate the intersection point using p(t) = p0 + dt
		iPt = { rayOrigin.x + (rayDirection.x * t),
				rayOrigin.y + (rayDirection.y * t),
				rayOrigin.z + (rayDirection.z * t) };

		return true;
	}

	bool triangleIntersectionTest(const Triangle &triangle, glm::vec3 &iPt)
	{
		glm::vec3 triangleA, triangleB, triangleC;
		getVertex(triangle, &triangleA, &triangleB, &triangleC);

		glm::vec3  BminusA, CminusA, CminusB, AminusC;
		getEdges(triangle, &BminusA, &CminusA, &CminusB, &AminusC);

		//populate the normal (B-A) x (C-A)
		glm::vec3 n = glm::cross(BminusA, CminusA);

		glm::vec3 normalised_n = glm::normalize(n);

		//check if ray and normal are parallel 
		double nDotd = glm::dot(normalised_n, rayDirection);
		if (std::abs(nDotd) < 1e-5)
			return false;

		//calculate implicit function cofficients
		double a = normalised_n.x;
		double b = normalised_n.y;
		double c = normalised_n.z;
		double d = -(a * triangleA.x) - (b * triangleA.y) - (c * triangleA.z);
		
		double lowerVal = a * rayDirection.x + b * rayDirection.y + c * rayDirection.z;
		
		if (lowerVal == 0)
			return false;
		
		glm::vec3 vectorOtoA = rayOrigin - triangleA;
		double t = -(glm::dot(vectorOtoA, normalised_n) / nDotd);

		if (t <= 0.001)
			return false;

		iPt = { rayOrigin.x + rayDirection.x * t,
				rayOrigin.y + rayDirection.y * t,
				rayOrigin.z + rayDirection.z * t };

		glm::vec3 vector_ipt_A = { iPt.x - triangleA.x , iPt.y - triangleA.y, iPt.z - triangleA.z };
		glm::vec3 vector_ipt_B = { iPt.x - triangleB.x , iPt.y - triangleB.y, iPt.z - triangleB.z };
		glm::vec3 vector_ipt_C = { iPt.x - triangleC.x , iPt.y - triangleC.y, iPt.z - triangleC.z };

		glm::vec3 cross_ipt_CA = glm::cross(AminusC, vector_ipt_C);
		glm::vec3 cross_ipt_CB = glm::cross(CminusB, vector_ipt_B);
		glm::vec3 cross_ipt_BA = glm::cross(BminusA, vector_ipt_A);;

		double val1 = glm::dot(cross_ipt_BA, n);
		double val2 = glm::dot(cross_ipt_CB, n);
		double val3 = glm::dot(cross_ipt_CA, n);

		if ((glm::dot(cross_ipt_BA, n) < 0) ||
			(glm::dot(cross_ipt_CB, n) < 0) ||
			(glm::dot(cross_ipt_CA, n) < 0))
			return false;

		return true;

	}
};


struct Color
{
	double r, g, b;
	Color() { r = 0.0;  g = 0.0; b = 0.0; }
	Color(double rI, double gI, double bI) {
		r = rI; g = gI; b = bI;
	}

};

Ray* generateRays(double x, double y, int num_rays)
{
	Ray* rays = new Ray[num_rays];

	if(num_rays == 1)
	{
		double screenX = (2 * ((x + 0.5) / WIDTH)) - 1;
		double screenY = (2 * ((y + 0.5) / HEIGHT)) - 1;
		double cameraX = aspectRatio * tanfov * screenX;
		double cameraY = tanfov * screenY;

		glm::vec3 genRay = { cameraX , cameraY , -1 };

		glm::vec3 dir = glm::normalize(genRay);

		glm::vec3 rayOrigin = { 0.0, 0.0, 0.0 };
		Ray generatedRay(rayOrigin, dir);

		rays[0] = generatedRay;
		return rays;
	}

	else {
		if (num_rays > 0)
		{
			//Ray 1
			double screenX = (2 * ((x + 0.25) / WIDTH)) - 1;
			double screenY = (2 * ((y + 0.25) / HEIGHT)) - 1;
			double cameraX = aspectRatio * tanfov * screenX;
			double cameraY = tanfov * screenY;
			glm::vec3 genRay = { cameraX , cameraY , -1 };
			glm::vec3 dir = glm::normalize(genRay);
			glm::vec3 rayOrigin = { 0.0, 0.0, 0.0 };
			Ray ray1(rayOrigin, dir);
			rays[0] = ray1;

			//Ray 2
			screenX = (2 * ((x + 0.50) / WIDTH)) - 1;
			screenY = (2 * ((y + 0.25) / HEIGHT)) - 1;
			cameraX = aspectRatio * tanfov * screenX;
			cameraY = tanfov * screenY;
			genRay = { cameraX , cameraY , -1 };
			dir = glm::normalize(genRay);
			rayOrigin = { 0.0, 0.0, 0.0 };
			Ray ray2(rayOrigin, dir);
			rays[1] = ray2;

			//Ray 3
			screenX = (2 * ((x + 0.25) / WIDTH)) - 1;
			screenY = (2 * ((y + 0.50) / HEIGHT)) - 1;
			cameraX = aspectRatio * tanfov * screenX;
			cameraY = tanfov * screenY;
			genRay = { cameraX , cameraY , -1 };
			dir = glm::normalize(genRay);
			rayOrigin = { 0.0, 0.0, 0.0 };
			Ray ray3(rayOrigin, dir);
			rays[2] = ray3;

			//Ray 4 
			screenX = (2 * ((x + 0.50) / WIDTH)) - 1;
			screenY = (2 * ((y + 0.50) / HEIGHT)) - 1;
			cameraX = aspectRatio * tanfov * screenX;
			cameraY = tanfov * screenY;
			genRay = { cameraX , cameraY , -1 };
			dir = glm::normalize(genRay);
			rayOrigin = { 0.0, 0.0, 0.0 };
			Ray ray4(rayOrigin, dir);
			rays[3] = ray4;

			//Ray 5 
			screenX = (2 * ((x + 0.75) / WIDTH)) - 1;
			screenY = (2 * ((y + 0.50) / HEIGHT)) - 1;
			cameraX = aspectRatio * tanfov * screenX;
			cameraY = tanfov * screenY;
			genRay = { cameraX , cameraY , -1 };
			dir = glm::normalize(genRay);
			rayOrigin = { 0.0, 0.0, 0.0 };
			Ray ray5(rayOrigin, dir);
			rays[4] = ray5;

			//Ray 6 
			screenX = (2 * ((x + 0.75) / WIDTH)) - 1;
			screenY = (2 * ((y + 0.75) / HEIGHT)) - 1;
			cameraX = aspectRatio * tanfov * screenX;
			cameraY = tanfov * screenY;
			genRay = { cameraX , cameraY , -1 };
			dir = glm::normalize(genRay);
			rayOrigin = { 0.0, 0.0, 0.0 };
			Ray ray6(rayOrigin, dir);
			rays[5] = ray6;

			return rays;
		}

		else
			return rays;
	}
}

Color PhongShading(glm::vec3 Lpos, glm::vec3 Ldir, double lColor[], char objectType, int objectIndex, glm::vec3 IntersectionPt)
{

#pragma region sphere phong shading
	if (objectType == 's')
	{
		Sphere sphere = spheres[objectIndex];
		glm::vec3 sphereCenter = glm::make_vec3(sphere.position);
		glm::vec3 n = glm::normalize(IntersectionPt - sphereCenter);

		double ldotn = glm::dot(Ldir, n);
		ldotn = clampValues(ldotn);

		glm::vec3 r = { 2 * ldotn * n.x - Ldir.x,
						2 * ldotn * n.y - Ldir.y,
						2 * ldotn * n.z - Ldir.z };

		r = glm::normalize(r);
		double rdotv = glm::dot(r, glm::normalize(-IntersectionPt));
		rdotv = clampValues(rdotv);

		double rSpecular = sphere.color_specular[0] * std::pow(rdotv, sphere.shininess);
		double gSpecular = sphere.color_specular[1] * std::pow(rdotv, sphere.shininess);
		double bSpecular = sphere.color_specular[2] * std::pow(rdotv, sphere.shininess);

		double rDiffuse = sphere.color_diffuse[0] * ldotn;
		double gDiffuse = sphere.color_diffuse[1] * ldotn;
		double bDiffuse = sphere.color_diffuse[2] * ldotn;

		double rI = lColor[0] * (rDiffuse + rSpecular);
		double gI = lColor[1] * (gDiffuse + gSpecular);
		double bI = lColor[2] * (bDiffuse + bSpecular);

		return Color(rI, gI, bI);
	}
#pragma endregion sphere phong shading

#pragma region triangle phong shading
	if (objectType == 't')
	{
		Triangle triangle = triangles[objectIndex];

		Vertex v0 = triangle.v[0];
		Vertex v1 = triangle.v[1];
		Vertex v2 = triangle.v[2];

		//glm::vec3 triangleA = { v0.position[0], v0.position[1], v0.position[2] };
		//glm::vec3 triangleB = { v1.position[0], v1.position[1], v1.position[2] };
		//glm::vec3 triangleC = { v2.position[0], v2.position[1], v2.position[2] };
		
		//glm::vec3 triangleA, triangleB, triangleC;
		//getVertex(triangle, &triangleA, &triangleB, &triangleC);

		//populate the normal (B-A) x (C-A)
		/*glm::vec3 diff1 = triangleB - triangleA;
		glm::vec3 diff2 = triangleC - triangleA;

		glm::vec3 n = glm::cross(diff1, diff2);
		*/
		//glm::vec3 normal_normalised = glm::normalize(n);

		double alpha, beta, gamma;
		computeBarrycentricCoords(triangle, IntersectionPt, &alpha, &beta);
		gamma = 1.0 - alpha - beta;

		glm::vec3 interpolatedN = { alpha * v0.normal[0] + beta * v1.normal[0] + gamma * v2.normal[0],
									alpha * v0.normal[1] + beta * v1.normal[1] + gamma * v2.normal[1],
									alpha * v0.normal[2] + beta * v1.normal[2] + gamma * v2.normal[2] };
		
		interpolatedN = glm::normalize(interpolatedN);

		double ldotn = glm::dot(Ldir, interpolatedN);
		ldotn = clampValues(ldotn);

		glm::vec3 r = { 2 * ldotn * interpolatedN.x - Ldir.x ,
						2 * ldotn * interpolatedN.y - Ldir.y ,
						2 * ldotn * interpolatedN.z - Ldir.z };

		r = glm::normalize(r);
		double rdotv = glm::dot(r, glm::normalize(-IntersectionPt));
		rdotv = clampValues(rdotv);

		double a = alpha * v0.shininess + beta * v1.shininess + gamma * v2.shininess;

		double rSpecular = (alpha * v0.color_specular[0] + beta * v1.color_specular[0] + gamma * v2.color_specular[0]) * std::pow(rdotv, a);
		double gSpecular = (alpha * v0.color_specular[1] + beta * v1.color_specular[1] + gamma * v2.color_specular[1]) * std::pow(rdotv, a);
		double bSpecular = (alpha * v0.color_specular[2] + beta * v1.color_specular[2] + gamma * v2.color_specular[2]) * std::pow(rdotv, a);

		double rDiffuse = (alpha * v0.color_diffuse[0] + beta * v1.color_diffuse[0] + gamma * v2.color_diffuse[0]) * ldotn;
		double gDiffuse = (alpha * v0.color_diffuse[1] + beta * v1.color_diffuse[1] + gamma * v2.color_diffuse[1]) * ldotn;
		double bDiffuse = (alpha * v0.color_diffuse[2] + beta * v1.color_diffuse[2] + gamma * v2.color_diffuse[2]) * ldotn;

		double rI = lColor[0] * (rDiffuse + rSpecular);
		double gI = lColor[1] * (gDiffuse + gSpecular);
		double bI = lColor[2] * (bDiffuse + bSpecular);

		//std::cout << rI << ' ' << gI << ' ' << bI << std::endl;

		return Color(rI, gI, bI);

	}

#pragma endregion triangle phong shading

}

bool shadowIntersectionTest(glm::vec3 Lpos, glm::vec3 Ldir, char testWith , glm::vec3 intersectionPt, const int objectIndex)
{
	bool flag_shadow = false;
	Ray shadowRay(intersectionPt, Ldir);

	if (testWith == 's')
	{
		//shadowRay.rayOrigin = intersectionPt;
		//shadowRay.rayDirection = Ldir;

		for (int i = 0; i < num_spheres; i++)
		{
			glm::vec3 shadowSphereIpt = { 0.0, 0.0, -1e-10 };
			if ((i != objectIndex) && shadowRay.sphereIntersectionTest(spheres[i], shadowSphereIpt))
			{
				//check the distance between intersection point and light source
				//(std::pow(v1.x, 2) + std::pow(v1.y, 2) + std::pow(v1.z, 2))
				double sphere_shadow = glm::length(shadowSphereIpt - intersectionPt);
				double light_sphere = glm::length(Lpos - intersectionPt);
				if (sphere_shadow < light_sphere)
				{
					flag_shadow = true;
					break;
				}
			}
		}

		for (int j = 0; j < num_triangles; j++)
		{
			glm::vec3 shadowTriangleIpt = { 0.0, 0.0, -1e-10 };
			if (shadowRay.triangleIntersectionTest(triangles[j], shadowTriangleIpt))
			{
				double triangle_shadow = glm::length(shadowTriangleIpt - intersectionPt);
				double light_triangle = glm::length(Lpos - intersectionPt);
				if (triangle_shadow < light_triangle)
				{
					flag_shadow = true;
					break;
				}
			}
		}

		return flag_shadow;
	}

	if (testWith == 't')
	{
		//shadowRay.rayOrigin = intersectionPt;
		//shadowRay.rayDirection = Ldir;

		for (int i = 0; i < num_spheres; i++)
		{
			glm::vec3 shadowSphereIpt = { 0.0, 0.0, EPSILON };
			if (shadowRay.sphereIntersectionTest(spheres[i], shadowSphereIpt))
			{
				double sphere_shadow = glm::length(shadowSphereIpt - intersectionPt);
				double light_sphere = glm::length(Lpos - intersectionPt);
				if (sphere_shadow < light_sphere)
				{
					flag_shadow = true;
					break;
				}
			}
		}

		for (int j = 0; j < num_triangles; j++)
		{
			glm::vec3 shadowTriangleIpt = { 0.0, 0.0, EPSILON };
			if ((j != objectIndex) && shadowRay.triangleIntersectionTest(triangles[j], shadowTriangleIpt))
			{
				double triangle_shadow = glm::length(shadowTriangleIpt - intersectionPt);
				double light_triangle = glm::length(Lpos - intersectionPt);
				if (triangle_shadow < light_triangle)
				{
					flag_shadow = true;
					break;
				}
			}
		}

		return flag_shadow;
	}

	return flag_shadow;
}

Color raytracer(Ray Primaryray)
{
	double nearestIntersection = -1e10;
	Color finColor = Color(1.0, 1.0, 1.0);

#pragma region ray-trace spheres 

	for (int i = 0; i < num_spheres; i++)
	{
		glm::vec3 iPt = {0.0, 0.0, -1e-10 };
		bool intersectsSphere = Primaryray.sphereIntersectionTest(spheres[i], iPt);
		bool isNearest = iPt.z > nearestIntersection ? true : false;
		if (intersectsSphere && isNearest)
		{
			for (int j = 0; j <num_lights; j++)
			{
				glm::vec3 currLightPos = {lights[j].position[0], lights[j].position[1], lights[j].position[2]};
				glm::vec3 lightDirection = glm::normalize(currLightPos - iPt);
				
				bool shadowTest = shadowIntersectionTest(currLightPos, lightDirection, 's', iPt, i);
				//std::cout << shadowTest << std::endl;
				if (shadowTest)
					finColor = Color(0.0, 0.0, 0.0);
				else
				{
					finColor = Color(0.0, 0.0, 0.0);
					Color phong = PhongShading(currLightPos, lightDirection, lights[j].color, 's', i, iPt);
					finColor.r = finColor.r + phong.r;
					finColor.g = finColor.g + phong.g;
					finColor.b = finColor.b + phong.b;
				}
			}

			nearestIntersection = iPt.z;
		}
	}
#pragma endregion ray-trace spheres 

#pragma region ray-trace triangles 

	//nearestIntersection = EPSILON;
	for (int i = 0; i < num_triangles; i++)
	{
		glm::vec3 iPt = {0.0, 0.0, -1e-10};
		bool intersectsTriangle = Primaryray.triangleIntersectionTest(triangles[i], iPt);
		bool isNearest = iPt.z > nearestIntersection ? true : false;
		if (intersectsTriangle && isNearest)
		{
			for (int j = 0; j < num_lights; j++)
			{
				glm::vec3 currLightPos = { lights[j].position[0], lights[j].position[1], lights[j].position[2] };
				glm::vec3 lightDirection = glm::normalize(currLightPos - iPt);

				bool shadowTest = shadowIntersectionTest(currLightPos, lightDirection, 't', iPt, i);
				//std::cout << shadowTest << std::endl;
				if (shadowTest)
					finColor = Color(0.0, 0.0, 0.0);
				else
				{
					finColor = Color(0.0, 0.0, 0.0);
					Color phong = PhongShading(currLightPos, lightDirection, lights[j].color, 't' , i , iPt);
					finColor.r = finColor.r + phong.r;
					finColor.g = finColor.g + phong.g;
					finColor.b = finColor.b + phong.b;
				}
			}
			nearestIntersection = iPt.z;
		}
	}
#pragma endregion ray-trace triangles 

	Color addAmbient = { ambient_light[0],ambient_light[1], ambient_light[2] };

	finColor.r = clampValues(finColor.r + addAmbient.r);
	finColor.g = clampValues(finColor.g + addAmbient.g);
	finColor.b = clampValues(finColor.b + addAmbient.b);

	//std::cout << finColor.r << ',' << finColor.g << ',' << finColor.b << std::endl;
	
	return finColor;
}

Color softshadowtracer(Ray Primaryray) {

	double nearestIntersection = -1e10;
	Color finColor(1.0, 1.0, 1.0);

	for (int i = 0; i < num_spheres; i++) {
		glm::vec3 iPt = { 0.0, 0.0, -1e-10 };
		bool intersectsSphere = Primaryray.sphereIntersectionTest(spheres[i], iPt);
		bool isNearest = iPt.z > nearestIntersection ? true : false;
		if (intersectsSphere && isNearest)
		{
			for (int j = 0; j < num_lights; j++)
			{

				Light sampledLight[20];
				for (int k = 0; k < 20; k++)
				{
					double c = (double)rand() / (RAND_MAX);
					double t = (double)(rand() % 360) * 0.1745;
					double p = (double)(rand() % 360) * 0.1745;

					double positions[3] = {
						lights[j].position[0] + (c / 10.0) * std::sin(t) * std::cos(p),
						lights[j].position[1] + (c / 10.0) * std::sin(t) * std::cos(p),
						lights[j].position[2] + (c / 10.0) * std::cos(t)
					};

					double col[3] = { lights[j].color[0] / 20.0, lights[j].color[1] / 20.0, lights[j].color[2] / 20.0 };

					Light sample;
					sample.position[0] = positions[0];
					sample.position[1] = positions[1];
					sample.position[2] = positions[2];
					sample.color[0] = col[0];
					sample.color[1] = col[1];
					sample.color[2] = col[2];

					sampledLight[k] = sample;
				}


				for (int l = 0; l < 20; l++)
				{
					glm::vec3 currLightPos = { sampledLight[l].position[0], sampledLight[l].position[1], sampledLight[l].position[2] };
					glm::vec3 lightDirection = glm::normalize(currLightPos - iPt);

					bool shadowTest = shadowIntersectionTest(currLightPos, lightDirection, 's', iPt, i);
					//std::cout << shadowTest << std::endl;
					if (shadowTest)
						finColor = Color(0.0, 0.0, 0.0);
					else
					{
						finColor = Color(0.0, 0.0, 0.0);
						Color phong = PhongShading(currLightPos, lightDirection, lights[j].color, 's', i, iPt);
						finColor.r = finColor.r + phong.r;
						finColor.g = finColor.g + phong.g;
						finColor.b = finColor.b + phong.b;
					}
				}
			}

			nearestIntersection = iPt.z;
		}
	}

	for (int i = 0; i < num_triangles; i++) {
		glm::vec3 iPt = { 0.0, 0.0, -1e-10 };
		bool intersectsTriangle = Primaryray.triangleIntersectionTest(triangles[i], iPt);
		bool isNearest = iPt.z > nearestIntersection ? true : false;
		if (intersectsTriangle && isNearest)
		{
			for (int j = 0; j < num_lights; j++)
			{
				Light sampledLight[20];
				for (int k = 0; k < 20; k++)
				{
					double c = (double)rand() / (RAND_MAX);
					double t = (double)(rand() % 360) * 0.1745;
					double p = (double)(rand() % 360) * 0.1745;

					double positions[3] = {
						lights[j].position[0] + (c / 10.0) * std::sin(t) * std::cos(p),
						lights[j].position[1] + (c / 10.0) * std::sin(t) * std::cos(p),
						lights[j].position[2] + (c / 10.0) * std::cos(t)
					};

					double col[3] = { lights[j].color[0] / 20.0, lights[j].color[1] / 20.0, lights[j].color[2] / 20.0 };

					Light sample;
					sample.position[0] = positions[0];
					sample.position[1] = positions[1];
					sample.position[2] = positions[2];
					sample.color[0] = col[0];
					sample.color[1] = col[1];
					sample.color[2] = col[2];

					sampledLight[k] = sample;
				}


				for (int l = 0; l < 20; l++)
				{
					glm::vec3 currLightPos = { sampledLight[l].position[0], sampledLight[l].position[1], sampledLight[l].position[2] };
					glm::vec3 lightDirection = glm::normalize(currLightPos - iPt);

					bool shadowTest = shadowIntersectionTest(currLightPos, lightDirection, 't', iPt, i);
					//std::cout << shadowTest << std::endl;
					if (shadowTest)
						finColor = Color(0.0, 0.0, 0.0);
					else
					{
						finColor = Color(0.0, 0.0, 0.0);
						Color phong = PhongShading(currLightPos, lightDirection, lights[j].color, 't', i, iPt);
						finColor.r = finColor.r + phong.r;
						finColor.g = finColor.g + phong.g;
						finColor.b = finColor.b + phong.b;
					}
				}
			}

			nearestIntersection = iPt.z;
		}
	}

	Color addAmbient = { ambient_light[0],ambient_light[1], ambient_light[2] };

	finColor.r = clampValues(finColor.r + addAmbient.r);
	finColor.g = clampValues(finColor.g + addAmbient.g);
	finColor.b = clampValues(finColor.b + addAmbient.b);

	return finColor;
}

//MODIFY THIS FUNCTION
void draw_scene()
{
  //a simple test output
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
		if (supersample) 
		{
			Color sampleColor;
			Ray* rays = generateRays(x, y, 6);
			for (int i = 0; i < 6; i++)
			{
				Ray primaryRay = rays[i];
				Color colorPixel;
				if (softshadows)
					colorPixel = softshadowtracer(primaryRay);
				else
					colorPixel	= raytracer(primaryRay);
				sampleColor.r = sampleColor.r + colorPixel.r;
				sampleColor.g = sampleColor.g + colorPixel.g;
				sampleColor.b = sampleColor.b + colorPixel.b;
			}
			sampleColor.r = sampleColor.r / 6.0;
			sampleColor.g = sampleColor.g / 6.0;
			sampleColor.b = sampleColor.b / 6.0;

			plot_pixel(x, y, sampleColor.r * 255, sampleColor.g * 255, sampleColor.b * 255);
		}
		else
		{
			Ray* rays = generateRays(x, y , 1);
			Ray primaryRay = rays[0];
			Color colorPixel;
			if (softshadows)
				colorPixel = softshadowtracer(primaryRay);
			else
				colorPixel= raytracer(primaryRay);
			plot_pixel(x, y, colorPixel.r * 255, colorPixel.g * 255, colorPixel.b * 255);
		}
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 5))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }

  if (argc == 5)
  {
	  mode = MODE_JPEG;
	  filename = argv[2];
	  char *antialias = argv[3];
	  if (*antialias == 'a' || *antialias == 'A')
	  {
		  supersample = true;
		  std::cout << "Antialising Enabled" << std::endl;
	  }
	  else
	  {
		  supersample = false;
		  std::cout << "Antialising Disabled" << std::endl;
	  }

	  char *s = argv[4];
	  if (*s == 's' || *s == 'S')
	  {
		  softshadows = true;
		  std::cout << "Soft shadows Enabled" << std::endl;
	  }
	  else
	  {
		  softshadows = false;
		  std::cout << "Soft shadows disabled" << std::endl;
	  }
  }

  if (argc == 4)
  {
	  mode = MODE_JPEG;
	  filename = argv[2];
	  char *antialias = argv[3];
	  if (*antialias == 'y' || *antialias == 'Y')
	  {
		  supersample = true;
		  std::cout << "Antialising Enabled" << std::endl;
	  }
	  else
	  {
		  supersample = false;
		  std::cout << "Antialising Disabled" << std::endl;
	  }
  }

  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
	std::cout << "Antialising/Soft Shadows disabled" << std::endl;
  }
  else if (argc == 2)
  {
	  mode = MODE_DISPLAY;
	  std::cout << "Antialising/Soft Shadows disabled" << std::endl;
  }

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

