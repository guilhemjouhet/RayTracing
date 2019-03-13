#define _CRT_SECURE_NO_WARNINGS 1
#include <random>
#include <Vector>
#include <algorithm>
#include <stdlib.h>
#include <omp.h>
#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#include "stb.h"
#include "stb_image_write.h"
#include "stb_image.h"
std::default_random_engine engine;
std::uniform_real_distribution <double> distrib(0, 1);


class Vector{
	//a class supporting 3D vectors
public :
	Vector(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z){};
	double x,y,z;

	//Calculating the norm 2
	double norm2() {
		return x * x + y * y + z * z;
	}

	//A method that normalises vectors
	void normalize(){
		double n = sqrt(norm2());
		x/=n;
		y/=n;
		z/=n;
	}

	//a transfomation that moves, rotates and scales a mesh
	void transforme() {
		double tempy,tempz;
		tempy = 15.*x;
		tempz = 15.*y ;
		x=-15.*z+10.;
		y = tempy;
		z = tempz;
	}
};

//VARIOUS OPERATIONS ON VECTORS
double dot(const Vector& a, const Vector& b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vector cross(const Vector& a, const Vector& b) {
	return Vector(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

Vector operator+(const Vector& a, const Vector& b){
	return Vector(a.x+b.x,a.y+b.y,a.z+b.z);
}
Vector operator-(const Vector& a, const Vector& b){
	return Vector(a.x-b.x,a.y-b.y,a.z-b.z);
}
Vector operator+(const Vector& a, const double b){
	return Vector(a.x+b,a.y+b,a.z+b);
}
Vector operator*(const Vector& a, const Vector& b) {
	return Vector(a.x * b.x, a.y * b.y, a.z * b.z);
}
Vector operator*(const double b,const Vector& a){
	return Vector(b*a.x,b*a.y,b*a.z);
}
Vector operator/(const Vector& a, const double b){
	return Vector(b/a.x,b/a.y,b/a.z);
}


class Ray {
	//a class that represents a ray of light.
	//C is the origin of the ray, U its direction
public:
	Ray(const Vector& C, const Vector& U) : C(C), U(U) {};
	Vector C;
	Vector U;
};

class Object {
	//A generic class to represent various types of object like spheres and meshes
public:
	//Object():{};
	virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t, Vector &albedo) = 0;
	Vector color;
	bool miroir;
	bool lumiere;
	bool transparent;

	void setColor(Vector &newcolor) {
		color = newcolor;
	}
};


class Sphere:public Object{
	//This class represents spheres
public:

	Sphere(const Vector& O, double R, const Vector& color, const bool miroir, const bool lumiere, const bool transparent = false) : O(O), R(R) {
		this->color=color;
		this->miroir = miroir;
		this->lumiere = lumiere;
		this->transparent = transparent;
	};
	Vector O;
	double R;

//detection of intersection between ray and spheres
	virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t, Vector &albedo) {
		double b = 2 * dot(r.U, r.C-O);
		double a = 1;
		double c = (r.C - O).norm2() - R * R;
		double delta = b * b - 4 * a*c;
		if (delta >= 0){
			double rac = sqrt(delta);
			double tp = (-b + rac) / (2*a);
			double tn = (-b - rac) / (2*a);
			if (tn > 0) {
				P = r.C + tn * r.U; //plus petite racine positive
				t = tn;
				albedo = color;
			}
			else {
				if (tp > 0) {
					P = r.C + tp * r.U; //sinon la seule positive
					t = tp;
					albedo = color;
				}
				else {
					return false;
				}
			}
			N = P - O; // Normale = rayon de la sphère
			N.normalize();
			return true;
		}
		else {
			return false;
		}
	}
};

class Triangle{
	//a triangle have 3 vertexes
public:
	Triangle(const Vector& A, const Vector& B, const Vector& C) : A(A),B(B),C(C){};
	Vector A;
	Vector B;
	Vector C;
	//intersection between ray and triangle
	virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t, double &alpha, double &beta, double &gamma) {
		//intersection avec le plan du triangle
		Vector Normale;
		Normale = cross(B - A, A - C);
		if (dot(r.U, Normale) > 0) Normale = Vector(0,0,0)-Normale;
		Normale.normalize();
		t = dot(A - r.C, Normale) / dot(r.U, Normale);
		if (t <= 0) {
			return false;
		}
		else {
			Vector Pinter = r.C + t * r.U;
			//on teste si l'intersection est dans le triangle ABC
			beta = (dot(Pinter - A, B - A)*(C - A).norm2() - dot(C - A, B - A)*dot(Pinter - A, C - A)) / ((B - A).norm2()*(C - A).norm2() - dot(B - A, C - A)*dot(B - A, C - A));
			gamma = ((B - A).norm2()*dot(Pinter - A, C - A) - dot(B - A, C - A)*dot(Pinter - A, B - A)) / ((B - A).norm2()*(C - A).norm2() - dot(B - A, C - A)* dot(B - A, C - A));
			alpha = 1 -gamma - beta;
			if ((beta > 0) && (gamma > 0) && (alpha > 0) && (beta < 1) && (gamma < 1) && (alpha < 1)) {
				P = Pinter;
				N = Normale;
				return true;
			}
			else { return false; }
		}
	}
};

class TriangleIndices {
	//class provided by Geometry class
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class Geometry : public Object {
	//A geometry represents a mesh
public:
	~Geometry() {}
	Geometry(const Vector& color, const bool miroir, const bool lumiere, const bool transparent = false) {
		this->color = color;
		this->miroir = miroir;
		this->lumiere = lumiere;
		this->transparent = transparent;
		//loading mesh
		this->readOBJ("C:\\Users\\guilh_000\\Documents\\Etudes\\ecl-3A\\INFO\\MOS Infographie\\Projet\\Beautiful Girl.obj");
		//initialisation
		for (int i = 0; i < vertices.size(); i++) {
			vertices[i].transforme();
		}
		for (int ii = 0; ii < normals.size(); ii++) {
			normals[ii].transforme();
		}
		double x, y, z;
		Xmin = vertices[0].x; Ymin = vertices[0].y; Zmin = vertices[0].z; Xmax = vertices[0].x; Ymax = vertices[0].y; Zmax = vertices[0].z;
		//déclaration des normales
		normalesBoite.push_back(Vector(-1, 0, 0)); normalesBoite.push_back(Vector(1, 0, 0)); normalesBoite.push_back(Vector(0,-1, 0)); normalesBoite.push_back(Vector(0,1, 0)); normalesBoite.push_back(Vector(0,0,-1)); normalesBoite.push_back(Vector(0,0,1));
		//on itère sur les sommets
		for (int i = 1; i < vertices.size(); i++) {
			x = vertices[i].x; y = vertices[i].y; z = vertices[i].z;
			if (x < Xmin)  Xmin = x;
			else if (x > Xmax) Xmax = x;
			if (y < Ymin)  Ymin = y;
			else if (y > Ymax) Ymax = y;
			if (z < Zmin)  Zmin = z;
			else if (z> Zmax) Zmax = z;
		}
		std::cout << "xmin " << Xmin << "\n"; std::cout << "xmax " << Xmax << "\n";
		std::cout << "ymin " << Ymin << "\n"; std::cout << "ymax " << Ymax << "\n";
		std::cout << "zmin " << Zmin << "\n"; std::cout << "zmax " << Zmax << "\n";
		//loading textures
		add_texture("visage.bmp");
		add_texture("cheveux.bmp");
		add_texture("corps.bmp");
		add_texture("pantalon.bmp");
		add_texture("accessoires.bmp");
		add_texture("mains.bmp");
	};

	void add_texture(const char* filename) {
		int w=0, h=0, c=0;
		textures.push_back(stbi_load(filename, &w, &h, &c,3));
		textures_W.push_back(w);
		textures_H.push_back(h);
	}


	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec.x, &vec.y, &vec.z, &col.x, &col.y, &col.z) == 6) {
					col.x = std::min(1., std::max(0., col.x));
					col.y = std::min(1., std::max(0., col.y));
					col.z = std::min(1., std::max(0., col.z));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				}
				else {
					sscanf(line, "v %lf %lf %lf\n", &vec.x, &vec.y, &vec.z);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec.x, &vec.y, &vec.z);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec.x, &vec.y);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				}
				else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					}
					else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						}
						else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					}
					else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						}
						else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							}
							else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								}
								else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);

	}

	double min(double a, double b, double c) {
		if (a <= b&&a <= c) return a;
		else if (b < a&&b <= c) return b;
		else return c;
	}
	double max(double a, double b, double c) {
		if (a >= b && a >= c) return a;
		else if (b> a&&b >= c) return b;
		else return c;
	}

	bool intersectbox(const Ray& r) {
		double txm, txM, tym, tyM, tzm, tzM;
		//intersections avec les faces de la boite englobante
		txm = dot(Vector(Xmin,0,0) - r.C, normalesBoite[0]) / dot(r.U, normalesBoite[0]);
		txM = dot(Vector(Xmax, 0, 0) - r.C, normalesBoite[1]) / dot(r.U, normalesBoite[1]);
		tym = dot(Vector(0, Ymin, 0) - r.C, normalesBoite[2]) / dot(r.U, normalesBoite[2]);
		tyM = dot(Vector(0,Ymax, 0) - r.C, normalesBoite[3]) / dot(r.U, normalesBoite[3]);
		tzm = dot(Vector(0,0,Zmin) - r.C, normalesBoite[4]) / dot(r.U, normalesBoite[4]);
		tzM = dot(Vector(0,0,Zmax) - r.C, normalesBoite[5]) / dot(r.U, normalesBoite[5]);

		double txmin = min(txm, txM, txM); double txMax = max(txm, txM, txM);
		double tymin = min(tym, tyM, tyM); double tyMax = max(tym, tyM, tyM);
		double tzmin = min(tzm, tzM, tzM); double tzMax = max(tzm, tzM, tzM);
		//Condition d’intersection avec la boite
		if (max(txmin, tymin, tzmin) < min(txMax, tyMax, tzMax)) {
			return true;
		}
		else return false;
	}

	//éléments du maillage
	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;

	std::vector<unsigned char*> textures;
	std::vector<int> textures_W;
	std::vector<int> textures_H;

	//boite englobante
	double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
	std::vector<Vector> normalesBoite;

	virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t, Vector &albedo) {
		double smallestt = 1E15;
		if (intersectbox(r)) {
			for (int indice = 0; indice < indices.size(); indice++) {
				Vector Plocal, Nlocal;
				Vector UV;
				double tloc = 0;
				//création des triangles
				Triangle triangle(vertices[indices[indice].vtxi], vertices[indices[indice].vtxj], vertices[indices[indice].vtxk]);
				double alpha, beta, gamma;
				//calcul des intersections
				bool inter = triangle.intersect(r, Plocal, Nlocal, tloc,alpha,beta,gamma);
				if (inter && tloc < smallestt) {
					smallestt = tloc;
					P = Plocal;
					Vector norma = normals[indices[indice].ni];
					Vector normb = normals[indices[indice].nj];
					Vector normc = normals[indices[indice].nk];
					//interpolation de la normale
					if (dot(norma, r.U) > 0) {
						norma = Vector(0, 0, 0) - norma;
					}
					if (dot(normb, r.U) > 0) {
						normb = Vector(0, 0, 0) - normb;
					}
					if (dot(normc, r.U) > 0) {
						normc = Vector(0, 0, 0) - normc;
					}
					N = alpha * norma + beta * normb + gamma * normc;
					N.normalize();
					if (dot(N, r.U) > 0) {
						N = Vector(0, 0, 0) - N;
					}
					/*if (indices[indice].group == 4) {
						albedo = Vector(255, 0, 0);
					}
					else if (indices[indice].group == 5) {
						albedo = Vector(0,255, 0);
					}
					else {
						albedo = Vector(0, 0, 0);
					}*/
					UV = alpha * uvs[indices[indice].uvi] + beta * uvs[indices[indice].uvj] + gamma * uvs[indices[indice].uvk];
					int x = fabs(fmod(UV.x, 1.))*(textures_W[indices[indice].group]-1);
					int y = fabs(fmod(UV.y, 1.))*(textures_H[indices[indice].group]-1);
					int adresse = 3*(y * textures_W[indices[indice].group] + x);
					albedo = Vector(textures[indices[indice].group][adresse], textures[indices[indice].group][adresse + 1], textures[indices[indice].group][adresse + 2]);
				}
			}
			if (smallestt < 1E14) {
				t = smallestt;
				return true;
			}
			else return false;
		}
		else return false;
	}
};

class Scene {
	//a scene have various object that will be represented in the image
public:
	Scene() {};
	std::vector<Object*> objects;

	void addObject(Object* s) {
		objects.push_back(s);
	}
//intersection between ray and all the scene's object
	bool intersect(const Ray& r, Vector &P, Vector &N, int &ii, double &tprime,Vector &albedo) {
		double smallestt = 1E15;
		//iteration over all the objects
		for (int i = 0; i < objects.size(); i++) {
			Vector Plocal, Nlocal, albedoloc;
			double t = 0;
			//calling the good intersection routine
			bool inter = objects[i]->intersect(r, Plocal, Nlocal,t,albedoloc);
			if (inter && t < smallestt) {
				smallestt = t;
				P = Plocal;
				N = Nlocal;
				ii = i;
				albedo = albedoloc;
			}
		}
		if (smallestt < 1E14) {
			tprime = smallestt;
			return true;
		}
		else {
			return false;
		}
	}

 
 //getcolor returns the color of a ray of light. Various rays are computed for each pixel
	//This method calls itself in order to get indirect light
	Vector getcolor(const Ray& ray, const Vector& L, int nrebond, int nChemin=1) {
		//here is the condition of stop
		if (nrebond == 0 || nChemin==0) { 
			return Vector(0, 0, 0); 
		}
		double I = 20000000;
		Vector P, N;//création des vecteurs d'intesection : point d'intersection et normale à la sphère
		int ii = -1;
		double t;
		Vector albedo;
		//searching for the first object intecepted by the ray
		if (intersect(ray, P, N, ii, t,albedo)) { //on trouve une sphère
			Object* SS = objects[ii];
			//miror case
			if (SS->miroir) {
				Vector R = ray.U - 2 * dot(ray.U, N)*N; //rayon réfléchit en P vers l'infini
				R.normalize();
				Ray rayon_lumiere(P + 0.0001*N, R); //idem
				Vector color = getcolor(rayon_lumiere, L, nrebond-1,nChemin);//appel récursif
				return color;
			}
			//transparent case
			else if (SS->transparent) {
				float nairnverre = 0.5;
				//tuning the direction of the normal of the sphere
				if (dot(ray.U, N) > 0) {
					N = operator-(Vector(0,0,0),N);
					nairnverre = 1 / nairnverre;
				}
				//Snell descartes law
				float det = 1 - nairnverre * nairnverre*(1 - dot(ray.U, N)*dot(ray.U, N));
				if (det < 0) {
					//reflected ray
					Vector R = ray.U - 2 * dot(ray.U, N)*N; //rayon réfléchit en P vers l'infini
					R.normalize();
					Ray rayon_lumiere(P + 0.0001*N, R); //idem
					Vector color = getcolor(rayon_lumiere, L, nrebond-1,nChemin);
					return color;

				}
				else {
					//transmitted ray
					Vector R = nairnverre * ray.U - (nairnverre*dot(ray.U, N) + sqrt(det))*N;
					R.normalize();
					Ray rayon_lumiere(P - 0.0001*N, R); //idem
					Vector color = getcolor(rayon_lumiere, L, nrebond-1,nChemin);
					return color;
				}
			}
			//Source of light
			else if (SS->lumiere) {
				Vector color = I/4./3.14/3.14/dynamic_cast<Sphere*>(SS)->R/ dynamic_cast<Sphere*>(SS)->R*albedo;
				return color;
			}
			
			else {//contribution directe
				Vector color(0, 0, 0);

				Vector PL = L - P;
				double distance = PL.norm2();
				PL.normalize();
				Ray rayon_lumiere(P + 0.00001*N, PL);
				
				Vector Pprime, Nprime;
				int objetprime;
				double tprime;
				Vector newalb;
				bool ombre = intersect(rayon_lumiere, Pprime, Nprime, objetprime, tprime,newalb);
				Vector PPprime = P - Pprime;
				double distanceprime = PPprime.norm2();
				if (ombre && distanceprime < distance && !objects[objetprime]->transparent && !objects[objetprime]->lumiere) {
					//on laisse color à 0 si on est dans l'ombre d'un autre object 
				}
				else if (objects[objetprime]->lumiere){
					// direct cas ponctuel
					color = (I / 3.14*std::max(0., dot(N, PL)) / (4 * 3.14*distance))*albedo;

					//direct cas étendu
					/*
					Vector OX = P - L;
					Vector OXprime = Pprime - L;
					color = I / 4. / 3.14 / 255.*abs(dot(ray.U, N)*dot(PPprime, Nprime)/ dot(OX, OXprime)) / distanceprime *spheres[0].color;
					int a = 1;*/
				}
				if (nChemin > 1) {//éclairage indirect for the first ray of the three : we split it in various rays
					for (int j = 0; j < nChemin; j++) {
						double r1 = distrib(engine);
						double r2 = distrib(engine);
						Vector R_aleatoire_local(cos(2 * 3.14*r1)*sqrt(1 - r2), sin(2 * 3.14*r1)*sqrt(1 - r2), sqrt(r2));
						Vector aleatoire(distrib(engine)-0.5, distrib(engine)-0.5, distrib(engine)-0.5);
						Vector tangent1 = cross(N, aleatoire);
						tangent1.normalize();
						Vector tangent2 = cross(N, tangent1);

						Vector R_aleatoire = R_aleatoire_local.z*N + R_aleatoire_local.x*tangent1 + R_aleatoire_local.y*tangent2;
						R_aleatoire.normalize();
						Ray rayon_aleatoire(P + 0.0001*N, R_aleatoire); //idem
						color = color + 1/ (double)nChemin/255.*getcolor(rayon_aleatoire, L, nrebond-1,1)*albedo;// on met nrebond=1 pour pas diviser le rayon lors des futurs rebonds
					}
				}
				else {
					Vector R = ray.U - 2 * dot(ray.U, N)*N; //rayon réfléchit en P vers l'infini
					R.normalize();
					Ray rayon_lumiere(P + 0.0001*N, R); //idem
					//on diminue le nombre de rebonds restants
					color = color +  1./255.*getcolor(rayon_lumiere, L, nrebond - 1, 1)*albedo;
				}
				return color;
			}
		}
		else {
			Vector color(0, 0, 0);
			return color;
		}
	}
};

int main() {
	//création de la scène et de ses objets
	Scene scene;
	Vector C(0, 0, 55);
	Vector L(-20,-20, 30);
	Vector colorL(255, 255, 255);
	Sphere Lum(L, 10, colorL, false, true);//lumiere
	scene.addObject(&Lum);
	int W =200;
	int H = 200;
	int R = 990;
	double fov = 35* 3.1415 / 180;
	double tanhalfov = -W / 2 / tan(fov / 2);
	std::vector<unsigned char> image(W*H * 3, 0);
	Vector O(0, -15, 0); // centre de la sphère
	Vector color(255, 0,0);
	Sphere S(O, 10, color, true,false); //miroir
	scene.addObject(&S);
	Vector OO(0, 15,0); // centre de la sphère transparente
	Sphere St(OO, 10, color, false,false,true); //Sphère transparente
	scene.addObject(&St);
	Vector colorb(180, 180, 180);
	Vector Ob(1000, 0, 0); // centre de la sphère du bas
	Sphere Sb(Ob, R, colorb, false,false);
	scene.addObject(&Sb);
	Vector colorh(0, 100,100);
	Vector Oh(-1200, 0, 0); // centre de la sphère du haut
	Sphere Sh(Oh, R, colorh, false, false);
	scene.addObject(&Sh);
	Vector colorder(255, 25,255);
	Vector Oder(0, 0, 1000); // centre de la sphère de derrière
	Sphere Sder(Oder, 940, colorder, false,false);
	scene.addObject(&Sder);
	Vector colorg(20, 150,30);
	Vector Og(0, 1000, 0); // centre de la sphère de gauche
	Sphere Sg(Og, 970, colorg, false,false);
	scene.addObject(&Sg);
	Vector colord(100, 10, 10);
	Vector Od(0, -1000, 0); // centre de la sphère de gauche
	Sphere Sd(Od, 970, colord, false, false);
	scene.addObject(&Sd);
	Vector colorf(0, 0,255);
	Vector Of(0, 0, -1200); // centre de la sphère du fond
	Sphere Sf(Of, 980, colorf, false,false);
	scene.addObject(&Sf);

	Vector AA(0, 15, 0);
	Vector BB(0,-15,0);
	Vector CC(-10, 0, -10);
	Vector colortri(0,255, 0);
	Geometry girl(colortri, false, false);
	
	scene.addObject(&girl);

	
	int nChemin =  5;//50
	int nrebond =  6;//10
	int nrays = 4;
	double distnette = 55;
	//Calcul de l'image
#pragma omp parallel for
	for (int i = 0; i < H; i++) { //on itère sur les lignes
		if (i % 10 == 0) std::cout << i / (double)H * 100 << "%\n";//affichage de l'avancement
		for (int j = 0; j < W; j++) {
			{//on itère sur les colones
				Vector color(0., 0., 0.);// on initialise la couleur du pixel au noir
				for (int k = 0; k < nrays; k++) {// On envoie des rayons sur ce pixel avec un petit angle aléatoir de déviation : Anti-aliasing
					//Angle aléatoire pour l'Anti-aliasing
					double r1 = distrib(engine);
					double r2 = distrib(engine);
					double R = 0.25*sqrt(-2 * log(r1));
					double dx = R * cos(2 * 3.14*r2);
					double dy = R * sin(2 * 3.14*r2);
					//Angle aléatoire pour le flou de profondeur
					double rc1 = distrib(engine);
					double rc2 = distrib(engine);
					double Rc = sqrt(-2 * log(rc1));
					double dcx = Rc * cos(2 * 3.14*rc2);
					double dcy = Rc * sin(2 * 3.14*rc2);
					//Création du vecteur directeur du rayon lumineux
					Vector U(i + 0.5 - H / 2 + dx, -j - 0.5 + dy + W / 2, tanhalfov); //calcul du vecteur directeur du rayon
					U.normalize();
					//Variation de ce vecteur pour le flou de profondeur
					Vector UU = distnette * U;
					Vector DC(dcx, dcy, 0);
					Vector origin = C + DC;
					Vector direction = UU - DC;
					direction.normalize();
					//création du rayon lumineux
					Ray ray(origin, direction);
					//calcul de la couleur
					Vector coloradd = scene.getcolor(ray, L, nrebond, nChemin);
					//ajout de la couleur calculé pour ce rayon
					color = color + 1 / (double)nrays*coloradd;
				}
				//Création de l'image
				image[(i*W + j) * 3 + 0] = std::min(255., std::pow(color.x, 0.45));
				image[(i*W + j) * 3 + 1] = std::min(255., std::pow(color.y, 0.45));
				image[(i*W + j) * 3 + 2] = std::min(255., std::pow(color.z, 0.45));
			}
		}
		stbi_write_png("image.png", W, H, 3, &image[0], 0);
		return 0;
	}
};