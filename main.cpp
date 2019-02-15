#define _CRT_SECURE_NO_WARNINGS 1
#include <random>
#include <Vector>
#include <algorithm>
#include <stdlib.h>
#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
//#define double PI = 3.1415
#include "stb_image_write.h"
std::default_random_engine engine;
std::uniform_real_distribution <double> distrib(0, 1);

class Vector{
public :
	Vector(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}; //this.x=x,...
	double x,y,z;
	double norm2() {
		return x * x + y * y + z * z;
	}
	void normalize(){
		double n = sqrt(norm2());
		x/=n;
		y/=n;
		z/=n;
	}
	void transforme() {
		double temp;
		temp = 15.*x;
		y = 15.*y ;
		x=-15.*z+10.;
		z = temp;
	}
};


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
public:
	Ray(const Vector& C, const Vector& U) : C(C), U(U) {};
	Vector C;
	Vector U;
};

class Object {
public:
	//Object():{};
	virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t) = 0;
	Vector color;
	bool miroir;
	bool lumiere;
	bool transparent;
};


class Sphere:public Object{
public:

	Sphere(const Vector& O, double R, const Vector& color, const bool miroir, const bool lumiere, const bool transparent = false) : O(O), R(R) {
		this->color=color;
		this->miroir = miroir;
		this->lumiere = lumiere;
		this->transparent = transparent;
	};
	Vector O;
	double R;

	virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t) {
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
			}
			else {
				if (tp > 0) {
					P = r.C + tp * r.U; //sinon la seule positive
					t = tp;
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

class Triangle :public Object{
public:
	Triangle(const Vector& A, const Vector& B, const Vector& C, const Vector& color, const bool miroir, const bool lumiere, const bool transparent = false) : A(A),B(B),C(C){
		this->color = color;
		this->miroir = miroir;
		this->lumiere = lumiere;
		this->transparent = transparent;
	};
	Vector A;
	Vector B;
	Vector C;

	virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t) {
		//intersection avec plan
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
			double beta = (dot(Pinter - A, B - A)*(C - A).norm2() - dot(C - A, B - A)*dot(Pinter - A, C - A)) / ((B - A).norm2()*(C - A).norm2() - dot(B - A, C - A)*dot(B - A, C - A));
			double gamma = ((B - A).norm2()*dot(Pinter - A, C - A) - dot(B - A, C - A)*dot(Pinter - A, B - A)) / ((B - A).norm2()*(C - A).norm2() - dot(B - A, C - A)* dot(B - A, C - A));
			double alpha = 1 - beta - gamma;
			if ((beta > 0) && (gamma > 0) && (alpha > 0) && (beta < 1) && (gamma < 1) && (alpha < 1)) {
				P = Pinter;
				N = Normale;
				return true;
			}
			else { return false; }
		}
	};
};

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class Geometry : public Object {
public:
	~Geometry() {}
	Geometry(const Vector& color, const bool miroir, const bool lumiere, const bool transparent = false) {
		this->color = color;
		this->miroir = miroir;
		this->lumiere = lumiere;
		this->transparent = transparent;
		this->readOBJ("C:\\Users\\guilh_000\\Documents\\Etudes\\ecl-3A\\INFO\\MOS Infographie\\Projet\\Beautiful Girl.obj");
		//détermination de la boite englobante
		//initialisation
		for (int i = 0; i < vertices.size(); i++) {
			vertices[i].transforme();
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
	};

	//0-10,0
	//*15
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
		txm = dot(Vector(Xmin,0,0) - r.C, normalesBoite[0]) / dot(r.U, normalesBoite[0]);
		txM = dot(Vector(Xmax, 0, 0) - r.C, normalesBoite[1]) / dot(r.U, normalesBoite[1]);
		tym = dot(Vector(0, Ymin, 0) - r.C, normalesBoite[2]) / dot(r.U, normalesBoite[2]);
		tyM = dot(Vector(0,Ymax, 0) - r.C, normalesBoite[3]) / dot(r.U, normalesBoite[3]);
		tzm = dot(Vector(0,0,Zmin) - r.C, normalesBoite[4]) / dot(r.U, normalesBoite[4]);
		tzM = dot(Vector(0,0,Zmax) - r.C, normalesBoite[5]) / dot(r.U, normalesBoite[5]);
		/*
		Vector Pxm = r.C + txm * r.U;
		Vector PxM = r.C + txM * r.U;
		Vector Pym = r.C + tym * r.U;
		Vector PyM = r.C + tyM * r.U;
		Vector Pzm = r.C + tzm * r.U;
		Vector PzM = r.C + tzM * r.U;
		*/
		double txmin = min(txm, txM, txM); double txMax = max(txm, txM, txM);
		double tymin = min(tym, tyM, tyM); double tyMax = max(tym, tyM, tyM);
		double tzmin = min(tzm, tzM, tzM); double tzMax = max(tzm, tzM, tzM);
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

	//boite englobante
	double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
	std::vector<Vector> normalesBoite;

	virtual bool intersect(const Ray& r, Vector &P, Vector &N, double &t) {
		double smallestt = 1E15;
		if (intersectbox(r)) {
			for (int indice = 0; indice < indices.size(); indice++) {
				Vector Plocal, Nlocal;
				double tloc = 0;
				Triangle triangle(vertices[indices[indice].vtxi], vertices[indices[indice].vtxj], vertices[indices[indice].vtxk], color, miroir, lumiere, transparent);
				bool inter = triangle.intersect(r, Plocal, Nlocal, tloc);
				if (inter && tloc < smallestt) {
					smallestt = tloc;
					P = Plocal;
					N = Nlocal;
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
public:
	Scene() {};
	std::vector<Object*> objects;

	void addObject(Object* s) {
		objects.push_back(s);
	}

	bool intersect(const Ray& r, Vector &P, Vector &N, int &ii, double &tprime) {
		double smallestt = 1E15;
		for (int i = 0; i < objects.size(); i++) {
			Vector Plocal, Nlocal;
			double t = 0;
			bool inter = objects[i]->intersect(r, Plocal, Nlocal,t);
			if (inter && t < smallestt) {
				smallestt = t;
				P = Plocal;
				N = Nlocal;
				ii = i;
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

	Vector getcolor(const Ray& ray, const Vector& L, int nrebond, int nChemin=1) {
		if (nrebond == 0 || nChemin==0) { 
			return Vector(0, 0, 0); 
		}
		double I = 20000000;
		Vector P, N;//création des vecteurs d'intesection : point d'intersection et normale à la sphère
		int ii = -1;
		double t;
		if (intersect(ray, P, N, ii, t)) { //on trouve une sphère
			Object* SS = objects[ii];
			if (SS->miroir) {
				Vector R = ray.U - 2 * dot(ray.U, N)*N; //rayon réfléchit en P vers l'infini
				R.normalize();
				Ray rayon_lumiere(P + 0.0001*N, R); //idem
				Vector color = getcolor(rayon_lumiere, L, nrebond-1,nChemin);
				return color;
			}
			else if (SS->transparent) {
				float nairnverre = 0.5;
				if (dot(ray.U, N) > 0) {
					N = operator-(Vector(0,0,0),N);
					nairnverre = 1 / nairnverre;
				}
				float det = 1 - nairnverre * nairnverre*(1 - dot(ray.U, N)*dot(ray.U, N));
				if (det < 0) {
					Vector R = ray.U - 2 * dot(ray.U, N)*N; //rayon réfléchit en P vers l'infini
					R.normalize();
					Ray rayon_lumiere(P + 0.0001*N, R); //idem
					Vector color = getcolor(rayon_lumiere, L, nrebond-1,nChemin);
					return color;

				}
				else {
					Vector R = nairnverre * ray.U - (nairnverre*dot(ray.U, N) + sqrt(det))*N;
					R.normalize();
					Ray rayon_lumiere(P - 0.0001*N, R); //idem
					Vector color = getcolor(rayon_lumiere, L, nrebond-1,nChemin);
					return color;
				}
			}
			else if (SS->lumiere) {
				Vector color = I/4./3.14/3.14/dynamic_cast<Sphere*>(SS)->R/ dynamic_cast<Sphere*>(SS)->R*SS->color;
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
				bool ombre = intersect(rayon_lumiere, Pprime, Nprime, objetprime, tprime);
				Vector PPprime = P - Pprime;
				double distanceprime = PPprime.norm2();
				if (ombre && distanceprime < distance && !objects[objetprime]->transparent && !objects[objetprime]->lumiere) {
					//on laisse color à 0
				}
				else if (objects[objetprime]->lumiere){
					// direct cas ponctuel
					color = (I / 3.14*std::max(0., dot(N, PL)) / (4 * 3.14*distance))*SS->color;

					//direct cas étendu
					/*
					Vector OX = P - L;
					Vector OXprime = Pprime - L;
					color = I / 4. / 3.14 / 255.*abs(dot(ray.U, N)*dot(PPprime, Nprime)/ dot(OX, OXprime)) / distanceprime *spheres[0].color;
					int a = 1;*/
				}

				//éclairage indirect

				if (nChemin > 1) {
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
						color = color + 1/ (double)nChemin/255.*getcolor(rayon_aleatoire, L, nrebond-1,1)*SS->color;// on met nrebond=1 pour pas passer en exponentiel
					}
				}
				else {
					Vector R = ray.U - 2 * dot(ray.U, N)*N; //rayon réfléchit en P vers l'infini
					R.normalize();
					Ray rayon_lumiere(P + 0.0001*N, R); //idem
					color = color +  1./255.*getcolor(rayon_lumiere, L, nrebond - 1, 1)*SS->color;
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
	Scene scene;
	Vector C(0, 0, 55);
	Vector L(-20,-20, 30);
	Vector colorL(255, 255, 255);
	Sphere Lum(L, 10, colorL, false, true);//lumiere
	scene.addObject(&Lum);
	//double I = 10000; // Intensité de la source
	int W = 1080;
	int H = 720;
	int R = 990;
	double fov = 65* 3.1415 / 180;
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
	Vector colortri(255, 0, 0);
	Triangle triangle1(AA, BB, CC, colortri,false,false);
	//scene.addObject(&triangle1);
	Geometry girl(colortri, false, false);
	
	scene.addObject(&girl);

	//Sphere SS(O, 10, color, false,false);
	int nChemin =  10;//50
	int nrebond =  7;//10
	int nrays = 4;
	double distnette = 55;
	#pragma omp parallel for schedule(dynamic,12)
	for (int i = 0; i < H; i++) { //pt pixel
		if (i%10==0) std::cout << i/(double)H*100 <<"%\n";
		for (int j = 0; j < W; j++) {
			Vector color(0., 0., 0.);
			for (int k = 0; k < nrays; k++) {
				double r1 = distrib(engine);
				double r2 = distrib(engine);
				double R = 0.25*sqrt(-2 * log(r1));
				double dx = R*cos(2 * 3.14*r2);
				double dy = R*sin(2 * 3.14*r2);
				double rc1 = distrib(engine);
				double rc2 = distrib(engine);
				double Rc = sqrt(-2 * log(rc1));
				double dcx = Rc * cos(2 * 3.14*rc2);
				double dcy = Rc * sin(2 * 3.14*rc2);

				dcx = 0;
				dcy = 0;

				Vector U(i + 0.5 - H / 2 + dx, -j - 0.5 + dy + W / 2, tanhalfov); //calcul du vecteur directeur du rayon
				U.normalize();
				Vector UU = distnette * U;
				Vector DC (dcx, dcy, 0);
				Vector origin = C + DC;
				Vector direction = UU-DC;
				direction.normalize();
				Ray ray(origin, direction);
				Vector coloradd = scene.getcolor( ray, L, nrebond, nChemin);
				color = color + 1 / (double)nrays*coloradd;
			}
			image[(i*W + j) * 3 + 0] = std::min(255., std::pow(color.x, 0.45));
			image[(i*W + j) * 3 + 1] = std::min(255., std::pow(color.y, 0.45));
			image[(i*W + j) * 3 + 2] = std::min(255., std::pow(color.z, 0.45));
			/*ICI*/
			/*
			Vector P, N;//création des vecteurs d'intesection : point d'intersection et normale à la sphère
			int ii = -1;
			double t;
			if (scene.intersect(ray, P, N, ii, t)) { //on trouve un objet
				Sphere SS = scene.spheres[ii];
				if (SS.miroir){ // c'est un miroir
					Vector PC = C - P;
					double distance1 = PC.norm2();
					Vector R = ray.U - 2 * dot(ray.U, N)*N; //rayon réfléchit en P vers l'infini
					R.normalize();
					Ray rayon_lumiere(P + 0.0001*N, R); //idem
					Vector Pprime, Nprime;
					int objetprime;
					double tprime;
					bool refletobjet = scene.intersect(rayon_lumiere, Pprime, Nprime, objetprime, tprime);
					if (refletobjet) {
						Vector PPprime = P - Pprime;
						double distance2 = PPprime.norm2();
						Vector PprimeL = L-Pprime;
						double distancelum = PprimeL.norm2()+distance1+distance2;
						Vector pixelColor = (I / 3.14*std::max(0., dot(Nprime, PprimeL)) / (4 * 3.14*distancelum))*scene.spheres[objetprime].color; // Luminosité fonction de la distance et de l'angle
						image[(i*W + j) * 3 + 0] = std::min(255., std::pow(pixelColor.x, 0.45));
						image[(i*W + j) * 3 + 1] = std::min(255., std::pow(pixelColor.y, 0.45));
						image[(i*W + j) * 3 + 2] = std::min(255., std::pow(pixelColor.z, 0.45));
					}
					else {
						image[(i*W + j) * 3 + 0] = 0;
						image[(i*W + j) * 3 + 1] = 0;
						image[(i*W + j) * 3 + 2] = 0;
					}
				}
				else {
					Vector PL = L-P;
					double distance = PL.norm2();
					PL.normalize();
					Ray rayon_lumiere(P + 0.00001*N, PL);
					Vector Pprime, Nprime;
					int objetprime;
					double tprime;
					bool ombre = scene.intersect(rayon_lumiere, Pprime, Nprime, objetprime, tprime);
					Vector PPprime = P - Pprime;
					double distanceprime = PPprime.norm2();
					if (ombre && distanceprime < distance) {
						image[(i*W + j) * 3 + 0] = 0;
						image[(i*W + j) * 3 + 1] = 0;
						image[(i*W + j) * 3 + 2] = 0;
					}
					else {
						Vector pixelColor = (I / 3.14*std::max(0., dot(N, PL)) / (4 * 3.14*distance))*SS.color; // Luminosité fonction de la distance et de l'angle
						image[(i*W + j) * 3 + 0] = std::min(255., std::pow(pixelColor.x, 0.45));
						image[(i*W + j) * 3 + 1] = std::min(255., std::pow(pixelColor.y, 0.45));
						image[(i*W + j) * 3 + 2] = std::min(255., std::pow(pixelColor.z, 0.45));
					}
				}

			}
			else {
				image[(i*W + j) * 3 + 0] = 0;
				image[(i*W + j) * 3 + 1] = 0;
				image[(i*W + j) * 3 + 2] = 0;

			}
		}
		}
		*/

		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);
	return 0;
};