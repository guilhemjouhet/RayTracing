#define _CRT_SECURE_NO_WARNINGS 1

#include <Vector>
#include <algorithm>
#include <stdlib.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
//#define double PI = 3.1415
#include "stb_image_write.h"

//R=ray.u-2*dot(ray.u-n)*N

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
};

Vector operator+(const Vector& a, const Vector& b){
	return Vector(a.x+b.x,a.y+b.y,a.z+b.z);
}

Vector operator-(const Vector& a, const Vector& b){
	return Vector(a.x-b.x,a.y-b.y,a.z-b.z);
}

Vector operator+(const Vector& a, const double b){
	return Vector(a.x+b,a.y+b,a.z+b);
}

Vector operator*(const double b,const Vector& a){
	return Vector(b*a.x,b*a.y,b*a.z);
}

Vector operator/(const Vector& a, const double b){
	return Vector(b/a.x,b/a.y,b/a.z);
}

double dot(const Vector& a, const Vector& b){
	return a.x*b.x+a.y*b.y+a.z*b.z;
}

class Ray {
public:
	Ray(const Vector& C, const Vector& U) : C(C), U(U) {};
	Vector C;
	Vector U;
};

class Sphere{
public:
	Sphere(const Vector& O, double R, const Vector& color,const bool miroir, const bool transparent=false) : O(O), R(R), color(color),miroir(miroir),transparent(transparent) {};
	Vector O;
	double R;
	Vector color;
	bool miroir;
	bool transparent;

	bool intersect(const Ray& r, Vector &P, Vector &N, double &t) {
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

class Scene {
public:
	Scene() {};
	std::vector<Sphere> spheres;

	void addSphere(const Sphere& s) {
		spheres.push_back(s);
	}

	bool intersect(const Ray& r, Vector &P, Vector &N, int &ii, double &tprime) {
		double smallestt = 1E15;
		for (int i = 0; i < spheres.size(); i++) {
			Vector Plocal, Nlocal;
			double t = 0;
			bool inter = spheres[i].intersect(r, Plocal, Nlocal, t);
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

	Vector getcolor(const Ray& ray, const Vector& L) {
		double I = 10000000;
		Vector P, N;//création des vecteurs d'intesection : point d'intersection et normale à la sphère
		int ii = -1;
		double t;
		if (intersect(ray, P, N, ii, t)) { //on trouve un objet
			Sphere SS = spheres[ii];

			if (SS.miroir) {
				Vector R = ray.U - 2 * dot(ray.U, N)*N; //rayon réfléchit en P vers l'infini
				R.normalize();
				Ray rayon_lumiere(P + 0.0001*N, R); //idem
				Vector color = getcolor(rayon_lumiere, L);
				return color;
			}
			else if (SS.transparent) {
				float nairnverre = 1 / 0.8;
				if (dot(ray.U, N) > 0) {
					N = Vector(0,0,0)-N;
					nairnverre = 1 / nairnverre;
				}
				float det = 1 - nairnverre * nairnverre*(1 - dot(ray.U, N)*dot(ray.U, N));
				if (det < 0) {
					Vector R = ray.U - 2 * dot(ray.U, N)*N; //rayon réfléchit en P vers l'infini
					R.normalize();
					Ray rayon_lumiere(P + 0.0001*N, R); //idem
					Vector color = getcolor(rayon_lumiere, L);
					return color;

				}
				else {
					Vector R = nairnverre * ray.U - (nairnverre*dot(ray.U, N) + sqrt(det))*N;
					R.normalize();
					Ray rayon_lumiere(P - 0.0001*N, R); //idem
					Vector color = getcolor(rayon_lumiere, L);
					return color;
				}
			}
			else {
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
				if (ombre && distanceprime < distance && !spheres[objetprime].transparent) {
					Vector color(0, 0, 0);
					return color;
				}
				else {
					Vector color = (I / 3.14*std::max(0., dot(N, PL)) / (4 * 3.14*distance))*SS.color; // Luminosité fonction de la distance et de l'angle
					return color;
				}
			}
		}
		else {
			Vector color(0, 0, 0);
			return color;
		}
	}
};

int main() {
	Vector C(0, 0, 55);
	Vector L(-10,20, 40);
	double I = 1000000; // Intensité de la source
	int W = 720;
	int H = 720;
	int R = 990;
	double fov = 80 * 3.1415 / 180;
	double tanhalfov = -W / 2 / tan(fov / 2);
	std::vector<unsigned char> image(W*H * 3, 0);
	Vector O(0, 15, 0); // centre de la sphère
	Vector color(255, 0, 0);
	Sphere S(O, 10, color, true); //Sphère centrée en O de rayon R=10
	Scene scene;
	scene.addSphere(S);
	Vector OO(0, -15, 0); // centre de la sphère transparente
	Sphere St(OO, 10, color, false,true); //Sphère transparente
	scene.addSphere(St);
	Vector colorb(0, 255, 0);
	Vector Ob(1000, 0, 0); // centre de la sphère du bas
	Sphere Sb(Ob, R, colorb, false);
	scene.addSphere(Sb);
	Vector colorder(0, 255, 255);
	Vector Oder(0, 0, 1000); // centre de la sphère de derrière
	Sphere Sder(Oder, 940, colorder, false);
	scene.addSphere(Sder);
	Vector colorg(255, 255, 0);
	Vector Og(0, 1000, 0); // centre de la sphère de gauche
	Sphere Sg(Og, 960, colorg, false);
	scene.addSphere(Sg);
	Vector colorf(255, 255, 255);
	Vector Of(0, 0, -1000); // centre de la sphère du fond
	Sphere Sf(Of, R, colorf, false);
	scene.addSphere(Sf);
	Sphere SS(O, 10, color, false);
	for (int i = 0; i < H; i++) { //pt pixel
		for (int j = 0; j < W; j++) {
			Vector U(i - H / 2, -j + W / 2, tanhalfov); //calcul du vecteur directeur du rayon
			U.normalize();
			Ray ray(C, U);
			Vector color = scene.getcolor(ray, L);
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