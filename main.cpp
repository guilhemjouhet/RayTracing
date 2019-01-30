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
Vector operator*(const Vector& a, const Vector& b) {
	return Vector(a.x * b.x, a.y * b.y, a.z * b.z);
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

Vector cross(const Vector& a, const Vector& b) {
	return Vector(a.x*b.y - a.y*b.x, a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z);
}

class Ray {
public:
	Ray(const Vector& C, const Vector& U) : C(C), U(U) {};
	Vector C;
	Vector U;
};

class Sphere{
public:
	Sphere(const Vector& O, double R, const Vector& color, const bool miroir, const bool lumiere, const bool transparent = false) : O(O), R(R), color(color), miroir(miroir),lumiere(lumiere), transparent(transparent) {};
	Vector O;
	double R;
	Vector color;
	bool miroir;
	bool lumiere;
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

	Vector getcolor(const Ray& ray, const Vector& L, int nrebond, int nChemin=1) {
		if (nrebond == 0 || nChemin==0) { 
			return Vector(0, 0, 0); 
		}
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
				Vector color = getcolor(rayon_lumiere, L, nrebond-1,nChemin);
				return color;
			}
			else if (SS.transparent) {
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
			else if (SS.lumiere) {
				Vector color = I/4./3.14/3.14/SS.R/SS.R*SS.color;
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
				if (ombre && distanceprime < distance && !spheres[objetprime].transparent && !spheres[objetprime].lumiere) {
					//on laisse color à 0
				}
				else if (spheres[objetprime].lumiere){
					// direct cas ponctuel
					color = (I / 3.14*std::max(0., dot(N, PL)) / (4 * 3.14*distance))*SS.color;

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
						color = color + 1/ (double)nChemin/255.*getcolor(rayon_aleatoire, L, nrebond-1,1)*SS.color;// on met nrebond=1 pour pas passer en exponentiel
					}
				}
				else {
					Vector R = ray.U - 2 * dot(ray.U, N)*N; //rayon réfléchit en P vers l'infini
					R.normalize();
					Ray rayon_lumiere(P + 0.0001*N, R); //idem
					color = color +  1./255.*getcolor(rayon_lumiere, L, nrebond - 1, 1)*SS.color;
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
	Vector L(-10,20, 40);
	Vector colorL(255, 255, 255);
	Sphere Lum(L, 10, colorL, false, true);//lumiere
	scene.addSphere(Lum);
	//double I = 10000; // Intensité de la source
	int W = 720;
	int H = 480;
	int R = 990;
	double fov = 65 * 3.1415 / 180;
	double tanhalfov = -W / 2 / tan(fov / 2);
	std::vector<unsigned char> image(W*H * 3, 0);
	Vector O(0, -15, 0); // centre de la sphère
	Vector color(255, 0,0);
	//Sphere S(O, 10, color, false, false);
	Sphere S(O, 10, color, true,false); //miroir
	scene.addSphere(S);
	Vector OO(0, 15,-20); // centre de la sphère transparente
	//Sphere St(OO, 10, color, false, false);
	Sphere St(OO, 10, color, false,false,true); //Sphère transparente
	scene.addSphere(St);
	Vector colorb(255, 10, 255);
	Vector Ob(1000, 0, 0); // centre de la sphère du bas
	Sphere Sb(Ob, R, colorb, false,false);
	scene.addSphere(Sb);
	Vector colorh(0, 100,100);
	Vector Oh(-1200, 0, 0); // centre de la sphère du haut
	Sphere Sh(Oh, R, colorh, false, false);
	scene.addSphere(Sh);
	Vector colorder(255, 25,255);
	Vector Oder(0, 0, 1000); // centre de la sphère de derrière
	Sphere Sder(Oder, 940, colorder, false,false);
	scene.addSphere(Sder);
	Vector colorg(20, 150,30);
	Vector Og(0, 1000, 0); // centre de la sphère de gauche
	Sphere Sg(Og, 970, colorg, false,false);
	scene.addSphere(Sg);
	Vector colord(100, 10, 10);
	Vector Od(0, -1000, 0); // centre de la sphère de gauche
	Sphere Sd(Od, 970, colord, false, false);
	scene.addSphere(Sd);
	Vector colorf(0, 0,255);
	Vector Of(0, 0, -1200); // centre de la sphère du fond
	Sphere Sf(Of, 980, colorf, false,false);
	scene.addSphere(Sf);
	//Sphere SS(O, 10, color, false,false);
	int nChemin =  50;
	int nrebond =  10;
	int nrays =  50;
	double distnette = 55;
	#pragma omp parallel for schedule(dynamic,12)
	for (int i = 0; i < H; i++) { //pt pixel
		if (i%50==0) std::cout << i/(double)H*100 <<"%\n";
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