#ifndef PHYSICS_HANDLER_H
#define PHYSICS_HANDLER_H
#include "Glew_Initialization.h"
#include "autodiff/forward/dual.hpp"
#include "FADBAD++/fadiff.h"

#define P0  101325 // initial pressure
#define T0  295. //initial temperature
#define R   8.31447 //gas const
   //Cp/Cv
#define m  29.e-3 //g/mol
  //for viscosity coefficient
#define I 1000 //amount of grids
//#define L 1.
//#define CFL 3.
//pistone acceleration


#define CFL 3

autodiff::dual H_z(autodiff::dual x, autodiff::dual y, autodiff::dual z) {
	return 1 + x + y + z + x * y + y * z + x * z + x * y * z + exp(x / y + y / z);/* expression for H_z */;
}

autodiff::dual fff(autodiff::dual x, autodiff::dual y, autodiff::dual z)
{
	

	return 1 + x + y + z + x * y + y * z + x * z + x * y * z + exp(x / y + y / z);
}

class PhysicsHandler
{
private:
	float Y = 1.4f;
	
public:
	float dx = 1e-3f;
	std::vector<float> new_Pressure;//
	std::vector<float> before_Pressure;//
	std::vector<float> new_Density;
	std::vector<float> before_Density;
	std::vector<float> new_Velocity;
	std::vector<float> before_Velocity;
	std::vector<float> new_X;
	std::vector<float> before_X;
	std::vector<double> Veloci;
	std::vector<double> magnetic;
	std::vector<double> Z_axis;
	std::vector<double> magnetic_fie;
	std::vector<double> current;
	std::vector<double> Electric;
	float k = 0;
	double dt, a0, pi, vis, A, B, C, meuu,dq;
	std::vector<float> XN, XB, uN, uB, PN, PB, rN, rB, eB ,eN, tB,tN;
	std::vector<float> x1, x2, x3, x4,xn;
	std::vector<float> u1, u2, u3, u4,un;
	std::vector<float> p1, p2, p3, p4,pn;
	std::vector<float> T1, T2, T3, T4, Tn;
	std::vector<float> pp1, pp2, pp3, pp4, ppp1,ppp2,ppp3,ppp4,longx;

	static void Step(Transform& transform)
	{


	}
	void reset()
	{
		XN.clear();
		PN.clear();
		PB.clear();
		rN.clear();
		XB.clear();
		uB.clear();
		uN.clear();
		rB.clear();
		eB.clear();
		eN.clear();
		for (int i = 0; i < I + 1; i++)
		{
			XN.push_back(0);
			PN.push_back(101325);
			PB.push_back(101325);
			rN.push_back(0);
			rB.push_back(0);
			XB.push_back(0);
			uB.push_back(0);
			uN.push_back(0);
			eB.push_back(0);
			eN.push_back(0);
		}
		XB.push_back(0);
		uB.push_back(0);
		uN.push_back(0);
		x1.clear();
		u1.clear();
		x2.clear();
		x3.clear();
		x4.clear();
		u2.clear();
		u3.clear();
		u4.clear();
		p1.clear();
		p2.clear();
		p4.clear();
		p3.clear();
		xn.clear();
		un.clear();
		pn.clear();
		T1.clear();
		T2.clear();
		T3.clear();
		T4.clear();
		Tn.clear();
		pp1.clear();
		pp2.clear();
		pp3.clear();
		pp4.clear();
		ppp1.clear();
		ppp2.clear();
		ppp3.clear();
		ppp4.clear();
		longx.clear();
	}
	PhysicsHandler()
	{
		XN.clear();
		PN.clear();
		PB.clear();
		rN.clear();
		XB.clear();
		uB.clear();
		uN.clear();
		tB.clear();
		tN.clear();
		rB.clear();
		eB.clear();
		eN.clear();
		for (int i = 0; i < I + 1; i++)
		{
			XN.push_back(0);
			PN.push_back(1.e5);
			PB.push_back(1.e5);
			rN.push_back(0);
			rB.push_back(0);
			XB.push_back(0);
			uB.push_back(0);
			uN.push_back(0);

			tB.push_back(295);
			tN.push_back(0);
			eB.push_back(0);
			eN.push_back(0);
		}
		XB.push_back(0);
		uB.push_back(0);
		uN.push_back(0);
		x1.clear();
		u1.clear();
		x2.clear();
		x3.clear();
		x4.clear();
		u2.clear();
		u3.clear();
		u4.clear();
		p1.clear();
		p2.clear();
		p4.clear();
		p3.clear();
		xn.clear();
		un.clear();
		pn.clear();
	}
	/*void Simulation()
	{
        
		new_Velocity[0] = aa * dt;
		int tM = 5;
		int t = 0;
		while (t <= dt)
		{
			for (int i = 0; i < I; i++)
			{
				before_Velocity[i] = new_Velocity[i];
				before_X[i] = new_Velocity[i];

			}
			for (int i = 0; i < I; i++)
			{

				new_X[i] = before_X[i] + before_Velocity[i] * dt;

			}
			for (int i = 0; i < I; i++)
			{
				if (new_X[i + 1] <= new_X[i])
				{
					std::cout << "i = %d\t" << i<<"\n";
					exit(0);
				}

				new_Density[i] = (before_Density[i] * (before_X[i + 1] - before_X[i])) / (new_X[i + 1] - new_X[i]);

			}
			for (int i = 0; i < I; i++)
			{
				new_Pressure[i] = before_Pressure[i] * powf(new_Density[i], gamma) / (powf(before_Density[i], gamma));
			}
			for (int i = 1; i < I; i++)
			{
				float temp = (new_X[i - 1] + new_X[i]) / 2 + (new_X[i + 1] + new_X[i]) / 2;
				new_Velocity[i] = before_Velocity[i] - (new_Pressure[i] - new_Pressure[i - 1]) / (((new_Density[i] + new_Density[i - 1]) / 2) * temp);
			}

			new_Velocity[0] = aa * dt * t;
			new_Velocity.back() = 0.0f;
			t++;


		}
	}*/

	void Simulation_2()
	{
		reset();
		
        float aa  = pisAcc;
		k = somek;
		dx = mydx; //x step 1mm
		a0 = sqrt(Y * R * T0 / m);
		dt =mydt;//dx/a0/CFL; //time step = 1.-6 [s]
		pi = 2. * asin(1.);

		//initial conditions begins
		for (int i = 0; i <= I; i++) {
			XN[i] = dx * i;
			uN[i] = 0.;
		}

		for (int i = 0; i < I; i++) {
			PN[i] = P0;
			rN[i] = P0 * m / R / T0; 
			eN[i] = PN[i] / rN[i] / (Y - 1.);
		}

		vis = rN[1] * dx * k; // viscosity
		//vis = 0;
		//initial conditions end

		//time loop begin
		for (int t = 1; t <= tmax; t++) {
			uN[0] = aa * dt * t;
			// switching times layers
			for (int i = 0; i <= I; i++) {
				XB[i] = XN[i];
				uB[i] = uN[i];
			}
			for (int i = 0; i < I; i++) {
				PB[i] = PN[i];
				rB[i] = rN[i];
				eB[i] = eN[i];
			} // 

			// New coordinates XN
			for (int i = 0; i < I; i++) {
				XN[i] = XB[i] + uB[i] * dt;
			}

			// New Density and New Pressure
			for (int i = 0; i < I; i++) {
				if (XN[i + 1] <= XN[i]) {
					printf("error t= %d, i= %d\n", t, i);
					exit(0);
				}
				rN[i] = rB[i] * (XB[i + 1] - XB[i]) / (XN[i + 1] - XN[i]);
				//PN[i] = pow(rN[i], Y) * PB[i] / pow(rB[i], Y);
				A = (uB[i + 1] - uB[i]) / (XN[i + 1] - XN[i]);
				eN[i] = eB[i] - dt * PB[i] * A / rN[i] + dt * vis * A * A / rN[i];
				PN[i] = eN[i] * rN[i] * (Y - 1.);
			}

			// New Velocity
			for (int i = 1; i < I; i++) {
				A = (XN[i + 1] + XN[i]) / 2. - (XN[i] + XN[i - 1]) / 2.;
				B = (uB[i + 1] - uB[i]) / (XB[i + 1] - XB[i]) - (uB[i] - uB[i - 1]) / (XB[i] - XB[i - 1]);
				C = (XB[i + 1] + XB[i]) / 2. - (XB[i] + XB[i - 1]) / 2.;
				uN[i] = uB[i] - dt * 2. * (PN[i] - PN[i - 1]) / (rN[i - 1] + rN[i]) / A + dt * vis * B * 2. / C / (rB[i - 1] + rB[i]);
			}

			

			if (t == 250) 
			{
				for (int i = 0; i < I; i++)
				{
					//fprintf(f1, "%g\t%g\t%g\n", XN[i], uN[i], PN[i] / P0);
					x1.push_back(XN[i]);
					u1.push_back(uN[i]);
					p1.push_back(PN[i] / P0);
					//std::cout << uN[i] << "\n";
				}
			}
			if (t == 500) {
				for (int i = 0; i < I; i++)
				{
					x2.push_back(XN[i]);
					u2.push_back(uN[i]);
					p2.push_back(PN[i] / P0);
				}
					//fprintf(f2, "%g\t%g\t%g\n", XN[i], uN[i], PN[i] / P0);
			}
			if (t == 750) {
				for (int i = 0; i < I; i++)
				{
					x3.push_back(XN[i]);
					u3.push_back(uN[i]);
					p3.push_back(PN[i] / P0);
				}
					//fprintf(f3, "%g\t%g\t%g\n", XN[i], uN[i], PN[i] / P0);
			}
			if (t == 1000) {
				for (int i = 0; i < I; i++)
				{
					x4.push_back(XN[i]);
					u4.push_back(uN[i]);
					p4.push_back(PN[i] / P0);
				}
					//fprintf(f4, "%g\t%g\t%g\n", XN[i], uN[i], PN[i] / P0);
			}
			if (t == myN * 1000)
			{
				for (int i = 0; i < I; i++)
				{
					xn.push_back(XN[i]);
					un.push_back(uN[i]);
					pn.push_back(PN[i] / P0);
				}

			}

			if (fmod(t, 100.) == 0.)
			{
				//std::cout<< "t= %d, t*dt= %g, XN[0]= %g\n" << t << "\t" << t * dt << "\t" << XN[0] << "\n";
			}
			

		} // time loop end
		
	}

	void Simulation_3()
	{
		reset();
		float aa = pisAcc;
		k = somek;
		dx = mydx; //x step 1mm
		a0 = sqrt(Y * R * T0 / m);
		dt = mydt;//dx/a0/CFL; //time step = 1.-6 [s]
		pi = 2. * asin(1.);

		//initial conditions begins
		for (int i = 0; i <= I; i++) {
			XN[i] = dx * i;
			uN[i] = 0.;
		}

		for (int i = 0; i < I; i++) {

			PN[i] =(i<I/2)? 3*P0:P0;

			rN[i] = P0 * m / R / T0;
			eN[i] = PN[i] / rN[i] / (Y - 1.);
		}

		vis = rN[1] * dx * k; // viscosity
		//vis = 0;
		//initial conditions end

		//time loop begin
		for (int t = 1; t <= tmax; t++) {
			uN[0] = aa * dt * t;
			// switching times layers
			for (int i = 0; i <= I; i++) {
				XB[i] = XN[i];
				uB[i] = uN[i];
			}
			for (int i = 0; i < I; i++) {
				PB[i] = PN[i];
				rB[i] = rN[i];
				eB[i] = eN[i];
			} // 

			// New coordinates XN
			for (int i = 0; i < I; i++) {
				XN[i] = XB[i] + uB[i] * dt;
			}

			// New Density and New Pressure
			for (int i = 0; i < I; i++) {
				if (XN[i + 1] <= XN[i]) {
					printf("error t= %d, i= %d\n", t, i);
					exit(0);
				}
				rN[i] = rB[i] * (XB[i + 1] - XB[i]) / (XN[i + 1] - XN[i]);
				//PN[i] = pow(rN[i], Y) * PB[i] / pow(rB[i], Y);
				A = (uB[i + 1] - uB[i]) / (XN[i + 1] - XN[i]);
				eN[i] = eB[i] - dt * PB[i] * A / rN[i] + dt * vis * A * A / rN[i];
				PN[i] = eN[i] * rN[i] * (Y - 1.);
			}

			// New Velocity
			for (int i = 1; i < I; i++) {
				A = (XN[i + 1] + XN[i]) / 2. - (XN[i] + XN[i - 1]) / 2.;
				B = (uB[i + 1] - uB[i]) / (XB[i + 1] - XB[i]) - (uB[i] - uB[i - 1]) / (XB[i] - XB[i - 1]);
				C = (XB[i + 1] + XB[i]) / 2. - (XB[i] + XB[i - 1]) / 2.;
				uN[i] = uB[i] - dt * 2. * (PN[i] - PN[i - 1]) / (rN[i - 1] + rN[i]) / A + dt * vis * B * 2. / C / (rB[i - 1] + rB[i]);
			}



			if (t == 250)
			{
				for (int i = 0; i < I; i++)
				{
					//fprintf(f1, "%g\t%g\t%g\n", XN[i], uN[i], PN[i] / P0);
					x1.push_back(XN[i]);
					u1.push_back(uN[i]);
					p1.push_back(PN[i] / P0);
					//std::cout << uN[i] << "\n";
					
				}
				
			}
			if (t == 500) {
				for (int i = 0; i < I; i++)
				{
					x2.push_back(XN[i]);
					u2.push_back(uN[i]);
					p2.push_back(PN[i] / P0);
				}
				//fprintf(f2, "%g\t%g\t%g\n", XN[i], uN[i], PN[i] / P0);
			}
			if (t == 750) {
				for (int i = 0; i < I; i++)
				{
					x3.push_back(XN[i]);
					u3.push_back(uN[i]);
					p3.push_back(PN[i] / P0);
				}
				//fprintf(f3, "%g\t%g\t%g\n", XN[i], uN[i], PN[i] / P0);
			}
			if (t == 1000) {
				for (int i = 0; i < I; i++)
				{
					x4.push_back(XN[i]);
					u4.push_back(uN[i]);
					p4.push_back(PN[i] / P0);
				}
				//fprintf(f4, "%g\t%g\t%g\n", XN[i], uN[i], PN[i] / P0);
			}
			if (t == myN * 1000)
			{
				for (int i = 0; i < I; i++)
				{
					xn.push_back(XN[i]);
					un.push_back(uN[i]);
					pn.push_back(PN[i] / P0);
				}

			}

			if (fmod(t, 100.) == 0.)
			{
				//std::cout<< "t= %d, t*dt= %g, XN[0]= %g\n" << t << "\t" << t * dt << "\t" << XN[0] << "\n";
			}


		} // time loop end

	}
	void Simulation_4()
	{
		reset();
		meuu = 6;
		float aa = pisAcc;
		k = somek;
		dx = mydx; //x step 1mm
		a0 = sqrt(Y * R * T0 / m);
		dt = mydt;//dx/a0/CFL; //time step = 1.-6 [s]
		pi = glm::pi<float>();
		dq = (P0 * m / R / T0) * dx;
		//initial conditions begins
		for (int i = 0; i <= I; i++) {
			XN[i] = dx * i;
			uN[i] = 0.;
		}

		for (int i = 0; i < I; i++) {

			PN[i] = (i < I / 2) ? 3 * P0 : P0;

			rN[i] = P0 * m / R / T0;
			eN[i] = PN[i] / rN[i] / (Y - 1.);
		}

		vis = rN[1] * dx * k; // viscosity
		//vis = 0;
		//initial conditions end

		//time loop begin
		for (int t = 1; t <= tmax; t++) {
			uN[0] = aa * dt * t;
			// switching times layers
			for (int i = 0; i <= I; i++) {
				XB[i] = XN[i];
				uB[i] = uN[i];
			}
			for (int i = 0; i < I; i++) {
				PB[i] = PN[i];
				rB[i] = rN[i];
				eB[i] = eN[i];

				tB[i] = T0;
			} // 

			// New coordinates XN
			for (int i = 0; i < I; i++) {
				XN[i] = XB[i] + uB[i] * dt;
			}

			


			// New Density and New Pressure
			for (int i = 0; i < I; i++) {
				if (XN[i + 1] <= XN[i]) {
					printf("error t= %d, i= %d\n", t, i);
					exit(0);
				}
				
				
				float d2e_dx2 = 0.0;
				if ( i > 0 && i < I - 1) {
					d2e_dx2 = (eB[i + 1] - 2.0 * eB[i] + eB[i - 1]) / (dx * dx);
				}
				float T = eN[i] * m * (Y - 1) / R;

				rN[i] = rB[i] * (XB[i + 1] - XB[i]) / (XN[i + 1] - XN[i]);
				//PN[i] = pow(rN[i], Y) * PB[i] / pow(rB[i], Y);
				A = (uB[i + 1] - uB[i]) / (XN[i + 1] - XN[i]);
				eN[i] = eB[i] - dt * PB[i] * A / rN[i] + dt * vis * A * A / rN[i];
				PN[i] = eN[i] * rN[i] * (Y - 1.);
				
			}

			// New Velocity
			for (int i = 1; i < I; i++) 
			{
				A = (XN[i + 1] + XN[i]) / 2. - (XN[i] + XN[i - 1]) / 2.;
				B = (uB[i + 1] - uB[i]) / (XB[i + 1] - XB[i]) - (uB[i] - uB[i - 1]) / (XB[i] - XB[i - 1]);
				C = (XB[i + 1] + XB[i]) / 2. - (XB[i] + XB[i - 1]) / 2.;
				uN[i] = uB[i] - dt * 2. * (PN[i] - PN[i - 1]) / (rN[i - 1] + rN[i]) / A + dt * vis * B * 2. / C / (rB[i - 1] + rB[i]);
			}



			if (t == 250)
			{
				for (int i = 0; i < I; i++)
				{
					//fprintf(f1, "%g\t%g\t%g\n", XN[i], uN[i], PN[i] / P0);*muee /(R*
					x1.push_back(XN[i]);
					u1.push_back(uN[i]);
					p1.push_back(PN[i] / P0);
					T1.push_back(p1[i]*m/(R* rN[i]));
					//std::cout << uN[i] << "\n";
				}

				float Mach = sqrt(((p1[5] / p1[I / 2])*(Y+1)+(Y-1))/(2*Y));
				float Mach2 = 2 + (Y - 1) * Mach * Mach / (2 * Y * Mach * Mach - (Y - 1));
				std::cout << "M1 : " << Mach<<"\t" << "M2 : " << Mach2 << "\n";
				float dist = ((0.73 - 0.62) * 1000 / (0.250)) / a0;

				float relative = (2 * Y) / (1 + Y) * dist * dist - (Y - 1) / (Y + 1);
				pp1.push_back(relative);
				ppp1.push_back(p1[5] /p1[I / 2] );
			
			}
			if (t == 500) {
				for (int i = 0; i < I; i++)
				{
					x2.push_back(XN[i]);
					u2.push_back(uN[i]);
					p2.push_back(PN[i] / P0);
					T2.push_back(p2[i] * m / (R * rN[i]));
					
				}
				//fprintf(f2, "%g\t%g\t%g\n", XN[i], uN[i], PN[i] / P0);
				float Mach = sqrt(((p2[5] / p2[I / 2]) * (Y + 1) + (Y - 1)) / (2 * Y));
				float Mach2 = 2 + (Y - 1) * Mach * Mach / (2 * Y * Mach * Mach - (Y - 1));
				std::cout << "M1 : " << Mach << "\t" << "M2 : " << Mach2 << "\n";
				float dist = ((0.73 - 0.62) * 1000 / (0.250))/a0;
			
				float relative = (2 * Y) / (1 + Y) * dist * dist - (Y - 1) / (Y + 1);
				//std::cout << p2[5]/ p2[I / 2] << "\n";
				//std::cout << relative << "\n";
				pp1.push_back(relative);
				ppp1.push_back(p2[5]/p2[I / 2] );


			}
			if (t == 750) {
				for (int i = 0; i < I; i++)
				{
					x3.push_back(XN[i]);
					u3.push_back(uN[i]);
					p3.push_back(PN[i] / P0);
					T3.push_back(p3[i] * m / (R * rN[i]));
				}
				//fprintf(f3, "%g\t%g\t%g\n", XN[i], uN[i], PN[i] / P0);
				float Mach = sqrt(((p3[5] / p3[I / 2]) * (Y + 1) + (Y - 1)) / (2 * Y));
				float Mach2 = 2 + (Y - 1) * Mach * Mach / (2 * Y * Mach * Mach - (Y - 1));
				std::cout << "M1 : " << Mach << "\t" << "M2 : " << Mach2 << "\n";
				float dist = ((0.73 - 0.62) * 1000 / (0.250)) / a0;

				float relative = (2 * Y) / (1 + Y) * dist * dist - (Y - 1) / (Y + 1);
				pp1.push_back(relative);
				ppp1.push_back(p3[5] / p3[I / 2]);
			}
			if (t == 1000) {
				for (int i = 0; i < I; i++)
				{
					x4.push_back(XN[i]);
					u4.push_back(uN[i]);
					p4.push_back(PN[i] / P0);
					T4.push_back(p4[i] * m / (R * rN[i]));
				}
				//fprintf(f4, "%g\t%g\t%g\n", XN[i], uN[i], PN[i] / P0);
				float Mach = sqrt(((p4[5] / p4[I / 2]) * (Y + 1) + (Y - 1)) / (2 * Y));
				float Mach2 = 2 + (Y - 1) * Mach * Mach / (2 * Y * Mach * Mach - (Y - 1));
				std::cout << "M1 : " << Mach << "\t" << "M2 : " << Mach2 << "\n\n";
				float dist = ((0.73 - 0.62) * 1000 / (0.250)) / a0;

				float relative = (2 * Y) / (1 + Y) * dist * dist - (Y - 1) / (Y + 1);
				pp1.push_back(relative);
				ppp1.push_back(p4[5] / p4[I / 2]);
				
				for (int i = 0; i < 4; i++)
				{
					longx.push_back(i);
				}
				std::cout << "p2/p1 simulated : " << ppp1[0] << "\t" << "p2/p1 calculated : " << pp1[0] << "\n";
				std::cout << "p2/p1 simulated : " << ppp1[1] << "\t" << "p2/p1 calculated : " << pp1[1] << "\n";
				std::cout << "p2/p1 simulated : " << ppp1[2] << "\t" << "p2/p1 calculated : " << pp1[2] << "\n";
				std::cout << "p2/p1 simulated : " << ppp1[3] << "\t" << "p2/p1 calculated : " << pp1[3] << "\n\n";




			}
			if (t == myN * 1000)
			{
				for (int i = 0; i < I; i++)
				{
					xn.push_back(XN[i]);
					un.push_back(uN[i]);
					pn.push_back(PN[i] / P0);
				}

			}
			
			if (fmod(t, 100.) == 0.)
			{
				//std::cout<< "t= %d, t*dt= %g, XN[0]= %g\n" << t << "\t" << t * dt << "\t" << XN[0] << "\n";
			}


		} // time loop end

	}
	void Poisseulle()
	{

		Z_axis.clear();
		Veloci.clear();
		magnetic.clear();
		magnetic_fie.clear();
		current.clear();
		Electric.clear();
		for (float i = -range; i < range+0.1; i = i+0.1)
		{
			Z_axis.push_back(i);
			Veloci.push_back(Vz(i));
			magnetic.push_back(Vzf(i));
			magnetic_fie.push_back(Hz(i));
			current.push_back(Jz(i));
			Electric.push_back(Ez(i));
		}
	}

	double Vzf(float z)
	{
		//return -powf(l_of_pipe,2)*grad_p()* (1 - powf(z/l_of_pipe, 2)) / (2 * 0.01f);
		return -powf(l_of_pipe, 2) * grad_p() * (1 - powf(z / l_of_pipe, 2)) / (2 * viscousity);
	}
	double Vx()
	{
		//return -powf(l_of_pipe,2)*grad_p()* (1 - powf(z/l_of_pipe, 2)) / (2 * 0.01f);
		float x = (cosh(Ha) - 1) / (Ha * sinh(Ha));
		return -powf(l_of_pipe, 2) * grad_p() *x / (viscousity);
	}
	double Hz(float z)
	{
		float x = (sinh(Ha * z / l_of_pipe) - (z / l_of_pipe) * sinh(Ha)) / (cosh(Ha) - 1);
		return Vx() * 4 * sqrt(sigma * viscousity) * glm::pi<float>() * x / speed_light;
	}
	double Vz(float z)
	{
		return Vx() * (cosh(Ha) - cosh(Ha * z / l_of_pipe)) / (cosh(Ha) - 1);
	}
	double grad_p()
	{
		return 8 * viscousity * l_of_pipe * flow_rate / (glm::pi<float>() *powf(pipe_raduis, 4));
	}
	double Jz(float z)
	{
		return Cz(z).y;
	}
	
	glm::dvec3 Mag_field(float z)
	{
		double mini_H = (speed_light / l_of_pipe) * sqrt(viscousity / sigma) * Ha;
		return glm::dvec3(Hz(z), 0.00001, mini_H);
	}
	

	glm::dvec3 Cz(float Z)
	{
		
		autodiff::dual x = Mag_field(Z).x;
		autodiff::dual y = Mag_field(Z).y;
		autodiff::dual z = Mag_field(Z).z;

		autodiff::dual curl_H_x_val = 0;
		autodiff::dual curl_H_y_val = autodiff::derivative(H_z, autodiff::wrt(y), autodiff::at(x, y, z));;
		autodiff::dual curl_H_z_val = 0;

		autodiff::dual u = fff(x, y, z);
		//std::cout << "Before derivative computation:" << std::endl;
		//std::cout << "x: " << x << ", y: " << y << ", z: " << z << std::endl;
		//std::cout << "fff: " << u.val << std::endl;

		//autodiff::dual dudx = autodiff::derivative(fff, autodiff::wrt(x), autodiff::at(x, y, z));
		//autodiff::dual dudy = autodiff::derivative(fff, autodiff::wrt(y), autodiff::at(x, y, z));
		//autodiff::dual dudz = autodiff::derivative(fff, autodiff::wrt(z), autodiff::at(x, y, z));



		//std::cout << "After derivative computation:" << std::endl;
		//std::cout << "dudx: " << dudx.val << ", dudy: " << dudy.val << ", dudz: " << dudz.val << std::endl;

		double temp = speed_light / (4 * glm::pi<float>());
		//double dudx_float = dudx.val;
		//double dudy_float = dudy.val;
		//double dudz_float = dudz.val;

		//double dudx_float = curl_H_x_val.val;
		//double dudy_float = curl_H_y_val.val;
		//double dudz_float = curl_H_z_val.val;


		
		
		//return glm::dvec3(temp * dudx_float, temp * dudy_float, temp * dudz_float);
		
		glm::dvec3 curl_H = glm::cross(Mag_field(Z), glm::dvec3(0, 0, 1));
		//glm::dvec3 curl_H = glm::dvec3(curl_H_x_val, curl_H_y_val, curl_H_z_val);
		return  ((double)speed_light / (4 * glm::pi<double>()) * curl_H);

	}
	double Ez(float z)
	{
		float mini_H = (speed_light / l_of_pipe) * sqrt(viscousity / sigma) * Ha;
		glm::dvec3 bbb = (Cz(z) - (glm::cross((1 / (double)speed_light) * glm::dvec3(Vz(z) / mini_H, 0, 0), Mag_field(z)))) / (double)sigma;
		return bbb.y;
	}
};
#endif