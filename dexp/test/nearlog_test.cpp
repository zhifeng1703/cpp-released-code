#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <chrono>
// #include <string>
#include <map>

#include <fstream>
#include <iostream>

#include "matOp.hpp"
#include "logmNear.hpp"
#include "skewSchFac.hpp"

#include "basics.hpp"

#define num_main_option 8

int main(int argc, char *argv[])
{
	string def_val[num_main_option] = {"10", "10", "10", "200", "1.0", "9527", "10", "../figures/dexp"};
	string options[num_main_option] = {"-dim", "-stack", "-minscale", "-maxiter", "-length", "-loop", "-seed", "-file"};
	char char_buff[100];

	auto paras = read_paras(argc, argv, num_main_option, def_val, options);

	INTE_TYPE dim_beg = 5;
	INTE_TYPE dim_end = stoi(paras["-dim"]);
	INTE_TYPE seed = stoi(paras["-seed"]);
	INTE_TYPE stack = stoi(paras["-stack"]);
	INTE_TYPE scale = stoi(paras["-minscale"]);
	INTE_TYPE queue = stoi(paras["-maxiter"]);

	REAL_TYPE length = stoi(paras["-length"]);
	REAL_TYPE radius = 0.0;
	REAL_TYPE dist = 0.0;

	INTE_TYPE stack_ind = 0;
	INTE_TYPE scale_ind = 0;
	INTE_TYPE queue_ind = 0;
	INTE_TYPE max_scale = 0;

	INTE_TYPE record_cnt = 0;
	INTE_TYPE loops = stoi(paras["-loop"]);
	for (INTE_TYPE outerloop = 0; outerloop < 2; outerloop++)
	{
		INTE_TYPE cur_loop = outerloop == 0 ? 1 : loops;
		if (outerloop)
			std::printf("\nPreheat finished! Experiments have started!\n");
		else
			std::printf("Preheating...");

		std::default_random_engine rng;
		rng.seed(seed);

		for (INTE_TYPE dim_ind = 0; dim_ind < dim_end - dim_beg + 1; dim_ind++)
		{
			INTE_TYPE n = dim_beg + dim_ind;
			INTE_TYPE m = n / 2;

			radius = 0;
			dist = 0;

			stack_ind = 0;
			scale_ind = 0;
			queue_ind = 0;
			max_scale = 0;

			if (outerloop)
				printf("Testing Matrices of dimension\t %lld \tx\t %lld...\n", n, n);

			auto nearlog = LogmSkewSymm(n);

			auto MatA = SkewSymmMat(n);
			auto MatM = SkewSymmMat(n);

			MatA.init_low_vec();
			MatM.init_low_vec();

			MatA.Rand(rng, -1.0, 1.0);
			MatA.low2upp();
			MatM.Rand(rng, -length, length);
			MatM.low2upp();

			auto MateA = ColMat<REAL_TYPE>(n, n);
			auto MateM = ColMat<REAL_TYPE>(n, n);

			auto MatQ = ColMat<REAL_TYPE>(n, n);

			auto StackScalar = ArrVec<INTE_TYPE>(queue);
			auto StackSOMatQ = ArrVec<ColMat<REAL_TYPE> *>(stack);
			auto StackSchVec = ArrVec<ColMat<REAL_TYPE> *>(stack);
			auto StackAngles = ArrVec<ArrVec<REAL_TYPE> *>(stack);

			for (auto ind = 0; ind < stack; ind++)
			{
				StackSOMatQ.v[ind] = new ColMat<REAL_TYPE>(n, n);
				StackSchVec.v[ind] = new ColMat<REAL_TYPE>(n, n);
				StackAngles.v[ind] = new ArrVec<REAL_TYPE>(m);
			}

			auto VecMateM = ArrVec<ColMat<REAL_TYPE> *>(scale);

			auto QueueRadius = ArrVec<REAL_TYPE>(queue);
			auto QueueSSMatA = ArrVec<ColMat<REAL_TYPE> *>(queue);
			auto QueueSOMatQ = ArrVec<ColMat<REAL_TYPE> *>(queue);
			auto QueueSchVec = ArrVec<ColMat<REAL_TYPE> *>(queue);
			auto QueueAngles = ArrVec<ArrVec<REAL_TYPE> *>(queue);

			for (auto ind = 0; ind < queue; ind++)
			{
				QueueSSMatA.v[ind] = new ColMat<REAL_TYPE>(n, n);
				QueueSOMatQ.v[ind] = new ColMat<REAL_TYPE>(n, n);
				QueueSchVec.v[ind] = new ColMat<REAL_TYPE>(n, n);
				QueueAngles.v[ind] = new ArrVec<REAL_TYPE>(m);
			}

			nearlog.SSF.Factor_SkewSymm(MatA.v, n);
			nearlog.SSF.Explict_Vector();
			nearlog.SSF.Exponential(MateA);

			radius = nearlog.SSF.Dist2ConjugateLocus();
			if (!outerloop)
				std::printf("\nRadius of the inscribed ball:\t %f\n", radius);

			QueueRadius[queue_ind] = radius;
			QueueSSMatA[queue_ind]->Assign(MatA);
			QueueSOMatQ[queue_ind]->Assign(MateA);
			QueueSchVec[queue_ind]->Assign(nearlog.SSF.R);
			QueueAngles[queue_ind]->Assign(nearlog.SSF.A);

			nearlog.SSF.Factor_SkewSymm(MatM.v, n);
			nearlog.SSF.Explict_Vector();
			for (auto ind = 0; ind < scale; ind++)
			{
				VecMateM.v[ind] = new ColMat<REAL_TYPE>(n, n);
				nearlog.SSF.Exponential(*VecMateM.v[ind]);
				for (auto angle_ind = 0; angle_ind < m; angle_ind++)
					nearlog.SSF.A.v[angle_ind] *= 0.5;
			}

			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, MateA.v, n, VecMateM.v[0]->v, n, 0.0, MatQ.v, n);
			StackSOMatQ[stack_ind]->Assign(MatQ.v, n);
			StackScalar[stack_ind] = 0;
			nearlog.SSF.Factor_SpecOrth(StackSOMatQ[stack_ind]->v, n);
			nearlog.SSF.Explict_Vector();
			StackSchVec[stack_ind]->Assign(nearlog.SSF.R);
			StackAngles[stack_ind]->Assign(nearlog.SSF.A);

			while (stack_ind >= 0)
			{
				nearlog._diff(QueueSSMatA[queue_ind]->v, n, StackSchVec[stack_ind]->v, n, StackAngles[stack_ind]->v);
				dist = nearlog._norm(nearlog.SSF);
				if (!outerloop)
					std::printf("Distance to the center:\t\t\t %f\n", dist);

				if (dist < QueueRadius[queue_ind])
				{
					queue_ind++;
					if (queue_ind >= queue)
					{
						std::cout << "Error: Nearlog not found within the given number of iterations!!!\n";
						for (auto ind = 0; ind < queue; ind++)
							QueueSSMatA[ind]->printf((std::snprintf(char_buff, sizeof(char_buff), "Iteration %d:\n", ind), char_buff));
						throw(1);
					}
					nearlog.SSF.R.Assign(*StackSchVec[stack_ind]);
					nearlog.SSF.A.Assign(*StackAngles[stack_ind]);
					radius = nearlog.SSF.Dist2ConjugateLocus();
					if (!outerloop)
						std::printf("\nRadius of the inscribed ball:\t %f\n", radius);

					QueueRadius[queue_ind] = radius;
					nearlog.SSF.GetSkewSymm(QueueSSMatA[queue_ind]->v, n);
					skewl2m(*QueueSSMatA[queue_ind]);
					QueueSOMatQ[queue_ind]->Assign(*StackSOMatQ[stack_ind]);
					QueueSchVec[queue_ind]->Assign(nearlog.SSF.R);
					QueueAngles[queue_ind]->Assign(nearlog.SSF.A);

					stack_ind--;
				}
				else
				{
					StackScalar[stack_ind] += 1;
					scale_ind = StackScalar[stack_ind];
					if (scale_ind > max_scale)
						max_scale = scale_ind;
					stack_ind++;
					if (stack_ind >= stack || scale_ind > scale)
					{
						std::cout << "Error: Stack Overflow in nearlog!!!\n";
						throw(1);
					}
					StackScalar[stack_ind] = scale_ind;
					cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, QueueSOMatQ[queue_ind]->v, n, VecMateM[scale_ind]->v, n, 0.0, StackSOMatQ[stack_ind]->v, n);

					nearlog.SSF.Factor_SpecOrth(StackSOMatQ[stack_ind]->v, n);
					nearlog.SSF.Explict_Vector();
					StackSchVec[stack_ind]->Assign(nearlog.SSF.R);
					StackAngles[stack_ind]->Assign(nearlog.SSF.A);
				}
			}

			printf("Total counts of logarithms:\t %lld. Maximum order of scaling:\t %lld.\n", queue_ind, max_scale);

			// for (auto ind = 0; ind < queue_ind + 1; ind++)
			//{
			//	QueueSSMatA[ind]->printf((std::snprintf(char_buff, sizeof(char_buff), "Iteration %d,\tSkew Symmetric A with distance bound\t%f:\n", ind, QueueRadius[ind]), char_buff));
			//	// QueueSOMatQ[ind]->printf((std::snprintf(char_buff, sizeof(char_buff), "Iteration %d,\tSpectial Orthogonal Q:\n", ind), char_buff));
			// }

			// Release array
			for (auto ind = 0; ind < stack; ind++)
			{
				delete StackSOMatQ.v[ind];
				delete StackSchVec.v[ind];
				delete StackAngles.v[ind];
			}

			for (auto ind = 0; ind < scale; ind++)
				delete VecMateM.v[ind];

			for (auto ind = 0; ind < queue; ind++)
			{
				delete QueueSSMatA.v[ind];
				delete QueueSOMatQ.v[ind];
				delete QueueSchVec.v[ind];
				delete QueueAngles.v[ind];
			}
		}
	}
	return 1;
};