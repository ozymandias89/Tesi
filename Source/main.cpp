/*
 @file    main.cc
 @author  Riccardo Zanella, riccardozanella89@gmail.com
 @version 1.0
 */

/* --*- C++ -*------x-----------------------------------------------------------
 *
 *
 * Description: This main resolve a linear relaxation of the problem in the form
 * 				min cT x
 * 				s.t Ax=b
 * 				x€R+
 *
 * 				then print the higher fractional variable and create a branch
 * 				with P1 and P2 new problem to resolve.
 * 				In the end create problem step 7.
 *
 * -----------------x-----------------------------------------------------------
 */


#include "solve.cpp"

int main(int argc, char const *argv[]) {

	if (argc < 2) { // you must insert file
		throw std::runtime_error("usage: ./main filename.txt");

	}

	ifstream myfile(argv[1], std::ios::in);

	load_problem(myfile);

	myfile.close();

	try {

		// --------------------------------------------------
		// 1. Initialization problem
		// --------------------------------------------------
		DECL_ENV(env);
		DECL_PROB(env, lp, "resolve problem RL");
		setupLP(env, lp);


		// --------------------------------------------------
		// 2. program (first part)
		// --------------------------------------------------
		solve (env, lp);

		int num_rows = (CPXgetnumrows(env, lp));


		// ---------------------------------------------------------
		// 3. if P_1 and P_2 have solution
		// ---------------------------------------------------------
		if (!flag_find) {
			// --------------------------------------------------
			// 4. Initialization of the second problem
			// --------------------------------------------------
			cout << "Now we have solved both subproblems and we have the min z. " << endl;
			cout << "Your guide to all the data is shown here: " << endl;

			DECL_ENV(env_dual);
			DECL_PROB(env_dual, lp_dual, "resolve second problem");
			setupSP(env_dual, lp_dual, num_rows);

			print_matrix();
			print_vect_b();
			print_cut_A();
			print_cut_b();

			print_u_variables();
			cout << endl;
			print_v_variables();
			cout << endl;
			cout << "gamma= " << gam << endl;
			cout << endl;
			cout << "min solution= " << min_sol << endl;
			cout << endl;
			cout << "index fractional variable (e_k)= " << k << endl;

			cout << "It was also set the second problem... check the correctness of the second model (second_problem.lp)." << endl;
			cout << "Warning: In the first primal problem (problem.lp) you don't evaluate the last line (last constraint).  " << endl;
			cout << " The last constraint, indeed, is trash." << endl;

			CHECKED_CPX_CALL(CPXwriteprob, env_dual, lp_dual, "../data/second_problem.lp", 0);

			CPXfreeprob(env_dual, &lp_dual);
			CPXcloseCPLEX(&env_dual);
		}

		// ---------------------------------------------------------
		// 5. free allocate memory
		// ---------------------------------------------------------


		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env);


		for (int i = 0; i < num_constraint; ++i) {
			delete[] A[i];
		}
		delete[] A;
		delete[] b;
		delete[] c;


			for (std::vector<double*>::const_iterator j = cut_A.begin();
					j != cut_A.end(); ++j) {
				delete *j;
			}


	} catch (std::exception& e) {
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
	}

	return 0;
}
