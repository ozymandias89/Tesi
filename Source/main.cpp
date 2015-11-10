/*
 @file    main.cc
 @author  Riccardo Zanella, riccardozanella89@gmail.com
 @version 2.0
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
#include "SecondProblem.cpp"

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

			print_matrix();

			SecondProblem* sec_prob = new SecondProblem();

			sec_prob->setupSP(env_dual, lp_dual);

			sec_prob->solve(env_dual, lp_dual, false);

			// --------------------------------------------------
			// 5. Evaluate vector r
			// --------------------------------------------------
			sec_prob->evaluate_rT();

			print_matrix();

			// --------------------------------------------------
			// 6. Cycle
			// --------------------------------------------------

			//while (Y!=Ytilde)
			sec_prob->step8_1(env_dual, lp_dual);

			sec_prob->step8_2(env_dual, lp_dual);


			sec_prob->solve(env_dual, lp_dual, true);



			// --------------------------------------------------
			// 6. Show dates
			// --------------------------------------------------
			cout << endl;
			print_vect_b();
			cout << endl;
			print_vect_c();
//			cout << endl;
//			print_u_variables();
//			cout << endl;
//			print_v_variables();
//			cout << endl;
			cout << "gamma= " << gam << endl;
			cout << endl;
			cout << "min solution= " << min_sol << endl;
			cout << endl;
			cout << "index fractional variable (e_k)= " << k << endl;

			CHECKED_CPX_CALL(CPXwriteprob, env_dual, lp_dual, "../data/second_problem.lp", 0);

			CPXfreeprob(env_dual, &lp_dual);
			CPXcloseCPLEX(&env_dual);
			free(sec_prob);
		}

		// ---------------------------------------------------------
		// 5. free allocate memory
		// ---------------------------------------------------------


		CPXfreeprob(env, &lp);
		CPXcloseCPLEX(&env);


	} catch (std::exception& e) {
		std::cout << ">>>EXCEPTION: " << e.what() << std::endl;
	}

	return 0;
}
