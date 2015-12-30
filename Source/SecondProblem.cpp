/*
 * SecondProblem.cpp
 *
 *  Created on: 04 nov 2015
 *      Author: riccardo
 */

#include "SecondProblem.h"

SecondProblem::SecondProblem(bool verbose) {
	this->violated_constraint = false;
	this->verbose = verbose;
	this->cost = c;
}

SecondProblem::~SecondProblem() {
	// TODO Auto-generated destructor stub
}

void SecondProblem::print_u() {

	cout << endl;
	for (unsigned int i = 0; i < u.size(); ++i)
		cout << "u_" << i + 1 << " " << " = " << u[i] << endl;

}

void SecondProblem::print_v() {

	cout << endl;
	for (unsigned int i = 0; i < v.size(); ++i)
		cout << "v_" << i + 1 << " = " << v[i] << endl;

}

void SecondProblem::print_r() {
	cout << endl;
	cout << "vector r " << endl;
	for (std::vector<double>::const_iterator j = rt.begin(); j != rt.end(); ++j)
		cout << *j << " ";
	cout << endl << endl;
}

void SecondProblem::print_c() {
	cout << endl;
	cout << "vector c " << endl;
	for (std::vector<double>::const_iterator j = cost.begin(); j != cost.end();
			++j)
		cout << *j << " ";
	cout << endl << endl;
}

void SecondProblem::print_a() {

	cout << endl;
	for (unsigned int i = 0; i < a.size(); ++i)
		cout << "a_" << i << " " << " = " << a[i] << endl;

}

void SecondProblem::print_u0() {

	cout << "u_0" << " = " << u0 << endl;
}

void SecondProblem::print_v0() {
	cout << endl;
	cout << "v0 " << " = " << v0 << endl;
}

void SecondProblem::print_beta() {
	cout << endl;
	cout << "beta " << " = " << beta << endl;

}

void SecondProblem::print_y_tilde() {

	cout << endl;
	cout << "y_tilde= ";
	for (unsigned int i = 0; i < y_tilde.size(); ++i)
		cout << y_tilde[i] << " ";

	cout << endl;
}

void SecondProblem::set_up_variable(CEnv env, Prob lp) {

	// ----------------------------------------------------
	// second problem order of variable: (u_0 u a b v_0 v)
	// ----------------------------------------------------

	if (verbose) {
		cout << endl;
		cout << "Initialization dual problem... " << endl;
	}
	// variables
	static const char* varType = NULL;
	double obj = 0;
	double lb = -CPX_INFBOUND;
	double ub = 0.0;

	// variable u_0
	snprintf(name, NAME_SIZE, "u_%i", 0);
	char* varName = (char*) (&name[0]);
	CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType, &varName);

	ub = CPX_INFBOUND;

	// variable u
	for (int i = 1; i <= num_constraint; i++) {
		snprintf(name, NAME_SIZE, "u_%i", i);
		varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);
	}

	// variables a
	for (int i = 0; i < N; i++) {
		snprintf(name, NAME_SIZE, "a_%i", i);
		varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);
	}

	// variables b
	snprintf(name, NAME_SIZE, "b");
	varName = (char*) (&name[0]);
	CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType, &varName);

	// variables v_0
	lb = 0.0;
	snprintf(name, NAME_SIZE, "v_%i", 0);
	varName = (char*) (&name[0]);
	CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType, &varName);

	lb = -CPX_INFBOUND;

	// variables v
	for (int i = 1; i <= num_constraint; i++) {
		snprintf(name, NAME_SIZE, "v_%i", i);
		varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);
	}

}

void SecondProblem::setupSP(CEnv env, Prob lp) {

	set_up_variable(env, lp);

	// constraints -A_T * u - e_k * u_0 + a
	std::vector<std::vector<double> > temp = A;

	//change sign A matrix
	change_sign_A();

	{
		std::vector<int> idx;
		std::vector<double> coef;

		// --------------------------------------------------
		//  -A_T * u
		// --------------------------------------------------
		for (int i = 0; i < N; i++) {
			char sense = 'G';
			int matbeg = 0;
			double rhs = 0;
			int nzcnt = 0;

			int iter = 0;
			int u = 1;

			while (iter < num_constraint) {

				if (A[iter][i] != 0) {
					idx.push_back(u);
					coef.push_back(A[iter][i]);
					nzcnt++;
				}
				iter++;
				u++;
			}

			// --------------------------------------------------
			//  -e_k * u0
			// --------------------------------------------------
			if (i == k) {
				idx.push_back(0);
				coef.push_back(-1);
				nzcnt++;
			}

			// --------------------------------------------------
			//  +a_i
			// --------------------------------------------------
			idx.push_back(num_constraint + 1 + i);
			coef.push_back(1);
			nzcnt++;

			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();
		}
	}

	int v_0;

	//constraint b_T * u + u_0 * gamma - b

	{
		std::vector<int> idx;
		std::vector<double> coef;

		char sense = 'G';
		int matbeg = 0;
		double rhs = 0;
		int nzcnt = 0;
		int u = 1;

		for (int i = 0; i < num_constraint; i++) {

			if (b[i] != 0) {
				idx.push_back(u);
				coef.push_back(b[i]);
				nzcnt++;
			}
			u++;

		}

		if (gam != 0) {
			idx.push_back(0);
			coef.push_back(gam);
			nzcnt++;
		}

		// --------------------------------------------------
		//  -b
		// --------------------------------------------------
		idx.push_back(num_constraint + 1 + N);
		coef.push_back(-1);
		nzcnt++;

		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
				&matbeg, &idx[0], &coef[0], 0, 0);

		idx.clear();
		coef.clear();
		v_0 = num_constraint + N + 2;

	}

	// constraints -A_T * v - e_k * v_0 + a

	{
		std::vector<int> idx;
		std::vector<double> coef;

		// --------------------------------------------------
		//  -A_T * v
		// --------------------------------------------------
		for (int i = 0; i < N; i++) {
			char sense = 'G';
			int matbeg = 0;
			double rhs = 0;
			int nzcnt = 0;

			int iter = 0;
			int v = v_0;
			v++;

			while (iter < num_constraint) {

				if (A[iter][i] != 0) {
					idx.push_back(v);
					coef.push_back(A[iter][i]);
					nzcnt++;
				}
				iter++;
				v++;
			}

			// --------------------------------------------------
			//  -e_k * v0
			// --------------------------------------------------
			if (i == k) {
				idx.push_back(v_0);
				coef.push_back(-1);
				nzcnt++;
			}

			// --------------------------------------------------
			//  +a_i
			// --------------------------------------------------
			idx.push_back(num_constraint + 1 + i);
			coef.push_back(1);
			nzcnt++;

			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();
		}
	}

	//constraint b_T * v + v_0 * (gamma+1) - b

	{
		std::vector<int> idx;
		std::vector<double> coef;

		char sense = 'G';
		int matbeg = 0;
		double rhs = 0;
		int nzcnt = 0;
		int v = v_0;
		v++;

		// --------------------------------------------------
		//  b_T * v
		// --------------------------------------------------
		for (int i = 0; i < num_constraint; i++) {

			if (b[i] != 0) {
				idx.push_back(v);
				coef.push_back(b[i]);
				nzcnt++;
			}
			v++;

		}

		// --------------------------------------------------
		//  (gamma+1) * v0
		// --------------------------------------------------
		idx.push_back(v_0);
		coef.push_back(gam + 1);
		nzcnt++;

		// --------------------------------------------------
		//  -b
		// --------------------------------------------------
		idx.push_back(num_constraint + N + 1);
		coef.push_back(-1);
		nzcnt++;

		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
				&matbeg, &idx[0], &coef[0], 0, 0);

		idx.clear();
		coef.clear();

	}

	// A_i* a + b_i* beta + u_i + v_i  =  A_i*c + b_i*z + u_i + v_i

	{

		std::vector<int> idx;
		std::vector<double> coef;

		//A_i*a
		for (int i = 0; i < num_constraint; i++) {

			double r = 0;
			int nzcnt = 0;

			for (int iter = 0; iter < N; iter++) {
				r += temp[i][iter] * cost[iter];
				if (temp[i][iter] != 0) {
					idx.push_back(num_constraint + 1 + iter);
					coef.push_back(temp[i][iter]);
					nzcnt++;
				}

			}
			r += b[i] * min_sol + dual_varVals_P1[i] + dual_varVals_P2[i];

			//b[i]*beta
			if (b[i] != 0) {
				idx.push_back(num_constraint + 1 + N);
				coef.push_back(b[i]);
				nzcnt++;
			}

			//u_i
			idx.push_back(1 + i);
			coef.push_back(1);
			nzcnt++;

			//v_i
			idx.push_back(num_constraint + 3 + N + i);
			coef.push_back(1);
			nzcnt++;

			char sense = 'E';
			int matbeg = 0;
			double rhs = r;

			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();
		}

	}

}

void SecondProblem::evaluate_rT() {

	double sum = 0;

	// --------------------------------------------------
	//evaluate coefficients for u (A matrix + b vector)
	// --------------------------------------------------
	for (int i = 0; i < num_constraint; i++) {

		sum = 0;
		for (int j = 0; j < N; j++)
			sum += A[i][j];

		sum += b[i];

		rt.push_back(sum);

	}

	// --------------------------------------------------
	//evaluate u_0
	// --------------------------------------------------
	sum = 0;
	sum = gam - 2;
	rt.push_back(sum);

	// --------------------------------------------------
	//duplicate vector for other constraints with gamma+1
	// --------------------------------------------------
	std::vector<double> temp = rt;
	temp[temp.size() - 1] = temp[temp.size() - 1] + 3;

	rt.insert(rt.end(), temp.begin(), temp.end());

	// --------------------------------------------------
	//evaluate coefficients for b
	// --------------------------------------------------
	sum = 0;
	sum = -2;
	rt.insert(rt.begin(), sum);

	// --------------------------------------------------
	//evaluate coefficients for a
	// --------------------------------------------------
	sum = 0;
	temp.clear();
	for (int i = 0; i < N; i++)
		temp.push_back(2.0);

	rt.insert(rt.begin(), temp.begin(), temp.end());

	if (verbose)
		print_r();

}

void SecondProblem::set_solution(CEnv env, Prob lp) {

	vector<double> varibles;

//	cout << "VARIABLES SECOND PROBLEM: " << endl;
	int cur_numcols = CPXgetnumcols(env, lp);

	varibles.clear();
	varibles.resize(cur_numcols);
	CHECKED_CPX_CALL(CPXgetx, env, lp, &varibles[0], 0, cur_numcols - 1);

	int surplus;
	status = CPXgetcolname(env, lp, NULL, NULL, 0, &surplus, 0,
			cur_numcols - 1);
	int cur_colnamespace = -surplus; // the space needed to save the names

	// allocate memory
	char** cur_colname = (char **) malloc(sizeof(char *) * cur_numcols);
	char* cur_colnamestore = (char *) malloc(cur_colnamespace);

	// get the names
	CPXgetcolname(env, lp, cur_colname, cur_colnamestore, cur_colnamespace,
			&surplus, 0, cur_numcols - 1);

	//  set variables
	u0 = varibles[0];

	u.clear();

	for (int i = 1; i <= num_constraint; i++)
		u.push_back(varibles[i]);

	a.clear();

	for (int i = num_constraint + 1; i < num_constraint + 1 + N; i++)
		a.push_back(varibles[i]);

	beta = varibles[num_constraint + N + 1];

	v0 = varibles[num_constraint + N + 2];

	v.clear();

	for (int i = num_constraint + N + 3; i < 2 * num_constraint + N + 3; i++)
		v.push_back(varibles[i]);

	//create y_tilde
	y_tilde.clear();
	y_tilde = a;
	y_tilde.push_back(beta);
	y_tilde.insert(y_tilde.end(), u.begin(), u.end());
	y_tilde.push_back(u0);
	y_tilde.insert(y_tilde.end(), v.begin(), v.end());
	y_tilde.push_back(v0);

	if (verbose) {
		print_y_bar();
		print_y_tilde();
		cout << endl << "VARIABLES SECOND PROBLEM: " << endl;
		print_u0();
		print_u();
		print_a();
		print_beta();
		print_v0();
		print_v();
	}
	//add y tilde in the R set
	R.insert(y_tilde);

	// free
	free(cur_colname);
	free(cur_colnamestore);

}

void SecondProblem::solve(CEnv env, Prob lp) {

	CHECKED_CPX_CALL(CPXlpopt, env, lp);

	int stat = CPXgetstat(env, lp);

	if (stat == CPX_STAT_UNBOUNDED)
		throw std::runtime_error("Second problem are unbounded");

	if (stat != CPX_STAT_INFEASIBLE) {
		set_solution(env, lp);
	} else {
		throw std::runtime_error("Second problem has conflict!!!!!");
	}

}

void SecondProblem::step8_1(CEnv env, Prob lp) {

	this->violated_constraint = false;
	std::vector<int> idx;
	std::vector<double> coef;
	int count_constraint = 0;
	char sense = 'E';
	int matbeg = 0;
	double rhs = 0;
	int nzcnt = 0;

	// --------------------------------------------------
	// Estimation A_T * u - e_k * u_0 + a = 0
	// --------------------------------------------------

	//  A_T * u
	double sum = 0;
	for (int j = 0; j < N; j++) {
		if (!(this->violated_constraint)) {
			sum = 0;
			for (int i = 0; i < num_constraint; i++) {

				if (A[i][j] != 0) {
					sum += A[i][j] * dual_varVals_P1[i];

				}

			}

			//  -e_k * u0
			if (j == k) {
				sum -= dual_varVals_P1[num_constraint];
			}

			//  +a_i
			sum += cost[j];

			//tolerance error
			if (sum < epsilon_8_1 && sum > -epsilon_8_1)
				sum = 0.0;

			// --------------------------------------------------
			//  print respect constraint
			// --------------------------------------------------

			if (sum == 0) {

				// --------------------------------------------------
				// add new constraint A_T * u - e_k * u_0 + a = 0
				// --------------------------------------------------
				nzcnt = 0;

				//  -A_T * u
				int iter = 0;
				int u = 1;

				while (iter < num_constraint) {

					if (A[iter][j] != 0) {
						idx.push_back(u);
						coef.push_back(A[iter][j]);
						nzcnt++;
					}
					iter++;
					u++;
				}

				//  -e_k * u0
				if (j == k) {
					idx.push_back(0);
					coef.push_back(-1);
					nzcnt++;
				}

				//  +a_i
				idx.push_back(num_constraint + 1 + j);
				coef.push_back(1);
				nzcnt++;

				// --------------------------------------------------
				// add new satisfy constraint if isn't in set
				// --------------------------------------------------
				if (satisfy_constraint_list.count(count_constraint) == 0) {
					satisfy_constraint_list.insert(count_constraint);

					CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs,
							&sense, &matbeg, &idx[0], &coef[0], 0, 0);
					if (verbose)
						cout << "The constraint number " << count_constraint
								<< " add " << endl << endl;
				}
				idx.clear();
				coef.clear();
			} else {
				if (satisfy_constraint_list.count(count_constraint) != 0) {
					this->violated_constraint = true;
					if (verbose) {
						cout << "The constraint number " << count_constraint
								<< " doesn't respects the equation " << "sum = "
								<< sum << endl << endl;
						cout << "Violated constraint that before was satisfied!"
								<< endl;
					}
//				throw std::runtime_error(
//					"Violated constraint that before was satisfied!");
				}
			}
			count_constraint++;
		}
	}

	if (!(this->violated_constraint)) {
		// --------------------------------------------------
		// Estimation b_T * u + u_0 * gamma - b
		// --------------------------------------------------

		//  b_T * u
		sum = 0;
		for (int i = 0; i < num_constraint; i++) {

			if (b[i] != 0) {
				sum += b[i] * dual_varVals_P1[i];

			}
		}

		//  u_0 * gamma
		sum += dual_varVals_P1[num_constraint] * gam;

		// --------------------------------------------------
		//  -b
		// --------------------------------------------------
		sum -= min_sol;

		// --------------------------------------------------
		//  print respect constraint
		// --------------------------------------------------

		//tolerance error
		if (sum < epsilon_8_1 && sum > -epsilon_8_1)
			sum = 0.0;

		if (sum == 0) {

			// --------------------------------------------------
			// add new constraint b_T * u + u_0 * gamma - b = 0
			// --------------------------------------------------

			nzcnt = 0;
			char sense = 'E';
			int matbeg = 0;
			double rhs = 0;
			int nzcnt = 0;
			int u = 1;

			for (int i = 0; i < num_constraint; i++) {

				if (b[i] != 0) {
					idx.push_back(u);
					coef.push_back(b[i]);
					nzcnt++;
				}
				u++;

			}

			if (gam != 0) {
				idx.push_back(0);
				coef.push_back(gam);
				nzcnt++;
			}

			//  -b
			idx.push_back(num_constraint + 1 + N);
			coef.push_back(-1);
			nzcnt++;

			// --------------------------------------------------
			// add new satisfy constraint if isn't in set
			// --------------------------------------------------
			if (satisfy_constraint_list.count(count_constraint) == 0) {
				satisfy_constraint_list.insert(count_constraint);

				CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
						&matbeg, &idx[0], &coef[0], 0, 0);
				if (verbose)
					cout << "The constraint number " << count_constraint
							<< " add " << endl << endl;
			}
			idx.clear();
			coef.clear();

		} else {
			if (satisfy_constraint_list.count(count_constraint) != 0) {
				this->violated_constraint = true;
				if (verbose) {
					cout << "The constraint number " << count_constraint
							<< " doesn't respects the equation " << "sum = "
							<< sum << endl << endl;
					cout << "Violated constraint that before was satisfied!"
							<< endl;
				}
//				throw std::runtime_error(
//						"Violated constraint that before was satisfied!");
			}
		}

		count_constraint++;
	}

	// --------------------------------------------------
	// Estimation A_T * v - e_k * u_v + a = 0
	// --------------------------------------------------

	//  A_T * v
	int v_0 = num_constraint + N + 2;
	sum = 0;
	for (int j = 0; j < N; j++) {
		if (!(this->violated_constraint)) {
			sum = 0;
			for (int i = 0; i < num_constraint; i++) {

				if (A[i][j] != 0) {
					sum += A[i][j] * dual_varVals_P2[i];
				}

			}

			//  -e_k * v0

			if (j == k) {
				sum -= dual_varVals_P2[num_constraint];
			}

			//  +a_i

			sum += cost[j];

			//tolerance error
			if (sum < epsilon_8_1 && sum > -epsilon_8_1)
				sum = 0.0;

			// --------------------------------------------------
			//  print respect constraint
			// --------------------------------------------------

			if (sum == 0) {
				// --------------------------------------------------
				// add new constraint A_T * v - e_k * v_0 + a = 0
				// --------------------------------------------------

				nzcnt = 0;
				int iter = 0;
				int v = v_0;
				v++;

				while (iter < num_constraint) {

					if (A[iter][j] != 0) {
						idx.push_back(v);
						coef.push_back(A[iter][j]);
						nzcnt++;
					}
					iter++;
					v++;
				}

				//  -e_k * v0

				if (j == k) {
					idx.push_back(v_0);
					coef.push_back(-1);
					nzcnt++;
				}

				//  +a_i

				idx.push_back(num_constraint + 1 + j);
				coef.push_back(1);
				nzcnt++;

				// --------------------------------------------------
				// add new satisfy constraint if isn't in set
				// --------------------------------------------------
				if (satisfy_constraint_list.count(count_constraint) == 0) {
					satisfy_constraint_list.insert(count_constraint);

					CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs,
							&sense, &matbeg, &idx[0], &coef[0], 0, 0);
					if (verbose)
						cout << "The constraint number " << count_constraint
								<< " add " << endl << endl;
				}
				idx.clear();
				coef.clear();

			} else {

				if (satisfy_constraint_list.count(count_constraint) != 0) {
					this->violated_constraint = true;
					if (verbose) {
						cout << "The constraint number " << count_constraint
								<< " doesn't respects the equation " << "sum = "
								<< sum << endl << endl;

						cout << "Violated constraint that before was satisfied!"
								<< endl;
					}
//				throw std::runtime_error(
//						"Violated constraint that before was satisfied!");
				}
			}
			count_constraint++;
		}
	}

	if (!(this->violated_constraint)) {
		// --------------------------------------------------
		// Estimation b_T * v + v_0 * (gamma+1) - b = 0
		// --------------------------------------------------

		//  b_T * v
		sum = 0;
		for (int i = 0; i < num_constraint; i++) {

			if (b[i] != 0) {
				sum += b[i] * dual_varVals_P2[i];
			}
		}

		//  v_0 * (gamma + 1)

		sum += dual_varVals_P2[num_constraint] * (gam + 1);

		//  -b
		sum -= min_sol;

		// --------------------------------------------------
		//  print respect constraint
		// --------------------------------------------------
		//tolerance error
		if (sum < epsilon_8_1 && sum > -epsilon_8_1)
			sum = 0.0;

		if (sum == 0) {

			// --------------------------------------------------
			// add new constraint b_T * v + v_0 * (gamma+1) - b = 0
			// --------------------------------------------------

			nzcnt = 0;
			int v = v_0;
			v++;

			//  b_T * v
			for (int i = 0; i < num_constraint; i++) {

				if (b[i] != 0) {
					idx.push_back(v);
					coef.push_back(b[i]);
					nzcnt++;
				}
				v++;

			}

			//  (gamma+1) * v0
			idx.push_back(v_0);
			coef.push_back(gam + 1);
			nzcnt++;

			//  -b
			idx.push_back(num_constraint + N + 1);
			coef.push_back(-1);
			nzcnt++;

			// --------------------------------------------------
			// add new satisfy constraint if isn't in set
			// --------------------------------------------------
			if (satisfy_constraint_list.count(count_constraint) == 0) {
				satisfy_constraint_list.insert(count_constraint);

				CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
						&matbeg, &idx[0], &coef[0], 0, 0);
				if (verbose)
					cout << "The constraint number " << count_constraint
							<< " add " << endl << endl;
			}
			idx.clear();
			coef.clear();
		} else {
			if (satisfy_constraint_list.count(count_constraint) != 0) {
				this->violated_constraint = true;
				if (verbose) {
					cout << "The constraint number " << count_constraint
							<< " doesn't respects the equation " << "sum = "
							<< sum << endl << endl;
					cout << "Violated constraint that before was satisfied!"
							<< endl;
				}

//				throw std::runtime_error(
//						"Violated constraint that before was satisfied!");
			}
		}

		count_constraint++;
	}

	if (!(this->violated_constraint)) {
		// --------------------------------------------------
		// Estimation -u_0>=0
		// --------------------------------------------------
		sum = 0;
		sum -= dual_varVals_P1.back();

		//tolerance error
		if (sum < epsilon_8_1 && sum > -epsilon_8_1)
			sum = 0.0;

		if (sum == 0) {
			nzcnt = 1;
			idx.push_back(0);
			coef.push_back(-1.0);

			// --------------------------------------------------
			// add new satisfy constraint if isn't in set
			// --------------------------------------------------
			if (satisfy_constraint_list.count(count_constraint) == 0) {
				satisfy_constraint_list.insert(count_constraint);

				CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
						&matbeg, &idx[0], &coef[0], 0, 0);
				if (verbose)
					cout << "The constraint number " << count_constraint
							<< " add " << endl << endl;
			}
			idx.clear();
			coef.clear();
		} else {
			if (satisfy_constraint_list.count(count_constraint) != 0) {
				this->violated_constraint = true;
				if (verbose) {
					cout << "The constraint number " << count_constraint
							<< " doesn't respects the equation " << "sum = "
							<< sum << endl << endl;
					cout << "Violated constraint that before was satisfied!"
							<< endl;
				}
//      			throw std::runtime_error(
//						"Violated constraint that before was satisfied!");
			}
		}

		count_constraint++;
	}

	if (!(this->violated_constraint)) {
		// --------------------------------------------------
		// Estimation v_0>=0
		// --------------------------------------------------
		sum = 0;
		sum += dual_varVals_P2.back();

		//tolerance error
		if (sum < epsilon_8_1 && sum > -epsilon_8_1)
			sum = 0.0;

		if (sum == 0) {
			nzcnt = 1;
			idx.push_back(v_0);
			coef.push_back(1.0);

			// --------------------------------------------------
			// add new satisfy constraint if isn't in set
			// --------------------------------------------------
			if (satisfy_constraint_list.count(count_constraint) == 0) {
				satisfy_constraint_list.insert(count_constraint);

				CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
						&matbeg, &idx[0], &coef[0], 0, 0);
				if (verbose)
					cout << "The constraint number " << count_constraint
							<< " add " << endl << endl;
			}
			idx.clear();
			coef.clear();
		} else {
			if (satisfy_constraint_list.count(count_constraint) != 0) {
				this->violated_constraint = true;
				if (verbose) {
					cout << "The constraint number " << count_constraint
							<< " doesn't respects the equation " << "sum = "
							<< sum << endl << endl;
					cout << "Violated constraint that before was satisfied!"
							<< endl;
				}
//				throw std::runtime_error(
//						"Violated constraint that before was satisfied!");
			}
		}
	}
}

void SecondProblem::step8_2(CEnv env, Prob lp) {

	std::vector<int> idx;
	std::vector<double> coef;

	if (verbose)
		cout << "r_T" << "  " << "y_bar" << endl;
	// --------------------------------------------------
	//evaluate right side r * y
	// --------------------------------------------------

	double r = 0;
	int i = 0;

	while (i < N) {
		r += rt[i] * cost[i];
		if (verbose)
			cout << rt[i] << "  " << cost[i] << endl;
		i++;
	}

	r += rt[i] * min_sol;
	if (verbose)
		cout << rt[i] << "  " << min_sol << endl;
	i++;

	int j = 0;
	while (i < N + 1 + num_constraint) {
		r += rt[i] * dual_varVals_P1[j];
		if (verbose)
			cout << rt[i] << "  " << dual_varVals_P1[j] << endl;

		i++;
		j++;
	}

	r += rt[i] * dual_varVals_P1[j];
	if (verbose)
		cout << rt[i] << "  " << dual_varVals_P1[j] << endl;

	i++;

	j = 0;
	while (i < N + 2 + num_constraint + num_constraint) {
		r += rt[i] * dual_varVals_P2[j];
		if (verbose)
			cout << rt[i] << "  " << dual_varVals_P2[j] << endl;

		i++;
		j++;
	}

	r += rt[i] * dual_varVals_P2[j];
	if (verbose)
		cout << rt[i] << "  " << dual_varVals_P2[j] << endl;

	// --------------------------------------------------
	// add constraint r_T * y = r_T * y
	// --------------------------------------------------
	char sense = 'E';
	int matbeg = 0;
	double rhs = r;
	int nzcnt = 0;

// a
	int p = 0;
	while (p < N) {
		idx.push_back(num_constraint + 1 + p);
		coef.push_back(rt[p]);
		nzcnt++;
		p++;
	}

	//b
	idx.push_back(num_constraint + 1 + N);
	coef.push_back(rt[p]);
	nzcnt++;
	p++;

	//u

	for (int iter = 1; iter <= num_constraint; iter++) {
		idx.push_back(iter);
		coef.push_back(rt[p]);
		nzcnt++;
		p++;
	}

	//u_0
	idx.push_back(0);
	coef.push_back(rt[p]);
	nzcnt++;
	p++;

	//v

	int v_0 = num_constraint + N + 2;
	int v = v_0;
	v++;
	for (int iter = 0; iter < num_constraint; iter++) {
		idx.push_back(v);
		coef.push_back(rt[p]);
		nzcnt++;
		p++;
		v++;
	}

	//v_0
	idx.push_back(v_0);
	coef.push_back(rt[p]);
	nzcnt++;

	CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense, &matbeg,
			&idx[0], &coef[0], 0, 0);

	idx.clear();
	coef.clear();

}

bool SecondProblem::y_tilde_EQ_y_bar() {

	bool equal = true;
	double difference;

	// --------------------------------------------------
	// 1. c-a
	// --------------------------------------------------

	for (unsigned int i = 0; i < a.size(); i++) {
		difference = cost[i] - a[i];

		//tolerance error
		if (difference < epsilon_8_3 && difference > -epsilon_8_3)
			difference = 0.0;

		if (difference != 0) {
			equal = false;
			return equal;
		}
	}

	// --------------------------------------------------
	// 2. z-b
	// --------------------------------------------------

	difference = min_sol - beta;

	//tolerance error
	if (difference < epsilon_8_3 && difference > -epsilon_8_3)
		difference = 0.0;

	if (difference != 0) {
		equal = false;
		return equal;
	}

	// --------------------------------------------------
	// 3. u-u
	// --------------------------------------------------

	for (unsigned int i = 0; i < u.size(); i++) {
		difference = dual_varVals_P1[i] - u[i];

		//tolerance error
		if (difference < epsilon_8_3 && difference > -epsilon_8_3)
			difference = 0.0;

		if (difference != 0) {
			equal = false;
			return equal;
		}
	}

	// --------------------------------------------------
	// 4. u_0-u_0
	// --------------------------------------------------

	difference = dual_varVals_P1.back() - u0;

	//tolerance error
	if (difference < epsilon_8_3 && difference > -epsilon_8_3)
		difference = 0.0;

	if (difference != 0) {
		equal = false;
		return equal;
	}

	// --------------------------------------------------
	// 5. v-v
	// --------------------------------------------------

	for (unsigned int i = 0; i < v.size(); i++) {
		difference = dual_varVals_P2[i] - v[i];

		//tolerance error
		if (difference < epsilon_8_3 && difference > -epsilon_8_3)
			difference = 0.0;

		if (difference != 0) {
			equal = false;
			return equal;
		}
	}

	// --------------------------------------------------
	// 6. v_0-v_0
	// --------------------------------------------------

	difference = dual_varVals_P2.back() - v0;

	//tolerance error
	if (difference < epsilon_8_3 && difference > -epsilon_8_3)
		difference = 0.0;

	if (difference != 0) {
		equal = false;
		return equal;
	}

	if (verbose)
		cout << endl;
	return equal;

}
