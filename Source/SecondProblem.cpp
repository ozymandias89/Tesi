/*
 * SecondProblem.cpp
 *
 *  Created on: 04 nov 2015
 *      Author: riccardo
 */

#include "SecondProblem.h"

SecondProblem::SecondProblem() {
	// TODO Auto-generated constructor stub

}

SecondProblem::~SecondProblem() {
	// TODO Auto-generated destructor stub
}

void SecondProblem::print_u() {

	cout << endl;
	for (unsigned int i = 0; i < u.size(); ++i)
		cout << "u_" << i+1 << " " << " = " << u[i] << endl;

}

void SecondProblem::print_v() {

	cout << endl;
	for (unsigned  int i = 0; i < v.size(); ++i)
		cout << "v_" << i+1 << " = " << v[i] << endl;

}

void SecondProblem::print_a() {
	cout << endl;
	for (unsigned  int i = 0; i < a.size(); ++i)
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

void SecondProblem::setupSP(CEnv env, Prob lp) {

	{
		cout << endl;
		cout << "DUAL PROBLEM: " << endl;
		// variables
		static const char* varType = NULL;
		double obj = 1.0;
		double lb = -CPX_INFBOUND;
		double ub = 0.0;

		// variable u_0
		snprintf(name, NAME_SIZE, "u_%i", 0);
		char* varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);

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
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);

		// variables v_0
		lb = 0.0;
		snprintf(name, NAME_SIZE, "v_%i", 0);
		varName = (char*) (&name[0]);
		CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
				&varName);

		lb = -CPX_INFBOUND;

		// variables v
		for (int i = 1; i <= num_constraint; i++) {
			snprintf(name, NAME_SIZE, "v_%i", i);
			varName = (char*) (&name[0]);
			CHECKED_CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, varType,
					&varName);
		}

	}

	// constraints -A_T * u - e_k * u_0 + a
	std::vector<std::vector<double> > temp = A;

	//change sign matrix A
	for (unsigned int i = 0; i < A.size(); i++) {
		for (unsigned int j = 0; j < A[i].size(); j++) {
			if (A[i][j] != 0)
				A[i][j] = -A[i][j];
		}
	}

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
		if (gam != 0) {
			idx.push_back(v_0);
			coef.push_back(gam + 1);
			nzcnt++;
		}

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
				r += temp[i][iter] * c[iter];
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

std::vector<double> SecondProblem::evaluate_rT() {

	std::vector<double> rt;
	double sum = 0;

	// --------------------------------------------------
	//evaluate coefficients for u (A matrix + b vector)
	// --------------------------------------------------
	for (int i = 0; i < num_constraint; i++) {

		sum = 0;
		for (int j = 0; j < N; j++)
			sum += A[i][j];

		sum +=b[i];

		rt.push_back(sum);

	}

	// --------------------------------------------------
	//evaluate u_0
	// --------------------------------------------------
	sum = 0;
	sum = gam -2;
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
	for (int i=0;i<N;i++)
		temp.push_back(2.0);

	rt.insert(rt.begin(), temp.begin(), temp.end());

	cout << endl;
	cout << "vector r " << endl;
	for (std::vector<double>::const_iterator j = rt.begin(); j != rt.end(); ++j)
		cout << *j << " ";
	cout << endl;

	return rt;
}

void SecondProblem::set_solution(CEnv env, Prob lp) {

	vector<double> varibles;

	cout << "VARIABLES SECOND PROBLEM: " << endl;
	int cur_numcols = CPXgetnumcols(env, lp);

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

	for (int i = 1; i <= num_constraint; i++)
		u.push_back(varibles[i]);

	for (int i = num_constraint + 1; i < num_constraint + 1 + N; i++)
		a.push_back(varibles[i]);

	beta = varibles[num_constraint + N + 1];

	v0 = varibles[num_constraint + N + 2];

	for (int i = num_constraint + N + 3; i < 2 * num_constraint + N + 3; i++)
		v.push_back(varibles[i]);

	// free
	free(cur_colname);
	free(cur_colnamestore);

}

void SecondProblem::solve(CEnv env, Prob lp) {

	CHECKED_CPX_CALL(CPXlpopt, env, lp);
	print_objval(env, lp);

	set_solution(env, lp);

	print_u0();
	print_u();
	print_a();
	print_beta();
	print_v0();
	print_v();


}

void SecondProblem::step8_1(CEnv env, Prob lp) {

	std::vector<int> idx;
	std::vector<double> coef;

	// --------------------------------------------------
	//  A_T * u
	// --------------------------------------------------
	char sense = 'E';
	int matbeg = 0;
	double rhs = 0;
	int nzcnt = 0;

	double sum = 0;
	for (int j = 0; j < N; j++) {
		sum = 0;
		for (int i = 0; i < num_constraint; i++) {

			if (A[i][j] != 0) {
				sum += A[i][j] * u[i];
				//cout << " A[i][j] " << A[i][j] << " u[i] " << u[i] << endl;
			}

		}
		// --------------------------------------------------
		//  -e_k * u0
		// --------------------------------------------------
		if (j == k) {
			sum -= u0;
			//cout << " u0 " << u0 << endl;
		}
		// --------------------------------------------------
		//  +a_i
		// --------------------------------------------------
		sum += a[j];
		//cout << "a[j] " << a[j] << endl;

		// --------------------------------------------------
		//  FAMMI VEDERE QUALI DI QUESTI VINCOLI SODDISFANO L'EQUAZIONE
		// --------------------------------------------------
		if (sum < std::numeric_limits<double>::epsilon()
				&& sum > -std::numeric_limits<double>::epsilon()) {
			sum = 0.0;
		}

		// --------------------------------------------------
		// add new constraint A_T * u - e_k * u_0 + a = 0
		// --------------------------------------------------
		cout << endl;
		cout << "STAMPA DI SUM " << sum << " vincolo numero " << j << endl;
		if (sum >= 0) {
			cout << "Il vincolo numero " << j << " soddisfa l'equazione"
					<< endl;

			//add new constraint
			nzcnt = 0;
			// --------------------------------------------------
			//  -A_T * u
			// --------------------------------------------------
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

			// --------------------------------------------------
			//  -e_k * u0
			// --------------------------------------------------
			if (j == k) {
				idx.push_back(0);
				coef.push_back(-1);
				nzcnt++;
			}

			// --------------------------------------------------
			//  +a_i
			// --------------------------------------------------
			idx.push_back(num_constraint + 1 + j);
			coef.push_back(1);
			nzcnt++;

			cout << "VINCOLO AGGIUNTO " << endl;
			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();
		}

	}

	// --------------------------------------------------
	//  b_T * u
	// --------------------------------------------------
	sum = 0;
	for (int i = 0; i < num_constraint; i++) {

		if (b[i] != 0) {
			sum += b[i] * u[i];
			//cout << b[i] << " " << u[i] << endl;
		}
	}

	// --------------------------------------------------
	//  u_0 * gamma
	// --------------------------------------------------
	sum += u0 * gam;
	//cout << u0 << " " << gam << endl;

	// --------------------------------------------------
	//  -b
	// --------------------------------------------------
	sum -= beta;
	//cout << beta << endl;

	// --------------------------------------------------
	//  FAMMI VEDERE QUALI DI QUESTI VINCOLI SODDISFANO L'EQUAZIONE
	// --------------------------------------------------
	if (sum < std::numeric_limits<double>::epsilon()
			&& sum > -std::numeric_limits<double>::epsilon()) {
		sum = 0.0;
	}

	// --------------------------------------------------
	// add new constraint b_T * u + u_0 * gamma - b = 0
	// --------------------------------------------------
	cout << endl;
	cout << "STAMPA DI SUM " << sum << endl;

	if (sum >= 0) {
		cout << "Il vincolo con beta soddisfa l'equazione " << endl;
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

	}

	// --------------------------------------------------
	//  A_T * v
	// --------------------------------------------------
	int v_0 = num_constraint + N + 2;
	sum = 0;
	for (int j = 0; j < N; j++) {
		sum = 0;
		for (int i = 0; i < num_constraint; i++) {

			if (A[i][j] != 0) {
				sum += A[i][j] * v[i];
				cout << " A[i][j] " << A[i][j] << " v[i] " << v[i] << endl;
			}

		}
		// --------------------------------------------------
		//  -e_k * v0
		// --------------------------------------------------
		if (j == k) {
			sum -= v0;
			cout << " v0 " << v0 << endl;
		}
		// --------------------------------------------------
		//  +a_i
		// --------------------------------------------------
		sum += a[j];
		cout << "a[j] " << a[j] << endl;

		// --------------------------------------------------
		//  FAMMI VEDERE QUALI DI QUESTI VINCOLI SODDISFANO L'EQUAZIONE
		// --------------------------------------------------
		if (sum < std::numeric_limits<double>::epsilon()
				&& sum > -std::numeric_limits<double>::epsilon()) {
			sum = 0.0;
		}

		cout << endl;
		cout << "STAMPA DI SUM " << sum << " vincolo numero " << j << endl;
		if (sum >= 0) {
			cout << "Il vincolo numero " << j << " soddisfa l'equazione"
					<< endl;

			// --------------------------------------------------
			//  A_T * v
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

			// --------------------------------------------------
			//  -e_k * v0
			// --------------------------------------------------
			if (j == k) {
				idx.push_back(v_0);
				coef.push_back(-1);
				nzcnt++;
			}

			// --------------------------------------------------
			//  +a_i
			// --------------------------------------------------
			idx.push_back(num_constraint + 1 + j);
			coef.push_back(1);
			nzcnt++;

			cout << "VINCOLO AGGIUNTO " << endl;
			CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
					&matbeg, &idx[0], &coef[0], 0, 0);

			idx.clear();
			coef.clear();

		}
	}

	// --------------------------------------------------
	//  b_T * v
	// --------------------------------------------------

	sum = 0;
	for (int i = 0; i < num_constraint; i++) {

		if (b[i] != 0) {
			sum += b[i] * v[i];
			//cout << b[i] << " " << v[i] << endl;
		}
	}

	// --------------------------------------------------
	//  v_0 * (gamma + 1)
	// --------------------------------------------------

	sum += v0 * (gam + 1);
	//cout << v0 << " " << gam << endl;

	// --------------------------------------------------
	//  -b
	// --------------------------------------------------
	sum -= beta;
	//cout << beta << endl;

	// --------------------------------------------------
	//  FAMMI VEDERE QUALI DI QUESTI VINCOLI SODDISFANO L'EQUAZIONE
	// --------------------------------------------------
	if (sum < std::numeric_limits<double>::epsilon()
			&& sum > -std::numeric_limits<double>::epsilon()) {
		sum = 0.0;
	}
	cout << endl;
	cout << "STAMPA DI SUM " << sum << endl;
	if (sum >= 0) {
		cout << "Il vincolo con beta soddisfa l'equazione " << endl;
		nzcnt = 0;
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
		if (gam != 0) {
			idx.push_back(v_0);
			coef.push_back(gam + 1);
			nzcnt++;
		}

		// --------------------------------------------------
		//  -b
		// --------------------------------------------------
		idx.push_back(num_constraint + N + 1);
		coef.push_back(-1);
		nzcnt++;

		cout << "VINCOLO AGGIUNTO " << endl;
		CHECKED_CPX_CALL(CPXaddrows, env, lp, 0, 1, nzcnt, &rhs, &sense,
				&matbeg, &idx[0], &coef[0], 0, 0);

		idx.clear();
		coef.clear();
	}
}
