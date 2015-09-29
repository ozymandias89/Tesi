/**
 * @file cpx_macro.h
 * CPLEX Helper Macros
 *
 */
//ddd
#ifndef CPX_MACRO_H
#define CPX_MACRO_H

#include <cstring>
#include <string>
#include <stdexcept>
#include <ilcplex/cplex.h>

#define STRINGIZE(something) STRINGIZE_HELPER(something)
#define STRINGIZE_HELPER(something) #something

/**
 * typedefs of basic Callable Library entities,
 * i.e. environment (Env) and problem pointers (Prob).
 */

typedef CPXENVptr Env;
typedef CPXCENVptr CEnv;
typedef CPXLPptr Prob;
typedef CPXCLPptr CProb;

/* CPLEX Error Status and Message Buffer */

extern int status;

const unsigned int BUF_SIZE = 4096;

extern char errmsg[BUF_SIZE];

/* Shortcut for declaring a CPLEX Env */
#define DECL_ENV(envPtr)\
Env envPtr = CPXopenCPLEX(&status);\
if (status) {\
    CPXgeterrorstring(NULL, status, errmsg);\
    int trailer = std::strlen(errmsg) - 1;\
    if (trailer >= 0)\
        errmsg[trailer] = '\0';\
    throw std::runtime_error(std::string(__FILE__) + ":" + STRINGIZE(__LINE__) + ": " + errmsg);\
}

/* Shortcut for declaring a CPLEX Problem */
#define DECL_PROB(cEnvPtr, lpPtr, problName)\
Prob lpPtr = CPXcreateprob(cEnvPtr, &status, std::string(problName).c_str());\
if (status) {\
    CPXgeterrorstring(NULL, status, errmsg);\
    int trailer = std::strlen(errmsg) - 1;\
    if (trailer >= 0)\
        errmsg[trailer] = '\0';\
    throw std::runtime_error(std::string(__FILE__) + ":" + STRINGIZE(__LINE__) + ": " + errmsg);\
}

/* Make a checked call to a CPLEX API function */
#define CHECKED_CPX_CALL(func, env, ...)\
do {\
    status = func(env, __VA_ARGS__);\
    if (status) {\
        CPXgeterrorstring(env, status, errmsg);\
        int trailer = std::strlen(errmsg) - 1;\
        if (trailer >= 0)\
            errmsg[trailer] = '\0';\
        throw std::runtime_error(std::string(__FILE__) + ":" + STRINGIZE(__LINE__) + ": " + errmsg);\
    }\
} while(false)

#endif /* CPX_MACRO_H */
