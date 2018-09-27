//-----------------------------------------------------------------------------
//   swmm5.h
//
//   Project: EPA SWMM5
//   Version: 5.1
//   Date:    03/24/14  (Build 5.1.001)
//   Author:  L. Rossman
//
//   Prototypes for SWMM5 functions exported to swmm5.dll.
//
//-----------------------------------------------------------------------------
#ifndef SWMM5_H
#define SWMM5_H

/*
// --- define WINDOWS

#undef WINDOWS
#ifdef _WIN32
  #define WINDOWS
#endif
#ifdef __WIN32__
  #define WINDOWS
#endif

// --- define DLLEXPORT

#ifdef WINDOWS
  #define DLLEXPORT __declspec(dllexport) __stdcall
#else
#define DLLEXPORT
#endif

// --- use "C" linkage for C++ programs

*/

#ifdef __cplusplus
extern "C" {
#endif

int swmm_run(char* f1, char* f2, char* f3);
int swmm_open(char* f1, char* f2, char* f3);
int swmm_start(int saveFlag);
int swmm_step(double* elapsedTime);
int swmm_end(void);
int swmm_report(void);
int swmm_getMassBalErr(float* runoffErr, float* flowErr,
                 float* qualErr);
int swmm_close(void);
int swmm_getVersion(void);

// Cosimulation getters
double  swmm_get( char* id, int attribute, int units );
double  swmm_get_from_input(char* filename, char *id, int attribute);
int  swmm_save_all(char* input_file, int object_type, int attribute);
// Cosimulation setters
int  swmm_modify_setting(char* id, double new_setting, double tstep);
int  swmm_modify_input(char* input_file, char *id, int attribute, double value);
int  swmm_save_results();

#ifdef __cplusplus 
}   // matches the linkage specification from above */ 
#endif

#endif
