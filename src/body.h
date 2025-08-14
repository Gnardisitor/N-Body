#ifndef BODY_H
#define BODY_H

#include <time.h>
#include <math.h>
#include "api.h"

/* BODY DEFINITION */

extern const double masses[9];
extern const char *names[9];
extern const char *id[9];

typedef struct {
	double mass;    // Mass
	double **hist;  // Previous vectors

	// Current vectors
	double x, y, z;
	double vx, vy, vz;
	double ax, ay, az;
} Body;

void init_body(Body *body, double mass, double *body_vectors, long unsigned int hist_size);

/* SYSTEM DEFINITION */

#define G       6.6743E-11
#define DIST    1.496E+11
#define TIME    86400.0
#define ACC	    4.989946524E-2 // (TIME ** 2) / DIST

typedef enum {
	EULER,
	VERLET,
	RK4,
	PEFRL
} Method;

typedef struct {
	Body *bodies;	// Array of bodies

	// Simulation parameters
	long unsigned int total_steps;
	double total_time;
	double step;
	unsigned int size;
	Method method;
} System;

void init_system(System *system, Method method, int year, unsigned int size, double step, double total_time);
void add_hist(System *system, long unsigned int t);
void set_acceleration(System *system);
void run(System *system);
void write_data(System *system, FILE *json);
void update_progress(long unsigned int t, long unsigned int total_t, clock_t *start_time);

void euler(System *system, long unsigned int t);
void verlet(System *system, long unsigned int t);
void f(System *system, double *y);
void rk4(System *system, long unsigned int t);
//void pefrl(System *system, long unsigned int t);	// Remvoved due to bugs with PEFRL method

#endif
