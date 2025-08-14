#include "body.h"

/* BODY DEFINITION */

const double masses[9] = {1.989E30, 3.301E23, 4.868E24, 5.972E24, 6.417E23, 1.898E27, 5.683E26, 8.681E25, 1.024E26};
const char *names[9] = {"sun", "mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune"};
const char *id[9] = {"010", "199", "299", "399", "499", "599", "699", "799", "899"};

void init_body(Body *body, double mass, double *body_vectors, long unsigned int hist_size) {
	// Initialize body parameters
	body->mass = mass;
	body->x = body_vectors[0];
	body->y = body_vectors[1];
	body->z = body_vectors[2];
	body->vx = body_vectors[3];
	body->vy = body_vectors[4];
	body->vz = body_vectors[5];

	// Allocate position history array
	body->hist = calloc(hist_size + 1, sizeof(double *));
	if (body->hist == NULL) {
		printf("ERROR: Could not allocate memory for body's history\n");
		exit(EXIT_FAILURE);
	}
	
	// Allocate memory for each history vector
	for (long unsigned int i = 0; i <= hist_size; i++) {
		body->hist[i] = calloc(6, sizeof(double));
		if (body->hist[i] == NULL) {
			printf("ERROR: Could not allocate memory for body's history vector\n");
			exit(EXIT_FAILURE);
		}
	}
}

/* SYSTEM DEFINITION */

void init_system(System *system, Method method, int year, unsigned int size, double step, double total_time) {
    // Initialize system parameters
    system->size = size + 1; // Add one for the Sun
    system->step = step;
    system->method = method;

    // Initialize total steps based on step size
    system->total_steps = (long unsigned int)(total_time / step);
	system->total_time = total_time;

    // Allocate memory for bodies
    system->bodies = calloc(system->size, sizeof(Body));
    if (system->bodies == NULL) {
        printf("ERROR: Could not allocate memory for bodies\n");
        exit(EXIT_FAILURE);
    }

	// Create curl object
	curl_global_init(CURL_GLOBAL_DEFAULT);
	CURL *curl = curl_easy_init();
	memory chunk;
	if (curl == NULL) {
		printf("ERROR: Could not initialize curl\n");
		exit(EXIT_FAILURE);
	}
	

	// Set curl options
	curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_memory);
	curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *)&chunk);

    // Initialize each body
	double *all_body_vectors[system->size];
	for (unsigned int i = 0; i < system->size; i++) {
		// Allocate memory for vectors
		all_body_vectors[i] = calloc(6, sizeof(double));
		if (all_body_vectors[i] == NULL) {
			printf("ERROR: Could not allocate memory for body vectors\n");
			exit(EXIT_FAILURE);
		}

		// Get body vectors from Horizons API and initialize
		get_vectors(curl, &chunk, all_body_vectors[i], year, id[i]);
		init_body(&system->bodies[i], masses[i], all_body_vectors[i], system->total_steps);
	}

	// Free all allocated memory
	for (unsigned int i = 0; i < system->size; i++) {
		free(all_body_vectors[i]);
	}

	// Show that initialization is complete
	printf("System initialized\n");
}

void add_hist(System *system, long unsigned int t) {
	// Add position to history
	for (unsigned int i = 0; i < system->size; i++) {
		(system->bodies[i]).hist[t][0] = (system->bodies[i]).x;
		(system->bodies[i]).hist[t][1] = (system->bodies[i]).y;
		(system->bodies[i]).hist[t][2] = (system->bodies[i]).z;
	}
}

void set_acceleration(System *system) {
    // Reset acceleration for all bodies to zero (no accumulation)
	for (unsigned int i = 0; i < system->size; i++) {
		system->bodies[i].ax = 0.0;
		system->bodies[i].ay = 0.0;
		system->bodies[i].az = 0.0;
	}

    // Compute acceleration for all pairs of bodies
	for (unsigned int i = 0; i < system->size - 1; i++) {
		for (unsigned int j = i + 1; j < system->size; j++) {
			// Compute delta vector
			double dx, dy, dz;
			dx = DIST * (system->bodies[j].x - system->bodies[i].x);
			dy = DIST * (system->bodies[j].y - system->bodies[i].y);
			dz = DIST * (system->bodies[j].z - system->bodies[i].z);

			// Square root of vector
			double r = sqrt((dx * dx) + (dy * dy) + (dz * dz));

			// Compute acceleration and assign to each body
			double mag = (ACC * G) / (r * r * r);
			double m1 = system->bodies[j].mass;
			double m2 = system->bodies[i].mass * -1.0;

			// Compute acceleration vector of first body
			system->bodies[i].ax += mag * m1 * dx;
			system->bodies[i].ay += mag * m1 * dy;
			system->bodies[i].az += mag * m1 * dz;

			// Compute acceleration vector of second body
			system->bodies[j].ax += mag * m2 * dx;
			system->bodies[j].ay += mag * m2 * dy;
			system->bodies[j].az += mag * m2 * dz;
		}
	}
}

void euler(System *system, long unsigned int t) {
	// Compute acceleration for all bodies
	set_acceleration(system);

	// Compute new velocity and position
	for (unsigned int i = 0; i < system->size; i++) {
		// Compute new velocity
		system->bodies[i].vx = system->bodies[i].vx + (system->step * system->bodies[i].ax);
		system->bodies[i].vy = system->bodies[i].vy + (system->step * system->bodies[i].ay);
		system->bodies[i].vz = system->bodies[i].vz + (system->step * system->bodies[i].az);

		// Compute new position
		system->bodies[i].x = system->bodies[i].x + (system->step * system->bodies[i].vx);
		system->bodies[i].y = system->bodies[i].y + (system->step * system->bodies[i].vy);
		system->bodies[i].z = system->bodies[i].z + (system->step * system->bodies[i].vz);
	}

	// Add position to history
	add_hist(system, t);
}

void verlet(System *system, long unsigned int t) {
	// Compute position at half-step
	for (unsigned int i = 0; i < system->size; i++) {
		system->bodies[i].x = system->bodies[i].x + (0.5 * system->step * system->bodies[i].vx);
		system->bodies[i].y = system->bodies[i].y + (0.5 * system->step * system->bodies[i].vy);
		system->bodies[i].z = system->bodies[i].z + (0.5 * system->step * system->bodies[i].vz);
	}

	// Compute acceleration for all bodies
	set_acceleration(system);

	// Compute next half-step
	for (unsigned int i = 0; i < system->size; i++) {
		// Compute new velocity
		system->bodies[i].vx = system->bodies[i].vx + (system->step * system->bodies[i].ax);
		system->bodies[i].vy = system->bodies[i].vy + (system->step * system->bodies[i].ay);
		system->bodies[i].vz = system->bodies[i].vz + (system->step * system->bodies[i].az);

		// Compute new position at the end of the step
		system->bodies[i].x = system->bodies[i].x + (0.5 * system->step * system->bodies[i].vx);
		system->bodies[i].y = system->bodies[i].y + (0.5 * system->step * system->bodies[i].vy);
		system->bodies[i].z = system->bodies[i].z + (0.5 * system->step * system->bodies[i].vz);
	}

	// Add position to history
	add_hist(system, t);
}

// RK4 variables
double *y_1, *y_2, *y_3, *y_4, *y_n, *tmp, *k1, *k2, *k3, *k4;

void f(System *system, double *y) {
	// Assign new position
	for (unsigned int i = 0; i < system->size; i++) {
		system->bodies[i].x = y[6 * i + 0];
		system->bodies[i].y = y[6 * i + 1];
		system->bodies[i].z = y[6 * i + 2];
	}

	// Compute acceleration for all bodies
	set_acceleration(system);

	// Set values for tmp array (acts as return)
	for (unsigned int i = 0; i < system->size; i++) {
		// Velocity
		tmp[6 * i + 0] = y[6 * i + 3];
		tmp[6 * i + 1] = y[6 * i + 4];
		tmp[6 * i + 2] = y[6 * i + 5];

		// Acceleration
		tmp[6 * i + 3] = system->bodies[i].ax;
		tmp[6 * i + 4] = system->bodies[i].ay;
		tmp[6 * i + 5] = system->bodies[i].az;
	}
}

void rk4(System *system, long unsigned int t) {
	// Assign values of y1
	for (unsigned int i = 0; i < system->size; i++) {
		y_1[6 * i + 0] = system->bodies[i].x;
		y_1[6 * i + 1] = system->bodies[i].y;
		y_1[6 * i + 2] = system->bodies[i].z;
		y_1[6 * i + 3] = system->bodies[i].vx;
		y_1[6 * i + 4] = system->bodies[i].vy;
		y_1[6 * i + 5] = system->bodies[i].vz;
	}

	// Compute k1
	for (unsigned int i = 0; i < system->size; i++) {
		f(system, y_1);
		for (int j = 0; j < 6; j++) k1[6 * i + j] = tmp[6 * i + j];
	}

	// Set y2 = y + step * k1 / 2 (using tmp)
	for (unsigned int i = 0; i < system->size; i++) {
		for (int j = 0; j < 6; j++) y_2[6 * i + j] = y_1[6 * i + j] + (0.5 * system->step * k1[6 * i + j]);
	}

	// Compute k2
	for (unsigned int i = 0; i < system->size; i++) {
		f(system, y_2);
		for (int j = 0; j < 6; j++) k2[6 * i + j] = tmp[6 * i + j];
	}

	// Set y3 = y + step * k2 / 2 (using tmp)
	for (unsigned int i = 0; i < system->size; i++) {
		for (int j = 0; j < 6; j++) y_3[6 * i + j] = y_1[6 * i + j] + (0.5 * system->step * k2[6 * i + j]);
	}

	// Compute k3
	for (unsigned int i = 0; i < system->size; i++) {
		f(system, y_3);
		for (int j = 0; j < 6; j++) k3[6 * i + j] = tmp[6 * i + j];
	}

	// Set y4 = y + step * k3 (using tmp)
	for (unsigned int i = 0; i < system->size; i++) {
		for (int j = 0; j < 6; j++) y_4[6 * i + j] = y_1[6 * i + j] + (system->step * k3[6 * i + j]);
	}

	// Compute k4
	for (unsigned int i = 0; i < system->size; i++) {
		f(system, y_4);
		for (int j = 0; j < 6; j++) k4[6 * i + j] = tmp[6 * i + j];
	}

	// Compute weighted average
	for (unsigned int i = 0; i < system->size; i++) {
		for (int j = 0; j < 6; j++) y_n[6 * i + j] = y_1[6 * i + j] + (system->step / 6.0) * (k1[6 * i + j] + (2.0 * k2[6 * i + j]) + (2.0 * k3[6 * i + j]) + k4[6 * i + j]);
	}

	// Set new position, velocity, and history
	for (unsigned int i = 0; i < system->size; i++) {
		// New position
		system->bodies[i].x = y_n[6 * i + 0];
		system->bodies[i].y = y_n[6 * i + 1];
		system->bodies[i].z = y_n[6 * i + 2];

		// New velocity
		system->bodies[i].vx = y_n[6 * i + 3];
		system->bodies[i].vy = y_n[6 * i + 4];
		system->bodies[i].vz = y_n[6 * i + 5];
	}

	// Add position to history
	add_hist(system, t);
}

void run(System *system) {
	// Variables for progress update
	clock_t start_time = clock();

	switch (system->method) {
		case EULER:
			for (long unsigned int t = 0; t < system->total_steps; t++) {
				update_progress(t, system->total_steps, &start_time);
				euler(system, t);
			}
			break;
		case VERLET:
			for (long unsigned int t = 0; t < system->total_steps; t++) {
				update_progress(t, system->total_steps, &start_time);
				verlet(system, t);
			}
			break;
		case RK4:
			// Initialize RK4 variables
			y_1 = calloc(system->size * 6, sizeof(double));
			y_2 = calloc(system->size * 6, sizeof(double));
			y_3 = calloc(system->size * 6, sizeof(double));
			y_4 = calloc(system->size * 6, sizeof(double));
			y_n = calloc(system->size * 6, sizeof(double));
			tmp = calloc(system->size * 6, sizeof(double));

			k1 = calloc(system->size * 6, sizeof(double));
			k2 = calloc(system->size * 6, sizeof(double));
			k3 = calloc(system->size * 6, sizeof(double));
			k4 = calloc(system->size * 6, sizeof(double));

			// Check for memory allocation errors
			if (y_1 == NULL || y_2 == NULL || y_3 == NULL || y_4 == NULL || y_n == NULL) {
				printf("ERROR: No more memory to allocate\n");
				exit(EXIT_FAILURE);
			}
			if (tmp == NULL || k1 == NULL || k2 == NULL || k3 == NULL || k4 == NULL) {
				printf("ERROR: No more memory to allocate\n");
				exit(EXIT_FAILURE);
			}

			// Run RK4 for each time step
			for (long unsigned int t = 0; t < system->total_steps; t++) {
				update_progress(t, system->total_steps, &start_time);
				rk4(system, t);
			}

			// Free RK4 variables
			free(y_1);
			free(y_2);
			free(y_3);
			free(y_4);
			free(y_n);
			free(tmp);
			free(k1);
			free(k2);
			free(k3);
			free(k4);
			break;
		/*
		case PEFRL:
			for (long unsigned int t = 0; t < system->total_steps; t++) {
				update_progress(t, system->total_steps, &start_time);
				pefrl(system, t);
			}
			break;
		*/
		default:
			printf("ERROR: Unknown method\n");
			exit(EXIT_FAILURE);
	}
}

void update_progress(long unsigned int t, long unsigned int total_t, clock_t *start_time) {
		// Show progress and estimate time left
        double percent = 100.0 * (double)(t + 1) / (double)total_t;
        clock_t current_time = clock();
        double elapsed = (double)(current_time - *start_time) / CLOCKS_PER_SEC;
        double estimated_total = elapsed / (percent / 100.0);
        double time_left = estimated_total - elapsed;

        // Print progress
        printf("\rRendering %.1f%% | Elapsed: %.1fs | Left: %.1fs | Total: %.1fs", percent, elapsed, time_left, estimated_total);
        fflush(stdout);
}

void write_data(System *system, FILE *json) {
	// Write to json file all results and extra details
	fprintf(json, "{\n\"n_steps\": %lu,\n\"step\": %lf,\n\"n\": %u,\n\"planets\": [\n", system->total_steps, system->step, system->size);
		for (unsigned int i = 0; i < system->size; i++) {
			fprintf(json, "\t{\n\t\t\"name\": \"%s\",\n\t\t\"pos\": [\n", names[i]);

			// Write all positions of body
			for (long unsigned int t = 0; t < system->total_steps; t++) {
				fprintf(json, "\t\t\t[%lf, %lf, %lf],\n", system->bodies[i].hist[t][0], system->bodies[i].hist[t][1], system->bodies[i].hist[t][2]);
			}

			// Write last position without comma
			fprintf(json, "\t\t\t[%lf, %lf, %lf]\n", system->bodies[i].hist[system->total_steps][0], system->bodies[i].hist[system->total_steps][1], system->bodies[i].hist[system->total_steps][2]);
			if (i < system->size - 1) fprintf(json, "\t\t]\n\t},\n");
			else fprintf(json, "\t\t]\n\t}\n");
		}
		fprintf(json, "]\n}");

	fclose(json);

	printf("\nFinished simulation\n");
}
