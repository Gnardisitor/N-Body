#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <curl/curl.h>

// Required constants
#define G		6.67e-11
#define DIST	1.496e11
#define TIME	86400.0
#define ACC		(TIME * TIME) / DIST

// Constants for bodies
const double masses[9] = {1.989e30, 3.301e23, 4.868e24, 5.972e24, 6.417e23, 1.898e27, 5.683e26, 8.681e25, 1.024e26};
const char *id[9] = {"010", "199", "299", "399", "499", "599", "699", "799", "899"};

// Required array variables
long unsigned int total_time;
unsigned int size;
double step;
char *algo;

// Curl variables
CURL *curl;
CURLcode result;
double body_vars[6];
int year;
char url[250];
char *pos[7];
char var[50];

// Required RK4 variables
double *y_1, *y_2, *y_3, *y_4, *y_n, *tmp, *k1, *k2, *k3, *k4;

// Create body struct
typedef struct {
	// Mass
	double mass;

	// Previous vectors
	double **hist;

	// Current vectors
	double x, y, z;
	double vx, vy, vz;
	double ax, ay, az;
} Body;

// Function to write data
struct MemoryStruct {
  char *memory;
  size_t size;
};

static size_t WriteMemoryCallback(void *contents, size_t size, size_t nmemb, void *userp) {
	size_t realsize = size * nmemb;
	struct MemoryStruct *mem = (struct MemoryStruct *)userp;
	char *ptr = realloc(mem->memory, mem->size + realsize + 1);
	if(!ptr) {
    	printf("not enough memory (realloc returned NULL)\n");
     return 0;
	}

	mem->memory = ptr;
	memcpy(&(mem->memory[mem->size]), contents, realsize);
	mem->size += realsize;
	mem->memory[mem->size] = 0;

	return realsize;
}

int get_body_vars(unsigned int body) {
	// Create curl object
	int a;
	curl = curl_easy_init();
	if (curl == NULL) {
		printf("An error has occured!\n");
		exit(1);
	}

	// Create chunk
	struct MemoryStruct chunk;
	chunk.memory = malloc(1);
	chunk.size = 0;

	// Get correct url
	//strcpy(url, "https://ssd.jpl.nasa.gov/api/horizons.api?format=text&COMMAND='399'&CENTER='@0'&EPHEM_TYPE='VECTOR'&VEC_TABLE='2'&OUT_UNITS='AU-D'&START_TIME='2000-01-01'&STOP_TIME='2000-01-02'&STEP_SIZE='2%20d'");
	sprintf(url, "https://ssd.jpl.nasa.gov/api/horizons.api?format=text&COMMAND='%s'&CENTER='@0'&EPHEM_TYPE='VECTOR'&VEC_TABLE='2'&OUT_UNITS='AU-D'&START_TIME='%d-01-01'&STOP_TIME='2000-01-02'&STEP_SIZE='2%%20d'", id[body], year);

	// Set curl operation
	curl_easy_setopt(curl, CURLOPT_URL, url);
	//curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
	curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteMemoryCallback);
	curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *)&chunk);

	result = curl_easy_perform(curl);
	if (result != CURLE_OK) {
		fprintf(stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror(result));
		return 1;
	}

	// Find all required positions in string
	pos[0] = strstr(chunk.memory, "X =");
	pos[1] = strstr(chunk.memory, "Y =");
	pos[2] = strstr(chunk.memory, "Z =");
	pos[3] = strstr(chunk.memory, "VX=");
	pos[4] = strstr(chunk.memory, "VY=");
	pos[5] = strstr(chunk.memory, "VZ=");
	pos[6] = strstr(chunk.memory, "$$EOE");

	for (int i = 0; i < 6; i++) {
		strncpy(var, (pos[i] + 3), pos[i + 1] - pos[i] + 3);
		body_vars[i] = atof(var);
	}

	curl_easy_cleanup(curl);
	free(chunk.memory);
	return 0;
}

// Create body array
Body *init_bodies(unsigned int wanted_size, unsigned int wanted_time, double wanted_step) {
	curl_global_init(CURL_GLOBAL_ALL);

	// Assign size of array
	if (wanted_size < 1 || wanted_size > 8) {
		printf("Size must be between 1 and 8 planets!\n");
		return NULL;
	}
	size = ++wanted_size;

	// Allocate memory for array
	Body *array = calloc(size, sizeof(Body));
	if (array == NULL) {
		printf("Array could not be created!\n");
		return NULL;
	}

	// Assign max total_time and step in steps
	if (wanted_time < 1 || wanted_step <= 0.0) {
		printf("Time must be larger than 1 day and step must be larger than 0!\n");
		free(array);
		return NULL;
	}
	step = wanted_step;
	total_time = (long unsigned int) (((double) wanted_time) / step);

	// Initialize each body
	for (unsigned int i = 0; i < size; i++) {
		array[i].mass = masses[i];
		array[i].hist = calloc(total_time, sizeof(double *));
		if (array[i].hist == NULL) {
			printf("History array could not be created!\n");
			return NULL;
		}

		// Assign initial positions and velocities
		if (get_body_vars(i) == 1) {
			printf("Could not get initial conditions for bodies!\n");
			return NULL;
		};
		array[i].x = body_vars[0];
		array[i].y = body_vars[1];
		array[i].z = body_vars[2];
		array[i].vx = body_vars[3];
		array[i].vy = body_vars[4];
		array[i].vz = body_vars[5];
	}

	curl_global_cleanup();
	return array;
}

// Free body array
void free_bodies(Body *array) {
	// Free history arrays
	for (unsigned int i = 0; i < size; i++) {
		// Free each position history array
		for (long unsigned int j = 0; j < total_time; j++) {
			free(array[i].hist[j]);
		}

		// Free history of each array
		free(array[i].hist);
	}

	// Free body array
	free(array);
}

void add_hist(Body *array, unsigned int index, unsigned int t) {
	// Create new double array
	array[index].hist[t] = calloc(3, sizeof(double));
	if (array[index].hist[t] == NULL) {
		printf("No more memory to allocate!\n");
		exit(1);
	}

	// Add position to history
	array[index].hist[t][0] = array[index].x;
	array[index].hist[t][1] = array[index].y;
	array[index].hist[t][2] = array[index].z;
}

// Compute acceleration for bodies
void set_acc(Body *array) {
	for (unsigned int i = 0; i < size; i++) {
		array[i].ax = 0.0;
		array[i].ay = 0.0;
		array[i].az = 0.0;
	}

	for (unsigned int i = 0; i < size - 1; i++) {
		for (long unsigned int j = i + 1; j < size; j++) {
			// Compute vector
			double dx, dy, dz;
			dx = DIST * (array[j].x - array[i].x);
			dy = DIST * (array[j].y - array[i].y);
			dz = DIST * (array[j].z - array[i].z);

			// Square root of vector
			double r= sqrt((dx * dx) + (dy * dy) + (dz * dz));

			// Compute acceleration and assign to each body
			double mag = (ACC * G) / (r * r * r);
			double m1 = array[j].mass;
			double m2 = array[i].mass * -1.0;

			// Compute acceleration vector of first body
			array[i].ax += mag * m1 * dx;
			array[i].ay += mag * m1 * dy;
			array[i].az += mag * m1 * dz;

			// Compute acceleration vector of second body
			array[j].ax += mag * m2 * dx;
			array[j].ay += mag * m2 * dy;
			array[j].az += mag * m2 * dz;
		}
	}
}

// Compute next step using Euler-Cromer
void euler(Body *array, unsigned int t) {
	// Compute acceleration for all bodies
	set_acc(array);

	// Compute new velocity and position
	for (unsigned int i = 0; i < size; i++) {
		// Compute new velocity
		array[i].vx = array[i].vx + (step * array[i].ax);
		array[i].vy = array[i].vy + (step * array[i].ay);
		array[i].vz = array[i].vz + (step * array[i].az);

		// Compute new position
		array[i].x = array[i].x + (step * array[i].vx);
		array[i].y = array[i].y + (step * array[i].vy);
		array[i].z = array[i].z + (step * array[i].vz);

		// Add position to history
		add_hist(array, i, t);
	}
}

// Compute next step using position Verlet
void verlet(Body *array, unsigned int t) {
	// Compute position at half-step
	for (unsigned int i = 0; i < size; i++) {
		array[i].x = array[i].x + (0.5 * step * array[i].vx);
		array[i].y = array[i].y + (0.5 * step * array[i].vy);
		array[i].z = array[i].z + (0.5 * step * array[i].vz);
	}

	// Compute acceleration for all bodies
	set_acc(array);

	// Compute next half-step
	for (unsigned int i = 0; i < size; i++) {
		// Compute new velocity
		array[i].vx = array[i].vx + (step * array[i].ax);
		array[i].vy = array[i].vy + (step * array[i].ay);
		array[i].vz = array[i].vz + (step * array[i].az);

		// Compute new position at the end of the step
		array[i].x = array[i].x + (0.5 * step * array[i].vx);
		array[i].y = array[i].y + (0.5 * step * array[i].vy);
		array[i].z = array[i].z + (0.5 * step * array[i].vz);

		// Add position to history
		add_hist(array, i, t);
	}
}

// Allocate memory for RK4 arrays
void init_rk4(void) {
	// Allocate memory for rows
	y_1 = calloc(size * 6, sizeof(double));
	y_2 = calloc(size * 6, sizeof(double));
	y_3 = calloc(size * 6, sizeof(double));
	y_4 = calloc(size * 6, sizeof(double));
	y_n = calloc(size * 6, sizeof(double));
	tmp = calloc(size * 6, sizeof(double));

	k1 = calloc(size * 6, sizeof(double));
	k2 = calloc(size * 6, sizeof(double));
	k3 = calloc(size * 6, sizeof(double));
	k4 = calloc(size * 6, sizeof(double));

	// Check for memory allocation errors
	if (y_1 == NULL || y_2 == NULL || y_3 == NULL || y_4 == NULL || y_n == NULL || tmp == NULL || k1 == NULL || k2 == NULL || k3 == NULL || k4 == NULL) {
		printf("No more memory to allocate!\n");
		exit(1);
	}
}

// Free allocated RK4 arrays
void free_rk4(void) {
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
}

void f(Body *array, double *y) {
	// Assign new position
	for (unsigned int i = 0; i < size; i++) {
		array[i].x = y[6 * i + 0];
		array[i].y = y[6 * i + 1];
		array[i].z = y[6 * i + 2];
	}

	// Compute acceleration for all bodies
	set_acc(array);

	// Set values for tmp array (acts as return)
	for (unsigned int i = 0; i < size; i++) {
		// Velocity
		tmp[6 * i + 0] = y[6 * i + 3];
		tmp[6 * i + 1] = y[6 * i + 4];
		tmp[6 * i + 2] = y[6 * i + 5];

		// Acceleration
		tmp[6 * i + 3] = array[i].ax;
		tmp[6 * i + 4] = array[i].ay;
		tmp[6 * i + 5] = array[i].az;
	}
}

void rk4(Body *array, unsigned int t) {
	// Assign values of y1
	for (unsigned int i = 0; i < size; i++) {
		y_1[6 * i + 0] = array[i].x;
		y_1[6 * i + 1] = array[i].y;
		y_1[6 * i + 2] = array[i].z;
		y_1[6 * i + 3] = array[i].vx;
		y_1[6 * i + 4] = array[i].vy;
		y_1[6 * i + 5] = array[i].vz;
	}

	// Compute k1
	for (unsigned int i = 0; i < size; i++) {
		f(array, y_1);
		for (int j = 0; j < 6; j++) k1[6 * i + j] = tmp[6 * i + j];
	}

	// Set y2 = y + step * k1 / 2 (using tmp)
	for (unsigned int i = 0; i < size; i++) {
		for (int j = 0; j < 6; j++) y_2[6 * i + j] = y_1[6 * i + j] + (0.5 * step * k1[6 * i + j]);
	}

	// Compute k2
	for (unsigned int i = 0; i < size; i++) {
		f(array, y_2);
		for (int j = 0; j < 6; j++) k2[6 * i + j] = tmp[6 * i + j];
	}

	// Set y3 = y + step * k2 / 2 (using tmp)
	for (unsigned int i = 0; i < size; i++) {
		for (int j = 0; j < 6; j++) y_3[6 * i + j] = y_1[6 * i + j] + (0.5 * step * k2[6 * i + j]);
	}

	// Compute k3
	for (unsigned int i = 0; i < size; i++) {
		f(array, y_3);
		for (int j = 0; j < 6; j++) k3[6 * i + j] = tmp[6 * i + j];
	}

	// Set y4 = y + step * k3 (using tmp)
	for (unsigned int i = 0; i < size; i++) {
		for (int j = 0; j < 6; j++) y_4[6 * i + j] = y_1[6 * i + j] + (step * k3[6 * i + j]);
	}

	// Compute k4
	for (unsigned int i = 0; i < size; i++) {
		f(array, y_4);
		for (int j = 0; j < 6; j++) k4[6 * i + j] = tmp[6 * i + j];
	}

	// Compute weighted average
	for (unsigned int i = 0; i < size; i++) {
		for (int j = 0; j < 6; j++) y_n[6 * i + j] = y_1[6 * i + j] + (step / 6.0) * (k1[6 * i + j] + (2.0 * k2[6 * i + j]) + (2.0 * k3[6 * i + j]) + k4[6 * i + j]);
	}

	// Set new position, velocity, and history
	for (unsigned int i = 0; i < size; i++) {
		// New position
		array[i].x = y_n[6 * i + 0];
		array[i].y = y_n[6 * i + 1];
		array[i].z = y_n[6 * i + 2];

		// New velocity
		array[i].vx = y_n[6 * i + 3];
		array[i].vy = y_n[6 * i + 4];
		array[i].vz = y_n[6 * i + 5];

		// Add position to history
		add_hist(array, i, t);
	}
}

// Main function
int main(int argc, char **argv) {
		// Check for correct argument count
		if (argc != 6) {
			printf("USAGE: ./nbody BODIES YEAR TIME STEP ALGORITHM\n");
			return 1;
		}
		// Create body array
		year = atoi(argv[2]);
		Body *bodies = init_bodies(atoi(argv[1]), atoi(argv[3]), atof(argv[4]));
		if (bodies == NULL) {
			printf("There was an error creating the body array!\n");
			return 1;
		}

		// Iterate using wanted algorithm
		if (strcmp(argv[5], "euler") == 0) {
			for (long unsigned int t = 0; t < total_time; t++) euler(bodies, t);
			algo = "euler.csv";
		}
		else if (strcmp(argv[5], "verlet") == 0) {
			for (long unsigned int t = 0; t < total_time; t++) verlet(bodies, t);
			algo = "verlet.csv";
		}
		else if (strcmp(argv[5], "rk4") == 0) {
			init_rk4();
			for (long unsigned int t = 0; t < total_time; t++) rk4(bodies, t);
			free_rk4();
			algo = "rk4.csv";
		}
		else {
			printf("Must have a valid algorithm for use!\n");
			free_bodies(bodies);
			return 1;
		}

		printf("Finished simulation!\n");

		// Create csv and write down all the
		FILE *csv = fopen(algo, "wt");
		if (csv == NULL) {
			printf("Could not open file!\n");
			return 1;
		}

		// Print one total_time step at a total_time
		for (long unsigned int t = 0; t < total_time; t++) {
			for (unsigned int i = 0; i < size; i++) {
				fprintf(csv, "%lf,%lf,%lf,", bodies[i].hist[t][0], bodies[i].hist[t][1], bodies[i].hist[t][2]);
			}
			fprintf(csv, "\n");
		}

		// Free all memory
		fclose(csv);
		free_bodies(bodies);

		return 0;
}
