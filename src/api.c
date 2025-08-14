#include "api.h"

size_t write_memory(void *contents, size_t size, size_t nmemb, void *userp) {
    // Allocate more memory
    size_t realsize = size * nmemb;
    memory *mem = (memory *)userp;
    char *ptr = realloc(mem->memory, mem->size + realsize + 1);
    if(ptr == NULL) {
        printf("ERROR: Cannot allocate more memory\n");
        return 0;
    }

    // Copy contents to memory
    mem->memory = ptr;
    memcpy(&(mem->memory[mem->size]), contents, realsize);
    mem->size += realsize;
    mem->memory[mem->size] = 0;

    return realsize;
}

void get_vectors(CURL *curl, memory *chunk, double *vectors, int year, const char *id) {
    // Allocate chunk memory
    chunk->memory = malloc(1);
    if (chunk->memory == NULL) {
        printf("ERROR: Could not allocate memory for curl chunk\n");
        exit(EXIT_FAILURE);
    }
    chunk->size = 0;

    // Initialize variables for parsing
    char url[250];
    char *pos[7];
    char var[50];

    // Get correct url and curl operations
    sprintf(url, "https://ssd.jpl.nasa.gov/api/horizons.api?format=text&COMMAND='%s'&CENTER='@0'&EPHEM_TYPE='VECTOR'&VEC_TABLE='2'&OUT_UNITS='AU-D'&START_TIME='%d-01-01'&STOP_TIME='%d-01-02'&STEP_SIZE='2%%20d'", id, year, year);
    curl_easy_setopt(curl, CURLOPT_URL, url);

    // perform the request and check result
    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
        printf("ERROR: Curl failed\n");
        free(chunk->memory);
        exit(EXIT_FAILURE);
    }

    // Find all required positions in string
    pos[0] = strstr(chunk->memory, "X =");
    pos[1] = strstr(chunk->memory, "Y =");
    pos[2] = strstr(chunk->memory, "Z =");
    pos[3] = strstr(chunk->memory, "VX=");
    pos[4] = strstr(chunk->memory, "VY=");
    pos[5] = strstr(chunk->memory, "VZ=");
    pos[6] = strstr(chunk->memory, "$$EOE");

    // Extract values from the positions
    for (int i = 0; i < 6; i++) {
        strncpy(var, (pos[i] + 3), pos[i + 1] - pos[i] + 3);
        vectors[i] = atof(var);
    }

    // Cleanup chunk
    free(chunk->memory);
}
